setGeneric("fromSample", def=function(obj,...) standardGeneric("fromSample"))
setGeneric("toGene", def=function(obj,...) standardGeneric("toGene"))

setGeneric("configure", def=function(obj,...) standardGeneric("configure"))

setGeneric("populate", def=function(obj,...) standardGeneric("populate"))

setGeneric("getMatrix", def=function(obj,...) standardGeneric("getMatrix"))
setGeneric("getAnnotation", def=function(obj,...) standardGeneric("getAnnotation"))
setGeneric("dataTypes", def=function(obj,...) standardGeneric("dataTypes"))
setGeneric("nodeName", def=function(obj,...) standardGeneric("nodeName"))
setGeneric("relNames", def=function(obj,...) standardGeneric("relNames"))
setGeneric("sampleEdge", def=function(obj,...) standardGeneric("sampleEdge"))
setGeneric("geneEdge", def=function(obj,...) standardGeneric("geneEdge"))

#' Class Helper Functions
#'
#' These functions are the preferred way to initialize the classes defined in \code{hwhelper}
#'
#' @name class_helpers
NULL

#' Neo4j Data for HitWalker2
#'
#' Basic class representing common components that all classes which load data into Neo4j for HitWalker2 should be derived from.
#'
#' @slot sample.edge.name The name of the edge going from sample to assay unit (e.g. SNP ID or probe ID).
#' @slot gene.edge.name The name of the edge going from assay unit to gene
#' @slot node.name The name of the node representation of the assay unit
#' @slot base.query A cypher template which provides the basic information about the assay at the gene level (see provided classes for examples).
#' @slot template.query A cypher template which forms the basis of aggregation-based queries (see provided classes for examples)
NeoData <- setClass(Class="NeoData", representation=list(sample.edge.name="character", gene.edge.name="character", node.name="character", base.query="character", template.query="character"))

#' Hit Parameter Data for HitWalker2
#'
#' Basic class representing common components for hit parameters that all classes which load data into Neo4j for HitWalker2 should be derived from.
#'
#' @slot default The default value used to determine if a datatype should be considered a 'hit' for a given sample
#' @slot direction Direction of the test
#' @slot range Valid range (vector) that the user can adjust the values towards
#' @slot display_name Name to be displayed to the user
HwHit <- setClass(Class="HwHit", representation=list(default="numeric", direction="character", range="numeric", display_name="character"))

setMethod("nodeName", signature("NeoData"), function(obj){
    return(obj@node.name)
})

#' @describeIn NeoData Extracts the name of the sample -> assay unit edge
#' @param obj An object of class \code{NeoData}
setMethod("sampleEdge", signature("NeoData"), function(obj){
    return(obj@sample.edge.name)
})

#' @describeIn NeoData Extracts the name of the  assay unit -> gene edge
setMethod("geneEdge", signature("NeoData"), function(obj){
    return(obj@gene.edge.name)
})

#' Subject-level information
#'
#' A class representing the subject and sample-level information in HitWalker2
#'
#' @slot subject.info, a \code{data.frame} containing subject level information.  The name of the first column will become the name associated with the subjects.
#' An alias column can be provided which will allow the end-user the ability to search on more than just the subject name.  It must be '&' delimited.
#' @slot subject.to.sample A \code{data.frame} containing the mapping from subject to sample, is best populated through \code{addSamples}.
#' @slot type The column to be used to describe the subject in terms of the study (e.g. 'Case', 'Wildtype', 'Tissue')
setClass(Class="Subject", representation=list(subject.info="data.frame", subject.to.sample="data.frame", type="character"))

#' HitWalker2 Configuration
#'
#' A class representing a HitWalker2 configuration
#'
#' @slot subject A Subject object to be used as the basis for this HitWalker2 database
#' @slot data.list Named list containing the experiemntal data
#' @slot data.types A list of the form: list(seeds=seed.vec,target='target') where the names in seed and target correspond to the names in data.list
#' @slot gene.model String naming the type of gene model to be used, currently only entrez is supported.
HW2Config <- setClass(Class="HW2Config", representation=list(subject="Subject", data.list="list", data.types="list", gene.model="character"))

#' @rdname class_helpers
#' @param subject An object of class \code{Subject}
#' @param gene.model String naming the type of gene model to be used, currently only entrez is supported.
#' @param data.types A list of the form: list(seeds=seed.vec,target='target') where the names in seed and target correspond to the names in data.list
makeHW2Config <- function(subject, gene.model=c("entrez", "ensembl"), data.types=NULL,...)
{
    gene.model <- match.arg(gene.model)
    
    data.list <- list(...)
    
    if (is.null(names(data.list)))
    {
        stop("ERROR: data.list needs to have names")
    }
    
    if (class(data.types) != "list" || is.null(data.types) == T || all(c('seeds', 'target') %in% names(data.types)) == F)
    {
        stop("ERROR: data.types needs to be a named list of the following form: list(seeds=seed.vec,target='target')")
    }
    
    data.names <- as.character(unlist(data.types))
    
    if (all(data.names %in% names(data.list)) == F)
    {
        stop("ERROR: All of the strings specified in data.types need to correspond to names supplied here")
    }
    
    return(new("HW2Config", subject=subject, data.list=data.list, data.types=data.types, gene.model=gene.model))
}

setGeneric("subjectName", def=function(obj,...) standardGeneric("subjectName"))

#' @describeIn HW2Config Extracts the subject names defined in a HW2Config object
#' @param obj An object of class \code{HW2Config}
setMethod("subjectName", signature("HW2Config"), function(obj){
    
    return(nodeName(obj@subject)) 
})

setMethod("dataTypes", signature("HW2Config"), function(obj){
    return(names(obj@data.list))  
})

setMethod("relNames", signature("HW2Config"), function(obj, data.type=NULL,  rel.type=c("from", "to")){
    
    if (missing(data.type) || is.null(data.type) || all(is.na(data.type)))
    {
        data.type <- names(obj@data.list)
    }
    
    rel.type <- match.arg(rel.type)
    
    return(sapply(obj@data.list[data.type], function(x) switch(rel.type, from=sampleEdge, to=geneEdge)(x)))
    
})

#' @describeIn HW2Config Populates a neo4j database.  The subject/sample info is populated first followed by the remaining entries.
#' @param skip If \code{neo.path} is specified, the neo4j-shell executable is expected at neo.path/bin/neo4j-shell.  Otherwise it is expected to be part of your path.
#' @param skip If \code{skip} is specified, it should be an integer vector indicating which of entries in the \code{data.list} slot should be skipped.
setMethod("populate", signature("HW2Config"), function(obj, neo.path=NULL, skip=NULL){
    
    if (missing(skip) || is.null(skip) || all(is.na(skip)))
    {
        skip <- integer(0)
    }else if (is.numeric(skip)){
        
        skip <- as.integer(skip)
        
    }else{
        stop("ERROR: skip needs to be an integer value or not specified at all")
    }
    
    obj.list <- obj@data.list
    
    gene.model <- obj@gene.model
    
    #populate the subject->sample info
    
    subject.info <- hw2.conf@subject@subject.info
    
    if (ncol(subject.info) > 1)
    {
        names(subject.info)[2:ncol(subject.info)] <- paste(names(subject.info)[1], names(subject.info)[2:ncol(subject.info)], sep=".")
    }
    
    load.neo4j(.data=subject.info, edge.name=NULL,commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")
    load.neo4j(.data=hw2.conf@subject@subject.to.sample, edge.name="DERIVED", commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&", merge.from=F)
    
    for(i in setdiff(seq_along(obj.list), skip))
    {
        message(paste("Starting:", i))
        fromSample(obj.list[[i]], neo.path=neo.path)
        toGene(obj.list[[i]], neo.path=neo.path, gene.model=gene.model)
    }
    
})


multi.gsub <- function(patterns, replacements, use.str)
{
    stopifnot(length(patterns) == length(replacements))

    for(i in 1:length(patterns))
    {
        use.str <- gsub(patterns[i], replacements[i], use.str, fixed=T)
    }
    
    return(use.str)
} 

#' @describeIn HW2Config Configures a HitWalker2 instance by substituting values into the template files defined in 'base.dir' and placing them into 'dest.dir'.
#' The defaults should suffice for Vagrant-based HitWalker2 instances.
#' @param base.dir Directory where the template files are located
#' @param dest.dir Directory where the complete template files should be placed.
#' @param make.graph.struct If TRUE, generate and stage a graph_struct file
setMethod("configure", signature("HW2Config"), function(obj, base.dir="/home/vagrant/HitWalker2/populate/",dest.dir="/home/vagrant/HitWalker2/network/", make.graph.struct=T){
    #copySubstitute() --which is part of Biobase
    
    if (file.exists(dest.dir) == F)
    {
        dir.create(dest.dir)
    }
    
    message("Creating config files from templates")
    
    src.files <- file.path(base.dir, c("tmpl_config.py", "tmpl_custom_functions.py"))
    
    #make sure the seeds are added to the configure file as a list
    obj@data.types$seeds <- as.list(obj@data.types$seeds)
    
    subj.name <- capwords(subjectName(obj))
    
    base.query <- paste0('MATCH (subject:',subj.name,')-[]->(sample) WITH subject, sample, CASE subject.alias WHEN null THEN [subject.name] ELSE [subject.name]+subject.alias END AS query_names WHERE ANY(x IN query_names WHERE x = "$$sample$$") WITH subject, sample ')
    
    prev.types <- character()
    
    for(i in dataTypes(obj)){
        base.query <- append(base.query, paste0('OPTIONAL MATCH (sample)-[:',relNames(obj, data.type=i, rel.type="from"),']-(res) WITH subject, sample',paste(c(" ", prev.types), collapse=" , "),', COUNT(res) AS ',i,' '))
        prev.types <- append(prev.types, i)
    }
    
    base.query <- append(base.query, paste0(' WHERE ', paste(paste0(dataTypes(obj), ' > 0'), collapse=" OR "), ' RETURN subject.name AS ',subj.name,', subject.',obj@subject@type, ' AS Type, sample.name AS Sample, ',
                                            paste(dataTypes(obj), collapse=" , "), ', CASE WHEN ', paste(paste0(unlist(obj@data.types), ' > 0'), collapse=" AND "), ' THEN 1 ELSE 0 END AS required_data'))
    
    #Add in the base queries for the seeds and target as well
    
    base.queries <- character()
    templ.queries <- character()
    hit.params <- ""
    #'conv_thresh':{'type':'numeric', 'default':1e-10, 'range':[0,1], 'comparison':'<', 'name':'Convergence Threshold'}
    for(i in dataTypes(obj))
    {
        if (inherits(obj@data.list[[i]], "HwHit"))
        {
            handler <- 'core.handle_gene_hits'
            use.params <- paste0("[['",tolower(i),"']]")
            cur.obj <- obj@data.list[[i]]
            hit.params <- append(hit.params, paste0("'",tolower(i),"':{'type':'numeric', 'default':",cur.obj@default, " , 'range':[",paste(cur.obj@range, collapse=","),"], 'comparison':'",cur.obj@direction,
                                                    "' , 'name':'",cur.obj@display_name,"'}"))
            
        }else{
            handler <- 'core.handle_gene_targets'
            use.params <- "None"
        }
        
        
        temp.query <- paste0(i,"= {'query':'", multi.gsub(c("$PAR_NAME$", "$DATA_NAME$", "$SUBJECT$"), c(paste0("{", tolower(i), "}"), i, subj.name), gsub("\n\\s+", " ", obj@data.list[[i]]@base.query, perl=T)), "', 'handler':",handler,", 'session_params':",use.params, "}")
        base.queries <- append(base.queries, temp.query)
        templ.queries <- append(templ.queries, paste0(i,"_tmpl = {'title':'$$ret_type$$s with ",i," hits for $$result$$','text':'",i,"', 'query':'", multi.gsub(c("$PAR_NAME$", "$DATA_NAME$", "$SUBJECT$"), c(paste0("{", tolower(i), "}"), i, subj.name), gsub("\n\\s+", " ", obj@data.list[[i]]@template.query, perl=T)), "', 'handler':None, 'session_params':",use.params, "}"))
    }
    
    sub.list <- list(DATA_TYPES=toJSON(obj@data.types), SUBJECT=subj.name, REL_QUERY_STR=paste(base.query, collapse="\n"), BASE_QUERIES=paste(base.queries, collapse="\n\n"),
                     TEMPLATE_QUERIES=paste(templ.queries, collapse="\n\n"), HIT_PARAMS=paste(hit.params, collapse=",\n"), USE_DATA=paste0("[", paste(paste0("'", dataTypes(obj) ,"'"), collapse=",") ,"]"))
    
    message("Staging config files")
    
    copySubstitute(src=src.files, dest=dest.dir, symbolValues=sub.list, symbolDelimiter="@", recursive=T)
    
    if (file.exists(file.path(dirname(dest.dir), "change_hw2_instance.sh")))
    {
        message("Setting up config files")
        cur.dir <- getwd()
        setwd(dirname(dest.dir))
        system("./change_hw2_instance.sh tmpl")
        setwd(cur.dir)
    }else{
        message("Cannot find 'change_hw2_instance.sh', skipping config file setup")
    }
    
    if (make.graph.struct){
        ccle.graph <- compute.graph.structure()
        make.graph.struct(ccle.graph, graph.struct.path="/var/www/hitwalker2_inst/static/network/data/graph_struct.json")
    }
    
})

setGeneric("addSamples<-", def=function(obj,..., value) standardGeneric("addSamples<-"))

#' @describeIn Subject Adds sample information to a \code{Subject} object.  This can either be done by passing
#' in a \code{data.frame} or by passing in an object which has a \code{sampleNames} method.  The \code{data.frame}
#' should contain a column referencing the subject, a 'sample' column and a 'type' column.  The 'type' column provides
#' a mechanism through which samples with the same name can be differentiated in terms of data type.  If an object is
#' supplied, the type column can be supplied as part of the method call (e.g. addSamples(subj, type="variant") <- object).
#' @param obj A object of class \code{Subject}
#' @param ... Additional argument list, currently used when \code{value} is not a \code{data.frame} to specify additional subject attributes such as 'type'.
#' @param value Either an object with a \code{sampleNames} method defined or a \code{data.frame} with at least a 'sample' and column named the same as the subject. 
setReplaceMethod("addSamples", signature("Subject"), function(obj, ..., value){
    
    if (is.data.frame(value) == F)
    {
        subj.name <- names(obj@subject.info)[1]
        
        use.names <- intersect(sampleNames(value), obj@subject.info[,1])
        
        call.list <- append(setNames(list(use.names, use.names), c(subj.name,"sample")), list(...))
        
        call.list$stringsAsFactors <- F
        
        value <- do.call(data.frame, call.list)
    }
    
    if (('sample' %in% names(value) && names(obj@subject.info)[1] %in% names(value)) == F)
    {
        stop("ERROR: the names of the supplied data.frame need to reference 'sample' as well as the subject node name")
    }
    
    if (nrow(obj@subject.to.sample) == 0)
    {
        obj@subject.to.sample <- value
        
    }else{
        
        if (all(names(obj@subject.to.sample) %in% names(value)) == F)
        {
            stop("ERROR: all the names of the current subject.to.sample need to be present in the supplied data.frame")
        }
        
        obj@subject.to.sample <- rbind(obj@subject.to.sample, value[,names(obj@subject.to.sample)])
    }
    
    validObject(obj)
    return(obj)
    
})

setMethod("nodeName", signature("Subject"), function(obj){
    return(names(obj@subject.info)[1])  
})

#' @rdname class_helpers
#' @param subject.info A \code{data.frame} with at least one column and row, where the first column name corresponds to how the subjects' should be named in HitWalker2.
#' @param subject.to.sample A \code{data.frame} describing the mapping from subject to sample, can be left as NULL here and filled in later using \code{addSamples}
#' @param type.col If specified, the column of \code{subject.info} that should be used as a description of the subject in the displays.  Defaults as 'N/A'.
Subject <- function(subject.info, subject.to.sample=NULL, type.col=NULL)
{
    if(missing(subject.info) || is.null(subject.info) || all(is.na(subject.info)))
    {
        stop("ERROR: subject.info needs to be supplied")
    }else if (class(subject.info) != "data.frame" || ncol(subject.info) < 1 || nrow(subject.info) == 0){
        stop("ERROR: subject.info needs to be a data.frame with at least one column and one row")
    }
    
    if (missing(type.col) || is.null(type.col) || all(is.na(type.col)))
    {
        subject.info$type <- "N/A"
        type <- "type"
        
    }else if (length(type.col) == 1 && is.character(type.col) == T && type.col %in% colnames(subject.info)){
        type <- type.col
    }else{
        stop("ERROR: type.col needs to be a single column indicating the 'type' of a subject (e.g. diagnosis  or type of cell line)")
    }
    
    if (missing(subject.to.sample) || is.null(subject.to.sample) || all(is.na(subject.to.sample)))
    {
        subject.to.sample <- data.frame()
    }
    
    return(new("Subject", subject.info=subject.info, subject.to.sample=subject.to.sample, type=type))
}


##Note here, that for dense datatypes like expression, don't actually return possible hits
#' Basic Representation for Expression Array Data
#'
#' A basic class for representing Affymetrix expression array data.
#'
#' @slot exprs, an \code{ExpressionSet} containing expression data.
HW2exprSet_class <- setClass(Class="HW2exprSet", representation=list(exprs="ExpressionSet"), contains=c("NeoData", "HwHit"),
         prototype=list(sample.edge.name="HAS_EXPRESSION", gene.edge.name="PS_MAPPED_TO", node.name="probeSet",
                        default=.75, direction=">", range=c(0,1), display_name="Expression (Hit) Threshold",
                        base.query='MATCH(n:$SUBJECT$)-[d:DERIVED]-()-[r:HAS_EXPRESSION]-()-[:PS_MAPPED_TO]-(m:EntrezID{name:{GENE}}) WHERE d.type = "Affy_Expression" AND HAS(r.score) AND n.name IN {SAMPLE}
                        AND r.score > $PAR_NAME$ RETURN m.name AS gene, n.name AS sample, "$DATA_NAME$" AS var, MAX(r.score) AS score, true AS is_hit',
                        
                        template.query='MATCH(subject:$SUBJECT$)-[d:DERIVED]-()-[r:HAS_EXPRESSION]-()-[:PS_MAPPED_TO]-(gene:EntrezID) WHERE d.type = "Affy_Expression" AND r.score > $PAR_NAME$ AND $$lower_coll_type$$.name IN {$$coll_type$$}
                        WITH $$lower_ret_type$$.name AS ret_type, COLLECT(DISTINCT $$lower_coll_type$$.name) AS use_coll WHERE LENGTH(use_coll) = {$$coll_type$$_length} RETURN ret_type'))

#' @rdname class_helpers
#' @param exprs An ExpressionSet
#' @param sample.edge.name The name of the relationship between the samples and the assay identifier e.g. probeset, drug or variant.
#' @param gene.edge.name The name of the relationship between the assay identifiers and genes
#' @param node.name The name that should be given to the assay units in the Neo4j database.
HW2exprSet <- function(exprs, sample.edge.name="HAS_EXPRESSION", gene.edge.name="PS_MAPPED_TO", node.name="probeSet"){
    
    return(new("HW2exprSet", exprs=exprs, sample.edge.name=sample.edge.name, gene.edge.name=gene.edge.name, node.name=node.name))
}

#' @describeIn HW2exprSet_class Implements the loading of sample to probeset data into Neo4j.  The sample names are derived from the column names and the
#' probeset names are derived from the rownames. The score attribute is derived from the values of the matrix.
#' @param obj An object of class \code{HW2exprSet}.
#' @param neo.path The optional path to a Neo4j database.
setMethod("fromSample", signature("HW2exprSet"), function(obj, neo.path=NULL){
    
    use.exprs <- exprs(obj@exprs)
    
    melt.use.exprs <- melt(use.exprs)
    names(melt.use.exprs) <- c("probeset", "sample", "score")
    
    melt.use.exprs$sample <- as.character(melt.use.exprs$sample)
    melt.use.exprs$probeset <- as.character(melt.use.exprs$probeset)
    
    #now load the sample -> probe mappings
    
    melt.use.exprs <-  melt.use.exprs[,c("sample", "probeset", "score")]
    names(melt.use.exprs)[2] <- nodeName(obj)
    
    load.neo4j(.data=melt.use.exprs, edge.name=sampleEdge(obj), commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&", unique.rels=F)
    
})

#' @describeIn HW2exprSet_class Implements the loading of the probeset to gene mapping data into Neo4j.  The mapping data is derived from the annotation
#' database specified using the \code{annotation} method for the \code{ExpressionSet}.
#' @param gene.model The type of gene model to utilize.  Currently only 'entrez' is supported.
setMethod("toGene", signature("HW2exprSet"), function(obj, neo.path=NULL,gene.model=c("entrez", "ensembl")){
    
    gene.model <- match.arg(gene.model)
    
    if (length(annotation(obj@exprs)) == 0 || require(annotation(obj@exprs), character.only=T) == F)
    {
        stop("ERROR: need to specify an annotation package ")
    }else{
        annotation.package <- annotation(obj@exprs)
    }
    
    #Then create a mapping file from probeset to gene
    
    gene.type <- switch(gene.model, entrez="ENTREZID", ensembl="ENSEMBL")
    
    probe.to.gene <- suppressWarnings(select(eval(parse(text=annotation.package)), keys=featureNames(obj@exprs), column=gene.type, keytypes="PROBEID"))
    
    #Discard for now those that do not map to either type of genes.
    
    probe.to.gene <- probe.to.gene[is.na(probe.to.gene[,gene.type]) == F,]
    
    #probesets that map to multiple genes are perhaps not that reliable either, so discard them as well for now...
    
    dup.probes <- probe.to.gene$PROBEID[duplicated(probe.to.gene$PROBEID)]
    
    probe.to.gene <- probe.to.gene[probe.to.gene$PROBEID %in% dup.probes == F,]
    names(probe.to.gene) <- c(nodeName(obj), switch(gene.model, entrez="entrezID", ensembl="ensembl"))
    
    #load the probe->gene mappings
    
    load.neo4j(.data=probe.to.gene, edge.name=geneEdge(obj), commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")
})


#' CCLE MAF Representation
#'
#' A basic class for representing the Cancer Cell Line Encyclopedia style Mutation Annotation Format files.
#'
#' @slot maf, a \code{data.frame} representation of the MAF file.
CCLEMaf_class <- setClass(Class="CCLEMaf", representation=list(maf="data.frame"), contains="NeoData",
         prototype=list(
                        base.query='MATCH (n:$SUBJECT$)-[d:DERIVED]-()-[r:HAS_DNASEQ]-(var)-[r2:IMPACTS]-(gene:EntrezID{name:{GENE}})-[:REFFERED_TO]-(symb) WHERE d.type = "DNASeq" AND n.name IN {SAMPLE}
                            RETURN var.name AS Variant_Position, r2.transcript AS Transcript, gene.name AS Gene, symb.name AS Symbol,
                            r.ref_counts as Ref_Counts, r.alt_counts AS Alt_Counts, REPLACE(RTRIM(REDUCE(str="",n IN var.dbsnp|str+n+" ")), " ", ";") AS dbSNP,
                            r2.variant_classification AS Variant_classification, r2.protein AS Protein_Change, 0 AS query_ind, 2 AS gene_ind, var.name + "_" + gene.name AS row_id, n.name AS Sample ',
                        
                        template.query='MATCH (subject:$SUBJECT$)-[d:DERIVED]-()-[r:HAS_DNASEQ]-(var)-[r2:IMPACTS]-(gene:EntrezID) WHERE d.type = "DNASeq" AND $$lower_coll_type$$.name IN {$$coll_type$$}
                        WITH $$lower_ret_type$$.name AS ret_type, COLLECT(DISTINCT $$lower_coll_type$$.name) AS use_coll WHERE LENGTH(use_coll) = {$$coll_type$$_length} RETURN ret_type'
                        ))
#' @rdname class_helpers
#' @param file.name The path to the CCLE MAF file.
readMAF.ccle <- function(file.name, sample.edge.name="HAS_DNASEQ", gene.edge.name="IMPACTS", node.name="variation")
{
    use.maf <- read.delim(file.name, sep="\t", stringsAsFactors=F)
    
    keep.maf <- use.maf[,c("Entrez_Gene_Id", "Genome_Change", "Variant_Classification", "Annotation_Transcript", "Transcript_Strand", "cDNA_Change", "Codon_Change", "Protein_Change",
                           "Tumor_Sample_Barcode", "Genome_Change", "Center", "Sequencer", "Alternative_allele_reads_count", "Reference_allele_reads_count", "dbSNP_RS", "dbSNP_Val_Status")]

    return(new("CCLEMaf", maf=keep.maf, sample.edge.name=sample.edge.name, gene.edge.name=gene.edge.name, node.name=node.name))
}

setGeneric("maf", def=function(obj,...) standardGeneric("maf"))
setMethod("maf", signature("CCLEMaf"), function(obj){
    return(obj@maf)
})

setMethod("sampleNames", signature("CCLEMaf"), function(object){
    return(unique(maf(object)$Tumor_Sample_Barcode))
})

#' @describeIn CCLEMaf_class Loads the sample to variant data into Neo4j.  The sample annotation is taken from the 'Tumor_Sample_Barcode' column and the variant annotation is taken from the 'Genome_Change' column.
#' The following additional columns are added as attributes as well:
#' \itemize{
#'   \item Center
#'   \item Sequencer
#'   \item Alternative_allele_reads_count
#'   \item Reference_allele_reads_count
#'   \item dbSNP_RS
#'   \item dbSNP_Val_Status
#' }
#' @param obj The optional path to a Neo4j database.
#' @param neo.path The optional path to a Neo4j database.
setMethod("fromSample", signature("CCLEMaf"), function(obj, neo.path=NULL){
    #first sample -> variant
    #the name here will be derived from the Genome_Change column as that provides potentially enough information to uniquely id a variant (indels might still be tricky...)
    #will keep missing values as "" for now
    
    sample.maf <- maf(obj)
    
    samp.maf <- sample.maf[,c("Tumor_Sample_Barcode", "Genome_Change", "Center", "Sequencer", "Alternative_allele_reads_count", "Reference_allele_reads_count", "dbSNP_RS", "dbSNP_Val_Status")]
    names(samp.maf) <- c("sample", nodeName(obj), "center", "sequencer", "alt_counts", "ref_counts", "variation.dbsnp", "variation.dbsnp_val_status")
    samp.maf$variation.dbsnp <- gsub(";", "&", samp.maf$variation.dbsnp)
    samp.maf$variation.dbsnp_val_status <- gsub(";", "&", samp.maf$variation.dbsnp_val_status)
    
    #also note there that things like presence in dbSNP or COSMIC etc could be used as a property in the Variation node--should add in Variant_Type here...
    
    load.neo4j(.data=samp.maf, edge.name=sampleEdge(obj), commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")
})

#' @describeIn CCLEMaf_class Loads the variant to gene data into Neo4j.  The mapping is produced using the 'Genome_Change' and 'Entrez_Gene_Id' columns.
#' Additionally, the following columns are also included as attributes:
#' \itemize{
#'      \item Variant_Classification
#'      \item Annotation_Transcript
#'      \item Transcript_Strand
#'      \item cDNA_Change
#'      \item Codon_Change
#'      \item Protein_Change
#' }
#' @param gene.model The type of gene model to utilize.  Currently only 'entrez' is supported.
setMethod("toGene", signature("CCLEMaf"), function(obj, neo.path=NULL, gene.model="entrez"){
     #then add in the Variation->EntrezID relationships
    
    if (gene.model != "entrez")
    {
        stop("ERROR: Only gene.model = 'entrez' is supported for MAF files")
    }
    
    cur.maf <- maf(obj)
    
    #only keep one row for each gene/variant
    var.gene.dta <- cur.maf[!duplicated(cur.maf[,c("Entrez_Gene_Id", "Genome_Change")]),]
    
    var.gene.dta <- var.gene.dta[,c("Genome_Change", "Entrez_Gene_Id", "Variant_Classification", "Annotation_Transcript", "Transcript_Strand", "cDNA_Change", "Codon_Change", "Protein_Change")]
    
    names(var.gene.dta) <- c(nodeName(obj), "entrezID","variant_classification", "transcript", "transcript_strand", "cdna", "codon", "protein")
    
    load.neo4j(.data=var.gene.dta, edge.name=geneEdge(obj), commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")  
    
})


#' Drug Data Representation
#'
#' A basic class for representing drug panel results for a set of samples
#'
#' @slot matrix, a matrix containing a sinlge summary value for each drug (rows) and each sample (columns).  
#' @slot mapping A \code{data.frame} containing the mapping from drug to gene.  It should contain a 'drug' column, a 'gene' column and a 'weight' column which indicates the confidence of the mapping.
DrugMatrix_class <- setClass(Class="DrugMatrix", representation=list(matrix="matrix", mapping="data.frame"), contains=c("NeoData", "HwHit"),
         prototype=list(base.query='MATCH (n:$SUBJECT$)-[d:DERIVED]-()-[r:HAS_DRUG_ASSAY]-(m)-[r2:ACTS_ON]-(o:EntrezID{name:{GENE}}) WHERE d.type = "Drug_Assay" AND n.name IN {SAMPLE}
                        WITH n, o, SUM(CASE WHEN r.score <= (m.median_ic50 / 5.0) THEN r2.weight ELSE -r2.weight END) AS effect_score
                        RETURN o.name AS gene, n.name AS sample, "$DATA_NAME$" AS var, effect_score AS score, effect_score > $PAR_NAME$ AS is_hit;',
                        
                        template.query='MATCH (subject:$SUBJECT$)-[d:DERIVED]-()-[r:HAS_DRUG_ASSAY]-(m)-[r2:ACTS_ON]-(gene:EntrezID) WHERE d.type = "Drug_Assay" WITH subject, gene,
                        SUM(CASE WHEN r.score <= (m.median_ic50 / 5.0) THEN r2.weight ELSE -r2.weight END) AS effect_score WHERE effect_score > $PAR_NAME$ AND
                        $$lower_coll_type$$.name IN {$$coll_type$$} WITH $$lower_ret_type$$.name AS ret_type, COLLECT(DISTINCT $$lower_coll_type$$.name) AS use_coll
                        WHERE LENGTH(use_coll) = {$$coll_type$$_length} RETURN ret_type',
                        
                        sample.edge.name="HAS_DRUG_ASSAY", gene.edge.name="ACTS_ON", node.name="drug",
                        
                        default=0, direction=">", range=c(-100, 100), display_name="GeneScore (Hit) Threshold"))

setMethod("getMatrix", signature("DrugMatrix"), function(obj)
          {
                return(obj@matrix)
          })

setMethod("getAnnotation", signature("DrugMatrix"), function(obj){
        return(obj@mapping)
})

setMethod("sampleNames", signature("DrugMatrix"), function(object){
    return(unique(colnames(getMatrix(object))))
})

#' @rdname class_helpers
#' @param mat A \code{matrix} of the form drug x sample with named rows and columns.
#' @param mapping A \code{data.frame} containing the mappings between drug and gene with at least column names 'drug' and 'gene'.
DrugMatrix <- function(mat, mapping,sample.edge.name="HAS_DRUG_ASSAY", gene.edge.name="ACTS_ON", node.name="drug"){
    
    if(missing(mat) || is.null(mat) || all(is.na(mat)) || class(mat) != "matrix")
    {
        stop("ERROR: need to supply a matrix for mat")
    }
    
    if (missing(mapping) || is.null(mapping) || all(is.na(mapping)) || class(mapping) != "data.frame")
    {
        stop("ERROR: need to supply a mapping data.frame containing the drug->gene mappings")
    }
    
    if (all(c("drug", "gene") %in% names(mapping)) == F)
    {
        stop("ERROR: the mapping data.frame needs to have columns for both drug and gene")
    }
    
    if (length(intersect(mapping$drug, rownames(mat))) == 0)
    {
        stop("ERROR: There is no overlap between the rownames of mat and mapping$drug.  Is the matrix of the form: drug x sample?")
    }
    
    return(new("DrugMatrix", matrix=mat, mapping=mapping, sample.edge.name=sample.edge.name, gene.edge.name=gene.edge.name, node.name=node.name))
}

#' @describeIn DrugMatrix_class The mapping from sample to drug taken from the column and rownames of the \code{matrix} slot.
#' The \code{score} attribute consists of the elements of the matrix.  In addition, a 'median_IC50' attribute is added to the
#' drug node.
#' @param obj The optional path to a Neo4j database.
#' @param neo.path The optional path to a Neo4j database. 
setMethod("fromSample", signature("DrugMatrix"), function(obj, neo.path=NULL){
    
    drug.mat <- getMatrix(obj)
    
    drug.dta <- melt(drug.mat)
    
    names(drug.dta) <- c("drug", "sample", "score")
    
    drug.dta$drug <- as.character(drug.dta[,"drug"])
    drug.dta$sample <- as.character(drug.dta$sample)
    
    #remove any na scores
    
    stopifnot(sum(is.na(drug.dta$score)) == sum(is.na(drug.dta)))
    
    drug.dta <- drug.dta[complete.cases(drug.dta),]
    
    #reorder the edges
    
    drug.dta <- drug.dta[,c("sample", "drug", "score")]
    
    #compute a median for the drugs
    drug.medians <- aggregate(score~drug, median, data=drug.dta)
    names(drug.medians) <- c(nodeName(obj), paste0(nodeName(obj),".median_ic50"))
    
    load.neo4j(.data=drug.medians, edge.name=NULL, commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")
    
    names(drug.dta)[2] <- nodeName(obj)
    
    load.neo4j(.data=drug.dta, edge.name=sampleEdge(obj), commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")
    
})

#' @describeIn DrugMatrix_class The mapping from drug to gene taken from the \code{mapping} slot.
#' @param gene.model The type of gene model to utilize.  Currently only 'entrez' is supported.
setMethod("toGene", signature("DrugMatrix"), function(obj, neo.path=NULL, gene.model=c("entrez", "ensembl")){
    
    drug.genes <- getAnnotation(obj)
    
    drug.genes <- drug.genes[complete.cases(drug.genes),]
    drug.genes <- drug.genes[,c("drug", "gene", "weight")]
    names(drug.genes) <- c(nodeName(obj), switch(gene.model, entrez="entrezID", ensembl="ensembl"), "weight")
    
    load.neo4j(.data=drug.genes, edge.name=geneEdge(obj), commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")
})
