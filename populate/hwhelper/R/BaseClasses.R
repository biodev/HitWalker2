setGeneric("fromSample", def=function(obj,...) standardGeneric("fromSample"))
setGeneric("toGene", def=function(obj,...) standardGeneric("toGene"))

setGeneric("configure", def=function(obj,...) standardGeneric("configure"))
setGeneric("populate", def=function(obj,...) standardGeneric("populate"))

setGeneric("getMatrix", def=function(obj,...) standardGeneric("getMatrix"))
setGeneric("getAnnotation", def=function(obj,...) standardGeneric("getAnnotation"))
setGeneric("dataTypes", def=function(obj,...) standardGeneric("dataTypes"))
setGeneric("relNames", def=function(obj,...) standardGeneric("relNames"))
setGeneric("sampleEdge", def=function(obj,...) standardGeneric("sampleEdge"))
setGeneric("geneEdge", def=function(obj,...) standardGeneric("geneEdge"))


#' Shared Generics
#'
#' Methods that perform similar tasks for currently defined classes
#'
#' @name shared_generics
NULL

#' @rdname shared_generics
#' @param obj Either an object of class \code{Subject} or an object inheriting from \code{NeoData}
#' @return The name of the subject or node respectively
setGeneric("nodeName", def=function(obj) standardGeneric("nodeName"))


#' Testing Generics
#'
#' These methods are used for testing the Neo4j database and defined cypher queries as well as the web front end
#'
#' @name test_helpers
NULL

#' @rdname test_helpers
#' @param obj An object of class \code{HW2Config} or derived from \code{NeoData}
#' @param ... Additional values to pass to the method such as the name of a datatype defined in \code{obj}
setGeneric("getFrequency", def=function(obj,...) standardGeneric("getFrequency"))

#' @rdname test_helpers
#' @param samples A vector of sample names corresponding to the class in \code{obj}.
#' @param limit.to.hits A logial value indicating whether the resulting data.frame should be limited to rows that contained hits.
#' @return For \code{findHits} and classes derived from \code{NeoData} a \code{data.frame} with indicating the 'Sample' and 'Gene' as well as whether or not it was determined
#' to be a hit ('IsHit').
setGeneric("findHits", def=function(obj,...) standardGeneric("findHits"))

#' @rdname test_helpers
setGeneric("subjectAttrs", def=function(obj,...) standardGeneric("subjectAttrs"))

#' @rdname test_helpers
setGeneric("subjectSubset", def=function(obj,...) standardGeneric("subjectSubset"))




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
#' @aliases fromSample toGene sampleEdge geneEdge
#' @slot sample.edge.name The name of the edge going from sample to assay unit (e.g. SNP ID or probe ID).
#' @slot gene.edge.name The name of the edge going from assay unit to gene
#' @slot node.name The name of the node representation of the assay unit
#' @slot base.query A cypher template which provides the basic information about the assay at the gene level (see provided classes for examples).
#' @slot template.query A cypher template which forms the basis of aggregation-based queries (see provided classes for examples)
NeoData <- setClass(Class="NeoData", representation=list(sample.edge.name="character", gene.edge.name="character", node.name="character", base.query="character", template.query="character"))

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
#' @aliases populate configure dataTypes relNames
#' @slot subject A Subject object to be used as the basis for this HitWalker2 database
#' @slot data.list Named list containing the experiemntal data
#' @slot data.types A list of the form: list(seeds=seed.vec,target='target') where the names in seed and target correspond to the names in data.list
#' @slot gene.model String naming the type of gene model to be used, currently only entrez is supported.
HW2Config <- setClass(Class="HW2Config", representation=list(subject="Subject", data.list="list", data.types="list", gene.model="character"))

#' @rdname class_helpers
#' @param subject An object of class \code{Subject}
#' @param gene.model String naming the type of gene model to be used, currently only entrez is supported.
#' @param data.types A list of the form: list(seeds=seed.vec,target='target') where the names in seed and target correspond to the names specified to the function.
#' @param ... A series of name=value pairs which are to be loaded into the Neo4j database. All of the values in \code{data.types} should correspond to these names though
#' it is not required that all the datatypes specified here correspond to those needed for prioritization.
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
#' @param neo.path If \code{neo.path} is specified, the neo4j-shell executable is expected at neo.path/bin/neo4j-shell.  Otherwise it is expected to be part of your path.
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
    
    subject.info <- obj@subject@subject.info
    
    if (ncol(subject.info) > 1)
    {
        names(subject.info)[2:ncol(subject.info)] <- paste(names(subject.info)[1], names(subject.info)[2:ncol(subject.info)], sep=".")
    }
    
    load.neo4j(.data=subject.info, edge.name=NULL,commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")
    load.neo4j(.data=obj@subject@subject.to.sample, edge.name="DERIVED", commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&", merge.from=F)
    
    for(i in setdiff(seq_along(obj.list), skip))
    {
        message(paste("Starting:", i))
        fromSample(obj.list[[i]], neo.path=neo.path)
        toGene(obj.list[[i]], neo.path=neo.path, gene.model=gene.model)
    }
    
})

#' @rdname test_helpers
#' @return For \code{subjectSubset}: A vector of subject names in principal to be part of a given metanode
setMethod("subjectSubset", signature("HW2Config"), function(obj, subset, subset_type=c("Subject", "Subject_Category")){
  
    subj <- obj@subject
    
    subset_type = match.arg(subset_type)
    
    if (subset_type == "Subject_Category"){
        
        which.subset <- apply(subj@subject.info[,-c(1:2)] == subset, 1, any)
        
        use.subjs <- subj@subject.info[which.subset,1]
        
    }else{
        
        use.subjs <- subset
    }
    
    return (use.subjs)
})

#' @rdname test_helpers
#' @return For \code{subjectAttrs}: A \code{data.frame} with a 'Type' column which indicates the variable name, a 'Value' that indicates a variable's value for the subjects in question and
#'a 'Count' column which indicates the number of subjects with the particular value.
setMethod("subjectAttrs", signature("HW2Config"), function(obj, subset, subset_type=c("Subject", "Subject_Category")){
    
    use.subjs <- subjectSubset(obj, subset, subset_type)
    
    sub.info <- obj@subject@subject.info
    
    sub.info <- sub.info[sub.info[,1] %in% use.subjs,]
    
    if (ncol(sub.info) > 2){
        
        melt.sub <- melt(measure.vars=names(sub.info)[3:ncol(sub.info)], data=sub.info[,3:ncol(sub.info), drop=F],as.is=T)
        
        melt.sub$count <- 1
        
        sum.tab <- aggregate(count~variable+value, sum, data=melt.sub)
        
        names(sum.tab) <- c("Type", "Value", "Count")
        
        sum.tab$Type <- as.character(sum.tab$Type)
        
        return(sum.tab)
        
    }else{
        return(data.frame(Type=character(), Value=character(), Count=character(), stringsAsFactors=F))
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
    
    subj.dta <- obj@subject@subject.info
    
    #as subj.name is capitalized above
    
    lc.subj.name <- subj.name
    
    substr(lc.subj.name, 1, 1) <- tolower(substr(lc.subj.name, 1, 1))
    
    use.cols <- subj.dta[,names(subj.dta) %in% c(lc.subj.name, "alias") == F]
    
    subj_atts <- "{"
    
    for(i in names(use.cols)){
        
        if (which(i == names(use.cols)) > 1){
            subj_atts <- paste(subj_atts, ",")
        }
        
        subj_atts <- paste(subj_atts, "'", i, "':",
                            paste("set([", paste(paste0("'", unique(use.cols[,i])  ,"'"), collapse=",") ,"])", sep=""),
                           sep="")
    }
    
    subj_atts <- paste(subj_atts, "}")
    
    sub.list <- list(DATA_TYPES=toJSON(obj@data.types), SUBJECT=subj.name, REL_QUERY_STR=paste(base.query, collapse="\n"), BASE_QUERIES=paste(base.queries, collapse="\n\n"),
                     SUBJECT_ATTRIBUTES=subj_atts,TEMPLATE_QUERIES=paste(templ.queries, collapse="\n\n"), HIT_PARAMS=paste(hit.params, collapse=",\n"), USE_DATA=paste0("[", paste(paste0("'", dataTypes(obj) ,"'"), collapse=",") ,"]"))
    
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
#' @aliases addSamples<-
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


GeneSet_class <- setClass(Class="GeneSet", contains=c("AnnotatedMatrix", "HwHit"),
                          prototype = list(base.query='MATCH (subject:$SUBJECT$)-[d:DERIVED]-(sample)-[r:HAS_GENE_SET]-(m)-[r2:ASSIGNED_TO]-(o:EntrezID{name:{GENE}}) WHERE d.type = "Gene_Set" AND subject.name IN {SAMPLE}
                                           RETURN o.name AS gene, subject.name AS sample, "$DATA_NAME$" AS var, r.score*r2.weight AS score, r.score*r2.weight > $PAR_NAME$ AS is_hit;',
                                           template.query='MATCH (subject:$SUBJECT$)-[d:DERIVED]-()-[r:HAS_GENE_SET]-(m)-[r2:ASSIGNED_TO]-(gene:EntrezID) WHERE d.type = "Gene_Set" AND r.score*r2.weight > $PAR_NAME$ AND
                                           $$lower_coll_type$$.name IN {$$coll_type$$} WITH $$lower_ret_type$$.name AS ret_type, COUNT(DISTINCT $$lower_coll_type$$.name) AS use_coll
                                           ORDER BY use_coll DESC RETURN ret_type, use_coll',
                                           sample.edge.name="HAS_GENE_SET", gene.edge.name="ASSIGNED_TO", node.name="suppliedSymbol",
                                           
                                           default=0, direction=">", range=c(0, 100), display_name="GeneSet (Hit) Threshold"))

GeneSet <- function(mat, mapping){
  
  if(missing(mat) || is.null(mat) || all(is.na(mat)) || class(mat) != "matrix")
  {
    stop("ERROR: need to supply a matrix for mat")
  }
  
  if (missing(mapping) || is.null(mapping) || all(is.na(mapping)) || class(mapping) != "data.frame")
  {
    stop("ERROR: need to supply a mapping data.frame containing the suppliedSymbol->gene mappings")
  }
  
  if (all(c("suppliedSymbol", "gene") %in% names(mapping)) == F)
  {
    stop("ERROR: the mapping data.frame needs to have columns for both suppliedSymbol and gene")
  }
  
  if (length(intersect(mapping$suppliedSymbol, rownames(mat))) == 0)
  {
    stop("ERROR: There is no overlap between the rownames of mat and mapping$drug.  Is the matrix of the form: suppliedSymbol x sample?")
  }
  
  return(new("GeneSet", matrix=mat, mapping=mapping))
}



#' Drug Data Representation
#'
#' A basic class for representing drug panel results for a set of samples
#'
#' @aliases getMatrix getAnnotation
#' @slot matrix, a matrix containing a sinlge summary value for each drug (rows) and each sample (columns).  
#' @slot mapping A \code{data.frame} containing the mapping from drug to gene.  It should contain a 'drug' column, a 'gene' column and a 'weight' column which indicates the confidence of the mapping.
DrugMatrix_class <- setClass(Class="DrugMatrix", contains=c("HwHit", "AnnotatedMatrix"),
         prototype=list(base.query='MATCH (subject:$SUBJECT$)-[d:DERIVED]-(sample)-[r:HAS_DRUG_ASSAY]-(m)-[r2:ACTS_ON]-(o:EntrezID{name:{GENE}}) WHERE d.type = "Drug_Assay" AND subject.name IN {SAMPLE}
                        WITH subject, o, SUM(CASE WHEN r.score <= (m.median_ic50 / 5.0) THEN r2.weight ELSE -r2.weight END) AS effect_score
                        RETURN o.name AS gene, subject.name AS sample, "$DATA_NAME$" AS var, effect_score AS score, effect_score > $PAR_NAME$ AS is_hit;',
                        
                        template.query='MATCH (subject:$SUBJECT$)-[d:DERIVED]-()-[r:HAS_DRUG_ASSAY]-(m)-[r2:ACTS_ON]-(gene:EntrezID) WHERE d.type = "Drug_Assay" WITH subject, gene,
                        SUM(CASE WHEN r.score <= (m.median_ic50 / 5.0) THEN r2.weight ELSE -r2.weight END) AS effect_score WHERE effect_score > $PAR_NAME$ AND
                        $$lower_coll_type$$.name IN {$$coll_type$$} WITH $$lower_ret_type$$.name AS ret_type, COUNT(DISTINCT $$lower_coll_type$$.name) AS use_coll
                        ORDER BY use_coll DESC RETURN ret_type, use_coll',
                        
                        sample.edge.name="HAS_DRUG_ASSAY", gene.edge.name="ACTS_ON", node.name="drug",
                        
                        default=0, direction=">", range=c(-100, 100), display_name="GeneScore (Hit) Threshold"))

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

# @describeIn DrugMatrix_class The mapping from drug to gene taken from the \code{mapping} slot.
# @param gene.model The type of gene model to utilize.  Currently only 'entrez' is supported.
# setMethod("toGene", signature("DrugMatrix"), function(obj, neo.path=NULL, gene.model=c("entrez", "ensembl")){
#     
#     drug.genes <- getAnnotation(obj)
#     
#     drug.genes <- drug.genes[complete.cases(drug.genes),]
#     drug.genes <- drug.genes[,c("drug", "gene", "weight")]
#     names(drug.genes) <- c(nodeName(obj), switch(gene.model, entrez="entrezID", ensembl="ensembl"), "weight")
#     
#     load.neo4j(.data=drug.genes, edge.name=geneEdge(obj), commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")
# })


#' @rdname test_helpers
#' @param subjects A vector of subject names or a single category name depending on \code{subject_types}.
#' @param genes A vector of gene names or a single pathway name depending on \code{gene_types}.
#' @param subject_types The type of value(s) specified in the \code{subjects} parameter.
#' @param gene_types type type of value(s) specified in the \code{genes} parameter.
#'
#' @return For \code{findHits},HW2Config: Returns a \code{data.frame} with the following columns: 'Subject', 'Gene', 'IsHit' and 'Datatype'.  The IsHit column is
#' a logical vector indicating whether or not there was seen to be a hit for the given Subject, Gene and Datatype.
setMethod("findHits", signature("HW2Config"), function(obj, subjects, genes, subject_types=c("Subject", "Subject_Category"), gene_types=c("Gene", "Pathway")){
    
    require(RNeo4j)
    
    subject_types = match.arg(subject_types)
    gene_types = match.arg(gene_types)
    
    #get the sample names for the involved subjects
    
    #should probably add this to the class definition and have it populated using addSamples<-
    sample.rel.names <- c(Expression="Affy_Expression", Variants="DNASeq", GeneScore="Drug_Assay")
    
    subjects <- subjectSubset(obj, subset=subjects, subset_type=subject_types)
    
    #get the entrez IDs and display names for the involved genes
    
    graph = startGraph("http://localhost:7474/db/data/")
    
    if (gene_types == "Gene"){
      gene.name <- cypher(graph, "MATCH (n:EntrezID)-[r:REFERRED_TO]-(m) RETURN n.name AS gene, m.name AS symbol")
      genes <- gene.name$gene[gene.name$symbol %in% genes]
      
    }else{
      genes <- cypher(graph, paste0('MATCH (path:Pathway)-[:PATHWAY_CONTAINS]->(gene) WHERE path.name="',genes,'" RETURN gene.name'))[,1]
    }
    
    do.call(rbind, lapply(names(hw2.conf@data.list), function(x){
        
        print (x)
        
        subset.samples <- obj@subject@subject.to.sample[(obj@subject@subject.to.sample$type == sample.rel.names[x]) & (obj@subject@subject.to.sample[,1] %in% subjects),]
        
        hit.dta <- findHits(hw2.conf@data.list[[x]], subset.samples$sample, genes, F)
        
        if (nrow(hit.dta) == 0){
            return(data.frame(Subject=character(), Gene=character(), IsHit=logical(), Datatype=character()))
        }
        
        hit.dta$Datatype <- x
        
        #translate samples -> subjects
        sample.subject <- subset.samples[,1]
        names(sample.subject) <- subset.samples$sample
        
        hit.dta$Subject <- as.character(sample.subject[as.character(hit.dta$Sample)])
        
        hit.dta <- hit.dta[,c("Subject", "Gene", "IsHit", "Datatype")]
        
        #if any of the samples are a hit, then the entire subject is...
        hit.dta <- hit.dta[!duplicated(hit.dta),]
        
        return(hit.dta)
    }))
    
})


#' @rdname test_helpers
#' @param datatype the datatype used for the summary as specified in the configuration object (e.g. Expression).
#' @param subset The subset of the data to compute the summary over can be by subject name, category or gene as specified in \code{subset_type}.
#' @param subset_type The type of values specified in \code{subset}.
#' @return For \code{getFrequency}: A \code{data.frame} with columns indicating the subset_type and a column indicating the Frequency as a character
#' value (e.g. 50%).
setMethod("getFrequency", signature("HW2Config"), function(obj, datatype, subset, subset_type=c("Subject", "Subject_Category", "Gene")){
    
    #should probably add this to the class definition and have it populated using addSamples<-
    sample.rel.names <- c(Expression="Affy_Expression", Variants="DNASeq", GeneScore="Drug_Assay")
    
    if (subset_type %in% c("Subject", "Subject_Category")){
        use.subjs <- subjectSubset(obj, subset, subset_type)
        subset_type = "Subject"
    }
    
   
    
    subset.samples <- obj@subject@subject.to.sample[(obj@subject@subject.to.sample$type == sample.rel.names[datatype]),]
    
    get.freq(obj@data.list[[datatype]], use.subjs, subset_type, subset.samples)
    
})

get.freq <- function(obj, subset, type, sample.mapping){
    
    if (type == "Subject"){
    
        sample.subset <- sample.mapping$sample[sample.mapping[,1] %in% subset]
        
        sum.gene <- findHits(obj, sample.subset)
        
        freq.by(sum.gene, "Sample", "Gene", sample.subset, "Genes", length(subset))
    
    }else if (type == "Gene"){
        
    }else{
        stop("ERROR: Expected either Subject or Gene as a type")
    }
    
}

freq.by <- function(dta, freq.name, agg.name, subset, type, subset.length){
    
    sub.dta <- dta[dta[,freq.name] %in% subset,]
    
    init.agg <- aggregate(as.formula(paste(freq.name, agg.name, sep="~")), function(x) length(unique(x)), data=sub.dta)
    
    init.freq <- as.data.frame(table(init.agg[,freq.name]))
    
    ord.freq <- order(init.freq$Var1, decreasing=T)
        
    init.freq$Frequency <- paste0(round((as.numeric(as.character(init.freq$Var1))/subset.length)*100), "%")
    
    ret.dta <- init.freq[ord.freq,c("Freq", "Frequency")]
    
    names(ret.dta) <- c(type, "Frequency")
    
    return(ret.dta)
}

#' @rdname test_helpers
setMethod("findHits", signature("HW2exprSet"), function(obj, samples, genes=NULL, limit.to.hits=T){
    
    annot.pack <- annotation(obj@exprs)
    
    if (grepl("\\.db$", annot.pack) == F){
        annot.pack <- paste0(annot.pack, ".db")
    }
    
    require(annot.pack, character.only=T)
    require(reshape2)
    
    mapping <- select(get(annot.pack), columns="ENTREZID", featureNames(obj@exprs))
        
    if (missing(genes) || is.null(genes) || all(is.na(genes))){
      
      sub.mapping <- mapping
      
    }else{
      
      sub.mapping <- mapping[mapping$ENTREZID %in% genes,]
    }
    
    if (nrow(sub.mapping) == 0){
      
      return(data.frame(Sample=character(0), Gene=character(0), IsHit=logical(0)))
      
    }else{
      common.samples <- intersect(samples, sampleNames(obj@exprs))
      
      probeset.dta <- melt(exprs(obj@exprs)[,common.samples,drop=F], as.is=T)
      
      probeset.gene <- merge(probeset.dta, sub.mapping, by.x="Var1", by.y="PROBEID", sort=F)
      
      sum.gene <- aggregate(value~Var2+ENTREZID, max, data=probeset.gene)
      
      names(sum.gene) <- c("Sample", "Gene", "value")
      
      if (limit.to.hits){
        
        sum.gene <- sum.gene[sum.gene$value > obj@default,]
        
        sum.gene <- sum.gene[,c("Sample", "Gene")]
        
      }else{
        
        sum.gene$IsHit <- sum.gene$value > obj@default
        sum.gene <- sum.gene[,c("Sample", "Gene", "IsHit")]
      }
      
      return(sum.gene)
    }
        
    
})

#' @rdname test_helpers
setMethod("findHits", signature("CCLEMaf"), function(obj, samples, genes=NULL, limit.to.hits=T){
    
    if (missing(genes) || is.null(genes) || all(is.na(genes))){
        gene.maf <- obj@maf    
    }else{
        gene.maf <- obj@maf[obj@maf$Entrez_Gene_Id %in% genes,]
    }
    
    sub.maf <- gene.maf[gene.maf$Tumor_Sample_Barcode %in% samples,c("Tumor_Sample_Barcode", "Entrez_Gene_Id")]
    
    names(sub.maf) <- c("Sample", "Gene")
    
    if (limit.to.hits==F){
        sub.maf$IsHit <- T
    }
    
    return(sub.maf)
    
})

#' @rdname test_helpers
setMethod("findHits", signature("DrugMatrix"), function(obj, samples, genes=NULL, limit.to.hits=T){
    
    require(reshape2)
    
    if (missing(genes) || is.null(genes) || all(is.na(genes))){
        use.mapping <- obj@mapping    
    }else{
        use.mapping <- obj@mapping[obj@mapping$gene %in% genes,]
    }
    
    drug.dta <- melt(obj@matrix, as.is=T)
    
    median.ic50 <- aggregate(value~Var1, median, data=drug.dta)
    names(median.ic50)[2] <- "median"
    
    drug.dta.merge <- merge(drug.dta, median.ic50, by="Var1")
    
    drug.dta.merge$is.hit <- with(drug.dta.merge, value <= (median*.2))
    
    drug.dta.genes <- merge(drug.dta.merge, use.mapping, by.x="Var1", by.y="drug")
    
    drug.dta.genes <- drug.dta.genes[as.character(drug.dta.genes$Var2) %in% samples,]
    
    if (nrow(drug.dta.genes) == 0){
        
        if (limit.to.hits){
            return(data.frame(Sample=character(), Gene=character()))
        }else{
            return(data.frame(Sample=character(), Gene=character(), IsHit=logical()))
        }
    }
    
    drug.dta.genes$score <- with(drug.dta.genes, ifelse(is.hit, weight, -weight))
    
    drug.sum <- aggregate(score~Var2+gene, sum, data=drug.dta.genes)
    
    names(drug.sum) <- c("Sample", "Gene", "score")
    
    if (limit.to.hits){
        use.hits <- drug.sum[drug.sum$score > obj@default,]
    
        use.hits <- use.hits[,c("Sample", "Gene")]
    }else{
     
        drug.sum$IsHit <- drug.sum$score > obj@default
        
        use.hits <- drug.sum[,c("Sample", "Gene", "IsHit")]
        
    }
    
    return(use.hits)
})
