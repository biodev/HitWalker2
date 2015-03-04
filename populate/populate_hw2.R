library(Biobase)
library(reshape2)
library(igraph)
library(rjson)

#importing pathway commons

#download .gmt file from http://www.pathwaycommons.org/pc2/downloads.html
#here it is in /Users/bottomly/Desktop/hitwalker2_paper/Pathway Commons.4.All.GSEA.gmt.gz

#to incorporate into hitwalker use the org.Hs.eg.db database to relate uniprot to entrez
#both the private and public versions will have reference to Entrez IDs so this will be useful...

setGeneric("fromSample", def=function(obj,...) standardGeneric("fromSample"))
setGeneric("toGene", def=function(obj,...) standardGeneric("toGene"))
setGeneric("getMatrix", def=function(obj,...) standardGeneric("getMatrix"))
setGeneric("getAnnotation", def=function(obj,...) standardGeneric("getAnnotation"))
setGeneric("configure", def=function(obj,...) standardGeneric("configure"))
setGeneric("populate", def=function(obj,...) standardGeneric("populate"))
setGeneric("dataTypes", def=function(obj,...) standardGeneric("dataTypes"))

setGeneric("nodeName", def=function(obj,...) standardGeneric("nodeName"))
setGeneric("relNames", def=function(obj,...) standardGeneric("relNames"))
setGeneric("sampleEdge", def=function(obj,...) standardGeneric("sampleEdge"))
setGeneric("geneEdge", def=function(obj,...) standardGeneric("geneEdge"))

setClass(Class="NeoData", representation=list(sample.edge.name="character", gene.edge.name="character", node.name="character"))

setMethod("nodeName", signature("NeoData"), function(obj){
    return(obj@node.name)
})

setMethod("sampleEdge", signature("NeoData"), function(obj){
    return(obj@sample.edge.name)
})

setMethod("geneEdge", signature("NeoData"), function(obj){
    return(obj@gene.edge.name)
})

setClass(Class="Subject", representation=list(subject.info="data.frame", subject.to.sample="data.frame", type="character"))

#data.types needs to be a list like: list(seeds=seed.vec,target='target')
#where the names in seeds and target need to correspond to the names of the ... arguments to makeHW2Config

setClass(Class="HW2Config", representation=list(subject="Subject", data.list="list", data.types="list", gene.model="character"))

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

setMethod("populate", signature("HW2Config"), function(obj, neo.path, skip=NULL){
    
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

setMethod("configure", signature("HW2Config"), function(obj, dest.dir="/Users/bottomly/Desktop/github_projects/HitWalker2/network/"){
    #copySubstitute() --which is part of Biobase
    
    if (file.exists(dest.dir) == F)
    {
        dir.create(dest.dir)
    }
    
    src.files <- paste0("/Users/bottomly/Desktop/github_projects/HitWalker2/populate/", c("tmpl_config.py", "tmpl_custom_functions.py"))
    
    #make sure the seeds are added to the configure file as a list
    obj@data.types$seeds <- as.list(obj@data.types$seeds)
    
    subj.name <- capwords(subjectName(obj))
    
    base.query <- paste0('MATCH (subject:',subj.name,')-[]->(sample) WHERE ANY(x IN [subject.name] + subject.alias WHERE x = "$$sample$$") WITH subject, sample ')
    
    prev.types <- character()
    
    for(i in dataTypes(obj)){
        base.query <- append(base.query, paste0('OPTIONAL MATCH (sample)-[:',relNames(obj, data.type=i, rel.type="from"),']-(res) WITH subject, sample',paste(c(" ", prev.types), collapse=" , "),', COUNT(res) AS ',i,' '))
        prev.types <- append(prev.types, i)
    }
    
    base.query <- append(base.query, paste0(' WHERE ', paste(paste0(dataTypes(obj), ' > 0'), collapse=" OR "), ' RETURN subject.name AS ',subj.name,', subject.',obj@subject@type, ' AS Type, sample.name AS Sample, ',
                                            paste(dataTypes(obj), collapse=" , "), ', CASE WHEN ', paste(paste0(obj@data.types$seeds, ' > 0'), collapse=" AND "), ' THEN 1 ELSE 0 END AS required_data'))
    
    #Add in the base queries for the seeds and target as well
    
    base.queries <- character()
    
    for(i in unlist(obj@data.types))
    {
        temp.query <- paste0(i,"= {'query':'", gsub("\n\\s+", " ", obj@data.list[[i]]@base.query, perl=T), "', 'handler':None, 'session_params':None", "}")
        base.queries <- append(base.queries, temp.query)
    }
    
    rep.query <- paste0("{'", obj@data.types$target, "':[{'query':'", gsub("\n\\s+", " ", obj@data.list[[obj@data.types$target]]@report.query, perl=T)  ,"', 'handler':core.handle_query_prior, 'session_params':None}]}")
    
    sub.list <- list(DATA_TYPES=toJSON(obj@data.types), SUBJECT=subj.name, REL_QUERY_STR=paste(base.query, collapse="\n"), BASE_QUERIES=paste(base.queries, collapse="\n\n"),
                     REPORT_QUERY=rep.query)
    
    copySubstitute(src=src.files, dest=dest.dir, symbolValues=sub.list, symbolDelimiter="@", recursive=T)
    
    #also might need the rjson package
    #write(toJSON())
    
})



setGeneric("addSamples<-", def=function(obj,..., value) standardGeneric("addSamples<-"))
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

#expression utils, affy for now...

setClass(Class="HW2exprSet", representation=list(exprs="ExpressionSet"), contains="NeoData")

HW2exprSet <- function(exprs, sample.edge.name="HAS_EXPRESSION", gene.edge.name="PS_MAPPED_TO", node.name="probeSet"){
    
    return(new("HW2exprSet", exprs=exprs, sample.edge.name=sample.edge.name, gene.edge.name=gene.edge.name, node.name=node.name))
}

setMethod("fromSample", signature("HW2exprSet"), function(obj, neo.path){
    
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

setMethod("toGene", signature("HW2exprSet"), function(obj, neo.path,gene.model=c("entrez", "ensembl")){
    
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

#MAF class utils

setClass(Class="CCLEMaf", representation=list(maf="data.frame", base.query="character", report.query="character"), contains="NeoData",
         prototype=list(base.query='MATCH (n:Sample)-[r:HAS_DNASEQ]-(var)-[r2:IMPACTS]-(gene) WHERE n.name IN {SAMPLE} RETURN n,r,m,r2,o',
                        report.query='MATCH (n:Sample)-[r:HAS_DNASEQ]-(var)-[r2:IMPACTS]-(gene)-[:REFFERED_TO]-(symb) WHERE n.name IN {name}
                            RETURN var.name AS Variant_Position, r2.transcript AS Transcript, gene.name AS Gene, symb.name AS Symbol,
                            r.ref_counts as Ref_Counts, r.alt_counts AS Alt_Counts, REPLACE(RTRIM(REDUCE(str="",n IN var.dbsnp|str+n+" ")), " ", ";") AS dbSNP,
                            r2.variant_classification AS Variant_classification, r2.protein AS Protein_Change, 0 AS query_ind, 2 AS gene_ind, var.name + "_" + gene.name AS row_id '))

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


setMethod("fromSample", signature("CCLEMaf"), function(obj, neo.path){
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

setMethod("toGene", signature("CCLEMaf"), function(obj, neo.path, gene.model="entrez"){
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

#Drug class utils

setClass(Class="DrugMatrix", representation=list(matrix="matrix", mapping="data.frame", base.query="character"), contains="NeoData",
         prototype=list(base.query='MATCH (n:Sample)-[r:HAS_DRUG_ASSAY]-(m)-[r2:ACTS_ON]-(o:EntrezID{name:{GENE}}) WHERE n.name IN {SAMPLE}
                        WITH n, o, SUM(CASE WHEN r.score <= (m.median_ic50 / 5.0) THEN r2.weight ELSE -r2.weight END) AS effect_score
                        RETURN o.name AS gene, n.name AS sample, "GeneScore" AS var, effect_score AS score, effect_score > 0 AS is_hit;'))

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

setMethod("fromSample", signature("DrugMatrix"), function(obj, neo.path){
    
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

setMethod("toGene", signature("DrugMatrix"), function(obj, neo.path, gene.model=c("entrez", "ensembl")){
    
    drug.genes <- getAnnotation(obj)
    
    drug.genes <- drug.genes[complete.cases(drug.genes),]
    drug.genes <- drug.genes[,c("drug", "gene", "weight")]
    names(drug.genes) <- c(nodeName(obj), switch(gene.model, entrez="entrezID", ensembl="ensembl"), "weight")
    
    load.neo4j(.data=drug.genes, edge.name=geneEdge(obj), commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")
})

read.pc.gmt <- function(filename, organism.code="9606")
{
    org.pack <- switch(make.names(organism.code), X9606="org.Hs.eg.db")
    
    require(org.pack, character.only=T)
    
    gmt.lines <- readLines(filename)
    gmt.split.lines <- strsplit(gmt.lines, "\t")
    
    path.org <- sapply(gmt.split.lines, function(x) strsplit(x[1], ":")[[1]][1])
    use.paths <- is.na(path.org) == F & path.org == organism.code
    
    use.path.lines <- gmt.split.lines[use.paths]
    #
    path.meta <- lapply(use.path.lines, function(x) c(sub(paste0(organism.code,":\\s+"), "", x[1]), regmatches(x[2], regexec("datasource:\\s+(\\w+);\\s+organism:\\s+(\\d+);\\s+id\\s+type:\\s+(\\w+)",x[2]))[[1]][-1]))
    
    path.meta.dta.list <- lapply(1:length(path.meta), function(x) cbind(matrix(path.meta[[x]], ncol=length(path.meta[[x]]), nrow=length(use.path.lines[[x]])-2, byrow=T), use.path.lines[[x]][-c(1:2)]))
    
    path.meta.dta <- data.frame(do.call("rbind", path.meta.dta.list), stringsAsFactors=F)
    names(path.meta.dta) <- c("pathway", "database", "organism", "id_type", "id")
    
    stopifnot(all(path.meta.dta$id_type == "uniprot") && all(path.meta.dta$organism == "9606"))
    
    #warnings due to multi-way relationships
    suppressWarnings(uni.to.ent <- select(eval(parse(text=org.pack)), keys=as.character(path.meta.dta$id), columns="ENTREZID", keytype="UNIPROT"))
    
    #merge them and return the result
    
    path.meta.dta.merge <- merge(path.meta.dta, uni.to.ent, by.x="id", by.y="UNIPROT", all=T, incomparables=NA)
    
    ret.dta <- path.meta.dta.merge[,c("pathway", "database", "id", "ENTREZID")]
    names(ret.dta)[3:4] <- c("uniprot", "entrezID")
    
    return(ret.dta)
}

#from toupper docs
capwords <- function(s, strict = FALSE) {
         cap <- function(s) paste(toupper(substring(s, 1, 1)),
                       {s <- substring(s, 2); if(strict) tolower(s) else s},
                                  sep = "", collapse = " " )
         sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
     }

make.property.str.node <- function(prop.names, dta, pos, array.delim)
    {
        base.str <- paste0("{name:line[",pos,"]")
        
        for(i in prop.names)
        {
            use.name <- sapply(strsplit(i, "\\."), function(x) x[-1])
            use.pos <- which(names(dta) == i)-1
            
            stopifnot(length(use.pos) == 1)
            
            if (is.integer(dta[,i]))
            {
                base.str <- paste0(base.str, ",",use.name,":toInt(line[",use.pos,"])")
            }else if (is.numeric(dta[,i]))
            {
                base.str <- paste0(base.str, ",",use.name,":toFloat(line[",use.pos,"])")
            }else{
                
                if (is.null(array.delim) == F && is.character(array.delim) == T && length(array.delim) == 1)
                {
                    if (any(grepl(array.delim, dta[,i])))
                    {
                        base.str <- paste0(base.str, ",",use.name,":split(line[",use.pos,"], '",array.delim,"')")
                    }else{
                        base.str <- paste0(base.str, ",",use.name,":line[",use.pos,"]")
                    }
                }else{
                    base.str <- paste0(base.str, ",",use.name,":line[",use.pos,"]")
                }
            }
        }
        
        base.str <- paste0(base.str, "}")
        
        return(base.str)
    }

get.or.create.constraints <- function(neo.path, node.names)
{
    #before executing, things will go A LOT faster if constraints/indexes are available, check first and if not, then create them
    
    found.schema <- system(paste(file.path(neo.path, "bin", "neo4j-shell"), "-c schema"), intern=T)
    
    constraint.section <- grep("Constraints", found.schema)
    
    if (length(constraint.section) > 0)
    {
        found.constraints <- found.schema[(constraint.section+1):length(found.schema)]
    
        present.constraints <- sapply(node.names, function(x)
               {
                    any(grepl(paste0("\\s+ON\\s+\\(",tolower(x),":",x,"\\) ASSERT\\s+", tolower(x), ".name IS UNIQUE"), found.constraints))
               })
    }else{
        #no constraints were found so set both notes to present.constraints=F
        present.constraints <- rep(F,2)
        names(present.constraints) <- node.names
    }
    
    if (all(present.constraints))
    {
        message("Found all necessary CONSTRAINTS, continuing...")
        missing.constraints <- character(0)
    }else{
        missing.constraints <- names(present.constraints)[present.constraints == F]
        message(paste("Missing CONSTRAINT(s) for:", paste(missing.constraints, collapse=","), "will build them first and continue with loading..."))
    }
    
    if (length(missing.constraints) > 0)
    {
        constr.str <- paste0(file.path(neo.path, "bin", "neo4j-shell"), " -c 'CREATE CONSTRAINT ON (",tolower(missing.constraints),":",missing.constraints,") ASSERT ",tolower(missing.constraints),".name IS UNIQUE;'")
    }else{
        constr.str <- character(0)
    }
    
    return (constr.str)
}

load.neo4j <- function(.data, edge.name=NULL, commit.size=1000, neo.path="/Users/bottomly/Desktop/resources/programs/neo4j-community-2.1.2", dry.run=F, array.delim="&", unique.rels=T, merge.from=T, merge.to=T)
{
    if (missing(edge.name) || is.null(edge.name))
    {
        load.neo4j.node(.data=.data, commit.size=commit.size, neo.path=neo.path, dry.run=dry.run, array.delim=array.delim)
        
    }else{
        load.neo4j.edge(.data=.data, edge.name=edge.name, commit.size=commit.size, neo.path=neo.path, dry.run=dry.run, array.delim=array.delim, unique.rels=unique.rels, merge.from=merge.from, merge.to=merge.to)
    }
}

load.neo4j.node <- function(.data, commit.size=1000, neo.path="/Users/bottomly/Desktop/resources/programs/neo4j-community-2.1.2", dry.run=F, array.delim="&")
{
    if (any(is.na(.data)))
    {
        stop("ERROR: .data currently cannot have NAs, please remove and try again...")
    }
    
    base.node.names <- names(.data)[1]
    
    node.names <- capwords(base.node.names)
    
    node.props.dta <- lapply(base.node.names, function(x)
                         {
                            return(names(.data)[grep(paste0(x,"\\."), names(.data))])
                         })
    
    use.temp <- tempfile()
    
    write.table(.data, sep="\t", file=use.temp, col.names=F, row.names=F, quote=F)
    
    cypher.temp <- tempfile()
    
    constr.str <- get.or.create.constraints(neo.path, node.names)
    
    cypher.stats <- c(
        paste0("USING PERIODIC COMMIT ",commit.size),
        paste0("LOAD CSV FROM 'file://", use.temp,"' AS line FIELDTERMINATOR '\t'"),
        paste0("MERGE (n:",node.names[1], make.property.str.node(node.props.dta[[1]], .data, 0, array.delim),");")
    )
    
    writeLines(cypher.stats, con=cypher.temp)
    
    if (dry.run == F)
    {
        for (i in constr.str)
        {
            system(i)
        }
        
        message("Loading into Neo4j...")
        system(paste(file.path(neo.path, "bin", "neo4j-shell"), "-file", cypher.temp))
    }else{
        
        for (i in constr.str)
        {
            message(i)
        }
        
        message(paste(file.path(neo.path, "bin", "neo4j-shell"), "-file", cypher.temp))
    }
}

#data.frame should be in the form:
#data.frame(node1, node2, node[12].property, edge property(no x'.'y just name))
load.neo4j.edge <- function(.data, edge.name=NULL, commit.size=1000, neo.path="/Users/bottomly/Desktop/resources/programs/neo4j-community-2.1.2", dry.run=F, array.delim="&", unique.rels=T, merge.from=T, merge.to=T)
{
    message("Preprocessing...")
    
    if (missing(edge.name) || is.null(edge.name) || is.na(edge.name))
    {
        stop("ERROR: Need to supply an edge name")
    }else if (is.character(edge.name) == F || length(edge.name) != 1)
    {
        stop("ERROR: edge.name needs to be a single string value")
    }
    
    if (ncol(.data) < 2 || is.null(names(.data)))
    {
        stop("ERROR: .data needs to have at least two columns and be named")
    }
    
    if (any(is.na(.data)))
    {
        stop("ERROR: .data currently cannot have NAs, please remove and try again...")
    }
    
    
    .make.property.str.edge <- function(prop.names, dta, array.delim)
    {
        if (length(prop.names) == 0)
        {
            return("")
        }
        
        base.str <- "{"
        
        for(i in prop.names)
        {
            if (which(prop.names == i) > 1)
            {
                base.str <- paste0(base.str, ",")
            }
            
            use.pos <- which(names(dta) == i)-1
            
            stopifnot(length(use.pos) == 1)
            
            if (is.integer(dta[,i]))
            {
                base.str <- paste0(base.str,i,":toInt(line[",use.pos,"])")
            }else if (is.numeric(dta[,i]))
            {
                base.str <- paste0(base.str,i,":toFloat(line[",use.pos,"])")
            }else{
                
                if (is.null(array.delim) == F && is.character(array.delim) == T && length(array.delim) == 1)
                {
                    if (any(grepl(array.delim, dta[,i])))
                    {   
                        base.str <- paste0(base.str,i,":split(line[",use.pos,"], '", array.delim,"')")
                    }else{
                        base.str <- paste0(base.str,i,":line[",use.pos,"]")
                    }
                }else{
                
                    base.str <- paste0(base.str,i,":line[",use.pos,"]")
                
                }
            }
        }
        
        base.str <- paste0(base.str, "}")
        
        return(base.str)
    }
    
    base.node.names <- names(.data)[1:2]
    
    node.names <- capwords(base.node.names)
    
    node.props.dta <- lapply(base.node.names, function(x)
                         {
                            return(names(.data)[grep(paste0(x,"\\."), names(.data))])
                         })
    
    edge.props <- setdiff(names(.data), c(base.node.names, unlist(node.props.dta)))
    
    use.temp <- tempfile()
    
    write.table(.data, sep="\t", file=use.temp, col.names=F, row.names=F, quote=F)
    
    cypher.temp <- tempfile()
    
    constr.str <- get.or.create.constraints(neo.path, node.names)
    
    rel.unique.str <- ifelse(unique.rels, "UNIQUE", "")
    
    if (merge.from){
        from_str = "MERGE"
    }else{
        from_str = "MATCH"
    }
    
    if (merge.to){
        to_str = "MERGE"
    }else{
        to_str = "MATCH"
    }
    
    cypher.stats <- c(
        paste0("USING PERIODIC COMMIT ",commit.size),
        paste0("LOAD CSV FROM 'file://", use.temp,"' AS line FIELDTERMINATOR '\t'"),
        paste0(from_str, " (n:",node.names[1], make.property.str.node(node.props.dta[[1]], .data, 0, array.delim),")"),
        paste0(to_str, " (m:",node.names[2],make.property.str.node(node.props.dta[[2]], .data, 1, array.delim),")"),
        paste0("CREATE ",rel.unique.str," (n)-[:",edge.name,.make.property.str.edge(edge.props, .data, array.delim),"]->(m);")
    )
    
    writeLines(cypher.stats, con=cypher.temp)
    
    if (dry.run == F)
    {
        for (i in constr.str)
        {
            system(i)
        }
        
        message("Loading into Neo4j...")
        system(paste(file.path(neo.path, "bin", "neo4j-shell"), "-file", cypher.temp))
    }else{
        
        for (i in constr.str)
        {
            message(i)
        }
        
        message(paste(file.path(neo.path, "bin", "neo4j-shell"), "-file", cypher.temp))
    }
}



sample.neo4j.graph <- function(graph_struct="/var/www/hitwalker_2_inst/graph_struct.json", start.node.label=NULL, num.start.nodes.sample=10, base.query="MATCH (n)-[r]-(m)", max.rels=100, dump.file="test.ndp", neo.path="/Users/bottomly/Desktop/resources/programs/neo4j-community-2.1.6")
{
    #read.graph_struct
    
    g.struct <- read.graph_struct(graph_struct)
    
    if (missing( start.node.label) || is.null( start.node.label) || all(is.na( start.node.label)))
    {
        V(g.struct)$name
    }
    
    all.edges <- get.edgelist(g.struct)
    
    previous.vars <- character(0)
    
    entire.statement <- ""
    
    for(i in 1:nrow(all.edges))
    {
        cur.e.label <- E(g.struct)[V(g.struct)[name == all.edges[i,1]] %->% V(g.struct)[name == all.edges[i,2]]]$label
        
        previous.vars <- union(previous.vars, c(tolower(all.edges[i,1]), tolower(all.edges[i,2]), tolower(cur.e.label)))
        
        entire.statement <- paste(entire.statement, paste0('MATCH (',tolower(all.edges[i,1]),':',all.edges[i,1],')-[',tolower(cur.e.label),':',cur.e.label,']-(',tolower(all.edges[i,2]),':',all.edges[i,2],') WITH ', paste(previous.vars, collapse=","), ' LIMIT ', max.rels))
    }
    
    entire.statement <- paste(entire.statement, "RETURN", paste(previous.vars, collapse=","))
    
    tmp.file <- tempfile()
    
    writeLines(paste0("DUMP", entire.statement, ";"), con=tmp.file)
    
    system(paste(file.path(neo.path, "bin", "neo4j-shell"), "-file", tmp.file,">", dump.file))
}

clean.neo4j.res <- function(result)
{
    if (length(result) == 6)
    {
        return(NULL)
    }else{
        
        use.res <- result[seq(from=4, to=length(result)-3)]
    
        #also get the header
        use.res <- append(result[2], use.res)
        
        use.res <- gsub("\"", "", gsub("|\\[|\\]|\\s", "", use.res), fixed=T)
        
        split.res <- sapply(strsplit(use.res, "\\|"), function(x) x[x!=""])
        
        if (class(split.res) == "character")
        {
            return(split.res[-1])
        }else if (class(split.res) == "matrix")
        {
            header <- split.res[,1]
            split.res <- split.res[,-1]
            ret.dta <- data.frame(t(split.res), stringsAsFactors=F)
            names(ret.dta) <- make.names(header)
            return(ret.dta)
        }else{
            stop("ERROR: Unexpected type for result")   
        }
    }
    
    
}

compute.graph.structure <- function(neo.path)
{
    message("Getting labels from DB")
    
    label.query <- "'MATCH (n) RETURN DISTINCT LABELS(n);'"
    
    labs <- system(paste(file.path(neo.path, "bin", "neo4j-shell"), "-c", label.query), intern=T)
    
    use.labs <- clean.neo4j.res(labs)
    
    lab.list <- lapply(use.labs, function(i)
           {
                message(paste("Starting", i))
                cur.query <- paste0('MATCH (n:', i ,')-[r]->(m) RETURN DISTINCT LABELS(n) AS from_node, LABELS(m) AS to_node, TYPE(r) AS type;')
        
                cur.res <- system(paste0(file.path(neo.path, "bin", "neo4j-shell"), " -c ","'",cur.query, "'"), intern=T)
                
                return(clean.neo4j.res(cur.res))
           })
    
    lab.dta <- do.call("rbind", lab.list)
    
    #make an igraph object
    
    return(graph.data.frame(lab.dta))
}

#neo.path=normalizePath("neo4j-community-2.1.6/")
make.graph.struct <- function(neo.graph, graph.struct.path="test_graph_struct.json")
{
    
    if (missing(graph.struct.path) || is.null(graph.struct.path) || all(is.na(graph.struct.path)) || (is.character(graph.struct.path) == F))
    {
       stop("ERROR: Need to supply a valid path for the graph_struct file")
    }
    
    if (missing(neo.graph) || is.null(neo.graph) || all(is.na(neo.graph)) || (class(neo.graph) != "igraph"))
    {
        stop("ERROR: neo.graph needs to be an igraph object")
    }
    
    lab.graph <- neo.graph
    
    un.neo.graph <- as.undirected(neo.graph, edge.attr.comb="first")

    graph.list <- lapply(get.adjedgelist(un.neo.graph), function(x)
                         {
                            unique.x <- unique(x)
                            
                            temp.list <- lapply(unique.x, function(y)
                                   {
                                        return(V(un.neo.graph)[inc(E(un.neo.graph)[y])]$name)
                                   })
                            
                            names(temp.list) <- E(un.neo.graph)[unique.x]$type
                            return(temp.list)
                         })
    
    pretty.graph.list <- mapply(function(x,y){
        
        lapply(x, function(z)
               {
                    if (length(z) > 1)
                    {
                        return(z[z != y])
                    }else{
                        return(z)
                    }
                    
               })
    }, graph.list, names(graph.list))
    
    write(toJSON(pretty.graph.list), file=graph.struct.path)
    
}

read.graph_struct  <- function(graph_struct="/var/www/hitwalker_2_inst/graph_struct.json")
{
    require(rjson)
    require(igraph)
    
    json.obj <- fromJSON(file=graph_struct)
    
    #warnings are about the rownames...
    json.dta <- suppressWarnings(data.frame(do.call(rbind, lapply(names(json.obj), function(x) {
        if (length(json.obj[[x]]) > 0)
        {
            temp.dta <- cbind(x, names(json.obj[[x]]), unlist(json.obj[[x]]))
        }else{
            return(NULL)
        }
        
    })), stringsAsFactors=F))
    names(json.dta) <- c("from", "label", "to")
    json.dta <- json.dta[!duplicated(json.dta$label),]
    json.graph <- graph.data.frame(json.dta[,c("from", "to", "label")], directed=F)
    return(json.graph)
    #plot(json.graph)
}

vcf.to.graph <- function()
{
    require(VariantAnnotation)
    requireNamespace(ensemblVEP)
    
    vcf.obj <- readVcf("test_out.vcf.new", genome="ncbi")
    
    vcf.csq <- ensemblVEP::parseCSQToGRanges(vcf.obj, VCFRowID=rownames(vcf.obj))
    
    annot.dta <- make.additional.flags(vcf.csq)
    
    #this will form the basic of Variation -> Transcript relationships (IMPACTS)

    #also needs to have a LabID -> Variation mapping (DNA_DIFF)
}


make.additional.flags <- function(annot.grange)
{
  
    colnames(mcols(annot.grange))[colnames(mcols(annot.grange)) == "HGNC"] <- "Symbol"
    mcols(annot.grange)$coding <- as.integer(ifelse(is.na(mcols(annot.grange)$CDS_position), FALSE, TRUE))
    synonymous <- c("synonymous_variant", "stop_retained_variant")
    non.synonymous <- c("stop_gained", "frameshift_variant", "stop_lost", "initiator_codon_variant", "inframe_insertion", "inframe_deletion", "missense_variant")
    syn.var <- in.csv.col(vec = mcols(annot.grange)$Consequence, search.vals = synonymous, delim.str = "&", match.func = any)
    non.syn.var <- in.csv.col(vec = mcols(annot.grange)$Consequence, search.vals = non.synonymous, delim.str = "&", match.func = any)
    both <- syn.var & non.syn.var
    syn.var[both] <- FALSE
    
    mcols(annot.grange)$Consequence <- gsub("&", ";", as.character(mcols(annot.grange)$Consequence))
    
    mcols(annot.grange)$Cons_cat <- "Other"
    mcols(annot.grange)$Cons_cat[syn.var] <- "Synonymous"
    mcols(annot.grange)$Cons_cat[non.syn.var] <- "NonSynonymous"
    
    mcols(annot.grange)$var_type <- mapply(function(ref.width, alt.width) {
        
        if (alt.width == ref.width && ref.width == 1) {
            return("SNV")
        }
        else if (alt.width == ref.width && ref.width > 1) {
            return("MNV")
        }
        else if (alt.width > ref.width) {
            return("INS")
        }
        else if (alt.width < ref.width) {
            return("DEL")
        }
        else {
            stop("ERROR: Unknown variant type")
        }
    }, width(mcols(annot.grange)$REF), width(mcols(annot.grange)$ALT))
    
    if('GMAF' %in% names(mcols(annot.grange)))
    {
        in.1kg <- is.na(mcols(annot.grange)$GMAF) == FALSE
        mcols(annot.grange)$in_1kg <- as.integer(in.1kg)
    }
    
    if ('Existing_variation' %in% names(mcols(annot.grange)) && all(is.na(mcols(annot.grange)$Existing_variation)) == F)
    {
        in.dbsnp <- grepl("rs\\d+", mcols(annot.grange)$Existing_variation)
        mcols(annot.grange)$in_dbsnp <- as.integer(in.dbsnp)
    }
    
    if ('CANONICAL' %in% names(mcols(annot.grange)))
    {
        stopifnot(levels(mcols(annot.grange)$CANONICAL) == "YES")
        mcols(annot.grange)$CANONICAL <- as.character(mcols(annot.grange)$CANONICAL)
        mcols(annot.grange)$CANONICAL[is.na(mcols(annot.grange)$CANONICAL)] <- "NO"
    }
   
    return(annot.grange)
}

in.csv.col <- function(vec, search.vals, delim.str=",", match.func=any)
{
    vec <- as.character(vec)
    
    in.vec <- sapply(strsplit(vec, delim.str), function(x)
                     {
                        match.func(x %in% search.vals) 
                     })
    return(in.vec)
}

#taken from older code and needs to be adapted and checked...
create.initial.graph <- function()
{
    library(igraph)

    string.graph <- read.delim("/Users/bottomly/Desktop/tyner_results/string_preparation/protein.links.detailed.v9.05.9606.ncol", sep="", header=FALSE, stringsAsFactors=FALSE)
    names(string.graph) <- c("from.node", "to.node", "weight")
    
    string.graph$weight <- string.graph$weight/1000
    
    use.graph <- graph.data.frame(string.graph, directed=FALSE)
    
    temp.sparse <- get.adjacency(graph=use.graph, sparse=TRUE, attr="weight", type="both")/2
    
    writeMM(temp.sparse, file="protein.links.detailed.v9.05.9606.mm.mtx")
    
    all(colnames(temp.sparse) == rownames(temp.sparse))#TRUE
    
    write.table(rownames(temp.sparse), file="protein.links.detailed.v9.05.9606.mm.names", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
}

