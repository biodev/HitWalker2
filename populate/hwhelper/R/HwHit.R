#' Hit Parameter Data for HitWalker2
#'
#' Basic class representing common components for hit parameters that all classes which load data into Neo4j for HitWalker2 should be derived from.
#'
#' @slot default The default value used to determine if a datatype should be considered a 'hit' for a given sample
#' @slot direction Direction of the test
#' @slot range Valid range (vector) that the user can adjust the values towards
#' @slot display_name Name to be displayed to the user
HwHit <- setClass(Class="HwHit", representation=list(default="numeric", direction="character", range="numeric", display_name="character"))

AnnotatedMatrix <- setClass(Class="AnnotatedMatrix", representation = list(matrix="matrix", mapping="data.frame"), contains="NeoData")

setMethod("getMatrix", signature("AnnotatedMatrix"), function(obj)
{
  return(obj@matrix)
})

setMethod("getAnnotation", signature("AnnotatedMatrix"), function(obj){
  return(obj@mapping)
})

setMethod("sampleNames", signature("AnnotatedMatrix"), function(object){
  return(unique(colnames(getMatrix(object))))
})

setMethod("toGene", signature("AnnotatedMatrix"), function(obj, neo.path=NULL, gene.model=c("entrez", "ensembl")){
  
  mat.genes <- getAnnotation(obj)
  
  mat.genes <- mat.genes[complete.cases(mat.genes),]
  
  if ("weight" %in% colnames(mat.genes) == F){
    mat.genes$weight <- 1
  }
  
  mat.genes <- mat.genes[,c(nodeName(obj), "gene", "weight")]
  names(mat.genes) <- c(nodeName(obj), switch(gene.model, entrez="entrezID", ensembl="ensembl"), "weight")
  
  load.neo4j(.data=mat.genes, edge.name=geneEdge(obj), commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")
})

setMethod("fromSample", signature("AnnotatedMatrix"), function(obj, neo.path=NULL){
  
  gene.mat <- getMatrix(obj)
  
  gene.dta <- melt(gene.mat)
  
  names(gene.dta) <- c(nodeName(obj), "sample", "score")
  
  for(i in 1:2){
    gene.dta[,i] <- as.character(gene.dta[,i])
  }
  
  #remove any na scores
  
  if (sum(is.na(gene.dta$score)) != sum(is.na(gene.dta))){
    stop("ERROR: NAs in other columns than the score column detected")
  }
  
  gene.dta <- gene.dta[complete.cases(gene.dta),]
  
  #reorder the edges
  
  gene.dta <- gene.dta[,c("sample", nodeName(obj), "score")]
  
  load.neo4j(.data=gene.dta, edge.name=sampleEdge(obj), commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")
})

##in progress, generalizing the templates so it is less work for the configurer (person)

fill.query.slots <- function(obj){
  
  sum.dir <- ifelse(obj@direction == ">", "MAX", "MIN")
  
  hit.crit <- ifelse(obj@is.dense, "true", paste0(sum.dir, "(r.score*r2.weight) > $PAR_NAME$"))
  
  type.name <- typeName(obj)
  
  if (length(obj@base.query) == 0){
    
    obj@base.query <- gsub("\\s*\n\\s*", " ", paste0('MATCH (subject:$SUBJECT$)-[d:DERIVED]-(sample)-[r:',sampleEdge(obj),']-()-[r2:',geneEdge(obj),']-(o:EntrezID{name:{GENE}}) WHERE d.type = "',type.name,'" AND 
     ANY(x IN [subject.name, sample.name] WHERE x IN {SAMPLE})
     RETURN o.name AS gene, subject.name AS sample, "$DATA_NAME$" AS var, ',sum.dir,'(r.score*r2.weight) AS score, ',hit.crit,' AS is_hit;'))
    
  }
  
  if (length(obj@template.query) == 0){
    
    obj@template.query <- gsub("\\s*\n\\s*", " ", paste0('MATCH(subject:$SUBJECT$)-[d:DERIVED]-()-[r:',sampleEdge(obj),']-()-[r2:',geneEdge(obj),']-(gene:EntrezID) WHERE d.type = "',type.name,'" 
    AND r.score*r2.weight > $PAR_NAME$ AND $$lower_coll_type$$.name IN {$$coll_type$$}
    WITH $$lower_ret_type$$.name AS ret_type, COUNT(DISTINCT $$lower_coll_type$$.name) AS use_coll ORDER BY use_coll DESC RETURN ret_type, use_coll'))
  
    }
  
  return(obj)
}

#tests for geneset and HW2exprSet : temp <- HW2exprSet(ExpressionSet())
#base.query='MATCH (subject:$SUBJECT$)-[d:DERIVED]-(sample)-[r:HAS_GENE_SET]-(m)-[r2:ASSIGNED_TO]-(o:EntrezID{name:{GENE}}) WHERE d.type = "Gene_Set" AND subject.name IN {SAMPLE}
#RETURN o.name AS gene, subject.name AS sample, "$DATA_NAME$" AS var, r.score*r2.weight AS score, r.score*r2.weight > $PAR_NAME$ AS is_hit;',
#template.query='MATCH (subject:$SUBJECT$)-[d:DERIVED]-()-[r:HAS_GENE_SET]-(m)-[r2:ASSIGNED_TO]-(gene:EntrezID) WHERE d.type = "Gene_Set" AND r.score*r2.weight > $PAR_NAME$ AND
#$$lower_coll_type$$.name IN {$$coll_type$$} WITH $$lower_ret_type$$.name AS ret_type, COUNT(DISTINCT $$lower_coll_type$$.name) AS use_coll
#ORDER BY use_coll DESC RETURN ret_type, use_coll',
#

##Note here, that for dense datatypes like expression, don't actually return possible hits
#' Basic Representation for Expression Array Data
#'
#' A basic class for representing Affymetrix expression array data.
#'
#' @slot exprs, an \code{ExpressionSet} containing expression data.
HW2exprSet_class <- setClass(Class="HW2exprSet", representation=list(exprs="ExpressionSet"), contains=c("NeoData", "HwHit"),
                             prototype=list(sample.edge.name="HAS_EXPRESSION", gene.edge.name="PS_MAPPED_TO", node.name="probeSet",
                                            default=2, direction=">", range=c(2, 10), display_name="Expression Z-Score Hit Threshold", is.dense=TRUE))

#' @rdname class_helpers
#' @param exprs An ExpressionSet
#' @param sample.edge.name The name of the relationship between the samples and the assay identifier e.g. probeset, drug or variant.
#' @param gene.edge.name The name of the relationship between the assay identifiers and genes
#' @param node.name The name that should be given to the assay units in the Neo4j database.
HW2exprSet <- function(exprs, sample.edge.name="HAS_EXPRESSION", gene.edge.name="PS_MAPPED_TO", node.name="probeSet"){
  
  return(fill.query.slots(new("HW2exprSet", exprs=exprs, sample.edge.name=sample.edge.name, gene.edge.name=gene.edge.name, node.name=node.name)))
}

#' @describeIn HW2exprSet_class Implements the loading of sample to probeset data into Neo4j.  The sample names are derived from the column names and the
#' probeset names are derived from the rownames. The score attribute is derived from the values of the matrix.  The expression values are converted to z scores
#' prior to loading and only those sample to gene relationships that exceed the specified lower range are kept.
#' @param obj An object of class \code{HW2exprSet}.
#' @param neo.path The optional path to a Neo4j database.
setMethod("fromSample", signature("HW2exprSet"), function(obj, neo.path=NULL){
  
  use.exprs <- exprs(obj@exprs)
  
  #by default we z transform the data
  
  scale.exprs <- t(scale(t(use.exprs)))
  
  melt.use.exprs <- melt(scale.exprs)
  names(melt.use.exprs) <- c("probeset", "sample", "score")
  
  melt.use.exprs$sample <- as.character(melt.use.exprs$sample)
  melt.use.exprs$probeset <- as.character(melt.use.exprs$probeset)
  
  use.expr <- paste("melt.use.exprs$score", obj@direction, min(obj@range))
  
  exprs.hits <- eval(parse(text=use.expr))
  
  melt.use.exprs <- melt.use.exprs[exprs.hits,]
  
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
  
  if (length(annotation(obj@exprs)) > 0){
    
    if (grepl("\\.db$", annotation(obj@exprs), perl=T)){
      
      annotation.package <- annotation(obj@exprs)
      
    }else{
      
      annotation.package <- paste0(annotation(obj@exprs), ".db")
    }
  }else{
    stop("ERROR: Need to specify an annotation package as part of the ExpressionSet")
  }
  
  require(annotation.package, character.only=T)
  
  #Then create a mapping file from probeset to gene
  
  gene.type <- switch(gene.model, entrez="ENTREZID", ensembl="ENSEMBL")
  
  probe.to.gene <- suppressWarnings(select(eval(parse(text=annotation.package)), keys=featureNames(obj@exprs), column=gene.type, keytypes="PROBEID"))
  
  #Discard for now those that do not map to either type of genes.
  
  probe.to.gene <- probe.to.gene[is.na(probe.to.gene[,gene.type]) == F,]
  
  #probesets that map to multiple genes are perhaps not that reliable either, so discard them as well for now...
  
  dup.probes <- probe.to.gene$PROBEID[duplicated(probe.to.gene$PROBEID)]
  
  probe.to.gene <- probe.to.gene[probe.to.gene$PROBEID %in% dup.probes == F,]
  names(probe.to.gene) <- c(nodeName(obj), switch(gene.model, entrez="entrezID", ensembl="ensembl"))
  
  probe.to.gene$weight <- 1
  
  #load the probe->gene mappings
  
  load.neo4j(.data=probe.to.gene, edge.name=geneEdge(obj), commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")
})

GeneSet_class <- setClass(Class="GeneSet", contains=c("AnnotatedMatrix", "HwHit"),
                          prototype = list(
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
  
  return(fill.query.slots(new("GeneSet", matrix=mat, mapping=mapping)))
}


#GeneScore and siRNA representation here

#' Drug Data Representation
#'
#' A basic class for representing drug panel results for a set of samples
#'
#' @aliases getMatrix getAnnotation
#' @slot matrix, a matrix containing a sinlge summary value for each drug (rows) and each sample (columns).  
#' @slot mapping A \code{data.frame} containing the mapping from drug to gene.  It should contain a 'drug' column, a 'gene' column and a 'weight' column which indicates the confidence of the mapping.
DrugMatrix_class <- setClass(Class="DrugMatrix", contains=c("HwHit", "AnnotatedMatrix"),
                             prototype=list(base.query='MATCH (subject:$SUBJECT$)-[d:DERIVED]-(sample)-[r:HAS_DRUG_ASSAY]-(m)-[r2:ACTS_ON]-(o:EntrezID{name:{GENE}}) WHERE d.type = "Drug_Assay" AND 
                                            ANY(x IN [subject.name, sample.name] WHERE x IN {SAMPLE})
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