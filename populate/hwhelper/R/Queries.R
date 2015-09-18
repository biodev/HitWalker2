#' GATK/Ensembl VEP VCF Representation
#' 
#' A class for representing general variant data stored in GATK flavored VCF files annotated with Ensembl VEP
#'

gatkvcf_class <- setClass(Class="VCFTable", representation=list(vcf.dta="data.frame", ensembl.to.entrez="data.frame"), contains="NeoData",
                          prototype = list(base.query='MATCH (subject:$SUBJECT$)-[d:DERIVED]-(sample)-[r:HAS_DNASEQ]-(var)-[r2:IMPACTS]-(gene:EntrezID{name:{GENE}})-[:REFERRED_TO]-(symb) WHERE d.type = "DNASeq" AND subject.name IN {SAMPLE}
                                          RETURN var.name AS Variant_Position, gene.name AS Gene, symb.name AS Symbol,
                                           r2.consequnce AS Consequence , REPLACE(RTRIM(REDUCE(str="",n IN var.Existing_Variation|str+n+" ")), " ", ";") AS Existing_Variation,
                                           r2.sift AS SIFT, r2.polyphen AS Polyphen, 0 AS query_ind, 1 AS gene_ind, var.name + "_" + gene.name AS row_id, subject.name AS Sample ',
                                           template.query='MATCH (subject:$SUBJECT$)-[d:DERIVED]-()-[r:HAS_DNASEQ]-(var)-[r2:IMPACTS]-(gene:EntrezID) WHERE d.type = "DNASeq" AND $$lower_coll_type$$.name IN {$$coll_type$$}
                                          WITH $$lower_ret_type$$.name AS ret_type, COUNT(DISTINCT $$lower_coll_type$$.name) AS use_coll ORDER BY use_coll DESC RETURN ret_type, use_coll',
                                           sample.edge.name="HAS_DNASEQ", 
                                           gene.edge.name="IMPACTS", 
                                           node.name="variation"))

#' @rdname class_helpers
#' @param vcf.file The path to a VCF file that has been annotated by Ensembl VEP.  If the 'PICK' column is defined, only those variants will be imported.
#' @param genome Optional genome version used in the VCF file.
#' @param ensembl.to.entrez needs to be a data.frame with a 'Gene' column indicating the Ensembl IDs and an entrezID column containing the EntrezIDs or NA if not available.
parse.vcf.vt <- function(vcf.file, genome="Unknown",ensembl.to.entrez=get.biomart.mapping()){
  
  require(VariantAnnotation)
  require(genotools)
  
  vcf.obj <- readVcf(file=vcf.file, genome=genome)
  
  vcf.dta <- vcf.to.table(vcf.obj)
  
  sub.vcf.dta <- vcf.dta[,c("seqnames", "start", "end", "REF", "ALT", "FILTER", "MQ0", "Gene", "Feature", "HGVSp", "HGVSc",
                            "Consequence", "Protein_position", "Amino_acids", "Existing_variation", "SYMBOL", "SIFT", 
                            "PolyPhen", "RefSeq", "samples", "genotypes")]
  
  return(new("VCFTable", vcf.dta=sub.vcf.dta, ensembl.to.entrez=ensembl.to.entrez, sample.edge.name=sample.edge.name, gene.edge.name=gene.edge.name, node.name=node.name))
}

setMethod("sampleNames", signature("VCFTable"), function(object){
  return(unique(object@vcf.dta$samples))
})

setMethod("fromSample", signature("VCFTable"), function(obj, neo.path=NULL){
  #first sample -> variant
  #the name here will be derived from the HGVSp column as that provides potentially enough information to uniquely id a variant as these should all be protein-coding variants
  #will keep missing values as "" for now
  
  sample.vcf <- obj@vcf.dta
  
  samp.vcf <- sample.vcf[,c("samples", "HGVSc", "FILTER", "MQ0", "seqnames", "start", "end", "REF", "ALT", "Existing_variation")]
  names(samp.vcf) <- c("sample", nodeName(obj), "FILTER", "MQ0", paste(nodeName(obj), c("seqnames", "start", "end", "REF", "ALT", "Existing_variation"), sep="."))
  
  samp.vcf <- factors.to.chars(samp.vcf)
  
  samp.vcf[is.na(samp.vcf)] <- ""
  
  load.neo4j(.data=samp.vcf, edge.name=sampleEdge(obj), commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")
})

setMethod("toGene", signature("VCFTable"), function(obj, neo.path=NULL, gene.model="entrez"){
  #then add in the Variation->EntrezID relationships
  
  if (gene.model != "entrez")
  {
    stop("ERROR: Only gene.model = 'entrez' is currently supported")
  }
  
  cur.vcf <- obj@vcf.dta
  
  #only keep one row for each gene/variant then convert to entrez IDs
  #note there are still potentially many entrez ids to one ensembl ID..
  
  ens.gene.dta <- cur.vcf[!duplicated(cur.vcf[,c("Gene", "HGVSp")]),]
  
  ens.gene.dta$Gene <- as.character(ens.gene.dta$Gene)
  obj@ensembl.to.entrez$Gene <- as.character(obj@ensembl.to.entrez$Gene)
  
  ens.gene.merge <- merge(ens.gene.dta, obj@ensembl.to.entrez, by="Gene", all=F, incomparables=NA, sort=F)
  
  ent.genes <- ens.gene.merge[is.na(ens.gene.merge$entrezID)==F,]
  
  ent.genes <- ent.genes[,c("HGVSc", "entrezID", "Consequence", "Protein_position", "Amino_acids", "SIFT", "PolyPhen")]
  
  names(ent.genes) <- c(nodeName(obj), "entrezID","consequence", "protein_position", "amino_acids", "sift", "polyphen")
  
  ent.genes <- factors.to.chars(ent.genes)
  
  ent.genes[is.na(ent.genes)] <- ""
  
  load.neo4j(.data=ent.genes, edge.name=geneEdge(obj), commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&", unique.rels=F)  
  
})



#' CCLE MAF Representation
#'
#' A basic class for representing the Cancer Cell Line Encyclopedia style Mutation Annotation Format files.
#'
#' @slot maf, a \code{data.frame} representation of the MAF file.
CCLEMaf_class <- setClass(Class="CCLEMaf", representation=list(maf="data.frame"), contains="NeoData",
                          prototype=list(
                            base.query='MATCH (subject:$SUBJECT$)-[d:DERIVED]-(sample)-[r:HAS_DNASEQ]-(var)-[r2:IMPACTS]-(gene:EntrezID{name:{GENE}})-[:REFERRED_TO]-(symb) WHERE d.type = "DNASeq" AND subject.name IN {SAMPLE}
                            RETURN var.name AS Variant_Position, r2.transcript AS Transcript, gene.name AS Gene, symb.name AS Symbol,
                            r.ref_counts as Ref_Counts, r.alt_counts AS Alt_Counts, REPLACE(RTRIM(REDUCE(str="",n IN var.dbsnp|str+n+" ")), " ", ";") AS dbSNP,
                            r2.variant_classification AS Variant_classification, r2.protein AS Protein_Change, 0 AS query_ind, 2 AS gene_ind, var.name + "_" + gene.name AS row_id, subject.name AS Sample ',
                            
                            template.query='MATCH (subject:$SUBJECT$)-[d:DERIVED]-()-[r:HAS_DNASEQ]-(var)-[r2:IMPACTS]-(gene:EntrezID) WHERE d.type = "DNASeq" AND $$lower_coll_type$$.name IN {$$coll_type$$}
                            WITH $$lower_ret_type$$.name AS ret_type, COUNT(DISTINCT $$lower_coll_type$$.name) AS use_coll ORDER BY use_coll DESC RETURN ret_type, use_coll'
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

