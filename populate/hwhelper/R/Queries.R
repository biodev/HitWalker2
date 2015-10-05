#' GATK/Ensembl VEP VCF Representation
#' 
#' A class for representing general variant data stored in GATK flavored VCF files annotated with Ensembl VEP
#'

gatkvcf_class <- setClass(Class="VCFTable", representation=list(ensembl.to.entrez="data.frame"), contains="DenseNeoData")

VCFTable <- function(vcf.dta, ensembl.to.entrez, node.name="variation", sample.edge.name="HAS_DNASEQ"){
  
  #add in some ancillary data:
  
  #cohort count/frequency in terms of samples
  coh.count <- aggregate(samples~HGVSc, length, data=vcf.dta)
  coh.count$cohort_freq <- coh.count$samples/length(unique(vcf.dta$samples))
  names(coh.count)[names(coh.count) == "samples"] <- "cohort_count"
  
  coh.count$HGVSc <- as.character(coh.count$HGVSc)
  vcf.dta$HGVSc <- as.character( vcf.dta$HGVSc)
  
  vcf.dta <- merge(vcf.dta, coh.count, by="HGVSc", all.x=T, all.y=F, incomparables=NA, sort=F)
  
  #some slicing and dicing of the Existing_variation info
  
  vcf.dta$in_dbsnp <- as.integer(grepl("rs\\d+", vcf.dta$Existing_variation))
  vcf.dta$in_esp <- as.integer(grepl("TMP_ESP_", vcf.dta$Existing_variation))
  vcf.dta$in_cosmic <- as.integer(grepl("COSM", vcf.dta$Existing_variation))
  
  #add in the entrez gene mappings (note that there will be instances of one ensembl to many entrez and vice versa that is not dealt with here...)
  
  actual.maps <- ensembl.to.entrez[complete.cases(ensembl.to.entrez),] 
  
  vcf.dta$Gene <- as.character(vcf.dta$Gene)
  
  samp.vcf <- merge(vcf.dta, actual.maps, by="Gene", all=F, incomparables=NA, sort=F)
  
  names(samp.vcf)[names(samp.vcf)=="Gene"] <- "Ensembl_Gene"
  names(samp.vcf)[names(samp.vcf)=="entrezID"] <- "gene"
  names(samp.vcf)[names(samp.vcf) == "HGVSc"] <- "name"
  names(samp.vcf)[names(samp.vcf) == "samples"] <- "sample"
  
  samp.vcf <- samp.vcf[,c("sample", "gene", "name", "FILTER", "MQ0", "seqnames", "start", "end", "REF", "ALT", "Existing_variation", "Consequence",
                         "HGVSp", "GMAF", "RefSeq", "SIFT", "PolyPhen", "cohort_count", "cohort_freq", "in_dbsnp", "in_esp")]
  
  samp.vcf <- factors.to.chars(samp.vcf)
  
  return(new("VCFTable", data=samp.vcf, ensembl.to.entrez=ensembl.to.entrez, node.name=node.name, sample.edge.name=sample.edge.name))
  
}


#' @rdname class_helpers
#' @param vcf.file The path to a VCF file that has been annotated by Ensembl VEP.  If the 'PICK' column is defined, only those variants will be imported.
#' @param genome Optional genome version used in the VCF file.
#' @param ensembl.to.entrez needs to be a data.frame with a 'Gene' column indicating the Ensembl IDs and an entrezID column containing the EntrezIDs or NA if not available.
parse.vcf.vt <- function(vcf.file, genome="Unknown",ensembl.to.entrez=get.biomart.mapping(), node.name="variation", sample.edge.name="HAS_DNASEQ"){
  
  require(VariantAnnotation)
  require(genotools)
  
  vcf.obj <- readVcf(file=vcf.file, genome=genome)
  
  vcf.dta <- vcf.to.table(vcf.obj)
  
  return(VCFTable(vcf.dta, ensembl.to.entrez))
  
}

#' @rdname class_helpers
#' @param file.name The path to the CCLE MAF file.
readMAF.ccle <- function(file.name, node.name="variation", sample.edge.name="HAS_DNASEQ")
{
  use.maf <- read.delim(file.name, sep="\t", stringsAsFactors=F)
  
  keep.maf <- use.maf[,c("Entrez_Gene_Id", "Genome_Change", "Variant_Classification", "Annotation_Transcript", "Transcript_Strand", "cDNA_Change", "Codon_Change", "Protein_Change",
                         "Tumor_Sample_Barcode", "Center", "Sequencer", "Alternative_allele_reads_count", "Reference_allele_reads_count", "dbSNP_RS", "dbSNP_Val_Status")]
  
  names(keep.maf)[c(1:2,9)] <- c("gene", "name", "sample")
  
  return(new("DenseNeoData", data=keep.maf, node.name=node.name, sample.edge.name=sample.edge.name))
}

#older stuff to be incorporated into docs etc

# #' @describeIn CCLEMaf_class Loads the sample to variant data into Neo4j.  The sample annotation is taken from the 'Tumor_Sample_Barcode' column and the variant annotation is taken from the 'Genome_Change' column.
# #' The following additional columns are added as attributes as well:
# #' \itemize{
# #'   \item Center
# #'   \item Sequencer
# #'   \item Alternative_allele_reads_count
# #'   \item Reference_allele_reads_count
# #'   \item dbSNP_RS
# #'   \item dbSNP_Val_Status
# #' }
# #' @param obj The optional path to a Neo4j database.
# #' @param neo.path The optional path to a Neo4j database.
# #' @describeIn CCLEMaf_class Loads the variant to gene data into Neo4j.  The mapping is produced using the 'Genome_Change' and 'Entrez_Gene_Id' columns.
# #' Additionally, the following columns are also included as attributes:
# #' \itemize{
# #'      \item Variant_Classification
# #'      \item Annotation_Transcript
# #'      \item Transcript_Strand
# #'      \item cDNA_Change
# #'      \item Codon_Change
# #'      \item Protein_Change
# #' }
# #' @param gene.model The type of gene model to utilize.  Currently only 'entrez' is supported.


# setMethod("sampleNames", signature("VCFTable"), function(object){
#   return(unique(object@vcf.dta$samples))
# })
# 
# setMethod("fromSample", signature("VCFTable"), function(obj, neo.path=NULL){
#   #first sample -> variant
#   #the name here will be derived from the HGVSp column as that provides potentially enough information to uniquely id a variant as these should all be protein-coding variants
#   #will keep missing values as "" for now
#   
#   sample.vcf <- obj@vcf.dta
#   
#   samp.vcf <- sample.vcf[,c("samples", "HGVSc", "FILTER", "MQ0", "seqnames", "start", "end", "REF", "ALT", "Existing_variation")]
#   names(samp.vcf) <- c("sample", nodeName(obj), "FILTER", "MQ0", paste(nodeName(obj), c("seqnames", "start", "end", "REF", "ALT", "Existing_variation"), sep="."))
#   
#   samp.vcf <- factors.to.chars(samp.vcf)
#   
#   samp.vcf[is.na(samp.vcf)] <- ""
#   
#   load.neo4j(.data=samp.vcf, edge.name=sampleEdge(obj), commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")
# })
# 
# setMethod("toGene", signature("VCFTable"), function(obj, neo.path=NULL, gene.model="entrez"){
#   #then add in the Variation->EntrezID relationships
#   
#   if (gene.model != "entrez")
#   {
#     stop("ERROR: Only gene.model = 'entrez' is currently supported")
#   }
#   
#   cur.vcf <- obj@vcf.dta
#   
#   #only keep one row for each gene/variant then convert to entrez IDs
#   #note there are still potentially many entrez ids to one ensembl ID..
#   
#   ens.gene.dta <- cur.vcf[!duplicated(cur.vcf[,c("Gene", "HGVSp")]),]
#   
#   ens.gene.dta$Gene <- as.character(ens.gene.dta$Gene)
#   obj@ensembl.to.entrez$Gene <- as.character(obj@ensembl.to.entrez$Gene)
#   
#   ens.gene.merge <- merge(ens.gene.dta, obj@ensembl.to.entrez, by="Gene", all=F, incomparables=NA, sort=F)
#   
#   ent.genes <- ens.gene.merge[is.na(ens.gene.merge$entrezID)==F,]
#   
#   ent.genes <- ent.genes[,c("HGVSc", "entrezID", "Consequence", "Protein_position", "Amino_acids", "SIFT", "PolyPhen")]
#   
#   names(ent.genes) <- c(nodeName(obj), "entrezID","consequence", "protein_position", "amino_acids", "sift", "polyphen")
#   
#   ent.genes <- factors.to.chars(ent.genes)
#   
#   ent.genes[is.na(ent.genes)] <- ""
#   
#   load.neo4j(.data=ent.genes, edge.name=geneEdge(obj), commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&", unique.rels=F)  
#   
# })


