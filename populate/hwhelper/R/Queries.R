#' GATK/Ensembl VEP VCF Representation
#' 
#' A class for representing general variant data stored in GATK flavored VCF files annotated with Ensembl VEP
#'

VCFTable <- function(vcf.dta,node.name="variation", sample.edge.name="HAS_DNASEQ"){
  
  #add in some ancillary data:
  
  #cohort count/frequency in terms of samples
  coh.count <- aggregate(sample~HGVSc, length, data=vcf.dta)
  coh.count$cohort_freq <- coh.count$sample/length(unique(vcf.dta$sample))
  names(coh.count)[names(coh.count) == "sample"] <- "cohort_count"
  
  coh.count$HGVSc <- as.character(coh.count$HGVSc)
  vcf.dta$HGVSc <- as.character( vcf.dta$HGVSc)
  
  vcf.dta <- merge(vcf.dta, coh.count, by="HGVSc", all.x=T, all.y=F, incomparables=NA, sort=F)
  
  #some slicing and dicing of the Existing_variation info
  
  vcf.dta$in_dbsnp <- as.integer(grepl("rs\\d+", vcf.dta$Existing_variation))
  vcf.dta$in_esp <- as.integer(grepl("TMP_ESP_", vcf.dta$Existing_variation))
  vcf.dta$in_cosmic <- as.integer(grepl("COSM", vcf.dta$Existing_variation))
  
  #add in the entrez gene mappings (note that there will be instances of one ensembl to many entrez and vice versa that is not dealt with here...)
  
  names(vcf.dta)[names(vcf.dta)=="Gene"] <- "gene"
  names(vcf.dta)[names(vcf.dta) == "HGVSc"] <- "name"
  
  vcf.dta <- vcf.dta[,c("sample", "gene", "name", "seqnames", "start", "end", "REF", "ALT", "Protein_position", "Amino_acids", "allele_count", "allele_reads", 
                        "total_reads", "MQ0", "FS", "MQ", "QD", "SB" ,"Existing_variation", "Consequence","HGVSp", "GMAF", "Feature", "SIFT", "PolyPhen", 
                        "cohort_count", "cohort_freq", "in_dbsnp", "in_esp", "in_cosmic")]
  
  vcf.dta <- factors.to.chars(vcf.dta)
  
  return(new("DenseNeoData", data=vcf.dta, node.name=node.name, sample.edge.name=sample.edge.name))
  
}

.csq.to.df <- function(gds, header){
  
  csq.dta.list <- seqGetData(gds, "annotation/info/CSQ")
  
  exp.inds <- unlist(lapply(seq_along(csq.dta.list$length), function(x) rep(x, csq.dta.list$length[x])))
  
  raw <- strsplit(csq.dta.list$data, "\\|")
  
  csq <- matrix(nrow = length(raw), ncol = length(header))
  
  for (i in 1:nrow(csq)) csq[i, 1:length(raw[[i]])] <- raw[[i]]
  colnames(csq) <- header
  
  var.ind <- index.gdsn(gds, "variant.id", silent = TRUE)
  variant.id <- read.gdsn(var.ind,simplify="none")
  
  csq <- data.frame(variant_id=variant.id[exp.inds],csq, stringsAsFactors = F)
  
  csq <- csq[order(csq$variant_id, decreasing=F),]
  
  return(csq)
}

#A slower modification of the same function from SeqVarTools that can deal with multi-alleleic vars
.getVariableLengthData <- function(gdsobj, var.name, use.names = TRUE){
  var.list <- seqApply(gdsobj, var.name, function(x) {
    x
  })
  
  cols <- max(sapply(var.list, ncol))
  rows <- unique(sapply(var.list, nrow))
  stopifnot(length(rows) == 1)
  
  ar.shape <- c(rows, cols, length(var.list))
  
  new.list <- lapply(var.list, function(x){
    cbind(x, matrix(NA, nrow=rows, ncol=cols-ncol(x)))
  })
  
  var <- array(unlist(new.list, use.names = FALSE), dim=ar.shape)
  
  var <- aperm(var, c(2, 1, 3))
  if (dim(var)[1] == 1) {
    var <- var[1, , ]
  }
  if (length(dim(var)) == 3) {
    dimnames(var) <- list(n = NULL, sample = NULL, variant = NULL)
  }
  else {
    dimnames(var) <- list(sample = NULL, variant = NULL)
  }
  if (use.names) 
    SeqVarTools:::.applyNames(gdsobj, var)
  else var
  
}

.info <- function(gds, rm.cols=c("CSQ", "@CSQ"), info.import=NULL){
  
  n <- index.gdsn(gds, "annotation/info", silent = TRUE)
  
  keep.cols <- setdiff(ls.gdsn(n), rm.cols)
  
  res.mat <- sapply(keep.cols, function(x){
    col.ind <- index.gdsn(n, x, silent = TRUE)
    read.gdsn(col.ind)
  })
  
  if(is.null(info.import) == F){
    
    for(i in setdiff(info.import[info.import %in% rm.cols == F], colnames(res.mat))){
      res.mat <- cbind(res.mat, NA)
      colnames(res.mat)[ncol(res.mat)] <- i
    }
    
  }
  
  return(res.mat)
}


#' @rdname class_helpers
#' @param vcfs The path to VCF file(s) that have been annotated by Ensembl VEP using the '--refseq' option.  
#' Our typical workflow additionally involves subsetting to protein impacting variants 
#' and keeping only consequences chosen via the allele-specific 'PICK' column.
#' @param keep.gds Should the intermediate GDS file be kept after the \code{data.frame} is generated?
make.vcf.table <- function(vcfs, info.import=c("FS", "MQ0", "MQ", "QD", "SB", "CSQ"), fmt.import="AD", gds.out="variant.gds", keep.gds=F){
  
  if(file.exists(gds.out)){
    
    message(paste("File:", gds.out, "exists, utilizing it for this analysis"))
    
    gds <- seqOpen(gds.out)
    
  }else{
    
    seqVCF2GDS(vcfs, gds.out,  info.import=info.import, fmt.import=fmt.import)
    
    gds <- seqOpen(gds.out)
    
    if (length(vcfs) > 1){
      var.ids <- seqGetData(gds, "variant.id")
      
      seqClose(gds)
      setVariantID(gds.out, seq_along(var.ids))
      gds <- seqOpen(gds.out)
    }
    
  }
  
  n <- index.gdsn(gds, "annotation/info/CSQ", silent = TRUE)
  hdr <- get.attr.gdsn(n)$Description
  nms <- unlist(strsplit(strsplit(hdr, "Format: ")[[1]][2], "\\|"))
  
  csq.res <- .csq.to.df(gds, nms)
  
  #add in the variant-level data
  
  which.biallelic <- nAlleles(gds)[csq.res$variant_id] == 2
  
  alt.list <- as(alt(gds)[csq.res$variant_id], "CharacterList")
  
  csq.res$ALT <- NA
  
  csq.res$ALT[which.biallelic] <- unlist(alt.list[which.biallelic])
  csq.res$ALT[which.biallelic == F] <- mapply(function(alt, num){
    return(alt[num])
  }, alt.list[which.biallelic == F], as.integer(csq.res$ALLELE_NUM)[which.biallelic == F])
  
  
  #add in the ref data
  csq.res$REF <- as.character(ref(gds)[csq.res$variant_id])
  
  csq.res <- csq.res[,colnames(csq.res) != "Allele"]
  
  #add in the range data
  
  var.grange <- as.data.frame(granges(gds))[csq.res$variant_id,]
  
  csq.res <- cbind(var.grange, csq.res)
  
  #add in the info data as well
  
  csq.res <- cbind(csq.res, .info(gds, info.import=info.import)[csq.res$variant_id,])
  
  #add in genotypes
  
  #actually a matrix...
  geno.ar <- getGenotype(gds)
  
  melt.geno <- melt(t(geno.ar))
  sub.melt.geno <- melt.geno[is.na(melt.geno$value) == F & melt.geno$value %in% c(".","0/0") == F,]
  
  #for each unique genotype allele, assign each sample by variant and ALLELE_NUM 
  
  split.genos <- strsplitAsListOfIntegerVectors(as.character(sub.melt.geno$value), sep="/")
  
  geno.mat <- t(sapply(split.genos, "["))
  
  sub.melt.geno <- cbind(sub.melt.geno, allele=geno.mat)
  
  allele.dta <- melt(sub.melt.geno, measure.vars=c("allele.1", "allele.2"), value.name="ALLELE_NUM")
  allele.dta <- allele.dta[allele.dta$ALLELE_NUM != 0,]
  allele.dta <- allele.dta[,names(allele.dta) %in% c("variable", "value") == F]
  allele.dta$sample <- as.character(allele.dta$sample)
  names(allele.dta)[1] <- "variant_id"
  
  al.num <- aggregate(ALLELE_NUM~variant_id+sample, function(x) ifelse(length(unique(x)) == 1 && length(x) == 2, 2, 1), data=allele.dta)
  names(al.num)[ncol(al.num)] <- "allele_count"
  
  allele.dta <- merge(allele.dta, al.num, by=c("variant_id", "sample"), all=F, sort=F)
  
  nd.alleles <- allele.dta[!duplicated(allele.dta),]
  
  stopifnot(is.integer(csq.res$variant_id) && is.character(csq.res$ALLELE_NUM))
  stopifnot(is.integer(nd.alleles$variant_id) && is.integer(nd.alleles$ALLELE_NUM))
  
  csq.res$ALLELE_NUM <- as.integer(csq.res$ALLELE_NUM)
  
  csq.allele <- merge(csq.res, nd.alleles, by=c("variant_id", "ALLELE_NUM"), all=F, sort=F)
  
  stopifnot(nrow(csq.allele) == nrow(nd.alleles))
  
  #get count data
  count.ar <- .getVariableLengthData(gds, "annotation/format/AD")
  
  melt.ar <- melt(count.ar, na.rm=T)
  
  total.count <- aggregate(value~sample+variant, sum, data=melt.ar)
  
  total.count$sample <- as.character(total.count$sample)
  
  csq.allele.t <- merge(csq.allele, total.count, by.x=c("variant_id", "sample"), by.y=c("variant", "sample"), all.x=T, all.y=F, sort=F)
  names(csq.allele.t)[ncol(csq.allele.t)] <- "total_reads"
  
  melt.ar$sample <- as.character(melt.ar$sample)
  melt.ar$n <- melt.ar$n-1
  
  fin.csq <- merge(csq.allele.t, melt.ar, by.x=c("variant_id", "ALLELE_NUM", "sample"), by.y=c("variant", "n", "sample"), all.x=T, all.y=F, sort=F)
  names(fin.csq)[ncol(fin.csq)] <- "allele_reads"
  
  seqClose(gds)
  
  if(keep.gds==F){
    file.remove(gds.out)
  }
  
  fin.csq <- fin.csq[,c("variant_id", "seqnames", "start", "end", "width", "REF", "ALT", "sample", "allele_count", "total_reads", "allele_reads", 
                        setdiff(info.import, "CSQ"), setdiff(nms, "Allele"))]
  
  fin.csq[is.na(fin.csq) | fin.csq == ""] <- NA
  
  return(fin.csq)
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


