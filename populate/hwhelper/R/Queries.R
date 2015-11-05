#taken from http://dec2014.archive.ensembl.org/info/genome/variation/predicted_data.html#consequences
ens.cons.tab <- "transcript_ablation 	A feature ablation whereby the deleted region includes a transcript feature 	SO:0001893 	Transcript ablation
splice_donor_variant 	A splice variant that changes the 2 base region at the 5' end of an intron 	SO:0001575 	Essential splice site
splice_acceptor_variant 	A splice variant that changes the 2 base region at the 3' end of an intron 	SO:0001574 	Essential splice site
stop_gained 	A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript 	SO:0001587 	Stop gained
frameshift_variant 	A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three 	SO:0001589 	Frameshift coding
stop_lost 	A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript 	SO:0001578 	Stop lost
initiator_codon_variant 	A codon variant that changes at least one base of the first codon of a transcript 	SO:0001582 	Non synonymous coding
transcript_amplification 	A feature amplification of a region containing a transcript 	SO:0001889 	Transcript amplification
inframe_insertion 	An inframe non synonymous variant that inserts bases into in the coding sequence 	SO:0001821 	Non synonymous coding
inframe_deletion 	An inframe non synonymous variant that deletes bases from the coding sequence 	SO:0001822 	Non synonymous coding
missense_variant 	A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved 	SO:0001583 	Non synonymous coding
splice_region_variant 	A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron 	SO:0001630 	Splice site
incomplete_terminal_codon_variant 	A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed 	SO:0001626 	Partial codon
stop_retained_variant 	A sequence variant where at least one base in the terminator codon is changed, but the terminator remains 	SO:0001567 	Synonymous coding
synonymous_variant 	A sequence variant where there is no resulting change to the encoded amino acid 	SO:0001819 	Synonymous coding
coding_sequence_variant 	A sequence variant that changes the coding sequence 	SO:0001580 	Coding unknown
mature_miRNA_variant 	A transcript variant located with the sequence of the mature miRNA 	SO:0001620 	Within mature miRNA
5_prime_UTR_variant 	A UTR variant of the 5' UTR 	SO:0001623 	5prime UTR
3_prime_UTR_variant 	A UTR variant of the 3' UTR 	SO:0001624 	3prime UTR
non_coding_transcript_exon_variant 	A sequence variant that changes non-coding exon sequence in a non-coding transcript 	SO:0001792 	Within non coding gene
intron_variant 	A transcript variant occurring within an intron 	SO:0001627 	Intronic
NMD_transcript_variant 	A variant in a transcript that is the target of NMD 	SO:0001621 	NMD transcript
non_coding_transcript_variant 	A transcript variant of a non coding RNA gene 	SO:0001619 	Within non coding gene
upstream_gene_variant 	A sequence variant located 5' of a gene 	SO:0001631 	Upstream
downstream_gene_variant 	A sequence variant located 3' of a gene 	SO:0001632 	Downstream
TFBS_ablation 	A feature ablation whereby the deleted region includes a transcription factor binding site 	SO:0001895 	Tfbs ablation
TFBS_amplification 	A feature amplification of a region containing a transcription factor binding site 	SO:0001892 	Tfbs amplification
TF_binding_site_variant 	A sequence variant located within a transcription factor binding site 	SO:0001782 	Regulatory region
regulatory_region_ablation 	A feature ablation whereby the deleted region includes a regulatory region 	SO:0001894 	Regulatory region ablation
regulatory_region_amplification 	A feature amplification of a region containing a regulatory region 	SO:0001891 	Regulatory region amplification
regulatory_region_variant 	A sequence variant located within a regulatory region 	SO:0001566 	Regulatory region
feature_elongation 	A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence 	SO:0001907 	Feature elongation
feature_truncation 	A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence 	SO:0001906 	Feature truncation
intergenic_variant 	A sequence variant located in the intergenic region, between genes 	SO:0001628 	Intergenic"

my.consequence.order <- function(){
  
  cons.lines <- strsplit(ens.cons.tab, "\\n")
  
  return(sapply(strsplit(cons.lines[[1]], "\\s+"), "[", 1))
}

my.short.prots <- c(Ala="A", Arg="R", Asn="N", Asp="D", Asx ="B", Cys="C", Glu="E", Gln ="Q", Glx="Z", Gly="G", His= "H", Ile ="I", Leu= "L",
                    Lys ="K", Met= "M", Phe ="F", Pro= "P", Ser= "S", Thr= "T", Trp= "W", Tyr= "Y", Val= "V", Xxx ="X", Ter ="*" )

.fix.protein.ids <- function(csq.df){
  
  #approach adapted from vcf2maf
  
  proteins <- as.character(csq.df$HGVSp)
  
  # Remove transcript ID from HGVS codon/protein changes, to make it easier on the eye
  proteins <- sub("^.*:", "", proteins, perl=T)
  
  # Remove the prefixed HGVSc code in HGVSp, if found
  proteins <- sapply(regmatches(proteins, regexec("\\(*p\\.\\S+\\)*",proteins)), "[", 1)
  
  proteins <- gsub("[\\(\\)]", "", proteins)
  
  # Create a shorter HGVS protein format using 1-letter codes
  
  for (i in names(my.short.prots)){
    proteins <- gsub(i, my.short.prots[i], proteins)
  }
  
  # Fix HGVSp_Short,for splice acceptor/donor variants
  
  splice.pos <- csq.df$Variant_Classification %in% c("splice_acceptor_variant", "splice_donor_variant")
  
  if (sum(splice.pos) > 0){
    
    c.pos <- as.numeric(sapply(regmatches(as.character(csq.df$HGVSc), regexec("c\\.(\\d+)", as.character(csq.df$HGVSc))), "[", 2))
    
    c.pos <- ifelse(is.na(c.pos) ==F & c.pos < 0, 1, c.pos)
    
    proteins <- ifelse((splice.pos == T) & (is.na(c.pos) ==F) , paste0("p.X", sprintf("%.0f", (c.pos + c.pos %% 3)/3), "_splice"), proteins)
  }
  
  # Fix HGVSp_Short for Silent mutations, so it mentions the amino-acid and position
  
  p.pos <- as.numeric(sapply(regmatches(as.character(csq.df$Protein_position), regexec("^(\\d+)\\/\\d+$", as.character(csq.df$Protein_position))), "[", 2))
  
  fin.prots <- ifelse(csq.df$Variant_Classification == "synonymous_variant", paste0("p.", csq.df$Amino_acids, p.pos,csq.df$Amino_acids) , proteins)
  
  csq.df$HGVSp_Short <- fin.prots
  
  return(csq.df)
  
}


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
  
  vcf.dta$Cons_Summary <- ifelse(vcf.dta$Variant_Classification %in% my.consequence.order()[1:12], "nonsynonymous", "synonymous")
  
  vcf.dta <- vcf.dta[,c("sample", "gene", "name", "seqnames", "start", "end", "REF", "ALT", "Protein_position", "Amino_acids", "allele_count", "allele_reads", 
                        "total_reads", "MQ0", "FS", "MQ", "QD", "SB","Existing_variation", "Variant_Classification", "Cons_Summary" ,"HGVSp", "GMAF", "Feature", "SIFT", "PolyPhen", 
                        "cohort_count", "cohort_freq", "in_dbsnp", "in_esp", "in_cosmic")]
  
  vcf.dta <- factors.to.chars(vcf.dta)
  
  return(new("DenseNeoData", data=vcf.dta, node.name=node.name, sample.edge.name=sample.edge.name))
  
}

.csq.to.df <- function(gds, header, pick.only=T){
  
  csq.dta.list <- seqGetData(gds, "annotation/info/CSQ")
  
  exp.inds <- unlist(lapply(seq_along(csq.dta.list$length), function(x) rep(x, csq.dta.list$length[x])))
  
  raw <- strsplit(csq.dta.list$data, "\\|")
  
  csq <- matrix(nrow = length(raw), ncol = length(header))
  
  for (i in 1:nrow(csq)) csq[i, 1:length(raw[[i]])] <- raw[[i]]
  colnames(csq) <- header
  
  var.ind <- index.gdsn(gds, "variant.id", silent = TRUE)
  variant.id <- read.gdsn(var.ind,simplify="none")
  
  csq <- data.frame(variant_id=variant.id[exp.inds],csq, stringsAsFactors = F)
  
  if (pick.only){
    csq <- csq[is.na(csq$PICK) == F & csq$PICK == 1,]
  }
  
  csq <- csq[order(csq$variant_id, decreasing=F),]
  
  return(csq)
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

default.protein.filters <- function(dta){
  
  cons.order <- my.consequence.order()
  
  ret.dta <- dta[(is.na(dta$HGVSp) == F) & (is.na(dta$Variant_Classification) == F) & (dta$Variant_Classification %in% cons.order[1:12] == T),]
  
  return(ret.dta)
}


#' @rdname class_helpers
#' @param vcfs The path to VCF file(s) that have been annotated by Ensembl VEP using the '--refseq' option.  
#' Our typical workflow additionally involves subsetting to protein impacting variants 
#' and keeping only consequences chosen via the allele-specific 'PICK' column.
#' @param keep.gds Should the intermediate GDS file be kept after the \code{data.frame} is generated?
make.vcf.table <- function(vcfs, info.import=c("FS", "MQ0", "MQ", "QD", "SB", "CSQ"), keep.samples=NULL, ignore.genotype=F, readcount.import="AD", gds.out="variant.gds", rm.zero.na=T, keep.gds=F, filter.by=default.protein.filters){
  
  if (length(readcount.import) %in% c(1, 2) == F){
    stop("ERROR: Expect readcount.import to be of length 1 or 2")
  }
  
  all.samples <- unique(unlist(lapply(vcfs, seqVCF.SampID)))
  
  if (missing(keep.samples) || is.null(keep.samples) || all(is.na(keep.samples))){
    keep.samples <- all.samples
  }else if ((is.character(keep.samples) == F) || (all(keep.samples %in% all.samples)) == F){
    stop("Unexpected input for 'keep.samples'")
  }
  
  if(file.exists(gds.out)){
    
    message(paste("File:", gds.out, "exists, utilizing it for this analysis"))
    
    gds <- seqOpen(gds.out)
    
  }else{
    
    seqVCF2GDS(vcfs, gds.out,  info.import=info.import, fmt.import=readcount.import)
    
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
  
  if (sum(which.biallelic==F) > 0){
    csq.res$ALT[which.biallelic == F] <- mapply(function(alt, num){
      return(alt[num])
    }, alt.list[which.biallelic == F], as.integer(csq.res$ALLELE_NUM)[which.biallelic == F])
  }
  
  #add in the ref data
  csq.res$REF <- as.character(ref(gds)[csq.res$variant_id])
  
  csq.res <- csq.res[,colnames(csq.res) != "Allele"]
  
  #add in the range data
  
  var.grange <- as.data.frame(granges(gds))[csq.res$variant_id,]
  
  csq.res <- cbind(var.grange, csq.res)
  
  #add in the info data as well
  
  if (any(info.import %in% names(csq.res))){
    message(paste("Removing CSQ column with specified info name", paste(intersect(info.import, names(csq.res)), collapse=",")))
    csq.res <- csq.res[,-which(names(csq.res) %in% info.import)]
  }
  
  csq.res <- cbind(csq.res, .info(gds, info.import=info.import)[csq.res$variant_id,])
  
  if (!ignore.genotype){
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
  
  }else if ((ignore.genotype == T) && (all(which.biallelic)==F)){
    
    stop("ERROR: Cannot use 'ignore.genotypes' when multi-allelic variants are present")
    
  }else{
    
    #add in the sample information, assume that each sample has each row...
    
    csq.allele <- do.call(rbind, lapply(seqGetData(gds, "sample.id"), function(x){
      csq.res$sample <- x
      csq.res
    }))
    
    csq.allele$allele_count <- 1
    csq.allele$ALLELE_NUM <- as.integer(csq.allele$ALLELE_NUM)
  }
  
  #get count data
  
  if (length(readcount.import) == 1){
    #assume it contains both counts
    
    count.ar <- getVariableLengthData(gds, paste0("annotation/format/", readcount.import))
    
    melt.ar <- melt(count.ar, na.rm=T)
    
  }else{
    #assume counts are provided in the form reference,alternative: readcount.import[1:2]
    
    reference <- getVariableLengthData(gds, paste0("annotation/format/", readcount.import[1]))
    
    melt.ref <- melt(reference, na.rm=T)
    melt.ref <- cbind(n=1, melt.ref)
    
    alt <- getVariableLengthData(gds, paste0("annotation/format/", readcount.import[2]))
    
    melt.alt <- melt(alt, na.rm=T)
    melt.alt <- cbind(n=2, melt.alt)
    
    melt.ar <- rbind(melt.ref, melt.alt)
    
  }
  
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
  
  fin.csq <- fin.csq[fin.csq$sample %in% keep.samples,]
  
  fin.csq <- fin.csq[,c("variant_id", "seqnames", "start", "end", "width", "REF", "ALT", "sample", "allele_count", "total_reads", "allele_reads", 
                        setdiff(info.import, "CSQ"), setdiff(nms, c("Allele", info.import)))]
  
  fin.csq[is.na(fin.csq) | fin.csq == ""] <- NA
  
  use.cons <- my.consequence.order()
  
  fin.csq$Variant_Classification <- sapply(strsplit(as.character(fin.csq$Consequence), "&"), function(x){
    
    if (length(x) > 1){
      
      new.x <- factor(x, levels=use.cons, ordered=T)
      return(as.character(sort(new.x)[1]))
    }else{
      return(x)
    }
    
  })
  
  message("applying post-summarization filters")
  
  fin.csq <- filter.by(fin.csq)
  
  if (rm.zero.na){
    
    message("removing rows with 0's in allele_reads")
  
    fin.csq <- fin.csq[is.na(fin.csq$allele_reads) == F & fin.csq$allele_reads > 0,]
  }
  
  fin.csq <- .fix.protein.ids(fin.csq)
  
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


