#Insulin/IGF pathway-protein kinase B signaling cascade (panther)
#[u'ALK', u'MAP3K1', u'HEPG2_LIVER', u'MAP3K13', u'MAP2K7', u'HSP90AA1', u'KRAS']
#[u'ALK', u'MAP3K1', u'MAP3K7', u'HEPG2_LIVER', u'MAP3K13', u'MAP3K14', u'MAP2K7']

#load("~/Desktop/hitwalker2_paper/ccle_conf.RData")
#temp <- findHits(hw2.conf, 'HEPG2_LIVER', c('ALK', 'MAP3K1', 'MAP3K7', 'HEPG2_LIVER', 'MAP3K13','MAP3K14', 'MAP2K7' ), 'Subject', 'Gene')
#temp.2 <- findHits(hw2.conf, 'HEPG2_LIVER', 'SHC1 events in EGFR signaling (reactome)', 'Subject', 'Pathway')
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

setMethod("findHits", signature("HW2exprSet"), function(obj, samples, genes=NULL, limit.to.hits=T){
    
    print ('expression method')
    
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