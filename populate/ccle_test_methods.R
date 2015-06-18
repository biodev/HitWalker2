#uses Rneo4j https://github.com/nicolewhite/RNeo4j--should add to hwhelper

#library(devtools)
#devtools::install_github("nicolewhite/RNeo4j")



#These generics should be added to hwhelper
setGeneric("subjectAttrs", def=function(obj,...) standardGeneric("subjectAttrs"))
setGeneric("subjectSubset", def=function(obj,...) standardGeneric("subjectSubset"))
setGeneric("findHits", def=function(obj,...) standardGeneric("findHits"))


#returns a vector of subject names in prinical to be part of a given metanode
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
    gene.name <- cypher(graph, "MATCH (n:EntrezID)-[r:REFFERED_TO]-(m) RETURN n.name AS gene, m.name AS symbol")
    
    stopinfot(gene_type == "Gene")
    
    genes <- gene.name$gene[gene.name$symbol %in% genes]
    
    do.call(rbind, lapply(names(hw2.conf@data.list), function(x){
        
        subset.samples <- obj@subject@subject.to.sample$sample[(obj@subject@subject.to.sample$type == sample.rel.names[i]) & (obj@subject@subject.to.sample[,1] %in% subjects)]
        
        hit.dta <- findHits(hw2.conf@data.list[[i]], subset.samples, genes)
        
        hit.dta$Datatype <- x
        
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

setMethod("findHits", signature("HW2exprSet"), function(obj, samples, genes=NULL){
    
    print ('expression method')
    
    annot.pack <- annotation(obj@exprs)
    
    if (grepl("\\.db$", annot.pack) == F){
        annot.pack <- paste0(annot.pack, ".db")
    }
    
    require(annot.pack, character.only=T)
    require(reshape2)
    
    #below is kind of a slow process so will save it and reuse if available
    
    #if (file.exists("saved_exprs.RData")){
    #    
    #    load("saved_exprs.RData")
    #    
    #}else{
    
        mapping <- select(get(annot.pack), columns="ENTREZID", featureNames(obj@exprs))
    
        if (missing(genes) || is.null(genes) || all(is.na(genes))){
        
            sub.mapping <- mapping
        
        }else{
            
            sub.mapping <- mapping[mapping$ENTREZID %in% genes,]
        }
        
        print (samples)
        common.samples <- intersect(samples, sampleNames(obj@exprs))
        print (common.samples)
        
        probeset.dta <- melt(exprs(obj@exprs)[,common.samples], as.is=T)
        
        probeset.gene <- merge(probeset.dta, sub.mapping, by.x="Var1", by.y="PROBEID", sort=F)
        
        sum.gene <- aggregate(value~Var2+ENTREZID, max, data=probeset.gene)
        
        #save(sum.gene, file="saved_exprs.RData")
    #}
    
    sum.gene <- sum.gene[sum.gene$value > obj@default,]
    
    names(sum.gene) <- c("Sample", "Gene")
    
    return(sum.gene)
})

setMethod("findHits", signature("CCLEMaf"), function(obj, samples, genes=NULL){
    
    if (missing(genes) || is.null(genes) || all(is.na(genes))){
        gene.maf <- obj@maf    
    }else{
        gene.maf <- obj@maf[obj@maf$Entrez_Gene_Id %in% genes,]
    }
    
    sub.maf <- gene.maf[gene.maf$Tumor_Sample_Barcode %in% samples,c("Tumor_Sample_Barcode", "Entrez_Gene_Id")]
    
    names(sub.maf) <- c("Sample", "Gene")
    
    return(sub.maf)
    
})

setMethod("findHits", signature("DrugMatrix"), function(obj, samples, genes=NULL){
    
    require(reshape2)
    
    drug.dta <- melt(obj@matrix, as.is=T)
    
    median.ic50 <- aggregate(value~Var1, median, data=drug.dta)
    names(median.ic50)[2] <- "median"
    
    drug.dta.merge <- merge(drug.dta, median.ic50, by="Var1")
    
    drug.dta.merge$is.hit <- with(drug.dta.merge, value <= (median*.2))
    
    drug.dta.genes <- merge(drug.dta.merge, obj@mapping, by.x="Var1", by.y="drug")
    
    drug.dta.genes$score <- with(drug.dta.genes, ifelse(is.hit, weight, -weight))
    
    drug.sum <- aggregate(score~Var2+gene, sum, data=drug.dta.genes)
    
    names(drug.sum) <- c("Sample", "Gene", "score")
    
    use.hits <- drug.sum[drug.sum$score > obj@default,]
    
    use.hits <- use.hits[,c("Sample", "Gene")]
    
    return(use.hits)
})

#setMethod("getFrequency", signature("CCLEMaf"), function(obj, subset, type, sample.mapping){
#    
#    print ('variant method')
#    
#    maf <- findHits(obj, subset)
#    
#    if (type == "Subject"){
#        
#        sample.subset <- sample.mapping$sample[sample.mapping[,1] %in% subset]
#        
#        #note that the frequency is of the total metanode size...
#        subset.length <- length(subset)
#        
#        freq.by(maf, "Sample", "Gene", sample.subset, "Genes", subset.length)
#        
#    }else if (type == "Gene"){
#        
#    }else{
#        stop("ERROR: Expected either Subject or Gene as a type")
#    }
#    
#})
#
#setMethod("getFrequency", signature("DrugMatrix"), function(obj, subset, type, sample.mapping){
#    
#    print ('drug method')
#    
#    gene.score <- findHits(obj, subset)
#    
#    if (type == "Subject"){
#    
#        sample.subset <- sample.mapping$sample[sample.mapping[,1] %in% subset]
#    
#        freq.by(gene.score, "subject", "gene", sample.subset, "Genes", length(subset))
#        
#    }else if (type == "Gene"){
#        
#        freq.by(gene.score, "gene", "subject", subset, "Subjects")
#        
#    }else{
#        stop("ERROR: Expected either Subject or Gene as a type")
#    }
#   
#})
#
#setMethod("getFrequency", signature("HW2exprSet"), function(obj, subset, type, sample.mapping){
#    
#    sum.gene <- findHits(obj, subset)
#    
#    if (type == "Subject"){
#    
#        sample.subset <- sample.mapping$sample[sample.mapping[,1] %in% subset]
#        freq.by(sum.gene, "Sample", "Gene", sample.subset, "Genes", length(subset))
#    
#    }else if (type == "Gene"){
#        
#    }else{
#        stop("ERROR: Expected either Subject or Gene as a type")
#    }
#})
