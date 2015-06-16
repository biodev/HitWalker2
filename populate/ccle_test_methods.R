setMethod("getFrequency", signature("HW2Config"), function(obj,datatype, subset, subset_type=c("Subject", "Subject_Category", "Gene")){
    
    #should probably add this to the class definition and have it populated using addSamples<-
    sample.rel.names <- c(Expression="Affy_Expression", Variants="DNASeq", GeneScore="Drug_Assay")
    
    subset_type = match.arg(subset_type)
    
    .get_freq_by_subset <- function(obj, datatype, subset, sample.mapping){
        
        subj <- obj@subject
        
        which.subset <- apply(subj@subject.info[,-c(1:2)] == subset, 1, any)
        
        use.subjs <- subj@subject.info[which.subset,1]
        
        getFrequency(obj@data.list[[datatype]],
                     use.subjs,
                     "Subject",
                     sample.mapping)
    }
    
    subset.samples <- obj@subject@subject.to.sample[(obj@subject@subject.to.sample$type == sample.rel.names[datatype]),]
    
    switch(subset_type,
           Subject=getFrequency(obj@data.list[[datatype]], subset.samples, subset_type, subset.samples),
           Gene=getFrequency(obj@data.list[[datatype]], subset.samples, subset_type, subset.samples),
           Subject_Category=.get_freq_by_subset(obj, datatype, subset, subset.samples)
           )
    
})

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

setMethod("getFrequency", signature("HW2exprSet"), function(obj, subset, type, sample.mapping){
    
    print ('expression method')
    annot.pack <- annotation(obj@exprs)
    
    if (grepl("\\.db$", annot.pack) == F){
        annot.pack <- paste0(annot.pack, ".db")
    }
    
    require(annot.pack, character.only=T)
    require(reshape2)
    
    mapping <- select(get(annot.pack), columns="ENTREZID", featureNames(obj@exprs))
    
    probeset.dta <- melt(exprs(obj@exprs), as.is=T)
    
    probeset.gene <- merge(probeset.dta, mapping, by.x="Var1", by.y="PROBEID", sort=F)
    
    sum.gene <- aggregate(value~Var2+ENTREZID, max, data=probeset.gene)
    
    sum.gene <- sum.gene[sum.gene$value > obj@default,]
    
    
    if (type == "Subject"){
    
        sample.subset <- sample.mapping$sample[sample.mapping[,1] %in% subset]
        freq.by(sum.gene, "Var2", "ENTREZID", sample.subset, "Genes", length(subset))
    
    }else if (type == "Gene"){
        
    }else{
        stop("ERROR: Expected either Subject or Gene as a type")
    }
})

setMethod("getFrequency", signature("CCLEMaf"), function(obj, subset, type, sample.mapping){
    
    print ('variant method')
    
    if (type == "Subject"){
        
        sample.subset <- sample.mapping$sample[sample.mapping[,1] %in% subset]
        
        #note that the frequency is of the total metanode size...
        subset.length <- length(subset)
        
        freq.by(obj@maf, "Tumor_Sample_Barcode", "Entrez_Gene_Id", sample.subset, "Genes", subset.length)
        
    }else if (type == "Gene"){
        
    }else{
        stop("ERROR: Expected either Subject or Gene as a type")
    }
    
})

drug.gene.score <- function(dm.obj){
    
    require(reshape2)
    
    drug.dta <- melt(dm.obj@matrix, as.is=T)
    
    median.ic50 <- aggregate(value~Var1, median, data=drug.dta)
    names(median.ic50)[2] <- "median"
    
    drug.dta.merge <- merge(drug.dta, median.ic50, by="Var1")
    
    drug.dta.merge$is.hit <- with(drug.dta.merge, value <= (median*.2))
    
    drug.dta.genes <- merge(drug.dta.merge, dm.obj@mapping, by.x="Var1", by.y="drug")
    
    drug.dta.genes$score <- with(drug.dta.genes, ifelse(is.hit, weight, -weight))
    
    drug.sum <- aggregate(score~Var2+gene, sum, data=drug.dta.genes)
    
    names(drug.sum) <- c("subject", "gene", "score")
    
    return(drug.sum)
}

setMethod("getFrequency", signature("DrugMatrix"), function(obj, subset, type, sample.mapping){
    
    print ('drug method')
    
    gene.score <- drug.gene.score(obj)
    
    gene.score <- gene.score[gene.score$score > obj@default,]
    
    if (type == "Subject"){
    
        sample.subset <- sample.mapping$sample[sample.mapping[,1] %in% subset]
    
        freq.by(gene.score, "subject", "gene", sample.subset, "Genes", length(subset))
        
    }else if (type == "Gene"){
        
        freq.by(gene.score, "gene", "subject", subset, "Subjects")
        
    }else{
        stop("ERROR: Expected either Subject or Gene as a type")
    }
   
})

