setMethod("getFrequency", signature("HW2Config"), function(obj,datatype, subset, subset_type=c("Subject", "Subject_Category", "Gene")){
    
    subset_type = match.arg(subset_type)
    
    .get_freq_by_subset <- function(obj, datatype, subset){
        
        subj <- obj@subject
        
        which.subset <- apply(subj@subject.info[,-c(1:2)] == subset, 1, any)
        
        getFrequency(obj@data.list[[datatype]], subj@subject.info[which.subset,1], "Subject")
    }
    
    switch(subset_type,
           Subject=getFrequency(obj@data.list[[datatype]], subset, subset_type),
           Gene=getFrequency(obj@data.list[[datatype]], subset, subset_type),
           Subject_Category=.get_freq_by_subset(obj, datatype, subset)
           )
    
})


setMethod("getFrequency", signature("HW2exprSet"), function(obj, subset, type){
    
    
    print('hello')
})

setMethod("getFrequency", signature("CCLEMaf"), function(obj, subset, type){
    
    if (type == "Subject"){
        
        sub.maf <- obj@maf[obj@maf$Tumor_Sample_Barcode %in% subset,]
        
        gene.agg <- aggregate(Tumor_Sample_Barcode~Entrez_Gene_Id, function(x) length(unique(x)) ,data=sub.maf)
        
        gene.freq <- as.data.frame(table(gene.agg$Tumor_Sample_Barcode))
        
        ord.freq <- order(gene.freq$Var1, decreasing=T)
        
        names(gene.freq) <- c("Frequency","Genes")
        
        gene.freq$Frequency <- paste0(round((as.numeric(as.character(gene.freq$Frequency))/length(subset))*100), "%")
        
        return(gene.freq[ord.freq, c("Genes", "Frequency")])
        
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

setMethod("getFrequency", signature("DrugMatrix"), function(obj, subset, type){
    
    gene.score <- drug.gene.score(obj)
    
    if (type == "Subject"){
        
        
        
        
    }else if (type == "Gene"){
        
        
    }else{
        stop("ERROR: Expected either Subject or Gene as a type")
    }
})

