#' @rdname test_helpers
setMethod("findHits", signature("HW2exprSet"), function(obj, samples, genes=NULL, limit.to.hits=T){
  
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

##' @rdname test_helpers
# setMethod("findHits", signature("CCLEMaf"), function(obj, samples, genes=NULL, limit.to.hits=T){
#   
#   if (missing(genes) || is.null(genes) || all(is.na(genes))){
#     gene.maf <- obj@maf    
#   }else{
#     gene.maf <- obj@maf[obj@maf$Entrez_Gene_Id %in% genes,]
#   }
#   
#   sub.maf <- gene.maf[gene.maf$Tumor_Sample_Barcode %in% samples,c("Tumor_Sample_Barcode", "Entrez_Gene_Id")]
#   
#   names(sub.maf) <- c("Sample", "Gene")
#   
#   if (limit.to.hits==F){
#     sub.maf$IsHit <- T
#   }
#   
#   return(sub.maf)
#   
# })

#' @rdname test_helpers
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
