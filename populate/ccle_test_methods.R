#uses Rneo4j https://github.com/nicolewhite/RNeo4j--should add to hwhelper

#library(devtools)
#devtools::install_github("nicolewhite/RNeo4j")



#These generics should be added to hwhelper
setGeneric("subjectAttrs", def=function(obj,...) standardGeneric("subjectAttrs"))
setGeneric("subjectSubset", def=function(obj,...) standardGeneric("subjectSubset"))
setGeneric("findHits", def=function(obj,...) standardGeneric("findHits"))


process_matrix_graph <- function(mm_file_base, string_conf){ 
    
    require(Matrix)
    require(igraph)
    
    init.mat <- readMM(paste0(mm_file_base, ".mtx"))
    
    new.graph.names <- read.delim(paste0(mm_file_base, ".names"), header=FALSE, stringsAsFactors=FALSE)
    
    colnames(init.mat) <- new.graph.names[,1]
    
    cur.graph <- graph.adjacency(init.mat,mode="directed", weighted="score")
    
    new.graph <- subgraph.edges(cur.graph, E(cur.graph)[ score > string_conf ], delete.vertices = TRUE)
    
    new.graph <- as.undirected(new.graph, mode="collapse", edge.attr.comb="max")
    
    return(new.graph)    
}

#sub_graph <- process_matrix_graph("/var/www/hitwalker2_inst/static/network/data/9606.protein.links.v9.1.mm", .95)
#get_direct_connections(sub_graph, "Signaling by EGFR in Cancer (reactome)", "Pathway")
get_direct_connections <- function(use_graph, node_set, node_types=c("Gene", "Pathway")){
    
    require(RNeo4j)
    
    node_types = match.arg(node_types)
    
    graph = startGraph("http://localhost:7474/db/data/")
    
    if (node_types == "Pathway"){
      
      genes <- cypher(graph, paste0('MATCH (path:Pathway)-[:PATHWAY_CONTAINS]->(gene) WHERE path.name="',node_set,'" RETURN gene.name'))[,1]
      
    }else{
      
      gene.symbs <- cypher(graph, paste0("MATCH (n)-[:REFERRED_TO]-(m)  n.name AS gene, m.name AS symbol"))
      
      sub.genes <- gene.symbs[gene.symbs$symbol %in% node_set,]
      
      stopifnot(nrow(sub.genes) == length(node_set))
      
      genes <- sub.genes$gene
    }
    
    string.name <- cypher(graph, "MATCH (n)-[:MAPPED_TO]-(o)-[:REFERRED_TO]-(m) RETURN n.name AS string, o.name AS gene, m.name AS symbol")
   
    node.comb <- data.frame(t(combn(genes, 2)), stringsAsFactors=F)
    
    names(node.comb) <- c("from", "to")
    
    should.keep <- sapply(1:nrow(node.comb), function(x){
        
        from.name <- string.name$string[string.name$gene == node.comb$from[x]]
        to.name <- string.name$string[string.name$gene == node.comb$to[x]]
        
        stopifnot(length(from.name) == 1 && length(to.name) == 1)
        
        are.connected(use_graph,V(use_graph)[name == from.name], V(use_graph)[name == to.name])
    })
   
   direct.cons <- node.comb[should.keep,]
   
   merged.from <- merge(direct.cons, string.name, by.x="from", by.y="gene")
   merged.ft <- merge(merged.from, string.name, by.x="to", by.y="gene")
   
   ret.edges <- merged.ft[,c("symbol.x", "symbol.y")]
   names(ret.edges) <- c("from", "to")
   
   ret.edges <- ret.edges[!duplicated(ret.edges),]
   
   return(ret.edges)
}

get_gene_connections <- function(use.graph, seeds, targs){
    
    #use.graph <- process_matrix_graph("/var/www/hitwalker2_inst/static/network/data/9606.protein.links.v9.1.mm", .4)
    #seeds =c('MAP2K7', 'ALK', 'HSP90AA1')
    #targs=c('MAP3K13', 'MAP3K1', 'KRAS')
    
    require(RNeo4j)
    
    graph = startGraph("http://localhost:7474/db/data/")
    string.name <- cypher(graph, "MATCH (n)-[:MAPPED_TO]-()-[:REFERRED_TO]-(m) RETURN n.name AS string, m.name AS symbol")
    
    gene.grid <- expand.grid(list(seeds=seeds, targs=targs), stringsAsFactors = FALSE)
    
    string.graph.dta <- do.call(rbind, lapply(1:nrow(gene.grid), function(i){
      
        var.names <- string.name$string[string.name$symbol == gene.grid$targs[i]]
        seed.names <- string.name$string[string.name$symbol == gene.grid$seeds[i]]
        
        if ((length(var.names) == 1 && length(seed.names) == 1) == F){
          stop(paste("ERROR: unexpected lengths:", paste(var.names, collapse=","), paste(seed.names, collapse=",")))
        }
        
        paths <-  get.all.shortest.paths(use.graph, from=V(use.graph)[name %in% var.names], to = V(use.graph)[ name %in% seed.names], mode = "all", weights=NULL)
        
        #just check the first as these should all be the same length...s
        if (length(paths$res) > 0 && length(paths$res[[1]]) < 4){
          
            which.path <- which.max(sapply(paths$res, function(x) sum(E(use.graph, path=x)$score)))
            
            edge.mat <- get.edges(use.graph, E(use.graph, path=paths$res[[which.path]]))
            
            return(data.frame(from=V(use.graph)$name[edge.mat[,1]], to=V(use.graph)$name[edge.mat[,2]], stringsAsFactors=F))
        
        }else{
            return(NULL)
        }
        
    }))
    
   #then add in all the direct connections between these nodes (and the initial set if not present)
   
   node_set <- union(unlist(string.graph.dta), string.name$string[string.name$symbol %in% c(gene.grid$seeds, gene.grid$targs)])
   
   node.comb <- data.frame(t(combn(node_set, 2)), stringsAsFactors=F)
   names(node.comb) <- c("from", "to")
   
   should.keep <- sapply(1:nrow(node.comb), function(x){
        are.connected(use.graph,V(use.graph)[name == node.comb[x,1]], V(use.graph)[name == node.comb[x,2]])
   })
   
   direct.cons <- node.comb[should.keep,]
   
   all.edges <- rbind(string.graph.dta, direct.cons)
   
   all.edges <- all.edges[!duplicated(all.edges),]
   
   merged.from <- merge(all.edges, string.name, by.x="from", by.y="string")
   merged.ft <- merge(merged.from, string.name, by.x="to", by.y="string")
   
   ret.edges <- merged.ft[,c("symbol.x", "symbol.y")]
   names(ret.edges) <- c("from", "to")
   
   return(ret.edges)
   
   #return the graph in the form: data.frame(from=, to=)
    
}

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

# structure(list(Subject = c("HEPG2_LIVER", "HEPG2_LIVER", "HEPG2_LIVER",
# "HEPG2_LIVER", "HEPG2_LIVER", "HEPG2_LIVER", "HEPG2_LIVER"), 
# Gene = c("4214", "238", "3320", "5609", "3845", "9175", "4214"), 
# IsHit = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), Datatype = c("Expression",
# "GeneScore", "GeneScore", "GeneScore", "Variants", "Variants",
# "Variants")), .Names = c("Subject", "Gene", "IsHit", "Datatype"))

encode_groups <- function(hit.dta, ids.to.symbs=F, prioritized_subject=NULL, group.by=c("None", "Subject", "Gene"), dense.datatype=NULL){
    
    require(RNeo4j)
    
    group.by <- match.arg(group.by)
  
    graph = startGraph("http://localhost:7474/db/data/")
    
    if (class(hit.dta) == "list"){
      hit.dta <- as.data.frame(hit.dta)
    }
    
    if (missing(dense.datatype) || is.null(dense.datatype) || all(is.na(dense.datatype))){
      dense.datatype <- character(0)
    }
    
    print(dense.datatype)
    
    hit.dta <- hit.dta[!(hit.dta$IsHit == F & hit.dta$Datatype %in% dense.datatype == T),]
    
    hit.dta$FixedDt <- ifelse(hit.dta$IsHit == T, paste0("Observed_", hit.dta$Datatype), paste0("Possible_",hit.dta$Datatype))
                              
    if ((missing(prioritized_subject) || is.null(prioritized_subject) || all(is.na(prioritized_subject)))==F){
      
      #if the Subject is equal to the specified prioritized_subject, the Datatype == "Variants" and the Gene is annotated in string
      #then set FixedDt equal to 'Ranked_Variants'
      
      string.name <- cypher(graph, "MATCH (n) WHERE (n)-[:MAPPED_TO]->() RETURN n.name AS gene")
      
      hit.dta$FixedDt[hit.dta$IsHit == T & hit.dta$Datatype == "Variants" & as.character(hit.dta$Subject) %in% prioritized_subject & as.character(hit.dta$Gene) %in% string.name$gene] <- "Ranked_Variants"
      
    }
    
    sum.dta <- aggregate(FixedDt~Subject + Gene, paste, collapse=",", data=hit.dta)
    
    if (ids.to.symbs){
      
      gene.name <- cypher(graph, "MATCH (n:EntrezID)-[r:REFERRED_TO]-(m) RETURN n.name AS Gene, m.name AS symbol")
      
      sum.dta$Gene <- as.character(sum.dta$Gene)
      
      sum.dta <- merge(sum.dta, gene.name, by="Gene", all.x=T, all.y=F, sort=F)
      
      sum.dta$Gene <- sum.dta$symbol
      sum.dta <- sum.dta[,-which(names(sum.dta) == "symbol")]
    }
    
    if (group.by != "None"){
      
      lo.cols <- paste(setdiff(names(sum.dta), group.by), collapse="+")
      
      sum.dta <- aggregate(as.formula(paste(group.by, lo.cols, sep="~")), c , data=sum.dta)
      
      #only those values (Genes/Subjects) which are not seen outside of the group rows should be collapsed to groups
      
      is.group <- sapply(sum.dta[,group.by], length)
      
      if (any(is.group > 1)){
        
        which.group <- which(is.group > 1)
        
        should.keep <- sapply(names(which.group), function(x){
            if (any(sum.dta[,group.by][[x]] %in% unlist(sum.dta[,group.by][setdiff(names(sum.dta[,group.by]), x)]))){
              return(F)
            }else{
              return(T)
            }
        })
        
        for(i in names(which.group)[should.keep]){
          sum.dta[[group.by]][[i]] <- paste0(group.by, " (", length(sum.dta[[group.by]][[i]]), ")")
        }
        
        for(i in names(which.group)[should.keep==F]){
          which.row <- which(names(sum.dta[,group.by])==i)
          new.dta <- sum.dta[which.row,]
          for(j in sum.dta[[group.by]][[i]][-1]){
            new.dta[,group.by] <- j
            sum.dta <- rbind(sum.dta, new.dta)
          }
        }
        
        sum.dta[,group.by] <- sapply(sum.dta[,group.by], "[", 1)
        
      }
    
    }
      
    return(sum.dta) 
}

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
