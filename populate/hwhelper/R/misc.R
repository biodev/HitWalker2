
#' Testing functions
#'
#' Miscellaneous functions helpful for checking the sanity of the HitWalker2 web interface.
#'
#' @name testing_functions
#' 
#' @return For \code{get_gene_connections} and \code{get_direct_connections}: A \code{data.frame} of gene symbols with 
#' 'from' and 'to' columns describing the gene interactions.
#'
#' 
NULL


#' @rdname testing_functions
#' 
#' @param mm_file_base Path to the matrix files used by HitWalker2, assumes that the matrix will have a '.mtx' suffix and that there
#' will be an accompanying '.names' file. 
#' @param string_conf Confidence threshold for the protein-protein interaction data.
#' @return For \code{process_matrix_graph}: An \code{igraph} object containing the thresholded graph.
#'
process_matrix_graph <- function(mm_file_base, string_conf){ 
  
  require(Matrix)
  
  init.mat <- readMM(paste0(mm_file_base, ".mtx"))
  
  new.graph.names <- read.delim(paste0(mm_file_base, ".names"), header=FALSE, stringsAsFactors=FALSE)
  
  colnames(init.mat) <- new.graph.names[,1]
  
  cur.graph <- graph.adjacency(init.mat,mode="directed", weighted="score")
  
  new.graph <- subgraph.edges(cur.graph, E(cur.graph)[ score > string_conf ], delete.vertices = TRUE)
  
  new.graph <- as.undirected(new.graph, mode="collapse", edge.attr.comb="max")
  
  return(new.graph)    
}

#' @rdname testing_functions
#'
#' @param use_graph \code{igraph} object containing a protein-protein interaction graph
#' @param node_set A vector of gene names or a pathway name as specified in \code{node_types}
#' @param node_types Either 'Gene' or 'Pathway' indicating what the supplied IDs refer to.
#'
##uses Rneo4j https://github.com/nicolewhite/RNeo4j--should add to hwhelper
#library(devtools)
#devtools::install_github("nicolewhite/RNeo4j")
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


#' @rdname testing_functions
#'
#' @param seeds The gene symbols indiating the 'seeds' or 'hits'.
#' @param targs The gene symbols which were prioritized.
#'
#' 
get_gene_connections <- function(use_graph, seeds, targs){
  
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
    
    paths <-  get.all.shortest.paths(use_graph, from=V(use_graph)[name %in% var.names], to = V(use_graph)[ name %in% seed.names], mode = "all", weights=NULL)
    
    #just check the first as these should all be the same length...s
    if (length(paths$res) > 0 && length(paths$res[[1]]) < 4){
      
      which.path <- which.max(sapply(paths$res, function(x) sum(E(use_graph, path=x)$score)))
      
      edge.mat <- get.edges(use_graph, E(use_graph, path=paths$res[[which.path]]))
      
      return(data.frame(from=V(use_graph)$name[edge.mat[,1]], to=V(use_graph)$name[edge.mat[,2]], stringsAsFactors=F))
      
    }else{
      return(NULL)
    }
    
  }))
  
  #then add in all the direct connections between these nodes (and the initial set if not present)
  
  node_set <- union(unlist(string.graph.dta), string.name$string[string.name$symbol %in% c(gene.grid$seeds, gene.grid$targs)])
  
  node.comb <- data.frame(t(combn(node_set, 2)), stringsAsFactors=F)
  names(node.comb) <- c("from", "to")
  
  should.keep <- sapply(1:nrow(node.comb), function(x){
    are.connected(use_graph,V(use_graph)[name == node.comb[x,1]], V(use_graph)[name == node.comb[x,2]])
  })
  
  direct.cons <- node.comb[should.keep,]
  
  all.edges <- rbind(string.graph.dta, direct.cons)
  
  all.edges <- all.edges[!duplicated(all.edges),]
  
  merged.from <- merge(all.edges, string.name, by.x="from", by.y="string")
  merged.ft <- merge(merged.from, string.name, by.x="to", by.y="string")
  
  ret.edges <- merged.ft[,c("symbol.x", "symbol.y")]
  names(ret.edges) <- c("from", "to")
  
  return(ret.edges)
  
}


#' @rdname testing_functions
#'
#' @param hit.dta A \code{data.frame} as derived from \code{findHits} which contains 'Subject', 'Gene', 'IsHit' and 'Datatype' columns.
#' @param ids.to.symbs Should the input IDs be converted to symbols by looking up their display name?
#' @param prioritized_subject Indicates which of the subjects in the 'Subjects' column was used for prioritiation, can be NULL. 
#' @param group.by Which variable to group by when aggregating.
#' @param dense.datatype For some datatypes such as expression, almost every gene can be considered a 'possible' hit.  This parameter should be
#' used if a datatype is limited to only observed hits.
#'
#' @return For \code{encode_groups} a \code{data.frame} with 'Subject', 'Gene' and 'FixedDt' columns where the latter column is a delimited  list of 
#' relationship types corresponding to the Subject(s) and Gene(s) in question.
#'
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
    
    sum.dta$gene.rel <- apply(sum.dta[,setdiff(names(sum.dta), group.by)], 1, function(x) paste(unlist(x), collapse="-"))
    
    new.dta <- aggregate(as.formula(paste("gene.rel", group.by, sep="~")), paste, collapse=";" , data=sum.dta)
    
    coll.dta <- aggregate(as.formula(paste(group.by, "gene.rel", sep="~")), function(x){
        
        if (length(x) == 1){
          return(x)
        }else{
          return(paste0(group.by, " (", length(x) ,")"))
        }
        
    }, data=new.dta)
    
    ret.dta <- do.call(rbind, lapply(1:nrow(coll.dta), function(x){
        split.rel <- strsplit(coll.dta[x,"gene.rel"], ";")[[1]]
        split.rb <- do.call(rbind, strsplit(split.rel, "-"))
        colnames(split.rb) <- c(setdiff(names(sum.dta), c(group.by , "FixedDt", "gene.rel")), "FixedDt")
        temp.df <- data.frame(coll.dta[x,group.by], split.rb, stringsAsFactors=F)
        names(temp.df)[1] <- group.by
        
        return(temp.df)
    }))
  
  return(ret.dta) 
}
}



read.pc.gmt <- function(filename, organism.code="9606")
{
    org.pack <- switch(make.names(organism.code), X9606="org.Hs.eg.db")
    
    if (!require(org.pack, character.only=T)){
        stop(paste("ERROR: Need to have package", org.pack))
    }
    
    gmt.lines <- readLines(filename)
    gmt.split.lines <- strsplit(gmt.lines, "\t")
    
    path.org <- sapply(gmt.split.lines, function(x) strsplit(x[1], ":")[[1]][1])
    use.paths <- is.na(path.org) == F & path.org == organism.code
    
    use.path.lines <- gmt.split.lines[use.paths]
    #
    path.meta <- lapply(use.path.lines, function(x) c(sub(paste0(organism.code,":\\s+"), "", x[1]), regmatches(x[2], regexec("datasource:\\s+(\\w+);\\s+organism:\\s+(\\d+);\\s+id\\s+type:\\s+(\\w+)",x[2]))[[1]][-1]))
    
    path.meta.dta.list <- lapply(1:length(path.meta), function(x) cbind(matrix(path.meta[[x]], ncol=length(path.meta[[x]]), nrow=length(use.path.lines[[x]])-2, byrow=T), use.path.lines[[x]][-c(1:2)]))
    
    path.meta.dta <- data.frame(do.call("rbind", path.meta.dta.list), stringsAsFactors=F)
    names(path.meta.dta) <- c("pathway", "database", "organism", "id_type", "id")
    
    stopifnot(all(path.meta.dta$id_type == "uniprot") && all(path.meta.dta$organism == "9606"))
    
    #warnings due to multi-way relationships
    suppressWarnings(uni.to.ent <- select(eval(parse(text=org.pack)), keys=as.character(path.meta.dta$id), columns="ENTREZID", keytype="UNIPROT"))
    
    #merge them and return the result
    
    path.meta.dta.merge <- merge(path.meta.dta, uni.to.ent, by.x="id", by.y="UNIPROT", all=T, incomparables=NA)
    
    ret.dta <- path.meta.dta.merge[,c("pathway", "database", "id", "ENTREZID")]
    names(ret.dta)[3:4] <- c("uniprot", "entrezID")
    
    return(ret.dta)
}

#from toupper docs
capwords <- function(s, strict = FALSE) {
         cap <- function(s) paste(toupper(substring(s, 1, 1)),
                       {s <- substring(s, 2); if(strict) tolower(s) else s},
                                  sep = "", collapse = " " )
         sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
     }
     
in.csv.col <- function(vec, search.vals, delim.str=",", match.func=any)
{
    vec <- as.character(vec)
    
    in.vec <- sapply(strsplit(vec, delim.str), function(x)
                     {
                        match.func(x %in% search.vals) 
                     })
    return(in.vec)
}

get.biomart.mapping <- function(host="feb2014.archive.ensembl.org"){
  
  require(biomaRt)
  
  sel.obj <- useMart(host='feb2014.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
  
  #get all the gene IDs
  
  temp <- getBM(mart=sel.obj, attributes="ensembl_gene_id")
  
  ret.dta <- select(sel.obj, keys=temp[,1], columns=c("ensembl_gene_id", "entrezgene"), keytype="ensembl_gene_id")

  names(ret.dta) <- c("Gene", "entrezID")
  
  return(ret.dta)
}

factors.to.chars <- function(dta){
  for (i in colnames(dta)){
    if (is.factor(dta[,i])){
      dta[,i] <- as.character(dta[,i])
    }
  }
  
  return(dta)
}

