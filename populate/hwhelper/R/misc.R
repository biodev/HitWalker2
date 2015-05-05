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
