read.pc.gmt <- function(filename, organism.code="9606")
{
    org.pack <- switch(make.names(organism.code), X9606="org.Hs.eg.db")
    
    require(org.pack, character.only=T)
    
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
     
vcf.to.graph <- function()
{
    require(VariantAnnotation)
    requireNamespace(ensemblVEP)
    
    vcf.obj <- readVcf("test_out.vcf.new", genome="ncbi")
    
    vcf.csq <- ensemblVEP::parseCSQToGRanges(vcf.obj, VCFRowID=rownames(vcf.obj))
    
    annot.dta <- make.additional.flags(vcf.csq)
    
    #this will form the basic of Variation -> Transcript relationships (IMPACTS)

    #also needs to have a LabID -> Variation mapping (DNA_DIFF)
}


make.additional.flags <- function(annot.grange)
{
  
    colnames(mcols(annot.grange))[colnames(mcols(annot.grange)) == "HGNC"] <- "Symbol"
    mcols(annot.grange)$coding <- as.integer(ifelse(is.na(mcols(annot.grange)$CDS_position), FALSE, TRUE))
    synonymous <- c("synonymous_variant", "stop_retained_variant")
    non.synonymous <- c("stop_gained", "frameshift_variant", "stop_lost", "initiator_codon_variant", "inframe_insertion", "inframe_deletion", "missense_variant")
    syn.var <- in.csv.col(vec = mcols(annot.grange)$Consequence, search.vals = synonymous, delim.str = "&", match.func = any)
    non.syn.var <- in.csv.col(vec = mcols(annot.grange)$Consequence, search.vals = non.synonymous, delim.str = "&", match.func = any)
    both <- syn.var & non.syn.var
    syn.var[both] <- FALSE
    
    mcols(annot.grange)$Consequence <- gsub("&", ";", as.character(mcols(annot.grange)$Consequence))
    
    mcols(annot.grange)$Cons_cat <- "Other"
    mcols(annot.grange)$Cons_cat[syn.var] <- "Synonymous"
    mcols(annot.grange)$Cons_cat[non.syn.var] <- "NonSynonymous"
    
    mcols(annot.grange)$var_type <- mapply(function(ref.width, alt.width) {
        
        if (alt.width == ref.width && ref.width == 1) {
            return("SNV")
        }
        else if (alt.width == ref.width && ref.width > 1) {
            return("MNV")
        }
        else if (alt.width > ref.width) {
            return("INS")
        }
        else if (alt.width < ref.width) {
            return("DEL")
        }
        else {
            stop("ERROR: Unknown variant type")
        }
    }, width(mcols(annot.grange)$REF), width(mcols(annot.grange)$ALT))
    
    if('GMAF' %in% names(mcols(annot.grange)))
    {
        in.1kg <- is.na(mcols(annot.grange)$GMAF) == FALSE
        mcols(annot.grange)$in_1kg <- as.integer(in.1kg)
    }
    
    if ('Existing_variation' %in% names(mcols(annot.grange)) && all(is.na(mcols(annot.grange)$Existing_variation)) == F)
    {
        in.dbsnp <- grepl("rs\\d+", mcols(annot.grange)$Existing_variation)
        mcols(annot.grange)$in_dbsnp <- as.integer(in.dbsnp)
    }
    
    if ('CANONICAL' %in% names(mcols(annot.grange)))
    {
        stopifnot(levels(mcols(annot.grange)$CANONICAL) == "YES")
        mcols(annot.grange)$CANONICAL <- as.character(mcols(annot.grange)$CANONICAL)
        mcols(annot.grange)$CANONICAL[is.na(mcols(annot.grange)$CANONICAL)] <- "NO"
    }
   
    return(annot.grange)
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
