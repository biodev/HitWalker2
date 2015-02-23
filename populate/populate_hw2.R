library(Biobase)
<<<<<<< HEAD
library(reshape2)
=======
>>>>>>> 2290ff2320a5a94b4b05caf9301379b994fb00ab

#importing pathway commons

#download .gmt file from http://www.pathwaycommons.org/pc2/downloads.html
#here it is in /Users/bottomly/Desktop/hitwalker2_paper/Pathway Commons.4.All.GSEA.gmt.gz

#to incorporate into hitwalker use the org.Hs.eg.db database to relate uniprot to entrez
#both the private and public versions will have reference to Entrez IDs so this will be useful...

setGeneric("fromSample", def=function(obj,...) standardGeneric("fromSample"))
setGeneric("toGene", def=function(obj,...) standardGeneric("toGene"))
<<<<<<< HEAD
setGeneric("getMatrix", def=function(obj,...) standardGeneric("getMatrix"))
=======
>>>>>>> 2290ff2320a5a94b4b05caf9301379b994fb00ab

setClass(Class="HW2Config", representation=list(data.list="list", data.types="list", gene.models="character", neo.path="character"))

#expression utils, affy for now...

setMethod("fromSample", signature("ExpressionSet"), function(obj, neo.path, to.node="probeSet", edge.name="HAS_EXPRESSION"){
    
<<<<<<< HEAD
    use.exprs <- exprs(obj)
    
    melt.use.exprs <- melt(use.exprs)
    names(melt.use.exprs) <- c("probeset", "sample", "score")
    
    melt.use.exprs$sample <- as.character(melt.use.exprs$sample)
    melt.use.exprs$probeset <- as.character(melt.use.exprs$probeset)
    
    #now load the sample -> probe mappings
    
    melt.use.exprs <-  melt.use.exprs[,c("sample", "probeset", "score")]
    names(melt.use.exprs)[2] <- c("probeSet")
    
    load.neo4j(.data=melt.use.exprs, edge.name=edge.name, commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&", unique.rels=F)
    
})

setMethod("toGene", signature("ExpressionSet"), function(obj, neo.path, from.node="probeSet", gene.model=c("entrez", "ensembl"), annotation.package=NULL, edge.name="PS_MAPPED_TO"){
    
    if (is.null(annotation.package) || all(is.na(annotation.package)))
    {
        stop("ERROR: annotation.package needs to be specified")
    }else{
        require(annotation.package, character.only=T)
    }
    
    #Then create a mapping file from probeset to gene
    
    gene.type <- switch(gene.model, entrez="ENTREZID", ensembl="ENSEMBL")
    
    probe.to.gene <- select(eval(parse(text=annotation.package)), keys=featureNames(obj), column=gene.type, keytypes="PROBEID")
    
    #Discard for now those that do not map to either type of genes.
    
    probe.to.gene <- probe.to.gene[is.na(probe.to.gene[,gene.type]) == F,]
    
    #probesets that map to multiple genes are perhaps not that reliable either, so discard them as well for now...
    
    dup.probes <- probe.to.gene$PROBEID[duplicated(probe.to.gene$PROBEID)]
    
    probe.to.gene <- probe.to.gene[probe.to.gene$PROBEID %in% dup.probes == F,]
    names(probe.to.gene) <- c("probeSet", switch(gene.model, entrez="entrezID", ensembl="ensembl"))
    
    #load the probe->gene mappings
    
    load.neo4j(.data=probe.to.gene, edge.name=edge.name, commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")
=======
})

setMethod("toGene", signature("ExpressionSet"), function(obj, neo.path, from.node="probeSet", gene.model=c("entrez", "ensembl"), annotation.package="", edge.name="PS_MAPPED_TO"){
    
>>>>>>> 2290ff2320a5a94b4b05caf9301379b994fb00ab
})

#MAF class utils

setClass(Class="CCLEMaf", representation=list(maf="data.frame"))

readMAF.ccle <- function(file.name)
{
    use.maf <- read.delim("file.name", sep="\t", stringsAsFactors=F)
    
    keep.maf <- use.maf[,c("Entrez_Gene_Id", "Genome_Change", "Variant_Classification", "Annotation_Transcript", "Transcript_Strand", "cDNA_Change", "Codon_Change", "Protein_Change",
                           "Tumor_Sample_Barcode", "Genome_Change", "Center", "Sequencer", "Alternative_allele_reads_count", "Reference_allele_reads_count", "dbSNP_RS", "dbSNP_Val_Status")]

    return(new("CCLEMaf", maf=keep.maf))
}

setGeneric("maf", def=function(obj,...) standardGeneric("maf"))
setMethod("maf", signature("CCLEMaf"), function(obj){
    return(obj@maf)
})


setMethod("fromSample", signature("CCLEMaf"), function(obj, neo.path, to.node="variation", edge.name="HAS_DNASEQ"){
    #first sample -> variant
    #the name here will be derived from the Genome_Change column as that provides potentially enough information to uniquely id a variant (indels might still be tricky...)
    #will keep missing values as "" for now
    
    sample.maf <- maf(obj)
    
    samp.maf <- sample.maf[,c("Tumor_Sample_Barcode", "Genome_Change", "Center", "Sequencer", "Alternative_allele_reads_count", "Reference_allele_reads_count", "dbSNP_RS", "dbSNP_Val_Status")]
    names(samp.maf) <- c("sample", to.node, "center", "sequencer", "alt_counts", "ref_counts", "variation.dbsnp", "variation.dbsnp_val_status")
    samp.ccle$variation.dbsnp <- gsub(";", "&", samp.ccle$variation.dbsnp)
    samp.ccle$variation.dbsnp_val_status <- gsub(";", "&", samp.ccle$variation.dbsnp_val_status)
    
    #also note there that things like presence in dbSNP or COSMIC etc could be used as a property in the Variation node--should add in Variant_Type here...
    
    load.neo4j(.data=samp.ccle, edge.name=edge.name, commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")
})

setMethod("toGene", signature("CCLEMaf"), function(obj, neo.path, from.node="variation", gene.model="entrez", edge.name="IMPACTS"){
     #then add in the Variation->EntrezID relationships
    
    if (gene.model != "entrez")
    {
        stop("ERROR: Only gene.model = 'entrez' is supported for MAF files")
    }
    
    cur.maf <- maf(obj)
    
    #only keep one row for each gene/variant
    var.gene.dta <- cur.maf[!duplicated(non.na.ccle[,c("Entrez_Gene_Id", "Genome_Change")]),]
    
    var.gene.dta <- var.gene.dta[,c("Genome_Change", "Entrez_Gene_Id", "Variant_Classification", "Annotation_Transcript", "Transcript_Strand", "cDNA_Change", "Codon_Change", "Protein_Change")]
    
    names(var.gene.dta) <- c(from.node, "entrezID","variant_classification", "transcript", "transcript_strand", "cdna", "codon", "protein")
    
    load.neo4j(.data=var.gene.dta, edge.name=edge.name, commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")  
    
})

#Drug class utils

setClass(Class="DrugMatrix", representation=list(matrix="matrix", mapping="data.frame"))

setMethod("getMatrix", signature("DrugMatrix"), function(obj)
          {
                return(obj@matrix)
          })

setMethod("fromSample", signature("DrugMatrix"), function(obj, neo.path, to.node="variation", edge.name="HAS_DNASEQ"){
    
    drug.mat <- getMatrix(obj)
    
    drug.dta <- melt(drug.mat)
    
    names(drug.dta) <- c("drug", "sample", "score")
    
    drug.dta$drug <- as.character(drug.dta$drug)
    drug.dta$sample <- as.character(drug.dta$sample)
    
    #remove any na scores
    
    stopifnot(sum(is.na(drug.dta$value)) == sum(is.na(drug.dta)))
    
    drug.dta <- drug.dta[complete.cases(drug.dta),]
    
    load.neo4j(.data=drug.dta, edge.name=edge.name, commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")
    
})

setMethod("toGene", signature("DrugMatrix") function(obj){
    
})


make.hw2.database <- function(obj, data.types, neo.path, gene.model=c("entrez", "ensembl"))
{
    
}

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

make.property.str.node <- function(prop.names, dta, pos, array.delim)
    {
        base.str <- paste0("{name:line[",pos,"]")
        
        for(i in prop.names)
        {
            use.name <- sapply(strsplit(i, "\\."), function(x) x[-1])
            use.pos <- which(names(dta) == i)-1
            
            stopifnot(length(use.pos) == 1)
            
            if (is.integer(dta[,i]))
            {
                base.str <- paste0(base.str, ",",use.name,":toInt(line[",use.pos,"])")
            }else if (is.numeric(dta[,i]))
            {
                base.str <- paste0(base.str, ",",use.name,":toFloat(line[",use.pos,"])")
            }else{
                
                if (is.null(array.delim) == F && is.character(array.delim) == T && length(array.delim) == 1)
                {
                    if (any(grepl(array.delim, dta[,i])))
                    {
                        base.str <- paste0(base.str, ",",use.name,":split(line[",use.pos,"], '",array.delim,"')")
                    }else{
                        base.str <- paste0(base.str, ",",use.name,":line[",use.pos,"]")
                    }
                }else{
                    base.str <- paste0(base.str, ",",use.name,":line[",use.pos,"]")
                }
            }
        }
        
        base.str <- paste0(base.str, "}")
        
        return(base.str)
    }

get.or.create.constraints <- function(neo.path, node.names)
{
    #before executing, things will go A LOT faster if constraints/indexes are available, check first and if not, then create them
    
    found.schema <- system(paste(file.path(neo.path, "bin", "neo4j-shell"), "-c schema"), intern=T)
    
    constraint.section <- grep("Constraints", found.schema)
    
    if (length(constraint.section) > 0)
    {
        found.constraints <- found.schema[(constraint.section+1):length(found.schema)]
    
        present.constraints <- sapply(node.names, function(x)
               {
                    any(grepl(paste0("\\s+ON\\s+\\(",tolower(x),":",x,"\\) ASSERT\\s+", tolower(x), ".name IS UNIQUE"), found.constraints))
               })
    }else{
        #no constraints were found so set both notes to present.constraints=F
        present.constraints <- rep(F,2)
        names(present.constraints) <- node.names
    }
    
    if (all(present.constraints))
    {
        message("Found all necessary CONSTRAINTS, continuing...")
        missing.constraints <- character(0)
    }else{
        missing.constraints <- names(present.constraints)[present.constraints == F]
        message(paste("Missing CONSTRAINT(s) for:", paste(missing.constraints, collapse=","), "will build them first and continue with loading..."))
    }
    
    if (length(missing.constraints) > 0)
    {
        constr.str <- paste0(file.path(neo.path, "bin", "neo4j-shell"), " -c 'CREATE CONSTRAINT ON (",tolower(missing.constraints),":",missing.constraints,") ASSERT ",tolower(missing.constraints),".name IS UNIQUE;'")
    }else{
        constr.str <- character(0)
    }
    
    return (constr.str)
}

load.neo4j <- function(.data, edge.name=NULL, commit.size=1000, neo.path="/Users/bottomly/Desktop/resources/programs/neo4j-community-2.1.2", dry.run=F, array.delim="&", unique.rels=T, merge.from=T, merge.to=T)
{
    if (missing(edge.name) || is.null(edge.name))
    {
        load.neo4j.node(.data=.data, commit.size=commit.size, neo.path=neo.path, dry.run=dry.run, array.delim=array.delim)
        
    }else{
        load.neo4j.edge(.data=.data, edge.name=edge.name, commit.size=commit.size, neo.path=neo.path, dry.run=dry.run, array.delim=array.delim, unique.rels=unique.rels, merge.from=merge.from, merge.to=merge.to)
    }
}

load.neo4j.node <- function(.data, commit.size=1000, neo.path="/Users/bottomly/Desktop/resources/programs/neo4j-community-2.1.2", dry.run=F, array.delim="&")
{
    if (any(is.na(.data)))
    {
        stop("ERROR: .data currently cannot have NAs, please remove and try again...")
    }
    
    base.node.names <- names(.data)[1]
    
    node.names <- capwords(base.node.names)
    
    node.props.dta <- lapply(base.node.names, function(x)
                         {
                            return(names(.data)[grep(paste0(x,"\\."), names(.data))])
                         })
    
    use.temp <- tempfile()
    
    write.table(.data, sep="\t", file=use.temp, col.names=F, row.names=F, quote=F)
    
    cypher.temp <- tempfile()
    
    constr.str <- get.or.create.constraints(neo.path, node.names)
    
    cypher.stats <- c(
        paste0("USING PERIODIC COMMIT ",commit.size),
        paste0("LOAD CSV FROM 'file://", use.temp,"' AS line FIELDTERMINATOR '\t'"),
        paste0("MERGE (n:",node.names[1], make.property.str.node(node.props.dta[[1]], .data, 0, array.delim),");")
    )
    
    writeLines(cypher.stats, con=cypher.temp)
    
    if (dry.run == F)
    {
        for (i in constr.str)
        {
            system(i)
        }
        
        message("Loading into Neo4j...")
        system(paste(file.path(neo.path, "bin", "neo4j-shell"), "-file", cypher.temp))
    }else{
        
        for (i in constr.str)
        {
            message(i)
        }
        
        message(paste(file.path(neo.path, "bin", "neo4j-shell"), "-file", cypher.temp))
    }
}

#data.frame should be in the form:
#data.frame(node1, node2, node[12].property, edge property(no x'.'y just name))
load.neo4j.edge <- function(.data, edge.name=NULL, commit.size=1000, neo.path="/Users/bottomly/Desktop/resources/programs/neo4j-community-2.1.2", dry.run=F, array.delim="&", unique.rels=T, merge.from=T, merge.to=T)
{
    message("Preprocessing...")
    
    if (missing(edge.name) || is.null(edge.name) || is.na(edge.name))
    {
        stop("ERROR: Need to supply an edge name")
    }else if (is.character(edge.name) == F || length(edge.name) != 1)
    {
        stop("ERROR: edge.name needs to be a single string value")
    }
    
    if (ncol(.data) < 2 || is.null(names(.data)))
    {
        stop("ERROR: .data needs to have at least two columns and be named")
    }
    
    if (any(is.na(.data)))
    {
        stop("ERROR: .data currently cannot have NAs, please remove and try again...")
    }
    
    
    .make.property.str.edge <- function(prop.names, dta, array.delim)
    {
        if (length(prop.names) == 0)
        {
            return("")
        }
        
        base.str <- "{"
        
        for(i in prop.names)
        {
            if (which(prop.names == i) > 1)
            {
                base.str <- paste0(base.str, ",")
            }
            
            use.pos <- which(names(dta) == i)-1
            
            stopifnot(length(use.pos) == 1)
            
            if (is.integer(dta[,i]))
            {
                base.str <- paste0(base.str,i,":toInt(line[",use.pos,"])")
            }else if (is.numeric(dta[,i]))
            {
                base.str <- paste0(base.str,i,":toFloat(line[",use.pos,"])")
            }else{
                
                if (is.null(array.delim) == F && is.character(array.delim) == T && length(array.delim) == 1)
                {
                    if (any(grepl(array.delim, dta[,i])))
                    {   
                        base.str <- paste0(base.str,i,":split(line[",use.pos,"], '", array.delim,"')")
                    }else{
                        base.str <- paste0(base.str,i,":line[",use.pos,"]")
                    }
                }else{
                
                    base.str <- paste0(base.str,i,":line[",use.pos,"]")
                
                }
            }
        }
        
        base.str <- paste0(base.str, "}")
        
        return(base.str)
    }
    
    base.node.names <- names(.data)[1:2]
    
    node.names <- capwords(base.node.names)
    
    node.props.dta <- lapply(base.node.names, function(x)
                         {
                            return(names(.data)[grep(paste0(x,"\\."), names(.data))])
                         })
    
    edge.props <- setdiff(names(.data), c(base.node.names, unlist(node.props.dta)))
    
    use.temp <- tempfile()
    
    write.table(.data, sep="\t", file=use.temp, col.names=F, row.names=F, quote=F)
    
    cypher.temp <- tempfile()
    
    constr.str <- get.or.create.constraints(neo.path, node.names)
    
    rel.unique.str <- ifelse(unique.rels, "UNIQUE", "")
    
    if (merge.from){
        from_str = "MERGE"
    }else{
        from_str = "MATCH"
    }
    
    if (merge.to){
        to_str = "MERGE"
    }else{
        to_str = "MATCH"
    }
    
    cypher.stats <- c(
        paste0("USING PERIODIC COMMIT ",commit.size),
        paste0("LOAD CSV FROM 'file://", use.temp,"' AS line FIELDTERMINATOR '\t'"),
        paste0(from_str, " (n:",node.names[1], make.property.str.node(node.props.dta[[1]], .data, 0, array.delim),")"),
        paste0(to_str, " (m:",node.names[2],make.property.str.node(node.props.dta[[2]], .data, 1, array.delim),")"),
        paste0("CREATE ",rel.unique.str," (n)-[:",edge.name,.make.property.str.edge(edge.props, .data, array.delim),"]->(m);")
    )
    
    writeLines(cypher.stats, con=cypher.temp)
    
    if (dry.run == F)
    {
        for (i in constr.str)
        {
            system(i)
        }
        
        message("Loading into Neo4j...")
        system(paste(file.path(neo.path, "bin", "neo4j-shell"), "-file", cypher.temp))
    }else{
        
        for (i in constr.str)
        {
            message(i)
        }
        
        message(paste(file.path(neo.path, "bin", "neo4j-shell"), "-file", cypher.temp))
    }
}



sample.neo4j.graph <- function(graph_struct="/var/www/hitwalker_2_inst/graph_struct.json", start.node.label=NULL, num.start.nodes.sample=10, base.query="MATCH (n)-[r]-(m)", max.rels=100, dump.file="test.ndp", neo.path="/Users/bottomly/Desktop/resources/programs/neo4j-community-2.1.6")
{
    #read.graph_struct
    
    g.struct <- read.graph_struct(graph_struct)
    
    if (missing( start.node.label) || is.null( start.node.label) || all(is.na( start.node.label)))
    {
        V(g.struct)$name
    }
    
    all.edges <- get.edgelist(g.struct)
    
    previous.vars <- character(0)
    
    entire.statement <- ""
    
    for(i in 1:nrow(all.edges))
    {
        cur.e.label <- E(g.struct)[V(g.struct)[name == all.edges[i,1]] %->% V(g.struct)[name == all.edges[i,2]]]$label
        
        previous.vars <- union(previous.vars, c(tolower(all.edges[i,1]), tolower(all.edges[i,2]), tolower(cur.e.label)))
        
        entire.statement <- paste(entire.statement, paste0('MATCH (',tolower(all.edges[i,1]),':',all.edges[i,1],')-[',tolower(cur.e.label),':',cur.e.label,']-(',tolower(all.edges[i,2]),':',all.edges[i,2],') WITH ', paste(previous.vars, collapse=","), ' LIMIT ', max.rels))
    }
    
    entire.statement <- paste(entire.statement, "RETURN", paste(previous.vars, collapse=","))
    
    tmp.file <- tempfile()
    
    writeLines(paste0("DUMP", entire.statement, ";"), con=tmp.file)
    
    system(paste(file.path(neo.path, "bin", "neo4j-shell"), "-file", tmp.file,">", dump.file))
}

read.graph_struct  <- function(graph_struct="/var/www/hitwalker_2_inst/graph_struct.json")
{
    require(rjson)
    require(igraph)
    
    json.obj <- fromJSON(file=graph_struct)
    
    #warnings are about the rownames...
    json.dta <- suppressWarnings(data.frame(do.call(rbind, lapply(names(json.obj), function(x) {
        if (length(json.obj[[x]]) > 0)
        {
            temp.dta <- cbind(x, names(json.obj[[x]]), unlist(json.obj[[x]]))
        }else{
            return(NULL)
        }
        
    })), stringsAsFactors=F))
    names(json.dta) <- c("from", "label", "to")
    json.dta <- json.dta[!duplicated(json.dta$label),]
    json.graph <- graph.data.frame(json.dta[,c("from", "to", "label")], directed=F)
    return(json.graph)
    #plot(json.graph)
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

