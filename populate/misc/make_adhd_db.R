#R-3.1.2

#Basic Genotype data

library(hwhelper)

BasicGenotypes <- function(geno.sample, geno.gene, sample.edge.name="HAS_GENOTYPE", gene.edge.name="RS_MAPPED_TO", node.name="snpID"){
    
    return(new("BasicGenotypes", geno.sample=geno.sample, sample.edge.name=sample.edge.name, geno.gene=geno.gene, gene.edge.name=gene.edge.name, node.name=node.name))
}

  #Needs to look like:
 #template.query='MATCH(subject:$SUBJECT$)-[d:DERIVED]-()-[r:HAS_EXPRESSION]-()-[:PS_MAPPED_TO]-(gene:EntrezID) WHERE d.type = "Affy_Expression" AND r.score > $PAR_NAME$ AND $$lower_coll_type$$.name IN {$$coll_type$$}
 #                       WITH $$lower_ret_type$$.name AS ret_type, COUNT(DISTINCT $$lower_coll_type$$.name) AS use_coll ORDER BY use_coll DESC RETURN ret_type, use_coll'))

setClass(Class="BasicGenotypes", representation=list(geno.sample="data.frame", geno.gene="data.frame"), contains="NeoData",
          prototype=list(
                        base.query='MATCH (n:$SUBJECT$)-[d:DERIVED]-(samp)-[r:HAS_GENOTYPE]-(var)-[r2:RS_MAPPED_TO]-(gene:EntrezID{name:{GENE}})-[:REFFERED_TO]-(symb) WHERE d.type = "Genotype" AND ANY(x IN [n.name, samp.name] WHERE x IN {SAMPLE})
                            RETURN DISTINCT var.name AS rsID, gene.name AS Gene, symb.name AS Symbol, var.genotype AS Genotype,
                            0 AS query_ind, 1 AS gene_ind, var.name + "_" + gene.name AS row_id, n.name AS Sample',
                        
                        template.query='MATCH (subject:$SUBJECT$)-[d:DERIVED]-()-[r:HAS_GENOTYPE]-(var)-[r2:RS_MAPPED_TO]-(gene:EntrezID) WHERE d.type = "Genotype" AND $$lower_coll_type$$.name IN {$$coll_type$$}
                        WITH $$lower_ret_type$$.name AS ret_type, COUNT(DISTINCT $$lower_coll_type$$.name) AS use_coll ORDER BY use_coll DESC RETURN ret_type, use_coll',
                        
                        sample.edge.name="HAS_GENOTYPE",
                        gene.edge.name="RS_MAPPED_TO",
                        node.name="snpID"
                        ))

GWASResult <- function(geno.gene, gene.edge.name="RS_MAPPED_TO", node.name="snpID"){
    return(new("GWASResult", geno.gene=geno.gene, gene.edge.name=gene.edge.name, node.name=node.name))
}

setClass(Class="GWASResult", contains=c("BasicGenotypes", "HwHit"),
         prototype=list(
                        
                        geno.sample=data.frame(),
                        
                        base.query='MATCH (n:$SUBJECT$)-[d:DERIVED]-(samp) WHERE ANY(x IN [n.name, samp.name] WHERE x IN {SAMPLE}) WITH n MATCH (snp)-[r:RS_MAPPED_TO]-(gene:EntrezID{name:{GENE}})
                        WHERE HAS(r.score) AND r.score < $PAR_NAME$ WITH n, gene, MIN(r.score) AS gwas_pvalue
                        RETURN gene.name AS gene, n.name AS sample, "$DATA_NAME$" AS var, gwas_pvalue AS score, gwas_pvalue < $PAR_NAME$ AS is_hit;',
                        
                        template.query='MATCH (snp)-[r:RS_MAPPED_TO]-(gene) WHERE HAS(r.score) AND r.score < $PAR_NAME$ WITH gene MATCH (subject:MergeID)-[d:DERIVED]-(s) WHERE
                        $$lower_coll_type$$.name IN {$$coll_type$$} WITH $$lower_ret_type$$.name AS ret_type, COUNT(DISTINCT $$lower_coll_type$$.name) AS use_coll
                        ORDER BY use_coll DESC RETURN ret_type, use_coll',
                        
                        default=.0001, direction="<", range=c(0, 1), display_name="GWAS Hit Threshold"))

setMethod("sampleNames", signature("BasicGenotypes"), function(object){
    return(unique(object@geno.sample$Sample))
})

setClass(Class="MethylResult", representation=list(methyl.gene="data.frame"),contains=c("NeoData", "HwHit"),
         prototype=list(
                        
                        sample.edge.name="HAS_EXPRESSION",
                        gene.edge.name="METHYL_MAPPED_TO",
                        node.name="methylID",
                        
                        geno.sample=data.frame(),
                        
                        base.query='MATCH (n:$SUBJECT$)-[d:DERIVED]-(samp) WHERE ANY(x IN [n.name, samp.name] WHERE x IN {SAMPLE}) WITH n MATCH (meth)-[r:METHYL_MAPPED_TO]-(gene:EntrezID{name:{GENE}})
                        WHERE HAS(r.score) AND r.score < $PAR_NAME$ WITH n, gene, MIN(r.score) AS methyl_pvalue
                        RETURN gene.name AS gene, n.name AS sample, "$DATA_NAME$" AS var, methyl_pvalue AS score, methyl_pvalue < $PAR_NAME$ AS is_hit;',
                        
                        template.query='MATCH (meth)-[r:METHYL_MAPPED_TO]-(gene:EntrezID) WHERE HAS(r.score) AND r.score < $PAR_NAME$ WITH gene MATCH (subject:MergeID)-[d:DERIVED]-(s) WHERE
                        $$lower_coll_type$$.name IN {$$coll_type$$} WITH $$lower_ret_type$$.name AS ret_type, COUNT(DISTINCT $$lower_coll_type$$.name) AS use_coll
                        ORDER BY use_coll DESC RETURN ret_type, use_coll',
                        
                        default=.001, direction="<", range=c(0, 1), display_name="Methylation Hit Threshold"))

setMethod("fromSample", signature("MethylResult"), function(obj, neo.path=NULL){})

setMethod("toGene", signature("MethylResult"), function(obj, neo.path=NULL, gene.model="entrez"){
    
    load.neo4j(.data=obj@methyl.gene, edge.name=geneEdge(obj), commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")    
})


setMethod("fromSample", signature("BasicGenotypes"), function(obj, neo.path=NULL){
    
    #also note there that things like presence in dbSNP or COSMIC etc could be used as a property in the Variation node--should add in Variant_Type here...
    
    #This is because we want to simply add the GWASResult class results in without having to add them per-sample
    if (nrow(obj@geno.sample) > 0)
    {
        load.neo4j(.data=obj@geno.sample, edge.name=sampleEdge(obj), commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")
    }
})


setMethod("toGene", signature("BasicGenotypes"), function(obj, neo.path=NULL, gene.model="entrez"){
    
    load.neo4j(.data=obj@geno.gene, edge.name=geneEdge(obj), commit.size=10000L, neo.path=neo.path, dry.run=F, array.delim="&")    
})

main <- function()
{
    
    #as these are hg18, need to lift them over
    pcg.res <- read.delim("pgc.adhd.full.2012-10.txt", sep="", header=T, stringsAsFactors=F)
    
    geno.res <- read.delim("final.coded.txt", sep="\t", header=T, stringsAsFactors=F)
    
    library(rtracklayer)
    
    #get locations etc for the genotyps
    
    mySession <- browserSession()
    genome(mySession) <- "hg19"
    
    geno.query <- ucscTableQuery(mySession, track="All SNPs(138)", table="snp138", names=geno.res$Name)
    
    geno.annot <- getTable(geno.query)
    geno.grange <- track(geno.query)
    
    all(geno.annot$name %in% geno.res$Name)#T
    
    nrow(geno.res) == nrow(geno.annot)#T
    
    library(GenomicFeatures)
    
    ref.gene <- loadDb("refGene_hg19.db")
    
    ens.genes <- genes(ref.gene)
    
    closest.genes <- nearest(geno.grange, ens.genes, select="all")
    
    distances <- distance(geno.grange[queryHits(closest.genes)], ens.genes[subjectHits(closest.genes)])
    
    geno.gene <- data.frame(snpID=geno.grange$name[queryHits(closest.genes)], entrezID=ens.genes$gene_id[subjectHits(closest.genes)], distance=distances, stringsAsFactors=F)
    
    library(reshape2)
    
    geno.sample <- melt(geno.res, id.vars="Name", measure.vars=colnames(geno.res)[grep("sample", colnames(geno.res))])
    geno.sample$Name <- as.character(geno.sample$Name)
    geno.sample$variable <- as.character(geno.sample$variable)
    geno.sample$value <- as.integer(geno.sample$value)
    names(geno.sample) <- c("snpID", "sample", "genotype")
    geno.sample <- geno.sample[,c("sample", "snpID", "genotype")]
    
   
    
    #The hits themselves may be able to be represented simply as SnpID->EntrezID simply attached to the samples
    
    chain <- import.chain("/Users/bottomly/Desktop/resources/sequences/hg18ToHg19.over.chain")
    
    pgc.grange <- with(pcg.res, GRanges(seqnames=hg18chr, ranges=IRanges(start=bp, width=1), strand="*", snpID=snpid, score=pval))
    seqlevelsStyle(pgc.grange) <- "UCSC"
    
    pgc.lo <- liftOver(pgc.grange, chain)
    
    pgc.hg19 <- unlist(pgc.lo)
    
    closest.genes <- nearest(pgc.hg19, ens.genes, select="all")
    
    distances <- distance(pgc.hg19[queryHits(closest.genes)], ens.genes[subjectHits(closest.genes)])
    
    pcg.gene <- data.frame(snpID=pgc.hg19$snpID[queryHits(closest.genes)], entrezID=ens.genes$gene_id[subjectHits(closest.genes)], distance=distances, score=pgc.hg19$score[queryHits(closest.genes)], stringsAsFactors=F)
    
    save(pcg.gene, geno.sample, geno.gene, file="temp_gwas_snps.RData")
    
    load("temp_gwas_snps.RData")
    
    met.geno <- GWASResult(geno.gene=pcg.gene)
    
    geno.sample <- geno.sample[is.na(geno.sample$genotype) == F & geno.sample$genotype != 0,]
    
    basic.geno <- BasicGenotypes(geno.sample=geno.sample, geno.gene=geno.gene)
    
    subj.info <- read.delim("passed_Dt_withSaliva_forHW2.txt", sep="\t", header=T, stringsAsFactors=F, skip=2)
    
    sub.subj <- subj.info[,c("Mergeid", "Sex.x", "ADHD_UPDATE_MD_Recoded", "Access_Adhd_Status.x", "Y1_Stim_Med")]
    names(sub.subj) <- c("mergeID", "Gender", "SubPheno", "Status", "Meds")
    
    sub.subj$Gender[is.na(sub.subj$Gender)] <- "UnkGender"
    sub.subj$Gender[sub.subj$Gender == 1] <- "Male"
    sub.subj$Gender[sub.subj$Gender == 2] <- "Female"
    
    #not sure about SubPheno--leave as integers
    
    sub.subj$SubPheno[is.na(sub.subj$SubPheno)==F] <- paste0("SubPheno", sub.subj$SubPheno[is.na(sub.subj$SubPheno)==F])
    sub.subj$SubPheno[is.na(sub.subj$SubPheno)] <- "SubPhenoUnk"
    
    sub.subj$Status <- ifelse(sub.subj$Status == 1, "Case", "Control")
    
    sub.subj$Meds <- ifelse(is.na(sub.subj$Meds), "NoMeds", "Meds")
    
    subj <- Subject(subject.info=sub.subj,type.col="Status")
    
    subj.geno.map <- read.delim("IDmapping.txt", sep="\t", header=T, stringsAsFactors=F)
    oragene <- read.delim("Oragene_Log_Abbr_3_28_2015_BW.txt", sep="\t", stringsAsFactors=F)
    
    subj.geno.map$fix.date <- sapply(strsplit(subj.geno.map$Date, "/"), function(x) paste(x[1], x[2], substr(x[3], 3,4), sep="/"))
    subj.geno.map$fix.samp <- substr(subj.geno.map$ID, 1, 5)
    
    temp <- merge(oragene, subj.geno.map, by.x=c("FamID", "DNA_Date_Collect"), by.y=c("fix.samp", "fix.date"), all.x=F, all.y=T, sort=F)
    
    test.temp <- temp[is.na(temp$Child_Parent)==F & temp$Child_Parent == "Child",]
    
    sum.temp <- aggregate(Mergeid~Date+FamID+Sample, function(x) paste(unique(x), collapse=","), data=test.temp)
    sum.temp[sum.temp$Date == "5/29/2010" & sum.temp$FamID == "10418","Mergeid"] <- "104181"
    
    add.dta <- data.frame(Date=c("11/6/2010", "12/19/2009", "4/23/2010"), FamID=c("10230", "10235", "10305"), Mergeid=c("102301", "102351", "103051"), stringsAsFactors=F)
    
    add.dta$Sample <- NA
    
    for(i in 1:length(add.dta$FamID))
    {
        add.dta[i,"Sample"] <- subj.geno.map[subj.geno.map$ID == add.dta[i,"FamID"],"Sample"]
    }
    
    sum.temp.fin <- rbind(sum.temp, add.dta[,names(sum.temp)])
    
    sum.temp.fin <- rbind(sum.temp.fin, data.frame(Date=NA, FamID=10394, Sample=74, Mergeid=103941, stringsAsFactors=F))
    
    setdiff(subj.geno.map$Sample, sum.temp.fin$Sample)#integer(0)
    
    #still need to add in this part
    
    sum.temp.fin <- sum.temp.fin[,c("Mergeid", "Sample")]
    names(sum.temp.fin) <- c("mergeID","sample")
    sum.temp.fin$sample <- paste0("sample", sum.temp.fin$sample)
    sum.temp.fin$type <- "Genotype"
    
    addSamples(subj) <- sum.temp.fin
    
    #also add in the methylation
    
    library(GenomicFeatures)
    library(rtracklayer)
    
    meth.1 <- read.delim("Ill_final_Dec.txt", sep="\t", header=T, stringsAsFactors=F)
    
    meth.1 <- meth.1[is.na(meth.1$MAPINFO)==F,]
    
    sub.meth <- meth.1[,c("gc.p.value", "GENOME_BUILD", "CHR", "MAPINFO", "NAME")]
    
    #lift over the b36 probes
    
    b36.probes <- with(sub.meth[sub.meth$GENOME_BUILD == 36,], GRanges(seqnames=CHR, ranges=IRanges(start=MAPINFO, width=1), strand="*", p.value=gc.p.value, methylID=NAME))
    seqlevelsStyle(b36.probes) <- "UCSC"
    
    chain <- import.chain("/Users/bottomly/Desktop/resources/sequences/hg18ToHg19.over.chain")
    
    b36.to37.probes <- liftOver(b36.probes, chain)
    b36.to37.probes.unl <- unlist(b36.to37.probes)
    
    #add to the b37 probes
    b37.probes <- with(sub.meth[sub.meth$GENOME_BUILD == 37,], GRanges(seqnames=CHR, ranges=IRanges(start=MAPINFO, width=1), strand="*", p.value=gc.p.value, methylID=NAME))
    seqlevelsStyle(b37.probes) <- "UCSC"
    
    all.probes <- append(b37.probes, b36.to37.probes.unl)
    
    ref.gene <- loadDb("refGene_hg19.db")
    
    ens.genes <- genes(ref.gene)
    
    closest.genes <- nearest(all.probes, ens.genes, select="all")
    
    distances <- distance(all.probes[queryHits(closest.genes)], ens.genes[subjectHits(closest.genes)])
    
    meth.genes <- data.frame(methylID=all.probes$methylID[queryHits(closest.genes)], entrezID=ens.genes$gene_id[subjectHits(closest.genes)], distance=distances, score=all.probes$p.value[queryHits(closest.genes)], stringsAsFactors=F)
    
    #add to the methylation class
    
    save(meth.genes, file="temp_meth_results.RData")
    
    load("temp_meth_results.RData")
    
    methyl.dta <- new("MethylResult", methyl.gene=meth.genes)
    
    #Still need to do...
    hw2.conf <- makeHW2Config(subject=subj, gene.model="entrez",
        data.types=list(seeds="GWASResult", target="Genotype"),
        Methylation=methyl.dta, Genotype=basic.geno, GWASResult=met.geno)
    
    populate(hw2.conf)
    
    ccle.graph <- compute.graph.structure()
    
    save(ccle.graph, file="neo_graph.RData")
    
    #from the adhd folder
    #library(igraph)
    #load("neo_graph.RData")
    #plot(ccle.graph)
    
    configure(hw2.conf)
    
}

get.gene.info <- function(){
    
    refgene.hg19 <- makeTranscriptDbFromUCSC(genome="hg19", tablename="refGene")
    
    saveDb(refgene.hg19, "refGene_hg19.db")
}