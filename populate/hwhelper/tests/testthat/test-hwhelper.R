snv.file <- "/Volumes/Macintosh_HD_4/tests/broad_examp_vcf/CEUTrio.HiSeq.WGS.b37.NA12878.refseq.filt.vcf"
indel.file <- "/Volumes/Macintosh_HD_4/tests/broad_examp_vcf/test_indel.vcf"

mutect.file <- "/Volumes/Macintosh_HD_4/tests/db_test/Pair_15-00208_AML_15-00209_Skin.mutect.c.0.annot.filt.vcf"
varscan.file <- "/Volumes/Macintosh_HD_4/tests/db_test/Pair_15-00208_AML_15-00209_Skin.varscan2.t.100.n.100.snp.annot.fixed.filt.vcf"


if(! (file.exists(snv.file) && file.exists(indel.file) && file.exists(mutect.file) && file.exists(varscan.file)) ){
  stop("tests cannot be carried out")
}

test.list <- lapply(c(snv.file, indel.file), function(x){
  
  table <- make.vcf.table(x, info.import=c("FS", "MQ0", "MQ", "QD", "SB", "CSQ"),
                 readcount.import="AD", gds.out="variant.gds", keep.gds=F)
  
  obj <- readVcf(x, genome="test")
  return(list(table=table, obj=obj))  
})

somatic.list <- list()
  
mut.table <- make.vcf.table(mutect.file, 
                                    info.import = c("MQ0", "SOMATIC", "CSQ"), keep.samples="15-00208_AML",
                                    ignore.genotype=T, readcount.import= c("AD"))
mut.obj <- readVcf(mutect.file, genome="test")

somatic.list <- append(somatic.list, list(list(table=mut.table, obj=mut.obj)))

var.table <-  make.vcf.table(varscan.file, 
                             info.import = c("SS", "SOMATIC", "CSQ"), keep.samples="TUMOR",
                             ignore.genotype=F, readcount.import= c("RD", "AD"))
  
var.obj <- readVcf(varscan.file, genome="test")

somatic.list <- append(somatic.list, list(list(table=var.table, obj=var.obj)))

##Tests specific to the make.vcf.table function

test_that("all the specified info.import columns minus CSQ should be present",{
  
  for(i in seq_along(test.list)){
    expect_true(all(c("FS", "MQ0", "MQ", "QD", "SB") %in% names(test.list[[i]]$table)))
  }
   
})

test_that("the csq data should be the same as that from ensemblVEP",{
  
  for(i in seq_along(test.list)){
    
    test.vcf <- test.list[[i]]$obj
    test <- test.list[[i]]$table
    
    csq.grange <- parseCSQToGRanges(test.vcf, VCFRowID=rownames(test.vcf))
    names(csq.grange) <- NULL
    csqs <- as.data.frame(csq.grange)
    names(csqs)[names(csqs) == "VCFRowID"] <- "variant_id"
    csqs$REF <- as.character(ref(test.vcf))[csqs$variant_id]
    csqs <- csqs[,names(csqs) != "Allele"]
    csqs$ALT <- mapply(function(x,y){
      
      return(x[y])
    }, as(alt(test.vcf), "CharacterList")[csqs$variant_id], csqs$ALLELE_NUM)
    
    csqs <- hwhelper:::factors.to.chars(csqs)
    
    test.merge <- merge(test, csqs, all=F, sort=F)
    
    
    expect_equal(nrow(test.merge), nrow(test))
    
  }
  
  
  
})

test_that("allele_count, sample, variant_id should make sense with respect to the genotype data", {
  
  for(j in seq_along(test.list)){
    
    test.vcf <- test.list[[j]]$obj
    test <- test.list[[j]]$table
    
    for(i in colnames(test.vcf)){
      
      temp.vec <- geno(test.vcf)$GT[,i]
      
      which.ref.het <- which(is.na(temp.vec) == F & grepl("0/[1-9]", temp.vec)  )
      
      expect_true(all(test$allele_count[test$sample == i & test$variant_id %in% which.ref.het] == 1))
      
      which.hom <- which(is.na(temp.vec) == F & grepl("([1-9])/\\1", temp.vec))
      
      expect_true(all(test$allele_count[test$sample == i & test$variant_id %in% which.hom] == 2))
      
      which.multi <- which(is.na(temp.vec) == F & grepl("[1-9]/[1-9]", temp.vec) == T & grepl("([1-9])/\\1", temp.vec) == F)
      
      expect_true(all(test$allele_count[test$sample == i & test$variant_id %in% which.multi] == 1))
      #these should also be duplicated, once per allele
      
      
      expect_true(sum(test$sample == i) == sum(test$sample == i & test$variant_id %in% c(which.ref.het, which.hom, which.multi)))
    }
    
  }
  
  
  
})

test_that("total_reads and allele_reads should make sense with the AD field", {
  
  for(j in seq_along(test.list)){
    
    test.vcf <- test.list[[j]]$obj
    test <- test.list[[j]]$table
    
    for(i in colnames(test.vcf)){
      
      temp.dta <- test[test$sample == i,]
      
      total.reads <- sapply(geno(test.vcf)$AD[,i], sum)
      
      expect_true(all(temp.dta$total_reads == total.reads[temp.dta$variant_id]))
      
      specific_allele_count <- mapply(function(x,y){
        
        x[y]
        
      }, geno(test.vcf)$AD[temp.dta$variant_id,i], temp.dta$ALLELE_NUM+1)
      
      expect_true(all(temp.dta$allele_reads == specific_allele_count))
    }
      
  }
  
})

#the somatic variants don't have allele counts defined for them for now...
test_that("SOMATIC allele_count, sample, variant_id should make sense with respect to the genotype data", {
  
  for(j in seq_along(somatic.list)){
    
    test.vcf <- somatic.list[[j]]$obj
    test <- somatic.list[[j]]$table
    
    i <- colnames(test.vcf)[1]
    
    temp.vec <- geno(test.vcf)$GT[,i]
    
    which.ref.het <- which(is.na(temp.vec) == F & grepl("0/[1-9]", temp.vec)  )
    
    expect_true(all(test$allele_count[test$sample == i & test$variant_id %in% which.ref.het] == 1))
    
    which.hom <- which(is.na(temp.vec) == F & grepl("([1-9])/\\1", temp.vec))
    
    expect_true(all(test$allele_count[test$sample == i & test$variant_id %in% which.hom] == 2))
    
    which.multi <- which(is.na(temp.vec) == F & grepl("[1-9]/[1-9]", temp.vec) == T & grepl("([1-9])/\\1", temp.vec) == F)
    
    expect_true(all(test$allele_count[test$sample == i & test$variant_id %in% which.multi] == 1))
    #these should also be duplicated, once per allele
    
    
    expect_true(sum(test$sample == i) == sum(test$sample == i & test$variant_id %in% c(which.ref.het, which.hom, which.multi)))
    
  }
  
})

test_that("SOMATIC total_reads and allele_reads should make sense with the AD field", {
  
  for(j in seq_along(somatic.list)){
    
    test.vcf <- somatic.list[[j]]$obj
    test <- somatic.list[[j]]$table
    
    i <- colnames(test.vcf)[1]
    
    temp.dta <- test[test$sample == i,]
    
    if ("RD" %in% names(geno(test.vcf))){
        
        total.reads <- sapply(geno(test.vcf)$AD[,i], sum) + sapply(geno(test.vcf)$RD[,i], sum)
        
    }else{
        total.reads <- sapply(geno(test.vcf)$AD[,i], sum)
    }
    
    expect_true(all(temp.dta$total_reads == total.reads[temp.dta$variant_id]))
    
    specific_allele_count <- mapply(function(x,y){
      
      x[y]
      
    }, geno(test.vcf)$AD[temp.dta$variant_id,i], temp.dta$ALLELE_NUM+1)
    
    expect_true(all(temp.dta$allele_reads == specific_allele_count))
      
  }
  
})

