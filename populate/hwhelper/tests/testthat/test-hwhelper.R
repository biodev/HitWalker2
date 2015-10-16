
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