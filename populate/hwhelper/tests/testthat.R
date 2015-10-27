library(testthat)
library(VariantAnnotation)
library(ensemblVEP)

###Example of how the below file(s) were processed:
# base.dir <- "/Volumes/Macintosh_HD_4/"
# 
# vcf_input <- "CEUTrio.HiSeq.WGS.b37.NA12878.vcf"
# vcf_output <- sub(".vcf",".refseq.vcf", vcf_input)
#   
# Sys.setenv(PERL5LIB="/Volumes/Macintosh_HD_4/vep/")
#   
# vep.runner <- paste0(base.dir, "vep/ensembl-tools-release-78/scripts/variant_effect_predictor/variant_effect_predictor.pl
#                        --fork 4 --species homo_sapiens --input_file ", vcf_input, " --format vcf --output_file ", vcf_output, 
#                        " --stats_text --cache --dir ", base.dir,"vep/vep --offline --everything --check_existing --total_length --allele_number
#                        --refseq --vcf --flag_pick_allele --pick_order biotype,rank,length,canonical,tsl ")
#   
# system(gsub("\n", " ", vep.runner))
#   
# 
# #maybe only do the HGVSp for HitWalker and company, not for general use...
# 
# 
# input_name <- "CEUTrio.HiSeq.WGS.b37.NA12878.refseq.vcf"
# output_name <- sub("refseq", "refseq.filt", input_name)
# 
# filter.runner <- paste0('perl /Volumes/Macintosh_HD_4/vep/ensembl-tools-release-78/scripts/variant_effect_predictor/filter_vep.pl --format vcf --input_file ',
#                         input_name, ' --output_file ',output_name, ' --only_matched --filter "PICK and HGVSp"')
# 
# system(gsub("\n", " ", filter.runner))

test_check("hwhelper")

