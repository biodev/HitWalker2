#!/bin/bash

Rscript -e 'library(roxygen2)' -e 'roxygenize("hwhelper")'
R CMD build hwhelper

use_file=(`ls hwhelper_*`)

ar_len=${#use_file[@]}

last_pos=$(($ar_len - 1))

R CMD check ${use_file[${last_pos}]}
