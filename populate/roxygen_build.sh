#!/bin/bash

if [ $# -lt 1 ] || [ $1 != "check" -a $1 != "install" ]
then
    echo "Need to specify either 'check' or 'install'"
    exit
fi

Rscript -e 'library(roxygen2)' -e 'roxygenize("hwhelper")'
R CMD build hwhelper

use_file=(`ls hwhelper_*`)

ar_len=${#use_file[@]}

last_pos=$(($ar_len - 1))

if [ $1 == "check" ]
then

    R CMD check ${use_file[${last_pos}]}

else
    R CMD INSTALL ${use_file[${last_pos}]}

fi
