#!/bin/bash

if [ $# -ne 1 ]
then

	echo "Need to supply the prefix of a config/custom_functions file"	
	exit 1
fi

if [ -f network/config.py ]
then
	rm network/config.py
fi

if [ -f network/custom_functions.py ]
then

	rm network/custom_functions.py
fi 

cp network/$1_config.py network/config.py
cp network/$1_custom_functions.py network/custom_functions.py
