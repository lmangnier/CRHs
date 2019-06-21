#!/bin/bash

cd /home/nash/Documents/juicebox/Juicer
 
head -2 TADs_Neu_$1_chr1/$1_blocks.bedpe > TADs_$1_all_blocks.bedpe

for dir in TADs_Neu_$1*
	do 

	if find ${dir} -mindepth 1 | read; then
	
		tail -n +3 -q ${dir}/$1_blocks.bedpe >> TADs_$1_all_blocks.bedpe
	else 
		echo "${dir} is empty" 
	fi
	done
