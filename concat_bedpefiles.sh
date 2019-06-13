#!/bin/bash

cd /home/nash/Documents/juicebox/Juicer
 
head -2 TADs_Neu_5000_chr1/5000_blocks.bedpe > TADs_5000_all_blocks.bedpe

for dir in TADs_Neu*
	do 

	if find ${dir} -mindepth 1 | read; then
	
		tail -n +3 -q ${dir}/5000_blocks.bedpe >> TADs_5000_all_blocks.bedpe
	else 
		echo "${dir} is empty" 
	fi
	done
