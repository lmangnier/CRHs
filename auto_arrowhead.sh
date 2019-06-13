#!/bin/bash

chrs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"

for chr in $chrs
	do 
		java -jar ./juicer_tools_1.11.04_jcuda.0.8.jar arrowhead -r 5000 -c "$chr"  ../Neu_inter_sample_372787143.hic ../TADs_Neu_5000_chr"$chr"

	done
