# Topologically associated domains (TADs) analysis

library(rtracklayer)
library(mosaic)

# This script assumes that the script HiC_V*.R has been run before

TADS50 <- import("TADs/Neu.50000_blocks.bed", format="bed")
TADS100 <- import("TADs/Neu(1).100000_blocks.bed", format="bed")
length(TADS50)
#[1] 2525
length(TADS100)
#[1] 1248
# Much fewer TADs of 100kb blocks
sum(width(reduce(TADS50)))
[1] 2316700000
sum(width(reduce(TADS100)))
[1] 2205200000
# Genomic coverage is about the same for 50kb and 100kb blocks

# Creating a GRanges object with the Epigenomics roadmap-based clusters

# Extract chromosome for each cluster in 1-22,X format from chr1-chr22,chrX format
chr.cluster = tapply(genes.enhancers.EGRM$chr,compo.epg$membership[genes.enhancers.EGRM$geneSymbol],function(vec) substring(vec[1],4))
cluster.GRanges = GRanges(seqnames = chr.cluster,ranges = IRanges(start.cluster,end.cluster))

# Number of TADs overlapped by each cluster
cluster.TADs50 <- findOverlaps(cluster.GRanges, TADS50)
table(countLnodeHits(cluster.TADs50))
#0   1   2   3   4   5   6 
#173 373 218  83  25   8   1 
round(tally(countLnodeHits(cluster.TADs50),format="percent"))
#0  1  2  3  4  5  6 
#20 42 25  9  3  1  0 

cluster.TADs100 <- findOverlaps(cluster.GRanges, TADS100)
table(countLnodeHits(cluster.TADs100))
#0   1   2   3   4   5 
#240 349 202  68  18   4 
round(tally(countLnodeHits(cluster.TADs100),format="percent"))
#0  1  2  3  4  5 
#27 40 23  8  2  0 

# Gene-enhancer Epigenomics Roadmap pairs overlap with TADs
vec.X = ifelse(genes.enhancers.EGRM$TSS<=genes.enhancers.EGRM$enhancerstop,genes.enhancers.EGRM$TSS, genes.enhancers.EGRM$enhancerstart)
vec.Y = ifelse(genes.enhancers.EGRM$TSS>genes.enhancers.EGRM$enhancerstop,genes.enhancers.EGRM$TES, genes.enhancers.EGRM$enhancerstop)
genes.enhancers.EGRM.GRanges=GRanges(substring(genes.enhancers.EGRM$chr,4),ranges = IRanges(vec.X,vec.Y))

# Number of TADs overlapped by each gene-enhancer pair
genes.enhancers.EGRM.TADs50 <- findOverlaps(genes.enhancers.EGRM.GRanges, TADS50)
gecTADs50 = countLnodeHits(genes.enhancers.EGRM.TADs50)
table(gecTADs50)
#0   1   2   3   4   5   6 
#256 560 327 134  29   8   1 
round(tally(gecTADs50,format="percent"))
#0  1  2  3  4  5  6 
#19 43 25 10  2  1  0 
round(tally(gecTADs50[gecTADs50>0],format="percent"))
#1  2  3  4  5  6 
#53 31 13  3  1  0

genes.enhancers.EGRM.TADs100 <- findOverlaps(genes.enhancers.EGRM.GRanges, TADS100)
gecTADs100 = countLnodeHits(genes.enhancers.EGRM.TADs100)
table(gecTADs100)
#0   1   2   3   4   5 
#363 520 291 111  25   5 
round(tally(gecTADs100,format="percent"))
#0  1  2  3  4  5 
#28 40 22  8  2  0 
round(tally(gecTADs100[gecTADs100>0],format="percent"))
#1  2  3  4  5 
#55 31 12  3  1 

# The distributions of number of TADs overlapped are very similar for the clusters
# and the gene-enhancer pairs, even though the clusters are longer than the pairs

