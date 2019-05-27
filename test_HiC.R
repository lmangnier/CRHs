#Hi-C_analysis

## Quick Description 
# his project gathers all the chunks of code in relation to Hi-C data analysis.
#Data were retrieved from Rajarajan et al 2018 (1) and available on Psychencode Synapse Platform.
#Our goal is to integrate the non-coding regions of the genome in rare association test variants, with respect to 3D genome contacts. 
#So, first of all, we start by annoting genes and enhancers which are present in contact locus from the Hi-C file. 
#Some statistical analysis were made to detect peak enrichment in Hi-C 3D contact data with HiCCUPS(2) software. 
#HiCCUPS implements methodology for detecting peaks provided from Rao et al 2014.

#ALL useful librairies for annotations of contacts regions and gene-enhancer clusters analysis
library(GenomicRanges)
library(IRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(rtracklayer)
library(karyoploteR)
library(ggplot2)
library(igraph)
library(FantomEnhancers.hg19)

setwd("/home/nash/Documents/Psychencode/data/")

clean.hic <- read.table("clean_master_hiccups_loops_nonsubsampled.txt",header = TRUE, sep="\t")
#As we're working on schizophrenia an bipolar troubles, we just retrieve data from neurons
clean.hic.neu <- clean.hic[clean.hic$Cell == "Neu",]

#Let's take a look on data
head(clean.hic.neu)


significant.contact <- function(fdr_bottom_left, fdr_donut, fdr_horizontal, fdr_vertical, threshold){
  #Function which determines the number of significant contacts in relation to a specific FDR threshold
  #The function returns a 0,1 vector
  ifelse(fdr_bottom_left <= threshold && fdr_donut <= threshold && fdr_horizontal <= threshold && fdr_vertical<= threshold, 1,0)
}

sum(mapply(significant.contact, clean.hic.neu$fdrBL, clean.hic.neu$fdrDonut, clean.hic.neu$fdrH, clean.hic.neu$fdrV, .10))

#100% of pixels provided from the file are enriched. with a 10% FDR thresold, this value may changed with a 5% threshold or 1% threshold.
#We will do the same analysis with different threshold's values
sum(mapply(significant.contact, clean.hic.neu$fdrBL, clean.hic.neu$fdrDonut, clean.hic.neu$fdrH, clean.hic.neu$fdrV, .05))
sum(mapply(significant.contact, clean.hic.neu$fdrBL, clean.hic.neu$fdrDonut, clean.hic.neu$fdrH, clean.hic.neu$fdrV, .01))
#The rule is extracted of Rao and al 2014.

#Creation of one column for each specific rule 
clean.hic.neu$FDR_0.1 <- mapply(significant.contact, clean.hic.neu$fdrBL, clean.hic.neu$fdrDonut, clean.hic.neu$fdrH, clean.hic.neu$fdrV, .10)
clean.hic.neu$FDR_0.05 <- mapply(significant.contact, clean.hic.neu$fdrBL, clean.hic.neu$fdrDonut, clean.hic.neu$fdrH, clean.hic.neu$fdrV, .05)
clean.hic.neu$FDR_0.01 <-mapply(significant.contact, clean.hic.neu$fdrBL, clean.hic.neu$fdrDonut, clean.hic.neu$fdrH, clean.hic.neu$fdrV, .01)

head(clean.hic.neu)

#Analysis of FDR's correlation
FDR.combn <- combn(colnames(clean.hic.neu[, c("fdrBL","fdrDonut","fdrH","fdrV")]),2)
lapply(1:ncol(FDR.combn), FUN = function(x) paste(FDR.combn[,x][1],",",FDR.combn[,x][2],": ",cor(clean.hic.neu[,FDR.combn[,x][1]],clean.hic.neu[,FDR.combn[,x][2]],method = "pearson")))
#There is no strong correlation between the fourth FDRs.

#Are there trans contacts ?
table(clean.hic.neu$chr1, clean.hic.neu$chr2)

#Now let's start by gene annotations
x <- clean.hic.neu[,c("chr1","x1","x2")]
y <- clean.hic.neu[,c("chr2","y1","y2")]

#Creation of GRanges from contact locus
#We'll serve of contact locus to detect gene-enhancer relations into the tissue
#I generate 4 files to make reverse analysis
GRanges.X <- GRanges(seqnames = x$chr1, ranges = IRanges(x$x1, x$x2))
GRanges.X1 <- GRanges(seqnames = x$chr1, ranges = IRanges(x$x1, x$x2))

GRanges.Y <- GRanges(seqnames = y$chr2, ranges = IRanges(y$y1, y$y2))
GRanges.Y1 <- GRanges(seqnames = y$chr2, ranges = IRanges(y$y1, y$y2))

head(GRanges.X);head(GRanges.Y)

#Unique coverage from contact locus for the both of files
paste("Total genome coverage : ",sum(width(reduce(GRanges.X))), "bp")
paste("Total genome coverage : ",sum(width(reduce(GRanges.Y))), "bp")

#Hg19 genes
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes.hg19 <- genes(txdb)

symbol <- select(org.Hs.eg.db,keys = genes.hg19$gene_id, columns = c("SYMBOL"), keytype = "ENTREZID")

genes.hg19$geneSymbol <- symbol$SYMBOL 
#We're taking the TSS and TES from genes for the post analysis
genes.hg19$TSS <- start(genes)
genes.hg19$TES <- end(genes)

head(genes.hg19)

#Number of genes overlaps on contacts regions
genes.overlaps.X <- subsetByOverlaps(genes.hg19, GRanges.X)

#Unique Gene Coverage by the overlapping genes 
sum(width(reduce(genes.overlaps.X)))
#Mean of width Gene coverage from overlapping genes
mean(width(genes.overlaps.X))

#density for overlapping genes all along chromosomes 
kp <- plotKaryotype(genome="hg19")
kpPlotDensity(kp, data=genes.overlaps.X)

#Here, we do the genes annotation from contact regions
m <- findOverlaps(genes, GRanges.X)

mcols(GRanges.X)$geneSymbol <- NA
mcols(GRanges.X)$TSS <- NA
mcols(GRanges.X)$TES <- NA

mcols(GRanges.X)[subjectHits(m), "geneSymbol"] <- mcols(genes.hg19)[queryHits(m), "geneSymbol"]
mcols(GRanges.X)[subjectHits(m), "TSS"] <- mcols(genes.hg19)[queryHits(m), "TSS"]
mcols(GRanges.X)[subjectHits(m), "TES"] <- mcols(genes.hg19)[queryHits(m), "TES"]

head(GRanges.X)

#Have enhancers for each genes, data are extracted from Fantom5 
load("neurons_enhancers.RData")

#Let's take a look on data 
head(neurons_enhancers)

expr <- as.data.frame(mcols(neurons_enhancers))

#3 different filter steps were done 
filter.expr <- rowSums(expr) > 0 #At least one  enhancer is expressed 
maj.enh <- function(x,y,z) ifelse((x > 0 && y > 0) || (x > 0 && z> 0) || (z > 0 && y > 0) || (x > 0 && y > 0 && z > 0), TRUE, FALSE) #At least two enhancers are expressed
consensus.enh <- function(x,y,z) ifelse(x > 0 && y > 0 && z > 0, TRUE, FALSE) # All of three enhancers are expressed

sum(filter.expr)
sum(mapply(maj.enh, expr$CNhs12338, expr$CNhs12726, expr$CNhs13815)) #For the first step, we choose this filter rule
sum(mapply(consensus.enh, expr$CNhs12338, expr$CNhs12726, expr$CNhs13815))

expr.enhanc <- neurons_enhancers[filter.expr]
expr.enhanc.3 <- neurons_enhancers[mapply(consensus.enh, expr$CNhs12338, expr$CNhs12726, expr$CNhs13815)]

#Expressed enhancer in at least two cell type
expr.enhanc.2 <- neurons_enhancers[mapply(maj.enh, expr$CNhs12338, expr$CNhs12726, expr$CNhs13815)]

#Annotations for enhancers
mcols(expr.enhanc.2)$enhancerstart <- start(expr.enhanc.2)
mcols(expr.enhanc.2)$enhancerstop <- end(expr.enhanc.2)

t <- findOverlaps(expr.enhanc.2, GRanges.Y)
mcols(GRanges.Y)$enhancerstart <- NA
mcols(GRanges.Y)$enhancerstop <- NA

mcols(GRanges.Y)[subjectHits(t),"enhancerstart"] <- mcols(expr.enhanc.2)[queryHits(t),"enhancerstart"]
mcols(GRanges.Y)[subjectHits(t),"enhancerstop"] <- mcols(expr.enhanc.2)[queryHits(t), "enhancerstop"]

head(GRanges.Y)

#Because of the order of bins is preserved  we index the 2 files for retrieving the genes associated with enhancers
names(GRanges.X) <- 1:length(GRanges.X)
names(GRanges.Y) <- 1:length(GRanges.Y)

#All index with a non null enhancer
names(GRanges.Y[!is.na(mcols(GRanges.Y)$enhancerstart)])
length(names(GRanges.Y[!is.na(mcols(GRanges.Y)$enhancerstart)]))
#180 enhancers are expressed on overlapping locus 

asso.genes <- GRanges.X[names(GRanges.Y[!is.na(mcols(GRanges.Y)$enhancerstart)])]
asso.genes <- asso.genes[!is.na(mcols(asso.genes)$geneSymbol)]
index.non.null.genes <- names(asso.genes)
asso.enhanc <- GRanges.Y[index.non.null.genes]

#The methodology seems to work
asso.genes;asso.enhanc

#We gonna do the same analysis but with genes on y file and enhancers on x file 
#This step has to be automated 

p <- findOverlaps(genes, GRanges.Y1)

mcols(GRanges.Y1)$geneSymbol <- NA
mcols(GRanges.Y1)$TSS <- NA
mcols(GRanges.Y1)$TES <- NA

mcols(GRanges.Y1)[subjectHits(p), "geneSymbol"] <- mcols(genes.hg19)[queryHits(p), "geneSymbol"]
mcols(GRanges.Y1)[subjectHits(p), "TSS"] <- mcols(genes.hg19)[queryHits(p), "TSS"]
mcols(GRanges.Y1)[subjectHits(p), "TES"] <- mcols(genes.hg19)[queryHits(p), "TES"]

head(GRanges.Y1)

q <- findOverlaps(expr.enhanc.2, GRanges.X1)
mcols(GRanges.X1)$enhancerstart <- NA
mcols(GRanges.X1)$enhancerstop <- NA

mcols(GRanges.X1)[subjectHits(q),"enhancerstart"] <- mcols(expr.enhanc.2)[queryHits(q),"enhancerstart"]
mcols(GRanges.X1)[subjectHits(q),"enhancerstop"] <- mcols(expr.enhanc.2)[queryHits(q), "enhancerstop"]

head(GRanges.X1)

names(GRanges.X1) <- 1:length(GRanges.X1)
names(GRanges.Y1) <- 1:length(GRanges.Y1)

asso.genes.1 <- GRanges.Y1[names(GRanges.X1[!is.na(mcols(GRanges.X1)$enhancerstart)])]
asso.genes.1 <- asso.genes.1[!is.na(mcols(asso.genes.1)$geneSymbol)]
index.non.null.genes.1 <- names(asso.genes.1)
asso.enhanc.1 <- GRanges.X1[index.non.null.genes.1]


#Pairs of gene-enhancer in overlapping regions
df.genes <- as.data.frame(mcols(asso.genes))
df.genes.1 <- as.data.frame(mcols(asso.genes.1))
df.genes <- rbind(df.genes, df.genes.1)

df.enhanc <- as.data.frame(mcols(asso.enhanc))
df.enhanc.1 <- as.data.frame(mcols(asso.enhanc.1))
df.enhanc <- rbind(df.enhanc, df.enhanc.1)

df.genes.enh <- cbind(df.genes, df.enhanc)

length(unique(df.genes.enh$geneSymbol))
length(unique(df.genes.enh$enhancerstart, df.genes.enh$enhancerstop))

#In previous work we work on MCF7 data extracted from Wu et al 2018 to analyze cluster of gene-enhancer based on their contacts
#We apply MCF7 methodology here...
#Non redundant Genome coverage by enhancers
enhancers <- IRanges(start = df.genes.enh$enhancerstart , end = df.genes.enh$enhancerstop)
sum(width(reduce(enhancers)))

df.genes.enh$enhancer <- paste(df.genes.enh$enhancerstart, df.genes.enh$enhancerstop)

agg.table.enhancers <- aggregate(geneSymbol~enhancer, data=df.genes.enh, unique, na.rm=TRUE)
agg.table.genes <- aggregate(enhancer~geneSymbol, data=df.genes.enh, unique, na.rm=TRUE)

table(sapply(agg.table.enhancers$geneSymbol, length))
table(sapply(agg.table.genes$enhancer, length))

barplot(table(sapply(agg.table.enhancers$geneSymbol, length)), main = "Distribution of number of genes associated by enhancer", xlab = "Number of genes associated", ylab = "Frequency")
barplot(table(sapply(agg.table.genes$enhancer, length)), main = "Distribution of number of enhancers associated by genes", xlab = "Number of enhancers associated", ylab = "Frequency")

#Overlapping analysis
test.overlapping <- function(TSS.gene,TES.gene, enhS, enhE) {
  if(TSS.gene > enhE){
    
    return("downstream")
  }
  else if (TES.gene < enhS){
    
    return("upstream")
  }
  else{
    
    return("overlapping")
  }
}

df.genes.enh$overlapping <- mapply(test.overlapping, df.genes.enh$TSS, df.genes.enh$TES, df.genes.enh$enhancerstart, df.genes.enh$enhancerstop)
head(df.genes.enh)

barplot(table(df.genes.enh$overlapping), main = "Barplot of enhancer position in relation to his gene")

#Quantile Coverage analysis by overlapping status 
df.genes.enh$width.enh <- df.genes.enh$enhancerstop - df.genes.enh$enhancerstart
quant.overlapping <- aggregate(width.enh ~ overlapping, data=df.genes.enh, FUN = "quantile" ,probs=c(0.25, 0.50,0.75))
quant.overlapping

hist(df.genes.enh$width.enh)

contact.func <- function(col.to.agg, data){
  library(data.table)
  tmp <- copy(data)
  
  #If column which wants to aggregate has a type list, one more conversion step is needed
  tmp$col1.bis <- as.character(lapply(tmp[,col.to.agg], FUN = function(x) paste(x, collapse = "-")))
  
  if(col.to.agg == "geneSymbol") {
    tmp.agg.qual <- aggregate(enhancer ~ col1.bis, data = tmp, length)
    tmp.agg.qual$N.genes <- as.numeric(lapply(strsplit(tmp.agg.qual$col1.bis, "-"), length))
    tmp.agg.quant <- as.data.frame(table(tmp.agg.qual$N.genes, tmp.agg.qual$enhancer))
    colnames(tmp.agg.quant) <- c("N.Genes", "N.Enhancers", "N")
    
    tmp.agg.quant <- tmp.agg.quant[!tmp.agg.quant$N == 0,]
    
  }
  
  else{
    tmp.agg.qual <- aggregate(geneSymbol ~ col1.bis, data = tmp, length)
    tmp.agg.qual$N.genes <- as.numeric(lapply(strsplit(tmp.agg.qual$col1.bis, "-"), length))
    tmp.agg.quant <- as.data.frame(table(tmp.agg.qual$N.genes, tmp.agg.qual$geneSymbol))
    colnames(tmp.agg.quant) <- c("N.Enhancers", "N.Genes", "N")
    
    tmp.agg.quant <- tmp.agg.quant[!tmp.agg.quant$N == 0,]
  }
  
  
  list("qualitative_asso"= tmp.agg.qual, "quantitative_asso"= tmp.agg.quant)
}

analysis.g_e <- contact.func("geneSymbol",agg.table.enhancers)
head(analysis.g_e$quantitative_asso)

analysis.e_g <- contact.func("enhancer", agg.table.genes)
head(analysis.e_g$quantitative_asso)

#Linear mapping between genes and enhancers
plot(df.genes.enh$TSS, df.genes.enh$enhancerstop)
#It seems that there is a perfectly linear relation between genes and enhancers

#Network analysis from gene-enhancer clusters 
genes <- matrix(df.genes.enh$geneSymbol, ncol = 1)
enh <- matrix(df.genes.enh$enhancer, ncol = 1)

nodes.g_e <- as.data.frame(rbind(genes,enh))
colnames(nodes.g_e) <- "id"

nodes.g_e[1:nrow(genes),"type"] <- "gene"
nodes.g_e[1:nrow(genes),"col"] <- "orange"

nodes.g_e[nrow(genes)+1:nrow(nodes.g_e),"type"] <- "enhancer"
nodes.g_e[nrow(genes)+1:nrow(nodes.g_e),"col"] <- "blue"

links.g_e <- df.genes.enh[,c("geneSymbol", "enhancer")]
links.g_e$weight <- 1

nodes.g_e <- unique(na.omit(nodes.g_e))

net <- graph_from_data_frame(d=links.g_e,vertices = nodes.g_e ,directed = F)
plot(net, vertex.label=NA, vertex.size = 2, margin=-.1,asp=.35, vertex.color = nodes.g_e$col, edge.color = "red",main="Clustering based on gene-enhancer contacts")
legend(x=0, y=-1.1, c("gene","enhancer"), pch=21,
        col="#777777", pt.bg=unique(nodes.g_e$col), pt.cex=2, cex=.8, bty="n", ncol=1)

