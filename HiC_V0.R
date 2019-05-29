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
library(biomaRt)
library(FantomEnhancers.hg19)

clean.hic <- read.table(file.choose(),header = TRUE, sep="\t")

#Let's take a look on initial data
head(clean.hic.neu)

#As we're working on schizophrenia an bipolar troubles, we just retrieve data from neurons
clean.hic.neu <- clean.hic[clean.hic$Cell == "Neu",]

#Some functions to facilitate analysis 

significant.contact <- function(fdr_bottom_left, fdr_donut, fdr_horizontal, fdr_vertical, threshold){
  #Function which determines the number of significant contacts in relation to a specific FDR threshold
  #The function returns a 0,1 vector
  ifelse(fdr_bottom_left <= threshold && fdr_donut <= threshold && fdr_horizontal <= threshold && fdr_vertical<= threshold, 1,0)
}

correlation.function <- function(df,vector.of.cols,method = "pearson"){
  #This function computes the correlation between a set of columns by pairs
  tmp.combn <- combn(colnames(df[, vector.of.cols]),2)
  lapply(1:ncol(tmp.combn), FUN = function(x) paste(tmp.combn[,x][1],",",tmp.combn[,x][2],": ",cor(df[,tmp.combn[,x][1]],df[,tmp.combn[,x][2]],method = method)))
}

create.reverse.GRanges <- function(x,y){
  #Creation of GRanges from contact locus
  #We'll serve of contact locus to detect gene-enhancer relations into the tissue
  #I generate 4 files to make reverse analysis
  
  GRanges.X <- GRanges(seqnames = x$chr1, ranges = IRanges(x$x1, x$x2))
  GRanges.X1 <- GRanges(seqnames = x$chr1, ranges = IRanges(x$x1, x$x2))
  
  GRanges.Y <- GRanges(seqnames = y$chr2, ranges = IRanges(y$y1, y$y2))
  GRanges.Y1 <- GRanges(seqnames = y$chr2, ranges = IRanges(y$y1, y$y2))
  
  list("X"= GRanges.X , "Y"= GRanges.Y, "X1" = GRanges.X1, "Y1" = GRanges.Y1)
}

genes.annotations <- function(genes, GRange.to.annotate){
  #Here, we do the genes annotation from contact regions
  overlaps.index <- findOverlaps(genes, GRange.to.annotate)
  
  #Creation of geneSymbol, TSS et TES for genes
  mcols(GRange.to.annotate)$geneSymbol <- NA
  mcols(GRange.to.annotate)$TSS <- NA
  mcols(GRange.to.annotate)$TES <- NA
  
  #Correspondance between contacts locus and genes which overlaps this latter
  mcols(GRange.to.annotate)[subjectHits(overlaps.index), "geneSymbol"] <- mcols(genes)[queryHits(overlaps.index), "geneSymbol"]
  mcols(GRange.to.annotate)[subjectHits(overlaps.index), "TSS"] <- mcols(genes)[queryHits(overlaps.index), "TSS"]
  mcols(GRange.to.annotate)[subjectHits(overlaps.index), "TES"] <- mcols(genes)[queryHits(overlaps.index), "TES"]
  
  GRange.annotated <- GRange.to.annotate
  return(GRange.annotated)
}

filter.enhancers <- function(enhancer.GRange, rule=c("one","consensus","total")){
  #This function takes in input a enhancer file from Fantom5
  #Retrieving of expression column
  expr <- as.data.frame(mcols(enhancer.GRange))
  
  if(rule == "one"){
    #At least one enhancer is expressed into replicat cells
    filter.one <- rowSums(expr) > 0 
    return(enhancer.GRange[filter.one])
  }
  #At least two enhancers are expressed
  else if(rule == "consensus"){
    filter.cons <- function(x,y,z) ifelse((x > 0 && y > 0) || (x > 0 && z> 0) || (z > 0 && y > 0) || (x > 0 && y > 0 && z > 0), TRUE, FALSE)
    return(enhancer.GRange[mapply(filter.cons, expr$CNhs12338, expr$CNhs12726, expr$CNhs13815)])
  }
  
  #All the enhancers are expressed
  else if(rule =="total"){
    filter.total <- function(x,y,z) ifelse(x > 0 && y > 0 && z > 0, TRUE, FALSE)
    return(enhancer.GRange[mapply(filter.total, expr$CNhs12338, expr$CNhs12726, expr$CNhs13815)])
  }
}

find.medoid <- function(x,y){
  return(median(seq(x,y)))
}

enhancers.annotations <- function(enhancers, GRange.to.annotate){
  
  #Annotations for enhancers
  mcols(enhancers)$enhancerstart <- start(enhancers)
  mcols(enhancers)$enhancerstop <- end(enhancers)
  
  overlaps.index <- findOverlaps(enhancers, GRange.to.annotate)
  
  mcols(GRange.to.annotate)$enhancerstart <- NA
  mcols(GRange.to.annotate)$enhancerstop <- NA
  
  mcols(GRange.to.annotate)[subjectHits(overlaps.index),"enhancerstart"] <- mcols(enhancers)[queryHits(overlaps.index),"enhancerstart"]
  mcols(GRange.to.annotate)[subjectHits(overlaps.index),"enhancerstop"] <- mcols(enhancers)[queryHits(overlaps.index), "enhancerstop"]
  
  GRange.annotated <- GRange.to.annotate
  return(GRange.annotated)
  
}

test.overlapping <- function(TSS.gene,TES.gene, enhS, enhE) {
  #Overlapping analysis
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
correlation.function(clean.hic.neu, c("fdrBL","fdrDonut","fdrH","fdrV"))
#There is no strong correlation between the fourth FDRs.

#Are there trans contacts ?
table(clean.hic.neu$chr1, clean.hic.neu$chr2)

#Now let's start by gene annotations
#I create two datasets come from .hic file
x <- clean.hic.neu[,c("chr1","x1","x2")]
y <- clean.hic.neu[,c("chr2","y1","y2")]

contacts.locus <- create.reverse.GRanges(x,y)

head(contacts.locus$X);head(contacts.locus$Y)

#Unique coverage from contact locus for the both of files
paste("Total genome coverage : ",sum(width(reduce(contacts.locus$X))), "bp")
paste("Total genome coverage : ",sum(width(reduce(contacts.locus$Y))), "bp")

#Hg19 genes
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes.hg19 <- genes(txdb)

symbol <- select(org.Hs.eg.db,keys = genes.hg19$gene_id, columns = c("SYMBOL"), keytype = "ENTREZID")

genes.hg19$geneSymbol <- symbol$SYMBOL 
#We're taking the TSS and TES from genes for the post analysis
genes.hg19$TSS <- start(genes.hg19)
genes.hg19$TES <- end(genes.hg19)

head(genes.hg19)

#Number of gene overlaps on contact regions
counts.overlaps.X <- countOverlaps(genes.hg19, contacts.locus$X)
counts.overlaps.Y <- countOverlaps(genes.hg19, contacts.locus$Y)

table(counts.overlaps.X)
table(counts.overlaps.Y)

genes.overlaps.X <- subsetByOverlaps(genes.hg19, contacts.locus$X)
genes.overlaps.Y <- subsetByOverlaps(genes.hg19, contacts.locus$Y)

#Unique Gene Coverage by the overlapping genes 
sum(width(reduce(genes.overlaps.X)))
sum(width(reduce(genes.overlaps.Y)))

#Mean of width Gene coverage from overlapping genes
mean(width(genes.overlaps.X))
mean(width(genes.overlaps.Y))

quantiles.overlapping.X <- quantile(width(genes.overlaps.X), probs=c(.25,.50,.75))
quantiles.overlapping.Y <- quantile(width(genes.overlaps.Y), probs=c(.25,.50,.75))
  
hist(width(genes.overlaps.X))
hist(width(genes.overlaps.Y))

#density for overlapping genes all along chromosomes 
kpX <- plotKaryotype(genome="hg19")
kpPlotDensity(kpX, data=genes.overlaps.X)

kpY <- plotKaryotype(genome="hg19")
kpPlotDensity(kpY, data=genes.overlaps.Y)

#Have enhancers for each genes, data are extracted from Fantom5 
load("neurons_enhancers.RData")
#Let's take a look on data 
head(neurons_enhancers)

genes.annotations.X <- genes.annotations(genes.hg19,contacts.locus$X)
genes.annotations.Y <- genes.annotations(genes.hg19,contacts.locus$Y1)

filter.enhancers.cons <- filter.enhancers(neurons_enhancers, "consensus")
filter.enhancers.one <- filter.enhancers(neurons_enhancers, "one")

enhancers.annotations.Y <- enhancers.annotations(filter.enhancers.one,contacts.locus$Y)
enhancers.annotations.X <- enhancers.annotations(filter.enhancers.one,contacts.locus$X1)

#Because of the order of bins is preserved  we index the 2 type files for retrieving the genes associated with enhancers
names(genes.annotations.X) <- 1:length(genes.annotations.X)
names(enhancers.annotations.Y) <- 1:length(enhancers.annotations.Y)
names(enhancers.annotations.X) <- 1:length(enhancers.annotations.X)
names(genes.annotations.Y) <- 1:length(genes.annotations.Y)

asso.genes.XY <- genes.annotations.X[names(enhancers.annotations.Y[!is.na(mcols(enhancers.annotations.Y)$enhancerstart)])]
asso.genes.XY <- asso.genes.XY[!is.na(mcols(asso.genes.XY)$geneSymbol)]
index.non.null.genes.XY <- names(asso.genes.XY)
asso.enhanc.XY <- enhancers.annotations.Y[index.non.null.genes.XY]

asso.genes.XY.1 <- genes.annotations.Y[names(enhancers.annotations.X[!is.na(mcols(enhancers.annotations.X)$enhancerstart)])]
asso.genes.XY.1 <- asso.genes.XY.1[!is.na(mcols(asso.genes.XY.1)$geneSymbol)]
index.non.null.genes.XY.1 <- names(asso.genes.XY.1)
asso.enhanc.XY.1 <- enhancers.annotations.X[index.non.null.genes.XY.1]

#Enhancer contacts in respect of contacts regions 
#Here, we're looking after clusters of enhancers because mostly of genes inside have the same function
#Epigenomic roadmap consortium, 2015
tt <- enhancers.annotations.X[names(enhancers.annotations.Y[!is.na(mcols(enhancers.annotations.Y)$enhancerstart)])]
tt <- tt[!is.na(mcols(tt)$enhancerstart)]
index.tt <- names(tt)
tt.1 <- enhancers.annotations.Y[index.tt]

#Pairs of gene-enhancer in overlapping regions
df.genes.1 <- as.data.frame(mcols(asso.genes.XY))
df.genes.2<- as.data.frame(mcols(asso.genes.XY.1))

df.genes <- rbind(df.genes.1, df.genes.2)

df.enhanc.1 <- as.data.frame(mcols(asso.enhanc.XY))
df.enhanc.2 <- as.data.frame(mcols(asso.enhanc.XY.1))

df.enhanc <- rbind(df.enhanc.1, df.enhanc.2)
df.enhanc.c <- cbind(df.enhanc.1, df.enhanc.2)
df.genes.enh <- cbind(df.genes, df.enhanc)

head(df.genes.enh)

length(unique(df.genes.enh$geneSymbol))
length(unique(df.genes.enh$enhancerstart, df.genes.enh$enhancerstop))

df.t <- as.data.frame(mcols(tt))
df.tt <- as.data.frame(mcols(tt.1))

enhancer.contacts <- cbind(df.t, df.tt)
head(enhancer.contacts)
colnames(enhancer.contacts) <- c("enhancerstart1","enhancerstop1","enhancerstart2","enhancerstop2")

enhancer.contacts$enhancer1 <- paste(enhancer.contacts$enhancerstart1, enhancer.contacts$enhancerstop1)
enhancer.contacts$enhancer2 <- paste(enhancer.contacts$enhancerstart2, enhancer.contacts$enhancerstop2)

enhancer.contacts$med1 <- mapply(find.medoid,enhancer.contacts$enhancerstart1, enhancer.contacts$enhancerstop1)
enhancer.contacts$med2 <- mapply(find.medoid,enhancer.contacts$enhancerstart2, enhancer.contacts$enhancerstop2)

#Quantiles of distance between two connected enhancers
quantile(abs(enhancer.contacts$med1 - enhancer.contacts$med2), probs=c(.25,.5,.75))

plot(enhancer.contacts$enhancerstart1, enhancer.contacts$enhancerstop2)

length(unique(enhancer.contacts$enhancerstart1, enhancer.contacts$enhancerstop1))
length(unique(enhancer.contacts$enhancerstart2, enhancer.contacts$enhancerstop2))

#In previous work we work on MCF7 data extracted from Wu et al 2018 to analyze cluster of gene-enhancer based on their contacts
#We apply MCF7 methodology here...
#Non redundant Genome coverage by enhancers
enhancers <- IRanges(start = df.genes.enh$enhancerstart , end = df.genes.enh$enhancerstop)


df.genes.enh$enhancer <- paste(df.genes.enh$enhancerstart, df.genes.enh$enhancerstop)
plot(df.genes.enh$enhancerstart, df.genes.enh$enhancerstop)
agg.table.enhancers <- aggregate(geneSymbol~enhancer, data=df.genes.enh, unique, na.rm=TRUE)
agg.table.genes <- aggregate(enhancer~geneSymbol, data=df.genes.enh, unique, na.rm=TRUE)

table(sapply(agg.table.enhancers$geneSymbol, length))
table(sapply(agg.table.genes$enhancer, length))

barplot(table(sapply(agg.table.enhancers$geneSymbol, length)), main = "Distribution of number of genes associated by enhancer", xlab = "Number of genes associated", ylab = "Frequency")
barplot(table(sapply(agg.table.genes$enhancer, length)), main = "Distribution of number of enhancers associated by genes", xlab = "Number of enhancers associated", ylab = "Frequency")

df.genes.enh$overlapping <- mapply(test.overlapping, df.genes.enh$TSS, df.genes.enh$TES, df.genes.enh$enhancerstart, df.genes.enh$enhancerstop)
head(df.genes.enh)

barplot(table(df.genes.enh$overlapping), main = "Barplot of enhancer position in relation to his gene")

#Quantile Coverage analysis by overlapping status 
df.genes.enh$width.enh <- df.genes.enh$enhancerstop - df.genes.enh$enhancerstart
quant.overlapping <- aggregate(width.enh ~ overlapping, data=df.genes.enh, FUN = "quantile" ,probs=c(0.25, 0.50,0.75))
quant.overlapping

hist(df.genes.enh$width.enh)

analysis.g_e <- contact.func("geneSymbol",agg.table.enhancers)
head(analysis.g_e$quantitative_asso)

analysis.e_g <- contact.func("enhancer", agg.table.genes)
head(analysis.e_g$quantitative_asso)

#Linear mapping between genes and enhancers
plot(df.genes.enh$TSS, df.genes.enh$enhancerstop)
#It seems that there is a perfectly linear relation between genes and enhancers

#Network analysis from gene-enhancer clusters 
#To DO: - bipartite nework
#       - Summary statistics on cluster behaviour

genes <- matrix(df.genes.enh$geneSymbol, ncol = 1)
enh <- matrix(df.genes.enh$enhancer, ncol = 1)
enh.c1 <- matrix(enhancer.contacts$enhancer1, ncol=1)
enh.c2 <- matrix(enhancer.contacts$enhancer2, ncol=1)

nodes.g_e <- as.data.frame(do.call("rbind",list(genes,enh, enh.c1, enh.c2)))
colnames(nodes.g_e) <- "id"
nodes.g_e$type <- ifelse(nodes.g_e$id%in%genes, "gene","enhancer")
nodes.g_e$col<- ifelse(nodes.g_e$type == "gene", "orange","blue")

nodes.g_e <- unique(nodes.g_e)

links.g_e <- df.genes.enh[,c("geneSymbol", "enhancer")]
colnames(links.g_e) <- c("from","to")
colnames(enh.c) <- c("from","to")

links.g_e <- rbind(links.g_e,enh.c)
links.g_e$weight <- 1
links.g_e

net <- graph_from_data_frame(d=links.g_e,vertices = nodes.g_e,directed = FALSE)

plot(net, vertex.label=NA, vertex.size = 1.5, margin=-.1,asp=.25, vertex.color = nodes.g_e$col, edge.color = "red",main="Clustering based on gene-enhancer contacts")
legend(x=0, y=-1.3, c("gene","enhancer"), pch=21,
        col="#777777", pt.bg=unique(nodes.g_e$col), pt.cex=2, cex=.8, bty="n", ncol=1)


#Most connected elements in cluster: diameter function returns the number of edges
diameter(net, directed = FALSE, weights = NA)

#So there are 5 nodes 
get_diameter(net, directed=FALSE, weights=NA)

#Average number of edges between two nodes
mean_distance(net,directed=FALSE)
#Weak graph density because of the high number of unique gene-enhancer relationship
edge_density(net)
#Probability of two nodes which are connected to the same node are connected with each other
transitivity(net)

#Here we have 389 distincts communities
components(net)
#Average number of components in communities
mean(components(net)$csize)

#Here we define a cluster which is full-connected between all of its elements
comm <-cluster_infomap(net)
plot(comm, net,vertex.label=NA, vertex.size = 2, margin=-.1,asp=.35,main="Gene-Enhancer Community Clusters based on their contacts")

#GO for genes in this cluster
#Because of same genes present in a same enhancer-only cluster have the same regulatory functions
#I filter on the communities which have more 2 elements inside and I make GO annoations


#Epigenomics Roadmap: 25 imputation states-model

epg_states_models <- import("enhancers/E081_25_imputed12marks_dense.bed", format="bed")
length(unique(epg_states_models))
