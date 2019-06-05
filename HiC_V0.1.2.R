#Hi-C_analysis

## Quick Description 
#This project gathers all the chunks of code in relation to Hi-C data analysis.
#Data were retrieved from Rajarajan et al 2018 (1) and available on Psychencode Synapse Platform.
#Our goal is to integrate the non-coding regions of the genome in rare association test variants, with respect to 3D genome contacts. 
#So, first of all, we start by annoting genes and enhancers which are present in contact locus from the Hi-C file. 
#Some statistical analysis were previous made to detect peak enrichment in Hi-C 3D contact data with HiCCUPS(2) software. 
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
#to install  FantomEnhancers.hg19 package, thanks to uncomment the following code line:
#devtools::install_github("charlesjb/fantomenhancers.hg19")
library(FantomEnhancers.hg19)

clean.hic <- read.table(file.choose(),header = TRUE, sep="\t")

#Let's take a look on initial data
#HiC Data had already been processed to determine the significant contacts
head(clean.hic)

#As we're working on schizophrenia and bipolar troubles, we just retrieve data from neurons
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
  mcols(GRange.to.annotate)$chr <- NA
  mcols(GRange.to.annotate)$geneSymbol <- NA
  mcols(GRange.to.annotate)$TSS <- NA
  mcols(GRange.to.annotate)$TES <- NA
  
  #Correspondance between contacts locus and genes which overlaps this latter
  mcols(GRange.to.annotate)[subjectHits(overlaps.index), "chr"] <- seqnames(genes[queryHits(overlaps.index)])
  mcols(GRange.to.annotate)[subjectHits(overlaps.index), "geneSymbol"] <- mcols(genes)[queryHits(overlaps.index), "geneSymbol"]
  mcols(GRange.to.annotate)[subjectHits(overlaps.index), "TSS"] <- mcols(genes)[queryHits(overlaps.index), "TSS"]
  mcols(GRange.to.annotate)[subjectHits(overlaps.index), "TES"] <- mcols(genes)[queryHits(overlaps.index), "TES"]
  
  GRange.annotated <- GRange.to.annotate
  return(GRange.annotated)
}

filter.enhancers <- function(enhancer.GRange, rule=c("one","consensus","total")){
  #This function takes in input a enhancer file from Fantom5
  #A specific filter on level expression is done 
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
  #This function annotates a GRange with enhancer file
  #The two files have to be GRanges 
  
  #Annotations for enhancers
  mcols(enhancers)$enhancerstart <- start(enhancers)
  mcols(enhancers)$enhancerstop <- end(enhancers)
  
  overlaps.index <- findOverlaps(enhancers, GRange.to.annotate)
  
  mcols(GRange.to.annotate)$chr <- NA
  mcols(GRange.to.annotate)$enhancerstart <- NA
  mcols(GRange.to.annotate)$enhancerstop <- NA
  
  mcols(GRange.to.annotate)[subjectHits(overlaps.index), "chr"] <- seqnames(enhancers[queryHits(overlaps.index)])
  mcols(GRange.to.annotate)[subjectHits(overlaps.index),"enhancerstart"] <- mcols(enhancers)[queryHits(overlaps.index),"enhancerstart"]
  mcols(GRange.to.annotate)[subjectHits(overlaps.index),"enhancerstop"] <- mcols(enhancers)[queryHits(overlaps.index), "enhancerstop"]
  
  GRange.annotated <- GRange.to.annotate
  return(GRange.annotated)
  
}

all.contact <- function(annotated.genes, annotated.enhancers, annotated.genes1, annotated.enhancers1){
  #IMPORTANT: The input files have to be indexed for well working of the function
  #The function returns 3 files: two annotated gene-enhancer file and 1 enhancer-enhancer file
  asso.genes <- annotated.genes[names(annotated.enhancers[!is.na(mcols(annotated.enhancers)$enhancerstart)])]
  asso.genes <- asso.genes[!is.na(mcols(asso.genes)$geneSymbol)]
  
  index.non.null.genes <- names(asso.genes)
  
  asso.enhanc <- annotated.enhancers[index.non.null.genes]
  
  asso.genes.1 <- annotated.genes1[names(annotated.enhancers1[!is.na(mcols(annotated.enhancers1)$enhancerstart)])]
  asso.genes.1 <- asso.genes.1[!is.na(mcols(asso.genes.1)$geneSymbol)]
  
  index.non.null.genes.1 <- names(asso.genes.1)
  
  asso.enhanc.1 <- annotated.enhancers1[index.non.null.genes.1]
  
  enhancer.enhancer <- annotated.enhancers[names(annotated.enhancers1[!is.na(mcols(annotated.enhancers1)$enhancerstart)])]
  enhancer.enhancer <- enhancer.enhancer[!is.na(mcols(enhancer.enhancer)$enhancerstart)]
  
  index.e_e <- names(enhancer.enhancer)
  
  enhancer.enhancer.1<- annotated.enhancers1[index.e_e]
  
  tmp.genes <- as.data.frame(mcols(asso.genes))
  tmp.genes.1 <- as.data.frame(mcols(asso.genes.1))
  
  df.genes <- rbind(tmp.genes, tmp.genes.1)
  
  tmp.enhanc <- as.data.frame(mcols(asso.enhanc))
  tmp.enhanc.1 <- as.data.frame(mcols(asso.enhanc.1))
  
  df.enhanc <- rbind(tmp.enhanc, tmp.enhanc.1)
  
  contact.genes.enh <- cbind(df.genes, df.enhanc)
  contact.genes.enh$enhancer <- paste(contact.genes.enh$enhancerstart, contact.genes.enh$enhancerstop)
  
  tmp.ee <- as.data.frame(mcols(enhancer.enhancer))
  tmp.ee.1 <- as.data.frame(mcols(enhancer.enhancer.1))
  
  contact.enh.enh <- cbind(tmp.ee, tmp.ee.1)
  colnames(contact.enh.enh) <- c("chr1","enhancerstart1","enhancerstop1","chr2","enhancerstart2","enhancerstop2")
  
  contact.enh.enh$enhancer1 <- paste(contact.enh.enh$enhancerstart1, contact.enh.enh$enhancerstop1)
  contact.enh.enh$enhancer2 <- paste(contact.enh.enh$enhancerstart2, contact.enh.enh$enhancerstop2)
  
  return(list("gene_enhancer"=contact.genes.enh,"enhancer_enhancer"= contact.enh.enh))
  
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

contact.func <- function(data, col.to.agg){
  #This function takes as inputs a DataFrame and a column to aggregate
  #The function returns a qualitative and quantitative results based on "full-linked" cluster
  #This definition allows that each element of a cluster has to be linked each other
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

create.igraph.matrix <-function(df.genes.enhancers, df.enhancers.contact){
  #Function which converts DataFrame of gene enhancer contacts and DataFrame of enhancer enhancer contacts into nodes and links objects
  #The function returns a list with nodes and links attributes 
  #The results can be directly used in igraph function to create network
  tmp.genes <- matrix(df.genes.enhancers$geneSymbol,ncol=1)
  tmp.enh <- matrix(df.genes.enhancers$enhancer, ncol =1)
  
  enh.c <- df.enhancers.contact[,c("enhancer1","enhancer2")]
  enh.c1 <- matrix(enh.c$enhancer1, ncol=1)
  enh.c2 <- matrix(enh.c$enhancer2, ncol=1)
  
  nodes <- as.data.frame(do.call("rbind",list(tmp.genes,tmp.enh, enh.c1, enh.c2)))
  colnames(nodes) <- "id"
  nodes$type <- ifelse(nodes$id%in%tmp.genes, "gene","enhancer")
  nodes$col<- ifelse(nodes$type == "gene", "orange","blue")
  nodes <- unique(nodes)
  
  
  links <- df.genes.enhancers[,c("geneSymbol", "enhancer")]
  colnames(links) <- c("from","to")
  colnames(enh.c) <- c("from","to")
  
  links <- rbind(links,enh.c)
  links$weight <- 1
  
  list("nodes"=nodes, "links"= links)
}

network.analysis <- function(network){
  #Most connected elements in cluster: diameter function returns the number of edges
  diameter <- diameter(network, directed = FALSE, weights = NA)
  #The most-connected sub-network of graph
  g.diameter <- get_diameter(network, directed=FALSE, weights=NA)
  #Average number of edges between two nodes
  md <- mean_distance(network,directed=FALSE)
  #Graph Density: Probability of two nodes which are linked on all possible connections in graph
  d <- edge_density(network)
  #Probability of two nodes which are connected to the same node are connected with each other
  t <- transitivity(network)
  return(list("diameter"=diameter,"get_diameter"=g.diameter, "mean_distance"=md, "density"=d, "transitivity"=t))
}

sum(mapply(significant.contact, clean.hic.neu$fdrBL, clean.hic.neu$fdrDonut, clean.hic.neu$fdrH, clean.hic.neu$fdrV, .10))
#100% of pixels provided from the file are enriched. with a 10% FDR thresold, this value may changed with a 5% threshold or 1% threshold.
#We will do the same analysis with different threshold's values
sum(mapply(significant.contact, clean.hic.neu$fdrBL, clean.hic.neu$fdrDonut, clean.hic.neu$fdrH, clean.hic.neu$fdrV, .05))
sum(mapply(significant.contact, clean.hic.neu$fdrBL, clean.hic.neu$fdrDonut, clean.hic.neu$fdrH, clean.hic.neu$fdrV, .01))
#The rule is extracted from Rao and al 2014.

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
all.equal(contacts.locus$X,contacts.locus$Y)

#Unique coverage from contact locus for the both of files
paste("Total genome coverage : ",sum(width(reduce(contacts.locus$X))), "bp")
paste("Total genome coverage : ",sum(width(reduce(contacts.locus$Y))), "bp")

#Hg19 genes used for genes annotations
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

#Unified genes overlapping
unified.genes <- intersect(genes.overlaps.X, genes.overlaps.Y)
sum(width(reduce(unified.genes)))
mean(width(reduce(unified.genes)))

quantiles.overlapping.genes <- quantile(width(unified.genes), probs=c(.25,.50,.75))
#25%    50%    75% 
#111891 242716 462379 

hist(width(unified.genes), main="Distribution of overlapping contact-locus gene width")
#It seems to have a most important part of overlapping genes which are little genes

#density for overlapping genes all along chromosomes 
kp.genes <- plotKaryotype(genome="hg19")
kpPlotDensity(kp.genes, data=unified.genes)

#Have enhancers for each genes, data are extracted from Fantom5
#Data are available on Github repo
load(file.choose())
#Let's take a look on data 
head(neurons_enhancers)

genes.annotations.X <- genes.annotations(genes.hg19,contacts.locus$X)
genes.annotations.Y <- genes.annotations(genes.hg19,contacts.locus$Y1)

#Here we use the lowest filter size
filter.enhancers.cons <- filter.enhancers(neurons_enhancers, "consensus")
filter.enhancers.one <- filter.enhancers(neurons_enhancers, "one")

enhancers.annotations.Y <- enhancers.annotations(filter.enhancers.one,contacts.locus$Y)
enhancers.annotations.X <- enhancers.annotations(filter.enhancers.one,contacts.locus$X1)

#Because of the order of bins is preserved  we index the 2 type files for retrieving the genes associated with enhancers
names(genes.annotations.X) <- 1:length(genes.annotations.X)
names(enhancers.annotations.Y) <- 1:length(enhancers.annotations.Y)
names(enhancers.annotations.X) <- 1:length(enhancers.annotations.X)
names(genes.annotations.Y) <- 1:length(genes.annotations.Y)

genes.enhancers.contact <- all.contact(genes.annotations.X,enhancers.annotations.Y,genes.annotations.Y,enhancers.annotations.X)$gene_enhancer
enhancers.enhancers.contact <- all.contact(genes.annotations.X,enhancers.annotations.Y,genes.annotations.Y,enhancers.annotations.X)$enhancer_enhancer

enhancers.enhancers.contact$med1 <- mapply(find.medoid,enhancers.enhancers.contact$enhancerstart1, enhancers.enhancers.contact$enhancerstop1)
enhancers.enhancers.contact$med2 <- mapply(find.medoid,enhancers.enhancers.contact$enhancerstart2,enhancers.enhancers.contact$enhancerstop2)

#Distribution Quantiles of distance between two connected enhancers
hist(abs(enhancers.enhancers.contact$med1 - enhancers.enhancers.contact$med2),main = "Distribution of proximity (bp) between two connected enhancers")
#It seems that two connected enhancers are mostly closed
quantile(abs(enhancers.enhancers.contact$med1 - enhancers.enhancers.contact$med2), probs=c(.25,.5,.75))

all.equal(enhancers.enhancers.contact$enhancer1, enhancers.enhancers.contact$enhancer2)
#[1] "75 string mismatches"

length(unique(enhancers.enhancers.contact$enhancer1))
#74
length(unique(enhancers.enhancers.contact$enhancer2))
#75

enhancers.FANTOM <- GRanges(seqnames = genes.enhancers.contact$chr,ranges=IRanges(start = genes.enhancers.contact$enhancerstart , end = genes.enhancers.contact$enhancerstop),strand="*")
genes.FANTOM <- GRanges(seqnames = genes.enhancers.contact$chr,ranges=IRanges(start = genes.enhancers.contact$TSS , end = genes.enhancers.contact$TES), strand="*")

#Linear mapping between genes and enhancers, here we filter by chromosomes
for(chr in unique(genes.enhancers.contact$chr)){
  tmp.chr <- genes.enhancers.contact[genes.enhancers.contact$chr==chr,]
  plot(tmp.chr$TSS, tmp.chr$enhancerstop, main=paste("Linear Mapping for ",chr), xlab="TSS", ylab="enhancerstop")
}

#Density superpostion between genes and enhancers from Fantom5
kp <- plotKaryotype()
kpPlotDensity(kp, enhancers.FANTOM, col="#FF00007F")
kpPlotDensity(kp, genes.FANTOM,col="#0000FF7F")

#Non redundant Genome coverage by enhancers
sum(width(reduce(enhancers.FANTOM)))
sum(width(reduce(genes.FANTOM)))

#Two network definitions are now possible to determine gene-enhancer connections. 
#The first one is that we call "full connected network", by definition in this latter all cluster elements are connected each others
#The second one is that we call "at least one connected network", by definition in this network each element is linked by at least one other element of network
#We choose to present the second network version
nodes.FANTOM <- create.igraph.matrix(genes.enhancers.contact, enhancers.enhancers.contact)$nodes
links.FANTOM <- create.igraph.matrix(genes.enhancers.contact, enhancers.enhancers.contact)$links

network.FANTOM <- graph_from_data_frame(d= links.FANTOM,vertices=nodes.FANTOM,directed = FALSE)
plot(network.FANTOM, vertex.label=NA, vertex.size = 1.5, margin=-.1,asp=.25, vertex.color = nodes.FANTOM$col, edge.color = "red",main="Clustering based on Gene-Enhancer contacts")
legend(x=0, y=-1.3, c("gene","enhancer"), pch=21,
       col="#777777", pt.bg=unique(nodes.FANTOM$col), pt.cex=2, cex=.8, bty="n", ncol=1)

network.analysis(network.FANTOM)

compo.FANTOM <- components(network.FANTOM)
#Number of unique clusters
compo.FANTOM$no
#[1] 389
mean(compo.FANTOM$csize)
#2.329049
table(compo.FANTOM$csize)
# 2   3   4   5   6 
#302  56  24   4   3 

#"Full-connected cluster" definition here: 
#The code is not executed

#agg.table.enhancers <- aggregate(geneSymbol~enhancer, data=genes.enhancers.contact, unique, na.rm=TRUE)
#agg.table.genes <- aggregate(enhancer~geneSymbol, data=genes.enhancers.contact, unique, na.rm=TRUE)

#table(sapply(agg.table.enhancers$geneSymbol, length))
#table(sapply(agg.table.genes$enhancer, length))

#barplot(table(sapply(agg.table.enhancers$geneSymbol, length)), main = "Distribution of number of genes associated by enhancer", xlab = "Number of genes associated", ylab = "Frequency")
#barplot(table(sapply(agg.table.genes$enhancer, length)), main = "Distribution of number of enhancers associated by genes", xlab = "Number of enhancers associated", ylab = "Frequency")

#genes.enhancers.contact$overlapping <- mapply(test.overlapping, genes.enhancers.contact$TSS, genes.enhancers.contact$TES, genes.enhancers.contact$enhancerstart, genes.enhancers.contact$enhancerstop)

#barplot(table(genes.enhancers.contact$overlapping), main = "Barplot of enhancer position in relation to his gene")

#Quantile Coverage analysis by overlapping status 
#genes.enhancers.contact$width.enh <- genes.enhancers.contact$enhancerstop - genes.enhancers.contact$enhancerstart
#quant.overlapping <- aggregate(width.enh ~ overlapping, data=genes.enhancers.contact, FUN = "quantile" ,probs=c(0.25, 0.50,0.75))
#quant.overlapping

#analysis.g_e <- contact.func("geneSymbol",agg.table.enhancers)
#head(analysis.g_e$quantitative_asso)

#analysis.e_g <- contact.func("enhancer", agg.table.genes)
#head(analysis.e_g$quantitative_asso)

#Because of lack of complexity with FANTOM annotations, we propose the same methodology as above but with
#Epigenomic Roadmap annotations.
#We choose this kind of annotations because Won et al 2016 used this methodology for annotating their Hi-C contacts files
#WE focus on Epigenomics Roadmap: 15 state-model
#The method of functionnal annotations based on modification of histone patterns is presented in Epigenomics Roadmap Consortium, 2015
#Data are available on Github repo

#File from https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/

enhancers.EGRM <- import(file.choose(), format="bed")
head(enhancers.EGRM)

#Distribution of different states present in the file 
distribution.states <-table(mcols(enhancers.EGRM)$name)
barplot(sort((distribution.states)/sum(distribution.states))*100,las =2)

#We focus our analysis only on enhancers
enhancers.EGRM.filter <- mcols(enhancers.EGRM)$name %in% c("6_EnhG","7_Enh")
enhancers.EGRM <- enhancers.EGRM[enhancers.EGRM.filter]

enhancers.EGRM.Y <- enhancers.annotations(enhancers.EGRM,contacts.locus$Y)
enhancers.EGRM.X <- enhancers.annotations(enhancers.EGRM,contacts.locus$X1)

names(enhancers.EGRM.Y) <- 1:length(enhancers.EGRM.Y)
names(enhancers.EGRM.X) <- 1:length(enhancers.EGRM.X)

genes.enhancers.EGRM <- all.contact(genes.annotations.X,enhancers.EGRM.Y,genes.annotations.Y,enhancers.EGRM.X)$gene_enhancer
enhancers.enhancers.EGRM <- all.contact(genes.annotations.X,enhancers.EGRM.Y,genes.annotations.Y,enhancers.EGRM.X)$enhancer_enhancer

enhancers.EGRM <- GRanges(seqnames = genes.enhancers.EGRM$chr,ranges=IRanges(start = genes.enhancers.EGRM$enhancerstart , end = genes.enhancers.EGRM$enhancerstop),strand="*")
genes.EGRM <- GRanges(seqnames = genes.enhancers.EGRM$chr,ranges=IRanges(start = genes.enhancers.EGRM$TSS , end = genes.enhancers.EGRM$TES), strand="*")
#Non redundant Genome coverage by enhancers and genes
sum(width(reduce(enhancers.EGRM)))
#[1] 1479800
sum(width(reduce(genes.EGRM)))
#[1] 188602968

#Linear mapping between genes and enhancers, here we filter by chromosomes
for(chr in unique(genes.enhancers.EGRM$chr)){
  pdf(paste0("/home/nash/Documents/Psychencode/output/linear_mapping/linear_mapping_",chr,".pdf"))
  tmp.chr <- genes.enhancers.EGRM[genes.enhancers.EGRM$chr==chr,]
  plot(tmp.chr$TSS, tmp.chr$enhancerstop, main=paste("Linear Mapping for ",chr), xlab="TSS", ylab="enhancerstop")
  dev.off()
}

#Density superpostion between genes and enhancers from Epigenomic Roadmap
col1 <- rgb(0,0,255, max=255, alpha=50)
col2 <- rgb(0,255,0, max=255, alpha=50)
kp <- plotKaryotype()
kpPlotDensity(kp, enhancers.EGRM, col=col1)
kpPlotDensity(kp, genes.EGRM,col=col2)

#Linear mapping on 500kbp window to determine proximity between genes based on enhancerstop
ex.chr12 <- genes.enhancers.EGRM[genes.enhancers.EGRM$chr=="chr12",]

vec.X <- ex.chr12[ex.chr12$TSS <= 5000000,]$TSS
vec.Y <- ex.chr12[ex.chr12$TSS <= 5000000,]$enhancerstop
#We plot only on top diagonal
lm <- lm(vec.Y~vec.X)
for(i in 1:length(lm$residuals)){
  if(lm$residuals[i] < 0){
    vec.Y[i] = vec.Y[i] +2*abs(lm$residuals[i])
  }
}
plot(vec.X, vec.Y,
     xlab="TSS",ylab="enhancerstop", main="Linear mapping on 5Mbp window to determine proximity between genes")
text(vec.X, vec.Y, ex.chr12[ex.chr12$TSS<=5000000,]$geneSymbol,cex=0.6,col="red")


#Create useful data to igraph analysis, because of the high number of data and to improve the graph quality we perform the vizualisation by chromosomes
for(chr in unique(genes.enhancers.EGRM$chr)){
  pdf(paste0("/home/nash/Documents/Psychencode/output/networks/cluster_network",chr,".pdf"))
  tmp.genes.enhancers <- genes.enhancers.EGRM[genes.enhancers.EGRM$chr==chr,]
  tmp.enhancers.enhancers <- enhancers.enhancers.EGRM[enhancers.enhancers.EGRM$chr1 == chr,]
  tmp.nodes <- create.igraph.matrix(tmp.genes.enhancers, tmp.enhancers.enhancers)$nodes
  tmp.links <- create.igraph.matrix(tmp.genes.enhancers, tmp.enhancers.enhancers)$links
  network <- graph_from_data_frame(d=tmp.links,vertices = tmp.nodes,directed = F)
  plot(network, vertex.label=NA, vertex.size = 2, margin=-.1,asp=.35, vertex.color = tmp.nodes$col, edge.color = "green",main=paste("Clustering based on Gene-Enhancer contacts: ", chr))
  legend(x=0, y=-1.3, c("gene","enhancer"), pch=21,
         col="#777777", pt.bg=unique(tmp.nodes$col), pt.cex=2, cex=.8, bty="n", ncol=1)
  dev.off()
}

nodes.EGRM <- create.igraph.matrix(genes.enhancers.EGRM, enhancers.enhancers.EGRM)$nodes
links.EGRM <- create.igraph.matrix(genes.enhancers.EGRM, enhancers.enhancers.EGRM)$links
network.EGRM <- graph_from_data_frame(d=links.EGRM,vertices = nodes.EGRM,directed = F)

network.analysis(network.EGRM)
compo.epg <- components(network.EGRM)
compo.epg$no

mean(compo.epg$csize)
#[1] 2.948108
table(compo.epg$csize)
# 2   3   4   5   6   7   8  10  13 
#433 236 186  38  17  10   3   1   1 

#agg.table.enhancers.epg <- aggregate(geneSymbol~enhancer, data=genes.enhancers.EGRM, unique, na.rm=TRUE)
#agg.table.genes.epg <- aggregate(enhancer~geneSymbol, data=genes.enhancers.EGRM, unique, na.rm=TRUE)

#table(sapply(agg.table.enhancers.epg$geneSymbol, length))
#table(sapply(agg.table.genes.epg$enhancer, length))

#barplot(table(sapply(agg.table.enhancers.epg$geneSymbol, length)), main = "Distribution of number of genes associated by enhancer", xlab = "Number of genes associated", ylab = "Frequency")
#barplot(table(sapply(agg.table.genes.epg$enhancer, length)), main = "Distribution of number of enhancers associated by genes", xlab = "Number of enhancers associated", ylab = "Frequency")

#genes.enhancers.EGRM$overlapping <- mapply(test.overlapping, genes.enhancers.EGRM$TSS, genes.enhancers.EGRM$TES, genes.enhancers.EGRM$enhancerstart, genes.enhancers.EGRM$enhancerstop)

#barplot(table(genes.enhancers.EGRM$overlapping), main = "Barplot of enhancer position in relation to his gene")

#Quantile Coverage analysis by overlapping status 
#genes.enhancers.EGRM$width.enh <- genes.enhancers.EGRM$enhancerstop - genes.enhancers.EGRM$enhancerstart
#quant.overlapping <- aggregate(width.enh ~ overlapping, data=dgenes.enhancers.EGRM FUN = "quantile" ,probs=c(0.25, 0.50,0.75))
#quant.overlapping

#analysis.g_e.epg <- contact.func("geneSymbol",agg.table.enhancers.epg)
#head(analysis.g_e.epg$quantitative_asso)

#analysis.e_g.epg <- contact.func("enhancer", agg.table.genes.epg)
#head(analysis.e_g.epg$quantitative_asso)

n.genes <- length(unique(genes.enhancers.EGRM$geneSymbol))
n.enhancers <- length(unique(genes.enhancers.EGRM$enhancer))

one.g.one.e <- table(compo.epg$csize)[1][[1]]
one.g.one.e / n.genes ; one.g.one.e/n.enhancers
#[1] 0.3775065 :Proportion of genes which are single linked with an enhancer
#[1] 0.3497577 : Proportion of enhancers which are single linked with a gene

#Statistics per cluster
gene.enhancer.clusters = tapply(names(compo.epg$membership),compo.epg$membership,function(vec,genes) table(vec%in%genes),genes=genes.enhancers.EGRM$geneSymbol)
gene.enhancer.clusters.mat = matrix(unlist(gene.enhancer.clusters),length(unlist(gene.enhancer.clusters))/2,2,byrow = T)

# Analysis by gene
tapply(gene.enhancer.clusters.mat[,1],gene.enhancer.clusters.mat[,2],summary)
# Analysis by enhancer
tapply(gene.enhancer.clusters.mat[,2],gene.enhancer.clusters.mat[,1],summary)

# Proportion of genes with a single enhancer
sum(gene.enhancer.clusters.mat[gene.enhancer.clusters.mat[,1]==1,2])/n.genes
#[1] 0.3591979
# Proportion of enhancers linked to a single gene
sum(gene.enhancer.clusters.mat[gene.enhancer.clusters.mat[,2]==1,1])/n.enhancers
#[1] 0.7714055

# Total cluster length
longueur.cluster  = tapply(pmax(genes.enhancers.EGRM$enhancerstop,genes.enhancers.EGRM$TES),
                           compo.epg$membership[genes.enhancers.EGRM$geneSymbol],max) - tapply(pmin(genes.enhancers.EGRM$enhancerstart,genes.enhancers.EGRM$TSS),
                                                                                               compo.epg$membership[genes.enhancers.EGRM$geneSymbol],min)
summary(longueur.cluster)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#42137  209249  370407  534965  657837 5789608 

#Extraction of the larger network and summary analysis
larger<- which.max(table(compo.epg$membership))
larger.node <- nodes.epg[compo.epg$membership == larger[[1]],]
larger.links <- links.EGRM[links.EGRM$from%in%larger.node$id | links.EGRM$to%in% larger.node$id,]
larger.graph = graph_from_data_frame(d=larger.links,directed=F,vertices=larger.node)
V(larger.graph)$label.cex = 0.5
V(larger.graph)$label.color = larger.node$col
plot(larger.graph,vertex.size = 1.5, margin=-.1,asp=.25,vertex.color = larger.node$col, edge.color = "red",main="Larger Network")