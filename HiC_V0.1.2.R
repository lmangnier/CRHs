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

#Some functions to facilitate analysis....
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
  
  GRanges.X <- GRanges(seqnames = x$chr1, ranges = IRanges(x$x1, x$x2), contacts = x$observed)
  GRanges.X1 <- GRanges(seqnames = x$chr1, ranges = IRanges(x$x1, x$x2),contacts = x$observed)
  
  GRanges.Y <- GRanges(seqnames = y$chr2, ranges = IRanges(y$y1, y$y2), contacts = y$observed)
  GRanges.Y1 <- GRanges(seqnames = y$chr2, ranges = IRanges(y$y1, y$y2), contacts = y$observed)
  
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
  #Because we focus our work on epigenetic regulation, we do not include the coregulation between 2 genes
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
  colnames(contact.enh.enh) <- c("contacts1","chr1","enhancerstart1","enhancerstop1","contacts2","chr2","enhancerstart2","enhancerstop2")
  
  contact.enh.enh$enhancer1 <- paste(contact.enh.enh$enhancerstart1, contact.enh.enh$enhancerstop1)
  contact.enh.enh$enhancer2 <- paste(contact.enh.enh$enhancerstart2, contact.enh.enh$enhancerstop2)
  
  return(list("gene_enhancer" = contact.genes.enh,"enhancer_enhancer" = contact.enh.enh))
  
}

create.igraph.matrix <-function(df.genes.enhancers, df.enhancers.contact){
  #Function which converts DataFrame of gene enhancer contacts and DataFrame of enhancer enhancer contacts into nodes and links objects
  #The function returns a list with nodes and links attributes 
  #The results can be directly used in igraph function to create network
  tmp.genes <- matrix(df.genes.enhancers$geneSymbol,ncol=1)
  tmp.enh <- matrix(df.genes.enhancers$enhancer, ncol =1)
  
  enh.c <- df.enhancers.contact[,c("enhancer1","enhancer2","contacts2")]
  enh.c1 <- matrix(enh.c$enhancer1, ncol=1)
  enh.c2 <- matrix(enh.c$enhancer2, ncol=1)
  
  nodes <- as.data.frame(do.call("rbind",list(tmp.genes,tmp.enh, enh.c1, enh.c2)))
  colnames(nodes) <- "id"
  nodes$type <- ifelse(nodes$id%in%tmp.genes, "gene","enhancer")
  nodes$col<- ifelse(nodes$type == "gene", "orange","blue")
  nodes <- unique(nodes)
  
  
  links <- df.genes.enhancers[,c("geneSymbol", "enhancer","contacts")]
  colnames(links) <- c("from","to","weight")
  colnames(enh.c) <- c("from","to","weight")
  
  links <- rbind(links,enh.c)
  
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

Numextract <- function(string){
  #Function of seqnames formating column for GRanges
  unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
}

########################################################### END of FUNCTIONS ############################################

########################################################### Start of CODE ###############################################
clean.hic <- read.table(file.choose(),header = TRUE, sep="\t")
#We retrieve data from neurons and Astro in two  differents dataframes
clean.hic.neu <- clean.hic[clean.hic$Cell == "Neu",]
clean.hic.astro <- clean.hic[clean.hic$Cell == "Astro"]

contacts_DNA_DNA <- clean.hic.neu[,c("chr1","x1", "y2", "observed")]
names(contacts_DNA_DNA) <- c("chr", "start", "end", "observed")

GRanges.DNA.DNA <- makeGRangesFromDataFrame(contacts_DNA_DNA, keep.extra.columns = TRUE)

#Let's take a look on initial data
#HiC Data had already been processed to determine the significant contacts
head(clean.hic.neu)
dim(clean.hic.neu)
#[1] 2471   25
sum(mapply(significant.contact, clean.hic.neu$fdrBL, clean.hic.neu$fdrDonut, clean.hic.neu$fdrH, clean.hic.neu$fdrV, .10))
#[1] 2471

#Analysis of FDR's correlation
correlation.function(clean.hic.neu, c("fdrBL","fdrDonut","fdrH","fdrV"))

#Creation of one column for possible each specific rule 
#Here a boolean condition on each FDR
clean.hic.neu$FDR_0.1 <- mapply(significant.contact, clean.hic.neu$fdrBL, clean.hic.neu$fdrDonut, clean.hic.neu$fdrH, clean.hic.neu$fdrV, .10)
head(clean.hic.neu)

#Are there trans contacts ?
table(clean.hic.neu$chr1, clean.hic.neu$chr2)

#Now let's start by gene annotations
#I create two datasets come from .hic file
x <- clean.hic.neu[,c("chr1","x1","x2", "observed")]
y <- clean.hic.neu[,c("chr2","y1","y2", "observed")]

contacts.locus <- create.reverse.GRanges(x,y)
head(contacts.locus$X);head(contacts.locus$Y)

mcols(contacts.locus$X)$medoidX <- start(contacts.locus$X) + (end(contacts.locus$X) - start(contacts.locus$X))/2
mcols(contacts.locus$X)$medoidY <- start(contacts.locus$Y) + (end(contacts.locus$Y) - start(contacts.locus$Y))/2
mcols(contacts.locus$X)$dist.medtomed <- mcols(contacts.locus$X)$medoidY - mcols(contacts.locus$X)$medoidX

distance_between_medoid <- mcols(contacts.locus$X)$medoidY - mcols(contacts.locus$X)$medoidX
summary(distance_between_medoid)

median.contacts.by.dist <-aggregate(contacts~dist.medtomed,data=as(contacts.locus$X, "data.frame"), mean)

#Correlation between medoid distance and number of contact
cor(median.contacts.by.dist$dist.medtomed, median.contacts.by.dist$contacts, method="pearson")
#[1] -0.1190844

#Histogramme of mean number of contacts by distance between medoid of contact loci
hist(median.contacts.by.dist$contacts)
#Goodness of fit :
shapiro.test(median.contacts.by.dist$contacts)
#p-value = 1.992e-08 : We reject Ho, data do not follow normal distribution

mcols(contacts.locus$X)$medoidX <- NULL
mcols(contacts.locus$X)$medoidY <- NULL
mcols(contacts.locus$X)$dist.medtomed <- NULL

#Unique coverage from contact locus for the both of files
sum(width(reduce(union(contacts.locus$X,contacts.locus$Y,ignore.strand=TRUE))))

#Hg19 genes used for genes annotations
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes.hg19 <- genes(txdb)

symbol <- select(org.Hs.eg.db,keys = genes.hg19$gene_id, columns = c("SYMBOL"), keytype = "ENTREZID")
genes.hg19$geneSymbol <- symbol$SYMBOL 
#We're taking the TSS and TES from genes for the post analysis
genes.hg19$TSS <- start(genes.hg19)
genes.hg19$TES <- end(genes.hg19)

#Number of gene overlaps on unique contact regions
common.region.X.Y <- union(contacts.locus$X,contacts.locus$Y,ignore.strand=TRUE)
counts.overlaps.X <- countOverlaps(genes.hg19, common.region.X.Y)
table(counts.overlaps.X!=0)
#FALSE  TRUE 
#20401  2655

sum(width(reduce(subsetByOverlaps(genes.hg19,common.region.X.Y))))
#374678512
summary(width(reduce(subsetByOverlaps(genes.hg19,common.region.X.Y))))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#54   22371   67562  147686  172588 5128602 
hist(width(reduce(subsetByOverlaps(genes.hg19,common.region.X.Y))), main="Distribution of overlapping contact-locus gene width")
#It seems to have a most important part of overlapping genes which are little genes

#density for overlapping genes all along chromosomes 
kp.genes <- plotKaryotype(genome="hg19")
kpPlotDensity(kp.genes, data=common.region.X.Y)

genes.annotations.X <- genes.annotations(genes.hg19,contacts.locus$X)
genes.annotations.Y <- genes.annotations(genes.hg19,contacts.locus$Y1)

#THIS STEP IS CRUCIAL!!!!!!
#Because of the order of bins is preserved  we index the 2 type files for retrieving the genes associated with enhancers
names(genes.annotations.X) <- 1:length(genes.annotations.X)
names(genes.annotations.Y) <- 1:length(genes.annotations.Y)

#Because of lack of complexity with FANTOM annotations, we propose the same methodology as above but with
#Epigenomic Roadmap annotations.
#We choose this kind of annotations because Won et al 2016 used this methodology for annotating their Hi-C contacts files
#WE focus on Epigenomics Roadmap: 15 state-model
#The method of functionnal annotations based on modification of histone patterns is presented in Epigenomics Roadmap Consortium, 2015
#Data are available on Github repo

#File from https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/

enhancers.EGRM <- import(file.choose(), format="bed")
head(enhancers.EGRM)

#Distribution of different states overlapped Hi-C contact loci present in the file 
overlapping.states <- subsetByOverlaps(enhancers.EGRM, union(contacts.locus$X, contacts.locus$Y))
sort(table(mcols(overlapping.states)$name),decreasing = TRUE)

#We focus our analysis only on enhancers
enhancers.EGRM.filter <- mcols(enhancers.EGRM)$name %in% c("6_EnhG","7_Enh","12_EnhBiv")
enhancers.EGRM <- enhancers.EGRM[enhancers.EGRM.filter]

enhancers.EGRM.Y <- enhancers.annotations(enhancers.EGRM,contacts.locus$Y)
enhancers.EGRM.X <- enhancers.annotations(enhancers.EGRM,contacts.locus$X1)

length(union(enhancers.EGRM.X, enhancers.EGRM.Y))

names(enhancers.EGRM.Y) <- 1:length(enhancers.EGRM.Y)
names(enhancers.EGRM.X) <- 1:length(enhancers.EGRM.X)

genes.enhancers.EGRM <- all.contact(genes.annotations.X,enhancers.EGRM.Y,genes.annotations.Y,enhancers.EGRM.X)$gene_enhancer
enhancers.enhancers.EGRM <- all.contact(genes.annotations.X,enhancers.EGRM.Y,genes.annotations.Y,enhancers.EGRM.X)$enhancer_enhancer

enhancers.EGRM <- GRanges(seqnames = genes.enhancers.EGRM$chr,ranges=IRanges(start = genes.enhancers.EGRM$enhancerstart , end = genes.enhancers.EGRM$enhancerstop),strand="*")
genes.EGRM <- GRanges(seqnames = genes.enhancers.EGRM$chr,ranges=IRanges(start = genes.enhancers.EGRM$TSS , end = genes.enhancers.EGRM$TES), strand="*")

#Non redundant Genome coverage by enhancers and genes

sum(width(reduce(enhancers.EGRM)))
#[1] 1666800
summary(width(reduce(enhancers.EGRM)))

sum(width(reduce(genes.EGRM)))
#[1] 214500987
summary(width(reduce(genes.EGRM)))

distance_gene_enhancer <- abs((start(enhancers.EGRM) + (end(enhancers.EGRM)-start(enhancers.EGRM))/2) - (start(genes.EGRM) + (end(genes.EGRM)-start(genes.EGRM))/2))
summary(distance_gene_enhancer)

#Density superposition between genes and enhancers from Epigenomic Roadmap
col1 <- rgb(0,0,255, max=255, alpha=50)
col2 <- rgb(0,255,0, max=255, alpha=50)
kp <- plotKaryotype()
kpPlotDensity(kp, enhancers.EGRM, col=col1)
kpPlotDensity(kp, genes.EGRM,col=col2)


#Create useful data to igraph analysis, because of the high number of data and to improve the graph quality we perform the vizualisation by chromosomes
for(chr in unique(genes.enhancers.EGRM$chr)){
  #pdf(paste0("/home/nash/Documents/Psychencode/output/networks/cluster_network",chr,".pdf"))
  tmp.genes.enhancers <- genes.enhancers.EGRM[genes.enhancers.EGRM$chr==chr,]
  tmp.enhancers.enhancers <- enhancers.enhancers.EGRM[enhancers.enhancers.EGRM$chr1 == chr,]
  tmp.nodes <- create.igraph.matrix(tmp.genes.enhancers, tmp.enhancers.enhancers)$nodes
  tmp.links <- create.igraph.matrix(tmp.genes.enhancers, tmp.enhancers.enhancers)$links
  network <- graph_from_data_frame(d=tmp.links,vertices = tmp.nodes,directed = F)
  
  #Weighting of network based on contact number between elements
  E(network)$width <- 1 + E(network)$weight/12
  
  plot(network, vertex.label=NA, vertex.size = 2, margin=-.1,asp=.35, vertex.color = tmp.nodes$col, edge.color = "green",main=paste("Clustering based on Gene-Enhancer contacts: ", chr))
  legend(x=0, y=-1.3, c("gene","enhancer"), pch=21,
         col="#777777", pt.bg=unique(tmp.nodes$col), pt.cex=2, cex=.8, bty="n", ncol=1)
  #dev.off()
}

nodes.EGRM <- create.igraph.matrix(genes.enhancers.EGRM, enhancers.enhancers.EGRM)$nodes
links.EGRM <- create.igraph.matrix(genes.enhancers.EGRM, enhancers.enhancers.EGRM)$links
network.EGRM <- graph_from_data_frame(d=links.EGRM,vertices = nodes.EGRM,directed = F)

network.analysis(network.EGRM)
compo.epg <- components(network.EGRM)

tmp.links <- links.EGRM[,]

for(c in compo.epg$membership){
  tmp.nms <- names(compo.epg$membership[compo.epg$membership==c])
  tmp.links[tmp.links$from%in%tmp.nms | tmp.links$to%in%tmp.nms,"n.cluster"] = c
}

mean.contacts.clusters <- aggregate(weight~n.cluster, data=tmp.links, mean)
hist(mean.contacts.clusters$weight)

summary(mean.contacts.clusters$weight)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#14.00   34.00   42.00   44.45   52.00  108.00 

#Linear plots between genes and enhancers, here we filter by chromosomes
for(chr in unique(genes.enhancers.EGRM$chr)){
  pdf(paste0("linear_plots/linear_plot_",chr,".pdf"))
  tmp.chr <- genes.enhancers.EGRM[genes.enhancers.EGRM$chr==chr,]
  #We plot only above top diagonal
  vec.X = ifelse(tmp.chr$TSS<=tmp.chr$enhancerstop,tmp.chr$TSS, tmp.chr$enhancerstop)
  vec.Y = ifelse(tmp.chr$TSS>tmp.chr$enhancerstop,tmp.chr$TSS, tmp.chr$enhancerstop)
  plot(vec.X, vec.Y, main=paste("Gene-enhancer proximity for ",chr), xlab="TSS", ylab="enhancerstop",col=compo.epg$membership[tmp.chr$geneSymbol])
  abline(0,1)
  dev.off()
}

mean(compo.epg$csize)
#[1] 3.046287
table(compo.epg$csize)
#2   3   4   5   6   7   8   9  10  13 
#449 269 228  50  17  12   8   1   2   1 

n.genes <- length(unique(genes.enhancers.EGRM$geneSymbol))
n.enhancers <- length(unique(genes.enhancers.EGRM$enhancer))

one.g.one.e <- table(compo.epg$csize)[1][[1]]
one.g.one.e / n.genes ; one.g.one.e/n.enhancers
#[1] 0.344589 :Proportion of genes which are single linked with an enhancer
#[1] 0.3159747 : Proportion of enhancers which are single linked with a gene

#Statistics per cluster
gene.enhancer.clusters = tapply(names(compo.epg$membership),compo.epg$membership,function(vec,genes) table(vec%in%genes),genes=genes.enhancers.EGRM$geneSymbol)
gene.enhancer.clusters.mat = matrix(unlist(gene.enhancer.clusters),length(unlist(gene.enhancer.clusters))/2,2,byrow = T)

# Analysis by gene
tapply(gene.enhancer.clusters.mat[,1],gene.enhancer.clusters.mat[,2],summary)
# Analysis by enhancer
tapply(gene.enhancer.clusters.mat[,2],gene.enhancer.clusters.mat[,1],summary)

# Proportion of genes with a single enhancer
sum(gene.enhancer.clusters.mat[gene.enhancer.clusters.mat[,1]==1,2])/n.genes
#[1] 0.3169609

# Proportion of enhancers linked to a single gene
sum(gene.enhancer.clusters.mat[gene.enhancer.clusters.mat[,2]==1,1])/n.enhancers
#[1] 0.7714055

# Total cluster length
end.cluster = tapply(pmax(genes.enhancers.EGRM$enhancerstop,genes.enhancers.EGRM$TES),
                       compo.epg$membership[genes.enhancers.EGRM$geneSymbol],max)
start.cluster = tapply(pmin(genes.enhancers.EGRM$enhancerstart,genes.enhancers.EGRM$TSS),
                       compo.epg$membership[genes.enhancers.EGRM$geneSymbol],min)

longueur.cluster  =  end.cluster - start.cluster

summary(longueur.cluster)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#42137  218445  384409  546059  694835 5789608 


#Link between size of cluster and number of contacts inside cluster
HiC.intensity.contacts <- function(links,df.genes.enhancers, df.enhancers.enhancers){
  i <- 1
  #Creation and formatting of output dataframe
  df.contacts.cluster.analysis <- data.frame(matrix(ncol=6))
  colnames(df.contacts.cluster.analysis) <- c("id.cluster", "n.contacts", "n.elements","n.genes","n.enhancers","span.cluster")
  
  for(id in unique(links$n.cluster)){
    
    df.contacts.cluster.analysis[i,"id.cluster"] = id
    df.contacts.cluster.analysis[i, "n.contacts"] = mean(tmp.links[tmp.links$n.cluster == id,]$weight)
    
    elements <- unique(c(links[links$n.cluster==id,"from"],links[tmp.links$n.cluster==id,"to"]))
    df.contacts.cluster.analysis[i,"n.elements"] <- length(elements)
    
    tmp.gene.enh <- df.genes.enhancers[genes.enhancers.EGRM$geneSymbol%in%elements | genes.enhancers.EGRM$enhancer%in%elements,]
    tmp.enh.enh <- df.enhancers.enhancers[enhancers.enhancers.EGRM$enhancer1%in%elements | enhancers.enhancers.EGRM$enhancer2%in%elements,]
    
    if(nrow(tmp.enh.enh) == 0){
      
      df.contacts.cluster.analysis[i, "n.genes"] <- length(unique(tmp.gene.enh$geneSymbol))
      df.contacts.cluster.analysis[i, "n.enhancers"] <- length(unique(tmp.gene.enh$enhancer))
      
      max.max <- max(tmp.gene.enh[,c("TSS","TES", "enhancerstart", "enhancerstop")])
      min.min <- min(tmp.gene.enh[,c("TSS","TES", "enhancerstart", "enhancerstop")])
      
      df.contacts.cluster.analysis[i,"span.cluster"] <-  max.max - min.min
    }
    
    else if(nrow(tmp.gene.enh) == 0){
      
      df.contacts.cluster.analysis[i, "n.genes"] <- 0
      df.contacts.cluster.analysis[i, "n.enhancers"] <- length(unique(c(tmp.gene.enh$enhancer,tmp.enh.enh$enhancer1, tmp.enh.enh$enhancer2)))
      
      max.max <- max(tmp.enh.enh[,c("enhancerstart1","enhancerstop1", "enhancerstart2", "enhancerstop2")])
      min.min <- min(tmp.enh.enh[,c("enhancerstart1","enhancerstop1", "enhancerstart2", "enhancerstop2")])
      
      df.contacts.cluster.analysis[i,"span.cluster"] <-  max.max - min.min
      
    }
    else{
      
      df.contacts.cluster.analysis[i, "n.genes"] <- length(unique(tmp.gene.enh$geneSymbol))
      df.contacts.cluster.analysis[i, "n.enhancers"] <- length(unique(c(tmp.gene.enh$enhancer,tmp.enh.enh$enhancer1, tmp.enh.enh$enhancer2)))
      
      max.gene.enh <- max(tmp.gene.enh[,c("TSS","TES", "enhancerstart", "enhancerstop")])
      max.enh.enh <- max(tmp.enh.enh[,c("enhancerstart1","enhancerstop1", "enhancerstart2", "enhancerstop2")])
      
      min.gene.enh <- min(tmp.gene.enh[,c("TSS","TES", "enhancerstart", "enhancerstop")])
      min.enh.enh <- min(tmp.enh.enh[,c("enhancerstart1","enhancerstop1", "enhancerstart2", "enhancerstop2")])
      
      max.max <- max(max.gene.enh, max.enh.enh)
      min.min <- min(min.gene.enh, min.enh.enh)
      
      df.contacts.cluster.analysis[i,"span.cluster"] <-  max.max - min.min
      
      
    }
  i <- i +1
  
  }
  return(df.contacts.cluster.analysis)
}

df.contacts.cluster.analysis <- HiC.intensity.contacts(tmp.links, genes.enhancers.EGRM, enhancers.enhancers.EGRM) 
nrow(df.contacts.cluster.analysis) == length(unique(tmp.links$n.cluster))
#TRUE

#With respect to the activity by contact methodology presented in Fulco et al preprint, we integrate DHS peaks data and Chip-Seq data for H3k27ac
#We retrieved data for bipolar neurons and try to integrate its in our cluster analysis. 
#files are Narrowpeak tables in .bed format
columns <- c("chrom", "start", "end", "name", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")

dnase_peaks <- read.table(file.choose(), header = FALSE, sep="\t", stringsAsFactors = FALSE, quote="")
chipseq_peaks <- read.table(file.choose(), header = FALSE, sep="\t", stringsAsFactors = FALSE, quote="")

colnames(dnase_peaks) <- columns ; colnames(chipseq_peaks) <- columns

