#to install  FantomEnhancers.hg19 package, thanks to uncomment the following code line:
#devtools::install_github("charlesjb/fantomenhancers.hg19")
library(FantomEnhancers.hg19)

#Have enhancers for each genes, data are extracted from Fantom5
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
#Data are available on Github repo
load(file.choose())
#Let's take a look on data 
head(neurons_enhancers)

#Here we use the lowest filter size
filter.enhancers.cons <- filter.enhancers(neurons_enhancers, "consensus")
filter.enhancers.one <- filter.enhancers(neurons_enhancers, "one")

enhancers.annotations.Y <- enhancers.annotations(filter.enhancers.one,contacts.locus$Y)
enhancers.annotations.X <- enhancers.annotations(filter.enhancers.one,contacts.locus$X1)

names(enhancers.annotations.Y) <- 1:length(enhancers.annotations.Y)
names(enhancers.annotations.X) <- 1:length(enhancers.annotations.X)

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