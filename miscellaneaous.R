library(GenomicRanges)
library(Rcpp)
library(HiTC)

sourceCpp("straw-master/R/straw-R.cpp")

#Vecteur de chromosomes pour l'automatisation de l'importation des fichiers hic via straw
chroms <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14",
            "chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")

import.hic <- function(file.hic, chrom1, chrom2, metric, resolution, normalization = "NONE"){
  #Fonction d'automatisation 
  print(paste(normalization, file.hic, chrom1, chrom2, metric, resolution, collapse = " "))
  straw_R(paste(normalization, file.hic, chrom1, chrom2, metric, resolution, collapse = " "))
}

for(chrom in chroms){
  assign(paste0("hic.",chrom),import.hic("Neu_inter_sample_372787143.hic", chrom, chrom, "BP", "5000", "VC"))
}

convert.hic.to.HTClist <- function(hic.file){
  #Voir documentation HiTC pour le format de fichier en inut de HiTC
  #output retourné: fichier .mat et fichier .bed
  bed <- hic.file[,c(1,2)]
  tmp <- cbind(unique(hic.file$x), seq(1,length(unique(hic.file$x)), 1 ))
  colnames(tmp) <- c("y", "index")
  
  tmp.merge <- merge(hic.22.t, tmp, by="y")
  
  colnames(tmp) <- c("x", "index")
  
  tmp.merge <- merge(tmp.merge, tmp, by="x")
  
  mat <- tmp.merge[, c("index.y","index.x","counts")]
  mat <- mat[with(mat, order(index.y,index.x)),]
  
  bed <- cbind("chr22",bed)
  bed$index <- seq(1,nrow(bed),1)
  bed$col5 <- 0 
  bed$col6 <- "*"
  
  write.table(mat,paste0(hic.file,".mat"), col.names = F, row.names = F, sep = "\t")
  write.table(bed,paste0(hic.file, ".bed"), col.names = F, row.names = F, sep = "\t")
}

#A automatiser....
tmp <- cbind(unique(hic.22.t$x), seq(1,length(unique(hic.22.t$x)), 1 ))
colnames(tmp) <- c("y", "index")
tmp.merge <- merge(hic.22.t, tmp, by="y")
colnames(tmp) <- c("x", "index")
tmp.merge <- merge(tmp.merge, tmp, by="x")

mat <- tmp.merge[, c("index.y","index.x","counts")]
mat <- mat[with(mat, order(index.y,index.x)),]

bed <- hic.22.t[,c(1,2)]
bed <- cbind("chr22",bed)
bed$index <- seq(1,nrow(bed),1)
bed$col5 <- 0 
bed$col6 <- "*"

write.table(mat, "test.mat", col.names = F, row.names = F, sep = "\t")
write.table(bed, "test.bed", col.names = F, row.names = F, sep = "\t")

#Modification du fichier via la commande shell: sed -i -e 's/"//g' fichier.bed

#Importation des fichiers dans le format HiTC
hic <- importC("test.mat", "test.bed")

#Binning des datas en vue de l'exportation au format my5C
hic.binning <- binningC(hic$chr22chr22, method = "mean", optimize.by = "memory")
#Via la commande suivante on peut récupérer un fichier au format .my5C avec lq version du génome souhaitée
export.my5C(hic.binning, "contacts_matrix", genome = "hg19")

#Visualisation avec Sushi 
library(Sushi)

#Fonctions de contröle:
isComplete(hic)
isPairwise(hic)
#Quels chromosomes ?
seqlevels(hic)

summary(hic$chr22chr22)# affiche les stqtistiques descriptives du fichier

#HeatMqp contact:
#Résolution sur 10kB 
heatmap <- HTClist(mclapply(hic, binningC,binsize=10000, bin.adjust=FALSE, method="sum", step=1, optimize = "memory"))
mapC(heatmap)

#Fragments de restriction 
resFrag <- getRestrictionFragmentsPerChromosome(resSite="AAGCTT", overhangs5=1,chromosomes="chr22",genomePack="BSgenome.Hsapiens.UCSC.hg19")
resFrag

map_hg18 <- NULL

#Sites de coupures
cutSites <- getAnnotatedRestrictionSites(resSite="AAGCTT", overhangs5=1,chromosomes="chr22",genomePack="BSgenome.Hsapiens.UCSC.hg19",wingc=200, mappability=map_hg18, winmap=500)
head(cutSites)

t_chr22annot <- setGenomicFeatures(t$chr22chr22, cutSites)
x_intervals(t_chr22annot);y_intervals(t_chr22annot)

hic_x.binned <- binningC(hic$chr22chr22, binsize=500000, method="median", step=3)
mapC(hic_x.binned,tracks=list(RefSeqGene=genes),maxrange=10,ti="chr22 contacts")


#La génération des graphiques descriptifs
CQC(file)

#Genome coverage by the larger network
tmp.cvg <- df.genes.enh.epg[df.genes.enh.epg$geneSymbol %in% larger.node$id | df.genes.enh.epg$enhancer %in% larger.node$id,]
genes.cvg <- IRanges(start = tmp.cvg$TSS, end = tmp.cvg$TES, names = tmp.cvg$geneSymbol)
sum(width(reduce(genes.cvg)))
#308560 bp covered by genes
enh.cvg <- IRanges(start = tmp.cvg$enhancerstart, end = tmp.cvg$enhancerstop, names = tmp.cvg$enhancer)
sum(width(reduce(enh.cvg)))
#17200 bp covered by genes

findOverlaps(genes.cvg, enh.cvg)
#Here we define community cluster for the larger network
comm.larger <- cluster_infomap(larger.graph)
plot(comm.larger, larger.graph,vertex.size = 1.5, margin=-.1,asp=.25,vertex.color = larger.node$col, edge.color = "red",main="Gene-Enhancer Community Clusters based on their contacts for the larger network")

#GO for larger network
geneS <- unique(tmp.cvg$geneSymbol)

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
searchAttributes(mart = ensembl, pattern = "GO")

results <- getBM(attributes = c('hgnc_symbol',"namespace_1003", "name_1006"),
                 filters = c('hgnc_symbol'), 
                 values = list(geneS),
                 mart = ensembl)

which.max(table(results[results$namespace_1003 == "molecular_function", ]$name_1006))
sort(table(results[results$namespace_1003 == "molecular_function", ]$name_1006), decreasing=T)

#Here we can define a compartment analysis to cluster genes between them based on hierarchical clustering, euclidean distance is used
clusters <- hclust(dist(window.500kb[,c("start.x", "end.y")]))
plot(clusters)
clusterCut <- cutree(clusters,3)
window.500kb$cluster <- clusterCut
window.500kb

plot(window.500kb$start.x, window.500kb$end.y, main="Linear mapping between genes on 500 kbp window size", 
     col=factor(window.500kb$cluster),xlab="TSS", ylab="enhancerstop")

hcluster.coverage.analysis <- function(df.cluster){
  #Function which analyses genomic coverage by cluster
  #Function returns the nuber of genes in cluster, genomic coverage and common genomic coverage, see above for more details
  clusters <- unique(df.cluster$cluster)
  for(clust in clusters){
    tmp.clusters <- df.cluster[df.cluster$cluster == clust,]
    unique.cluster.gene <- unique(tmp.clusters$geneSymbol)
    tmp.enhancer <- GRanges(seqnames=tmp.clusters$seqnames.y, ranges = IRanges(start=tmp.clusters$start.y,end=tmp.clusters$end.y))
    tmp.gene <- GRanges(seqnames=tmp.clusters$seqnames.x, ranges = IRanges(start=tmp.clusters$start.x, end=tmp.clusters$end.x))
    tmp.union <- union(tmp.enhancer, tmp.gene, ignore.strand=TRUE)
    tmp.intersect <- intersect(tmp.enhancer, tmp.gene, ignore.strand=TRUE)
    print(paste("Number of genes in cluster n°",clust,": ", length(unique.cluster.gene)))
    print(paste("Cluster n°: ",clust, "---> Genomic coverage: ", sum(width(reduce(tmp.union)))))
    print(paste("Cluster n°: ",clust, "---> Common Genomic coverage: ", sum(width(reduce(tmp.intersect)))))
  }
}

hcluster.coverage.analysis(window.500kb)
#[1] "Number of genes in cluster n° 1 :  2"
#[1] "Cluster n°:  1 ---> Genomic coverage:  54252"
#[1] "Cluster n°:  1 ---> Common Genomic coverage:  3599"

#[1] "Number of genes in cluster n° 2 :  5"
#[1] "Cluster n°:  2 ---> Genomic coverage:  173020"
#[1] "Cluster n°:  2 ---> Common Genomic coverage:  56655"

#[1] "Number of genes in cluster n° 3 :  5"
#[1] "Cluster n°:  3 ---> Genomic coverage:  135262"
#[1] "Cluster n°:  3 ---> Common Genomic coverage:  22131"

#Linear plot on 500kbp window to determine proximity between genes based on enhancerstop
ex.chr12 <- genes.enhancers.EGRM[genes.enhancers.EGRM$chr=="chr12",]

#We plot only above top diagonal
tmp.X <- ex.chr12[ex.chr12$TSS <= 5000000,]$TSS
tmp.Y <- ex.chr12[ex.chr12$TSS <= 5000000,]$enhancerstop
vec.X <- ifelse(tmp.X<=tmp.Y,tmp.X,tmp.Y)
vec.Y <- ifelse(tmp.X>tmp.Y,tmp.X,tmp.Y)
pdf("linear_plots/linear_plot_example_neurons.pdf")
plot(vec.X, vec.Y,ylim=c(min(vec.X),max(vec.Y)),
     xlab="TSS",ylab="enhancerstop", main="Gene-enhancer proximity on a 5Mb window",
     col=compo.epg$membership[ex.chr12[ex.chr12$TSS<=5000000,]$geneSymbol])
#text(vec.X, vec.Y, ex.chr12[ex.chr12$TSS<=5000000,]$geneSymbol,cex=0.6,col="red")
abline(0,1)
dev.off()

#Extraction of the larger network and summary analysis
larger<- which.max(table(compo.epg$membership))
larger.node <- nodes.epg[compo.epg$membership == larger[[1]],]
larger.links <- links.EGRM[links.EGRM$from%in%larger.node$id | links.EGRM$to%in% larger.node$id,]
larger.graph = graph_from_data_frame(d=larger.links,directed=F,vertices=larger.node)
V(larger.graph)$label.cex = 0.5
V(larger.graph)$label.color = larger.node$col
plot(larger.graph,vertex.size = 1.5, margin=-.1,asp=.25,vertex.color = larger.node$col, edge.color = "red",main="Larger Network")

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

genes.enhancers.EGRM.TADs5 <- findOverlaps(genes.enhancers.EGRM.GRanges, TADS5)
gecTADs5 = countLnodeHits(genes.enhancers.EGRM.TADs5)
table(gecTADs5)
# 0   1   2   3   4   7 
#458 719 116  15   5   2 
round(tally(gecTADs5,format="percent"))
#0  1  2  3  4  7 
#35 55  9  1  0  0 
round(tally(gecTADs5[gecTADs5>0],format="percent"))
#1  2  3  4  7 
#84 14  2  1  0
#17% of enhancer-promoter pairs span 2 or more TADs, which is the same as results seens in other tissues

