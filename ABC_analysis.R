#After launching the ABC score program from Fulco and al 2019, we proced at some descriptive analysis
library(rtracklayer)
library(MASS)
library(biomaRt)
library(BSgenome.Hsapiens.UCSC.hg19)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(igraph)

setwd("/home/nash/Documents/recherche/data/ABC_impl/")

#files are available from master branch
Enhancers_Pred_txt <- read.table("Predictions/EnhancerPredictions.txt", header=TRUE, sep="\t")
Gene_Pred <- read.table("Predictions/GenePredictionStats.txt", header = TRUE)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes.hg19 <- genes(txdb)

symbol <- select(org.Hs.eg.db,keys = genes.hg19$gene_id, columns = c("SYMBOL"), keytype = "ENTREZID")
genes.hg19$geneSymbol <- symbol$SYMBOL 

#Extension to integrate gene promoters (+-500bp dowstream/upstream)
start(genes.hg19) <- start(genes.hg19) - 500
end(genes.hg19) <- end(genes.hg19) + 500


mean(Gene_Pred$nDistalEnhancersPredicted)
#2.11 : Nombre d'enhancers associés moyen par gène 
summary(Gene_Pred$nDistalEnhancersPredicted)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   1.000   2.000   2.117   3.000  19.000 

#Analysis between Expressed and non-Expressed genes
table(Gene_Pred$geneIsExpressed)
#False True
#9761   14728

#Analyse séparée entre gène exprimé et gène non exprimé
NEG <- Gene_Pred[Gene_Pred$geneIsExpressed=="False",]
EG <- Gene_Pred[Gene_Pred$geneIsExpressed=="True",]

MNEG <- NEG[NEG$nDistalEnhancersPredicted == max(NEG$nDistalEnhancersPredicted),]
#Gène non-exprimé le plus associé : LOC101927827 --> 19 enhancers prédits

#Most connected expressed gene
MEG <- EG[EG$nDistalEnhancersPredicted == max(EG$nDistalEnhancersPredicted),]
#Gène exprimé le plus associé : MIR9-3HG --> 14 enhancers prédits 

#Différence de distribution entre les gènes exprimés et non-exprimés
#Because of non normal distribution of Number of predicted enhancers, we do a Wilcoxon test to compare median
wilcox.test(NEG$nDistalEnhancersPredicted, EG$nDistalEnhancersPredicted)
# p-value <  2.2e-16 --> We reject H0

#One analysis for each kind of genes --> expressed vs non-expressed ?

#Correlation between number of candidat elements and number of connected elements ? 
plot(Gene_Pred$nDistalEnhancersPredicted,Gene_Pred$nEnhancersConsidered)
cor(Gene_Pred$nEnhancersConsidered, Gene_Pred$nDistalEnhancersPredicted)

cor(NEG$nEnhancersConsidered, NEG$nDistalEnhancersPredicted)
cor(EG$nEnhancersConsidered, EG$nDistalEnhancersPredicted)

#Modèle entre le nombre d'enhancers prédits associés comme fonction du nombre d'éléments candidats et du fait que le gène soit exprimé (ajout de terme d'intéraction)
#Nature discrète de la variable (dénombrement) --> fitter un modèle de poisson avec lien Log
model_P <- glm(nDistalEnhancersPredicted ~ nEnhancersConsidered + geneIsExpressed + nEnhancersConsidered*geneIsExpressed, data = Gene_Pred,
             family = poisson(link="log"))
summary(model_P)

#Test sur la surdispersion --> Modèle binomial négatif
model_BN <- glm.nb(nDistalEnhancersPredicted ~ nEnhancersConsidered + geneIsExpressed + nEnhancersConsidered*geneIsExpressed, data = Gene_Pred)
summary(model_BN)

l0 <- logLik(model_P)
l1 <- logLik(model_BN)

lik_stat <- 2*(l1-l0)
0.5 * (1-pchisq(lik_stat,1)) < .05 #TRUE --> on garde le modèle binomial négatif avec lien log

#Interprétation du modèle: 
#Significativité statistique --> effet de taille (24489 gènes considérés)
dim(Gene_Pred)

exp(coef(summary(model_BN))[,1]) #--> geneIsExpressed présente une mesure d'effet intéressante
#Le fait que le gène soit exprimé entraine une une diminution de 44% du nombre d'enhancers associés prédits. 
confint(model_BN)
exp(-0.6417731948);exp(-0.493517668) #--> IC 95%: [0.5263583;0.6104752]

#Nombre moyen de gènes associés par enhancers
mean(table(Enhancers_Pred_txt$name)) #--> 1.80

n_unique_enhancers <- length(unique(Enhancers_Pred_txt$name))
#16085
Enhancers_GRanges <- makeGRangesFromDataFrame(Enhancers_Pred_txt, keep.extra.columns = TRUE)

sum(width(reduce(Enhancers_GRanges)))
#15209201 bp

#Import the scz2.prs.txt SNP file 
snps <- read.table("/home/nash/Documents/recherche/data/scz2.prs.txt", header = TRUE, sep="\t")

#Significance filter (0.10)
signi_snps <- snps[snps$p <= .10,]
signi_id <- signi_snps$snpid

nsigni_snps <- snps[snps$p >.10,]
nsigni_id <- nsigni_snps$snpid

nrow(signi_snps)
#35895 significant SNPs with .10 threshold

#This following steps could be replaced by snps_nsigni.csv and snps_signi.csv from master branch
#snp <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

#snps_pos_signi <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
#                  filters = "snp_filter",
#                  values = signi_id,
#                  mart = snp)

#snps_pos_nsigni <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
#                        filters = "snp_filter",
#                        values = nsigni_id,
#                         mart = snp)

snps_pos_signi <- read.table("snps_signi.csv", header=TRUE, sep=";")
snps_pos_nsigni <- read.table("snps_nsigni.csv", header=TRUE, sep=";")

#To avoid issues with biomaRt, snps coordinates were exported in csv format
snps_pos_signi <- snps_pos_signi[snps_pos_signi$chr_name %in% 1:22 | snps_pos_signi$chr_name == "X",]
snps_pos_signi$chr_name <- paste0("chr",snps_pos_signi$chr_name)

snps_pos_nsigni <- snps_pos_nsigni[snps_pos_nsigni$chr_name %in% 1:22 | snps_pos_nsigni$chr_name == "X",]
snps_pos_nsigni$chr_name <- paste0("chr",snps_pos_nsigni$chr_name)

snps_pos_GRanges_signi <- GRanges(seqnames = snps_pos_signi$chr_name, ranges = IRanges(start = snps_pos_signi$chrom_start , end = snps_pos_signi$chrom_start, names = snps_pos_signi$refsnp_id))
snps_pos_GRanges_nsigni <- GRanges(seqnames = snps_pos_nsigni$chr_name, ranges = IRanges(start = snps_pos_nsigni$chrom_start , end = snps_pos_nsigni$chrom_start, names = snps_pos_nsigni$refsnp_id))

snps_pos_GRanges_signi$p <- 0
snps_pos_GRanges_nsigni$p <- 0

signi_snps_f <- signi_snps[signi_snps$snpid %in% names(snps_pos_GRanges_signi),]
nsigni_snps_f <- nsigni_snps[nsigni_snps$snpid %in% names(snps_pos_GRanges_nsigni),]

#To do --> vectorize this part
for(snp in signi_snps_f$snpid){
  
  p <- signi_snps_f[signi_snps_f$snpid == snp,"p"]
  snps_pos_GRanges_signi[snp]$p <- p 
 
}

for(snp in nsigni_snps_f$snpid){
  
  p <- nsigni_snps_f[nsigni_snps_f$snpid == snp,"p"]
  snps_pos_GRanges_nsigni[snp]$p <- p 
  
}

#Here I filter on genes which are present in output ABC Score file 
genes.hg19.filter <- genes.hg19[genes.hg19$geneSymbol %in% Enhancers_GRanges$TargetGene]
genes.hg19[!genes.hg19$geneSymbol %in% Enhancers_GRanges$TargetGene]

#Genes with LOC start names are genes of uncertain function

length(unique(Enhancers_GRanges$TargetGene)) - length(unique(genes.hg19.filter$geneSymbol))
#1390 --> number of lost genes 

#Number of SNPs present in clusters elements
#Present in Enhancers
sum(countOverlaps(Enhancers_GRanges,snps_pos_GRanges_signi)); sum(countOverlaps(snps_pos_GRanges_nsigni, Enhancers_GRanges))
#[1] 435;[1] 731

signi_Snps_in_clusters <- subsetByOverlaps(snps_pos_GRanges_signi,Enhancers_GRanges)
nsigni_Snps_in_clusters <- subsetByOverlaps(snps_pos_GRanges_nsigni,Enhancers_GRanges)

#Present in gene Body
sum(countOverlaps(genes.hg19.filter,snps_pos_GRanges_signi )); sum(countOverlaps(genes.hg19.filter,snps_pos_GRanges_nsigni))
#[1] 9993;[1] 18186

Enhancers_GRanges$nSNP_signi <- countOverlaps(Enhancers_GRanges,snps_pos_GRanges_signi)
Enhancers_GRanges$nSNP_nsigni <- countOverlaps(Enhancers_GRanges,snps_pos_GRanges_nsigni)

genes.hg19.filter$nSNP_signi <- countOverlaps(genes.hg19.filter,snps_pos_GRanges_signi)
genes.hg19.filter$nSNP_nsigni <- countOverlaps(genes.hg19.filter,snps_pos_GRanges_nsigni)

nSNP_signi_clust <- sum(Enhancers_GRanges$nSNP_signi) + sum(genes.hg19.filter$nSNP_signi)
nSNP_nsigni_clust <- sum(Enhancers_GRanges$nSNP_nsigni) + sum(genes.hg19.filter$nSNP_nsigni)

nSNP_signi_hclust <- length(snps_pos_GRanges_signi) - nSNP_signi_clust
nSNP_nsigni_hclust <- length(snps_pos_GRanges_nsigni)- nSNP_nsigni_clust

enrichment.SNPs <- matrix(c(nSNP_signi_clust, nSNP_signi_hclust,nSNP_nsigni_clust,nSNP_nsigni_hclust),nrow = 2,ncol = 2)
#    [,1]  [,2]
#[1,] 10428 18917
#[2,] 25116 47259

#Naive enrichment
#Fisher exact test 
fisher.test(enrichment.SNPs)
#p-value = 0.01179 --> SNP enrichment in clusters

Enhancers_Pred_full_txt <- read.table("/home/nash/Documents/recherche/data/ABC_impl/Predictions/EnhancerPredictionsAllPutative.txt", sep="\t", header = TRUE)
genes_not_clusters <- genes.hg19[!genes.hg19$geneSymbol%in% Enhancers_GRanges$TargetGene]

#Here we suppose that a candidat element is an enhancer
#ENhancers which ain't in clusters are simply the candidat elements with a ABC score under .02
not_significant_candidate_elements <- Enhancers_Pred_full_txt[Enhancers_Pred_full_txt$ABC.Score < 0.02,]

#For efficient we filter only on unique enhancers and make some filters about NAs
unique_candidat <- not_significant_candidate_elements[!duplicated(not_significant_candidate_elements$name),]
unique_candidat <- unique_candidat[!is.na(unique_candidat$start)& !is.na(unique_candidat$end),]

#Grange conversion
enhancers_not_clusters <- GRanges(seqnames=unique_candidat$chr, ranges = IRanges(unique_candidat$start,unique_candidat$end))

#Same steps as previously
sum(countOverlaps(genes_not_clusters, snps_pos_GRanges_signi));sum(countOverlaps(genes_not_clusters, snps_pos_GRanges_nsigni))
sum(countOverlaps(enhancers_not_clusters, snps_pos_GRanges_signi));sum(countOverlaps(enhancers_not_clusters, snps_pos_GRanges_nsigni))

nSNP_signi_hclust_ge <- sum(countOverlaps(genes_not_clusters, snps_pos_GRanges_signi)) + sum(countOverlaps(enhancers_not_clusters, snps_pos_GRanges_signi))
nSNP_nsigni_hclust_ge <- sum(countOverlaps(genes_not_clusters, snps_pos_GRanges_nsigni)) + sum(countOverlaps(enhancers_not_clusters, snps_pos_GRanges_nsigni))

enrichment.SNPs.ge <- matrix(c(nSNP_signi_clust, nSNP_signi_hclust_ge,nSNP_nsigni_clust,nSNP_nsigni_hclust_ge ),nrow = 2,ncol = 2)
#[,1]  [,2]
#[1,] 10428 18917
#[2,]  9350 17407
fisher.test(enrichment.SNPs.ge)
#p-value = 0.1444 --> When we consider genes and enhancers which aren't in clusters, there ain't SNP enrichment in ABC score clusters

library(karyoploteR)

kp <- plotKaryotype()
kpPoints(kp,data=signi_Snps_in_clusters, y=-log(signi_Snps_in_clusters$p), r0=0, r1=0.5,cex = 0.3 ,col="blue")
kpPoints(kp,data=nsigni_Snps_in_clusters, y=-log(nsigni_Snps_in_clusters$p), r0=0, r1=0.5,cex = 0.3 ,col="green")
kpPlotDensity(kp,data=genes.hg19.filter)

df.genes.enhancers <- data.frame(seqnames(Enhancers_GRanges),Enhancers_GRanges$name, Enhancers_GRanges$TargetGene)

nodes <- unique(as.data.frame(rbind(as.matrix(df.genes.enhancers$Enhancers_GRanges.name), as.matrix(df.genes.enhancers$Enhancers_GRanges.TargetGene))))
nodes$typ <- ifelse(nodes$V1 %in% df.genes.enhancers$Enhancers_GRanges.TargetGene, "gene", "enhancer")
nodes$col <- ifelse(nodes$typ == "gene", "blue", "orange")
network <- graph_from_data_frame(d=df.genes.enhancers[,c(2,3)], vertices = nodes, directed = F)

plot(network,vertex.label=NA, vertex.size = 0.5, margin=-.35,asp=.35, vertex.color = nodes$col, edge.color = "green")

get_diameter(network)
mean_distance(network)#4.19
edge_density(network)#0
transitivity(network)#0

compo_ABC <- components(network)
table(compo_ABC$csize)
# 2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23 
#1068  785  479  289  170  151   98   65   60   62   63   46   29   27   29   40   22   24   23   24   10   23 

#24   25   26   27   28   29   30   31   32   33   34   35   36   37   38   39   40   41   42   43   44   45 
#15   12    9   10   10    3   10   10    8    9    7    7    7    5    5    7    8    2    5    3    3    3 

#46   47   48   50   51   52   53   54   55   56   57   59   60   62   63   64   65   66   68   71   73   76 
#4    4    3    1    3    1    1    3    4    2    4    1    3    1    1    2    1    1    2    2    1    1 

#80   81   82   87   89   91   98  105  123  124  125  126  156  157 
#1    1    1    2    2    1    1    1    1    1    1    1    1    1 
summary(compo_ABC$csize)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.000   2.000   4.000   7.312   7.000 157.000 

n.enhancers <- length(unique(df.genes.enhancers$Enhancers_GRanges.name))
n.genes <- length(unique(df.genes.enhancers$Enhancers_GRanges.TargetGene))

monogamous_relationship <- table(compo_ABC$csize)[1][[1]]

monogamous_relationship / n.enhancers # 6.6%
monogamous_relationship / n.genes #9%

#Graphics by chr for cluster analysis
for(chr in unique(df.genes.enhancers$seqnames.Enhancers_GRanges.)){
  
  tmp.genes.enhancers <- df.genes.enhancers[df.genes.enhancers$seqnames.Enhancers_GRanges.==chr,]
  
  nodes <- unique(as.data.frame(rbind(as.matrix(tmp.genes.enhancers$Enhancers_GRanges.name), as.matrix(tmp.genes.enhancers$Enhancers_GRanges.TargetGene))))
  nodes$typ <- ifelse(nodes$V1 %in% tmp.genes.enhancers$Enhancers_GRanges.TargetGene, "gene", "enhancer")
  nodes$col <- ifelse(nodes$typ == "gene", "blue", "orange")
  network <- graph_from_data_frame(d=tmp.genes.enhancers[,c(2,3)], vertices = nodes, directed = F)
  
  plot(network, vertex.label=NA, vertex.size = 2, margin=-.1,asp=.35, vertex.color = nodes$col, edge.color = "green",main=paste("Clustering based on Gene-Enhancer contacts: ", chr))
  legend(x=0, y=-1.3, c("gene","enhancer"), pch=21,
         col="#777777", pt.bg=unique(nodes$col), pt.cex=2, cex=.8, bty="n", ncol=1)
}


for(gene in Enhancers_Pred_txt$TargetGene){
  
  start_gene <- start(genes.hg19.filter[genes.hg19.filter$geneSymbol==gene])
  end_gene <- end(genes.hg19.filter[genes.hg19.filter$geneSymbol==gene])
  
  if(length(start_gene) !=0 & length(end_gene) !=0){
    
    Enhancers_Pred_txt[Enhancers_Pred_txt$TargetGene == gene, "TSS"] <- start_gene
    Enhancers_Pred_txt[Enhancers_Pred_txt$TargetGene == gene, "TES"] <- end_gene
    
  }
  
}

Enhancers_Pred_txt$CellType <- NULL
Enhancers_Pred_txt <- na.omit(Enhancers_Pred_txt)

start.cluster = tapply(pmin(Enhancers_Pred_txt$start,Enhancers_Pred_txt$TSS),
                       compo_ABC$membership[Enhancers_Pred_txt$TargetGene],min)

end.cluster = tapply(pmax(Enhancers_Pred_txt$end,Enhancers_Pred_txt$TES),
                     compo_ABC$membership[Enhancers_Pred_txt$TargetGene],max)

TADS10 <- import(file.choose(), format="bed")

df.TADS10 <- as(TADS10, "data.frame")
df.TADS10$seqnames <- paste0("chr", df.TADS10$seqnames)
df.TADS10 <- df.TADS10[with(df.TADS10, order(seqnames,start)),]
df.TADS10$index.TADs <- seq(1, nrow(df.TADS10), 1)

for (i in unique(df.TADS10$seqnames)) {
  n <- rownames(df.TADS10[df.TADS10$seqnames == i,])
  
  df.TADS10[n,"index.TADs.chr"] <- seq(1, length(n), 1)
}

clusterTADoverlapAnalysis = function(cluster.TADs,TADs)
{
  #Clusters which overlap several TADs
  #Determines the structure of overlapped TADs by cluster
  #and tallies number of clusters with each structure 
  names.several.TADs <- names(table(queryHits(cluster.TADs))[table(queryHits(cluster.TADs))!=1])
  clusters.several.TADs <- cluster.TADs[queryHits(cluster.TADs) %in% names.several.TADs]
  
  distinct.TADs <- 0
  closer.distinct.TADs <- 0 
  nested.TADs <- 0 
  overlaps.TADs <- 0
  
  for(i in unique(queryHits(clusters.several.TADs))){
    tmp.subjectHits <- subjectHits(clusters.several.TADs[queryHits(clusters.several.TADs) == i])
    tmp.TADs <- TADs[tmp.subjectHits]
    
    n <- length(tmp.TADs)
    
    unique.span <- sum(width(reduce(tmp.TADs)))
    tt.span <- sum(width(tmp.TADs))
    w <- width(tmp.TADs[which.max(width(tmp.TADs))])
    
    # If the genome length uniquely spanned by the TADs equals the 
    # sum of the TADs lengths then all TADs are distincts
    if(unique.span == tt.span){
      distinct.TADs =  distinct.TADs + 1
      # If last TAD for this cluster is n-1 positions from first TAD in global index
      # then all TADs are consecutive
      if(mcols(tmp.TADs)$index.TADs[n] == mcols(tmp.TADs)$index.TADs[1] + n -1){
        closer.distinct.TADs = closer.distinct.TADs +1 
      }
    }
    # Else if the genome length uniquely spanned by the TADs equals the
    # width of the largest TAD, then all other TADs are nested in that largest TAD
    else if(unique.span == w){
      nested.TADs = nested.TADs + 1
    }
    # Else we conclude that we have some other form of overlap
    else{
      overlaps.TADs = overlaps.TADs + 1
    }
    
  }
  
  list(cluster.distinct.TADs=distinct.TADs,cluster.closer.distinct.TADs=closer.distinct.TADs,cluster.overlaps.TADs=overlaps.TADs,cluster.nested.TADs=nested.TADs)
}

TADS10 <- GRanges(df.TADS10)
sum(width(reduce(TADS10)))

chr.cluster = tapply(Enhancers_Pred_txt$chr,compo_ABC$membership[Enhancers_Pred_txt$TargetGene],function(vec) substring(vec[1],4))
cluster.GRanges = GRanges(seqnames = paste0("chr",chr.cluster),ranges = IRanges(start.cluster,end.cluster))

cluster.TADS10 <- findOverlaps(cluster.GRanges, TADS10)
table(countLnodeHits(cluster.TADS10))

#0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26 
#231 488 223  78  29  15  14  11   6   5   7   9   5   6   5   4   7   7   3   3   6   7  12   6   3   5   5 

#27  28  29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53 
#10   9   7   3   2  10  10   3  10   4   9   6   3   5   2   3  10   3   8   8   7  10   6   7   8   5   5 

#54  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71  72  73  74  75  76  77  78  79  80 
#8   3   4  10   9  10  14   9   4   8  11   5   5   8   9   8   7   2   5   7   8  11   9   2   7   7   8 

#81  82  83  84  85  86  87  88  89  90  91  92  93  94  95  96  97  98  99 100 101 102 103 104 105 106 107 
#7   7   4  10  10   6   8   7   4   9  10   7  13   4   5  11   8   8   8   5  11  14  19  10  17  16  25 

#108 109 110 111 112 113 114 115 116 117 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 
#8   6   3   7   9   7   6   4   5   5   7   3   4   6   3  12   3   5   3   2   3   4   6   7   8   7   1 

#136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 
#3   6   4   4   6   6   5   7   3  10   1 140   1   2   4   3   1   5   3   4   7   2   4   3   4   4   4 

#163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 184 185 186 187 188 189 190 
#2   6   2   2   3   4   4   5   5   3   5   3   3   6   4   3   4   2   4   3   4   1   6   3   4   5   2 

#191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 214 215 216 217 218 
#4   8   4   3   7   4   1   7   3   7   2   2   3   6   2   4   2   1   4   6   4   4   7   1   6   6   4 

#219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 
#3   4   5   1   2   2   6   2   5   3   4   3   3   1   2   5   2   3   2   4   3   5   2   1   3   7   3 

#246 247 248 249 250 251 252 253 254 256 258 259 260 261 262 263 264 265 266 267 268 269 271 273 280 282 285 
#7   2   1   2   1   2   2   2   3   3   7   3   3   3   3   2   2   1  14   8   6   1   1   2   4   1   1 

#287 288 289 290 293 294 298 311 312 315 317 321 322 324 325 328 330 337 356 395 
#1   2   1   1   1   1   2   1   1   1   1   1   1   1   1   1   1   1   1   1

vec.X = ifelse(Enhancers_Pred_txt$TSS<=Enhancers_Pred_txt$end,Enhancers_Pred_txt$TSS, Enhancers_Pred_txt$start)
vec.Y = ifelse(Enhancers_Pred_txt$TSS>Enhancers_Pred_txt$end,Enhancers_Pred_txt$TES, Enhancers_Pred_txt$end)

genes.enhancers.EGRM.GRanges=GRanges(seqnames = Enhancers_Pred_txt$chr,ranges = IRanges(vec.X,vec.Y))

cluster.TADs10.overlap = clusterTADoverlapAnalysis(cluster.TADS10,TADS10)
cluster.TADs10.overlap
#$cluster.distinct.TADs
#[1] 269

#$cluster.closer.distinct.TADs
#[1] 263

#$cluster.overlaps.TADs
#[1] 1562
