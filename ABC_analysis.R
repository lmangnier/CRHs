#After launching the ABC score program from Fulco and al 2019, we proced at some descriptive analysis
library(rtracklayer)
library(MASS)
library(biomaRt)
library(BSgenome.Hsapiens.UCSC.hg19)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(igraph)

setwd("/home/nash/Documents/recherche/data/ABC_impl/")
Enhancers_Pred <- import("Predictions/EnhancerPredictions.bedpe", format="bedpe")
Enhancers_Pred_txt <- read.table("Predictions/EnhancerPredictions.txt", header=TRUE, sep="\t")
Gene_Pred <- read.table("Predictions/GenePredictionStats.txt", header = TRUE)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes.hg19 <- genes(txdb)

symbol <- select(org.Hs.eg.db,keys = genes.hg19$gene_id, columns = c("SYMBOL"), keytype = "ENTREZID")
genes.hg19$geneSymbol <- symbol$SYMBOL 

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
snps <- read.table("/home/nash/Documents/recherche/data/scz2.prs.txt", header = TRUE, sep="\t")

signi_snps <- snps[snps$p <= .10,]
signi_id <- signi_snps$snpid
nsigni_snps <- snps[snps$p >.10,]
nsigni_id <- nsigni_snps$snpid

nrow(signi_snps)
#35895 snps significatifs au seuil de 10%
snp <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
snp <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")


snps_pos_signi <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
                  filters = "snp_filter",
                  values = signi_id,
                  mart = snp)

snps_pos_nsigni <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
                         filters = "snp_filter",
                         values = nsigni_id,
                         mart = snp)
#To avoid issues with biomaRt, snps coordinates were exported in csv format
snps_pos_signi <- snps_pos_signi[snps_pos_signi$chr_name %in% 1:22 | snps_pos_signi$chr_name == "X",]
snps_pos_signi$chr_name <- paste0("chr",snps_pos_signi$chr_name)

snps_pos_nsigni <- snps_pos_nsigni[snps_pos_nsigni$chr_name %in% 1:22 | snps_pos_nsigni$chr_name == "X",]
snps_pos_nsigni$chr_name <- paste0("chr",snps_pos_nsigni$chr_name)

snps_pos_GRanges_signi <- GRanges(seqnames = snps_pos_signi$chr_name, ranges = IRanges(start = snps_pos_signi$chrom_start , end = snps_pos_signi$chrom_start, names = snps_pos_signi$refsnp_id))
snps_pos_GRanges_nsigni <- GRanges(seqnames = snps_pos_nsigni$chr_name, ranges = IRanges(start = snps_pos_nsigni$chrom_start , end = snps_pos_nsigni$chrom_start, names = snps_pos_nsigni$refsnp_id))

genes.hg19.filter <- genes.hg19[genes.hg19$geneSymbol %in% Enhancers_GRanges$TargetGene]

length(unique(Enhancers_GRanges$TargetGene))-length(unique(genes.hg19.filter$geneSymbol))
#1390

#Number of SNPs present in clusters elements
#Present in Enhancers
sum(countOverlaps(Enhancers_GRanges,snps_pos_GRanges_signi)); sum(countOverlaps(snps_pos_GRanges_nsigni, Enhancers_GRanges))
#[1] 435;[1] 731

#Present in gene Body
sum(countOverlaps(genes.hg19.filter,snps_pos_GRanges_signi )); sum(countOverlaps(genes.hg19.filter,snps_pos_GRanges_nsigni))
#[1] 9858; [1] 17916

Enhancers_GRanges$nSNP_signi <- countOverlaps(Enhancers_GRanges,snps_pos_GRanges_signi)
Enhancers_GRanges$nSNP_nsigni <- countOverlaps(Enhancers_GRanges,snps_pos_GRanges_nsigni)

genes.hg19.filter$nSNP_signi <- countOverlaps(genes.hg19.filter,snps_pos_GRanges_signi)
genes.hg19.filter$nSNP_nsigni <- countOverlaps(genes.hg19.filter,snps_pos_GRanges_nsigni)

nSNP_signi_clust <- sum(Enhancers_GRanges$nSNP_signi) + sum(genes.hg19.filter$nSNP_signi)
nSNP_nsigni_clust <- sum(Enhancers_GRanges$nSNP_nsigni) + sum(genes.hg19.filter$nSNP_nsigni)

nSNP_signi_hclust <- length(snps_pos_GRanges_signi) - nSNP_signi_clust
nSNP_nsigni_hclust <- length(snps_pos_GRanges_nsigni)- nSNP_nsigni_clust

enrichment.SNPs <- matrix(c(nSNP_signi_clust, nSNP_signi_hclust,nSNP_nsigni_clust,nSNP_nsigni_hclust ),nrow = 2,ncol = 2)

#Fisher exact test 
fisher.test(enrichment.SNPs)
#p-value = 0.008701 --> SNP enrichment in clusters

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
