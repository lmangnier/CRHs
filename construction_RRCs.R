library(igraph)
library(GenomicRanges)
library(rtracklayer)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(plyr)
library(GenomicFeatures)
setwd("/home/loic/Documents/HiC/data/export_3Dfeatures/NEU/")

EnhancersPred = read.table("EnhancerPredictions.txt", header=T)
#chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22  chr3  chr4  chr5  chr6  chr7  chr8  chr9  chrX 
#7119  3014  4074  3690  1461  2218  2275  2933  4314  1091  4249  5365  2289   867  1952  3933  2663  3554  3294  3770  2683  2918  3000

#Suppression du chrX des analyses
EnhancersPred.WX = EnhancersPred[EnhancersPred$chr!="chrX",]
EnhancersPred.WX$chr = droplevels(EnhancersPred.WX$chr)

EnhancersPred.WX$typeOf = sub("\\|.*", "", EnhancersPred.WX$name)

#Nombre d'elements geniques et inter-geniques
table(unique(EnhancersPred.WX[,c("name", "typeOf")])[,"typeOf"])
# genic intergenic 
#16980      13819 

genes=unique(EnhancersPred.WX$TargetGene)
enhancers=unique(EnhancersPred.WX$name)

length(genes);length(enhancers)
#[1] 15345
#[1] 30799

#Pour faire la distinction entre les elements presents a l'interieur des introns du gene versus ceux presents dans le gene on recupere
#le promoter du gene
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes.hg19 <- genes(txdb)
symbol <- select(org.Hs.eg.db,keys = genes.hg19$gene_id, columns = c("SYMBOL"), keytype = "ENTREZID")
genes.hg19$geneSymbol <- symbol$SYMBOL

genesInRRCs = genes.hg19[genes.hg19$geneSymbol%in%genes]
promoters = promoters(genesInRRCs)

df.promoters = data.frame(promoters)
df.promoters = df.promoters[,c("seqnames", "start", "end", "geneSymbol")]
colnames(df.promoters) = c("chr", "startProm", "endProm", "TargetGene")

enhancers.promoters = merge(EnhancersPred.WX, df.promoters, by="TargetGene")

length(unique(enhancers.promoters$TargetGene))
#13871
length(unique(enhancers.promoters$name))
#29089

table(unique(enhancers.promoters[,c("name", "typeOf")])[,"typeOf"])
# genic intergenic 
#16100      12989 

summary(as.numeric(table(enhancers.promoters$name)))
summary(as.numeric(table(enhancers.promoters$TargetGene)))
GRanges.Promoters.RRCs = GRanges(seqnames = enhancers.promoters$chr.y, ranges = IRanges(start=enhancers.promoters$startProm, end=enhancers.promoters$endProm))
GRanges.Enhancers.RRCs = GRanges(seqnames = enhancers.promoters$chr.x, ranges = IRanges(start=enhancers.promoters$start, end=enhancers.promoters$end))

GRanges.Enhancers.Prom.RRCs = Pairs(GRanges.Enhancers.RRCs,GRanges.Promoters.RRCs)

unique.Promoters.RRCs = unique(GRanges.Promoters.RRCs)
unique.Enhancers.RRCs = unique(GRanges.Enhancers.RRCs)

length(findOverlaps(unique.Promoters.RRCs,unique.Enhancers.RRCs))/length(GRanges.Enhancers.Prom.RRCs)
#2% des promoters des genes chevauchent les enhancers identifies

#Fichier des RRCs sur la base du score-ABC contenant l'ensemble des contacts genes-enhancers
all.contacts = read.table("EnhancerPredictionsAllPutative.txt", header=T)
all.contacts.WX = all.contacts[all.contacts$chr!="chrX", ]
all.contacts.WX$chr = droplevels(all.contacts.WX$chr)

all.contacts.nRRCs = all.contacts.WX[all.contacts.WX$ABC.Score<0.012,c("chr", "start", "end", "name", "class", "TargetGene")]

rm(all.contacts)

#Recuperation des enhancers non inclus dans les RRCs
#Fulco et al., 2019 ont demontre que le Score-ABC ne performe pas bien pour les contacts promoters-promoters
enhancers.nRRCs = all.contacts.nRRCs[!(all.contacts.nRRCs$name%in%enhancers)&(all.contacts.nRRCs$class!="promoter"),]
enhancers.nRRCs = na.omit(enhancers.nRRCs)

enhancers.promoters.nRRCs = merge(enhancers.nRRCs, df.promoters, by="TargetGene")

GRanges.Promoters.nRRCs = GRanges(seqnames=enhancers.promoters.nRRCs$chr.y, ranges=IRanges(start=enhancers.promoters.nRRCs$startProm, end=enhancers.promoters.nRRCs$endProm))
GRanges.Enhancers.nRRCs = GRanges(seqnames=enhancers.promoters.nRRCs$chr.x, ranges=IRanges(start=enhancers.promoters.nRRCs$start, end=enhancers.promoters.nRRCs$end))


GRanges.Enhancers.Prom.nRRCs = Pairs(GRanges.Enhancers.nRRCs,GRanges.Promoters.nRRCs)
###################################################################################################

sum(width(reduce(first(GRanges.Enhancers.Prom.RRCs))))
# 25 557 347: couverture des enhancers geniques et intergeniques confondus
summary(width(reduce(first(GRanges.Enhancers.Prom.RRCs))))

sum(width(reduce(second(GRanges.Enhancers.Prom.RRCs))))
# 29 551 760: couverture des promoters
summary(width(reduce(second(GRanges.Enhancers.Prom.RRCs))))


#Distance entre le debut du promoteur et le debut du enhancer, sans tenir compte du type de enhancer
dist.prom.enh = enhancers.promoters$startProm-enhancers.promoters$start
dist.prom.enh.gen = enhancers.promoters[enhancers.promoters$typeOf=="genic", "startProm"] - enhancers.promoters[enhancers.promoters$typeOf=="genic", "start"]
dist.prom.enh.inter = enhancers.promoters[enhancers.promoters$typeOf=="intergenic", "startProm"] - enhancers.promoters[enhancers.promoters$typeOf=="intergenic", "start"]

summary(dist.prom.enh)
#   Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-4988075   -44348    -1012     -763    39730 86965203 

summary(abs(dist.prom.enh))
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#3    12165    42022   214735   202885 86965203 

summary(dist.prom.enh.gen)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-4963888   -34524     -670    -1940    33104 85376637 

summary(abs(dist.prom.enh.gen))
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#5    11337    33886   174250   158709 85376637 

summary(dist.prom.enh.inter)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-4988075   -61290    -1120      642    51239 86965203 

summary(abs(dist.prom.enh.inter))
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#3    13498    56138   263101   264058 86965203 

#Nombre de contacts aval et amont: 
n.aval = sum(ifelse(enhancers.promoters$startProm<enhancers.promoters$start,1,0))
n.amont= sum(ifelse(enhancers.promoters$startProm>enhancers.promoters$start,1,0))

n.aval/(n.amont+n.aval)
#Les enhancers presentent-ils un profil d'association similaire aux genes ?
asso.enhancers = as.numeric(table(EnhancersPred.WX$name))
asso.genes = as.numeric(table(EnhancersPred.WX$TargetGene))

summary(asso.enhancers)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   1.000   1.000   2.198   2.000  71.000

summary(asso.genes)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   2.000   4.000   4.386   6.000  26.000 

sd(asso.enhancers); sd(asso.genes)

e = data.frame(cbind("enhancers", asso.enhancers))
g = data.frame(cbind("genes", asso.genes))
colnames(e) = c("type", "n")
colnames(g) = colnames(e)
eg = rbind(e,g) 
eg$n = as.numeric(as.character(eg$n))

#Ici j'utilise une anova et sa contrepartie non-parametrique pour savoir si les resultats obtenus offrent la meme conclusion
anova(lm(n~type, eg)) #difference de moyennes entre les groupes
kruskal.test(n~type,eg) #difference entre les rangs moyens dans les deux groupes

#Parce que les deux tests offrent des resultats significatifs au seuil de 5%, on s'interesse au sens de l'effet, a savoir: 
#Est-ce que les genes sont plus associes que les enhancers ? 
#Sens de l'effet: 
t.test(asso.genes,asso.enhancers, alternative = "greater")
#Les genes sont en moyenne plus associes avec les enhancers

#Creation des paires promoters-enhancers
start.Pair = sapply(GRanges.Enhancers.Prom.RRCs, function(x) min(start(first(x)), start(second(x))))
end.Pair = sapply(GRanges.Enhancers.Prom.RRCs, function(x) max(end(first(x)), end(second(x))))

GRanges.Pair.Prom.Enh = GRanges(seqnames = seqnames(first(GRanges.Enhancers.Prom.RRCs)), ranges=IRanges(start.Pair, end.Pair, names=paste0("Pair",1:length(GRanges.Enhancers.Prom.RRCs))))

##############################################################################################
##################################TADs########################################################
##############################################################################################

##############################Directionality Index############################################
TADs.DI = read.table("DI/finaldomaincalls_allchrs", header = F, sep="\t")
TADs.DI = na.omit(TADs.DI)
TADs.DI.WX = TADs.DI[TADs.DI$V1 != "chr23",]

GRanges.TADs.DI = GRanges(seqnames = TADs.DI.WX$V1, ranges = IRanges(start=TADs.DI.WX$V2, end=TADs.DI.WX$V3, names=paste0("TAD", 1:nrow(TADs.DI.WX))))

sum(width(GRanges.TADs.DI))
#2 547 799 532

length(GRanges.TADs.DI)
#3309

summary(width(GRanges.TADs.DI))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#50001  330001  570001  769961  980001 4810001

export(Pairs(GRanges.TADs.DI,GRanges.TADs.DI), "TADs.DI.bedpe", format="bedpe")

#############################Insulation#########################################################
boundaries.Ins = read.table("INS/output_INS/allchrs_dense_annotated.is400001.ids100001.insulation.boundaries", header=T)
chroms = c(paste0("chr", 1:22), "chrX")

for(chr in chroms){
  boundaries.Ins[grepl(paste0(chr,":"),boundaries.Ins$header, fixed=T),"chr"] = chr
}

boundaries.Ins.WX = boundaries.Ins[boundaries.Ins$chr!="chrX",]

GRanges.boundaries = GRanges(seqnames = boundaries.Ins.WX$chr, ranges=IRanges(start=boundaries.Ins.WX$start,end=boundaries.Ins.WX$end))
GRanges.boundaries.noverlap = reduce(GRanges.boundaries)

split.GRanges.boundaries.noverlap = split(GRanges.boundaries.noverlap, seqnames(GRanges.boundaries.noverlap))

first.by.chr = lapply(split.GRanges.boundaries.noverlap, function(x) {
  GRanges(seqnames = seqnames(x[1]), ranges = IRanges(start = 0, end = start(x[1])))})

GRanges.TADs.INS = lapply(split.GRanges.boundaries.noverlap, function(x) {
  TADs.by.chr = GRanges()
  for(bound in 2:length(x)){
  TADs.by.chr[bound - 1] = GRanges(seqnames = seqnames(x[bound-1]), ranges = IRanges(end(x[bound-1]), start(x[bound])))
  }
  TADs.by.chr
})

grl = GRangesList(unlist(as(GRanges.TADs.INS,"GRangesList")),unlist(as(first.by.chr,"GRangesList")))
TADs.INS = unlist(grl)

length(TADs.INS)
#5211

sum(width(TADs.INS))
#2 493 485 211

summary(width(TADs.INS))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#10001   210001   410001   478504   590001 21130001 

export(Pairs(TADs.INS,TADs.INS),"TADs.INS.bedpe")

#######################################Arrowhead#############################################
#Importation des TADs generes par Arrowhead (sans chrX)
TADs.Rao = import("10000_blocks.bedpe", format="bedpe")

length(first(TADs.Rao))
#4252

sum(width(first(TADs.Rao)))
#1 701 990 000

summary(width(first(TADs.Rao)))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#120000  200000  280000  400280  410000 5470000

sd(width(GRanges.TADs.DI));sd(width(TADs.INS));sd(width(first(TADs.Rao)))

hist(width(TADs.INS))
hist(width(GRanges.TADs.DI))
hist(width(first(TADs.Rao)))

#Test d'adequation des distributions de couverture des TADs
ks.test(width(TADs.INS), width(GRanges.TADs.DI))
ks.test(width(TADs.INS), width(first(TADs.Rao)))
ks.test(width(first(TADs.Rao)), width(GRanges.TADs.DI))

#TADs consensus
OverlappingINSRao = subsetByOverlaps(TADs.INS,TADs.Rao)
OverlappingINSDI = subsetByOverlaps(TADs.INS,GRanges.TADs.DI)
OverlappingRaoDI = subsetByOverlaps(first(TADs.Rao),GRanges.TADs.DI)
#3866 TADs qui se chevauchent entre le score d'Insulation et methode de Rao
#5109 TADs qui se chevauchent entre le score d'Insulation et DI
#4235 TADs qui se chevauchent entre Rao et DI

OverlappingINSRaoDI = subsetByOverlaps(subsetByOverlaps(TADs.INS,TADs.Rao),GRanges.TADs.DI)
#3841 TADs qui se chevauchent entre les 3 definitions (procedure a ameliorer pour prendre le plus grand TAD)


##############################################################################################
################################RRCs & TADs###################################################
##############################################################################################
#Chevauchement des paires genes-enhancers sur les TADs

chisq.PairesTADS= matrix(0, ncol=2, nrow=3)
overlapsPairsDI = table(countOverlaps(GRanges.Pair.Prom.Enh , GRanges.TADs.DI))
#0     1     2     3     4     5     6     7     8     9    10    11    12   103   104 
#2684 50471  7048  1663   495   148    81    39    16     3     3     2     1     1     3 
overlapsPairsDI[-1]/sum(overlapsPairsDI[-1])

#Pourcentage des Paires chevauchant au moins 1 TAD
#           1            2            3            4            5            6            7            8            9           10 
#8.415480e-01 1.175176e-01 2.772868e-02 8.253577e-03 2.467736e-03 1.350585e-03 6.502818e-04 2.667823e-04 5.002168e-05 5.002168e-05 
#11           12          103          104 
#3.334778e-05 1.667389e-05 1.667389e-05 5.002168e-05 

DI.1=overlapsPairsDI[-1][1]
sup1DI = sum(overlapsPairsDI[-1][2:length(overlapsPairsDI[-1])])

chisq.PairesTADS[1,] = c(DI.1, sup1DI)

overlapsPairsINS = table(countOverlaps(GRanges.Pair.Prom.Enh, TADs.INS))
#0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    17   151   152   153   154 
#8525 44493  5725  2152   913   396   196   110    55    40    13     8    15     8     4     1     1     1     1     1  

overlapsPairsINS[-1]/sum(overlapsPairsINS[-1])
# 1            2            3            4            5            6            7            8            9           10 
#8.219201e-01 1.057580e-01 3.975394e-02 1.686587e-02 7.315316e-03 3.620712e-03 2.032032e-03 1.016016e-03 7.389208e-04 2.401493e-04 
#11           12           13           14           17          151          152          153          154 
#1.477842e-04 2.770953e-04 1.477842e-04 7.389208e-05 1.847302e-05 1.847302e-05 1.847302e-05 1.847302e-05 1.847302e-05 


INS.1=overlapsPairsINS[-1][1]
sup1INS= sum(overlapsPairsINS[-1][2:length(overlapsPairsINS[-1])])

chisq.PairesTADS[2,] = c(INS.1, sup1INS)

overlapsPairsRao = table(countOverlaps(GRanges.Pair.Prom.Enh , TADs.Rao))
# 0     1     2     3     4     5     6     7     8     9    10    11    13   100   101 
#17554 33093  8406  2192   788   351   159    64    30     9     6     1     1     3     1 

overlapsPairsRao[-1]/sum(overlapsPairsRao[-1])
#1            2            3            4            5            6            7            8            9           10 
#7.337043e-01 1.863693e-01 4.859879e-02 1.747073e-02 7.782015e-03 3.525186e-03 1.418943e-03 6.651295e-04 1.995388e-04 1.330259e-04 
#11           13          100          101 
#2.217098e-05 2.217098e-05 6.651295e-05 2.217098e-05 

Rao.1=overlapsPairsRao[-1][1]
sup1Rao= sum(overlapsPairsRao[-1][2:length(overlapsPairsRao[-1])])

chisq.PairesTADS[3,] = c(Rao.1, sup1Rao)

#H0: Independance du nombre de chevauchement des paires de genes-enhancers independamment de la methode de creation des TADs
#H1: La methode de creation des TADs impacte le nombre de chevauchement des paires genes-enhancers 
chisq.test(chisq.PairesTADS, correct = F)
#Dependance du nombre de chevauchement des paires genes-enhancers par rapport a la methode de creation des TADs

#Construction des graphs:

nodes = enhancers.promoters[,c("name", "TargetGene")]
vertices = data.frame(rbind(matrix(unique(enhancers.promoters$name), ncol=1), matrix(unique(enhancers.promoters$TargetGene), ncol=1)))

colnames(vertices) = c("name")
vertices$type = vertices[,"name"]%in%enhancers
vertices$col = ifelse(vertices[,"name"]%in%enhancers, "blue", "red")
graph = graph_from_data_frame(d=nodes, directed = F, vertices=vertices)

ego.graph = make_ego_graph(graph)
subgraphs = decompose(graph)

structure.relations.prom = lapply(subgraphs, function(x) {
  if(sum(names(V(x)) %in% unique(enhancers.promoters$TargetGene))==1 & length(V(x))>2){
    return(1)
  }
})

prom.Nenhancers = lapply(subgraphs, function(x) {
  if(sum(names(V(x)) %in% unique(enhancers.promoters$TargetGene))==1 & length(V(x))>2){
    return(x)
  }
})


sum(unlist(structure.relations.prom))
#461 RRCs ou 1 gene est relie a plusieurs enhancers 

prom.Nenhancers = list.clean(prom.Nenhancers,fun = is.null, recursive = F)
plot(prom.Nenhancers[[460]])

compo.graph = components(graph)

compo.graph$no
#1633

#Nombre moyen d'elements par RRCs
summary(as.numeric(table(compo.graph$csize)))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   1.000   2.000   9.226   4.000 252.000 

#metriques de complexite:
length(table(compo.graph$csize))
#177 structures de RRCs differentes 

n.prom = length(unique(enhancers.promoters$name))
n.enhancers = length(unique(enhancers.promoters$TargetGene))

one.prom.one.enh = table(compo.graph$csize)[1][[1]]

one.prom.one.enh/n.prom #0.08 % de promoteurs monogames
one.prom.one.enh/n.enhancers #1% d'enhancers monogames

prom.enhancer.clusters = tapply(names(compo.graph$membership),compo.graph$membership,function(vec,genes) table(vec%in%genes),genes=enhancers.promoters$TargetGene)
prom.enhancer.clusters.mat = matrix(unlist(prom.enhancer.clusters),length(unlist(prom.enhancer.clusters))/2,2,byrow = T)

#Analyse par gene 
tapply(prom.enhancer.clusters.mat[,1],prom.enhancer.clusters.mat[,2],summary)
#Analyse par enhancer
tapply(prom.enhancer.clusters.mat[,2],prom.enhancer.clusters.mat[,1],summary)

#Proportion des relation 1-1-n
sum(prom.enhancer.clusters.mat[prom.enhancer.clusters.mat[,1]==1,2])/n.prom
#0.01344151

sum(prom.enhancer.clusters.mat[prom.enhancer.clusters.mat[,2]==1,1])/n.enhancers
#0.1412299

df.membership = data.frame(compo.graph$membership)
df.membership$name = rownames(df.membership)
df.membership = df.membership[df.membership$name%in%nodes$name,]

enhancers.promoters = merge(enhancers.promoters, df.membership, by="name")
enhancers.promoters$minStart = apply(enhancers.promoters[,c("start", "startProm")], 1, min)
enhancers.promoters$maxEnd = apply(enhancers.promoters[,c("end", "endProm")], 1, max)

start.cluster = aggregate(minStart~compo.graph.membership, enhancers.promoters, min)
end.cluster = aggregate(maxEnd~compo.graph.membership, enhancers.promoters, max)

df.start.end = unique(merge(merge(start.cluster, end.cluster, by="compo.graph.membership"), enhancers.promoters[,c("chr.x", "compo.graph.membership")], by="compo.graph.membership"))


GRanges.cluster = GRanges(seqnames = df.start.end$chr, ranges = IRanges(start=df.start.end$minStart, end=df.start.end$maxEnd, names=df.start.end$compo.graph.membership))


chisq.RRCsTADS= matrix(0, ncol=2, nrow=3)

#Chevauchement de TADs par les RRCs

overlapsRRCsDI = table(countOverlaps(GRanges.cluster , GRanges.TADs.DI))
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   21   22   24   30   42  109 
#56 1044  167   70   45   45   37   30   31   22   20    9   13    9    8    6    3    4    2    4    1    3    1    1    1    1 

overlapsRRCsDI[-1]/sum(overlapsRRCsDI[-1])
#1            2            3            4            5            6            7            8            9           10 
#0.6620164870 0.1058972733 0.0443880786 0.0285351934 0.0285351934 0.0234622701 0.0190234623 0.0196575777 0.0139505390 0.0126823082 
#11           12           13           14           15           16           17           18           19           21 
#0.0057070387 0.0082435003 0.0057070387 0.0050729233 0.0038046925 0.0019023462 0.0025364616 0.0012682308 0.0025364616 0.0006341154 
#22           24           30           42          109 
#0.0019023462 0.0006341154 0.0006341154 0.0006341154 0.0006341154 

DI.1=overlapsRRCsDI[-1][1]
sup1DI = sum(overlapsRRCsDI[-1][2:length(overlapsRRCsDI[-1])])

overlapsRRCsINS = table(countOverlaps(GRanges.cluster , TADs.INS))
#  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  23  24  25  26  27  28  29  30  32  33  34  35  43 
#135 942 130  71  49  40  31  24  28  25  24  17  17  14  11  10   5   8  11   4   5   1   4   3   1   4   2   3   2   1   2   2   2   2   1 
#68 162 
#1   1 

overlapsRRCsINS[-1]/sum(overlapsRRCsINS[-1])
#1            2            3            4            5            6            7            8            9           10 
#0.6288384513 0.0867823765 0.0473965287 0.0327102804 0.0267022697 0.0206942590 0.0160213618 0.0186915888 0.0166889186 0.0160213618 
#11           12           13           14           15           16           17           18           19           20 
#0.0113484646 0.0113484646 0.0093457944 0.0073431242 0.0066755674 0.0033377837 0.0053404539 0.0073431242 0.0026702270 0.0033377837 
#21           23           24           25           26           27           28           29           30           32 
#0.0006675567 0.0026702270 0.0020026702 0.0006675567 0.0026702270 0.0013351135 0.0020026702 0.0013351135 0.0006675567 0.0013351135 
#33           34           35           43           68          162 
#0.0013351135 0.0013351135 0.0013351135 0.0006675567 0.0006675567 0.0006675567 

INS.1=overlapsRRCsINS[-1][1]
sup1INS = sum(overlapsRRCsINS[-1][2:length(overlapsRRCsINS[-1])])

overlapsRRCsRao = table(countOverlaps(GRanges.cluster , TADs.Rao))
#  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  29  31  33  51  58 107 
#348 682 192  79  56  42  28  33  28  20  14  21  10   9   9  10  10   9   4   4   1   3   3   4   1   3   2   2   1   1   1   1   1   1 

overlapsRRCsRao[-1]/sum(overlapsRRCsRao[-1])
#1            2            3            4            5            6            7            8            9           10 
#0.5307392996 0.1494163424 0.0614785992 0.0435797665 0.0326848249 0.0217898833 0.0256809339 0.0217898833 0.0155642023 0.0108949416 
#11           12           13           14           15           16           17           18           19           20 
#0.0163424125 0.0077821012 0.0070038911 0.0070038911 0.0077821012 0.0077821012 0.0070038911 0.0031128405 0.0031128405 0.0007782101 
#21           22           23           24           25           26           27           29           31           33 
#0.0023346304 0.0023346304 0.0031128405 0.0007782101 0.0023346304 0.0015564202 0.0015564202 0.0007782101 0.0007782101 0.0007782101 
#51           58          107 
#0.0007782101 0.0007782101 0.0007782101 

Rao.1=overlapsRRCsRao[-1][1]
sup1Rao = sum(overlapsRRCsRao[-1][2:length(overlapsRRCsRao[-1])])

chisq.RRCsTADS[1,]= c(DI.1, sup1DI)
chisq.RRCsTADS[2,]= c(INS.1, sup1INS)
chisq.RRCsTADS[3,]= c(Rao.1, sup1Rao)

chisq.test(chisq.RRCsTADS, correct = F)

#Proportion des RRCs incluent totalement dans un TAD construit sur la base du score d'insulation
OverlapsRRCINS = findOverlaps(GRanges.cluster, TADs.INS)

uniqueOverlapsRRCINS = names(table(queryHits(OverlapsRRCINS))[table(queryHits(OverlapsRRCINS)) == 1])
subsetOverlapsRRCINS = OverlapsRRCINS[queryHits(OverlapsRRCINS)%in%uniqueOverlapsRRCINS]

RRCsPleinementInclusINS = sum(ifelse(start(GRanges.cluster[queryHits(subsetOverlapsRRCINS)]) >= start(TADs.INS[subjectHits(subsetOverlapsRRCINS)]) & end(GRanges.cluster[queryHits(subsetOverlapsRRCINS)])<=end(TADs.INS[subjectHits(subsetOverlapsRRCINS)]), 1,0))
#Proportion des RRCs chevauchant un seul TAD inclus totalement dans un TAD
RRCsPleinementInclusINS/overlapsRRCsINS[2]
#Proportion des RRCs chevauchant au moins un TAD inclus dans un seul TAD
RRCsPleinementInclusINS/sum(overlapsRRCsINS[-1])

#Proportion des RRCs incluent totalement dans un TAD construit sur la base de DI
OverlapsRRCDI = findOverlaps(GRanges.cluster, GRanges.TADs.DI)

uniqueOverlapsRRCDI = names(table(queryHits(OverlapsRRCDI))[table(queryHits(OverlapsRRCDI)) == 1])
subsetOverlapsRRCDI = OverlapsRRCDI[queryHits(OverlapsRRCDI)%in%uniqueOverlapsRRCDI]

RRCsPleinementInclusDI = sum(ifelse(start(GRanges.cluster[queryHits(subsetOverlapsRRCDI)]) >= start(GRanges.TADs.DI[subjectHits(subsetOverlapsRRCDI)]) & end(GRanges.cluster[queryHits(subsetOverlapsRRCDI)])<=end(GRanges.TADs.DI[subjectHits(subsetOverlapsRRCDI)]), 1,0))
#Proportion des RRCs chevauchant un seul TAD inclus totalement dans un TAD
RRCsPleinementInclusDI/overlapsRRCsDI[2]
#Proportion des RRCs chevauchant au moins un TAD inclus dans un seul TAD
RRCsPleinementInclusDI/sum(overlapsRRCsDI[-1])

#Proportion des RRCs incluent totalement dans un TAD construit sur la base de DI
OverlapsRRCRao = findOverlaps(GRanges.cluster, TADs.Rao)

uniqueOverlapsRRCRao = names(table(queryHits(OverlapsRRCRao))[table(queryHits(OverlapsRRCRao)) == 1])
subsetOverlapsRRCRao = OverlapsRRCRao[queryHits(OverlapsRRCRao)%in%uniqueOverlapsRRCRao]

RRCsPleinementInclusRao = sum(ifelse(start(GRanges.cluster[queryHits(subsetOverlapsRRCRao)]) >= start(TADs.Rao[subjectHits(subsetOverlapsRRCRao)]) & end(GRanges.cluster[queryHits(subsetOverlapsRRCRao)])<=end(TADs.Rao[subjectHits(subsetOverlapsRRCRao)]), 1,0))

#Proportion des RRCs chevauchant un seul TAD inclus totalement dans un TAD
RRCsPleinementInclusRao/overlapsRRCsRao[2]
#Proportion des RRCs chevauchant au moins un TAD inclus dans un seul TAD
RRCsPleinementInclusRao/sum(overlapsRRCsRao[-1])


#Lien entre RRCs chevauchant les TADs et taille des TADs
df.nRRCs.lengthRao = data.frame(matrix(0, ncol=3, nrow=length(TADs.Rao)))
df.nRRCs.lengthINS = data.frame(matrix(0, ncol=3, nrow =length(TADs.INS)))
df.nRRCs.lengthDI = data.frame(matrix(0, ncol=3, nrow =length(GRanges.TADs.DI)))

ORaoRRCs = countOverlaps(TADs.Rao, GRanges.cluster)
ODIRRCs = countOverlaps(GRanges.TADs.DI, GRanges.cluster)
OINSRRCs = countOverlaps(TADs.INS, GRanges.cluster)

for(o in 1:length(ORaoRRCs)){
  
  df.nRRCs.lengthRao[o,] = c("Rao",width(first(TADs.Rao[o])),ORaoRRCs[o])
}

for(o in 1:length(ODIRRCs)){
  
  df.nRRCs.lengthDI[o,] = c("DI",width(GRanges.TADs.DI[o]),ODIRRCs[o])
}

for(o in 1:length(OINSRRCs)){
  
  df.nRRCs.lengthINS[o,] = c("INS",width(TADs.INS[o]),OINSRRCs[o])
}

df.nRRCs.length = data.frame(rbind( df.nRRCs.lengthRao,rbind(df.nRRCs.lengthDI, df.nRRCs.lengthINS)))
df.nRRCs.length$X2 = as.numeric(df.nRRCs.length$X2)
df.nRRCs.length$X3 = as.numeric(df.nRRCs.length$X3)

cor.test(df.nRRCs.length$X2, df.nRRCs.length$X3, method="spearman")

###############################################################################################
##################################RRCs & Compartiments#########################################
###############################################################################################

PC.compartiments.500Kb = read.table("A_B/output/PC/allchrs_PC.500Kb.txt", header=F)
colnames(PC.compartiments.500Kb) = c("chr", "bin", "PC1", "PC2", "PC3")
PC.compartiments.500Kb.WX = PC.compartiments.500Kb[PC.compartiments.500Kb$chr!="chrX",]

PC.compartiments.1Mb = read.table("A_B/output/PC/allchrs_PC.1Mb.txt", header=F)

colnames(PC.compartiments.1Mb) = colnames(PC.compartiments.500Kb)
PC.compartiments.1Mb.WX = PC.compartiments.1Mb[PC.compartiments.1Mb$chr!="chrX", ]

list.PC.1Mb.compartiments = list()
list.PC.500Kb.compartiments = list()

for(chr in chroms[-length(chroms)]){
  PC.compartiments.1Mb.tmp = PC.compartiments.1Mb.WX[PC.compartiments.1Mb.WX$chr ==chr,]
  
  windows.1Mb = seq(0, 1000000*(nrow(PC.compartiments.1Mb.tmp)), 1000000)
  start.stop.1Mb = data.frame(cbind(windows.1Mb[-length(windows.1Mb)], windows.1Mb[-1]))
  colnames(start.stop.1Mb) = c("start", "stop")
  
  
  PC.compartiments.500Kb.tmp = PC.compartiments.500Kb.WX[PC.compartiments.500Kb.WX$chr==chr,]
  
  windows.500Kb = seq(0, 500000*(nrow(PC.compartiments.500Kb.tmp)), 500000)
  start.stop.500Kb = data.frame(cbind(windows.500Kb[-length(windows.500Kb)], windows.500Kb[-1]))
  colnames(start.stop.500Kb) = c("start", "stop")
  
  list.PC.1Mb.compartiments[[chr]] = cbind(PC.compartiments.1Mb.tmp, start.stop.1Mb)
  list.PC.500Kb.compartiments[[chr]] = cbind(PC.compartiments.500Kb.tmp, start.stop.500Kb)
}

PC.compartiments.1Mb.analyse = do.call(rbind, list.PC.1Mb.compartiments)
PC.compartiments.500Kb.analyse = do.call(rbind, list.PC.500Kb.compartiments)

#Le signe des vecteurs propres est arbitraire, multiplication par -1
PC.compartiments.1Mb.analyse[,3:5] = PC.compartiments.1Mb.analyse[,3:5]*-1
PC.compartiments.500Kb.analyse[,3:5] = PC.compartiments.500Kb.analyse[,3:5]*-1

GRanges.PC.1Mb = GRanges(seqnames = PC.compartiments.1Mb.analyse$chr, ranges=IRanges(start = PC.compartiments.1Mb.analyse$start,end =PC.compartiments.1Mb.analyse$stop,names=paste0("bin", 1:nrow(PC.compartiments.1Mb.analyse))),
                     PC1=PC.compartiments.1Mb.analyse$PC1, PC2=PC.compartiments.1Mb.analyse$PC2, PC3=PC.compartiments.1Mb.analyse$PC3)

GRanges.PC.500Kb = GRanges(seqnames = PC.compartiments.500Kb.analyse$chr, ranges=IRanges(start = PC.compartiments.500Kb.analyse$start,end =PC.compartiments.500Kb.analyse$stop,names=paste0("bin", 1:nrow(PC.compartiments.500Kb.analyse))),
                         PC1=PC.compartiments.500Kb.analyse$PC1, PC2=PC.compartiments.500Kb.analyse$PC2, PC3=PC.compartiments.500Kb.analyse$PC3)

Overlaps.RRCs.Compartiments.1Mb = findOverlaps(GRanges.PC.1Mb, GRanges.cluster)
Overlaps.RRCs.Compartiments.500Kb = findOverlaps(GRanges.PC.500Kb, GRanges.cluster)

GRanges.PC.1Mb$ngenes = countOverlaps(GRanges.PC.1Mb, genes.hg19)
GRanges.PC.500Kb$ngenes = countOverlaps(GRanges.PC.500Kb, genes.hg19)

df.PC.1Mb = data.frame(GRanges.PC.1Mb)
df.PC.500Kb = data.frame(GRanges.PC.500Kb)

#Compartementalisation, l'analyse est faite en lien avec la densite de genes, d'enhancers et de RRCs

#Correlation de Spearman entre la densite de genes et la valeur des CP 
cor(df.PC.1Mb$PC1, df.PC.1Mb$ngenes, method="spearman")
cor(df.PC.1Mb$PC2, df.PC.1Mb$ngenes, method="spearman")
cor(df.PC.1Mb$PC3, df.PC.1Mb$ngenes, method="spearman")

#1Mb: Sur l'ensemble des autosomes, la premiere composante traduit la compartementalisation

cor(df.PC.500Kb$PC1, df.PC.500Kb$ngenes, method="spearman")
cor(df.PC.500Kb$PC2, df.PC.500Kb$ngenes, method="spearman")
cor(df.PC.500Kb$PC3, df.PC.500Kb$ngenes, method="spearman")


#Analyse par chromosome
cor.PC1.genes = ddply(df.PC.1Mb, .(seqnames), summarise, "corr" = cor(PC1, ngenes, method = "spearman"))
cor.PC2.genes = ddply(df.PC.1Mb, .(seqnames), summarise, "corr" = cor(PC2, ngenes, method = "spearman"))
cor.PC3.genes = ddply(df.PC.1Mb, .(seqnames), summarise, "corr" = cor(PC3, ngenes, method = "spearman"))

cor.densite.genes.PC = merge(merge(cor.PC1.genes, cor.PC2.genes, by="seqnames"),cor.PC3.genes, by="seqnames")
colnames(cor.densite.genes.PC) = c("chr", "COR.PC1.DENSITE", "COR.PC2.DENSITE", "COR.PC3.DENSITE")

cor.densite.genes.PC

#Distinction entre compartiment A et compartiment B: 
#Utilisation de la mediane ou de la moyenne pour distinguer les compartiments A des compartiments B
densite.moyenne.genes = mean(df.PC.1Mb$ngenes)
densite.mediane.genes = median(df.PC.1Mb$ngenes)

GRanges.PC.1Mb$Compartiment = ifelse(GRanges.PC.1Mb$ngenes>densite.mediane.genes, "A", "B")
df.PC.1Mb$Compartiment = ifelse(df.PC.1Mb$ngenes>densite.mediane.genes, "A", "B")


GRanges.PC.1Mb$nprom = countOverlaps(GRanges.PC.1Mb, unique(second(GRanges.Enhancers.Prom.RRCs)))
GRanges.PC.1Mb$nenh = countOverlaps(GRanges.PC.1Mb, unique(first(GRanges.Enhancers.Prom.RRCs)))

sum(width(GRanges.PC.1Mb))
#2 734 002 734

length(GRanges.PC.1Mb)

length(GRanges.PC.1Mb[GRanges.PC.1Mb$Compartiment=="A"])
#1345
sum(width(GRanges.PC.1Mb[GRanges.PC.1Mb$Compartiment=="A"]))
#1 345 001 345

length(GRanges.PC.1Mb[GRanges.PC.1Mb$Compartiment=="B"])
#1389
sum(width(GRanges.PC.1Mb[GRanges.PC.1Mb$Compartiment=="B"]))
#1 389 001 389

Overlap.RRCs.Comp = findOverlaps(GRanges.cluster, GRanges.PC.1Mb)
nOverlapsRRCsCOMP = countOverlaps(GRanges.cluster, GRanges.PC.1Mb)

table.nOverlapsRRCsCOMP = table(nOverlapsRRCsCOMP)
table.nOverlapsRRCsCOMP[-1] / sum(table.nOverlapsRRCsCOMP[-1])

sum(countOverlaps(GRanges.cluster, GRanges.PC.1Mb[GRanges.PC.1Mb$Compartiment=="A"]))
sum(countOverlaps(GRanges.cluster, GRanges.PC.1Mb[GRanges.PC.1Mb$Compartiment=="B"]))

GRanges.PC.1Mb$nRRCs = countOverlaps(GRanges.PC.1Mb, GRanges.cluster)

sum(GRanges.PC.1Mb$nRRCs[GRanges.PC.1Mb$Compartiment=="A"]) / sum(GRanges.PC.1Mb$nRRCs)

#RRCs inclus dans un seul compartiment: 
RRCsInCOMP=Overlap.RRCs.Comp[queryHits(Overlap.RRCs.Comp) %in%names(nOverlapsRRCsCOMP[nOverlapsRRCsCOMP==1])]
nRRCsInCOMP = 0

for(q in unique(queryHits(Overlap.RRCs.Comp))){
  s=subjectHits(Overlap.RRCs.Comp[queryHits(Overlap.RRCs.Comp)==q])
  if(min(start(GRanges.cluster[q])) >= min(start(GRanges.PC.1Mb[s])) & max(end(GRanges.cluster[q])) <= max(end(GRanges.PC.1Mb[s]))){
    
    nRRCsInCOMP = nRRCsInCOMP + 1
  }
}

nRRCsInCOMP/length(Overlap.RRCs.Comp)

#Nombre d'elements impliques dans les RRCs inclus dans les compartiments A
query.Overlap.RRCs.Comp = queryHits(Overlap.RRCs.Comp)
subject.Overlap.RRCs.Comp = subjectHits(Overlap.RRCs.Comp)

list.RRCs.Comp = list()

for(q in unique(query.Overlap.RRCs.Comp)){
  genes.enhancers.tmp = enhancers.promoters[enhancers.promoters$compo.graph.membership==q,c("chr.x","name","start","end","TargetGene", "startProm", "endProm")]
  
  GRanges.prom = GRanges(seqnames=genes.enhancers.tmp$chr, ranges=IRanges(start=genes.enhancers.tmp$startProm,end=genes.enhancers.tmp$endProm,names=genes.enhancers.tmp$TargetGene))
  GRanges.enh = GRanges(seqnames=genes.enhancers.tmp$chr, ranges=IRanges(start=genes.enhancers.tmp$start,end=genes.enhancers.tmp$end,names=genes.enhancers.tmp$name))
  
  subset.subjects = subjectHits(Overlap.RRCs.Comp[queryHits(Overlap.RRCs.Comp) == q])
  
  GRanges.subset = GRanges.PC.1Mb[subset.subjects]
  
  GRanges.subset$npromRRCsOverlap = countOverlaps(GRanges.subset,unique(GRanges.prom))
  GRanges.subset$nenhancersRRCsOverlap = countOverlaps(GRanges.subset,unique(GRanges.enh))
  
  list.RRCs.Comp[[q]] = GRanges.subset
}
 
proportion.prom.A = lapply(list.RRCs.Comp, function(x){
  if(length(x[x$Compartiment=="A"])>0) {
    sum(x[x$Compartiment=="A"]$npromRRCsOverlap)/ sum(x$npromRRCsOverlap)}
  }
)

proportion.enhancers.A = lapply(list.RRCs.Comp, function(x){
  if(length(x[x$Compartiment=="A"])>0) {
    sum(x[x$Compartiment=="A"]$nenhancersRRCsOverlap)/ sum(x$nenhancersRRCsOverlap)}
}
)

#Lorsque la moyenne est utilisee
# 0.8859238
# 0.8642006

#Utilisation de la mediane
mean(unlist(proportion.prom.A))
#0.9374591
mean(unlist(proportion.enhancers.A))
#0.9163816

export(Pairs(GRanges.PC.1Mb[GRanges.PC.1Mb$Compartiment=="A"],GRanges.PC.1Mb[GRanges.PC.1Mb$Compartiment=="A"]), "compartiments_1Mb.bedpe", format = "bedpe")

###############################################################################################
##################################RRCs & FIREs#################################################
###############################################################################################
FIREs = read.table("FIREs/allchrs_FIREs.matrix", header=T)
FIREs.WX = FIREs[FIREs$chr!="chrX",]

all.FIREs = read.table("FIREs/allchrs_complete_processed.matrix", header=T)
all.FIREs.WX = all.FIREs[all.FIREs$chr!="chrX",]

GRanges.FIREs = GRanges(seqnames = FIREs.WX$chr, ranges = IRanges(FIREs.WX$start, FIREs.WX$stop), ScoreFire = FIREs.WX$ScoreFire)
GRanges.allFIREs = GRanges(seqname = all.FIREs.WX$chr, ranges =IRanges(all.FIREs.WX$start, all.FIREs.WX$stop), pvalues=all.FIREs.WX$pvalues)

Overlaps.Enhancers.FIREs = findOverlaps(unique.Enhancers.RRCs, GRanges.FIREs)
Overlaps.Promoters.FIREs = findOverlaps(unique.Promoters.RRCs, GRanges.FIREs)

mcols(unique.Enhancers.RRCs)["ScoreFire"] = aggregate(GRanges.FIREs, Overlaps.Enhancers.FIREs, ScoreFire=mean(ScoreFire))$ScoreFire
mcols(unique.Promoters.RRCs)["ScoreFire"] = aggregate(GRanges.FIREs, Overlaps.Promoters.FIREs, ScoreFire=mean(ScoreFire))$ScoreFire


summary(FIREs.WX$ScoreFire)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.00287  0.44767  1.53162  1.14550  1.71861 76.60118

summary(mcols(unique.Enhancers.RRCs)$ScoreFire)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.01071 0.49048 1.63214 1.33534 1.80037 5.26158

summary(mcols(unique.Promoters.RRCs)$ScoreFire)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.01587 0.56306 1.65644 1.43711 1.80713 5.26158 


OEnhFIREs = table(countOverlaps(unique.Enhancers.RRCs, GRanges.FIREs))
OEnhFIREs/sum(OEnhFIREs)

OPromFIREs = table(countOverlaps(unique.Promoters.RRCs, GRanges.FIREs))
OPromFIREs/sum(OPromFIREs)
                                       
export(Pairs(GRanges.FIREs, GRanges.FIREs), "FIREs.bedpe", format="bedpe")

#Les enhancers sont-ils enrichis en FIREs?
GRanges.enhancers$FIREs.nsigni = countOverlaps(GRanges.enhancers, GRanges.allFIREs[GRanges.allFIREs$pvalues>0.05])
GRanges.enhancers$FIREs.signi = countOverlaps(GRanges.enhancers, GRanges.allFIREs[GRanges.allFIREs$pvalues<=0.05])

GRanges.enhancers.nRRCs$FIREs.nsigni = countOverlaps(GRanges.enhancers.nRRCs, GRanges.allFIREs[GRanges.allFIREs$pvalues>0.05])
GRanges.enhancers.nRRCs$FIREs.signi = countOverlaps(GRanges.enhancers.nRRCs, GRanges.allFIREs[GRanges.allFIREs$pvalues<=0.05])

table(GRanges.enhancers.nRRCs$FIREs.signi)
enrichissement.FIREs = matrix(c(sum(GRanges.enhancers$FIREs.signi), sum(GRanges.enhancers.nRRCs$FIREs.signi), sum(GRanges.enhancers$FIREs.nsigni), sum(GRanges.enhancers.nRRCs$FIREs.nsigni)), 
                              ncol=2, nrow=2, byrow=T)

#H0: Le fait d'etre FIRE est independant du fait d'etre integre dans les RRCs ou non (RC =1)
#H1: On observe plus de FIREs dans les elements inclus dans les RRCs versus les elements non inclus (RC > 1)
fisher.test(enrichissement.FIREs, alternative = "greater")

#Les enhancers intergeniques sont-ils plus enrichis que les enhancers geniques?
enrichissement.geniqueVSintergenique = matrix(c(sum(GRanges.enhancers[GRanges.enhancers$typeOf=="intergenic"]$FIREs.signi), sum(GRanges.enhancers[GRanges.enhancers$typeOf=="genic"]$FIREs.signi), sum(GRanges.enhancers[GRanges.enhancers$typeOf=="intergenic"]$FIREs.nsigni), sum(GRanges.enhancers[GRanges.enhancers$typeOf=="genic"]$FIREs.nsigni)), 
                              ncol=2, nrow=2, byrow=T)

(enrichissement.geniqueVSintergenique[1,1] * enrichissement.geniqueVSintergenique[2,2]) / (enrichissement.geniqueVSintergenique[1,2] * enrichissement.geniqueVSintergenique[2,1])
#1.14471

fisher.test(enrichissement.geniqueVSintergenique, alternative = "greater")
#L'enrichissement global des enhancers presents dans les RRCs ne semble pas venir du fait que les enhancers soient presents dans le corps du gene

###############################################################################################
#################################RRCs et DI####################################################
###############################################################################################

DI = read.table("DI/all_chrs_dense_annotated.matrix.DI", header=F)
colnames(DI) = c("chr", "start", "end", "DI")
DI$chr = paste0("chr", DI$chr)

DI.WX = DI[DI$chr!="chr23",]
GRanges.DI = GRanges(seqnames = DI.WX$chr, ranges=IRanges(start=DI.WX$start, end=DI.WX$end), DI= DI.WX$DI)

table(countOverlaps(unique.Enhancers.RRCs, GRanges.DI)) / sum(table(countOverlaps(unique.Enhancers.RRCs, GRanges.DI)))
table(countOverlaps(unique.Promoters.RRCs, GRanges.DI)) / sum(table(countOverlaps(unique.Promoters.RRCs, GRanges.DI)))


OEnhDI = findOverlaps(unique.Enhancers.RRCs, GRanges.DI)
OPromDI = findOverlaps(unique.Promoters.RRCs, GRanges.DI)

mcols(unique.Enhancers.RRCs)["DI"] = aggregate(GRanges.DI, OEnhDI, DI = mean(DI))$DI

summary(unique.Enhancers.RRCs$DI)
#Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-1483.3440   -29.3256     0.0001    -1.6459    27.5029  1355.6690 

mcols(unique.Promoters.RRCs)["DI"] = aggregate(GRanges.DI, OPromDI, DI = mean(DI))$DI

summary(unique.Promoters.RRCs$DI)
#  Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-932.4027  -21.5652    0.0109   -0.6515   21.8448  970.8814 

hist(unique.Enhancers.RRCs$DI)
hist(unique.Promoters.RRCs$DI)

t.test(unique.Enhancers.RRCs$DI, unique.Promoters.RRCs$DI)
#Les promoters et les enhancers n'ont pas un profil different en termes de DI

######################################################################################
#########################RRCs et INS##################################################
######################################################################################

INS = read.table("INS/output_INS/allchrs.insulation", header = T)

INS$insulationScore = as.numeric(as.character(INS$insulationScore))
INS$start = as.numeric(as.character(INS$start))
INS$end = as.numeric(as.character(INS$end))

INS = na.omit(INS)
INS$normalizedScore = (INS$insulationScore - min(INS$insulationScore) )/ (max(INS$insulationScore) - min(INS$insulationScore))
INS$chr=gsub(".*[|]([^.]+)[:].*", "\\1", INS$header)

INS.WX = INS[INS$chr!="chrX",]

GRanges.INS = GRanges(seqnames =  INS.WX$chr, ranges = IRanges(start = INS.WX$start, end = INS.WX$end), normalizedScore = INS.WX$normalizedScore)

OEnhINS = findOverlaps(unique.Enhancers.RRCs, GRanges.INS)
OPromINS = findOverlaps(unique.Promoters.RRCs, GRanges.INS)

mcols(unique.Enhancers.RRCs)["normalizedScore"] = aggregate(GRanges.INS,OEnhINS, normalizedScore=mean(normalizedScore))$normalizedScore
mcols(unique.Promoters.RRCs)["normalizedScore"] = aggregate(GRanges.INS,OPromINS, normalizedScore=mean(normalizedScore))$normalizedScore

summary(unique.Enhancers.RRCs$normalizedScore)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.4440  0.7767  0.8093  0.8051  0.8384  0.9437

summary(unique.Promoters.RRCs$normalizedScore)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.5174  0.7668  0.7998  0.7958  0.8292  0.9397

hist(unique.Enhancers.RRCs$normalizedScore)
hist(unique.Promoters.RRCs$normalizedScore)

t.test(unique.Enhancers.RRCs$normalizedScore, unique.Promoters.RRCs$normalizedScore)
#Diffence entre enhancers et promoters au regard de l'INS

#Lien caracteristiques 3D et RRCs:

df.Promoters.RRCs = data.frame(unique.Promoters.RRCs)
df.Promoters.RRCs = df.Promoters.RRCs[c("seqnames", "start", "end", "normalizedScore", "DI", "ScoreFire")]
colnames(df.Promoters.RRCs) = c("chr.x", "startProm", "endProm","normalizedScore.PROM", "DI.PROM", "ScoreFire.PROM")


df.Enhancers.RRCs = data.frame(unique.Enhancers.RRCs)
df.Enhancers.RRCs = df.Enhancers.RRCs[c("seqnames", "start", "end", "normalizedScore", "DI", "ScoreFire")]
colnames(df.Enhancers.RRCs) = c("chr.x", "start", "end","normalizedScore.Enh", "DI.Enh", "ScoreFire.Enh")


enhancers.promoters = merge(enhancers.promoters, df.Promoters.RRCs, by=c("chr.x", "startProm", "endProm"))
enhancers.promoters =merge(enhancers.promoters, df.Enhancers.RRCs, by=c("chr.x", "start", "end"))

head(enhancers.promoters)

meanRRCs.ScoreABC = aggregate(ABC.Score~compo.graph.membership, enhancers.promoters, mean)

#INS
meanRRCsPROM.INS = aggregate(normalizedScore.PROM~compo.graph.membership, enhancers.promoters, mean)
meanRRCsEnh.INS = aggregate(normalizedScore.Enh~compo.graph.membership, enhancers.promoters, mean)

meanRRCs.INS = merge(meanRRCsEnh.INS, meanRRCsPROM.INS)
cor(meanRRCs.INS$normalizedScore.Enh, meanRRCs.INS$normalizedScore.PROM)

meanRRCs.INS$meanINS = apply(meanRRCs.INS[,c("normalizedScore.Enh","normalizedScore.PROM")], 1, mean)

#DI
meanRRCsPROM.DI = aggregate(DI.PROM~compo.graph.membership, enhancers.promoters, mean)
meanRRCsEnh.DI = aggregate(DI.Enh~compo.graph.membership, enhancers.promoters, mean)

meanRRCs.DI = merge(meanRRCsEnh.DI, meanRRCsPROM.DI)
cor(meanRRCs.DI$DI.Enh, meanRRCs.DI$DI.PROM)

meanRRCs.DI$meanDI = apply(meanRRCs.DI[,c("DI.Enh","DI.PROM")], 1, mean)

#Score-FIRE
meanRRCsPROM.FIRE = aggregate(ScoreFire.PROM~compo.graph.membership, enhancers.promoters, mean)
meanRRCsEnh.FIRE = aggregate(ScoreFire.Enh~compo.graph.membership, enhancers.promoters, mean)

meanRRCs.FIRE = merge(meanRRCsEnh.FIRE, meanRRCsPROM.FIRE)
cor(meanRRCs.FIRE$ScoreFire.Enh, meanRRCs.FIRE$ScoreFire.PROM)

meanRRCs.FIRE$meanFIRE = apply(meanRRCs.FIRE[,c("ScoreFire.Enh","ScoreFire.PROM")], 1, mean)


ThreeD.RRCs = merge(meanRRCs.FIRE[,c("compo.graph.membership","meanFIRE")],merge(meanRRCs.DI[,c("compo.graph.membership","meanDI")], meanRRCs.INS[,c("compo.graph.membership","meanINS")], by="compo.graph.membership"), by="compo.graph.membership")

Full.3D.RRCs = merge(ThreeD.RRCs, meanRRCs.ScoreABC, by="compo.graph.membership" )

nByRRCs = lapply(unique(enhancers.promoters$compo.graph.membership), function(x) length(unique(enhancers.promoters[enhancers.promoters$compo.graph.membership==x,"TargetGene"])) + length(unique(enhancers.promoters[enhancers.promoters$compo.graph.membership==x,"name"])))

df.nByRRCs = data.frame(matrix(c(1:length(nByRRCs), unlist(nByRRCs)), ncol=2))
colnames(df.nByRRCs) = c("compo.graph.membership", "nElements")

Full.3D.RRCs = merge(Full.3D.RRCs,df.nByRRCs, by="compo.graph.membership")
rownames(Full.3D.RRCs) = Full.3D.RRCs$compo.graph.membership
Full.3D.RRCs$compo.graph.membership = NULL

heatmap(cor(Full.3D.RRCs))
t.test(Full.3D.RRCs$meanFIRE, GRanges.FIREs$ScoreFire, alternative = "greater")
t.test(Full.3D.RRCs$meanDI, GRanges.DI$DI)
t.test(Full.3D.RRCs$meanINS, GRanges.INS$normalizedScore, alternative = "greater")

#Analyse Graphique des caracteristiques 3D et du score-ABC au niveau individuel
#La moyenne sur le genome est affichee a titre indicatif
par(mfrow = c(1,1))
plot(Full.3D.RRCs$meanFIRE)
abline(h=mean(GRanges.FIREs$ScoreFire),lwd=2,col="blue")

plot(Full.3D.RRCs$meanDI)
abline(h=mean(GRanges.DI$DI),lwd=2, col="blue")

plot(Full.3D.RRCs$meanINS)
abline(h=mean(GRanges.INS$normalizedScore), lwd=2, col="blue")

plot(Full.3D.RRCs$ABC.Score)

#Analyse Graphique des caracteristiques 3D en fonction de la structure des RRCs en termes de nombres d'elements
boxplot(meanFIRE~nElements,Full.3D.RRCs, las=2)
abline(h=mean(GRanges.FIREs$ScoreFire))
boxplot(meanDI~nElements,Full.3D.RRCs, las=2)
abline(h=mean(GRanges.DI$DI))
boxplot(meanINS~nElements,Full.3D.RRCs, las=2)
abline(h=mean(GRanges.INS$normalizedScore))
boxplot(ABC.Score~nElements,Full.3D.RRCs, las=2)
