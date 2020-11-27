library(igraph)
library(GenomicRanges)
library(rtracklayer)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(plyr)
library(GenomicFeatures)
library(coin)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(ggpubr)
library(lme4)
library(PerformanceAnalytics)
setwd("/home/loic/Documents/HiC/data/export_3Dfeatures/NEU/")

enhancers.promoters = process.ABC("EnhancerPredictions.txt")
GRanges.Enhancers.Prom.RRCs = Pairs.Enh.Prom.ABC(enhancers.promoters)

summary(as.numeric(table(enhancers.promoters$name)))
summary(as.numeric(table(enhancers.promoters$TargetGene)))

unique.Promoters.RRCs = unique(second(GRanges.Enhancers.Prom.RRCs))
unique.Enhancers.RRCs = unique(first(GRanges.Enhancers.Prom.RRCs))

length(findOverlaps(unique.Promoters.RRCs,unique.Enhancers.RRCs))/length(GRanges.Enhancers.Prom.RRCs)
#2% des promoters des genes chevauchent les enhancers identifies
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes.hg19 <- genes(txdb)
symbol <- select(org.Hs.eg.db,keys = genes.hg19$gene_id, columns = c("SYMBOL"), keytype = "ENTREZID")
genes.hg19$geneSymbol <- symbol$SYMBOL

DNAseActivity.NEU = import("/home/loic/Documents/HiC/data/export_3Dfeatures/NEU/NEU_DNAse.macs2_peaks.narrowPeak.sorted", format = "narrowpeak")
DNAseActivity.NEU$signalValue = log2(DNAseActivity.NEU$signalValue+1)

H3K27acActivity.NEU = import("/home/loic/Documents/HiC/data/export_3Dfeatures/NEU/NEU_H3k27ac_sorted_macs2_peaks.narrowPeak.sorted", format = "narrowpeak")
H3K27acActivity.NEU$signalValue = log2(H3K27acActivity.NEU$signalValue+1)

H3K4me1Activity.NEU = import("/home/loic/Documents/HiC/data/H3K4me1.peaks.macs2_peaks.narrowPeak", format = "narrowpeak")
H3K4me1Activity.NEU$signalValue = log2(H3K4me1Activity.NEU$signalValue +1)

H3K4me3.NEU.REP1 = import("/home/loic/Documents/HiC/data/NEU_H3k4me3_REP1.bed.gz", format = "narrowpeak")
H3K4me3.NEU.REP2 = import("/home/loic/Documents/HiC/data/NEU_H3k4me3_REP2.bed.gz", format = "narrowpeak")
H3K4me3.NEU.REP3 = import("/home/loic/Documents/HiC/data/NEU_H3k4me3_REP3.bed.gz", format = "narrowpeak")

grl.H3K4me3.NEU = GRangesList(H3K4me3.NEU.REP1, H3K4me3.NEU.REP2, H3K4me3.NEU.REP3)

H3K4me3Activity.NEU = unique(do.call("c", as(grl.H3K4me3.NEU, "GRangesList")))
H3K4me3Activity.NEU$signalValue = log2(H3K4me3Activity.NEU$signalValue+1)

H3K27me3.NEU.REP1 = import("/home/loic/Documents/HiC/data/export_3Dfeatures/NEU/H3K27me3_Rep1.bed.gz", format = "narrowpeak")
H3K27me3.NEU.REP2 = import("/home/loic/Documents/HiC/data/export_3Dfeatures/NEU/H3K27me3_Rep2.bed.gz", format = "narrowpeak")
H3K27me3.NEU.REP3 = import("/home/loic/Documents/HiC/data/export_3Dfeatures/NEU/H3K27me3_Rep3.bed.gz", format = "narrowpeak")

grl.H3K27me3.NEU = GRangesList(H3K27me3.NEU.REP1, H3K27me3.NEU.REP2,H3K27me3.NEU.REP3)

H3K27me3Activity.NEU = unique(do.call("c", as(grl.H3K27me3.NEU, "GRangesList")))
H3K27me3Activity.NEU$signalValue = log2(H3K27me3Activity.NEU$signalValue+1)

CTCF.NEU.REP1 = import("NEU_REP1_CTCF.bed.gz", format = "narrowpeak")
CTCF.NEU.REP2 =   import("NEU_REP2_CTCF.bed.gz", format = "narrowpeak")
grL.CTCF.NEU = GRangesList(CTCF.NEU.REP1, CTCF.NEU.REP2)

CTCFActivity.NEU = unique(do.call("c", as(grL.CTCF.NEU, "GRangesList")))
CTCFActivity.NEU$signalValue = log2(CTCFActivity.NEU$signalValue+1)

p300.NEU.REP1 = import("ENCFF114UVX.bed.gz", format = "narrowpeak")
p300.NEU.REP2 = import("ENCFF287VNA.bed.gz", format = "narrowpeak")
p300.NEU.REP3 = import("ENCFF777EIJ.bed.gz", format = "narrowpeak")
grL.p300.NEU = GRangesList(p300.NEU.REP1, p300.NEU.REP2, p300.NEU.REP3)

p300Activity.NEU = unique(do.call("c", as(grL.p300.NEU, "GRangesList")))
p300Activity.NEU$signalValue = log2(p300Activity.NEU$signalValue+1)

#Gene Expression
expression.NEU = read.table("GSE142670_countdata_20M_neurons.txt", header=T)
expression.NEU$GENEID = rownames(expression.NEU)
expression.NEU[,-ncol(expression.NEU)] = log(expression.NEU[,-ncol(expression.NEU)]+1)
correspondance = ensembldb::select(EnsDb.Hsapiens.v86, keys=rownames(expression.NEU), keytype = "GENEID", columns=c("SYMBOL", "GENEID"))

expression.NEU = merge(expression.NEU, correspondance, by="GENEID")
expression.NEU$geneSymbol =  expression.NEU$SYMBOL
expression.NEU$SYMBOL = NULL

expression.NEU$median.Expr = apply(expression.NEU[,-c(1,ncol(expression.NEU))],1, median)
expression.genes = makeGRangesFromDataFrame(merge(data.frame(genes.hg19), expression.NEU[,c("geneSymbol", "median.Expr")], by="geneSymbol"), keep.extra.columns = T)

###################################################################################################
#Distance entre le debut du promoteur et le debut du enhancer, sans tenir compte du type de enhancer
summary(distance.between.Pairs(GRanges.Enhancers.Prom.RRCs, absolute = F))
#   Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-4988075   -44348    -1012     -763    39730 86965203 

summary(distance.between.Pairs(GRanges.Enhancers.Prom.RRCs, absolute = T))
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#3    12165    42022   214735   202885 86965203 


#Nombre de contacts aval et amont: 
n.aval = sum(ifelse(enhancers.promoters$startProm<enhancers.promoters$start,1,0))
n.amont= sum(ifelse(enhancers.promoters$startProm>enhancers.promoters$start,1,0))

n.aval/(n.amont+n.aval)

#Creation des paires promoters-enhancers
GRanges.Pair.ABC = coverage.By.Pair(GRanges.Enhancers.Prom.RRCs)

#Construction des graphs:
graph.ABC = AL1C.Crn.ABC(enhancers.promoters)
compo.graph.ABC = components(graph.ABC)

decompose(graph.ABC)
compo.graph.ABC$no
#1633


#Nombre moyen d'elements par RRCs
summary(as.numeric(compo.graph.ABC$csize))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.00    3.00    6.00   26.31   24.00  506.00 

#metriques de complexite:
length(table(compo.graph.ABC$csize))
#177 structures de RRCs differentes 
table(compo.graph.ABC$csize)/sum(table(compo.graph.ABC$csize))
complexity.ABC = Complexity.Crn(graph.ABC,enhancers.promoters$name,enhancers.promoters$TargetGene ,extract.central.genes=T)

GRanges.cluster.ABC = coverage.By.Crn(graph.ABC, enhancers.promoters)$clusters
enhancers.promoters = coverage.By.Crn(graph.ABC, enhancers.promoters)$data.frame

summary(width(GRanges.Pair.ABC))
summary(width(GRanges.cluster.ABC))

t.test(width(GRanges.cluster.ABC), width(GRanges.Pair.ABC), alternative = "greater")
wilcox.test(width(GRanges.cluster.ABC), width(GRanges.Pair.ABC))

###############################################################################################
##################################RRCs & Compartiments#########################################
###############################################################################################
t = define.active.compartments("A_B/output/PC/allchrs_PC.1Mb.txt", genes=genes.hg19)

t$ngenes = countOverlaps(t, genes.hg19)
# t$PC1 = t$PC1*-1

Overlaps.Comp.DNAse = findOverlaps(t, DNAseActivity.NEU)
subset.Comp.DNAse = DNAseActivity.NEU[subjectHits(Overlaps.Comp.DNAse)]

mcols(t)[unique(queryHits(Overlaps.Comp.DNAse)), "DNAse"] = aggregate(subset.Comp.DNAse$signalValue, list(queryHits(Overlaps.Comp.DNAse)), mean, na.rm=T)$x

Overlaps.Comp.expr = findOverlaps(t,expression.genes)
subset.Comp.expr = expression.genes[subjectHits(Overlaps.Comp.expr)]

mcols(t)[unique(queryHits(Overlaps.Comp.expr)), "expr"] = aggregate(subset.Comp.expr$median.Expr, list(queryHits(Overlaps.Comp.expr)), mean, na.rm=T)$x

t$Compartment = ifelse(t$PC1>0, "A","B")

t.A = data.frame(mcols(t[t$Compartment=="A"]))
t.A[is.na(t.A)] = 0

cor(t.A$PC1, t.A$ngenes, method="spearman")
cor(t.A$PC1, t.A$expr, method="spearman")
cor(t.A$PC1, t.A$DNAse, method="spearman")


consecutive.Compartment = rle(t$Compartment)

CompartmentsByChr = split(t, seqnames(t))
CompartmentsByChr$chrX =NULL

q = lapply(CompartmentsByChr, function(x){
  l = vector(mode="numeric")
  consecutive.Compartment = rle(x$Compartment)
  for(i in seq_along(consecutive.Compartment$lengths)){
    l = c(l,rep(i, consecutive.Compartment$lengths[i]))
  }
  x$group = l
  d = data.frame(x)
  d
})

r = lapply(q, function(x) {
  
  start.C = aggregate(start~group,x,min)
  end.C = aggregate(end~group,x,max)
  ngenes.C = aggregate(ngenes~group,x,sum)
  expr.C = aggregate(expr~group,x,mean)
  DNAse.C = aggregate(DNAse~group,x,mean)
  PC1.C = aggregate(PC1~group,x,sum)
  C = aggregate(Compartment~group,x,unique)
  cbind(unique(x$seqnames),merge(C ,merge(expr.C,merge(DNAse.C,merge(start.C,merge(end.C,merge(ngenes.C,PC1.C, by="group"), by="group"), by="group"), by="group"), by="group"),by="group"))
  #cbind(unique(x$seqnames),merge(H3K4me1.C,merge(H3K27ac.C,merge(C ,merge(expr.C,merge(DNAse.C,merge(start.C,merge(end.C,merge(ngenes.C,PC1.C, by="group"), by="group"), by="group"), by="group"), by="group"),by="group"), by="group"),by="group"))
})

a = do.call(rbind, r)
GRanges.Compartment.agg = GRanges(seqnames=a$`unique(x$seqnames)`, ranges=IRanges(start=a$start,end=a$end), ngenes=a$ngenes, expr=a$expr,DNAse=a$DNAse, PC1=a$PC1, Compartment=a$Compartment)

A = GRanges.Compartment.agg[GRanges.Compartment.agg$Compartment=="A"]

cor(A$ngenes, A$PC1,method = "spearman")
cor(A$expr, A$PC1, method="spearman")
cor(A$DNAse, A$PC1, method="spearman")

#Lien entre paires et compartiments
#On cherche a connaitre si 2 elements inclus dans une paire sont dans des compartiments A (AA), compartiments B (BB) ou un element dans un
#compartiment A et l'autre dans un compartiment B (AB)

table((countOverlaps(GRanges.Pair.ABC, GRanges.Compartment.agg)))/sum(table((countOverlaps(GRanges.Pair.ABC, GRanges.Compartment.agg))))
table((countOverlaps(GRanges.Pair.Rao, GRanges.Compartment.agg)))/sum(table((countOverlaps(GRanges.Pair.Rao, GRanges.Compartment.agg))))
table((countOverlaps(GRanges.Pair.DNAse, GRanges.Compartment.agg)))/sum(table((countOverlaps(GRanges.Pair.DNAse, GRanges.Compartment.agg))))

k = Pairs.with.Compartments(GRanges.Pair.ABC,GRanges.Compartment.agg)
l = Pairs.with.Compartments(GRanges.Pair.Rao,GRanges.Compartment.agg)
m = Pairs.with.Compartments(GRanges.Pair.DNAse,GRanges.Compartment.agg)
w = prop.table(t(rbind(k,l,m)), 2)
colnames(w) = c("ABC", "Rao", "DNAse")

barplot(w*100,col=c(rgb(0,0,1,0.3), rgb(0.80,0,0.2),rgb(0,1,0,0.3)), xlab="Annotation Type", ylab = "%", xlim=c(0,4.5), main="Repartition of between-compartment contacts for pairs of elements")
legend("topright", c("AA", "AB","BB"), col=c(rgb(0,0,1,0.3), rgb(0.80,0,0.2),rgb(0,1,0,0.3)), lwd=10)

table((countOverlaps(GRanges.cluster.ABC, GRanges.Compartment.agg)))/sum(table((countOverlaps(GRanges.cluster.ABC, GRanges.Compartment.agg))))
table((countOverlaps(GRanges.cluster.Rao, GRanges.Compartment.agg)))/sum(table((countOverlaps(GRanges.cluster.Rao, GRanges.Compartment.agg))))
table((countOverlaps(GRanges.cluster.DNAse, GRanges.Compartment.agg)))/sum(table((countOverlaps(GRanges.cluster.DNAse, GRanges.Compartment.agg))))

kC = Pairs.with.Compartments(GRanges.cluster.ABC,GRanges.Compartment.agg)
lC = Pairs.with.Compartments(GRanges.cluster.Rao,GRanges.Compartment.agg)
mC = Pairs.with.Compartments(GRanges.cluster.DNAse,GRanges.Compartment.agg)

wC = prop.table(t(rbind(kC,lC,mC)), 2)
colnames(wC) = c("ABC", "Rao", "DNAse")

barplot(wC*100,col=c(rgb(0,0,1,0.3), rgb(0.80,0,0.2),rgb(0,1,0,0.3)), xlab="Annotation Type", ylab = "%", xlim=c(0,4.5), main="Repartition of between-compartment contacts for RRCs")
legend("topright", c("AA", "AB","BB"), col=c(rgb(0,0,1,0.3), rgb(0.80,0,0.2),rgb(0,1,0,0.3)), lwd=10)

GRanges.cluster.ABC$nelements = as.numeric(table(enhancers.promoters$membership))
GRanges.cluster.DNAse$nelements = as.numeric(table(df.DNAse$compo.DNAse.membership))
GRanges.cluster.Rao$nelements = as.numeric(table(df.Rao$compo.Rao.membership))

f.ABC = findOverlaps(GRanges.Compartment.agg,GRanges.cluster.ABC)
f.Rao = findOverlaps(GRanges.Compartment.agg,GRanges.cluster.Rao)
f.DNAse = findOverlaps(GRanges.Compartment.agg, GRanges.cluster.DNAse)
  
keep_singles <- function(v){
  v[!(v %in% v[duplicated(v)])] 
}

nd.subject.ABC = keep_singles(subjectHits(f.ABC))
nd.query.ABC = queryHits(f.ABC[subjectHits(f.ABC)%in%nd.subject.ABC])

query.Comp.RRCs.ABC = GRanges.cluster.ABC[nd.subject.ABC]

nd.subject.Rao = keep_singles(subjectHits(f.Rao))
nd.query.Rao = queryHits(f.Rao[subjectHits(f.Rao)%in%nd.subject.Rao])

query.Comp.RRCs.Rao = GRanges.cluster.Rao[nd.subject.Rao]

nd.subject.DNAse = keep_singles(subjectHits(f.DNAse))
nd.query.DNAse = queryHits(f.DNAse[subjectHits(f.DNAse)%in%nd.subject.DNAse])

query.Comp.RRCs.DNAse = GRanges.cluster.DNAse[nd.subject.DNAse]


mcols(GRanges.Compartment.agg)[unique(nd.query.ABC), "mean_nelements_ABC"] = aggregate(query.Comp.RRCs.ABC$nelements, list(nd.query.ABC), mean, na.rm=T)$x
mcols(GRanges.Compartment.agg)[unique(nd.query.Rao), "mean_nelements.Rao"] = aggregate(query.Comp.RRCs.Rao$nelements, list(nd.query.Rao), mean, na.rm=T)$x
mcols(GRanges.Compartment.agg)[unique(nd.query.DNAse), "mean_nelements_DNAse"] = aggregate(query.Comp.RRCs.DNAse$nelements, list(nd.query.DNAse), mean, na.rm=T)$x

mean(mcols(GRanges.Compartment.agg)[GRanges.Compartment.agg$Compartment=="A","mean_nelements_ABC"], na.rm=T)
mean(mcols(GRanges.Compartment.agg)[GRanges.Compartment.agg$Compartment=="B","mean_nelements_ABC"], na.rm=T)

mean(mcols(GRanges.Compartment.agg)[GRanges.Compartment.agg$Compartment=="A","mean_nelements.Rao"], na.rm=T)
mean(mcols(GRanges.Compartment.agg)[GRanges.Compartment.agg$Compartment=="B","mean_nelements.Rao"], na.rm=T)

mean(mcols(GRanges.Compartment.agg)[GRanges.Compartment.agg$Compartment=="A","mean_nelements_DNAse"], na.rm=T)
mean(mcols(GRanges.Compartment.agg)[GRanges.Compartment.agg$Compartment=="B","mean_nelements_DNAse"], na.rm=T)

RRCs.Compartments.Complexity = data.frame(mcols(GRanges.Compartment.agg))
my_comparisons=list(c("A","B")) 

graph.A = ggviolin(RRCs.Compartments.Complexity, x="Compartment", y="mean_nelements_ABC", fill="Compartment",ylab="mean Elements",palette = c("#00AFBB", "#E7B800"),add = "boxplot", add.params = list(fill = "white"))+
           stat_compare_means(comparisons = my_comparisons, label = "p.signif") + stat_compare_means(label.y = 100, label.x=1.35) 
graph.B = ggviolin(RRCs.Compartments.Complexity, x="Compartment", y="mean_nelements.Rao", fill="Compartment",ylab="mean Elements",palette = c("#00AFBB", "#E7B800"),add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + stat_compare_means(label.y = 20, label.x=1.35)
graph.C = ggviolin(RRCs.Compartments.Complexity, x="Compartment", y="mean_nelements_DNAse", fill="Compartment",ylab="mean Elements",palette = c("#00AFBB", "#E7B800"),add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + stat_compare_means(label.y = 25, label.x=1.35) 

ggarrange(graph.A, graph.B, graph.C, 
          labels = c("ABC", "Rao", "Rao+DNAse"),
          ncol = 2, nrow = 2)

###############################################################################################
##################################RRCs & FIREs#################################################
###############################################################################################

FIREs = read.table("/home/loic/Documents/HiC/code/FIREs/FIREs_NEU.txt", header=T)
superFIREs = read.table("/home/loic/Documents/HiC/code/FIREs/NEU_SuperFIREs.txt", header=T)
FIREs.signi = FIREs[FIREs$NEU_indicator==1,]

GRanges.FIREs = GRanges(seqnames = FIREs$chr, ranges=IRanges(start=FIREs$start, end=FIREs$end, names = paste0("FIRE",1:nrow(FIREs))),ScoreFire=FIREs$NEU_neg_ln_pval)
GRanges.FIREs.signi = GRanges(seqnames = FIREs.signi$chr, ranges=IRanges(start=FIREs.signi$start, end=FIREs.signi$end, names = paste0("FIRE",1:nrow(FIREs.signi))), x =FIREs.signi$NEU_neg_ln_pval)
GRanges.superFIREs = GRanges(seqnames = superFIREs$chr, ranges=IRanges(start=superFIREs$start, end=superFIREs$end, names = paste0("superFIRE",1:nrow(superFIREs))))


unique.Enhancers.RRCs = annotate.3D.Features(GRanges.FIREs, unique.Enhancers.RRCs, kind="FIRE", aggregateFunction = "mean")
unique.Promoters.RRCs = annotate.3D.Features(GRanges.FIREs, unique.Promoters.RRCs, kind="FIRE", aggregateFunction = "mean")

distanceRRCs.ABC.FIREs = mcols(distanceToNearest(unique.Enhancers.RRCs,GRanges.FIREs.signi))$distance
distanceRRCs.Rao.FIREs = mcols(distanceToNearest(second(Pairs.prom.regulatory.Rao),GRanges.FIREs.signi))$distance
distanceRRCs.DNAse.FIREs = mcols(distanceToNearest(first(Pairs.prom.regulatory.DNAse),GRanges.FIREs.signi))$distance

boxplot(distanceRRCs.ABC.FIREs[distanceRRCs.ABC.FIREs<1000000],distanceRRCs.Rao.FIREs[distanceRRCs.Rao.FIREs<1000000],distanceRRCs.DNAse.FIREs[distanceRRCs.DNAse.FIREs<1000000])
boxplot(distanceRRCs.Rao.FIREs)
boxplot(distanceRRCs.DNAse.FIREs)

wilcox.test(distanceRRCs.ABC.FIREs,distanceRRCs.Rao.FIREs)
wilcox.test(distanceRRCs.ABC.FIREs,distanceRRCs.DNAse.FIREs)
wilcox.test(distanceRRCs.Rao.FIREs,distanceRRCs.DNAse.FIREs)

table(countOverlaps(unique.Enhancers.RRCs, GRanges.FIREs.signi))/sum(table(countOverlaps(unique.Enhancers.RRCs, GRanges.FIREs.signi)))
table(countOverlaps(second(Pairs.prom.regulatory.Rao), GRanges.FIREs.signi))/sum(table(countOverlaps(second(Pairs.prom.regulatory.Rao), GRanges.FIREs.signi)))
table(countOverlaps(first(Pairs.prom.regulatory.DNAse), GRanges.FIREs.signi))/sum(table(countOverlaps(first(Pairs.prom.regulatory.DNAse), GRanges.FIREs.signi)))

###############################################################################################
#################################RRCs et DI####################################################
###############################################################################################

DI = read.table("DI/all_chrs_dense_annotated.matrix.DI", header=F)
colnames(DI) = c("chr", "start", "end", "DI")
DI$chr = paste0("chr", DI$chr)

DI.WX = DI[DI$chr!="chr23",]
GRanges.DI = GRanges(seqnames = DI.WX$chr, ranges=IRanges(start=DI.WX$start, end=DI.WX$end), DI= DI.WX$DI)

unique.Enhancers.RRCs = annotate.3D.Features(GRanges.DI, unique.Enhancers.RRCs, kind = "DI", aggregateFunction = "mean")
unique.Promoters.RRCs = annotate.3D.Features(GRanges.DI, unique.Promoters.RRCs, kind = "DI", aggregateFunction = "mean")

summary(unique.Enhancers.RRCs$mean_DI)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000    6.079   28.400   66.335   85.013 1483.344 

summary(unique.Promoters.RRCs$mean_DI)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   5.205  21.757  51.480  63.595 970.881 


table(countOverlaps(unique.Enhancers.RRCs, GRanges.DI)) / sum(table(countOverlaps(unique.Enhancers.RRCs, GRanges.DI)))
table(countOverlaps(unique.Promoters.RRCs, GRanges.DI)) / sum(table(countOverlaps(unique.Promoters.RRCs, GRanges.DI)))

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

unique.Enhancers.RRCs = annotate.3D.Features(GRanges.INS, unique.Enhancers.RRCs, kind = "INS", aggregateFunction = "mean")
unique.Promoters.RRCs = annotate.3D.Features(GRanges.INS, unique.Promoters.RRCs, kind = "INS", aggregateFunction = "mean")

summary(unique.Enhancers.RRCs$mean_INS)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.4440  0.7767  0.8092  0.8051  0.8383  0.9437     334 

summary(unique.Promoters.RRCs$mean_INS)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.5174  0.7668  0.7998  0.7957  0.8291  0.9397     181 

#Lien caracteristiques 3D et RRCs:
#Premier niveau d'analyse: Elements inclus dans les RRCs (e.g Enhancers et Promoters)
#Profil de chevauchement:
table(countOverlaps(unique.Enhancers.RRCs, GRanges.INS))
table(countOverlaps(unique.Promoters.RRCs, GRanges.INS))

table(countOverlaps(unique.Enhancers.RRCs, GRanges.DI))
table(countOverlaps(unique.Promoters.RRCs, GRanges.DI))

table(countOverlaps(unique.Enhancers.RRCs, GRanges.FIREs))
table(countOverlaps(unique.Promoters.RRCs, GRanges.FIREs))

#Correlation entre l'activite de la DNAse/H3k27ac/CTCF et les caracteristiques 3D:
###############################################################################################################

unique.Enhancers.RRCs = annotate.Activity(CTCFActivity.NEU, unique.Enhancers.RRCs, kind = "CTCF")
unique.Enhancers.RRCs = annotate.Activity(DNAseActivity.NEU, unique.Enhancers.RRCs, kind="DNAse")
unique.Enhancers.RRCs = annotate.Activity(H3K27acActivity.NEU, unique.Enhancers.RRCs, kind="H3K27ac")
unique.Enhancers.RRCs = annotate.Activity(H3K4me1Activity.NEU, unique.Enhancers.RRCs, kind="H3K4me1")
unique.Enhancers.RRCs =  annotate.Activity(H3K4me3Activity.NEU, unique.Enhancers.RRCs, kind="H3K4me3")
unique.Enhancers.RRCs = annotate.Activity(p300Activity.NEU, unique.Enhancers.RRCs, kind="p300")

FeaturesPlusActivity.Enhancers = data.frame(mcols(unique.Enhancers.RRCs))

unique.Promoters.RRCs = annotate.Activity(CTCFActivity.NEU, unique.Promoters.RRCs, kind = "CTCF")
unique.Promoters.RRCs = annotate.Activity(DNAseActivity.NEU, unique.Promoters.RRCs, kind="DNAse")
unique.Promoters.RRCs = annotate.Activity(H3K27acActivity.NEU, unique.Promoters.RRCs, kind="H3K27ac")
unique.Promoters.RRCs = annotate.Activity(H3K4me1Activity.NEU, unique.Promoters.RRCs, kind="H3K4me1")
unique.Promoters.RRCs =  annotate.Activity(H3K4me3Activity.NEU, unique.Promoters.RRCs, kind="H3K4me3")
unique.Promoters.RRCs = annotate.Activity(p300Activity.NEU, unique.Promoters.RRCs, kind="p300")

FeaturesPlusActivity.Promoters = data.frame(mcols(unique.Promoters.RRCs))
FeaturesPlusActivity.Promoters.tmp = FeaturesPlusActivity.Promoters

FeaturesPlusActivity.Promoters.tmp$geneSymbol = rownames(FeaturesPlusActivity.Promoters.tmp)
FeaturesPlusActivity.Promoters.tmp = merge(FeaturesPlusActivity.Promoters.tmp, expression.NEU[,c("geneSymbol","median.Expr")], by="geneSymbol")
FeaturesPlusActivity.Promoters.tmp$geneSymbol = NULL


relevantCols.Cor.Regul.ABC = FeaturesPlusActivity.Enhancers[,-1]
colnames(relevantCols.Cor.Regul.ABC) = c("Score_FIRE", "DI", "INS", "CTCF", "DNAse", "H3K27ac", "H3K4me1", "H3K4me3","p300", "H3K27me3")
relevantCols.Cor.Prom.ABC = FeaturesPlusActivity.Promoters.tmp[,-1]
colnames(relevantCols.Cor.Prom.ABC) = c("Score_FIRE", "DI", "INS", "CTCF", "DNAse", "H3K27ac", "H3K4me1", "H3K4me3","p300", "H3K27me3", "Expression")

PerformanceAnalytics::chart.Correlation(relevantCols.Cor.Regul.ABC)
PerformanceAnalytics::chart.Correlation(relevantCols.Cor.Prom.ABC)



FeaturesPlusActivity.Promoters.Rao.tmp = prom.Rao[,-1]

FeaturesPlusActivity.Promoters.Rao.tmp = merge(FeaturesPlusActivity.Promoters.Rao.tmp, expression.NEU[,c("geneSymbol","median.Expr")], by="geneSymbol")
FeaturesPlusActivity.Promoters.Rao.tmp$geneSymbol = NULL

relevantCols.Cor.Regul.Rao = regul.Rao[,-c(1,2)]
colnames(relevantCols.Cor.Regul.Rao) = c("DI","Score_FIRE" , "INS", "CTCF", "H3K27ac","DNAse", "H3K4me1", "H3K4me3","p300", "H3K27me3")
colnames(FeaturesPlusActivity.Promoters.Rao.tmp) = c("DI","Score_FIRE" , "INS", "CTCF", "H3K27ac","DNAse", "H3K4me1", "H3K4me3","p300", "H3K27me3", "Expression")


PerformanceAnalytics::chart.Correlation(relevantCols.Cor.Regul.Rao[,c("Score_FIRE","DI","INS", "CTCF", "H3K27ac", "DNAse", "H3K4me1", "H3K4me3", "p300","H3K27me3")])
PerformanceAnalytics::chart.Correlation(FeaturesPlusActivity.Promoters.Rao.tmp[,c("Score_FIRE","DI","INS", "CTCF", "H3K27ac", "DNAse", "H3K4me1", "H3K4me3", "p300","H3K27me3","Expression")])

FeaturesPlusActivity.Promoters.DNAse.tmp = prom.DNAse[,-1]

FeaturesPlusActivity.Promoters.DNAse.tmp = merge(FeaturesPlusActivity.Promoters.DNAse.tmp, expression.NEU[,c("geneSymbol","median.Expr")], by="geneSymbol")
FeaturesPlusActivity.Promoters.DNAse.tmp$geneSymbol = NULL

relevantCols.Cor.Regul.DNAse = regul.DNAse[,-c(1,2)]
colnames(relevantCols.Cor.Regul.DNAse) = c("DI","Score_FIRE" , "INS", "CTCF", "H3K27ac","DNAse", "H3K4me1", "H3K4me3","p300", "H3K27me3")
colnames(FeaturesPlusActivity.Promoters.DNAse.tmp) = c("DI","Score_FIRE" , "INS", "CTCF", "H3K27ac","DNAse", "H3K4me1", "H3K4me3","p300", "H3K27me3", "Expression")


PerformanceAnalytics::chart.Correlation(relevantCols.Cor.Regul.DNAse[,c("Score_FIRE","DI","INS", "CTCF", "H3K27ac", "DNAse", "H3K4me1", "H3K4me3", "p300","H3K27me3")])
PerformanceAnalytics::chart.Correlation(FeaturesPlusActivity.Promoters.DNAse.tmp[,c("Score_FIRE","DI","INS", "CTCF", "H3K27ac", "DNAse", "H3K4me1", "H3K4me3", "p300","H3K27me3","Expression")])

#Correlation au niveau des RRCs:

#Ici on considere la relation entre les caracteristiques 3D pour les elements inclus dans les RRCs
Enhancers.3D = unique.Enhancers.RRCs[!is.na(mcols(unique.Enhancers.RRCs)$mean_ScoreFire)&!is.na(mcols(unique.Enhancers.RRCs)$mean_normalizedScore)&!is.na(mcols(unique.Enhancers.RRCs)$mean_DI)]
Promoters.3D = unique.Promoters.RRCs[!is.na(mcols(unique.Promoters.RRCs)$mean_ScoreFire)&!is.na(mcols(unique.Promoters.RRCs)$mean_normalizedScore)&!is.na(mcols(unique.Promoters.RRCs)$mean_DI)]

df.Promoters.RRCs = data.frame(unique.Promoters.RRCs)
df.Promoters.RRCs = df.Promoters.RRCs[,c("seqnames", "start", "end","ABC.Score" ,"mean_INS", "mean_DI", "mean_ScoreFire", "mean_CTCF", "mean_DNAse","mean_H3K27ac", "mean_H3K4me1", "mean_p300", "mean_H3K4me3")]
colnames(df.Promoters.RRCs) = c("chr.x", "startProm", "endProm","ABC.Score.PROM","INS.PROM", "DI.PROM", "ScoreFire.PROM","mean_CTCF_PROM", "mean_DNAse_PROM","mean_H3K27ac_PROM", "mean_H3K4me1_PROM", "mean_p300_PROM", "mean_H3K4me3_PROM")


df.Enhancers.RRCs = data.frame(unique.Enhancers.RRCs)
df.Enhancers.RRCs = df.Enhancers.RRCs[,c("seqnames", "start", "end", "ABC.Score","mean_INS", "mean_DI", "mean_ScoreFire","mean_CTCF", "mean_DNAse","mean_H3K27ac", "mean_H3K4me1", "mean_p300","mean_H3K4me3")]
colnames(df.Enhancers.RRCs) = c("chr.x", "start", "end","ABC.Score.ENH","INS.ENH", "DI.ENH", "ScoreFire.ENH","mean_CTCF_ENH", "mean_DNAse_ENH","mean_H3K27ac_ENH", "mean_H3K4me1_ENH", "mean_p300_ENH","mean_H3K4me3_ENH")


enhancers.promoters = merge(enhancers.promoters, df.Promoters.RRCs, by=c("chr.x", "startProm", "endProm"))
enhancers.promoters = merge(enhancers.promoters, df.Enhancers.RRCs, by=c("chr.x", "start", "end"))

meanRRCs.ScoreABC = aggregate(ABC.Score~membership, enhancers.promoters, mean, na.rm=T)

#INS
meanRRCsPROM.INS = aggregate(INS.PROM~membership, enhancers.promoters, mean, na.rm=T)
meanRRCsEnh.INS = aggregate(INS.ENH~membership, enhancers.promoters, mean, na.rm=T)

meanRRCs.INS = merge(meanRRCsPROM.INS,meanRRCsEnh.INS)
cor(meanRRCs.INS$INS.PROM, meanRRCs.INS$INS.ENH)
#[1] 0.9170433

meanRRCs.INS$meanINS = apply(meanRRCs.INS[,c("INS.ENH","INS.PROM")], 1, mean, na.rm=T)

#DI
meanRRCsPROM.DI = aggregate(DI.PROM~membership, enhancers.promoters, mean, na.rm=T)
meanRRCsEnh.DI = aggregate(DI.ENH~membership, enhancers.promoters, mean, na.rm=T)

meanRRCs.DI = merge(meanRRCsEnh.DI, meanRRCsPROM.DI)

cor(meanRRCs.DI$DI.ENH, meanRRCs.DI$DI.PROM)
#0.7446607

meanRRCs.DI$meanDI = apply(meanRRCs.DI[,c("DI.ENH","DI.PROM")], 1, mean, na.rm=T)

#Score-FIRE
meanRRCsPROM.FIRE = aggregate(ScoreFire.PROM~membership, enhancers.promoters, mean, na.rm=T)
meanRRCsEnh.FIRE = aggregate(ScoreFire.ENH~membership, enhancers.promoters, mean, na.rm=T)

meanRRCs.FIRE = merge(meanRRCsEnh.FIRE, meanRRCsPROM.FIRE)

cor(meanRRCs.FIRE$ScoreFire.ENH, meanRRCs.FIRE$ScoreFire.PROM)
#0.699876

meanRRCs.FIRE$meanFIRE = apply(meanRRCs.FIRE[,c("ScoreFire.ENH","ScoreFire.PROM")], 1, mean, na.rm=T)

#DNAse
meanRRCsPROM.DNAse = aggregate(mean_DNAse_PROM~membership, enhancers.promoters, mean, na.rm=T)
meanRRCsEnh.DNAse = aggregate(mean_DNAse_ENH~membership, enhancers.promoters, mean, na.rm=T)

meanRRCs.DNAse = merge(meanRRCsEnh.DNAse, meanRRCsPROM.DNAse)

cor(meanRRCs.DNAse$mean_DNAse_ENH, meanRRCs.DNAse$mean_DNAse_PROM)
#0.1562011

meanRRCs.DNAse$meanDNAse = apply(meanRRCs.DNAse[,c("mean_DNAse_ENH","mean_DNAse_PROM")], 1, mean, na.rm=T)

#CTCF
meanRRCsPROM.CTCF = aggregate(mean_CTCF_PROM~membership, enhancers.promoters, mean, na.rm=T)
meanRRCsEnh.CTCF = aggregate(mean_CTCF_ENH~membership, enhancers.promoters, mean, na.rm=T)

meanRRCs.CTCF = merge(meanRRCsPROM.CTCF, meanRRCsEnh.CTCF)

cor(meanRRCs.CTCF$mean_CTCF_PROM,meanRRCs.CTCF$mean_CTCF_ENH)
#0.0932113

meanRRCs.CTCF$meanCTCF = apply(meanRRCs.CTCF[,c("mean_CTCF_PROM","mean_CTCF_ENH")], 1, mean, na.rm=T)

#H3K27ac
meanRRCsPROM.H3K27ac = aggregate(mean_H3K27ac_PROM~membership, enhancers.promoters, mean, na.rm=T)
meanRRCsEnh.H3K27ac = aggregate(mean_H3K27ac_ENH~membership, enhancers.promoters, mean, na.rm=T)

meanRRCs.H3K27ac = merge(meanRRCsPROM.H3K27ac, meanRRCsEnh.H3K27ac)

cor(meanRRCs.H3K27ac$mean_H3K27ac_PROM,meanRRCs.H3K27ac$mean_H3K27ac_ENH)
#0.1453321

meanRRCs.H3K27ac$meanH3K27ac= apply(meanRRCs.H3K27ac[,c("mean_H3K27ac_PROM","mean_H3K27ac_ENH")], 1, mean, na.rm=T)

#H3K27me3
meanRRCsPROM.H3K27me3 = aggregate(mean_H3K27me3_PROM~membership, enhancers.promoters, mean, na.rm=T)
meanRRCsEnh.H3K27me3 = aggregate(mean_H3K27me3_ENH~membership, enhancers.promoters, mean, na.rm=T)

meanRRCs.H3K27me3 = merge(meanRRCsPROM.H3K27me3, meanRRCsEnh.H3K27me3)

cor(meanRRCs.H3K27me3$mean_H3K27me3_PROM,meanRRCs.H3K27me3$mean_H3K27me3_ENH)
#0.31

meanRRCs.H3K27me3$meanH3K27me3= apply(meanRRCs.H3K27me3[,c("mean_H3K27me3_PROM","mean_H3K27me3_ENH")], 1, mean, na.rm=T)


#H3K4me1
meanRRCsPROM.H3K4me1 = aggregate(mean_H3K4me1_PROM~membership, enhancers.promoters, mean, na.rm=T)
meanRRCsEnh.H3K4me1 = aggregate(mean_H3K4me1_ENH~membership, enhancers.promoters, mean, na.rm=T)

meanRRCs.H3K4me1 = merge(meanRRCsPROM.H3K4me1, meanRRCsEnh.H3K4me1)

cor(meanRRCs.H3K4me1$mean_H3K4me1_PROM,meanRRCs.H3K4me1$mean_H3K4me1_ENH)
#1] 0.1782799

meanRRCs.H3K4me1$meanH3K4me1= apply(meanRRCs.H3K4me1[,c("mean_H3K4me1_PROM","mean_H3K4me1_ENH")], 1, mean, na.rm=T)

#p300
meanRRCsPROM.p300 = aggregate(mean_p300_PROM~membership, enhancers.promoters, mean, na.rm=T)
meanRRCsEnh.p300 = aggregate(mean_p300_ENH~membership, enhancers.promoters, mean, na.rm=T)

meanRRCs.p300 = merge(meanRRCsPROM.p300, meanRRCsEnh.p300)

cor(meanRRCs.p300$mean_p300_PROM,meanRRCs.p300$mean_p300_ENH)
#[1] 0.06732281

meanRRCs.p300$meanp300= apply(meanRRCs.p300[,c("mean_p300_PROM","mean_p300_ENH")], 1, mean, na.rm=T)


#H3K4me3
meanRRCsPROM.H3K4me3 = aggregate(mean_H3K4me3_PROM~membership, enhancers.promoters, mean, na.rm=T)
meanRRCsEnh.H3K4me3 = aggregate(mean_H3K4me3_ENH~membership, enhancers.promoters, mean, na.rm=T)

meanRRCs.H3K4me3 = merge(meanRRCsPROM.H3K4me3, meanRRCsEnh.H3K4me3)

cor(meanRRCs.H3K4me3$mean_H3K4me3_PROM, meanRRCs.H3K4me3$mean_H3K4me3_ENH)
#[1] 0.2437115

meanRRCs.H3K4me3$meanH3K4me3= apply(meanRRCs.H3K4me3[,c("mean_H3K4me3_PROM","mean_H3K4me3_ENH")], 1, mean, na.rm=T)

#Expression 
expression.genes.RRCs = data.frame(mcols(expression.genes))
expression.genes.RRCs = expression.genes.RRCs[,c("geneSymbol", "median.Expr")]
colnames(expression.genes.RRCs) = c("TargetGene", "Expression")

enhancers.promoters = merge(enhancers.promoters, expression.genes.RRCs, by="TargetGene")

Full3D = merge(meanRRCs.FIRE[,c("membership","meanFIRE")],merge(meanRRCs.DI[,c("membership","meanDI")], meanRRCs.INS[,c("membership","meanINS")], by="membership"), by="membership")
FullAct = merge(meanRRCs.DNAse[,c("membership","meanDNAse")],merge(meanRRCs.H3K27ac[,c("membership","meanH3K27ac")], meanRRCs.CTCF[,c("membership","meanCTCF")], by="membership"), by="membership")
FullAct = merge(FullAct,merge(meanRRCs.Expre,merge(meanRRCs.H3K4me3[,c("membership", "meanH3K4me3")],merge(meanRRCs.H3K4me1[,c("membership","meanH3K4me1")],meanRRCs.p300[,c("membership", "meanp300")], by="membership"), by = "membership"), by="membership"),by="membership")

Full3DAct= merge(merge(Full3D, FullAct, by="membership"), meanRRCs.ScoreABC, by="membership" )
rownames(Full3DAct) = Full3DAct$membership
Full3DAct$membership = NULL


relevantCols.Cor.RRCs.ABC = Full3DAct[,-c(ncol(Full3DAct))] 
colnames(relevantCols.Cor.RRCs.ABC) = c("FIRE","DI","INS","DNAse","H3K27ac","CTCF","Expression","H3K4me3","H3K4me1","p300")

colnames(Full3DAct.DNAse) = colnames(relevantCols.Cor.RRCs.ABC)
colnames(Full3DAct.Rao) = colnames(relevantCols.Cor.RRCs.ABC)


PerformanceAnalytics::chart.Correlation(relevantCols.Cor.RRCs.ABC, method="pearson")
PerformanceAnalytics::chart.Correlation(Full3DAct.Rao, method="pearson")
PerformanceAnalytics::chart.Correlation(Full3DAct.DNAse, method="pearson")

######################################################################

######################################################################
library(lme4)

unique.Enhancers.RRCs = annotate.Activity(H3K27me3Activity.NEU, unique.Enhancers.RRCs, kind="H3K27me3")
unique.Promoters.RRCs = annotate.Activity(H3K27me3Activity.NEU, unique.Promoters.RRCs, kind="H3K27me3")

H3K27me3.enh = data.frame(unique.Enhancers.RRCs)[,c("seqnames", "start", "end","mean_H3K27me3")]
colnames(H3K27me3.enh) = c("chr.x", "start", "end", "mean_H3K27me3_ENH")
H3K27me3.prom = data.frame(unique.Promoters.RRCs)[,c("seqnames", "start", "end","mean_H3K27me3")]
colnames(H3K27me3.prom) = c("chr.x", "startProm", "endProm", "mean_H3K27me3_PROM")

enhancers.promoters = merge(enhancers.promoters, H3K27me3.prom, by=c("chr.x", "startProm", "endProm"))
enhancers.promoters = merge(enhancers.promoters, H3K27me3.enh, by=c("chr.x", "start", "end"))


enhancers.promoters$mean.Pair.FIREs = apply(enhancers.promoters[,c("ScoreFire.PROM","ScoreFire.ENH")], 1, mean, na.rm=T)
enhancers.promoters$mean.Pair.INS = apply(enhancers.promoters[,c("INS.PROM","INS.ENH")], 1, mean, na.rm=T)
enhancers.promoters$mean.Pair.DI = apply(enhancers.promoters[,c("DI.PROM","DI.ENH")], 1, mean, na.rm =T)
enhancers.promoters$mean.Pair.H3K27ac = apply(enhancers.promoters[,c("mean_H3K27ac_PROM","mean_H3K27ac_ENH")], 1, mean, na.rm=T)
enhancers.promoters$mean.Pair.CTCF = apply(enhancers.promoters[,c("mean_CTCF_PROM","mean_CTCF_ENH")], 1, mean, na.rm=T)
enhancers.promoters$mean.Pair.DNAse = apply(enhancers.promoters[,c("mean_DNAse_PROM","mean_DNAse_ENH")], 1, mean, na.rm=T)
enhancers.promoters$mean.Pair.H3K4me1 = apply(enhancers.promoters[,c("mean_H3K4me1_PROM","mean_H3K4me1_ENH")], 1, mean, na.rm=T)
enhancers.promoters$mean.Pair.H3K4me3 = apply(enhancers.promoters[,c("mean_H3K4me3_PROM","mean_H3K4me3_ENH")], 1, mean, na.rm=T)
enhancers.promoters$mean.Pair.p300 = apply(enhancers.promoters[,c("mean_p300_PROM","mean_p300_ENH")], 1, mean, na.rm=T)
enhancers.promoters$mean.Pair.H3K27me3 = apply(enhancers.promoters[,c("mean_H3K27me3_PROM","mean_H3K27me3_ENH")], 1, mean, na.rm=T)

#Restriction sur les RRCs avec plus d'une connexion:

enhancers.promoters.subset = enhancers.promoters[ave(1:nrow(enhancers.promoters), enhancers.promoters$membership, FUN = length)!=1,]

random.model = function(col2fit, level, data){
  f = as.formula(paste(col2fit,paste0("1|",level), sep="~"))
  return(lmer(f, data, REML = T))
}
ICC = function(model,nested=F){
  vc = data.frame(VarCorr(model))
  if(nested==T){
    return((vc$vcov[1]+vc$vcov[2])/sum(vc$vcov))
  }
  else{
    return(vc$vcov[1]/sum(vc$vcov))
  }
  
}


col2investiguate = c("mean.Pair.FIREs", "mean.Pair.INS", "mean.Pair.DI","mean.Pair.H3K27ac","mean.Pair.CTCF","mean.Pair.DNAse","mean.Pair.H3K4me1","mean.Pair.H3K4me3","mean.Pair.p300")
ICC.ABC = lapply(col2investiguate, 
       function(x) ICC(random.model(x,"membership",enhancers.promoters.subset)))

ABC.ICC=data.frame(cbind(col2investiguate, do.call(rbind, ICC.ABC)))
colnames(ABC.ICC) = c("Metric", "ICC")
ABC.ICC$Method = "ABC"

#On restreint l'analyse sur les RRCs avec au moins deux genes par RRC


ICC(expr.ABC.model)

plot(AL2Genes.ABC[AL2Genes.ABC$membership%in%1:100,"membership"], AL2Genes.ABC[AL2Genes.ABC$membership%in%1:100,"Expression"])
lines(aggregate(Expression~membership,AL2Genes.ABC[AL2Genes.ABC$membership%in%1:100,], mean)$membership,aggregate(Expression~membership,AL2Genes.ABC[AL2Genes.ABC$membership%in%1:100,], mean)$Expression, col="red", lwd=2)
#####################################################################################################

hist(GRanges.FIREs$ScoreFire[GRanges.FIREs$ScoreFire<6])
abline(v=mean(Full3DAct$meanFIRE))

hist(log2(DNAseActivity.NEU$signalValue))
abline(v=mean(Full3DAct$meanDNAse))

hist(log2(H3K27acActivity.NEU$signalValue))
abline(v=mean(Full3DAct$meanH3K27ac))

hist(log2(H3K4me1Activity.NEU$signalValue))
abline(v=mean(Full3DAct$meanH3K4me1))

hist(log2(p300Activity.NEU$signalValue))
abline(v=mean(Full3DAct$meanp300))

hist(log2(H3K4me3Activity.NEU$signalValue))
abline(v=mean(Full3DAct$meanH3K4me3))
################################################################################
#Nombre de genes par RRCs
agg.ABC = aggregate(TargetGene~membership, enhancers.promoters, function(x) length(unique(x)))


#Filtre sur les RRCs avec 1 Gene et ceux ou l'on observe plus de 2 genes
Max1Genes.ABC = enhancers.promoters[enhancers.promoters$membership%in% agg.ABC[agg.ABC$TargetGene==1,"membership"],]
AL2Genes.ABC = enhancers.promoters[enhancers.promoters$membership%in% agg.ABC[agg.ABC$TargetGene>1,"membership"],]

nGenes = aggregate(TargetGene~membership, AL2Genes.ABC, function(x) length(unique(x)))
nEnh = aggregate(name~membership, AL2Genes.ABC, function(x) length(unique(x)))

max.Expression.RRCs = aggregate(Expression~membership, AL2Genes.ABC, function(x) quantile(x, probs = 0.90))
max.H3K27ac.RRCs = aggregate(mean.Pair.H3K27ac~membership, AL2Genes.ABC, function(x) quantile(x, probs = 0.90))
max.CTCF.RRCs = aggregate(mean.Pair.CTCF~membership, AL2Genes.ABC, function(x) quantile(x, probs = 0.90))
max.DNAse.RRCs = aggregate(mean.Pair.DNAse~membership, AL2Genes.ABC, function(x) quantile(x, probs = 0.90))
max.H3K4me1.RRCs = aggregate(mean.Pair.H3K4me1~membership, AL2Genes.ABC, function(x) quantile(x, probs = 0.90))
max.H3K4me3.RRCs = aggregate(mean.Pair.H3K4me3~membership, AL2Genes.ABC, function(x) quantile(x, probs = 0.90))
max.H3K27me3.RRCs = aggregate(mean.Pair.H3K27me3~membership, AL2Genes.ABC, function(x) quantile(x, probs = 0.90))
max.p300.RRCs = aggregate(mean.Pair.p300~membership, AL2Genes.ABC, function(x) quantile(x, probs = 0.90))

#Relation entre nombre de genes et d'elements de regulation et expression maximale dans le RRC
Expression.nElements.RRCs = merge(max.p300.RRCs,merge(max.H3K27me3.RRCs,merge(max.H3K4me3.RRCs,merge(max.H3K4me1.RRCs,merge(max.DNAse.RRCs,merge(max.H3K27ac.RRCs,merge(max.CTCF.RRCs,merge(nGenes, merge(max.Expression.RRCs,nEnh, by="membership"), by="membership"), by="membership"), by="membership"), 
                                  by="membership"),by="membership"), by="membership"),by="membership"),by="membership")

o1 = ggscatter(Expression.nElements.RRCs[Expression.nElements.RRCs$TargetGene<=50,], x="TargetGene", y="Expression", color = "red", add="loess") + stat_cor(label.x=20, label.y = 5, method="spearman",cor.coef.name = "rho")
o1 = ggpar(o1, xlab="#Promoters in RRCs", ylab="Expression")
o2 = ggscatter(Expression.nElements.RRCs[Expression.nElements.RRCs$TargetGene<=50,], x="TargetGene", y="mean.Pair.CTCF", color = "blue", add="loess") + stat_cor(label.x=20, label.y = 5,method="spearman",cor.coef.name = "rho")
o2 = ggpar(o2, xlab="#Promoters in RRCs", ylab="CTCF")
o3 = ggscatter(Expression.nElements.RRCs[Expression.nElements.RRCs$TargetGene<=50,], x="TargetGene", y="mean.Pair.H3K27ac", color = "green", add="loess") + stat_cor(label.x=20, label.y = 4,method="spearman",cor.coef.name = "rho")
o3 = ggpar(o3, xlab="#Promoters in RRCs", ylab="H3K27ac")
o4 = ggscatter(Expression.nElements.RRCs[Expression.nElements.RRCs$TargetGene<=50,], x="TargetGene", y="mean.Pair.DNAse", color = "black", add="loess") + stat_cor(label.x=20, label.y = 4,method="spearman",cor.coef.name = "rho")
o4 = ggpar(o4, xlab="#Promoters in RRCs", ylab="DNAse")
o5 = ggscatter(Expression.nElements.RRCs[Expression.nElements.RRCs$TargetGene<=50,], x="TargetGene", y="mean.Pair.H3K4me1", color = "grey", add="loess") + stat_cor(label.x=20, label.y = 4,method="spearman",cor.coef.name = "rho")
o5 = ggpar(o5, xlab="#Promoters in RRCs", ylab="H3K4me1")
o6 = ggscatter(Expression.nElements.RRCs[Expression.nElements.RRCs$TargetGene<=50,], x="TargetGene", y="mean.Pair.H3K4me3", color = "yellow", add="loess") + stat_cor(label.x=20, label.y = 4,method="spearman",cor.coef.name = "rho")
o6 = ggpar(o6, xlab="#Promoters in RRCs", ylab="H3K4me3")
o7 = ggscatter(Expression.nElements.RRCs[Expression.nElements.RRCs$TargetGene<=50,], x="TargetGene", y="mean.Pair.H3K27me3", color = "purple", add="loess") + stat_cor(label.x=20, label.y = 4,method="spearman",cor.coef.name = "rho")
o7 = ggpar(o7, xlab="#Promoters in RRCs", ylab="H3K27me3")
o8 = ggscatter(Expression.nElements.RRCs[Expression.nElements.RRCs$TargetGene<=50,], x="TargetGene", y="mean.Pair.p300", color = "cyan", add="loess") + stat_cor(label.x=20, label.y = 4,method="spearman",cor.coef.name = "rho")
o8 = ggpar(o8, xlab="#Promoters in RRCs", ylab="p300")
ggarrange(o1,o2,o3,o4,o5,o6,o7,o8, ncol=2, nrow=4, labels = c("A","B","C","D","E","F","G", "H"))

o1 = ggscatter(Expression.nElements.RRCs[Expression.nElements.RRCs$name<=50,], x="name", y="Expression", color = "red", add="loess") + stat_cor(label.x=20, label.y = 4,method="spearman",cor.coef.name = "rho")
o1 = ggpar(o1, xlab="#Regulatory_Elements in RRCs", ylab="Expression")
o2 = ggscatter(Expression.nElements.RRCs[Expression.nElements.RRCs$name<=50,], x="name", y="mean.Pair.CTCF", color = "blue", add="loess") + stat_cor(label.x=20, label.y = 5,method="spearman",cor.coef.name = "rho")
o2 = ggpar(o2, xlab="#Regulatory_Elements in RRCs", ylab="CTCF")
o3 = ggscatter(Expression.nElements.RRCs[Expression.nElements.RRCs$name<=50,], x="name", y="mean.Pair.H3K27ac", color = "green", add="loess") + stat_cor(label.x=20, label.y = 4,method="spearman",cor.coef.name = "rho")
o3 = ggpar(o3, xlab="#Regulatory_Elements in RRCs", ylab="H3K27ac")
o4 = ggscatter(Expression.nElements.RRCs[Expression.nElements.RRCs$name<=50,], x="name", y="mean.Pair.DNAse", color = "black", add="loess") + stat_cor(label.x=20, label.y = 4,method="spearman",cor.coef.name = "rho")
o4 = ggpar(o4, xlab="#Regulatory_Elements in RRCs", ylab="DNAse")
o5 = ggscatter(Expression.nElements.RRCs[Expression.nElements.RRCs$name<=50,], x="name", y="mean.Pair.H3K4me1", color = "grey", add="loess") + stat_cor(label.x=20, label.y = 4,method="spearman",cor.coef.name = "rho")
o5 = ggpar(o5, xlab="#Regulatory_Elements in RRCs", ylab="H3K4me1")
o6 = ggscatter(Expression.nElements.RRCs[Expression.nElements.RRCs$name<=50,], x="name", y="mean.Pair.H3K4me3", color = "yellow", add="loess") + stat_cor(label.x=20, label.y = 4,method="spearman",cor.coef.name = "rho")
o6 = ggpar(o6, xlab="#Regulatory_Elements in RRCs", ylab="H3K4me3")
o7 = ggscatter(Expression.nElements.RRCs[Expression.nElements.RRCs$name<=50,], x="name", y="mean.Pair.H3K27me3", color = "purple", add="loess") + stat_cor(label.x=20, label.y = 4,method="spearman",cor.coef.name = "rho")
o7 = ggpar(o7, xlab="#Regulatory_Elements in RRCs", ylab="H3K27me3")
o8 = ggscatter(Expression.nElements.RRCs[Expression.nElements.RRCs$name<=50,], x="name", y="mean.Pair.p300", color = "cyan", add="loess") + stat_cor(label.x=20, label.y = 4,method="spearman",cor.coef.name = "rho")
o8 = ggpar(o8, xlab="#Regulatory_Elements in RRCs", ylab="p300")
ggarrange(o1,o2,o3,o4,o5,o6,o7,o8, ncol=2, nrow=4, labels = c("A","B","C","D","E","F","G", "H"))

#200 Genes les plus exprimes:
top200Expr.ABC = order(unique(enhancers.promoters[,c("TargetGene","Expression")])[,"Expression"], decreasing = T)[1:200]
membership.mostexpressed.ABC = unique(enhancers.promoters[top200Expr.ABC,"membership"])
SizeG.RRCs.mostexpressed = aggregate(TargetGene~membership,enhancers.promoters[enhancers.promoters$membership%in%membership.mostexpressed.ABC,], function(x) length(unique(x)))

SizeG.RRCs.ALL = aggregate(TargetGene~membership,enhancers.promoters[!enhancers.promoters$membership%in%membership.mostexpressed.ABC,], function(x) length(unique(x)))

mean(SizeG.RRCs.mostexpressed$TargetGene);mean(SizeG.RRCs.ALL$TargetGene)
wilcox.test(SizeG.RRCs.mostexpressed$TargetGene,SizeG.RRCs.ALL$TargetGene)
t.test(SizeG.RRCs.mostexpressed$TargetGene,SizeG.RRCs.ALL$TargetGene, alternative = "greater")


SizeG.RRCs.mostexpressed$MostExpressed = "YES"
SizeG.RRCs.ALL$MostExpressed = "NO"

NGenes.mostExpressed = rbind(SizeG.RRCs.mostexpressed, SizeG.RRCs.ALL)
NGenes.mostExpressed$type="Promoter"

SizeE.RRCs.mostexpressed = aggregate(name~membership,enhancers.promoters[enhancers.promoters$membership%in%membership.mostexpressed.ABC,], function(x) length(unique(x)))
SizeE.RRCs.ALL = aggregate(name~membership,enhancers.promoters[!enhancers.promoters$membership%in%membership.mostexpressed.ABC,], function(x) length(unique(x)))
wilcox.test(SizeE.RRCs.mostexpressed$name,SizeE.RRCs.ALL$name)
t.test(SizeE.RRCs.mostexpressed$name,SizeE.RRCs.ALL$name, alternative = "greater")

SizeE.RRCs.mostexpressed$MostExpressed = "YES"
SizeE.RRCs.ALL$MostExpressed = "NO"

NRegul.mostExpressed = rbind(SizeE.RRCs.mostexpressed, SizeE.RRCs.ALL)
NRegul.mostExpressed$type="Regulatory"
colnames(NGenes.mostExpressed) = c("membership", "N", "MostExpressed", "type")
colnames(NRegul.mostExpressed) = c("membership", "N", "MostExpressed", "type")

mostExpressed.ALL = rbind(NGenes.mostExpressed,NRegul.mostExpressed)

p = ggboxplot(mostExpressed.ALL,x="type", y="N", color="MostExpressed", add="jitter")
p = ggpar(p, legend.title = "Expression Status", ylab="#Elements in RRC", xlab="Element Type")


relationships.mostExpressed.Genes.ABC = enhancers.promoters[enhancers.promoters$TargetGene%in%enhancers.promoters[top200Expr.ABC,"TargetGene"],"TargetGene"]
relationships.mostExpressed.Genes.ABC = droplevels(relationships.mostExpressed.Genes.ABC)

summary(as.numeric(table(relationships.mostExpressed.Genes.ABC)))

relationships.ALL.Genes.ABC = enhancers.promoters[!enhancers.promoters$TargetGene%in%names(relationships.mostExpressed.Genes.ABC),"TargetGene"]
relationships.ALL.Genes.ABC = droplevels(relationships.ALL.Genes.ABC)

wer = data.frame((cbind("YES",table(relationships.mostExpressed.Genes.ABC))))
wfr = data.frame((cbind("NO",table(relationships.ALL.Genes.ABC))))

all.relationships.Genes.ABC = rbind(wer, wfr)
all.relationships.Genes.ABC$X0 = rownames(all.relationships.Genes.ABC)
colnames(all.relationships.Genes.ABC) = c("MostExpressed", "N","TargetGene")
rownames(all.relationships.Genes.ABC) = NULL
all.relationships.Genes.ABC$N = as.numeric(as.character(all.relationships.Genes.ABC$N))
q = ggboxplot(all.relationships.Genes.ABC,x="MostExpressed", y="N", color="MostExpressed", add="jitter")
q = ggpar(q, legend="",ylab="#Connected Elements to gene", xlab="Expression Status")


ggarrange(p+stat_compare_means(method = "t.test", aes(group=MostExpressed), label="p.format") , q+stat_compare_means(method = "t.test", label="p.format",ref.group = "NO", method.args=list(alternative="greater")),
          labels = c("A","B"), ncol=1, nrow=2)
summary(as.numeric(table(relationships.ALL.Genes.ABC)))


wilcox.test(as.numeric(table(relationships.mostExpressed.Genes.ABC)),as.numeric(table(relationships.ALL.Genes.ABC)))
t.test(as.numeric(table(relationships.mostExpressed.Genes.ABC)),as.numeric(table(relationships.ALL.Genes.ABC)), alternative = "greater")

#######################################################################
#Rao

first(Pairs.prom.regulatory.Rao) = annotate.Activity(H3K27me3Activity.NEU, first(Pairs.prom.regulatory.Rao), kind = "H3K27me3")
second(Pairs.prom.regulatory.Rao) = annotate.Activity(H3K27me3Activity.NEU, second(Pairs.prom.regulatory.Rao), kind = "H3K27me3")

H3K27me3.enh.Rao = data.frame(second(Pairs.prom.regulatory.Rao))[,c("seqnames", "start", "end","mean_H3K27me3")]
colnames(H3K27me3.enh.Rao) = c("chr.y", "start.y", "end.y", "mean_H3K27me3.y")
H3K27me3.prom.Rao = data.frame(first(Pairs.prom.regulatory.Rao))[,c("seqnames", "start", "end","mean_H3K27me3")]
colnames(H3K27me3.prom.Rao) = c("chr.x", "start.x", "end.x", "mean_H3K27me3.x")

df.Rao = merge(df.Rao, H3K27me3.prom.Rao, by=c("chr.x", "start.x", "end.x"))
df.Rao = merge(df.Rao, H3K27me3.enh.Rao, by=c("chr.y", "start.y", "end.y"))


nGenes.Rao = aggregate(geneSymbol.x~compo.Rao.membership, df.Rao, function(x) length(unique(x)))
Max1Genes.Rao = df.Rao[df.Rao$compo.Rao.membership %in% nGenes.Rao[nGenes.Rao$geneSymbol.x==1,"compo.Rao.membership"],]
AL2Genes.Rao = df.Rao[df.Rao$compo.Rao.membership %in% nGenes.Rao[nGenes.Rao$geneSymbol.x>1,"compo.Rao.membership"],]

MaxExpression.Rao = aggregate(Expression~compo.Rao.membership,AL2Genes.Rao, function(x) quantile(x, probs = 0.90))
max.H3K27ac.Rao = aggregate(mean.Pair.H3K27ac~compo.Rao.membership, AL2Genes.Rao, function(x) quantile(x, probs = 0.90))
max.CTCF.Rao = aggregate(mean.Pair.CTCF~compo.Rao.membership, AL2Genes.Rao, function(x) quantile(x, probs = 0.90))
max.DNAse.Rao = aggregate(mean.Pair.DNAse~compo.Rao.membership, AL2Genes.Rao, function(x) quantile(x, probs = 0.90))
max.H3K4me1.Rao = aggregate(mean.Pair.H3K4me1~compo.Rao.membership, AL2Genes.Rao, function(x) quantile(x, probs = 0.90))
max.H3K4me3.Rao = aggregate(mean.Pair.H3K4me3~compo.Rao.membership, AL2Genes.Rao, function(x) quantile(x, probs = 0.90))
max.H3K27me3.Rao = aggregate(mean.Pair.H3K27me3~compo.Rao.membership, AL2Genes.Rao, function(x) quantile(x, probs = 0.90))
max.p300.Rao = aggregate(mean.Pair.p300~compo.Rao.membership, AL2Genes.Rao, function(x) quantile(x, probs = 0.90))


nGenes.membership.Rao = aggregate(geneSymbol.x~compo.Rao.membership, AL2Genes.Rao, function(x) length(unique(x)))
nEnhan.membership.Rao = aggregate(name~compo.Rao.membership, AL2Genes.Rao, function(x) length(unique(x)))

Expression.nElements.RRCs.Rao = merge(max.p300.Rao,merge(max.H3K27me3.Rao,merge(max.H3K4me3.Rao,merge(max.H3K4me1.Rao,merge(max.DNAse.Rao,merge(max.CTCF.Rao,merge(max.H3K27ac.Rao,merge(nGenes.membership.Rao, merge(MaxExpression.Rao,nEnhan.membership.Rao, by="compo.Rao.membership"), by="compo.Rao.membership"),
                                                                                                                                                                   by="compo.Rao.membership"),by="compo.Rao.membership"),by="compo.Rao.membership"),by="compo.Rao.membership"),by="compo.Rao.membership"),by="compo.Rao.membership"),by="compo.Rao.membership")
Expression.nElements.RRCs.Rao$N = ifelse(Expression.nElements.RRCs.Rao$geneSymbol.x>4, "5+", as.character(Expression.nElements.RRCs.Rao$geneSymbol.x))
Expression.nElements.RRCs.Rao$N_name = ifelse(Expression.nElements.RRCs.Rao$name>4, "5+", as.character(Expression.nElements.RRCs.Rao$name))
Expression.nElements.RRCs.Rao$N_name = factor(Expression.nElements.RRCs.Rao$N_name, levels=c("1", "2","3","4","5+"))


o1.Rao = ggscatter(Expression.nElements.RRCs.Rao[Expression.nElements.RRCs.Rao$geneSymbol.x<20,], x="geneSymbol.x", y="Expression", color = "red", add="loess") + stat_cor(method="spearman",label.x=5, label.y = 10, cor.coef.name = "rho")
o1.Rao = ggpar(o1.Rao, xlab="#Promoters in RRCs", ylab="Expression")
o2.Rao = ggscatter(Expression.nElements.RRCs.Rao[Expression.nElements.RRCs.Rao$geneSymbol.x<20,], x="geneSymbol.x", y="mean.Pair.CTCF", color = "blue", add="loess") + stat_cor(method="spearman",label.x=5, label.y = 10, cor.coef.name = "rho")
o2.Rao = ggpar(o2.Rao, xlab="#Promoters in RRCs", ylab="CTCF")
o3.Rao = ggscatter(Expression.nElements.RRCs.Rao[Expression.nElements.RRCs.Rao$geneSymbol.x<20,], x="geneSymbol.x", y="mean.Pair.H3K27ac", color = "green", add="loess") + stat_cor(method="spearman",label.x=5, label.y = 3, cor.coef.name = "rho")
o3.Rao = ggpar(o3.Rao, xlab="#Promoters in RRCs", ylab="H3K27ac")
o4.Rao = ggscatter(Expression.nElements.RRCs.Rao[Expression.nElements.RRCs.Rao$geneSymbol.x<20,], x="geneSymbol.x", y="mean.Pair.DNAse", color = "black", add="loess") + stat_cor(method="spearman",label.x=5, label.y = 3, cor.coef.name = "rho")
o4.Rao = ggpar(o4.Rao, xlab="#Promoters in RRCs", ylab="DNAse")
o5.Rao = ggscatter(Expression.nElements.RRCs.Rao[Expression.nElements.RRCs.Rao$geneSymbol.x<20,], x="geneSymbol.x", y="mean.Pair.H3K4me1", color = "grey", add="loess") + stat_cor(method="spearman",label.x=5, label.y = 3, cor.coef.name = "rho")
o5.Rao = ggpar(o5.Rao, xlab="#Promoters in RRCs", ylab="H3K4me1")
o6.Rao = ggscatter(Expression.nElements.RRCs.Rao[Expression.nElements.RRCs.Rao$geneSymbol.x<20,], x="geneSymbol.x", y="mean.Pair.H3K4me3", color = "yellow", add="loess") + stat_cor(method="spearman",label.x=5, label.y = 3, cor.coef.name = "rho")
o6.Rao = ggpar(o6.Rao, xlab="#Promoters in RRCs", ylab="H3K4me3")
o7.Rao = ggscatter(Expression.nElements.RRCs.Rao[Expression.nElements.RRCs.Rao$geneSymbol.x<20,], x="geneSymbol.x", y="mean.Pair.H3K27me3", color = "purple", add="loess") + stat_cor(method="spearman",label.x=5, label.y = 7, cor.coef.name = "rho")
o7.Rao = ggpar(o7.Rao, xlab="#Promoters in RRCs", ylab="H3K27me3")
o8.Rao = ggscatter(Expression.nElements.RRCs.Rao[Expression.nElements.RRCs.Rao$geneSymbol.x<20,], x="geneSymbol.x", y="mean.Pair.p300", color = "cyan", add="loess") + stat_cor(method="spearman",label.x=5, label.y = 7, cor.coef.name = "rho")
o8.Rao = ggpar(o8.Rao, xlab="#Promoters in RRCs", ylab="p300")
ggarrange(o1.Rao,o2.Rao,o3.Rao,o4.Rao,o5.Rao,o6.Rao,o7.Rao,o8.Rao, ncol=2, nrow=4, labels = c("A","B","C","D","E","F","G", "H"))

o1.Rao = ggscatter(Expression.nElements.RRCs.Rao, x="name", y="Expression", color = "red", add="loess") + stat_cor(method="spearman",label.x=10, label.y = 10, cor.coef.name = "rho")
o1.Rao = ggpar(o1.Rao, xlab="#Regulatory-Elements in RRCs", ylab="Expression")
o2.Rao = ggscatter(Expression.nElements.RRCs.Rao, x="name", y="mean.Pair.CTCF", color = "blue", add="loess") + stat_cor(method="spearman",label.x=10, label.y = 10, cor.coef.name = "rho")
o2.Rao = ggpar(o2.Rao, xlab="#Regulatory-Elements in RRCs", ylab="CTCF")
o3.Rao = ggscatter(Expression.nElements.RRCs.Rao, x="name", y="mean.Pair.H3K27ac", color = "green", add="loess") + stat_cor(method="spearman",label.x=10, label.y = 4, cor.coef.name = "rho")
o3.Rao = ggpar(o3.Rao, xlab="#Regulatory-Elements in RRCs", ylab="H3K27ac")
o4.Rao = ggscatter(Expression.nElements.RRCs.Rao, x="name", y="mean.Pair.DNAse", color = "black", add="loess") + stat_cor(method="spearman",label.x=10, label.y = 4, cor.coef.name = "rho")
o4.Rao = ggpar(o4.Rao, xlab="#Regulatory-Elements in RRCs", ylab="DNAse")
o5.Rao = ggscatter(Expression.nElements.RRCs.Rao, x="name", y="mean.Pair.H3K4me1", color = "grey", add="loess") + stat_cor(method="spearman",label.x=10, label.y = 4, cor.coef.name = "rho")
o5.Rao = ggpar(o5.Rao, xlab="#Regulatory-Elements in RRCs", ylab="H3K4me1")
o6.Rao = ggscatter(Expression.nElements.RRCs.Rao, x="name", y="mean.Pair.H3K4me3", color = "yellow", add="loess") + stat_cor(method="spearman",label.x=10, label.y = 4, cor.coef.name = "rho")
o6.Rao = ggpar(o6.Rao, xlab="#Regulatory-Elements in RRCs", ylab="H3K4me3")
o7.Rao = ggscatter(Expression.nElements.RRCs.Rao, x="name", y="mean.Pair.H3K27me3", color = "purple", add="loess") + stat_cor(method="spearman",label.x=10, label.y = 9, cor.coef.name = "rho")
o7.Rao = ggpar(o7.Rao, xlab="#Regulatory-Elements in RRCs", ylab="H3K27me3")
o8.Rao = ggscatter(Expression.nElements.RRCs.Rao, x="name", y="mean.Pair.p300", color = "cyan", add="loess") + stat_cor(method="spearman",label.x=10, label.y = 8, cor.coef.name = "rho")
o8.Rao = ggpar(o8.Rao, xlab="#Regulatory-Elements in RRCs", ylab="p300")
ggarrange(o1.Rao,o2.Rao,o3.Rao,o4.Rao,o5.Rao,o6.Rao,o7.Rao,o8.Rao, ncol=2, nrow=4, labels = c("A","B","C","D","E","F","G", "H"))





r = ggboxplot(Expression.nElements.RRCs.Rao,x="N", y="Expression", color="N", add="jitter") +
  rotate_x_text(angle = 45)+geom_hline(yintercept=mean(Expression.nElements.RRCs.Rao$Expression), linetype=2)
r = ggpar(r, legend="", xlab="#promoters in RRC")


s = ggboxplot(Expression.nElements.RRCs.Rao,x="N_name", y="Expression", color="N_name", add="jitter") +
  rotate_x_text(angle = 45)+geom_hline(yintercept=mean(Expression.nElements.RRCs.Rao$Expression), linetype=2)
s = ggpar(s, legend="", xlab="#regulatory-elements in RRC")

ggarrange(r +  stat_compare_means(method = "anova", label.y = 10) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                     ref.group = "2", method.args=list(alternative="greater")),s +  stat_compare_means(method = "anova", label.y = 10)+stat_compare_means(label = "p.signif", method = "t.test",
                                                                           ref.group = "1", method.args=list(alternative="greater")), labels=c("A","B"), nrow=2, ncol=1)

Expression.nElements.Rao = merge(max.p300.Rao,merge(max.H3K27me3.Rao,merge(max.H3K4me3.Rao,merge(max.H3K4me1.Rao,merge(max.DNAse.Rao,merge(max.H3K27ac.Rao,merge(max.CTCF.Rao,merge(nGenes.Rao, merge(max.Expression.Rao,nEnhan.membership.Rao, by="compo.Rao.membership"), by="compo.Rao.membership"), by="compo.Rao.membership"), by="compo.Rao.membership"), 
                                                                                                                            by="compo.Rao.membership"),by="compo.Rao.membership"), by="compo.Rao.membership"),by="compo.Rao.membership"),by="compo.Rao.membership")

Expression.nElements.Rao$N_name = ifelse(Expression.nElements.Rao$name>4, "5+", as.character(Expression.nElements.Rao$name))
Expression.nElements.Rao$N_name = factor(Expression.nElements.Rao$N_name, levels=c("1", "2","3","4","5+"))

Expression.nElements.Rao$N = ifelse(Expression.nElements.Rao$geneSymbol.x>4, "5+", as.character(Expression.nElements.Rao$geneSymbol.x))
Expression.nElements.Rao$N = factor(Expression.nElements.Rao$N, levels=c("1", "2","3","4","5+"))

oo1 = ggboxplot(Expression.nElements.Rao, x = "N", y="mean.Pair.H3K27ac", color="N") 
oo1 = ggpar(oo1, legend = "", xlab="#promoters", ylab="H3K27ac", ylim=c(0,4))
oo1=oo1 + stat_compare_means(method = "anova", label.y = 3) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                             ref.group = "2", method.args=list(alternative="greater"))
oo2 = ggboxplot(Expression.nElements.Rao, x = "N", y="mean.Pair.H3K27me3", color="N") 
oo2 = ggpar(oo2, legend = "", xlab="#promoters", ylab="H3K27me3",ylim=c(0,7))
oo2=oo2 + stat_compare_means(method = "anova", label.y = 6) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                              ref.group = "2", method.args=list(alternative="greater"))
oo3 = ggboxplot(Expression.nElements.Rao, x = "N", y="mean.Pair.H3K4me1", color="N") 
oo3 = ggpar(oo3, legend = "", xlab="#promoters", ylab="H3K4me1",ylim=c(0,4))
oo3=oo3 + stat_compare_means(method = "anova", label.y = 3) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                              ref.group = "2", method.args=list(alternative="greater"))
oo4 = ggboxplot(Expression.nElements.Rao, x = "N", y="mean.Pair.H3K4me3", color="N") 
oo4 = ggpar(oo4, legend = "", xlab="#promoters", ylab="H3K4me3",ylim=c(0,5))
oo4 = oo4 + stat_compare_means(method = "anova", label.y = 3) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                              ref.group = "2", method.args=list(alternative="greater"))
oo5 = ggboxplot(Expression.nElements.Rao, x = "N", y="mean.Pair.CTCF", color="N") 
oo5 = ggpar(oo5, legend = "", xlab="#promoters", ylab="CTCF",ylim=c(0,12))
oo5 =oo5 + stat_compare_means(method = "anova", label.y = 9) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                              ref.group = "2", method.args=list(alternative="greater"))
oo6 = ggboxplot(Expression.nElements.Rao, x = "N", y="mean.Pair.DNAse", color="N") 
oo6 = ggpar(oo6, legend = "", xlab="#promoters", ylab="DNAse",ylim=c(0,4))
oo6= oo6 + stat_compare_means(method = "anova", label.y = 3) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                              ref.group = "2", method.args=list(alternative="greater"))

oo7 = ggboxplot(Expression.nElements.Rao, x = "N", y="mean.Pair.p300", color="N") 
oo7 = ggpar(oo7, legend = "", xlab="#promoters", ylab="p300",ylim=c(0,7))
oo7=oo7 + stat_compare_means(method = "anova", label.y = 6) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                              ref.group = "2", method.args=list(alternative="greater"))

ggarrange(oo1,oo2,oo3,oo4,oo5,oo6,oo7, labels=c("A","B","C","D","E","F","G"), nrow=4, ncol=2)

oo1 = ggboxplot(Expression.nElements.Rao, x = "N_name", y="mean.Pair.H3K27ac", color="N_name") 
oo1 = ggpar(oo1, legend = "", xlab="#Regulatory_elements", ylab="H3K27ac", ylim=c(0,4))
oo1=oo1 + stat_compare_means(method = "anova", label.y = 3) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                ref.group = "1", method.args=list(alternative="greater"))
oo2 = ggboxplot(Expression.nElements.Rao, x = "N_name", y="mean.Pair.H3K27me3", color="N_name") 
oo2 = ggpar(oo2, legend = "", xlab="#Regulatory_elements", ylab="H3K27me3",ylim=c(0,7))
oo2=oo2 + stat_compare_means(method = "anova", label.y = 6) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                ref.group = "1", method.args=list(alternative="greater"))
oo3 = ggboxplot(Expression.nElements.Rao, x = "N_name", y="mean.Pair.H3K4me1", color="N_name") 
oo3 = ggpar(oo3, legend = "", xlab="#Regulatory_elements", ylab="H3K4me1",ylim=c(0,4))
oo3=oo3 + stat_compare_means(method = "anova", label.y = 3) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                ref.group = "1", method.args=list(alternative="greater"))
oo4 = ggboxplot(Expression.nElements.Rao, x = "N_name", y="mean.Pair.H3K4me3", color="N_name") 
oo4 = ggpar(oo4, legend = "", xlab="#Regulatory_elements", ylab="H3K4me3",ylim=c(0,5))
oo4 = oo4 + stat_compare_means(method = "anova", label.y = 3) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                  ref.group = "1", method.args=list(alternative="greater"))
oo5 = ggboxplot(Expression.nElements.Rao, x = "N_name", y="mean.Pair.CTCF", color="N_name") 
oo5 = ggpar(oo5, legend = "", xlab="#Regulatory_elements", ylab="CTCF",ylim=c(0,12))
oo5 =oo5 + stat_compare_means(method = "anova", label.y = 9) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                 ref.group = "1", method.args=list(alternative="greater"))
oo6 = ggboxplot(Expression.nElements.Rao, x = "N_name", y="mean.Pair.DNAse", color="N_name") 
oo6 = ggpar(oo6, legend = "", xlab="#Regulatory_elements", ylab="DNAse",ylim=c(0,4))
oo6= oo6 + stat_compare_means(method = "anova", label.y = 3) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                 ref.group = "1", method.args=list(alternative="greater"))

oo7 = ggboxplot(Expression.nElements.Rao, x = "N_name", y="mean.Pair.p300", color="N_name") 
oo7 = ggpar(oo7, legend = "", xlab="#Regulatory_elements", ylab="p300",ylim=c(0,7))
oo7=oo7 + stat_compare_means(method = "anova", label.y = 6) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                ref.group = "1", method.args=list(alternative="greater"))

ggarrange(oo1,oo2,oo3,oo4,oo5,oo6,oo7, labels=c("A","B","C","D","E","F","G"), nrow=4, ncol=2)





summary(aov(Expression.nElements.RRCs.Rao$Expression~Expression.nElements.RRCs.Rao$N))
summary(aov(Expression.nElements.RRCs.Rao$Expression~Expression.nElements.RRCs.Rao$N_name))

top200Expr.Rao = order(unique(df.Rao[,c("geneSymbol.x","Expression")])[,"Expression"], decreasing = T)[1:200]
membership.mostexpressed.Rao = unique(df.Rao[top200Expr.Rao,"compo.Rao.membership"])

SizeG.RRCs.mostexpressed.Rao = aggregate(geneSymbol.x~compo.Rao.membership,df.Rao[df.Rao$compo.Rao.membership%in%membership.mostexpressed.Rao,], function(x) length(unique(x)))
SizeG.RRCs.ALL.Rao = aggregate(geneSymbol.x~compo.Rao.membership,df.Rao[!df.Rao$compo.Rao.membership%in%membership.mostexpressed.Rao,],function(x) length(unique(x)))

mean(SizeG.RRCs.mostexpressed.Rao$geneSymbol.x);mean(SizeG.RRCs.ALL.Rao$geneSymbol.x)
wilcox.test(SizeG.RRCs.mostexpressed.Rao$geneSymbol.x,SizeG.RRCs.ALL.Rao$geneSymbol.x)
t.test(SizeG.RRCs.mostexpressed.Rao$geneSymbol.x,SizeG.RRCs.ALL.Rao$geneSymbol.x, alternative = "greater")

SizeG.RRCs.mostexpressed.Rao$MostExpressed = "YES"
SizeG.RRCs.ALL.Rao$MostExpressed = "NO"

NGenes.mostExpressed.Rao = rbind(SizeG.RRCs.mostexpressed.Rao, SizeG.RRCs.ALL.Rao)
NGenes.mostExpressed.Rao$type="Promoter"

SizeE.RRCs.mostexpressed.Rao = aggregate(name~compo.Rao.membership,df.Rao[df.Rao$compo.Rao.membership%in%membership.mostexpressed.Rao,], function(x) length(unique(x)))
SizeE.RRCs.ALL.Rao = aggregate(name~compo.Rao.membership,df.Rao[!df.Rao$compo.Rao.membership%in%membership.mostexpressed.Rao,], function(x) length(unique(x)))
wilcox.test(SizeE.RRCs.mostexpressed.Rao$name,SizeE.RRCs.ALL.Rao$name)
t.test(SizeE.RRCs.mostexpressed.Rao$name,SizeE.RRCs.ALL.Rao$name, alternative = "greater")

SizeE.RRCs.mostexpressed.Rao$MostExpressed = "YES"
SizeE.RRCs.ALL.Rao$MostExpressed = "NO"

NRegul.mostExpressed.Rao = rbind(SizeE.RRCs.mostexpressed.Rao, SizeE.RRCs.ALL.Rao)
NRegul.mostExpressed.Rao$type="Regulatory"
colnames(NGenes.mostExpressed.Rao) = c("membership", "N", "MostExpressed", "type")
colnames(NRegul.mostExpressed.Rao) = c("membership", "N", "MostExpressed", "type")

mostExpressed.ALL.Rao = rbind(NGenes.mostExpressed.Rao,NRegul.mostExpressed.Rao)

v = ggboxplot(mostExpressed.ALL.Rao,x="type", y="N", color="MostExpressed", add="jitter")
v = ggpar(v,legend.title="Expression Status", ylab="#Elements in RRC", xlab="Element Type")


relationships.mostExpressed.Genes.Rao = df.Rao[df.Rao$geneSymbol.x%in%df.Rao[top200Expr.Rao,"geneSymbol.x"],"geneSymbol.x"]
relationships.mostExpressed.Genes.Rao = droplevels(relationships.mostExpressed.Genes.Rao)

summary(as.numeric(table(relationships.mostExpressed.Genes.Rao)))

relationships.ALL.Genes.Rao = df.Rao[!df.Rao$geneSymbol.x%in%names(relationships.mostExpressed.Genes.Rao),"geneSymbol.x"]

wer.Rao = data.frame((cbind("YES",table(relationships.mostExpressed.Genes.Rao))))
wfr.Rao = data.frame((cbind("NO",table(relationships.ALL.Genes.Rao))))

all.relationships.Genes.Rao = rbind(wer.Rao, wfr.Rao)
all.relationships.Genes.Rao$X0 = rownames(all.relationships.Genes.Rao)
colnames(all.relationships.Genes.Rao) = c("MostExpressed", "N","TargetGene")
rownames(all.relationships.Genes.Rao) = NULL
all.relationships.Genes.Rao$N = as.numeric(as.character(all.relationships.Genes.Rao$N))
i = ggboxplot(all.relationships.Genes.Rao,x="MostExpressed", y="N", color="MostExpressed", add="jitter")
i = ggpar(i, legend="", ylab="#Connected Elements to gene", xlab="Expression Status")


ggarrange(v+stat_compare_means(method = "t.test", aes(group=MostExpressed), label="p.format", label.y = 50)  ,i+stat_compare_means(method = "t.test", label="p.format") , labels=c("A","B"), nrow=2, ncol=1)

first(Pairs.prom.regulatory.DNAse) = annotate.Activity(H3K27me3Activity.NEU, first(Pairs.prom.regulatory.DNAse), kind = "H3K27me3")
second(Pairs.prom.regulatory.DNAse) = annotate.Activity(H3K27me3Activity.NEU, second(Pairs.prom.regulatory.DNAse), kind = "H3K27me3")

H3K27me3.enh.DNAse = data.frame(first(Pairs.prom.regulatory.DNAse))[,c("seqnames", "start", "end","mean_H3K27me3")]
colnames(H3K27me3.enh.DNAse) = c("chr.x", "startPeak", "endPeak", "mean_H3K27me3.x")
H3K27me3.prom.DNAse = data.frame(second(Pairs.prom.regulatory.DNAse))[,c("seqnames", "start", "end","mean_H3K27me3")]
colnames(H3K27me3.prom.DNAse) = c("chr.y", "start.y", "end.y", "mean_H3K27me3.y")

df.DNAse = merge(df.DNAse, H3K27me3.prom.DNAse, by=c("chr.y", "start.y", "end.y"))
df.DNAse = merge(df.DNAse, H3K27me3.enh.DNAse, by=c("chr.x", "startPeak", "endPeak"))

MaxExpression.DNAse = aggregate(Expression~compo.DNAse.membership,df.DNAse, max)

nGenes.DNAse = aggregate(geneSymbol.y~compo.DNAse.membership, df.DNAse, function(x) length(unique(x)))
Max1Genes.DNAse = df.DNAse[df.DNAse$compo.DNAse.membership %in% nGenes.DNAse[nGenes.DNAse$geneSymbol.y==1,"compo.DNAse.membership"],]
AL2Genes.DNAse = df.DNAse[df.DNAse$compo.DNAse.membership %in% nGenes.DNAse[nGenes.DNAse$geneSymbol.y>1,"compo.DNAse.membership"],]

MaxExpression.DNAse = aggregate(Expression~compo.DNAse.membership,AL2Genes.DNAse, function(x) quantile(x, probs = 0.9))
max.H3K27ac.DNAse = aggregate(mean.Pair.H3K27ac~compo.DNAse.membership, AL2Genes.DNAse, function(x) quantile(x, probs = 0.9))
max.CTCF.DNAse = aggregate(mean.Pair.CTCF~compo.DNAse.membership, AL2Genes.DNAse, function(x) quantile(x, probs = 0.9))
max.DNAse.DNAse = aggregate(mean.Pair.DNAse~compo.DNAse.membership, AL2Genes.DNAse, function(x) quantile(x, probs = 0.9))
max.H3K4me1.DNAse = aggregate(mean.Pair.H3K4me1~compo.DNAse.membership, AL2Genes.DNAse, function(x) quantile(x, probs = 0.9))
max.H3K4me3.DNAse = aggregate(mean.Pair.H3K4me3~compo.DNAse.membership, AL2Genes.DNAse, function(x) quantile(x, probs = 0.9))
max.H3K27me3.DNAse = aggregate(mean.Pair.H3K27me3~compo.DNAse.membership, AL2Genes.DNAse, function(x) quantile(x, probs = 0.9))
max.p300.DNAse = aggregate(mean.Pair.p300~compo.DNAse.membership, AL2Genes.DNAse, function(x) quantile(x, probs = 0.9))

nGenes.membership.DNAse = aggregate(geneSymbol.y~compo.DNAse.membership, AL2Genes.DNAse, function(x) length(unique(x)))
nEnhan.membership.DNAse = aggregate(name~compo.DNAse.membership, AL2Genes.DNAse, function(x) length(unique(x)))


Expression.nElements.RRCs.DNAse = merge(max.p300.DNAse,merge(max.H3K27me3.DNAse,merge(max.H3K4me3.DNAse,merge(max.H3K4me1.DNAse,merge(max.DNAse.DNAse,merge(max.CTCF.DNAse,merge(max.H3K27ac.DNAse,merge(nGenes.membership.DNAse, merge(MaxExpression.DNAse,nEnhan.membership.DNAse, by="compo.DNAse.membership"), by="compo.DNAse.membership"),
                                                                                                                                                                   by="compo.DNAse.membership"),by="compo.DNAse.membership"),by="compo.DNAse.membership"),by="compo.DNAse.membership"),by="compo.DNAse.membership"),by="compo.DNAse.membership"),by="compo.DNAse.membership")

Expression.nElements.RRCs.DNAse$N = ifelse(Expression.nElements.RRCs.DNAse$geneSymbol.y>4, "5+", as.character(Expression.nElements.RRCs.DNAse$geneSymbol.y))
Expression.nElements.RRCs.DNAse$N_name = ifelse(Expression.nElements.RRCs.DNAse$name>4, "5+", as.character(Expression.nElements.RRCs.DNAse$name))

Expression.nElements.RRCs.DNAse$N =factor(Expression.nElements.RRCs.DNAse$N, levels=c("1","2","3","4","5+"))
Expression.nElements.RRCs.DNAse$N_name =factor(Expression.nElements.RRCs.DNAse$N_name, levels=c("1","2","3","4","5+"))


o1.DNAse = ggscatter(Expression.nElements.RRCs.DNAse, x="geneSymbol.y", y="Expression", color = "red", add="loess") + stat_cor(method="spearman",label.x=6, label.y = 4, cor.coef.name = "rho")
o1.DNAse = ggpar(o1.DNAse, xlab="#Promoters in RRCs", ylab="Expression")
o2.DNAse = ggscatter(Expression.nElements.RRCs.DNAse, x="geneSymbol.y", y="mean.Pair.CTCF", color = "blue", add="loess") + stat_cor(method="spearman",label.x=6, label.y = 5, cor.coef.name = "rho")
o2.DNAse = ggpar(o2.DNAse, xlab="#Promoters in RRCs", ylab="CTCF")
o3.DNAse = ggscatter(Expression.nElements.RRCs.DNAse, x="geneSymbol.y", y="mean.Pair.H3K27ac", color = "green", add="loess") + stat_cor(method="spearman",label.x=6, label.y = 3, cor.coef.name = "rho")
o3.DNAse = ggpar(o3.DNAse, xlab="#Promoters in RRCs", ylab="H3K27ac")
o4.DNAse = ggscatter(Expression.nElements.RRCs.DNAse, x="geneSymbol.y", y="mean.Pair.DNAse", color = "black", add="loess") + stat_cor(method="spearman",label.x=6, label.y = 3, cor.coef.name = "rho")
o4.DNAse = ggpar(o4.DNAse, xlab="#Promoters in RRCs", ylab="DNAse")
o5.DNAse = ggscatter(Expression.nElements.RRCs.DNAse, x="geneSymbol.y", y="mean.Pair.H3K4me1", color = "grey", add="loess") + stat_cor(method="spearman",label.x=6, label.y = 3, cor.coef.name = "rho")
o5.DNAse = ggpar(o5.DNAse, xlab="#Promoters in RRCs", ylab="H3K4me1")
o6.DNAse = ggscatter(Expression.nElements.RRCs.DNAse, x="geneSymbol.y", y="mean.Pair.H3K4me3", color = "yellow", add="loess") + stat_cor(method="spearman",label.x=6, label.y = 3, cor.coef.name = "rho")
o6.DNAse = ggpar(o6.DNAse, xlab="#Promoters in RRCs", ylab="H3K4me3")
o7.DNAse = ggscatter(Expression.nElements.RRCs.DNAse, x="geneSymbol.y", y="mean.Pair.H3K27me3", color = "purple", add="loess") + stat_cor(method="spearman",label.x=6, label.y = 3, cor.coef.name = "rho")
o7.DNAse = ggpar(o7.DNAse, xlab="#Promoters in RRCs", ylab="H3K27me3")
o8.DNAse = ggscatter(Expression.nElements.RRCs.DNAse, x="geneSymbol.y", y="mean.Pair.p300", color = "cyan", add="loess") + stat_cor(method="spearman",label.x=6, label.y = 7, cor.coef.name = "rho")
o8.DNAse = ggpar(o8.DNAse, xlab="#Promoters in RRCs", ylab="p300")
ggarrange(o1.DNAse,o2.DNAse,o3.DNAse,o4.DNAse,o5.DNAse,o6.DNAse,o7.DNAse,o8.DNAse, ncol=2, nrow=4, labels = c("A","B","C","D","E","F","G", "H"))

o1.DNAse = ggscatter(Expression.nElements.RRCs.DNAse[Expression.nElements.RRCs.DNAse$name<30,], x="name", y="Expression", color = "red", add="loess") + stat_cor(method="spearman",label.x=10 ,label.y = 5, cor.coef.name = "rho")
o1.DNAse = ggpar(o1.DNAse, xlab="#Regulatory-Elements in RRCs", ylab="Expression")
o2.DNAse = ggscatter(Expression.nElements.RRCs.DNAse[Expression.nElements.RRCs.DNAse$name<30,], x="name", y="mean.Pair.CTCF", color = "blue", add="loess") + stat_cor(method="spearman",label.x=10, label.y = 5, cor.coef.name = "rho")
o2.DNAse = ggpar(o2.DNAse, xlab="#Regulatory-Elements in RRCs", ylab="CTCF")
o3.DNAse = ggscatter(Expression.nElements.RRCs.DNAse[Expression.nElements.RRCs.DNAse$name<30,], x="name", y="mean.Pair.H3K27ac", color = "green", add="loess") + stat_cor(method="spearman",label.x=10, label.y = 4, cor.coef.name = "rho")
o3.DNAse = ggpar(o3.DNAse, xlab="#Regulatory-Elements in RRCs", ylab="H3K27ac")
o4.DNAse = ggscatter(Expression.nElements.RRCs.DNAse[Expression.nElements.RRCs.DNAse$name<30,], x="name", y="mean.Pair.DNAse", color = "black", add="loess") + stat_cor(method="spearman",label.x=10, label.y = 4, cor.coef.name = "rho")
o4.DNAse = ggpar(o4.DNAse, xlab="#Regulatory-Elements in RRCs", ylab="DNAse")
o5.DNAse = ggscatter(Expression.nElements.RRCs.DNAse[Expression.nElements.RRCs.DNAse$name<30,], x="name", y="mean.Pair.H3K4me1", color = "grey", add="loess") + stat_cor(method="spearman",label.x=10, label.y = 4, cor.coef.name = "rho")
o5.DNAse = ggpar(o5.DNAse, xlab="#Regulatory-Elements in RRCs", ylab="H3K4me1")
o6.DNAse = ggscatter(Expression.nElements.RRCs.DNAse[Expression.nElements.RRCs.DNAse$name<30,], x="name", y="mean.Pair.H3K4me3", color = "yellow", add="loess") + stat_cor(method="spearman",label.x=10, label.y = 7, cor.coef.name = "rho")
o6.DNAse = ggpar(o6.DNAse, xlab="#Regulatory-Elements in RRCs", ylab="H3K4me3")
o7.DNAse = ggscatter(Expression.nElements.RRCs.DNAse[Expression.nElements.RRCs.DNAse$name<30,], x="name", y="mean.Pair.H3K27me3", color = "purple", add="loess") + stat_cor(method="spearman",label.x=10, label.y = 4, cor.coef.name = "rho")
o7.DNAse = ggpar(o7.DNAse, xlab="#Regulatory-Elements in RRCs", ylab="H3K27me3")
o8.DNAse = ggscatter(Expression.nElements.RRCs.DNAse[Expression.nElements.RRCs.DNAse$name<30,], x="name", y="mean.Pair.p300", color = "cyan", add="loess") + stat_cor(method="spearman",label.x=10, label.y = 7, cor.coef.name = "rho")
o8.DNAse = ggpar(o8.DNAse, xlab="#Regulatory-Elements in RRCs", ylab="p300")
ggarrange(o1.DNAse,o2.DNAse,o3.DNAse,o4.DNAse,o5.DNAse,o6.DNAse,o7.DNAse,o8.DNAse, ncol=2, nrow=4, labels = c("A","B","C","D","E","F","G", "H"))


t = ggboxplot(Expression.nElements.RRCs.DNAse,x="N", y="Expression", color="N", add="jitter") +
  rotate_x_text(angle = 45)+geom_hline(yintercept=mean(Expression.nElements.RRCs.DNAse$Expression), linetype=2)
t = ggpar(t, legend="", xlab="#promoters in RRCs", ylab="Expression")


u = ggboxplot(Expression.nElements.RRCs.DNAse,x="N_name", y="Expression", color="N_name", add="jitter") +
  rotate_x_text(angle = 45)+geom_hline(yintercept=mean(Expression.nElements.RRCs.DNAse$Expression), linetype=2)
u = ggpar(u, legend="", xlab="#regulatory-elements in RRCs", ylab="Expression")


ggarrange(t +  stat_compare_means(method = "anova", label.y = 10) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                      ref.group = "2", method.args=list(alternative="greater")),u +  stat_compare_means(method = "anova", label.y = 10)+stat_compare_means(label = "p.signif", method = "t.test",
                                                                           ref.group = "1", method.args=list(alternative="greater")), labels=c("A","B"), nrow=2, ncol=1)


Expression.nElements.DNAse = merge(max.p300.DNAse,merge(max.H3K27me3.DNAse,merge(max.H3K4me3.DNAse,merge(max.H3K4me1.DNAse,merge(max.DNAse.DNAse,merge(max.H3K27ac.DNAse,merge(max.CTCF.DNAse,merge(nGenes.DNAse, merge(max.Expression.DNAse,nEnhan.membership.DNAse, by="compo.DNAse.membership"), by="compo.DNAse.membership"), by="compo.DNAse.membership"), by="compo.DNAse.membership"), 
                                                                                                                       by="compo.DNAse.membership"),by="compo.DNAse.membership"), by="compo.DNAse.membership"),by="compo.DNAse.membership"),by="compo.DNAse.membership")

Expression.nElements.DNAse$N_name = ifelse(Expression.nElements.DNAse$name>4, "5+", as.character(Expression.nElements.DNAse$name))
Expression.nElements.DNAse$N_name = factor(Expression.nElements.DNAse$N_name, levels=c("1", "2","3","4","5+"))

Expression.nElements.DNAse$N = ifelse(Expression.nElements.DNAse$geneSymbol.y>4, "5+", as.character(Expression.nElements.DNAse$geneSymbol.y))
Expression.nElements.DNAse$N = factor(Expression.nElements.DNAse$N, levels=c("1", "2","3","4","5+"))

oo1 = ggboxplot(Expression.nElements.DNAse, x = "N", y="mean.Pair.H3K27ac", color="N") 
oo1 = ggpar(oo1, legend = "", xlab="#promoters", ylab="H3K27ac")
oo1=oo1 + stat_compare_means(method = "anova", label.y = 2.8) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                ref.group = "2", method.args=list(alternative="greater"))
oo2 = ggboxplot(Expression.nElements.DNAse, x = "N", y="mean.Pair.H3K27me3", color="N") 
oo2 = ggpar(oo2, legend = "", xlab="#promoters", ylab="H3K27me3", ylim = c(0,4.2))
oo2=oo2 + stat_compare_means(method = "anova", label.y = 4) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                ref.group = "2", method.args=list(alternative="greater"))
oo3 = ggboxplot(Expression.nElements.DNAse, x = "N", y="mean.Pair.H3K4me1", color="N") 
oo3 = ggpar(oo3, legend = "", xlab="#promoters", ylab="H3K4me1")
oo3=oo3 + stat_compare_means(method = "anova", label.y = 2.8) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                ref.group = "2", method.args=list(alternative="greater"))
oo4 = ggboxplot(Expression.nElements.DNAse, x = "N", y="mean.Pair.H3K4me3", color="N") 
oo4 = ggpar(oo4, legend = "", xlab="#promoters", ylab="H3K4me3")
oo4 = oo4 + stat_compare_means(method = "anova", label.y = 6) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                  ref.group = "2", method.args=list(alternative="greater"))
oo5 = ggboxplot(Expression.nElements.DNAse, x = "N", y="mean.Pair.CTCF", color="N") 
oo5 = ggpar(oo5, legend = "", xlab="#promoters", ylab="CTCF")
oo5 =oo5 + stat_compare_means(method = "anova", label.y = 10) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                 ref.group = "2", method.args=list(alternative="greater"))
oo6 = ggboxplot(Expression.nElements.DNAse, x = "N", y="mean.Pair.DNAse", color="N") 
oo6 = ggpar(oo6, legend = "", xlab="#promoters", ylab="DNAse")
oo6= oo6 + stat_compare_means(method = "anova", label.y = 3.3) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                 ref.group = "2", method.args=list(alternative="greater"))

oo7 = ggboxplot(Expression.nElements.DNAse, x = "N", y="mean.Pair.p300", color="N") 
oo7 = ggpar(oo7, legend = "", xlab="#promoters", ylab="p300")
oo7=oo7 + stat_compare_means(method = "anova", label.y = 7.5) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                ref.group = "2", method.args=list(alternative="greater"))

ggarrange(oo1,oo2,oo3,oo4,oo5,oo6,oo7, labels=c("A","B","C","D","E","F","G"), nrow=4, ncol=2)

oo1 = ggboxplot(Expression.nElements.DNAse, x = "N_name", y="mean.Pair.H3K27ac", color="N_name") 
oo1 = ggpar(oo1, legend = "", xlab="#Regulatory_elements", ylab="H3K27ac", ylim=c(0,3.2))
oo1=oo1 + stat_compare_means(method = "anova", label.y = 2.8) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                  ref.group = "1", method.args=list(alternative="greater"))
oo2 = ggboxplot(Expression.nElements.DNAse, x = "N_name", y="mean.Pair.H3K27me3", color="N_name") 
oo2 = ggpar(oo2, legend = "", xlab="#Regulatory_elements", ylab="H3K27me3", ylim = c(0,4.2))
oo2=oo2 + stat_compare_means(method = "anova", label.y = 4) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                ref.group = "1", method.args=list(alternative="greater"))
oo3 = ggboxplot(Expression.nElements.DNAse, x = "N_name", y="mean.Pair.H3K4me1", color="N_name") 
oo3 = ggpar(oo3, legend = "", xlab="#Regulatory_elements", ylab="H3K4me1",ylim=c(0,4))
oo3=oo3 + stat_compare_means(method = "anova", label.y = 2.8) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                  ref.group = "1", method.args=list(alternative="greater"))
oo4 = ggboxplot(Expression.nElements.DNAse, x = "N_name", y="mean.Pair.H3K4me3", color="N_name") 
oo4 = ggpar(oo4, legend = "", xlab="#Regulatory_elements", ylab="H3K4me3",ylim=c(0,7))
oo4 = oo4 + stat_compare_means(method = "anova", label.y = 6) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                  ref.group = "1", method.args=list(alternative="greater"))
oo5 = ggboxplot(Expression.nElements.DNAse, x = "N_name", y="mean.Pair.CTCF", color="N_name") 
oo5 = ggpar(oo5, legend = "", xlab="#Regulatory_elements", ylab="CTCF",ylim=c(0,12))
oo5 =oo5 + stat_compare_means(method = "anova", label.y = 10) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                  ref.group = "1", method.args=list(alternative="greater"))
oo6 = ggboxplot(Expression.nElements.DNAse, x = "N_name", y="mean.Pair.DNAse", color="N_name") 
oo6 = ggpar(oo6, legend = "", xlab="#Regulatory_elements", ylab="DNAse",ylim=c(0,4))
oo6= oo6 + stat_compare_means(method = "anova", label.y = 3.3) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                   ref.group = "1", method.args=list(alternative="greater"))

oo7 = ggboxplot(Expression.nElements.DNAse, x = "N_name", y="mean.Pair.p300", color="N_name") 
oo7 = ggpar(oo7, legend = "", xlab="#Regulatory_elements", ylab="p300",ylim=c(0,9))
oo7=oo7 + stat_compare_means(method = "anova", label.y = 7.5) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                  ref.group = "1", method.args=list(alternative="greater"))

ggarrange(oo1,oo2,oo3,oo4,oo5,oo6,oo7, labels=c("A","B","C","D","E","F","G"), nrow=4, ncol=2)



top200Expr.DNAse = order(unique(df.DNAse[,c("geneSymbol.y","Expression")])[,"Expression"], decreasing = T)[1:200]
membership.mostexpressed.DNAse = unique(df.DNAse[top200Expr.DNAse,"compo.DNAse.membership"])

SizeG.RRCs.mostexpressed.DNAse = aggregate(geneSymbol.y~compo.DNAse.membership,df.DNAse[df.DNAse$compo.DNAse.membership%in%membership.mostexpressed.DNAse,], function(x) length(unique(x)))
SizeG.RRCs.ALL.DNAse = aggregate(geneSymbol.y~compo.DNAse.membership,df.DNAse[!df.DNAse$compo.DNAse.membership%in%membership.mostexpressed.DNAse,],function(x) length(unique(x)))

mean(SizeG.RRCs.mostexpressed.DNAse$geneSymbol.y);mean(SizeG.RRCs.ALL.DNAse$geneSymbol.y)
wilcox.test(SizeG.RRCs.mostexpressed.DNAse$geneSymbol.y,SizeG.RRCs.ALL.DNAse$geneSymbol.y)
t.test(SizeG.RRCs.mostexpressed.DNAse$geneSymbol.y,SizeG.RRCs.ALL.DNAse$geneSymbol.y, alternative = "greater")

SizeG.RRCs.mostexpressed.DNAse$MostExpressed = "YES"
SizeG.RRCs.ALL.DNAse$MostExpressed = "NO"

NGenes.mostExpressed.DNAse = rbind(SizeG.RRCs.mostexpressed.DNAse, SizeG.RRCs.ALL.DNAse)
NGenes.mostExpressed.DNAse$type="Promoter"

SizeE.RRCs.mostexpressed.DNAse = aggregate(name~compo.DNAse.membership,df.DNAse[df.DNAse$compo.DNAse.membership%in%membership.mostexpressed.DNAse,], function(x) length(unique(x)))
SizeE.RRCs.ALL.DNAse = aggregate(name~compo.DNAse.membership,df.DNAse[!df.DNAse$compo.DNAse.membership%in%membership.mostexpressed.DNAse,], function(x) length(unique(x)))
wilcox.test(SizeE.RRCs.mostexpressed.DNAse$name,SizeE.RRCs.ALL.DNAse$name)
t.test(SizeE.RRCs.mostexpressed.DNAse$name,SizeE.RRCs.ALL.DNAse$name, alternative = "greater")

SizeE.RRCs.mostexpressed.DNAse$MostExpressed = "YES"
SizeE.RRCs.ALL.DNAse$MostExpressed = "NO"

NRegul.mostExpressed.DNAse = rbind(SizeE.RRCs.mostexpressed.DNAse, SizeE.RRCs.ALL.DNAse)
NRegul.mostExpressed.DNAse$type="Regulatory"
colnames(NGenes.mostExpressed.DNAse) = c("membership", "N", "MostExpressed", "type")
colnames(NRegul.mostExpressed.DNAse) = c("membership", "N", "MostExpressed", "type")

mostExpressed.ALL.DNAse = rbind(NGenes.mostExpressed.DNAse,NRegul.mostExpressed.DNAse)

w = ggboxplot(mostExpressed.ALL.DNAse,x="type", y="N", color="MostExpressed", add="jitter")
w = ggpar(w, legend.title="Expression Status",ylab="#Elements in RRC", xlab="Element Type")


relationships.mostExpressed.Genes.DNAse = df.DNAse[df.DNAse$geneSymbol.y%in%df.DNAse[top200Expr.DNAse,"geneSymbol.y"],"geneSymbol.y"]

summary(as.numeric(table(relationships.mostExpressed.Genes.DNAse)))

relationships.ALL.Genes.DNAse = df.DNAse[!df.DNAse$geneSymbol.y%in%names(relationships.mostExpressed.Genes.DNAse),"geneSymbol.y"]

wer.DNAse = data.frame((cbind("YES",table(relationships.mostExpressed.Genes.DNAse))))
wfr.DNAse = data.frame((cbind("NO",table(relationships.ALL.Genes.DNAse))))

all.relationships.Genes.DNAse = rbind(wer.DNAse, wfr.DNAse)
all.relationships.Genes.DNAse$X0 = rownames(all.relationships.Genes.DNAse)
colnames(all.relationships.Genes.DNAse) = c("MostExpressed", "N","TargetGene")
rownames(all.relationships.Genes.DNAse) = NULL
all.relationships.Genes.DNAse$N = as.numeric(as.character(all.relationships.Genes.DNAse$N))
j = ggboxplot(all.relationships.Genes.DNAse,x="MostExpressed", y="N", color="MostExpressed", add="jitter")
j = ggpar(j,legend="",ylab="#Connected Elements to gene", xlab="Expression Status")

ggarrange(w+stat_compare_means(method = "t.test", aes(group=MostExpressed), label="p.format", label.y=30)  ,j+stat_compare_means(method = "t.test", label="p.format") ,
          labels=c("A", "B"), nrow = 2, ncol=1)


#ICC pour chaque structure de RRCs considerees: 

model.global.ABC = lmer(Expression~(1|membership),Expression.Gene.membership.ABC, REML=T)
ICC(model.global.ABC, nested = F)
model.global.Rao = lmer(Expression~(1|compo.Rao.membership),Expression.Gene.membership.Rao, REML=T)
ICC(model.global.Rao, nested = F)
model.global.DNAse = lmer(Expression~(1|compo.DNAse.membership),Expression.Gene.membership.DNAse, REML=T)
ICC(model.global.DNAse, nested = F)

nGenesbyRRCs.ABC = aggregate(TargetGene~membership,Expression.Gene.membership.ABC, length)
nGenesbyRRCs.Rao = aggregate(geneSymbol.x~compo.Rao.membership,Expression.Gene.membership.Rao, length)
nGenesbyRRCs.DNAse = aggregate(geneSymbol.y~compo.DNAse.membership,Expression.Gene.membership.DNAse, length)

ICC.ABC = sapply(sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:60], function(i){
  subset.membership.ABC = nGenesbyRRCs.ABC[nGenesbyRRCs.ABC$TargetGene<=i,"membership"]
  model = lmer(Expression~(1|membership),Expression.Gene.membership.ABC[Expression.Gene.membership.ABC$membership%in%subset.membership.ABC,], REML=T)
  return(ICC(model,nested = F))
})

ICC.Rao = sapply(sort(unique(nGenesbyRRCs.Rao$geneSymbol.x)), function(i){
  subset.membership.Rao = nGenesbyRRCs.Rao[nGenesbyRRCs.Rao$geneSymbol.x<=i,"compo.Rao.membership"]
  model = lmer(Expression~(1|compo.Rao.membership),Expression.Gene.membership.Rao[Expression.Gene.membership.Rao$compo.Rao.membership%in%subset.membership.Rao,], REML=T)
  return(ICC(model,nested = F))
})

ICC.DNAse = sapply(sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)), function(i){
  subset.membership.DNAse = nGenesbyRRCs.DNAse[nGenesbyRRCs.DNAse$geneSymbol.y<=i,"compo.DNAse.membership"]
  model = lmer(Expression~(1|compo.DNAse.membership),Expression.Gene.membership.DNAse[Expression.Gene.membership.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,], REML=T)
  return(ICC(model,nested = F))
})

df.ICC.ABC = cbind("ABC",sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:60],data.frame(ICC.ABC))
df.ICC.Rao = cbind("Rao",sort(unique(nGenesbyRRCs.Rao$geneSymbol.x)),data.frame(ICC.Rao))
df.ICC.DNAse = cbind("DNAse",sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)),data.frame(ICC.DNAse))

colnames(df.ICC.ABC) = c("Method", "N", "ICC")
colnames(df.ICC.Rao)=colnames(df.ICC.ABC)
colnames(df.ICC.DNAse)=colnames(df.ICC.ABC)

df.ICC.all = rbind(df.ICC.ABC,df.ICC.Rao,df.ICC.DNAse)

a1 = ggline(df.ICC.all, x = "N", y="ICC", color="Method")
a1 = ggpar(a1, xlab="#promoters in RRC", title="Expression")

ICC.ABC.CTCF = sapply(sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:60], function(i){
  subset.membership.ABC = nGenesbyRRCs.ABC[nGenesbyRRCs.ABC$TargetGene<=i,"membership"]
  model = lmer(mean.Pair.CTCF~(1|membership),enhancers.promoters[enhancers.promoters$membership%in%subset.membership.ABC,], REML=T)
  return(ICC(model,nested = F))
})

ICC.Rao.CTCF = sapply(sort(unique(nGenesbyRRCs.Rao$geneSymbol.x)), function(i){
  subset.membership.Rao = nGenesbyRRCs.Rao[nGenesbyRRCs.Rao$geneSymbol.x<=i,"compo.Rao.membership"]
  model = lmer(mean.Pair.CTCF~(1|compo.Rao.membership),df.Rao[df.Rao$compo.Rao.membership%in%subset.membership.Rao,], REML=T)
  return(ICC(model,nested = F))
})

ICC.DNAse.CTCF = sapply(sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)), function(i){
  subset.membership.DNAse = nGenesbyRRCs.DNAse[nGenesbyRRCs.DNAse$geneSymbol.y<=i,"compo.DNAse.membership"]
  model = lmer(mean.Pair.CTCF~(1|compo.DNAse.membership),df.DNAse[df.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,], REML=T)
  return(ICC(model,nested = F))
})

df.ICC.ABC.CTCF = cbind("ABC",sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:60],data.frame(ICC.ABC.CTCF))
df.ICC.Rao.CTCF = cbind("Rao",sort(unique(nGenesbyRRCs.Rao$geneSymbol.x)),data.frame(ICC.Rao.CTCF))
df.ICC.DNAse.CTCF = cbind("DNAse",sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)),data.frame(ICC.DNAse.CTCF))

colnames(df.ICC.ABC.CTCF) = c("Method", "N", "ICC")
colnames(df.ICC.Rao.CTCF)=colnames(df.ICC.ABC.CTCF)
colnames(df.ICC.DNAse.CTCF)=colnames(df.ICC.ABC.CTCF)

df.ICC.all.CTCF = rbind(df.ICC.ABC.CTCF,df.ICC.Rao.CTCF,df.ICC.DNAse.CTCF)

a2 = ggline(df.ICC.all.CTCF, x = "N", y="ICC", color="Method")
a2 = ggpar(a2, xlab="#promoters in RRC", title="CTCF Activity")

ICC.ABC.H3K27ac = sapply(sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:60], function(i){
  subset.membership.ABC = nGenesbyRRCs.ABC[nGenesbyRRCs.ABC$TargetGene<=i,"membership"]
  model = lmer(mean.Pair.H3K27ac~(1|membership),enhancers.promoters[enhancers.promoters$membership%in%subset.membership.ABC,], REML=T)
  return(ICC(model,nested = F))
})

ICC.Rao.H3K27ac = sapply(sort(unique(nGenesbyRRCs.Rao$geneSymbol.x)), function(i){
  subset.membership.Rao = nGenesbyRRCs.Rao[nGenesbyRRCs.Rao$geneSymbol.x<=i,"compo.Rao.membership"]
  model = lmer(mean.Pair.H3K27ac~(1|compo.Rao.membership),df.Rao[df.Rao$compo.Rao.membership%in%subset.membership.Rao,], REML=T)
  return(ICC(model,nested = F))
})

ICC.DNAse.H3K27ac = sapply(sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)), function(i){
  subset.membership.DNAse = nGenesbyRRCs.DNAse[nGenesbyRRCs.DNAse$geneSymbol.y<=i,"compo.DNAse.membership"]
  model = lmer(mean.Pair.H3K27ac~(1|compo.DNAse.membership),df.DNAse[df.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,], REML=T)
  return(ICC(model,nested = F))
})

df.ICC.ABC.H3K27ac = cbind("ABC",sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:60],data.frame(ICC.ABC.H3K27ac))
df.ICC.Rao.H3K27ac = cbind("Rao",sort(unique(nGenesbyRRCs.Rao$geneSymbol.x)),data.frame(ICC.Rao.H3K27ac))
df.ICC.DNAse.H3K27ac = cbind("DNAse",sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)),data.frame(ICC.DNAse.H3K27ac))

colnames(df.ICC.ABC.H3K27ac) = c("Method", "N", "ICC")
colnames(df.ICC.Rao.H3K27ac)=colnames(df.ICC.ABC.H3K27ac)
colnames(df.ICC.DNAse.H3K27ac)=colnames(df.ICC.ABC.H3K27ac)

df.ICC.all.H3K27ac = rbind(df.ICC.ABC.H3K27ac,df.ICC.Rao.H3K27ac,df.ICC.DNAse.H3K27ac)

a3 = ggline(df.ICC.all.H3K27ac, x = "N", y="ICC", color="Method")
a3 = ggpar(a3, xlab="#promoters in RRC", title="H3K27ac Activity")

ICC.ABC.H3K4me3 = sapply(sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:60], function(i){
  subset.membership.ABC = nGenesbyRRCs.ABC[nGenesbyRRCs.ABC$TargetGene<=i,"membership"]
  model = lmer(mean.Pair.H3K4me3~(1|membership),enhancers.promoters[enhancers.promoters$membership%in%subset.membership.ABC,], REML=T)
  return(ICC(model,nested = F))
})

ICC.Rao.H3K4me3 = sapply(sort(unique(nGenesbyRRCs.Rao$geneSymbol.x)), function(i){
  subset.membership.Rao = nGenesbyRRCs.Rao[nGenesbyRRCs.Rao$geneSymbol.x<=i,"compo.Rao.membership"]
  model = lmer(mean.Pair.H3K4me3~(1|compo.Rao.membership),df.Rao[df.Rao$compo.Rao.membership%in%subset.membership.Rao,], REML=T)
  return(ICC(model,nested = F))
})

ICC.DNAse.H3K4me3 = sapply(sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)), function(i){
  subset.membership.DNAse = nGenesbyRRCs.DNAse[nGenesbyRRCs.DNAse$geneSymbol.y<=i,"compo.DNAse.membership"]
  model = lmer(mean.Pair.H3K4me3~(1|compo.DNAse.membership),df.DNAse[df.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,], REML=T)
  return(ICC(model,nested = F))
})

df.ICC.ABC.H3K4me3 = cbind("ABC",sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:60],data.frame(ICC.ABC.H3K4me3))
df.ICC.Rao.H3K4me3 = cbind("Rao",sort(unique(nGenesbyRRCs.Rao$geneSymbol.x)),data.frame(ICC.Rao.H3K4me3))
df.ICC.DNAse.H3K4me3 = cbind("DNAse",sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)),data.frame(ICC.DNAse.H3K4me3))

colnames(df.ICC.ABC.H3K4me3) = c("Method", "N", "ICC")
colnames(df.ICC.Rao.H3K4me3)=colnames(df.ICC.ABC.H3K4me3)
colnames(df.ICC.DNAse.H3K4me3)=colnames(df.ICC.ABC.H3K4me3)

df.ICC.all.H3K4me3 = rbind(df.ICC.ABC.H3K4me3,df.ICC.Rao.H3K4me3,df.ICC.DNAse.H3K4me3)
a4 = ggline(df.ICC.all.H3K4me3, x = "N", y="ICC", color="Method")
a4 = ggpar(a4, xlab="#promoters in RRC", title="H3K4me3 Activity")

ICC.ABC.H3K4me1 = sapply(sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:60], function(i){
  subset.membership.ABC = nGenesbyRRCs.ABC[nGenesbyRRCs.ABC$TargetGene<=i,"membership"]
  model = lmer(mean.Pair.H3K4me1~(1|membership),enhancers.promoters[enhancers.promoters$membership%in%subset.membership.ABC,], REML=T)
  return(ICC(model,nested = F))
})

ICC.Rao.H3K4me1 = sapply(sort(unique(nGenesbyRRCs.Rao$geneSymbol.x)), function(i){
  subset.membership.Rao = nGenesbyRRCs.Rao[nGenesbyRRCs.Rao$geneSymbol.x<=i,"compo.Rao.membership"]
  model = lmer(mean.Pair.H3K4me1~(1|compo.Rao.membership),df.Rao[df.Rao$compo.Rao.membership%in%subset.membership.Rao,], REML=T)
  return(ICC(model,nested = F))
})

ICC.DNAse.H3K4me1 = sapply(sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)), function(i){
  subset.membership.DNAse = nGenesbyRRCs.DNAse[nGenesbyRRCs.DNAse$geneSymbol.y<=i,"compo.DNAse.membership"]
  model = lmer(mean.Pair.H3K4me1~(1|compo.DNAse.membership),df.DNAse[df.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,], REML=T)
  return(ICC(model,nested = F))
})

df.ICC.ABC.H3K4me1 = cbind("ABC",sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:60],data.frame(ICC.ABC.H3K4me1))
df.ICC.Rao.H3K4me1 = cbind("Rao",sort(unique(nGenesbyRRCs.Rao$geneSymbol.x)),data.frame(ICC.Rao.H3K4me1))
df.ICC.DNAse.H3K4me1 = cbind("DNAse",sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)),data.frame(ICC.DNAse.H3K4me1))

colnames(df.ICC.ABC.H3K4me1) = c("Method", "N", "ICC")
colnames(df.ICC.Rao.H3K4me1)=colnames(df.ICC.ABC.H3K4me1)
colnames(df.ICC.DNAse.H3K4me1)=colnames(df.ICC.ABC.H3K4me1)

df.ICC.all.H3K4me1 = rbind(df.ICC.ABC.H3K4me1,df.ICC.Rao.H3K4me1,df.ICC.DNAse.H3K4me1)
a5 = ggline(df.ICC.all.H3K4me1, x = "N", y="ICC", color="Method")
a5 = ggpar(a5, xlab="#promoters in RRC", title="H3K4me1 Activity")

ICC.ABC.H3K27me3 = sapply(sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:60], function(i){
  subset.membership.ABC = nGenesbyRRCs.ABC[nGenesbyRRCs.ABC$TargetGene<=i,"membership"]
  model = lmer(mean.Pair.H3K27me3~(1|membership),enhancers.promoters[enhancers.promoters$membership%in%subset.membership.ABC,], REML=T)
  return(ICC(model,nested = F))
})

ICC.Rao.H3K27me3 = sapply(sort(unique(nGenesbyRRCs.Rao$geneSymbol.x)), function(i){
  subset.membership.Rao = nGenesbyRRCs.Rao[nGenesbyRRCs.Rao$geneSymbol.x<=i,"compo.Rao.membership"]
  model = lmer(mean.Pair.H3K27me3~(1|compo.Rao.membership),df.Rao[df.Rao$compo.Rao.membership%in%subset.membership.Rao,], REML=T)
  return(ICC(model,nested = F))
})

ICC.DNAse.H3K27me3 = sapply(sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)), function(i){
  subset.membership.DNAse = nGenesbyRRCs.DNAse[nGenesbyRRCs.DNAse$geneSymbol.y<=i,"compo.DNAse.membership"]
  model = lmer(mean.Pair.H3K27me3~(1|compo.DNAse.membership),df.DNAse[df.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,], REML=T)
  return(ICC(model,nested = F))
})

df.ICC.ABC.H3K27me3 = cbind("ABC",sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:60],data.frame(ICC.ABC.H3K27me3))
df.ICC.Rao.H3K27me3 = cbind("Rao",sort(unique(nGenesbyRRCs.Rao$geneSymbol.x)),data.frame(ICC.Rao.H3K27me3))
df.ICC.DNAse.H3K27me3 = cbind("DNAse",sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)),data.frame(ICC.DNAse.H3K27me3))

colnames(df.ICC.ABC.H3K27me3) = c("Method", "N", "ICC")
colnames(df.ICC.Rao.H3K27me3)=colnames(df.ICC.ABC.H3K27me3)
colnames(df.ICC.DNAse.H3K27me3)=colnames(df.ICC.ABC.H3K27me3)

df.ICC.all.H3K27me3 = rbind(df.ICC.ABC.H3K27me3,df.ICC.Rao.H3K27me3,df.ICC.DNAse.H3K27me3)
a6 = ggline(df.ICC.all.H3K27me3, x = "N", y="ICC", color="Method")
a6 = ggpar(a6, xlab="#promoters in RRC", title="H3K27me3 Activity")

ICC.ABC.DNAse = sapply(sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:60], function(i){
  subset.membership.ABC = nGenesbyRRCs.ABC[nGenesbyRRCs.ABC$TargetGene<=i,"membership"]
  model = lmer(mean.Pair.DNAse~(1|membership),enhancers.promoters[enhancers.promoters$membership%in%subset.membership.ABC,], REML=T)
  return(ICC(model,nested = F))
})

ICC.Rao.DNAse = sapply(sort(unique(nGenesbyRRCs.Rao$geneSymbol.x)), function(i){
  subset.membership.Rao = nGenesbyRRCs.Rao[nGenesbyRRCs.Rao$geneSymbol.x<=i,"compo.Rao.membership"]
  model = lmer(mean.Pair.DNAse~(1|compo.Rao.membership),df.Rao[df.Rao$compo.Rao.membership%in%subset.membership.Rao,], REML=T)
  return(ICC(model,nested = F))
})

ICC.DNAse.DNAse = sapply(sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)), function(i){
  subset.membership.DNAse = nGenesbyRRCs.DNAse[nGenesbyRRCs.DNAse$geneSymbol.y<=i,"compo.DNAse.membership"]
  model = lmer(mean.Pair.DNAse~(1|compo.DNAse.membership),df.DNAse[df.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,], REML=T)
  return(ICC(model,nested = F))
})

df.ICC.ABC.DNAse = cbind("ABC",sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:60],data.frame(ICC.ABC.DNAse))
df.ICC.Rao.DNAse = cbind("Rao",sort(unique(nGenesbyRRCs.Rao$geneSymbol.x)),data.frame(ICC.Rao.DNAse))
df.ICC.DNAse.DNAse = cbind("DNAse",sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)),data.frame(ICC.DNAse.DNAse))

colnames(df.ICC.ABC.DNAse) = c("Method", "N", "ICC")
colnames(df.ICC.Rao.DNAse)=colnames(df.ICC.ABC.DNAse)
colnames(df.ICC.DNAse.DNAse)=colnames(df.ICC.ABC.DNAse)

df.ICC.all.DNAse = rbind(df.ICC.ABC.DNAse,df.ICC.Rao.DNAse,df.ICC.DNAse.DNAse)
a7 = ggline(df.ICC.all.DNAse, x = "N", y="ICC", color="Method")
a7 = ggpar(a7, xlab="#promoters in RRC", title="DNAse Activity")

ICC.ABC.p300 = sapply(sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:60], function(i){
  subset.membership.ABC = nGenesbyRRCs.ABC[nGenesbyRRCs.ABC$TargetGene<=i,"membership"]
  model = lmer(mean.Pair.p300~(1|membership),enhancers.promoters[enhancers.promoters$membership%in%subset.membership.ABC,], REML=T)
  return(ICC(model,nested = F))
})

ICC.Rao.p300 = sapply(sort(unique(nGenesbyRRCs.Rao$geneSymbol.x)), function(i){
  subset.membership.Rao = nGenesbyRRCs.Rao[nGenesbyRRCs.Rao$geneSymbol.x<=i,"compo.Rao.membership"]
  model = lmer(mean.Pair.p300~(1|compo.Rao.membership),df.Rao[df.Rao$compo.Rao.membership%in%subset.membership.Rao,], REML=T)
  return(ICC(model,nested = F))
})

ICC.DNAse.p300 = sapply(sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)), function(i){
  subset.membership.DNAse = nGenesbyRRCs.DNAse[nGenesbyRRCs.DNAse$geneSymbol.y<=i,"compo.DNAse.membership"]
  model = lmer(mean.Pair.p300~(1|compo.DNAse.membership),df.DNAse[df.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,], REML=T)
  return(ICC(model,nested = F))
})

df.ICC.ABC.p300 = cbind("ABC",sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:60],data.frame(ICC.ABC.p300))
df.ICC.Rao.p300 = cbind("Rao",sort(unique(nGenesbyRRCs.Rao$geneSymbol.x)),data.frame(ICC.Rao.p300))
df.ICC.DNAse.p300 = cbind("DNAse",sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)),data.frame(ICC.DNAse.p300))

colnames(df.ICC.ABC.p300) = c("Method", "N", "ICC")
colnames(df.ICC.Rao.p300)=colnames(df.ICC.ABC.p300)
colnames(df.ICC.DNAse.p300)=colnames(df.ICC.ABC.p300)

df.ICC.all.p300 = rbind(df.ICC.ABC.p300,df.ICC.Rao.p300,df.ICC.DNAse.p300)
a8 = ggline(df.ICC.all.p300, x = "N", y="ICC", color="Method")
a8 = ggpar(a7, xlab="#promoters in RRC", title="p300 Activity")

a1
library(ggpubr)
ggarrange(a2,a3,a4,a5,a6,a7,a8, labels=c("A","B","C","D","E","F","G"), nrow=3, ncol=3)
#########################################################################################


Expression.ABC.membership = aggregate(Expression~membership+TargetGene, AL2Genes.ABC, mean)

mean.expression.ABC.membership = aggregate(Expression~membership, Expression.ABC.membership, mean)

mean.expression.ABC.membership$DiffMean = mean.expression.ABC.membership$Expression - mean(Expression.ABC.membership$Expression, na.rm=T)
jhyt = merge(nGenesbyRRCs.ABC,mean.expression.ABC.membership, by="membership")
Var.expression.ABC.membership = aggregate(Expression~membership, Expression.ABC.membership, var)
jhut = merge(nGenesbyRRCs.ABC,Var.expression.ABC.membership, by="membership")

jhyt$TargetGene_Retyp = ifelse(jhyt$TargetGene>14, "15+", as.character(jhyt$TargetGene))
jhyt$TargetGene_Retyp = factor(jhyt$TargetGene_Retyp, levels=c( "2", "3","4","5","6","7","8", "9", "10", "11", "12","13","14","15+"))
jhyt$DiffMean_Squared = jhyt$DiffMean^2

cor.test(jhyt$TargetGene, jhyt$DiffMean^2, method="spearman")

jhut$TargetGene_Retyp = ifelse(jhut$TargetGene>14, "15+", as.character(jhut$TargetGene))
jhut$TargetGene_Retyp = factor(jhut$TargetGene_Retyp, levels=c( "2", "3","4","5","6","7","8", "9", "10", "11", "12","13","14","15+"))

cor.test(jhut$TargetGene, jhut$Expression, method="spearman")

opi = ggboxplot(jhyt, x="TargetGene_Retyp", y ="DiffMean_Squared", color = "TargetGene_Retyp")
opi = ggpar(opi, legend = "", xlab="#promoters in RRCs", ylab="Difference between RRC mean and Population mean")
opi +stat_compare_means(method = "anova", label.y = 20) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                           ref.group = "2", method.args=list(alternative="less"))

opi.var = ggboxplot(jhut, x="TargetGene_Retyp", y ="Expression", color = "TargetGene_Retyp")
opi.var = ggpar(opi.var, legend="", xlab="#promoters in RRCs", ylab="Intra-RRC Expression Variance")
opi.var +stat_compare_means(method = "anova", label.y = 20) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                            ref.group = "2", method.args=list(alternative="greater"))
######################################
#CTCF
mean.CTCF.ABC.membership = aggregate(mean.Pair.CTCF~membership, AL2Genes.ABC, mean)

mean.CTCF.ABC.membership$DiffMean = mean.CTCF.ABC.membership$mean.Pair.CTCF - mean(AL2Genes.ABC$mean.Pair.CTCF, na.rm=T)
jhyt = merge(nGenesbyRRCs.ABC,mean.CTCF.ABC.membership, by="membership")

Var.CTCF.ABC.membership = aggregate(mean.Pair.CTCF~membership, AL2Genes.ABC, var)
jhut = merge(nGenesbyRRCs.ABC,Var.CTCF.ABC.membership, by="membership")

jhyt$TargetGene_Retyp = ifelse(jhyt$TargetGene>14, "15+", as.character(jhyt$TargetGene))
jhyt$TargetGene_Retyp = factor(jhyt$TargetGene_Retyp, levels=c( "2", "3","4","5","6","7","8", "9", "10", "11", "12","13","14","15+"))
jhyt$DiffMean_Squared = jhyt$DiffMean^2

cor.test(jhyt$TargetGene, jhyt$DiffMean_Squared, method="spearman")

jhut$TargetGene_Retyp = ifelse(jhut$TargetGene>14, "15+", as.character(jhut$TargetGene))
jhut$TargetGene_Retyp = factor(jhut$TargetGene_Retyp, levels=c( "2", "3","4","5","6","7","8", "9", "10", "11", "12","13","14","15+"))

cor.test(jhut$TargetGene, jhut$mean.Pair.CTCF, method="spearman")

opi = ggboxplot(jhyt, x="TargetGene_Retyp", y ="DiffMean_Squared", color = "TargetGene_Retyp")
opi = ggpar(opi, legend = "", xlab="#promoters in RRCs", ylab="Difference between RRC mean and Population mean")
opi +stat_compare_means(method = "anova", label.y = 20) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                            ref.group = "2", method.args=list(alternative="less"))

opi.var = ggboxplot(jhut, x="TargetGene_Retyp", y ="mean.Pair.CTCF", color = "TargetGene_Retyp")
opi.var = ggpar(opi.var, legend="", xlab="#promoters in RRCs", ylab="Intra-RRC CTCF Variance")
opi.var +stat_compare_means(method = "anova", label.y = 20) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                ref.group = "2", method.args=list(alternative="greater"))

##############################################################
#H3K27ac

mean.H3K27ac.ABC.membership = aggregate(mean.Pair.H3K27ac~membership, AL2Genes.ABC, mean)

mean.H3K27ac.ABC.membership$DiffMean = mean.H3K27ac.ABC.membership$mean.Pair.H3K27ac - mean(AL2Genes.ABC$mean.Pair.H3K27ac, na.rm=T)
jhyt = merge(nGenesbyRRCs.ABC,mean.H3K27ac.ABC.membership, by="membership")

Var.H3K27ac.ABC.membership = aggregate(mean.Pair.H3K27ac~membership, AL2Genes.ABC, var)
jhut = merge(nGenesbyRRCs.ABC,Var.H3K27ac.ABC.membership, by="membership")

jhyt$TargetGene_Retyp = ifelse(jhyt$TargetGene>14, "15+", as.character(jhyt$TargetGene))
jhyt$TargetGene_Retyp = factor(jhyt$TargetGene_Retyp, levels=c( "2", "3","4","5","6","7","8", "9", "10", "11", "12","13","14","15+"))
jhyt$DiffMean_Squared = jhyt$DiffMean^2

cor.test(jhyt$TargetGene, jhyt$DiffMean_Squared, method="spearman")

jhut$TargetGene_Retyp = ifelse(jhut$TargetGene>14, "15+", as.character(jhut$TargetGene))
jhut$TargetGene_Retyp = factor(jhut$TargetGene_Retyp, levels=c( "2", "3","4","5","6","7","8", "9", "10", "11", "12","13","14","15+"))

cor.test(jhut$TargetGene, jhut$mean.Pair.H3K27ac, method="spearman")

opi = ggboxplot(jhyt, x="TargetGene_Retyp", y ="DiffMean_Squared", color = "TargetGene_Retyp")
opi = ggpar(opi, legend = "", xlab="#promoters in RRCs", ylab="Difference between RRC mean and Population mean")
opi +stat_compare_means(method = "anova", label.y = 0.5) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                            ref.group = "2", method.args=list(alternative="less"))

opi.var = ggboxplot(jhut, x="TargetGene_Retyp", y ="mean.Pair.H3K27ac", color = "TargetGene_Retyp")
opi.var = ggpar(opi.var, legend="", xlab="#promoters in RRCs", ylab="Intra-RRC H3K27ac Variance")
opi.var +stat_compare_means(method = "anova", label.y = 0.5) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                ref.group = "2", method.args=list(alternative="greater"))

#########################################################################################
#H3K4me1
mean.H3K4me1.ABC.membership = aggregate(mean.Pair.H3K4me1~membership, AL2Genes.ABC, mean)

mean.H3K4me1.ABC.membership$DiffMean = mean.H3K4me1.ABC.membership$mean.Pair.H3K4me1 - mean(AL2Genes.ABC$mean.Pair.H3K4me1, na.rm=T)
jhyt = merge(nGenesbyRRCs.ABC,mean.H3K4me1.ABC.membership, by="membership")

Var.H3K4me1.ABC.membership = aggregate(mean.Pair.H3K4me1~membership, AL2Genes.ABC, var)
jhut = merge(nGenesbyRRCs.ABC,Var.H3K4me1.ABC.membership, by="membership")

jhyt$TargetGene_Retyp = ifelse(jhyt$TargetGene>14, "15+", as.character(jhyt$TargetGene))
jhyt$TargetGene_Retyp = factor(jhyt$TargetGene_Retyp, levels=c( "2", "3","4","5","6","7","8", "9", "10", "11", "12","13","14","15+"))
jhyt$DiffMean_Squared = jhyt$DiffMean^2

cor.test(jhyt$TargetGene, jhyt$DiffMean_Squared, method="spearman")

jhut$TargetGene_Retyp = ifelse(jhut$TargetGene>14, "15+", as.character(jhut$TargetGene))
jhut$TargetGene_Retyp = factor(jhut$TargetGene_Retyp, levels=c( "2", "3","4","5","6","7","8", "9", "10", "11", "12","13","14","15+"))

cor.test(jhut$TargetGene, jhut$mean.Pair.H3K4me1, method="spearman")

opi = ggboxplot(jhyt, x="TargetGene_Retyp", y ="DiffMean_Squared", color = "TargetGene_Retyp")
opi = ggpar(opi, legend = "", xlab="#promoters in RRCs", ylab="Difference between RRC mean and Population mean")
opi +stat_compare_means(method = "anova", label.y = 0.5) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                             ref.group = "2", method.args=list(alternative="less"))

opi.var = ggboxplot(jhut, x="TargetGene_Retyp", y ="mean.Pair.H3K4me1", color = "TargetGene_Retyp")
opi.var = ggpar(opi.var, legend="", xlab="#promoters in RRCs", ylab="Intra-RRC H3K4me1 Variance")
opi.var +stat_compare_means(method = "anova", label.y = 0.5) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                 ref.group = "2", method.args=list(alternative="greater"))
#########################################################################################
#H3K4me3
mean.H3K4me3.ABC.membership = aggregate(mean.Pair.H3K4me3~membership, AL2Genes.ABC, mean)

mean.H3K4me3.ABC.membership$DiffMean = mean.H3K4me3.ABC.membership$mean.Pair.H3K4me3 - mean(AL2Genes.ABC$mean.Pair.H3K4me3, na.rm=T)
jhyt = merge(nGenesbyRRCs.ABC,mean.H3K4me3.ABC.membership, by="membership")

Var.H3K4me3.ABC.membership = aggregate(mean.Pair.H3K4me3~membership, AL2Genes.ABC, var)
jhut = merge(nGenesbyRRCs.ABC,Var.H3K4me3.ABC.membership, by="membership")

jhyt$TargetGene_Retyp = ifelse(jhyt$TargetGene>14, "15+", as.character(jhyt$TargetGene))
jhyt$TargetGene_Retyp = factor(jhyt$TargetGene_Retyp, levels=c( "2", "3","4","5","6","7","8", "9", "10", "11", "12","13","14","15+"))
jhyt$DiffMean_Squared = jhyt$DiffMean^2

cor.test(jhyt$TargetGene, jhyt$DiffMean_Squared, method="spearman")

jhut$TargetGene_Retyp = ifelse(jhut$TargetGene>14, "15+", as.character(jhut$TargetGene))
jhut$TargetGene_Retyp = factor(jhut$TargetGene_Retyp, levels=c( "2", "3","4","5","6","7","8", "9", "10", "11", "12","13","14","15+"))

cor.test(jhut$TargetGene, jhut$mean.Pair.H3K4me3, method="spearman")

opi = ggboxplot(jhyt, x="TargetGene_Retyp", y ="DiffMean_Squared", color = "TargetGene_Retyp")
opi = ggpar(opi, legend = "", xlab="#promoters in RRCs", ylab="Difference between RRC mean and Population mean")
opi +stat_compare_means(method = "anova", label.y = 3) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                             ref.group = "2", method.args=list(alternative="less"))

opi.var = ggboxplot(jhut, x="TargetGene_Retyp", y ="mean.Pair.H3K4me3", color = "TargetGene_Retyp")
opi.var = ggpar(opi.var, legend="", xlab="#promoters in RRCs", ylab="Intra-RRC H3K4me3 Variance")
opi.var +stat_compare_means(method = "anova", label.y = 3) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                 ref.group = "2", method.args=list(alternative="greater"))
#########################################################################################
#H3K27me3
mean.H3K27me3.ABC.membership = aggregate(mean.Pair.H3K27me3~membership, AL2Genes.ABC, mean)

mean.H3K27me3.ABC.membership$DiffMean = mean.H3K27me3.ABC.membership$mean.Pair.H3K27me3 - mean(AL2Genes.ABC$mean.Pair.H3K27me3, na.rm=T)
jhyt = merge(nGenesbyRRCs.ABC,mean.H3K27me3.ABC.membership, by="membership")

Var.H3K27me3.ABC.membership = aggregate(mean.Pair.H3K27me3~membership, AL2Genes.ABC, var)
jhut = merge(nGenesbyRRCs.ABC,Var.H3K27me3.ABC.membership, by="membership")

jhyt$TargetGene_Retyp = ifelse(jhyt$TargetGene>14, "15+", as.character(jhyt$TargetGene))
jhyt$TargetGene_Retyp = factor(jhyt$TargetGene_Retyp, levels=c( "2", "3","4","5","6","7","8", "9", "10", "11", "12","13","14","15+"))
jhyt$DiffMean_Squared = jhyt$DiffMean^2

cor.test(jhyt$TargetGene, jhyt$DiffMean_Squared, method="spearman")

jhut$TargetGene_Retyp = ifelse(jhut$TargetGene>14, "15+", as.character(jhut$TargetGene))
jhut$TargetGene_Retyp = factor(jhut$TargetGene_Retyp, levels=c( "2", "3","4","5","6","7","8", "9", "10", "11", "12","13","14","15+"))

cor.test(jhut$TargetGene, jhut$mean.Pair.H3K27me3, method="spearman")

opi = ggboxplot(jhyt, x="TargetGene_Retyp", y ="DiffMean_Squared", color = "TargetGene_Retyp")
opi = ggpar(opi, legend = "", xlab="#promoters in RRCs", ylab="Difference between RRC mean and Population mean")
opi +stat_compare_means(method = "anova", label.y = 0.5) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                             ref.group = "2", method.args=list(alternative="less"))

opi.var = ggboxplot(jhut, x="TargetGene_Retyp", y ="mean.Pair.H3K27me3", color = "TargetGene_Retyp")
opi.var = ggpar(opi.var, legend="", xlab="#promoters in RRCs", ylab="Intra-RRC H3K27me3 Variance")
opi.var +stat_compare_means(method = "anova", label.y = 0.5) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                 ref.group = "2", method.args=list(alternative="greater"))
#########################################################################################
#DNAse
mean.DNAse.ABC.membership = aggregate(mean.Pair.DNAse~membership, AL2Genes.ABC, mean)

mean.DNAse.ABC.membership$DiffMean = mean.DNAse.ABC.membership$mean.Pair.DNAse - mean(AL2Genes.ABC$mean.Pair.DNAse, na.rm=T)
jhyt = merge(nGenesbyRRCs.ABC,mean.DNAse.ABC.membership, by="membership")

Var.DNAse.ABC.membership = aggregate(mean.Pair.DNAse~membership, AL2Genes.ABC, var)
jhut = merge(nGenesbyRRCs.ABC,Var.DNAse.ABC.membership, by="membership")

jhyt$TargetGene_Retyp = ifelse(jhyt$TargetGene>14, "15+", as.character(jhyt$TargetGene))
jhyt$TargetGene_Retyp = factor(jhyt$TargetGene_Retyp, levels=c( "2", "3","4","5","6","7","8", "9", "10", "11", "12","13","14","15+"))
jhyt$DiffMean_Squared = jhyt$DiffMean^2

cor.test(jhyt$TargetGene, jhyt$DiffMean_Squared, method="spearman")

jhut$TargetGene_Retyp = ifelse(jhut$TargetGene>14, "15+", as.character(jhut$TargetGene))
jhut$TargetGene_Retyp = factor(jhut$TargetGene_Retyp, levels=c( "2", "3","4","5","6","7","8", "9", "10", "11", "12","13","14","15+"))

cor.test(jhut$TargetGene, jhut$mean.Pair.DNAse, method="spearman")

opi = ggboxplot(jhyt, x="TargetGene_Retyp", y ="DiffMean_Squared", color = "TargetGene_Retyp")
opi = ggpar(opi, legend = "", xlab="#promoters in RRCs", ylab="Difference between RRC mean and Population mean")
opi +stat_compare_means(method = "anova", label.y = 0.5) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                             ref.group = "2", method.args=list(alternative="less"))

opi.var = ggboxplot(jhut, x="TargetGene_Retyp", y ="mean.Pair.DNAse", color = "TargetGene_Retyp")
opi.var = ggpar(opi.var, legend="", xlab="#promoters in RRCs", ylab="Intra-RRC DNAse Variance")
opi.var +stat_compare_means(method = "anova", label.y = 0.5) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                 ref.group = "2", method.args=list(alternative="greater"))
#########################################################################################
#p300
mean.p300.ABC.membership = aggregate(mean.Pair.p300~membership, AL2Genes.ABC, mean)

mean.p300.ABC.membership$DiffMean = mean.p300.ABC.membership$mean.Pair.p300 - mean(AL2Genes.ABC$mean.Pair.p300, na.rm=T)
jhyt = merge(nGenesbyRRCs.ABC,mean.p300.ABC.membership, by="membership")

Var.p300.ABC.membership = aggregate(mean.Pair.p300~membership, AL2Genes.ABC, var)
jhut = merge(nGenesbyRRCs.ABC,Var.p300.ABC.membership, by="membership")

jhyt$TargetGene_Retyp = ifelse(jhyt$TargetGene>14, "15+", as.character(jhyt$TargetGene))
jhyt$TargetGene_Retyp = factor(jhyt$TargetGene_Retyp, levels=c( "2", "3","4","5","6","7","8", "9", "10", "11", "12","13","14","15+"))
jhyt$DiffMean_Squared = jhyt$DiffMean^2

cor.test(jhyt$TargetGene, jhyt$DiffMean_Squared, method="spearman")

jhut$TargetGene_Retyp = ifelse(jhut$TargetGene>14, "15+", as.character(jhut$TargetGene))
jhut$TargetGene_Retyp = factor(jhut$TargetGene_Retyp, levels=c( "2", "3","4","5","6","7","8", "9", "10", "11", "12","13","14","15+"))

cor.test(jhut$TargetGene, jhut$mean.Pair.p300, method="spearman")

opi = ggboxplot(jhyt, x="TargetGene_Retyp", y ="DiffMean_Squared", color = "TargetGene_Retyp")
opi = ggpar(opi, legend = "", xlab="Cluster Size", ylab="Difference between RRC mean and Population mean")
opi +stat_compare_means(method = "anova", label.y = 2) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                             ref.group = "2", method.args=list(alternative="less"))

opi.var = ggboxplot(jhut, x="TargetGene_Retyp", y ="mean.Pair.p300", color = "TargetGene_Retyp")
opi.var = ggpar(opi.var, legend="", xlab="Cluster Size", ylab="Intra-RRC p300 Variance")
opi.var +stat_compare_means(method = "anova", label.y = 2) +stat_compare_means(label = "p.signif", method = "t.test",
                                                                                 ref.group = "2", method.args=list(alternative="greater"))

################################################3
#Detection de Phenomenes de compensation entre marques d'histones contradictoires 
#Consideration des marques d'histone H3K27me3/H3K27ac/H3K4me1/H3K4me3/CTCF

max.Expression.RRCs = aggregate(Expression~membership, AL2Genes.ABC, max)
max.H3K27ac.RRCs = aggregate(mean.Pair.H3K27ac~membership, AL2Genes.ABC, max)
max.CTCF.RRCs = aggregate(mean.Pair.CTCF~membership, AL2Genes.ABC, max)
max.DNAse.RRCs = aggregate(mean.Pair.DNAse~membership, AL2Genes.ABC, max)
max.H3K4me1.RRCs = aggregate(mean.Pair.H3K4me1~membership, AL2Genes.ABC, max)
max.H3K4me3.RRCs = aggregate(mean.Pair.H3K4me3~membership, AL2Genes.ABC, max)
max.H3K27me3.RRCs = aggregate(mean.Pair.H3K27me3~membership, AL2Genes.ABC, max)
max.p300.RRCs = aggregate(mean.Pair.p300~membership, AL2Genes.ABC, max)

Expr.H3K27me3.ABC = merge(max.Expression.RRCs, max.H3K27me3.RRCs, by="membership")
H3K27ac.H3K27me3.ABC = merge(max.H3K27ac.RRCs, max.H3K27me3.RRCs, by="membership")
H3K4me1.H3K27me3.ABC = merge(max.H3K4me1.RRCs, max.H3K27me3.RRCs, by="membership")
H3K4me3.H3K27me3.ABC = merge(max.H3K4me3.RRCs, max.H3K27me3.RRCs, by="membership")
H3K27ac.CTCF.ABC = merge(max.H3K27ac.RRCs, max.CTCF.RRCs,by="membership")
H3K27me3.CTCF.ABC = merge(max.H3K27me3.RRCs, max.CTCF.RRCs,by="membership")
H3K4me1.CTCF.ABC = merge(max.H3K4me1.RRCs,max.CTCF.RRCs, by="membership")
H3K4me3.CTCF.ABC = merge(max.H3K4me3.RRCs,max.CTCF.RRCs, by="membership")
p300.CTCF.ABC = merge(max.p300.RRCs,max.CTCF.RRCs, by="membership")
H3K27me3.p300.ABC = merge(max.H3K27me3.RRCs,max.p300.RRCs, by="membership")

comp = sapply(sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:25], function(i){
  subset.membership.ABC = nGenesbyRRCs.ABC[nGenesbyRRCs.ABC$TargetGene<=i,"membership"]
  return(cor(Expr.H3K27me3.ABC[Expr.H3K27me3.ABC$membership%in%subset.membership.ABC,"mean.Pair.H3K27me3"],Expr.H3K27me3.ABC[Expr.H3K27me3.ABC$membership%in%subset.membership.ABC,"Expression"]))
  
})

comp= data.frame(comp)
colnames(comp) = "COR"
comp$Variables = "Expression_H3K27me3"
comp$N = sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:25]


comp1 = sapply(sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:25], function(i){
  subset.membership.ABC = nGenesbyRRCs.ABC[nGenesbyRRCs.ABC$TargetGene<=i,"membership"]
  return(cor(H3K27ac.H3K27me3.ABC[H3K27ac.H3K27me3.ABC$membership%in%subset.membership.ABC,"mean.Pair.H3K27me3"],H3K27ac.H3K27me3.ABC[H3K27ac.H3K27me3.ABC$membership%in%subset.membership.ABC,"mean.Pair.H3K27ac"]))
  
})

comp1= data.frame(comp1)
colnames(comp1) = "COR"
comp1$Variables = "H3K27ac_H3K27me3"
comp1$N = sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:25]

comp2 = sapply(sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:25], function(i){
  subset.membership.ABC = nGenesbyRRCs.ABC[nGenesbyRRCs.ABC$TargetGene<=i,"membership"]
  return(cor(H3K4me1.H3K27me3.ABC[H3K4me1.H3K27me3.ABC$membership%in%subset.membership.ABC,"mean.Pair.H3K27me3"],H3K4me1.H3K27me3.ABC[H3K4me1.H3K27me3.ABC$membership%in%subset.membership.ABC,"mean.Pair.H3K4me1"]))
  
})
comp2= data.frame(comp2)
colnames(comp2) = "COR"
comp2$Variables = "H3K4me1_H3K27me3"
comp2$N = sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:25]

comp3 = sapply(sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:25], function(i){
  subset.membership.ABC = nGenesbyRRCs.ABC[nGenesbyRRCs.ABC$TargetGene<=i,"membership"]
  return(cor(H3K4me3.H3K27me3.ABC[H3K4me3.H3K27me3.ABC$membership%in%subset.membership.ABC,"mean.Pair.H3K27me3"],H3K4me3.H3K27me3.ABC[H3K4me3.H3K27me3.ABC$membership%in%subset.membership.ABC,"mean.Pair.H3K4me3"]))
  
})
comp3= data.frame(comp3)
colnames(comp3) = "COR"
comp3$Variables = "H3K4me3_H3K27me3"
comp3$N = sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:25]

comp4 = sapply(sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:25], function(i){
  subset.membership.ABC = nGenesbyRRCs.ABC[nGenesbyRRCs.ABC$TargetGene<=i,"membership"]
  return(cor(H3K27ac.CTCF.ABC[H3K27ac.CTCF.ABC$membership%in%subset.membership.ABC,"mean.Pair.CTCF"],H3K27ac.CTCF.ABC[H3K27ac.CTCF.ABC$membership%in%subset.membership.ABC,"mean.Pair.H3K27ac"]))
  
})
comp4= data.frame(comp4)
colnames(comp4) = "COR"
comp4$Variables = "CTCF_H3K27ac"
comp4$N = sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:25]

comp5 = sapply(sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:25], function(i){
  subset.membership.ABC = nGenesbyRRCs.ABC[nGenesbyRRCs.ABC$TargetGene<=i,"membership"]
  return(cor(H3K27me3.CTCF.ABC[H3K27me3.CTCF.ABC$membership%in%subset.membership.ABC,"mean.Pair.H3K27me3"],H3K27me3.CTCF.ABC[H3K27me3.CTCF.ABC$membership%in%subset.membership.ABC,"mean.Pair.CTCF"]))
  
})
comp5= data.frame(comp5)
colnames(comp5) = "COR"
comp5$Variables = "CTCF_H3K27me3"
comp5$N = sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:25]

comp6 = sapply(sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:25], function(i){
  subset.membership.ABC = nGenesbyRRCs.ABC[nGenesbyRRCs.ABC$TargetGene<=i,"membership"]
  return(cor(H3K4me1.CTCF.ABC[H3K4me1.CTCF.ABC$membership%in%subset.membership.ABC,"mean.Pair.H3K4me1"],H3K4me1.CTCF.ABC[H3K4me1.CTCF.ABC$membership%in%subset.membership.ABC,"mean.Pair.CTCF"]))
  
})
comp6= data.frame(comp6)
colnames(comp6) = "COR"
comp6$Variables = "CTCF_H3K4me1"
comp6$N = sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:25]

comp7 = sapply(sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:25], function(i){
  subset.membership.ABC = nGenesbyRRCs.ABC[nGenesbyRRCs.ABC$TargetGene<=i,"membership"]
  return(cor(H3K4me3.CTCF.ABC[H3K4me3.CTCF.ABC$membership%in%subset.membership.ABC,"mean.Pair.H3K4me3"],H3K4me3.CTCF.ABC[H3K4me3.CTCF.ABC$membership%in%subset.membership.ABC,"mean.Pair.CTCF"]))
  
})
comp7= data.frame(comp7)
colnames(comp7) = "COR"
comp7$Variables = "CTCF_H3K4me3"
comp7$N = sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:25]

comp8 = sapply(sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:25], function(i){
  subset.membership.ABC = nGenesbyRRCs.ABC[nGenesbyRRCs.ABC$TargetGene<=i,"membership"]
  return(cor(p300.CTCF.ABC[p300.CTCF.ABC$membership%in%subset.membership.ABC,"mean.Pair.p300"],p300.CTCF.ABC[p300.CTCF.ABC$membership%in%subset.membership.ABC,"mean.Pair.CTCF"]))
  
})

comp8= data.frame(comp8)
colnames(comp8) = "COR"
comp8$Variables = "CTCF_p300"
comp8$N = sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:25]

comp9 = sapply(sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:25], function(i){
  subset.membership.ABC = nGenesbyRRCs.ABC[nGenesbyRRCs.ABC$TargetGene<=i,"membership"]
  return(cor(H3K27me3.p300.ABC[H3K27me3.p300.ABC$membership%in%subset.membership.ABC,"mean.Pair.p300"],H3K27me3.p300.ABC[H3K27me3.p300.ABC$membership%in%subset.membership.ABC,"mean.Pair.H3K27me3"]))
  
})

comp9= data.frame(comp9)
colnames(comp9) = "COR"
comp9$Variables = "H3K27me3_p300"
comp9$N = sort(unique(nGenesbyRRCs.ABC$TargetGene))[1:25]


comp.tt = rbind(comp,comp1,comp2,comp3,comp4,comp9,comp5,comp6,comp7,comp8)

ggl = ggline(comp.tt, x="N", y="COR", color="Variables", plot_type = "l")
ggl = ggpar(ggl, legend.title = "Association", xlab="#promoters in RRC", ylab="Correlation")
ggl

max.Expression.Rao = aggregate(Expression~compo.Rao.membership, AL2Genes.Rao, max)
max.H3K27ac.Rao = aggregate(mean.Pair.H3K27ac~compo.Rao.membership, AL2Genes.Rao, max)
max.CTCF.Rao = aggregate(mean.Pair.CTCF~compo.Rao.membership, AL2Genes.Rao, max)
max.DNAse.Rao = aggregate(mean.Pair.DNAse~compo.Rao.membership, AL2Genes.Rao, max)
max.H3K4me1.Rao = aggregate(mean.Pair.H3K4me1~compo.Rao.membership, AL2Genes.Rao, max)
max.H3K4me3.Rao = aggregate(mean.Pair.H3K4me3~compo.Rao.membership, AL2Genes.Rao, max)
max.H3K27me3.Rao = aggregate(mean.Pair.H3K27me3~compo.Rao.membership, AL2Genes.Rao, max)
max.p300.Rao = aggregate(mean.Pair.p300~compo.Rao.membership, AL2Genes.Rao, max)

Expre.H3K27me3.Rao = merge(max.Expression.Rao, max.H3K27me3.Rao, by="compo.Rao.membership")
H3K27ac.H3K27me3.Rao = merge(max.H3K27ac.Rao, max.H3K27me3.Rao, by="compo.Rao.membership")
H3K4me1.H3K27me3.Rao = merge(max.H3K4me1.Rao, max.H3K27me3.Rao, by="compo.Rao.membership")
H3K4me3.H3K27me3.Rao = merge(max.H3K4me3.Rao, max.H3K27me3.Rao, by="compo.Rao.membership")
H3K27ac.CTCF.Rao = merge(max.H3K27ac.Rao, max.CTCF.Rao,by="compo.Rao.membership")
H3K27me3.CTCF.Rao = merge(max.H3K27me3.Rao, max.CTCF.Rao,by="compo.Rao.membership")
H3K4me1.CTCF.Rao = merge(max.H3K4me1.Rao,max.CTCF.Rao, by="compo.Rao.membership")
H3K4me3.CTCF.Rao = merge(max.H3K4me3.Rao,max.CTCF.Rao, by="compo.Rao.membership")
p300.CTCF.Rao = merge(max.p300.Rao,max.CTCF.Rao, by="compo.Rao.membership")
H3K27me3.p300.Rao = merge(max.H3K27me3.Rao,max.p300.Rao, by="compo.Rao.membership")

comp.Rao = sapply(sort(unique(nGenesbyRRCs.Rao$geneSymbol.x))[1:6], function(i){
  subset.membership.Rao = nGenesbyRRCs.Rao[nGenesbyRRCs.Rao$geneSymbol.x<=i,"compo.Rao.membership"]
  return(cor(Expre.H3K27me3.Rao[Expre.H3K27me3.Rao$compo.Rao.membership%in%subset.membership.Rao,"mean.Pair.H3K27me3"],Expre.H3K27me3.Rao[Expre.H3K27me3.Rao$compo.Rao.membership%in%subset.membership.Rao,"Expression"]))
  
})

comp.Rao= data.frame(comp.Rao)
colnames(comp.Rao) = "COR"
comp.Rao$Variables = "Expression_H3K27me3"
comp.Rao$N = sort(unique(nGenesbyRRCs.Rao$geneSymbol.x))[1:6]

comp1.Rao = sapply(sort(unique(nGenesbyRRCs.Rao$geneSymbol.x))[1:6], function(i){
  subset.membership.Rao = nGenesbyRRCs.Rao[nGenesbyRRCs.Rao$geneSymbol.x<=i,"compo.Rao.membership"]
  return(cor(H3K27ac.H3K27me3.Rao[H3K27ac.H3K27me3.Rao$compo.Rao.membership%in%subset.membership.Rao,"mean.Pair.H3K27me3"],H3K27ac.H3K27me3.Rao[H3K27ac.H3K27me3.Rao$compo.Rao.membership%in%subset.membership.Rao,"mean.Pair.H3K27ac"]))
  
})

comp1.Rao= data.frame(comp1.Rao)
colnames(comp1.Rao) = "COR"
comp1.Rao$Variables = "H3K27ac_H3K27me3"
comp1.Rao$N = sort(unique(nGenesbyRRCs.Rao$geneSymbol.x))[1:6]

comp2.Rao = sapply(sort(unique(nGenesbyRRCs.Rao$geneSymbol.x))[1:6], function(i){
  subset.membership.Rao = nGenesbyRRCs.Rao[nGenesbyRRCs.Rao$geneSymbol.x<=i,"compo.Rao.membership"]
  return(cor(H3K4me1.H3K27me3.Rao[H3K4me1.H3K27me3.Rao$compo.Rao.membership%in%subset.membership.Rao,"mean.Pair.H3K27me3"],H3K4me1.H3K27me3.Rao[H3K4me1.H3K27me3.Rao$compo.Rao.membership%in%subset.membership.Rao,"mean.Pair.H3K4me1"]))
  
})

comp2.Rao= data.frame(comp2.Rao)
colnames(comp2.Rao) = "COR"
comp2.Rao$Variables = "H3K4me1_H3K27me3"
comp2.Rao$N = sort(unique(nGenesbyRRCs.Rao$geneSymbol.x))[1:6]

comp3.Rao = sapply(sort(unique(nGenesbyRRCs.Rao$geneSymbol.x))[1:6], function(i){
  subset.membership.Rao = nGenesbyRRCs.Rao[nGenesbyRRCs.Rao$geneSymbol.x<=i,"compo.Rao.membership"]
  return(cor(H3K4me3.H3K27me3.Rao[H3K4me3.H3K27me3.Rao$compo.Rao.membership%in%subset.membership.Rao,"mean.Pair.H3K27me3"],H3K4me3.H3K27me3.Rao[H3K4me3.H3K27me3.Rao$compo.Rao.membership%in%subset.membership.Rao,"mean.Pair.H3K4me3"]))
  
})
comp3.Rao= data.frame(comp3.Rao)
colnames(comp3.Rao) = "COR"
comp3.Rao$Variables = "H3K4me3_H3K27me3"
comp3.Rao$N = sort(unique(nGenesbyRRCs.Rao$geneSymbol.x))[1:6]

comp4.Rao = sapply(sort(unique(nGenesbyRRCs.Rao$geneSymbol.x))[1:6], function(i){
  subset.membership.Rao = nGenesbyRRCs.Rao[nGenesbyRRCs.Rao$geneSymbol.x<=i,"compo.Rao.membership"]
  return(cor(H3K27ac.CTCF.Rao[H3K27ac.CTCF.Rao$compo.Rao.membership%in%subset.membership.Rao,"mean.Pair.CTCF"],H3K27ac.CTCF.Rao[H3K27ac.CTCF.Rao$compo.Rao.membership%in%subset.membership.Rao,"mean.Pair.H3K27ac"]))
  
})
comp4.Rao= data.frame(comp4.Rao)
colnames(comp4.Rao) = "COR"
comp4.Rao$Variables = "CTCF_H3K27ac"
comp4.Rao$N = sort(unique(nGenesbyRRCs.Rao$geneSymbol.x))[1:6]

comp5.Rao = sapply(sort(unique(nGenesbyRRCs.Rao$geneSymbol.x))[1:6], function(i){
  subset.membership.Rao = nGenesbyRRCs.Rao[nGenesbyRRCs.Rao$geneSymbol.x<=i,"compo.Rao.membership"]
  return(cor(H3K27me3.CTCF.Rao[H3K27me3.CTCF.Rao$compo.Rao.membership%in%subset.membership.Rao,"mean.Pair.H3K27me3"],H3K27me3.CTCF.Rao[H3K27me3.CTCF.Rao$compo.Rao.membership%in%subset.membership.Rao,"mean.Pair.CTCF"]))
  
})
comp5.Rao= data.frame(comp5.Rao)
colnames(comp5.Rao) = "COR"
comp5.Rao$Variables = "CTCF_H3K27me3"
comp5.Rao$N = sort(unique(nGenesbyRRCs.Rao$geneSymbol.x))[1:6]

comp6.Rao = sapply(sort(unique(nGenesbyRRCs.Rao$geneSymbol.x))[1:6], function(i){
  subset.membership.Rao = nGenesbyRRCs.Rao[nGenesbyRRCs.Rao$geneSymbol.x<=i,"compo.Rao.membership"]
  return(cor(H3K4me1.CTCF.Rao[H3K4me1.CTCF.Rao$compo.Rao.membership%in%subset.membership.Rao,"mean.Pair.H3K4me1"],H3K4me1.CTCF.Rao[H3K4me1.CTCF.Rao$compo.Rao.membership%in%subset.membership.Rao,"mean.Pair.CTCF"]))
  
})
comp6.Rao= data.frame(comp6.Rao)
colnames(comp6.Rao) = "COR"
comp6.Rao$Variables = "CTCF_H3K4me1"
comp6.Rao$N = sort(unique(nGenesbyRRCs.Rao$geneSymbol.x))[1:6]

comp7.Rao = sapply(sort(unique(nGenesbyRRCs.Rao$geneSymbol.x))[1:6], function(i){
  subset.membership.Rao = nGenesbyRRCs.Rao[nGenesbyRRCs.Rao$geneSymbol.x<=i,"compo.Rao.membership"]
  return(cor(H3K4me3.CTCF.Rao[H3K4me3.CTCF.Rao$compo.Rao.membership%in%subset.membership.Rao,"mean.Pair.H3K4me3"],H3K4me3.CTCF.Rao[H3K4me3.CTCF.Rao$compo.Rao.membership%in%subset.membership.Rao,"mean.Pair.CTCF"]))
  
})
comp7.Rao= data.frame(comp7.Rao)
colnames(comp7.Rao) = "COR"
comp7.Rao$Variables = "CTCF_H3K4me3"
comp7.Rao$N = sort(unique(nGenesbyRRCs.Rao$geneSymbol.x))[1:6]

comp8.Rao = sapply(sort(unique(nGenesbyRRCs.Rao$geneSymbol.x))[1:6], function(i){
  subset.membership.Rao = nGenesbyRRCs.Rao[nGenesbyRRCs.Rao$geneSymbol.x<=i,"compo.Rao.membership"]
  return(cor(p300.CTCF.Rao[p300.CTCF.Rao$compo.Rao.membership%in%subset.membership.Rao,"mean.Pair.p300"],p300.CTCF.Rao[p300.CTCF.Rao$compo.Rao.membership%in%subset.membership.Rao,"mean.Pair.CTCF"]))
  
})

comp8.Rao= data.frame(comp8.Rao)
colnames(comp8.Rao) = "COR"
comp8.Rao$Variables = "CTCF_p300"
comp8.Rao$N = sort(unique(nGenesbyRRCs.Rao$geneSymbol.x))[1:6]

comp9.Rao = sapply(sort(unique(nGenesbyRRCs.Rao$geneSymbol.x))[1:6], function(i){
  subset.membership.Rao = nGenesbyRRCs.Rao[nGenesbyRRCs.Rao$geneSymbol.x<=i,"compo.Rao.membership"]
  return(cor(H3K27me3.p300.Rao[H3K27me3.p300.Rao$compo.Rao.membership%in%subset.membership.Rao,"mean.Pair.p300"],H3K27me3.p300.Rao[H3K27me3.p300.Rao$compo.Rao.membership%in%subset.membership.Rao,"mean.Pair.H3K27me3"]))
  
})

comp9.Rao= data.frame(comp9.Rao)
colnames(comp9.Rao) = "COR"
comp9.Rao$Variables = "H3K27me3_p300"
comp9.Rao$N = sort(unique(nGenesbyRRCs.Rao$geneSymbol.x))[1:6]


comp.tt.Rao = rbind(comp.Rao, comp1.Rao,comp2.Rao,comp3.Rao,comp4.Rao,comp9.Rao,comp5.Rao,comp6.Rao,comp7.Rao, comp8.Rao)
ggl.Rao = ggline(comp.tt.Rao, x="N",y="COR",plot_type="l",color="Variables")
ggl.Rao = ggpar(ggl.Rao, legend.title = "Association", xlab="#promoters in RRC", ylab="Correlation")


max.Expression.DNAse = aggregate(Expression~compo.DNAse.membership, AL2Genes.DNAse, max)
max.H3K27ac.DNAse = aggregate(mean.Pair.H3K27ac~compo.DNAse.membership, AL2Genes.DNAse, max)
max.CTCF.DNAse = aggregate(mean.Pair.CTCF~compo.DNAse.membership, AL2Genes.DNAse, max)
max.DNAse.DNAse = aggregate(mean.Pair.DNAse~compo.DNAse.membership, AL2Genes.DNAse, max)
max.H3K4me1.DNAse = aggregate(mean.Pair.H3K4me1~compo.DNAse.membership, AL2Genes.DNAse, max)
max.H3K4me3.DNAse = aggregate(mean.Pair.H3K4me3~compo.DNAse.membership, AL2Genes.DNAse, max)
max.H3K27me3.DNAse = aggregate(mean.Pair.H3K27me3~compo.DNAse.membership, AL2Genes.DNAse, max)
max.p300.DNAse = aggregate(mean.Pair.p300~compo.DNAse.membership, AL2Genes.DNAse, max)

Expr.H3K27me3.DNAse = merge(max.Expression.DNAse, max.H3K27me3.DNAse, by="compo.DNAse.membership")
H3K27ac.H3K27me3.DNAse = merge(max.H3K27ac.DNAse, max.H3K27me3.DNAse, by="compo.DNAse.membership")
H3K4me1.H3K27me3.DNAse = merge(max.H3K4me1.DNAse, max.H3K27me3.DNAse, by="compo.DNAse.membership")
H3K4me3.H3K27me3.DNAse = merge(max.H3K4me3.DNAse, max.H3K27me3.DNAse, by="compo.DNAse.membership")
H3K27ac.CTCF.DNAse = merge(max.H3K27ac.DNAse, max.CTCF.DNAse,by="compo.DNAse.membership")
H3K27me3.CTCF.DNAse = merge(max.H3K27me3.DNAse, max.CTCF.DNAse,by="compo.DNAse.membership")
H3K4me1.CTCF.DNAse = merge(max.H3K4me1.DNAse,max.CTCF.DNAse, by="compo.DNAse.membership")
H3K4me3.CTCF.DNAse = merge(max.H3K4me3.DNAse,max.CTCF.DNAse, by="compo.DNAse.membership")
p300.CTCF.DNAse = merge(max.p300.DNAse,max.CTCF.DNAse, by="compo.DNAse.membership")
H3K27me3.p300.DNAse = merge(max.H3K27me3.DNAse,max.p300.DNAse, by="compo.DNAse.membership")


comp.DNAse = sapply(sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)), function(i){
  subset.membership.DNAse = nGenesbyRRCs.DNAse[nGenesbyRRCs.DNAse$geneSymbol.y<=i,"compo.DNAse.membership"]
  return(cor(Expr.H3K27me3.DNAse[Expr.H3K27me3.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,"mean.Pair.H3K27me3"],Expr.H3K27me3.DNAse[Expr.H3K27me3.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,"Expression"]))
  
})

comp.DNAse= data.frame(comp.DNAse)
colnames(comp.DNAse) = "COR"
comp.DNAse$Variables = "Expression_H3K27me3"
comp.DNAse$N = sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y))


comp1.DNAse = sapply(sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)), function(i){
  subset.membership.DNAse = nGenesbyRRCs.DNAse[nGenesbyRRCs.DNAse$geneSymbol.y<=i,"compo.DNAse.membership"]
  return(cor(H3K27ac.H3K27me3.DNAse[H3K27ac.H3K27me3.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,"mean.Pair.H3K27me3"],H3K27ac.H3K27me3.DNAse[H3K27ac.H3K27me3.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,"mean.Pair.H3K27ac"]))
  
})

comp1.DNAse= data.frame(comp1.DNAse)
colnames(comp1.DNAse) = "COR"
comp1.DNAse$Variables = "H3K27ac_H3K27me3"
comp1.DNAse$N = sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y))

comp2.DNAse = sapply(sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)), function(i){
  subset.membership.DNAse = nGenesbyRRCs.DNAse[nGenesbyRRCs.DNAse$geneSymbol.y<=i,"compo.DNAse.membership"]
  return(cor(H3K4me1.H3K27me3.DNAse[H3K4me1.H3K27me3.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,"mean.Pair.H3K27me3"],H3K4me1.H3K27me3.DNAse[H3K4me1.H3K27me3.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,"mean.Pair.H3K4me1"]))
  
})

comp2.DNAse= data.frame(comp2.DNAse)
colnames(comp2.DNAse) = "COR"
comp2.DNAse$Variables = "H3K4me1_H3K27me3"
comp2.DNAse$N = sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y))

comp3.DNAse = sapply(sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)), function(i){
  subset.membership.DNAse = nGenesbyRRCs.DNAse[nGenesbyRRCs.DNAse$geneSymbol.y<=i,"compo.DNAse.membership"]
  return(cor(H3K4me3.H3K27me3.DNAse[H3K4me3.H3K27me3.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,"mean.Pair.H3K27me3"],H3K4me3.H3K27me3.DNAse[H3K4me3.H3K27me3.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,"mean.Pair.H3K4me3"]))
  
})
comp3.DNAse= data.frame(comp3.DNAse)
colnames(comp3.DNAse) = "COR"
comp3.DNAse$Variables = "H3K4me3_H3K27me3"
comp3.DNAse$N = sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y))

comp4.DNAse = sapply(sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)), function(i){
  subset.membership.DNAse = nGenesbyRRCs.DNAse[nGenesbyRRCs.DNAse$geneSymbol.y<=i,"compo.DNAse.membership"]
  return(cor(H3K27ac.CTCF.DNAse[H3K27ac.CTCF.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,"mean.Pair.CTCF"],H3K27ac.CTCF.DNAse[H3K27ac.CTCF.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,"mean.Pair.H3K27ac"]))
  
})
comp4.DNAse= data.frame(comp4.DNAse)
colnames(comp4.DNAse) = "COR"
comp4.DNAse$Variables = "CTCF_H3K27ac"
comp4.DNAse$N = sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y))

comp5.DNAse = sapply(sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)), function(i){
  subset.membership.DNAse = nGenesbyRRCs.DNAse[nGenesbyRRCs.DNAse$geneSymbol.y<=i,"compo.DNAse.membership"]
  return(cor(H3K27me3.CTCF.DNAse[H3K27me3.CTCF.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,"mean.Pair.H3K27me3"],H3K27me3.CTCF.DNAse[H3K27me3.CTCF.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,"mean.Pair.CTCF"]))
  
})
comp5.DNAse= data.frame(comp5.DNAse)
colnames(comp5.DNAse) = "COR"
comp5.DNAse$Variables = "CTCF_H3K27me3"
comp5.DNAse$N = sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y))

comp6.DNAse = sapply(sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)), function(i){
  subset.membership.DNAse = nGenesbyRRCs.DNAse[nGenesbyRRCs.DNAse$geneSymbol.y<=i,"compo.DNAse.membership"]
  return(cor(H3K4me1.CTCF.DNAse[H3K4me1.CTCF.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,"mean.Pair.H3K4me1"],H3K4me1.CTCF.DNAse[H3K4me1.CTCF.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,"mean.Pair.CTCF"]))
  
})
comp6.DNAse= data.frame(comp6.DNAse)
colnames(comp6.DNAse) = "COR"
comp6.DNAse$Variables = "CTCF_H3K4me1"
comp6.DNAse$N = sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y))

comp7.DNAse = sapply(sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)), function(i){
  subset.membership.DNAse = nGenesbyRRCs.DNAse[nGenesbyRRCs.DNAse$geneSymbol.y<=i,"compo.DNAse.membership"]
  return(cor(H3K4me3.CTCF.DNAse[H3K4me3.CTCF.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,"mean.Pair.H3K4me3"],H3K4me3.CTCF.DNAse[H3K4me3.CTCF.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,"mean.Pair.CTCF"]))
  
})
comp7.DNAse= data.frame(comp7.DNAse)
colnames(comp7.DNAse) = "COR"
comp7.DNAse$Variables = "CTCF_H3K4me3"
comp7.DNAse$N = sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y))

comp8.DNAse = sapply(sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)), function(i){
  subset.membership.DNAse = nGenesbyRRCs.DNAse[nGenesbyRRCs.DNAse$geneSymbol.y<=i,"compo.DNAse.membership"]
  return(cor(p300.CTCF.DNAse[p300.CTCF.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,"mean.Pair.p300"],p300.CTCF.DNAse[p300.CTCF.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,"mean.Pair.CTCF"]))
  
})

comp8.DNAse= data.frame(comp8.DNAse)
colnames(comp8.DNAse) = "COR"
comp8.DNAse$Variables = "CTCF_p300"
comp8.DNAse$N = sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y))

comp9.DNAse = sapply(sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y)), function(i){
  subset.membership.DNAse = nGenesbyRRCs.DNAse[nGenesbyRRCs.DNAse$geneSymbol.y<=i,"compo.DNAse.membership"]
  return(cor(H3K27me3.p300.DNAse[H3K27me3.p300.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,"mean.Pair.p300"],H3K27me3.p300.DNAse[H3K27me3.p300.DNAse$compo.DNAse.membership%in%subset.membership.DNAse,"mean.Pair.H3K27me3"]))
  
})

comp9.DNAse= data.frame(comp9.DNAse)
colnames(comp9.DNAse) = "COR"
comp9.DNAse$Variables = "H3K27me3_p300"
comp9.DNAse$N = sort(unique(nGenesbyRRCs.DNAse$geneSymbol.y))

comp.tt.DNAse = rbind(comp.DNAse, comp1.DNAse,comp2.DNAse,comp3.DNAse,comp9.DNAse,comp4.DNAse,comp5.DNAse,comp6.DNAse,comp7.DNAse,comp8.DNAse)
ggl.DNAse = ggline(comp.tt.DNAse, x="N",y="COR",plot_type="l",color="Variables")
ggl.DNAse = ggpar(ggl.DNAse, legend.title = "Association", xlab="#promoters in RRC", ylab="Correlation")


ggarrange(ggl,ggl.Rao, ggl.DNAse, labels=c("A","B","C"), nrow=3, ncol=1)
