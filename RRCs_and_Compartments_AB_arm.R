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
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
source("HiCrn.R")

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes.hg19 <- genes(txdb)
symbol <- select(org.Hs.eg.db,keys = genes.hg19$gene_id, columns = c("SYMBOL"), keytype = "ENTREZID")
genes.hg19$geneSymbol <- symbol$SYMBOL

###############################################################################################
##################################RRCs & Compartiments#########################################
###############################################################################################
ta = define.active.compartments.arms("PC_oe/allchrarms_PC_oe.500Kb.txt", resolution=500000, genes=genes.hg19)
sum(abs(ta$COR.PC1.GDENSITE)>abs(ta$COR.PC2.GDENSITE) & abs(ta$COR.PC1.GDENSITE)>abs(ta$COR.PC3.GDENSITE))
[1] 22

ta.100Kb = define.active.compartments.arms("PC_oe/allchrarms_PC_oe.100Kb.txt", resolution=100000, genes=genes.hg19)
sum(abs(ta.100Kb$COR.PC1.GDENSITE)>abs(ta.100Kb$COR.PC2.GDENSITE) & abs(ta.100Kb$COR.PC1.GDENSITE)>abs(ta.100Kb$COR.PC3.GDENSITE))
[1] 22

talog = define.active.compartments.arms("PC_logoe/allchrarms_PC_logoe.500Kb.txt", resolution=500000, genes=genes.hg19)
sum(abs(talog$COR.PC1.GDENSITE)>abs(talog$COR.PC2.GDENSITE) & abs(talog$COR.PC1.GDENSITE)>abs(talog$COR.PC3.GDENSITE))
[1] 21

talogc = define.active.compartments.arms("PC_logoe/allchrarms_PC_logoe_clipped.500Kb.txt", resolution=500000, genes=genes.hg19)
sum(abs(talogc$COR.PC1.GDENSITE)>abs(talogc$COR.PC2.GDENSITE) & abs(talogc$COR.PC1.GDENSITE)>abs(talogc$COR.PC3.GDENSITE))
[1] 24

# Calcul de corrélation avec le contenu GC

GRanges.500Kb.PCarms = .make.GRanges.compartments.arms("PC_oe/allchrarms_PC_oe.500Kb.txt",resolution=500000)
GRanges.500Kb.PCarms$gc = compute_GCcontent(GRanges.500Kb.PCarms,genome=getBSgenome("hg19"))
tagc = define.active.compartments.arms.GC(GRanges.500Kb.PCarms)
sum(abs(tagc$COR.PC1.GC)>abs(tagc$COR.PC2.GC) & abs(tagc$COR.PC1.GC)>abs(tagc$COR.PC3.GC))
[1] 26

GRanges.100Kb.PCarms = .make.GRanges.compartments.arms("PC_oe/allchrarms_PC_oe.100Kb.txt",resolution=100000)
GRanges.100Kb.PCarms$gc = compute_GCcontent(GRanges.100Kb.PCarms,genome=getBSgenome("hg19"))
tagc.100Kb = define.active.compartments.arms.GC(GRanges.100Kb.PCarms)
sum(abs(tagc.100Kb$COR.PC1.GC)>abs(tagc.100Kb$COR.PC2.GC) & abs(tagc.100Kb$COR.PC1.GC)>abs(tagc.100Kb$COR.PC3.GC))
[1] 26

GRanges.500Kb.PCarms.log = .make.GRanges.compartments.arms("PC_logoe/allchrarms_PC_logoe.500Kb.txt",resolution=500000)
GRanges.500Kb.PCarms.log$gc = GRanges.500Kb.PCarms$gc 
taloggc = define.active.compartments.arms.GC(GRanges.500Kb.PCarms.log)
sum(abs(taloggc$COR.PC1.GC)>abs(taloggc$COR.PC2.GC) & abs(taloggc$COR.PC1.GC)>abs(taloggc$COR.PC3.GC))
[1] 26

GRanges.100Kb.PCarms.log = .make.GRanges.compartments.arms("PC_logoe/allchrarms_PC_logoe.100Kb.txt",resolution=100000)
GRanges.100Kb.PCarms.log$gc = GRanges.100Kb.PCarms$gc
taloggc.100Kb = define.active.compartments.arms.GC(GRanges.100Kb.PCarms.log)
sum(abs(taloggc.100Kb$COR.PC1.GC)>abs(taloggc.100Kb$COR.PC2.GC) & abs(taloggc.100Kb$COR.PC1.GC)>abs(taloggc.100Kb$COR.PC3.GC))
[1] 21

GRanges.500Kb.PCarms.logc = .make.GRanges.compartments.arms("PC_logoe/allchrarms_PC_logoe_clipped.500Kb.txt",resolution=500000)
GRanges.500Kb.PCarms.logc$gc = GRanges.500Kb.PCarms$gc 
taloggcc = define.active.compartments.arms.GC(GRanges.500Kb.PCarms.logc)
sum(abs(taloggcc$COR.PC1.GC)>abs(taloggcc$COR.PC2.GC) & abs(taloggcc$COR.PC1.GC)>abs(taloggcc$COR.PC3.GC))
[1] 30

GRanges.100Kb.PCarms.logc = .make.GRanges.compartments.arms("PC_logoe/allchrarms_PC_logoe_clipped.100Kb.txt",resolution=100000)
GRanges.100Kb.PCarms.logc$gc = GRanges.100Kb.PCarms$gc
taloggcc.100Kb = define.active.compartments.arms.GC(GRanges.100Kb.PCarms.logc)
sum(abs(taloggcc.100Kb$COR.PC1.GC)>abs(taloggcc.100Kb$COR.PC2.GC) & abs(taloggcc.100Kb$COR.PC1.GC)>abs(taloggcc.100Kb$COR.PC3.GC))
[1] 21


t$ngenes = countOverlaps(t, genes.hg19)

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
