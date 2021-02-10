library(igraph)
library(GenomicRanges)
library(rtracklayer)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(plyr)
library(GenomicFeatures)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(lme4)
library(gee)
library(ggplot2)
library(mgcv)
library(betareg)
library(broom)
library(knitr)
library(kableExtra)
library(sandwich)
library(mice)

set.seed(1258)
setwd("/home/loic/Documents/HiC/data/4Script/NEU")

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes.hg19 <- genes(txdb)
symbol <- select(org.Hs.eg.db,keys = genes.hg19$gene_id, columns = c("SYMBOL"), keytype = "ENTREZID")
genes.hg19$geneSymbol <- symbol$SYMBOL

regul.promoters.ABC = process.ABC("ABC/EnhancerPredictions.txt")
# all.putative.ABC = read.table("ABC/EnhancerPredictionsAllPutative.txt.gz", header = T)
# all.putative.ABC = all.putative.ABC[all.putative.ABC$chr!="chrX",]
# 
# all.putative.ABC.candidates = all.putative.ABC[!all.putative.ABC$name%in%regul.promoters.ABC$name&all.putative.ABC$class!="promoter",]

hic.loops = read.table("allloops_HiCCUPS/enriched_pixels_10000.bedpe", header=T)

peaks.NEU = import("Epigenetic/NEU_DNAse.macs2_peaks.narrowPeak.sorted.candidateRegions.bed", format="bed")
peaks.NEU = peaks.NEU[seqnames(peaks.NEU)%in%paste0("chr", 1:22)]

#Epigenetic data
#DNAse: accessbility peaks 
DNAseActivity.NEU = import("Epigenetic/NEU_DNAse.macs2_peaks.narrowPeak.sorted", format="narrowpeak")
DNAseActivity.NEU$signalValue = log2(DNAseActivity.NEU$signalValue+1)

#H3K27ac: peaks
H3K27acActivity.NEU = import("Epigenetic/NEU_H3k27ac_sorted_macs2_peaks.narrowPeak.sorted", format = "narrowpeak")
H3K27acActivity.NEU$signalValue = log2(H3K27acActivity.NEU$signalValue+1)

#H3K4me1: peaks
H3K4me1Activity.NEU = import("Epigenetic/NEU_H3K4me1.peaks.macs2_peaks.narrowPeak", format = "narrowpeak")
H3K4me1Activity.NEU$signalValue = log2(H3K4me1Activity.NEU$signalValue +1)

H3K4me3.NEU.REP1 = import("Epigenetic/NEU_H3k4me3_REP1.bed.gz", format = "narrowpeak")
H3K4me3.NEU.REP2 = import("Epigenetic/NEU_H3k4me3_REP2.bed.gz", format = "narrowpeak")
H3K4me3.NEU.REP3 = import("Epigenetic/NEU_H3k4me3_REP3.bed.gz", format = "narrowpeak")

grl.H3K4me3.NEU = GRangesList(H3K4me3.NEU.REP1, H3K4me3.NEU.REP2, H3K4me3.NEU.REP3)

#H3K4me3: peaks
H3K4me3Activity.NEU = unique(do.call("c", as(grl.H3K4me3.NEU, "GRangesList")))
H3K4me3Activity.NEU$signalValue = log2(H3K4me3Activity.NEU$signalValue+1)

#H3K27me3: peaks
H3K27me3.NEU.REP1 = import("Epigenetic/NEU_H3K27me3_Rep1.bed.gz", format = "narrowpeak")
H3K27me3.NEU.REP2 = import("Epigenetic/NEU_H3K27me3_Rep2.bed.gz", format = "narrowpeak")
H3K27me3.NEU.REP3 = import("Epigenetic/NEU_H3K27me3_Rep3.bed.gz", format = "narrowpeak")

grl.H3K27me3.NEU = GRangesList(H3K27me3.NEU.REP1, H3K27me3.NEU.REP2,H3K27me3.NEU.REP3)

H3K27me3Activity.NEU = unique(do.call("c", as(grl.H3K27me3.NEU, "GRangesList")))
H3K27me3Activity.NEU$signalValue = log2(H3K27me3Activity.NEU$signalValue+1)

#CTCF
CTCF.NEU.REP1 = import("Epigenetic/NEU_REP1_CTCF.bed.gz", format = "narrowpeak")
CTCF.NEU.REP2 =   import("Epigenetic/NEU_REP2_CTCF.bed.gz", format = "narrowpeak")
grL.CTCF.NEU = GRangesList(CTCF.NEU.REP1, CTCF.NEU.REP2)

CTCFActivity.NEU = unique(do.call("c", as(grL.CTCF.NEU, "GRangesList")))
CTCFActivity.NEU$signalValue = log2(CTCFActivity.NEU$signalValue+1)

#p300
p300.NEU.REP1 = import("Epigenetic/NEU_p300_REP1.bed.gz", format = "narrowpeak")
p300.NEU.REP2 = import("Epigenetic/NEU_p300_REP1.bed.gz", format = "narrowpeak")
p300.NEU.REP3 = import("Epigenetic/NEU_p300_REP1.bed.gz", format = "narrowpeak")
grL.p300.NEU = GRangesList(p300.NEU.REP1, p300.NEU.REP2, p300.NEU.REP3)

p300Activity.NEU = unique(do.call("c", as(grL.p300.NEU, "GRangesList")))
p300Activity.NEU$signalValue = log2(p300Activity.NEU$signalValue+1)

#Gene Expression
expression.NEU = read.table("Expression/Expression_neurons.txt",header=T)
expression.NEU$GENEID = rownames(expression.NEU)
expression.NEU[,-ncol(expression.NEU)] = log(expression.NEU[,-ncol(expression.NEU)]+1)
correspondance = ensembldb::select(EnsDb.Hsapiens.v86, keys=rownames(expression.NEU), keytype = "GENEID", columns=c("SYMBOL", "GENEID"))

expression.NEU = merge(expression.NEU, correspondance, by="GENEID")
expression.NEU$geneSymbol =  expression.NEU$SYMBOL
expression.NEU$SYMBOL = NULL

expression.NEU$median.Expr = apply(expression.NEU[,-c(1,ncol(expression.NEU))],1, median)
expression.genes = makeGRangesFromDataFrame(merge(data.frame(genes.hg19), expression.NEU[,c("geneSymbol", "median.Expr")], by="geneSymbol"), keep.extra.columns = T)

#3D features
#FIREs
FIREs = read.table("3D_features/FIREs_NEU.txt", header=T)
superFIREs = read.table("3D_features//NEU_SuperFIREs.txt", header=T)
FIREs.signi = FIREs[FIREs$NEU_indicator==1,]

GRanges.FIREs = GRanges(seqnames = FIREs$chr, ranges=IRanges(start=FIREs$start, end=FIREs$end, names = paste0("FIRE",1:nrow(FIREs))),ScoreFire=FIREs$NEU_neg_ln_pval)
GRanges.FIREs.signi = GRanges(seqnames = FIREs.signi$chr, ranges=IRanges(start=FIREs.signi$start, end=FIREs.signi$end, names = paste0("FIRE",1:nrow(FIREs.signi))), x =FIREs.signi$NEU_neg_ln_pval)
GRanges.superFIREs = GRanges(seqnames = superFIREs$chr, ranges=IRanges(start=superFIREs$start, end=superFIREs$end, names = paste0("superFIRE",1:nrow(superFIREs))))

#DI
DI = read.table("3D_features/all_chrs_dense_annotated.matrix.DI", header=F)
colnames(DI) = c("chr", "start", "end", "DI")
DI$chr = paste0("chr", DI$chr)

DI.WX = DI[DI$chr!="chr23",]
GRanges.DI = GRanges(seqnames = DI.WX$chr, ranges=IRanges(start=DI.WX$start, end=DI.WX$end), DI= DI.WX$DI)

#INS
INS = read.table("3D_features/allchrs.insulation", header = T)

INS$insulationScore = as.numeric(as.character(INS$insulationScore))
INS$start = as.numeric(as.character(INS$start))
INS$end = as.numeric(as.character(INS$end))

INS = na.omit(INS)
INS$normalizedScore = (INS$insulationScore - min(INS$insulationScore) )/ (max(INS$insulationScore) - min(INS$insulationScore))
INS$chr=gsub(".*[|]([^.]+)[:].*", "\\1", INS$header)

INS.WX = INS[INS$chr!="chrX",]

GRanges.INS = GRanges(seqnames =  INS.WX$chr, ranges = IRanges(start = INS.WX$start, end = INS.WX$end), normalizedScore = INS.WX$normalizedScore)


#Start of Script
GRanges.loops.bin1 = GRanges(seqnames = hic.loops$chr1, ranges=IRanges(start=hic.loops$x1, end=hic.loops$x2))
GRanges.loops.bin2 = GRanges(seqnames = hic.loops$chr2, ranges=IRanges(start=hic.loops$y1, end=hic.loops$y2))

Pairs.bins = Pairs(GRanges.loops.bin1, GRanges.loops.bin2,fdrBL = hic.loops$fdrBL ,fdrDonut = hic.loops$fdrDonut,fdrV=hic.loops$fdrV,fdrH = hic.loops$fdrH)
positive.contact = Pairs.bins[mcols(Pairs.bins)$fdrBL <= 0.15&mcols(Pairs.bins)$fdrDonut <=0.15&mcols(Pairs.bins)$fdrV <= 0.15&mcols(Pairs.bins)$fdrH <=0.15]


ABC.Pairs = Pairs.Enh.Prom.ABC(regul.promoters.ABC)
Rao.Pairs = process.Rao.DNAse(positive.contact, method = "Rao")
DNAse.Pairs = process.Rao.DNAse(positive.contact,peaks.NEU, method = "DNAse")


list.Prom.Regul = list("ABC" = list("Prom"= first(ABC.Pairs),"Regul"= second(ABC.Pairs)), "Rao" = list("Prom"= first(Rao.Pairs), "Regul"= second(Rao.Pairs)), 
                       "DNAse" = list("Prom"= first(DNAse.Pairs), "Regul"= second(DNAse.Pairs)))

list.Prom.Regul.annotated = lapply(1:length(list.Prom.Regul), function(x) lapply(1:length(list.Prom.Regul[[x]]), function(y) {
  list.Prom.Regul[[x]][[y]] = annotate.3D.Features(GRanges.DI, list.Prom.Regul[[x]][[y]], kind = "DI",aggregateFunction ="mean")
  list.Prom.Regul[[x]][[y]] = annotate.3D.Features(GRanges.FIREs, list.Prom.Regul[[x]][[y]], kind = "FIRE",aggregateFunction ="mean")
  list.Prom.Regul[[x]][[y]] = annotate.3D.Features(GRanges.INS, list.Prom.Regul[[x]][[y]], kind = "INS",aggregateFunction ="mean")
  
  list.Prom.Regul[[x]][[y]] = annotate.Activity(CTCFActivity.NEU, list.Prom.Regul[[x]][[y]], kind = "CTCF")
  list.Prom.Regul[[x]][[y]] = annotate.Activity(H3K27acActivity.NEU, list.Prom.Regul[[x]][[y]], kind = "H3K27ac")
  list.Prom.Regul[[x]][[y]] = annotate.Activity(DNAseActivity.NEU, list.Prom.Regul[[x]][[y]], kind = "DNAse")
  list.Prom.Regul[[x]][[y]] = annotate.Activity(H3K4me1Activity.NEU, list.Prom.Regul[[x]][[y]], kind = "H3K4me1")
  list.Prom.Regul[[x]][[y]] = annotate.Activity(H3K4me3Activity.NEU, list.Prom.Regul[[x]][[y]], kind = "H3K4me3")
  list.Prom.Regul[[x]][[y]] = annotate.Activity(p300Activity.NEU, list.Prom.Regul[[x]][[y]], kind = "p300")
  list.Prom.Regul[[x]][[y]] = annotate.Activity(H3K27me3Activity.NEU, list.Prom.Regul[[x]][[y]], kind = "H3K27me3")
  
}

))

ABC.Pairs = Pairs(list.Prom.Regul.annotated[[1]][[1]],list.Prom.Regul.annotated[[1]][[2]])
Rao.Pairs = Pairs(list.Prom.Regul.annotated[[2]][[1]],list.Prom.Regul.annotated[[2]][[2]])
DNAse.Pairs = Pairs(list.Prom.Regul.annotated[[3]][[1]],list.Prom.Regul.annotated[[3]][[2]])

regul.promoters.ABC = annotate.dataframe(ABC.Pairs, method = "ABC", df.ABC = regul.promoters.ABC)
df.Rao = annotate.dataframe(Rao.Pairs, method="Rao")
df.DNase = annotate.dataframe(DNAse.Pairs, method="DNAse")

#Unique elements by annotation methods 

unique.Promoters.ABC = unique(second(ABC.Pairs))
unique.Regul.ABC = unique(first(ABC.Pairs))
unique.Promoters.Rao = unique(first(Rao.Pairs))
unique.Regul.Rao = unique(second(Rao.Pairs))
unique.Promoters.DNAse = unique(first(DNAse.Pairs))
unique.Regul.DNAse =  unique(second(DNAse.Pairs))

#Individual correlations 
#ABC
individual.correlation(unique.Regul.ABC, annotation = "ABC", method="pearson",check.3D.circularity = TRUE)
individual.correlation(unique.Promoters.ABC, annotation = "ABC", method="pearson",check.3D.circularity = TRUE)
#Rao
individual.correlation(unique.Regul.Rao, annotation = "Rao", method="pearson")
individual.correlation(unique.Promoters.Rao, annotation = "Rao", method="pearson")
#DNAse
individual.correlation(unique.Regul.DNAse, annotation = "DNAse", method="pearson")
individual.correlation(unique.Promoters.DNAse, annotation = "DNAse", method="pearson")

#Enrichment in FIREs 
GRanges.FIREs.signi = GRanges.FIREs[GRanges.FIREs$ScoreFire>=3]
GRanges.FIREs.nsigni = GRanges.FIREs[GRanges.FIREs$ScoreFire<3]


#Enrichment of regulatory elements in CRNs 
Candidates.Signi.H3K27ac = import("Encode_crns/Encode_H3K27ac_signi.txt", format="bed")
Candidates.NSigni.H3K27ac = import("Encode_crns/Encode_H3K27ac_nsigni.txt", format="bed")

Candidates.Signi.H3K4me3 = import("Encode_crns/Encode_H3K4me3_signi.txt", format="bed")
Candidates.NSigni.H3K4me3 = import("Encode_crns/Encode_H3K4me3_nsigni.txt", format="bed")

Candidates.Signi.CTCF = import("Encode_crns/Encode_CTCF_signi.txt", format="bed")
Candidates.NSigni.CTCF = import("Encode_crns/Encode_CTCF_nsigni.txt", format="bed")

Candidates.Signi.DNAse = import("Encode_crns/Encode_DNAse_signi.txt", format="bed")
Candidates.NSigni.DNAse = import("Encode_crns/Encode_DNAse_nsigni.txt", format="bed")

Candidates.Signi.ALL = import("Encode_crns/Encode_ALL_signi.txt", format="bed")
Candidates.NSigni.ALL = import("Encode_crns/Encode_ALL_nsigni.txt", format="bed")

#Comparable set building for each annotation method
regul.candidates.ABC = make.comparable.set(unique.Regul.ABC,method="ABC", element="regulatory", DNAse = peaks.NEU)
promoters.candidates.ABC = make.comparable.set(unique.Promoters.ABC,method="ABC", element="promoter")
regul.candidates.Rao = make.comparable.set(unique.Regul.Rao,method="Rao", element="regulatory", contact =positive.contact )
promoters.candidates.Rao = make.comparable.set(unique.Promoters.Rao,method="Rao", element="promoter", contact = positive.contact)
regul.candidates.DNAse = make.comparable.set(unique.Promoters.DNAse,method="DNAse", element="regulatory",contact = positive.contact, DNAse = peaks.NEU)
promoters.candidates.DNAse = make.comparable.set(unique.Promoters.DNAse,method="DNAse", element="promoter", contact = positive.contact, DNAse = peaks.NEU)

OR.FIREs.prom.ABC = enrichments.analysis(unique.Promoters.ABC,promoters.candidates.ABC,GRanges.FIREs.signi,GRanges.FIREs.nsigni)
OR.FIREs.prom.Rao = enrichments.analysis(unique.Promoters.Rao,promoters.candidates.Rao,GRanges.FIREs.signi,GRanges.FIREs.nsigni)
OR.FIREs.prom.DNAse = enrichments.analysis(unique.Promoters.DNAse,promoters.candidates.DNAse,GRanges.FIREs.signi,GRanges.FIREs.nsigni)

OR.FIREs.regul.ABC = enrichments.analysis(unique.Regul.ABC,regul.candidates.ABC,GRanges.FIREs.signi,GRanges.FIREs.nsigni)
OR.FIREs.regul.Rao = enrichments.analysis(unique.Regul.Rao,regul.candidates.Rao,GRanges.FIREs.signi,GRanges.FIREs.nsigni)
OR.FIREs.regul.DNAse = enrichments.analysis(unique.Regul.DNAse,regul.candidates.DNAse,GRanges.FIREs.signi,GRanges.FIREs.nsigni)

FIREs.prom = plot.OR(OR.FIREs.prom.ABC,OR.FIREs.prom.Rao,OR.FIREs.prom.DNAse)
FIREs.regul = plot.OR(OR.FIREs.regul.ABC,OR.FIREs.regul.Rao,OR.FIREs.regul.DNAse)

ggpubr::ggarrange(FIREs.prom, FIREs.regul, labels = c("A","B"), ncol = 1, nrow = 2)

OR.H3K4me3.regul.ABC = enrichments.analysis(unique.Regul.ABC,regul.candidates.ABC,Candidates.Signi.H3K4me3,Candidates.NSigni.H3K4me3)
OR.CTCF.regul.ABC = enrichments.analysis(unique.Regul.ABC,regul.candidates.ABC,Candidates.Signi.CTCF,Candidates.NSigni.CTCF)
OR.ALL.regul.ABC = enrichments.analysis(unique.Regul.ABC,regul.candidates.ABC,Candidates.Signi.ALL,Candidates.NSigni.ALL)


df.OR.ALL.ABC = data.frame(OR.H3K4me3.regul.ABC$estimate[[1]], OR.H3K4me3.regul.ABC$conf.int[1],OR.H3K4me3.regul.ABC$conf.int[2])
df.OR.CTCF.ABC = data.frame(OR.CTCF.regul.ABC$estimate[[1]],OR.CTCF.regul.ABC$conf.int[1], OR.CTCF.regul.ABC$conf.int[2])
df.OR.H3K4me3.ABC = data.frame(OR.H3K4me3.regul.ABC$estimate[[1]], OR.H3K4me3.regul.ABC$conf.int[1],OR.H3K4me3.regul.ABC$conf.int[2])

colnames(df.OR.ALL.ABC) = c("OR", "CI_lower", "CI_upper")

colnames(df.OR.CTCF.ABC) = colnames(df.OR.ALL.ABC)
colnames(df.OR.H3K4me3.ABC) = colnames(df.OR.ALL.ABC)


df.OR.ALL.ABC$Element ="ALL"
df.OR.CTCF.ABC$Element = "CTCF"
df.OR.H3K4me3.ABC$Element="H3K4me3"


ors.ABC.regul = rbind(df.OR.ALL.ABC,df.OR.CTCF.ABC,df.OR.H3K4me3.ABC)


OR.H3K27ac.regul.Rao = enrichments.analysis(unique.Regul.Rao,regul.candidates.Rao,Candidates.Signi.H3K27ac,Candidates.NSigni.H3K27ac)
OR.CTCF.regul.Rao = enrichments.analysis(unique.Regul.Rao,regul.candidates.Rao,Candidates.Signi.CTCF,Candidates.NSigni.CTCF)
OR.H3K4me3.regul.Rao = enrichments.analysis(unique.Regul.Rao,regul.candidates.Rao,Candidates.Signi.H3K4me3,Candidates.NSigni.H3K4me3)
OR.DNAse.regul.Rao = enrichments.analysis(unique.Regul.Rao,regul.candidates.Rao,Candidates.Signi.DNAse,Candidates.NSigni.DNAse)
OR.ALL.regul.Rao = enrichments.analysis(unique.Regul.Rao,regul.candidates.Rao,Candidates.Signi.ALL,Candidates.NSigni.ALL)


df.OR.ALL.Rao = data.frame(OR.H3K4me3.regul.Rao$estimate[[1]], OR.H3K4me3.regul.Rao$conf.int[1],OR.H3K4me3.regul.Rao$conf.int[2])
df.OR.CTCF.Rao = data.frame(OR.CTCF.regul.Rao$estimate[[1]],OR.CTCF.regul.Rao$conf.int[1], OR.CTCF.regul.Rao$conf.int[2])
df.OR.H3K4me3.Rao = data.frame(OR.H3K4me3.regul.Rao$estimate[[1]], OR.H3K4me3.regul.Rao$conf.int[1],OR.H3K4me3.regul.Rao$conf.int[2])
df.OR.H3K27ac.Rao = data.frame(OR.H3K27ac.regul.Rao$estimate[[1]],OR.H3K27ac.regul.Rao$conf.int[1], OR.H3K27ac.regul.Rao$conf.int[2])
df.OR.DNAse.Rao = data.frame(OR.DNAse.regul.Rao$estimate[[1]], OR.DNAse.regul.Rao$conf.int[1],OR.DNAse.regul.Rao$conf.int[2])

colnames(df.OR.ALL.Rao) = c("OR", "CI_lower", "CI_upper")

colnames(df.OR.CTCF.Rao) = colnames(df.OR.ALL.Rao)
colnames(df.OR.H3K4me3.Rao) = colnames(df.OR.ALL.Rao)
colnames(df.OR.H3K27ac.Rao) = colnames(df.OR.ALL.Rao)
colnames(df.OR.DNAse.Rao) = colnames(df.OR.ALL.Rao)

df.OR.ALL.Rao$Element ="ALL"
df.OR.CTCF.Rao$Element = "CTCF"
df.OR.H3K4me3.Rao$Element="H3K4me3"
df.OR.H3K27ac.Rao$Element = "H3K27ac"
df.OR.DNAse.Rao$Element="DNAse"


ors.Rao.regul = rbind(df.OR.ALL.Rao,df.OR.CTCF.Rao,df.OR.H3K4me3.Rao,df.OR.H3K27ac.Rao,df.OR.DNAse.Rao)

OR.H3K27ac.regul.DNAse = enrichments.analysis(unique.Regul.DNAse,regul.candidates.DNAse,Candidates.Signi.H3K27ac,Candidates.NSigni.H3K27ac)
OR.CTCF.regul.DNAse = enrichments.analysis(unique.Regul.DNAse,regul.candidates.DNAse,Candidates.Signi.CTCF,Candidates.NSigni.CTCF)
OR.H3K4me3.regul.DNAse = enrichments.analysis(unique.Regul.DNAse,regul.candidates.DNAse,Candidates.Signi.H3K4me3,Candidates.NSigni.H3K4me3)
OR.ALL.regul.DNAse = enrichments.analysis(unique.Regul.DNAse,regul.candidates.DNAse,Candidates.Signi.ALL,Candidates.NSigni.ALL)

df.OR.ALL.DNAse = data.frame(OR.H3K4me3.regul.DNAse$estimate[[1]], OR.H3K4me3.regul.DNAse$conf.int[1],OR.H3K4me3.regul.DNAse$conf.int[2])
df.OR.CTCF.DNAse = data.frame(OR.CTCF.regul.DNAse$estimate[[1]],OR.CTCF.regul.DNAse$conf.int[1], OR.CTCF.regul.DNAse$conf.int[2])
df.OR.H3K4me3.DNAse = data.frame(OR.H3K4me3.regul.DNAse$estimate[[1]], OR.H3K4me3.regul.DNAse$conf.int[1],OR.H3K4me3.regul.DNAse$conf.int[2])
df.OR.H3K27ac.DNAse = data.frame(OR.H3K27ac.regul.DNAse$estimate[[1]],OR.H3K27ac.regul.DNAse$conf.int[1], OR.H3K27ac.regul.DNAse$conf.int[2])


colnames(df.OR.ALL.DNAse) = c("OR", "CI_lower", "CI_upper")

colnames(df.OR.CTCF.DNAse) = colnames(df.OR.ALL.DNAse)
colnames(df.OR.H3K4me3.DNAse) = colnames(df.OR.ALL.DNAse)
colnames(df.OR.H3K27ac.DNAse) = colnames(df.OR.ALL.DNAse)


df.OR.ALL.DNAse$Element ="ALL"
df.OR.CTCF.DNAse$Element = "CTCF"
df.OR.H3K4me3.DNAse$Element="H3K4me3"
df.OR.H3K27ac.DNAse$Element = "H3K27ac"


ors.DNAse.regul = rbind(df.OR.ALL.DNAse,df.OR.CTCF.DNAse,df.OR.H3K4me3.DNAse,df.OR.H3K27ac.DNAse)

ors.ABC.regul$Method = "ABC"
ors.Rao.regul$Method = "Rao"
ors.DNAse.regul$Method = "DNAse"

ors.ALL.regul = rbind(ors.DNAse.regul, ors.ABC.regul, ors.Rao.regul)
ggplot(ors.ALL.regul, aes(x=OR, y=Element)) + geom_pointrange(aes(xmin=CI_lower, xmax=CI_upper, color=Method, group=Method),  position=position_dodge(width=0.35)) +  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed")+
  ylab("Encode Elements")


OR.H3K27ac.Promoters.ABC = enrichments.analysis(unique.Promoters.ABC,promoters.candidates.ABC,Candidates.Signi.H3K27ac,Candidates.NSigni.H3K27ac)
OR.CTCF.Promoters.ABC = enrichments.analysis(unique.Promoters.ABC,promoters.candidates.ABC,Candidates.Signi.CTCF,Candidates.NSigni.CTCF)
OR.H3K4me3.Promoters.ABC = enrichments.analysis(unique.Promoters.ABC,promoters.candidates.ABC,Candidates.Signi.H3K4me3,Candidates.NSigni.H3K4me3)
OR.DNAse.Promoters.ABC = enrichments.analysis(unique.Promoters.ABC,promoters.candidates.ABC,Candidates.Signi.DNAse,Candidates.NSigni.DNAse)
OR.ALL.Promoters.ABC = enrichments.analysis(unique.Promoters.ABC,promoters.candidates.ABC,Candidates.Signi.ALL,Candidates.NSigni.ALL)

df.OR.ALL.ABC = data.frame(OR.H3K4me3.Promoters.ABC$estimate[[1]], OR.H3K4me3.Promoters.ABC$conf.int[1],OR.H3K4me3.Promoters.ABC$conf.int[2])
df.OR.CTCF.ABC = data.frame(OR.CTCF.Promoters.ABC$estimate[[1]],OR.CTCF.Promoters.ABC$conf.int[1], OR.CTCF.Promoters.ABC$conf.int[2])
df.OR.H3K4me3.ABC = data.frame(OR.H3K4me3.Promoters.ABC$estimate[[1]], OR.H3K4me3.Promoters.ABC$conf.int[1],OR.H3K4me3.Promoters.ABC$conf.int[2])
df.OR.H3K27ac.ABC = data.frame(OR.H3K27ac.Promoters.ABC$estimate[[1]],OR.H3K27ac.Promoters.ABC$conf.int[1], OR.H3K27ac.Promoters.ABC$conf.int[2])
df.OR.DNAse.ABC = data.frame(OR.DNAse.Promoters.ABC$estimate[[1]], OR.DNAse.Promoters.ABC$conf.int[1],OR.DNAse.Promoters.ABC$conf.int[2])

colnames(df.OR.ALL.ABC) = c("OR", "CI_lower", "CI_upper")

colnames(df.OR.CTCF.ABC) = colnames(df.OR.ALL.ABC)
colnames(df.OR.H3K4me3.ABC) = colnames(df.OR.ALL.ABC)
colnames(df.OR.H3K27ac.ABC) = colnames(df.OR.ALL.ABC)
colnames(df.OR.DNAse.ABC) = colnames(df.OR.ALL.ABC)

df.OR.ALL.ABC$Element ="ALL"
df.OR.CTCF.ABC$Element = "CTCF"
df.OR.H3K4me3.ABC$Element="H3K4me3"
df.OR.H3K27ac.ABC$Element = "H3K27ac"
df.OR.DNAse.ABC$Element="DNAse"


ors.ABC.Promoters = rbind(df.OR.ALL.ABC,df.OR.CTCF.ABC,df.OR.H3K4me3.ABC,df.OR.H3K27ac.ABC,df.OR.DNAse.ABC)

OR.H3K27ac.Promoters.Rao = enrichments.analysis(unique.Promoters.Rao,promoters.candidates.Rao,Candidates.Signi.H3K27ac,Candidates.NSigni.H3K27ac)
OR.CTCF.Promoters.Rao = enrichments.analysis(unique.Promoters.Rao,promoters.candidates.Rao,Candidates.Signi.CTCF,Candidates.NSigni.CTCF)
OR.H3K4me3.Promoters.Rao = enrichments.analysis(unique.Promoters.Rao,promoters.candidates.Rao,Candidates.Signi.H3K4me3,Candidates.NSigni.H3K4me3)
OR.DNAse.Promoters.Rao = enrichments.analysis(unique.Promoters.Rao,promoters.candidates.Rao,Candidates.Signi.DNAse,Candidates.NSigni.DNAse)
OR.ALL.Promoters.Rao = enrichments.analysis(unique.Promoters.Rao,promoters.candidates.Rao,Candidates.Signi.ALL,Candidates.NSigni.ALL)


df.OR.ALL.Rao = data.frame(OR.H3K4me3.Promoters.Rao$estimate[[1]], OR.H3K4me3.Promoters.Rao$conf.int[1],OR.H3K4me3.Promoters.Rao$conf.int[2])
df.OR.CTCF.Rao = data.frame(OR.CTCF.Promoters.Rao$estimate[[1]],OR.CTCF.Promoters.Rao$conf.int[1], OR.CTCF.Promoters.Rao$conf.int[2])
df.OR.H3K4me3.Rao = data.frame(OR.H3K4me3.Promoters.Rao$estimate[[1]], OR.H3K4me3.Promoters.Rao$conf.int[1],OR.H3K4me3.Promoters.Rao$conf.int[2])
df.OR.H3K27ac.Rao = data.frame(OR.H3K27ac.Promoters.Rao$estimate[[1]],OR.H3K27ac.Promoters.Rao$conf.int[1], OR.H3K27ac.Promoters.Rao$conf.int[2])
df.OR.DNAse.Rao = data.frame(OR.DNAse.Promoters.Rao$estimate[[1]], OR.DNAse.Promoters.Rao$conf.int[1],OR.DNAse.Promoters.Rao$conf.int[2])

colnames(df.OR.ALL.Rao) = c("OR", "CI_lower", "CI_upper")

colnames(df.OR.CTCF.Rao) = colnames(df.OR.ALL.Rao)
colnames(df.OR.H3K4me3.Rao) = colnames(df.OR.ALL.Rao)
colnames(df.OR.H3K27ac.Rao) = colnames(df.OR.ALL.Rao)
colnames(df.OR.DNAse.Rao) = colnames(df.OR.ALL.Rao)

df.OR.ALL.Rao$Element ="ALL"
df.OR.CTCF.Rao$Element = "CTCF"
df.OR.H3K4me3.Rao$Element="H3K4me3"
df.OR.H3K27ac.Rao$Element = "H3K27ac"
df.OR.DNAse.Rao$Element="DNAse"


ors.Rao.Promoters = rbind(df.OR.ALL.Rao,df.OR.CTCF.Rao,df.OR.H3K4me3.Rao,df.OR.H3K27ac.Rao,df.OR.DNAse.Rao)

OR.H3K27ac.Promoters.DNAse = enrichments.analysis(unique.Promoters.DNAse,promoters.candidates.DNAse,Candidates.Signi.H3K27ac,Candidates.NSigni.H3K27ac)
OR.CTCF.Promoters.DNAse = enrichments.analysis(unique.Promoters.DNAse,promoters.candidates.DNAse,Candidates.Signi.CTCF,Candidates.NSigni.CTCF)
OR.H3K4me3.Promoters.DNAse = enrichments.analysis(unique.Promoters.DNAse,promoters.candidates.DNAse,Candidates.Signi.H3K4me3,Candidates.NSigni.H3K4me3)
OR.DNAse.Promoters.DNAse = enrichments.analysis(unique.Promoters.DNAse,promoters.candidates.DNAse,Candidates.Signi.DNAse,Candidates.NSigni.DNAse)
OR.ALL.Promoters.DNAse = enrichments.analysis(unique.Promoters.DNAse,promoters.candidates.DNAse,Candidates.Signi.ALL,Candidates.NSigni.ALL)

df.OR.ALL.DNAse = data.frame(OR.H3K4me3.Promoters.DNAse$estimate[[1]], OR.H3K4me3.Promoters.DNAse$conf.int[1],OR.H3K4me3.Promoters.DNAse$conf.int[2])
df.OR.CTCF.DNAse = data.frame(OR.CTCF.Promoters.DNAse$estimate[[1]],OR.CTCF.Promoters.DNAse$conf.int[1], OR.CTCF.Promoters.DNAse$conf.int[2])
df.OR.H3K4me3.DNAse = data.frame(OR.H3K4me3.Promoters.DNAse$estimate[[1]], OR.H3K4me3.Promoters.DNAse$conf.int[1],OR.H3K4me3.Promoters.DNAse$conf.int[2])
df.OR.H3K27ac.DNAse = data.frame(OR.H3K27ac.Promoters.DNAse$estimate[[1]],OR.H3K27ac.Promoters.DNAse$conf.int[1], OR.H3K27ac.Promoters.DNAse$conf.int[2])
df.OR.DNAse.DNAse = data.frame(OR.DNAse.Promoters.DNAse$estimate[[1]], OR.DNAse.Promoters.DNAse$conf.int[1],OR.DNAse.Promoters.DNAse$conf.int[2])

colnames(df.OR.ALL.DNAse) = c("OR", "CI_lower", "CI_upper")

colnames(df.OR.CTCF.DNAse) = colnames(df.OR.ALL.DNAse)
colnames(df.OR.H3K4me3.DNAse) = colnames(df.OR.ALL.DNAse)
colnames(df.OR.H3K27ac.DNAse) = colnames(df.OR.ALL.DNAse)
colnames(df.OR.DNAse.DNAse) = colnames(df.OR.ALL.DNAse)

df.OR.ALL.DNAse$Element ="ALL"
df.OR.CTCF.DNAse$Element = "CTCF"
df.OR.H3K4me3.DNAse$Element="H3K4me3"
df.OR.H3K27ac.DNAse$Element = "H3K27ac"
df.OR.DNAse.DNAse$Element="DNAse"


ors.DNAse.Promoters = rbind(df.OR.ALL.DNAse,df.OR.CTCF.DNAse,df.OR.H3K4me3.DNAse,df.OR.H3K27ac.DNAse,df.OR.DNAse.DNAse)

ors.ABC.Promoters$Method = "ABC"
ors.Rao.Promoters$Method = "Rao"
ors.DNAse.Promoters$Method = "DNAse"

ors.ALL.Promoters = rbind(ors.DNAse.Promoters, ors.ABC.Promoters, ors.Rao.Promoters)
ggplot(ors.ALL.Promoters, aes(x=OR, y=Element)) + geom_pointrange(aes(xmin=CI_lower, xmax=CI_upper, color=Method, group=Method),  position=position_dodge(width=0.35)) +  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed")+
  ylab("Encode Elements")

#Essential genes
#What's the CRN typology for essential genes
essential.genes = read.csv2("../essential_genes.csv", sep=",")
# essential.genes.DEG = data.table::fread("../degannotation-e.dat")
# 
# unique.essential.genes.DEG = unique(data.frame(essential.genes.DEG[essential.genes.DEG$Refseq=="Homo sapiens","Gene_Ref"]))
unique.essential.genes = unique(essential.genes[, "GENE"])

regul.promoters.ABC$essentialGenes = ifelse(regul.promoters.ABC$TargetGene%in%unique.essential.genes,1,0)
df.Rao$essentialGenes = ifelse(df.Rao$geneSymbol.x%in%unique.essential.genes,1,0)
df.DNase$essentialGenes = ifelse(df.DNase$geneSymbol.x%in%unique.essential.genes,1,0)

#Networks creation
graph_from_ABC = create.Crn(regul.promoters.ABC, method="ABC")
compo.ABC = components(graph_from_ABC)

summary(compo.ABC$csize)

graph_from_Rao = create.Crn(df.Rao, method="Rao")
compo.Rao = components(graph_from_Rao)

summary(compo.Rao$csize)

graph_from_DNAse = create.Crn(df.DNase, method="DNAse")
compo.DNAse = components(graph_from_DNAse)

summary(compo.DNAse$csize)
#Membership adding for subsequent analysis 
regul.promoters.ABC = add.membership(graph_from_ABC, regul.promoters.ABC)
df.Rao = add.membership(graph_from_Rao, df.Rao)
df.DNase = add.membership(graph_from_DNAse, df.DNase)

#Add expression for genes 
expression.genes.RRCs = data.frame(mcols(expression.genes))
expression.genes.RRCs = expression.genes.RRCs[,c("geneSymbol", "median.Expr")]
colnames(expression.genes.RRCs) = c("TargetGene", "Expression")

regul.promoters.ABC = merge(regul.promoters.ABC, expression.genes.RRCs[!duplicated(expression.genes.RRCs$TargetGene),], by="TargetGene", all.x=T)
colnames(expression.genes.RRCs) = c("geneSymbol.x", "Expression")
df.Rao = merge(df.Rao, expression.genes.RRCs[!duplicated(expression.genes.RRCs$geneSymbol.x),], by="geneSymbol.x", all.x=T)
df.DNase = merge(df.DNase, expression.genes.RRCs[!duplicated(expression.genes.RRCs$geneSymbol.x),], by="geneSymbol.x", all.x=T)

summary(as.numeric(table(regul.promoters.ABC$TargetGene)))
summary(as.numeric(table(regul.promoters.ABC$name)))

hist(as.numeric(table(regul.promoters.ABC$TargetGene)))
hist(as.numeric(table(regul.promoters.ABC$name)))

summary(as.numeric(table(df.Rao$geneSymbol.x)))
summary(as.numeric(table(df.Rao$name)))

summary(as.numeric(table(df.DNase$geneSymbol.x)))
summary(as.numeric(table(df.DNase$name)))

wilcox.test(as.numeric(table(regul.promoters.ABC$TargetGene)), as.numeric(table(regul.promoters.ABC$name)), alternative = "greater")
wilcox.test(as.numeric(table(df.Rao$geneSymbol.x)), as.numeric(table(df.Rao$name)), alternative = "greater")
wilcox.test(as.numeric(table(df.DNase$geneSymbol.x)), as.numeric(table(df.DNase$name)), alternative = "greater")

#18-Chromatin States for E007
chromatin.States = import("Epigenetic/E007_18_core_K27ac_dense.bed.gz", format="bed")
unique(chromatin.States$name)

chromatin.States$meta_states = ifelse(chromatin.States$name %in% c("10_EnhA2", "9_EnhA1", "12_ZNF/Rpts", "1_TssA","2_TssFlnk","3_TssFlnkU","4_TssFlnkD", "5_Tx","7_EnhG1", "8_EnhG2"), "Active", ifelse(chromatin.States$name %in% c("6_TxWk" , "15_EnhBiv", "14_TssBiv"),"Weakly_Active", "InactiveOrRepressor" ))
#######################################################################
#Essentials genes 
length(unique(regul.promoters.ABC[regul.promoters.ABC$essentialGenes==1,"TargetGene"]))/length(unique(regul.promoters.ABC$TargetGene))
length(unique(df.Rao[df.Rao$essentialGenes==1,"geneSymbol.x"]))/length(unique(df.Rao$geneSymbol.x))
length(unique(df.DNase[df.DNase$essentialGenes==1,"geneSymbol.x"]))/length(unique(df.DNase$geneSymbol.x))

decompose.ABC = decompose(graph_from_ABC)
decompose.Rao = decompose(graph_from_Rao)
decompose.DNAse = decompose(graph_from_DNAse)

asso.complexity.ABC = add.Complexity.RRCs(regul.promoters.ABC, method="ABC", decompose.ABC)
asso.complexity.Rao = add.Complexity.RRCs(df.Rao, method="Rao", decompose.Rao)
asso.complexity.DNAse = add.Complexity.RRCs(df.DNase, method="DNAse", decompose.DNAse)


asso.complexity.ABC.CS = CS.by.CRN(unique.Regul.ABC, chromatin.States, regul.promoters.ABC,asso.complexity.ABC,method="ABC")
asso.complexity.DNAse.CS = CS.by.CRN(unique.Regul.DNAse, chromatin.States, df.DNase,asso.complexity.DNAse,method="DNAse")
asso.complexity.Rao.CS = CS.by.CRN(unique.Regul.Rao, chromatin.States, df.Rao,asso.complexity.Rao,method="Rao")

asso.complexity.ABC.tmp = asso.complexity.ABC
asso.complexity.Rao.tmp = asso.complexity.Rao
asso.complexity.DNAse.tmp = asso.complexity.DNAse

test.CS.ABC = CS.by.CRN(unique.Regul.ABC, chromatin.States, regul.promoters.ABC,asso.complexity.ABC,method="ABC")
test.CS.Rao = CS.by.CRN(unique.Regul.Rao, chromatin.States, df.Rao,asso.complexity.Rao,method="Rao")
test.CS.DNAse = CS.by.CRN(unique.Regul.DNAse, chromatin.States, df.DNase,asso.complexity.DNAse,method="DNAse")

dt = data.table::rbindlist(
  lapply(test.CS.ABC, function(x) data.table::data.table(t(x))),
  fill = TRUE
)
colSums(!is.na(dt))

g = data.frame(which(!is.na(dt), arr.ind = T))
table(g$col)
agg.g = aggregate(col~row, g, unique)
t.CS.ABC = table(do.call(rbind,lapply(agg.g$col, function(x) paste( unlist(x), collapse=';') )))
sum((sort(t.CS.ABC, decreasing = T)/sum(t.CS.ABC))*100)


dt.Rao = data.table::rbindlist(
  lapply(test.CS.Rao, function(x) data.table::data.table(t(x))),
  fill = TRUE
)
colSums(!is.na(dt.Rao))

g.Rao= data.frame(which(!is.na(dt.Rao), arr.ind = T))
table(g.Rao$col)
agg.g.Rao = aggregate(col~row, g.Rao, unique)
t.CS.Rao = table(do.call(rbind,lapply(agg.g.Rao$col, function(x) paste( unlist(x), collapse=';') )))
(sort(t.CS.Rao, decreasing = T)/sum(t.CS.Rao))*100



dt.DNAse = data.table::rbindlist(
  lapply(test.CS.DNAse, function(x) data.table::data.table(t(x))),
  fill = TRUE
)
colSums(!is.na(dt.DNAse))

g.DNAse = data.frame(which(!is.na(dt.DNAse), arr.ind = T))
table(g.DNAse$col)
agg.g.DNAse = aggregate(col~row, g.DNAse, unique)

t.CS.DNAse = table(do.call(rbind,lapply(agg.g.DNAse$col, function(x) paste( unlist(x), collapse=';') )))

(sort(t.CS.DNAse, decreasing = T)/sum(t.CS.DNAse))*100

length((sort(table(do.call(rbind,lapply(agg.g$col, function(x) paste( unlist(x), collapse=';') ))), decreasing = T)/sum(table(do.call(rbind,lapply(agg.g$col, function(x) paste( unlist(x), collapse=';') )))))*100)
length((sort(table(do.call(rbind,lapply(agg.g.Rao$col, function(x) paste( unlist(x), collapse=';') ))), decreasing = T)/sum(table(do.call(rbind,lapply(agg.g.Rao$col, function(x) paste( unlist(x), collapse=';') )))))*100)
length((sort(table(do.call(rbind,lapply(agg.g.DNAse$col, function(x) paste( unlist(x), collapse=';') ))), decreasing = T)/sum(table(do.call(rbind,lapply(agg.g.DNAse$col, function(x) paste( unlist(x), collapse=';') )))))*100)


asso.complexity.ABC.tmp$nCS = rowSums(!is.na(dt))
asso.complexity.Rao.tmp$nCS = rowSums(!is.na(dt.Rao))
asso.complexity.DNAse.tmp$nCS = rowSums(!is.na(dt.DNAse))



sum(t.CS.Rao[which(unlist(lapply(strsplit(names(t.CS.Rao), ";", fixed=T), length))<=2)]/sum(t.CS.Rao))
t.CS.ABC[which(unlist(lapply(strsplit(names(t.CS.ABC), ";", fixed=T), length))<=2)]/sum(t.CS.ABC)
sum(t.CS.DNAse[which(unlist(lapply(strsplit(names(t.CS.DNAse), ";", fixed=T), length))<=2)]/sum(t.CS.DNAse))

barplot(sort(apply(dt, 2, function(x) sum(!is.na(x))/nrow(dt)), decreasing=T))
barplot(sort(apply(dt.Rao, 2, function(x) sum(!is.na(x))/nrow(dt.Rao)),decreasing = T))
barplot(sort(apply(dt.DNAse, 2, function(x) sum(!is.na(x))/nrow(dt.DNAse)), decreasing = T))

iop = apply(dt, 1, function(x) names(which(!is.na(x)))) 
iop.subset = iop[sapply(iop, function(x) length(x)>=3)]

length(which(lapply(iop,function(x) all(c("18_Quies", "11_EnhWk", "6_TxWk") %in% x ))==T))/length(iop)

iop = apply(dt.Rao, 1, function(x) names(which(!is.na(x)))) 
iop.subset = iop[sapply(iop, function(x) length(x)>=3)]
length(which(lapply(iop,function(x) all(c("18_Quies", "11_EnhWk", "6_TxWk") %in% x ))==T))/length(iop)
iop = apply(dt.DNAse, 1, function(x) names(which(!is.na(x)))) 
iop.subset = iop[sapply(iop, function(x) length(x)>=3)]
length(which(lapply(iop.subset,function(x) all(c("18_Quies", "11_EnhWk", "6_TxWk") %in% x ))==T))/length(iop.subset)

i=ggplot(na.omit(asso.complexity.Rao.CS), aes(x=reorder(CS,CS,
                                          function(x)-length(x))))+
  geom_bar(stat="count", width=0.7, fill="white",aes(color=CS))+xlab("CHROMATIN STATE")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


j=ggplot(na.omit(asso.complexity.Rao.CS), aes(CS, max.prop.CS))+geom_boxplot(aes(color=CS))+xlab("CHROMATIN STATE")+ylab("CHROMATIN STATE PROPORTION")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
k=ggplot(na.omit(asso.complexity.Rao.CS[asso.complexity.Rao.CS$complexity<150,]), aes(CS, complexity))+geom_boxplot(aes(color=CS))+xlab("CHROMATIN STATE")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
l=ggplot(na.omit(asso.complexity.Rao.CS[asso.complexity.Rao.CS$complexity<150,]), aes(complexity, max.prop.CS))+geom_point(aes(color=CS))+geom_smooth(method="gam",formula = y ~ s(x, k = 10))+ylab("CHROMATIN STATE PROPORTION")+
  annotate(geom="text", x=50, y=0.75, label=paste(expression(rho),":",round(cor(na.omit(asso.complexity.Rao.CS)$complexity, na.omit(asso.complexity.Rao.CS)$max.prop.CS, method = "spearman"),2)),
           color="black", parse=T)


multiplot(i,j,k,l,cols=2)


gam.prop.CS.ABC = gam(max.prop.CS~s(complexity,k=100),data=asso.complexity.ABC.CS, method = "REML")
summary(gam.prop.CS.ABC)

kable(tidy(gam.prop.CS),"html")%>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>% 
  save_kable(file = "./gam_complexity_prop.png")

gam.prop.CS.Rao = gam(max.prop.CS~s(complexity,k=10),data=asso.complexity.Rao.CS, method = "REML")
summary(gam.prop.CS.Rao)

gam.prop.CS.DNAse = gam(max.prop.CS~s(complexity,k=10),data=asso.complexity.DNAse.CS, method = "REML")
summary(gam.prop.CS.DNAse)



betareg.CS.prop = betareg(hey.ho~CS,jhgt[jhgt$hey.ho>0&jhgt$hey.ho<1,], link="logit")
summary(betareg.CS.prop)

kable(tidy(betareg.CS.prop),"html")%>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>% 
  save_kable(file = "./betaReg_CS_prop.png")

asso.complexity.DNAse.CS$CS = factor(asso.complexity.DNAse.CS$CS)
asso.complexity.ABC.CS$CS = factor(asso.complexity.ABC.CS$CS)
asso.complexity.Rao.CS$CS = factor(asso.complexity.Rao.CS$CS)

asso.complexity.ABC.CS$CS = relevel(asso.complexity.ABC.CS$CS, ref="11_EnhWk")
asso.complexity.Rao.CS$CS = relevel(asso.complexity.Rao.CS$CS, ref="18_Quies")
asso.complexity.DNAse.CS$CS = relevel(asso.complexity.DNAse.CS$CS, ref="18_Quies")

betareg.CS.prop.Rao = betareg(max.prop.CS~CS,asso.complexity.Rao.CS[asso.complexity.Rao.CS$max.prop.CS>0&asso.complexity.Rao.CS$max.prop.CS<1,], link="logit")
summary(betareg.CS.prop.Rao)

kable(tidy(betareg.CS.prop.Rao),"html")%>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>% 
  save_kable(file = "./betaReg_CS_prop.Rao.png")


betareg.CS.prop.DNAse = betareg(max.prop.CS~CS,asso.complexity.DNAse.CS[asso.complexity.DNAse.CS$max.prop.CS>0&asso.complexity.DNAse.CS$max.prop.CS<1,], link="logit")
summary(betareg.CS.prop.DNAse)

kable(tidy(betareg.CS.prop.DNAse),"html")%>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>% 
  save_kable(file = "./betaReg_CS_prop.DNAse.png")

linearreg.CS.compl = lm(log(complexity)~CS,jhgt)
summary(linearreg.CS.compl)
kable(tidy(linearreg.CS.compl),"html")%>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>% 
  save_kable(file = "./linearReg_prop_complexity.png")

linearreg.CS.compl.Rao = lm(log(complexity)~CS,asso.complexity.Rao.CS)
summary(linearreg.CS.compl.Rao)

kable(tidy(linearreg.CS.compl.Rao),"html")%>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>% 
  save_kable(file = "./linearReg_prop_complexity.Rao.png")


linearreg.CS.compl.DNAse = lm(log(complexity)~CS,asso.complexity.DNAse.CS)
summary(linearreg.CS.compl.DNAse)

kable(tidy(linearreg.CS.compl.DNAse),"html")%>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>% 
  save_kable(file = "./linearReg_prop_complexity.DNAse.png")


#Are the CRNs enriched in essential genes ? 
enrichments.essential.genes.ABC = matrix(c(sum(names(unique.Promoters.ABC)%in%unique.essential.genes),sum(promoters.candidates.ABC$geneSymbol%in%unique.essential.genes),sum(!names(unique.Promoters.ABC)%in%unique.essential.genes),sum(!promoters.candidates.ABC$geneSymbol%in%unique.essential.genes)), ncol=2,nrow=2)
enrichments.essential.genes.Rao = matrix(c(sum(unique.Promoters.Rao$geneSymbol%in%unique.essential.genes),sum(promoters.candidates.Rao$geneSymbol%in%unique.essential.genes),sum(!unique.Promoters.Rao$geneSymbol%in%unique.essential.genes),sum(!promoters.candidates.Rao$geneSymbol%in%unique.essential.genes)), ncol=2,nrow=2)
enrichments.essential.genes.DNAse = matrix(c(sum(unique.Promoters.DNAse$geneSymbol%in%unique.essential.genes),sum(promoters.candidates.DNAse$geneSymbol%in%unique.essential.genes),sum(!unique.Promoters.DNAse$geneSymbol%in%unique.essential.genes),sum(!promoters.candidates.DNAse$geneSymbol%in%unique.essential.genes)), ncol=2,nrow=2)

fisher.test(enrichments.essential.genes.ABC)
fisher.test(enrichments.essential.genes.Rao)
fisher.test(enrichments.essential.genes.DNAse)
plot.OR(fisher.test(enrichments.essential.genes.ABC),fisher.test(enrichments.essential.genes.Rao),fisher.test(enrichments.essential.genes.DNAse))

#Are the essential genes more connected with regulatory elements ? 
summary(as.numeric(table(regul.promoters.ABC[regul.promoters.ABC$essentialGenes==1,"TargetGene"])[table(regul.promoters.ABC[regul.promoters.ABC$essentialGenes==1,"TargetGene"])>0]))
summary(as.numeric(table(regul.promoters.ABC[regul.promoters.ABC$essentialGenes==0,"TargetGene"])[table(regul.promoters.ABC[regul.promoters.ABC$essentialGenes==0,"TargetGene"])>0]))

summary(as.numeric(table(df.Rao[df.Rao$essentialGenes==1,"geneSymbol.x"])[table(df.Rao[df.Rao$essentialGenes==1,"geneSymbol.x"])>0]))
summary(as.numeric(table(df.Rao[df.Rao$essentialGenes==0,"geneSymbol.x"])[table(df.Rao[df.Rao$essentialGenes==0,"geneSymbol.x"])>0]))

summary(as.numeric(table(df.DNase[df.DNase$essentialGenes==1,"geneSymbol.x"])[table(df.DNase[df.DNase$essentialGenes==1,"geneSymbol.x"])>0]))
summary(as.numeric(table(df.DNase[df.DNase$essentialGenes==0,"geneSymbol.x"])[table(df.DNase[df.DNase$essentialGenes==0,"geneSymbol.x"])>0]))

Signal.By.RRCs.ABC = Signal.By.RRCs(regul.promoters.ABC, method = "ABC", aggfun = "90th")
Signal.By.RRCs.Rao = Signal.By.RRCs(df.Rao, method = "Rao", aggfun = "90th")
Signal.By.RRCs.DNAse = Signal.By.RRCs(df.DNase, method = "DNAse", aggfun = "90th")

Expression.CRN.ABC = aggregate(Expression~membership, regul.promoters.ABC, function(x) quantile(x, probs=0.90, na.rm=T))
NinetyPercent.Expression.EssentialGenes.ABC = aggregate(Expression~membership,regul.promoters.ABC[regul.promoters.ABC$essentialGenes==1,c("membership", "essentialGenes", "Expression")], function(x) quantile(x, probs=0.90, na.rm=T))
asso.complexity.ABC = merge(asso.complexity.ABC,merge(Signal.By.RRCs.ABC,merge(Expression.CRN.ABC,NinetyPercent.Expression.EssentialGenes.ABC, by="membership", all.x=T),by="membership", all.x=T),by="membership", all.x=T)

Expression.CRN.Rao = aggregate(Expression~membership, df.Rao, function(x) quantile(x, probs=0.90, na.rm=T))
NinetyPercent.Expression.EssentialGenes.Rao = aggregate(Expression~membership,df.Rao[df.Rao$essentialGenes==1,c("membership", "essentialGenes", "Expression")], function(x) quantile(x, probs=0.90, na.rm=T))
asso.complexity.Rao = merge(asso.complexity.Rao,merge(Signal.By.RRCs.Rao,merge(Expression.CRN.Rao,NinetyPercent.Expression.EssentialGenes.Rao, by="membership", all.x=T),by="membership", all.x=T),by="membership", all.x=T)


Expression.CRN.DNAse = aggregate(Expression~membership, df.DNase, function(x) quantile(x, probs=0.90, na.rm=T))
NinetyPercent.Expression.EssentialGenes.DNAse = aggregate(Expression~membership,df.DNase[df.DNase$essentialGenes==1,c("membership", "essentialGenes", "Expression")], function(x) quantile(x, probs=0.90, na.rm=T))
asso.complexity.DNAse = merge(asso.complexity.DNAse,merge(Signal.By.RRCs.DNAse,merge(Expression.CRN.DNAse,NinetyPercent.Expression.EssentialGenes.DNAse, by="membership", all.x=T),by="membership", all.x=T),by="membership", all.x=T)

asso.complexity.ABC$Expression.y[is.na(asso.complexity.ABC$Expression.y)] = 0
asso.complexity.Rao$Expression.y[is.na(asso.complexity.Rao$Expression.y)] = 0
asso.complexity.DNAse$Expression.y[is.na(asso.complexity.DNAse$Expression.y)] = 0

asso.complexity.ABC$essentielGenes.YESorNO = as.factor(asso.complexity.ABC$essentialGenes>0)
asso.complexity.Rao$essentielGenes.YESorNO = as.factor(asso.complexity.Rao$essentialGenes>0)
asso.complexity.DNAse$essentielGenes.YESorNO = as.factor(asso.complexity.DNAse$essentialGenes>0)

covariates = c("mean_CTCF.PROM",      "mean_CTCF.ENH",       "Expression.x"     ,   "mean_DNAse.PROM","mean_DNAse.ENH",     
               "mean_H3K27ac.PROM",   "mean_H3K27ac.ENH" ,   "mean_p300.PROM"   ,   "mean_p300.ENH"    ,   "mean_H3K4me3.PROM",  
               "mean_H3K4me3.ENH",    "mean_H3K4me1.PROM" ,  "mean_H3K4me1.ENH"  ,  "ScoreFire.PROM"      ,"ScoreFire.ENH",      
               "INS.PROM",            "INS.ENH"    ,         "DI.PROM"    ,         "DI.ENH"          ,    "mean_H3K27me3.PROM",
               "mean_H3K27me3.ENH")

asso.complexity.ABC.FULL$completeCase = ifelse(rowSums(is.na(asso.complexity.ABC.FULL[,covariates]))==0, 1 ,ifelse(rowSums(is.na(asso.complexity.ABC.FULL[,covariates]))==2,1, ifelse(
  rowSums(is.na(asso.complexity.ABC.FULL[,covariates]))==6, 2, ifelse(rowSums(is.na(asso.complexity.ABC.FULL[,covariates]))==8,3, ifelse(
    rowSums(is.na(asso.complexity.ABC.FULL[,covariates]))==10,4, ifelse(rowSums(is.na(asso.complexity.ABC.FULL[,covariates]))==12,5, ifelse(rowSums(is.na(asso.complexity.ABC.FULL[,covariates]))==14,6, ifelse(rowSums(is.na(asso.complexity.ABC.FULL[,covariates]))==16,7, ifelse(rowSums(is.na(asso.complexity.ABC.FULL[,covariates]))==18,8, ifelse(rowSums(is.na(asso.complexity.ABC.FULL[,covariates]))==19,9, ifelse(rowSums(is.na(asso.complexity.ABC.FULL[,covariates]))==21,10,11)))))))))))

multi.model = multinom(completeCase~log(complexity),asso.complexity.ABC.FULL)
weight.multi = fitted(multi.model)



# asso.complexity.ABC.FULL$completeCase = as.numeric(apply(asso.complexity.ABC.FULL[,c("mean_CTCF.PROM",      "mean_CTCF.ENH",       "Expression.x"     ,   "mean_DNAse.PROM","mean_DNAse.ENH",     
#                                                                                       "mean_H3K27ac.PROM",   "mean_H3K27ac.ENH" ,   "mean_p300.PROM"   ,   "mean_p300.ENH"    ,   "mean_H3K4me3.PROM",  
#                                                                                       "mean_H3K4me3.ENH",    "mean_H3K4me1.PROM" ,  "mean_H3K4me1.ENH"  ,  "ScoreFire.PROM"      ,"ScoreFire.ENH",      
#                                                                                       "INS.PROM",            "INS.ENH"    ,         "DI.PROM"    ,         "DI.ENH"          ,    "mean_H3K27me3.PROM",
#                                                                                       "mean_H3K27me3.ENH",   "Expression.y")], 1,function(x) !any(is.na(x))))

estimate.weights = glm(completeCase~log(complexity), data=asso.complexity.ABC.FULL, family = binomial(link="logit"))
weights = predict(estimate.weights, type="response")




model.3D.epige = lm(log(complexity)~ mean_CTCF.PROM+mean_CTCF.ENH+mean_DNAse.PROM+mean_DNAse.ENH+mean_H3K27ac.PROM+mean_H3K27ac.ENH+mean_p300.PROM+mean_p300.ENH+mean_H3K4me3.PROM+  
                mean_H3K4me3.ENH+mean_H3K4me1.PROM+mean_H3K4me1.ENH+ScoreFire.PROM+ScoreFire.ENH+      
                INS.PROM+INS.ENH+DI.PROM+DI.ENH+mean_H3K27me3.PROM+ 
                mean_H3K27me3.ENH+essentielGenes.YESorNO, asso.complexity.ABC.FULL, weights = ifelse(asso.complexity.ABC.FULL$completeCase==1, asso.complexity.ABC.FULL$completeCase/weight.multi,0))
summary(model.3D.epige);AIC(model.3D.epige)
plot(model.3D.epige)




#imputation for missing values:
regul.promoters.ABC.imput = imput.individual.elements(unique.Promoters.ABC,unique.Regul.ABC, regul.promoters.ABC, method="ABC")$all
df.Rao.imput = imput.individual.elements(unique.Promoters.Rao,unique.Regul.Rao, df.Rao, method="Rao")$all
df.DNase.imput = imput.individual.elements(unique.Promoters.DNAse,unique.Regul.DNAse, df.DNase, method="DNAse")$all

Signal.By.RRCs.ABC = Signal.By.RRCs(regul.promoters.ABC.imput, method="ABC", aggfun = "90th")
Signal.By.RRCs.Rao = Signal.By.RRCs(df.Rao.imput, method="Rao", aggfun = "90th")
Signal.By.RRCs.DNAse = Signal.By.RRCs(df.DNase.imput, method="DNAse", aggfun = "90th")

Expression.CRN.ABC.imput = aggregate(Expression~membership, regul.promoters.ABC.imput, function(x) quantile(x, probs=0.90, na.rm=T))
Expression.CRN.Rao.imput = aggregate(Expression~membership, df.Rao.imput, function(x) quantile(x, probs=0.90, na.rm=T))
Expression.CRN.DNAse.imput = aggregate(Expression~membership, df.DNase.imput, function(x) quantile(x, probs=0.90, na.rm=T))

NinetyPercent.Expression.EssentialGenes.ABC.imput = aggregate(Expression~membership,regul.promoters.ABC.imput[regul.promoters.ABC.imput$essentialGenes==1,c("membership", "essentialGenes", "Expression")], function(x) quantile(x, probs=0.90, na.rm=T))
NinetyPercent.Expression.EssentialGenes.Rao.imput = aggregate(Expression~membership,df.Rao.imput[df.Rao.imput$essentialGenes==1,c("membership", "essentialGenes", "Expression")], function(x) quantile(x, probs=0.90, na.rm=T))
NinetyPercent.Expression.EssentialGenes.DNAse.imput = aggregate(Expression~membership,df.DNase.imput[df.DNase.imput$essentialGenes==1,c("membership", "essentialGenes", "Expression")], function(x) quantile(x, probs=0.90, na.rm=T))



asso.complexity.ABC.FULL.imput=merge(merge(Signal.By.RRCs.ABC,Expression.CRN.ABC.imput , by="membership", all.x=T), NinetyPercent.Expression.EssentialGenes.ABC.imput, by="membership", all.x=T)
asso.complexity.ABC.FULL.imput$Expression.y[is.na(asso.complexity.ABC.FULL.imput$Expression.y)] = 0
asso.complexity.ABC.FULL.imput = merge(asso.complexity.ABC.FULL.imput,asso.complexity.ABC[,c("membership", "essentialGenes", "complexity", "TargetGene", "prop.EssentialGenes")], by="membership", all.x=T)
asso.complexity.ABC.FULL.imput$essentielGenes.YESorNO = as.factor(asso.complexity.ABC.FULL.imput$essentialGenes>0)

model.imput.ABC = lm(log(complexity)~ mean_CTCF.PROM+mean_CTCF.ENH+mean_DNAse.PROM+mean_DNAse.ENH+mean_H3K27ac.PROM+mean_H3K27ac.ENH+mean_p300.PROM+mean_p300.ENH+mean_H3K4me3.PROM+  
     mean_H3K4me3.ENH+mean_H3K4me1.PROM+mean_H3K4me1.ENH+ScoreFire.PROM+ScoreFire.ENH+      
     INS.PROM+INS.ENH+DI.PROM+DI.ENH+mean_H3K27me3.PROM+ 
     mean_H3K27me3.ENH+essentielGenes.YESorNO, asso.complexity.ABC.FULL.imput)


summary(model.imput.ABC)

kable(tidy(model.imput.ABC),"html")%>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>% 
  save_kable(file = "./linearReg_all_complexity.ABC.png")


asso.complexity.Rao.FULL.imput=merge(merge(Signal.By.RRCs.Rao,Expression.CRN.Rao.imput , by="membership", all.x=T), NinetyPercent.Expression.EssentialGenes.Rao.imput, by="membership", all.x=T)
asso.complexity.Rao.FULL.imput$Expression.y[is.na(asso.complexity.Rao.FULL.imput$Expression.y)] = 0

asso.complexity.Rao.FULL.imput = merge(asso.complexity.Rao.FULL.imput,asso.complexity.Rao[,c("membership", "essentialGenes", "complexity", "geneSymbol.x", "prop.EssentialGenes")], by="membership", all.x=T)
asso.complexity.Rao.FULL.imput$essentielGenes.YESorNO = as.factor(asso.complexity.Rao.FULL.imput$essentialGenes>0)


model.imput.Rao = lm(log(complexity)~ mean_CTCF.x+mean_CTCF.y+mean_DNAse.x+mean_DNAse.y+mean_H3K27ac.x+mean_H3K27ac.y+mean_p300.x+mean_p300.y+mean_H3K4me3.x+  
                       mean_H3K4me3.y+mean_H3K4me1.x+mean_H3K4me1.y+mean_ScoreFire.x+mean_ScoreFire.y+      
                       mean_INS.x+mean_INS.y+mean_DI.x+mean_DI.y+mean_H3K27me3.x+ 
                       mean_H3K27me3.y+essentielGenes.YESorNO, asso.complexity.Rao.FULL.imput)

plot(model.imput.Rao)
summary(model.imput.Rao)

kable(tidy(model.imput.Rao),"html")%>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>% 
  save_kable(file = "./linearReg_all_complexity.Rao.png")

asso.complexity.DNAse.FULL.imput=merge(merge(Signal.By.RRCs.DNAse,Expression.CRN.DNAse.imput , by="membership", all.x=T), NinetyPercent.Expression.EssentialGenes.DNAse.imput, by="membership", all.x=T)
asso.complexity.DNAse.FULL.imput$Expression.y[is.na(asso.complexity.DNAse.FULL.imput$Expression.y)] = 0

asso.complexity.DNAse.FULL.imput = merge(asso.complexity.DNAse.FULL.imput,asso.complexity.DNAse[,c("membership", "essentialGenes", "complexity", "geneSymbol.x", "prop.EssentialGenes")], by="membership", all.x=T)
asso.complexity.DNAse.FULL.imput$essentielGenes.YESorNO = as.factor(asso.complexity.DNAse.FULL.imput$essentialGenes>0)


model.imput.DNAse = lm(log(complexity)~ mean_CTCF.x+mean_CTCF.y+mean_DNAse.x+mean_DNAse.y+mean_H3K27ac.x+mean_H3K27ac.y+mean_p300.x+mean_p300.y+mean_H3K4me3.x+  
                       mean_H3K4me3.y+mean_H3K4me1.x+mean_H3K4me1.y+mean_ScoreFire.x+mean_ScoreFire.y+      
                       mean_INS.x+mean_INS.y+mean_DI.x+mean_DI.y+mean_H3K27me3.x+ 
                       mean_H3K27me3.y+essentielGenes.YESorNO, asso.complexity.DNAse.FULL.imput)

plot(model.imput.DNAse)
summary(model.imput.DNAse)

kable(tidy(model.imput.DNAse),"html")%>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>% 
  save_kable(file = "./linearReg_all_complexity.DNAse.png")


model.ABC = CI.regression(model.imput.ABC, "ABC")
model.Rao = CI.regression(model.imput.Rao, "Rao")
model.DNAse = CI.regression(model.imput.DNAse, "DNAse")

model.Rao$name = model.ABC$name
model.DNAse$name = model.ABC$name

#CI with robust variance
ggplot(rbind(model.ABC,model.Rao, model.DNAse), aes(x=Estimate, y=name, color=Method)) + 
  geom_pointrange(aes(xmin=CI.l, xmax=CI.u), position=position_dodge(width=0.5),size=0.2)


model.expression.essential.genes.ABC = gam(Expression.y~s(complexity), data=asso.complexity.ABC.FULL.imput[!is.na(asso.complexity.ABC.FULL.imput$Expression.y),], method="REML")
summary(model.expression.essential.genes.ABC)
cor(asso.complexity.ABC.FULL.imput[!is.na(asso.complexity.ABC.FULL.imput$Expression.y),"complexity"],asso.complexity.ABC.FULL.imput[!is.na(asso.complexity.ABC.FULL.imput$Expression.y),"Expression.y"], method = "spearman")


model.expression.essential.genes.Rao = gam(Expression.y~s(complexity), data=asso.complexity.Rao.FULL.imput[!is.na(asso.complexity.Rao.FULL.imput$Expression.y),], method="REML")
summary(model.expression.essential.genes.Rao)
cor(asso.complexity.Rao.FULL.imput[!is.na(asso.complexity.Rao.FULL.imput$Expression.y),"complexity"],asso.complexity.Rao.FULL.imput[!is.na(asso.complexity.Rao.FULL.imput$Expression.y),"Expression.y"], method = "spearman")


model.expression.essential.genes.DNAse = gam(Expression.y~s(complexity), data=asso.complexity.DNAse.FULL.imput[!is.na(asso.complexity.DNAse.FULL.imput$Expression.y),], method="REML")
summary(model.expression.essential.genes.DNAse)
cor(asso.complexity.DNAse.FULL.imput[!is.na(asso.complexity.DNAse.FULL.imput$Expression.y),"complexity"],asso.complexity.DNAse.FULL.imput[!is.na(asso.complexity.DNAse.FULL.imput$Expression.y),"Expression.y"], method = "spearman")


plot.model.expression.essential.genes.ABC = ggplot(na.omit(asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$Expression.y>0,]), aes(complexity, Expression.y))+geom_point()+xlab("complexity")+ylab("90th percentile of essential genes expression")+geom_smooth(method="gam",formula = y ~ s(x))+
  annotate(geom="text", x=1000, y=10, label=paste(expression(rho),":",round(cor(asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$Expression.y>0,"complexity"],asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$Expression.y>0,"Expression.y"], method = "spearman"),2)),
           color="black", parse=T)
plot.model.expression.essential.genes.Rao = ggplot(na.omit(asso.complexity.Rao.FULL.imput[asso.complexity.Rao.FULL.imput$Expression.y>0,]), aes(complexity, Expression.y))+geom_point()+xlab("complexity")+ylab("90th percentile of essential genes expression")+geom_smooth(method="gam",formula = y ~ s(x))+
  annotate(geom="text", x=50, y=10, label=paste(expression(rho),":",round(cor(asso.complexity.Rao.FULL.imput[asso.complexity.Rao.FULL.imput$Expression.y>0,"complexity"],asso.complexity.Rao.FULL.imput[asso.complexity.Rao.FULL.imput$Expression.y>0,"Expression.y"], method = "spearman"),2)),
           color="black", parse=T)
plot.model.expression.essential.genes.DNAse  = ggplot(na.omit(asso.complexity.DNAse.FULL.imput[asso.complexity.DNAse.FULL.imput$Expression.y>0,]), aes(complexity, Expression.y))+geom_point()+xlab("complexity")+ylab("90th percentile of essential genes expression")+geom_smooth(method="gam",formula = y ~ s(x))+
  annotate(geom="text", x=50, y=10, label=paste(expression(rho),":",round(cor(asso.complexity.DNAse.FULL.imput[asso.complexity.DNAse.FULL.imput$Expression.y>0,"complexity"],asso.complexity.DNAse.FULL.imput[asso.complexity.DNAse.FULL.imput$Expression.y>0,"Expression.y"], method = "spearman"),2)),
           color="black", parse=T)

multiplot(plot.model.expression.essential.genes.ABC,plot.model.expression.essential.genes.Rao,plot.model.expression.essential.genes.DNAse)

#Split Enhancers in large CRNs vs enhancers in little CRNs:

split.CRNs.regul = function(df.by.CRN, GRange,df.all, threshold=0, method=c("ABC","Rao", "DNAse")){
  if(method=="ABC"){
    membership.large.CRN.ABC = df.by.CRN[df.by.CRN$TargetGene >=threshold, "membership"]
    regul.large.CRN.ABC = GRange[names(GRange) %in%unique(df.all[df.all$membership%in%membership.large.CRN.ABC, "name" ])]
    regul.small.CRN.ABC = GRange[!names(GRange) %in%unique(df.all[df.all$membership%in%membership.large.CRN.ABC, "name" ])]
    
    df.small = data.frame("CHR"=seqnames(regul.small.CRN.ABC), "start"=start(regul.small.CRN.ABC), "end"=end(regul.small.CRN.ABC))
    df.large = data.frame("CHR"=seqnames(regul.large.CRN.ABC), "start"=start(regul.large.CRN.ABC), "end"=end(regul.large.CRN.ABC))
    return(list("small"=regul.small.CRN.ABC,"large"=regul.large.CRN.ABC))
  }
  
  else{
    membership.large.CRN = df.by.CRN[df.by.CRN$geneSymbol.x >=threshold, "membership"]
    
    names(GRange) = paste0(seqnames(GRange),":", start(GRange),":", end(GRange))
    
    regul.large.CRN = GRange[names(GRange) %in%unique(df.all[df.all$membership%in%membership.large.CRN, "name" ])]
    regul.small.CRN = GRange[!names(GRange) %in%unique(df.all[df.all$membership%in%membership.large.CRN, "name" ])]
    df.small = data.frame("CHR"=seqnames(regul.small.CRN), "start"=start(regul.small.CRN), "end"=end(regul.small.CRN))
    df.large = data.frame("CHR"=seqnames(regul.large.CRN), "start"=start(regul.large.CRN), "end"=end(regul.large.CRN))
    return(list("small"=regul.small.CRN,"large"=regul.large.CRN))
  }
}

split.CRNs.prom = function(df.by.CRN, GRange,df.all, threshold=0, method=c("ABC","Rao", "DNAse")){
  if(method=="ABC"){
    membership.medium.CRN.ABC = df.by.CRN[df.by.CRN$TargetGene > threshold & df.by.CRN$TargetGene <=25, "membership"]
    membership.small.CRN.ABC = df.by.CRN[df.by.CRN$TargetGene <= threshold , "membership"]
    membership.large.CRN.ABC = df.by.CRN[df.by.CRN$TargetGene >25, "membership"]
    
    regul.large.CRN.ABC = GRange[names(GRange) %in%unique(df.all[df.all$membership%in%membership.large.CRN.ABC, "TargetGene" ])]
    regul.small.CRN.ABC = GRange[names(GRange) %in%unique(df.all[df.all$membership%in%membership.small.CRN.ABC, "TargetGene" ])]
    regul.medium.CRN.ABC = GRange[names(GRange) %in%unique(df.all[df.all$membership%in%membership.medium.CRN.ABC, "TargetGene" ])]
    
    GeneSet.L = data.frame(matrix(names(regul.large.CRN.ABC), ncol=1))
    GeneCoord.L = cbind(GeneSet.L, data.frame(regul.large.CRN.ABC)[,1:3])
    
    GeneSet.S = data.frame(matrix(names(regul.small.CRN.ABC), ncol=1))
    GeneCoord.S = cbind(GeneSet.S, data.frame(regul.small.CRN.ABC)[,1:3])
    
    GeneSet.M = data.frame(matrix(names(regul.medium.CRN.ABC), ncol=1))
    GeneCoord.M = cbind(GeneSet.M, data.frame(regul.medium.CRN.ABC)[,1:3])
    
    colnames(GeneCoord.L) = c("GENE", "CHR", "START", "END")
    colnames(GeneCoord.S) = c("GENE", "CHR", "START", "END")
    colnames(GeneCoord.M) = c("GENE", "CHR", "START", "END")
    
    return(list("large" = list("GSet" = na.omit(GeneSet.L), "GCoor" = na.omit(GeneCoord.L)), "small" = list("GSet" = na.omit(GeneSet.S), "GCoor" = na.omit(GeneCoord.S)), "medium" = list("GSet" = na.omit(GeneSet.M), "GCoor" = na.omit(GeneCoord.M))) )
    
    # df.small = data.frame("CHR"=seqnames(regul.small.CRN.ABC), "start"=start(regul.small.CRN.ABC), "end"=end(regul.small.CRN.ABC))
    # df.large = data.frame("CHR"=seqnames(regul.large.CRN.ABC), "start"=start(regul.large.CRN.ABC), "end"=end(regul.large.CRN.ABC))
    # return(list("small"=df.small,"large"=df.large))
  }
  
  else{
    membership.large.CRN = df.by.CRN[df.by.CRN$geneSymbol.x >=threshold, "membership"]
    
    names(GRange) = GRange$geneSymbol
    
    regul.large.CRN = GRange[names(GRange) %in%unique(df.all[df.all$membership%in%membership.large.CRN, "geneSymbol.x" ])]
    regul.small.CRN = GRange[!names(GRange) %in%unique(df.all[df.all$membership%in%membership.large.CRN, "geneSymbol.x" ])]
    
    GeneSet.L = data.frame(matrix(names(regul.large.CRN), ncol=1))
    GeneCoord.L = cbind(GeneSet.L, data.frame(regul.large.CRN)[,1:3])
    
    GeneSet.S = data.frame(matrix(names(regul.small.CRN), ncol=1))
    GeneCoord.S = cbind(GeneSet.S, data.frame(regul.small.CRN)[,1:3])
    
    
    colnames(GeneCoord.L) = c("GENE", "CHR", "START", "END")
    colnames(GeneCoord.S) = c("GENE", "CHR", "START", "END")
    
    return(list("large" = list("GSet" = na.omit(GeneSet.L), "GCoor" = na.omit(GeneCoord.L)), "small" = list("GSet" = na.omit(GeneSet.S), "GCoor" = na.omit(GeneCoord.S))))
    
    
    # df.small = data.frame("CHR"=seqnames(regul.small.CRN), "start"=start(regul.small.CRN), "end"=end(regul.small.CRN))
    # df.large = data.frame("CHR"=seqnames(regul.large.CRN), "start"=start(regul.large.CRN), "end"=end(regul.large.CRN))
    # return(list("small"=df.small,"large"=df.large))
  }
}




split.DNAse.regul = split.CRNs.regul(asso.complexity.DNAse.FULL.imput, unique.Regul.DNAse, df.DNase, threshold = 3, method="DNAse")
split.Rao.regul = split.CRNs.regul(asso.complexity.Rao.FULL.imput, unique.Regul.Rao, df.Rao, threshold = 3, method="Rao")

split.DNAse.prom = split.CRNs.prom(asso.complexity.DNAse.FULL.imput, unique.Promoters.DNAse, df.DNase, threshold = 3, method="DNAse")
split.Rao.prom = split.CRNs.prom(asso.complexity.Rao.FULL.imput, unique.Promoters.Rao, df.Rao, threshold = 3, method="Rao")
split.ABC.prom = split.CRNs.prom(asso.complexity.ABC.FULL.imput, unique.Promoters.ABC, regul.promoters.ABC, threshold = 3, method="ABC")


membership.medium.CRN.ABC = asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$TargetGene >3 & asso.complexity.ABC.FULL.imput$TargetGene <=25, "membership"]
membership.small.CRN.ABC = asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$TargetGene <=3 , "membership"]
membership.large.CRN.ABC = asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$TargetGene >25, "membership"]

regul.medium.CRN.ABC = unique.Regul.ABC[names(unique.Regul.ABC) %in%unique(regul.promoters.ABC[regul.promoters.ABC$membership%in%membership.medium.CRN.ABC, "name" ])]
regul.small.CRN.ABC = unique.Regul.ABC[names(unique.Regul.ABC) %in%unique(regul.promoters.ABC[regul.promoters.ABC$membership%in%membership.small.CRN.ABC, "name" ])]
regul.large.CRN.ABC= unique.Regul.ABC[names(unique.Regul.ABC) %in%unique(regul.promoters.ABC[regul.promoters.ABC$membership%in%membership.large.CRN.ABC, "name" ])]

df.small = data.frame("CHR"=seqnames(regul.small.CRN.ABC), "start"=start(regul.small.CRN.ABC), "end"=end(regul.small.CRN.ABC))
df.large = data.frame("CHR"=seqnames(regul.large.CRN.ABC), "start"=start(regul.large.CRN.ABC), "end"=end(regul.large.CRN.ABC))
df.medium =  data.frame("CHR"=seqnames(regul.medium.CRN.ABC), "start"=start(regul.medium.CRN.ABC), "end"=end(regul.medium.CRN.ABC))

write.table(df.large, "Large_ENH_ABC.bed", sep="\t", row.names = F)
write.table(df.small, "small_ENH_ABC.bed", sep="\t", row.names = F)
write.table(df.medium, "medium_ENH_ABC.bed", sep="\t", row.names = F)


membership.medium.CRN.ABC = asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$TargetGene >3 & asso.complexity.ABC.FULL.imput$TargetGene <=25, "membership"]
membership.small.CRN.ABC = asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$TargetGene <=3 , "membership"]
membership.large.CRN.ABC = asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$TargetGene >25, "membership"]

regul.medium.CRN.ABC = unique.Regul.ABC[names(unique.Regul.ABC) %in%unique(regul.promoters.ABC[regul.promoters.ABC$membership%in%membership.medium.CRN.ABC, "name" ])]
regul.small.CRN.ABC = unique.Regul.ABC[names(unique.Regul.ABC) %in%unique(regul.promoters.ABC[regul.promoters.ABC$membership%in%membership.small.CRN.ABC, "name" ])]
regul.large.CRN.ABC= unique.Regul.ABC[names(unique.Regul.ABC) %in%unique(regul.promoters.ABC[regul.promoters.ABC$membership%in%membership.large.CRN.ABC, "name" ])]

df.small = data.frame("CHR"=seqnames(regul.small.CRN.ABC), "start"=start(regul.small.CRN.ABC), "end"=end(regul.small.CRN.ABC))
df.large = data.frame("CHR"=seqnames(regul.large.CRN.ABC), "start"=start(regul.large.CRN.ABC), "end"=end(regul.large.CRN.ABC))
df.medium =  data.frame("CHR"=seqnames(regul.medium.CRN.ABC), "start"=start(regul.medium.CRN.ABC), "end"=end(regul.medium.CRN.ABC))



write.table(split.DNAse.regul$large, "Large_ENH_DNAse.bed", sep="\t", row.names = F)
write.table(split.DNAse.regul$small, "small_ENH_DNAse.bed", sep="\t", row.names = F)
write.table(split.Rao.regul$large, "Large_ENH_Rao.bed", sep="\t", row.names = F)
write.table(split.Rao.regul$small, "small_ENH_Rao.bed", sep="\t", row.names = F)

write.table(split.DNAse.prom$large$GSet, "NEU_DNAse.LGeneSet", sep="\t", row.names = F, col.names = F)
write.table(split.DNAse.prom$large$GCoor, "NEU_DNAse_LGenecoord.txt", sep="\t", row.names = F, col.names = T)

write.table(split.DNAse.prom$small$GSet, "NEU_DNAse.SGeneSet", sep="\t", row.names = F, col.names = F)
write.table(split.DNAse.prom$small$GCoor, "NEU_DNAse_SGenecoord.txt", sep="\t", row.names = F, col.names = T)


write.table(split.Rao.prom$large$GSet, "NEU_Rao.LGeneSet", sep="\t", row.names = F, col.names = F)
write.table(split.Rao.prom$large$GCoor, "NEU_Rao_LGenecoord.txt", sep="\t", row.names = F, col.names = T)

write.table(split.Rao.prom$small$GSet, "NEU_Rao.SGeneSet", sep="\t", row.names = F, col.names = F)
write.table(split.Rao.prom$small$GCoor, "NEU_Rao_SGenecoord.txt", sep="\t", row.names = F, col.names = T)

write.table(split.ABC.prom$large$GSet, "NEU_ABC.LGeneSet", sep="\t", row.names = F, col.names = F)
write.table(split.ABC.prom$large$GCoor, "NEU_ABC_LGenecoord.txt", sep="\t", row.names = F, col.names = T)

write.table(split.ABC.prom$small$GSet, "NEU_ABC.SGeneSet", sep="\t", row.names = F, col.names = F)
write.table(split.ABC.prom$small$GCoor, "NEU_ABC_SGenecoord.txt", sep="\t", row.names = F, col.names = T)

write.table(split.ABC.prom$medium$GSet, "NEU_ABC.MGeneSet", sep="\t", row.names = F, col.names = F)
write.table(split.ABC.prom$medium$GCoor, "NEU_ABC_MGenecoord.txt", sep="\t", row.names = F, col.names = T)


#Do Promoter caracteristics influence the relationship behaviour with enhancers ? 
NRegul.Prom.ABC = data.frame(table(regul.promoters.ABC.imput$TargetGene))
colnames(NRegul.Prom.ABC) = c("TargetGene", "NRegul")
NRegul.Prom.Rao = data.frame(table(df.Rao.imput$geneSymbol.x))
colnames(NRegul.Prom.Rao) = c("geneSymbol.x", "NRegul")
NRegul.Prom.DNAse = data.frame(table(df.DNase.imput$geneSymbol.x))
colnames(NRegul.Prom.DNAse) = c("geneSymbol.x", "NRegul")

df.asso.NRegul.PROM.ABC = unique(merge(regul.promoters.ABC.imput, NRegul.Prom.ABC, by="TargetGene", all.x=T)[,c("TargetGene", "NRegul", "DI.PROM", "ScoreFire.PROM","INS.PROM", "mean_CTCF.PROM","mean_H3K27ac.PROM","mean_DNAse.PROM","mean_H3K4me1.PROM","mean_H3K4me3.PROM","mean_p300.PROM","mean_H3K27me3.PROM","essentialGenes")])
df.asso.NRegul.PROM.Rao = unique(merge(df.Rao.imput, NRegul.Prom.Rao, by="geneSymbol.x", all.x=T)[,c("geneSymbol.x", "NRegul", "mean_DI.x", "mean_ScoreFire.x","mean_INS.x", "mean_CTCF.x","mean_H3K27ac.x","mean_DNAse.x","mean_H3K4me1.x","mean_H3K4me3.x","mean_p300.x","mean_H3K27me3.x","essentialGenes")])
df.asso.NRegul.PROM.DNAse = unique(merge(df.DNase.imput, NRegul.Prom.DNAse, by="geneSymbol.x", all.x=T)[,c("geneSymbol.x", "NRegul", "mean_DI.x", "mean_ScoreFire.x","mean_INS.x", "mean_CTCF.x","mean_H3K27ac.x","mean_DNAse.x","mean_H3K4me1.x","mean_H3K4me3.x","mean_p300.x","mean_H3K27me3.x","essentialGenes")])

asso.NRegul.PROM.ABC = glm(NRegul~DI.PROM+ScoreFire.PROM+INS.PROM+mean_CTCF.PROM+mean_H3K27ac.PROM+mean_DNAse.PROM+mean_H3K4me1.PROM+
                                 mean_H3K4me3.PROM+mean_p300.PROM+mean_H3K27me3.PROM+essentialGenes, data=df.asso.NRegul.PROM.ABC, family = poisson(link="log"))
summary(asso.NRegul.PROM.ABC)

asso.NRegul.PROM.Rao = glm(NRegul~mean_DI.x+mean_ScoreFire.x+mean_INS.x+mean_CTCF.x+mean_H3K27ac.x+mean_DNAse.x+mean_H3K4me1.x+
                               mean_H3K4me3.x+mean_p300.x+mean_H3K27me3.x+essentialGenes, data=df.asso.NRegul.PROM.Rao, family = poisson(link="log"))
summary(asso.NRegul.PROM.Rao)

asso.NRegul.PROM.DNAse = glm(NRegul~mean_DI.x+mean_ScoreFire.x+mean_INS.x+mean_CTCF.x+mean_H3K27ac.x+mean_DNAse.x+mean_H3K4me1.x+
                             mean_H3K4me3.x+mean_p300.x+mean_H3K27me3.x+essentialGenes, data=df.asso.NRegul.PROM.DNAse, family = poisson(link="log"))
summary(asso.NRegul.PROM.DNAse)

#Integration of TADs for correction of auto-correlation in residuals 
# TADs.Ins = read.table("/home/loic/Documents/HiC/data/4Script/NEU/TADs/TADs.INS.bedpe", header=T)
# GRanges.TADs = GRanges(seqnames=TADs.Ins$chrom1, ranges=IRanges(start=TADs.Ins$start1, end=TADs.Ins$end1, names=1:nrow(TADs.Ins)))
# 
# regul.promoters.ABC$minStart = apply(regul.promoters.ABC[,c("start", "startProm")], 1, min)
# regul.promoters.ABC$maxEnd = apply(regul.promoters.ABC[,c("end", "endProm")], 1, max)
# 
# start.cluster = aggregate(minStart~membership, regul.promoters.ABC, min)
# end.cluster = aggregate(maxEnd~membership, regul.promoters.ABC, max)
# 
# df.start.end = unique(merge(merge(start.cluster, end.cluster, by="membership"), regul.promoters.ABC[,c("chr.x", "membership")], by="membership"))
# 
# GRanges.cluster = GRanges(seqnames = df.start.end$chr, ranges = IRanges(start=df.start.end$minStart, end=df.start.end$maxEnd, names=df.start.end$compo.graph.membership ) )
# 
# member2TADs = data.frame(findOverlaps(GRanges.cluster, GRanges.TADs)[queryHits(findOverlaps(GRanges.cluster, GRanges.TADs))%in%which(countOverlaps(GRanges.cluster, GRanges.TADs)==1)])
# colnames(member2TADs) = c("membership", "TADs")
# 
# asso.complexity.ABC.FULL.imput = merge(asso.complexity.ABC.FULL.imput, member2TADs, by="membership", all.x=T)
# 
# asso.complexity.ABC.FULL.withTADs = asso.complexity.ABC.FULL.imput[!is.na(asso.complexity.ABC.FULL.imput$TADs), ]
# asso.complexity.ABC.FULL.withTADs$essentielGenes.YESorNO = as.numeric(asso.complexity.ABC.FULL.withTADs$essentialGenes>0)
# 
# 
# library(gee)
# osub = as.numeric(order(asso.complexity.ABC.FULL.withTADs$TADs))
# asso.complexity.ABC.FULL.withTADs = asso.complexity.ABC.FULL.withTADs[osub,]
# 
# gee.model.imput = gee(log(complexity)~ mean_CTCF.PROM+mean_CTCF.ENH+mean_DNAse.PROM+mean_DNAse.ENH+mean_H3K27ac.PROM+mean_H3K27ac.ENH+mean_p300.PROM+mean_p300.ENH+mean_H3K4me3.PROM+  
#       mean_H3K4me3.ENH+mean_H3K4me1.PROM+mean_H3K4me1.ENH+ScoreFire.PROM+ScoreFire.ENH+      
#       INS.PROM+INS.ENH+DI.PROM+DI.ENH+mean_H3K27me3.PROM+ 
#       mean_H3K27me3.ENH+essentielGenes.YESorNO,data=asso.complexity.ABC.FULL.imput,id=membership,family = gaussian, corstr = "exchangeable")
# 
# summary(gee.model.imput)
# p.value.gee = 2*pnorm(abs(gee.model.imput$coefficients/sqrt(diag(gee.model.imput$robust.variance))), lower.tail=F)
# 
# gee.model.imput.withTADs = gee(log(complexity)~ mean_CTCF.PROM+mean_CTCF.ENH+mean_DNAse.PROM+mean_DNAse.ENH+mean_H3K27ac.PROM+mean_H3K27ac.ENH+mean_p300.PROM+mean_p300.ENH+mean_H3K4me3.PROM+  
#                         mean_H3K4me3.ENH+mean_H3K4me1.PROM+mean_H3K4me1.ENH+ScoreFire.PROM+ScoreFire.ENH+      
#                         INS.PROM+INS.ENH+DI.PROM+DI.ENH+mean_H3K27me3.PROM+ 
#                         mean_H3K27me3.ENH+essentielGenes.YESorNO,data=asso.complexity.ABC.FULL.withTADs,id=TADs,family = gaussian, corstr = "exchangeable")
# 
# summary(gee.model.imput.withTADs)
# p.value.gee.TADs = 2*pnorm(abs(gee.model.imput.withTADs$coefficients/sqrt(diag(gee.model.imput.withTADs$robust.variance))), lower.tail=F)


#Can we modeleize the statut of genes: essential versus non-essential
boxplot(Expression~essentialGenes,unique.EssentialGenes.CRN.ABC)
summary(aov(Expression~essentialGenes,unique.EssentialGenes.CRN.ABC))
t.test(unique.EssentialGenes.CRN.ABC[unique.EssentialGenes.CRN.ABC$essentialGenes==1,"Expression"], unique.EssentialGenes.CRN.ABC[unique.EssentialGenes.CRN.ABC$essentialGenes==0,"Expression"], alternative = "greater")

logit.model.genes.status = glm(essentialGenes~Expression+DI.ENH+ScoreFire.ENH+INS.ENH+mean_CTCF.ENH+mean_H3K27ac.ENH+mean_DNAse.ENH+mean_H3K4me1.ENH+mean_H3K4me3.ENH+mean_p300.ENH+mean_H3K27me3.ENH+
                    DI.PROM+ScoreFire.PROM+INS.PROM+mean_CTCF.PROM+mean_H3K27ac.PROM+mean_DNAse.PROM+mean_H3K4me1.PROM+mean_H3K4me3.PROM+mean_p300.PROM+mean_H3K27me3.PROM,regul.promoters.ABC.imput , family = binomial(link="logit"))
summary(logit.model.genes.status)

null.model.genes.status = glm(essentialGenes~1, regul.promoters.ABC.imput , family = binomial(link="logit"))
summary(null.model.genes.status)

1-deviance(logit.model.genes.status)/deviance(null.model.genes.status)

#SCZ SNPs analysis 
SCZ3.all = read.table("../Genetic/SCZ3.tsv.gz", header=T, fill = T)
SCZ3.all = na.omit(SCZ3.all)

SCZ3.all.clumped = read.table("../Genetic/PGC3_SCZ_wave3_public.clumped.v2.tsv", header=T, fill = T)
GRanges.snps.SCZ3.clumped = GRanges(seqnames=paste0("chr",SCZ3.all.clumped$CHR), ranges=IRanges(start=SCZ3.all.clumped$BP, end=SCZ3.all.clumped$BP, names=SCZ3.all.clumped$SNP), pval=SCZ3.all.clumped$P)
commom.variants = read.table("allchrs_filter.bim", header=F)
commom.variants$V1 = paste0("chr",commom.variants$V1)

GRanges.common = GRanges(seqnames=commom.variants$V1, ranges=IRanges(start=commom.variants$V4, end=commom.variants$V4, names = commom.variants$V2))

background.beh.overlap = sum(countOverlaps(GRanges.common, unique.Regul.ABC))/length(GRanges.common)


threshold = c(0.10,0.05,0.025,0.01,0.001,0.0001,0.00001,0.000001, 0.0000001)

subset.SNPs.by.signi = lapply(threshold, function(x) GRanges.snps.SCZ3.clumped[GRanges.snps.SCZ3.clumped$pval<=x])
subset.SNPs.by.nsigni = lapply(threshold, function(x) GRanges.snps.SCZ3.clumped[GRanges.snps.SCZ3.clumped$pval>x])

overlap.beh.ABC = sapply(1:length(threshold), function(x) sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.signi[[x]])) / (sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.signi[[x]]))+sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.nsigni[[x]]))))
overlap.beh.Rao =sapply(1:length(threshold), function(x) sum(countOverlaps(unique.Regul.Rao, subset.SNPs.by.signi[[x]])) / (sum(countOverlaps(unique.Regul.Rao, subset.SNPs.by.signi[[x]]))+sum(countOverlaps(unique.Regul.Rao, subset.SNPs.by.nsigni[[x]]))))
overlap.beh.DNAse =sapply(1:length(threshold), function(x) sum(countOverlaps(unique.Regul.DNAse, subset.SNPs.by.signi[[x]])) / (sum(countOverlaps(unique.Regul.DNAse, subset.SNPs.by.signi[[x]]))+sum(countOverlaps(unique.Regul.DNAse, subset.SNPs.by.nsigni[[x]]))))

htr = rbind(matrix(overlap.beh.ABC, ncol=1),matrix(overlap.beh.Rao, ncol=1),matrix(overlap.beh.DNAse, ncol=1))
htr = cbind(htr,matrix(as.vector(sapply(c("ABC", "Rao", "DNAse"), function (x) rep(x,9))),ncol=1))
htr = data.frame(cbind(matrix(rep(threshold, 3),ncol=1), htr))

colnames(htr) = c("Pvalue","Prop", "Method")

htr$Pvalue = as.numeric(as.character(htr$Pvalue))
htr$Prop = as.numeric(as.character(htr$Prop))
ggplot(htr, aes(x=Pvalue, y=Prop, color=Method))+geom_line()+geom_hline(yintercept = background.beh.overlap,linetype = "dashed")+ylab("SNP overlapping proportion")+xlab("Pvalue")

#Enrichment between regulatory-elements and accessible regions + rest of genome
#Deletion of promoter or gene regions 

SNPs.enrichment.ABC = lapply(1:length(threshold), function(x)
  {
  signi.Regul.ABC = sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.signi[[x]]))
  nsigni.Regul.ABC = sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.nsigni[[x]]))
  
  signi.Regul.candidates = sum(countOverlaps(regul.candidates.ABC, subset.SNPs.by.signi[[x]]))
  nsigni.Regul.candidates = sum(countOverlaps(regul.candidates.ABC, subset.SNPs.by.nsigni[[x]]))
  
  
  signi.Prom.ABC = sum(countOverlaps(unique.Promoters.ABC, subset.SNPs.by.signi[[x]]))
  nsigni.Prom.ABC = sum(countOverlaps(unique.Promoters.ABC, subset.SNPs.by.nsigni[[x]]))
  
  signi.Prom.candidates = sum(countOverlaps(promoters.candidates.ABC, subset.SNPs.by.signi[[x]]))
  nsigni.Prom.candidates = sum(countOverlaps(promoters.candidates.ABC, subset.SNPs.by.nsigni[[x]]))
  
  e.regul = fisher.test(matrix(c(signi.Regul.ABC, signi.Regul.candidates, nsigni.Regul.ABC,nsigni.Regul.candidates), ncol=2), alternative="two.sided")$estimate
  ci.regul = fisher.test(matrix(c(signi.Regul.ABC, signi.Regul.candidates, nsigni.Regul.ABC,nsigni.Regul.candidates), ncol=2), alternative="two.sided")$conf.int
  e.prom =  fisher.test(matrix(c(signi.Prom.ABC, signi.Prom.candidates, nsigni.Prom.ABC,nsigni.Prom.candidates), ncol=2), alternative="two.sided")$estimate
  ci.prom =  fisher.test(matrix(c(signi.Prom.ABC, signi.Prom.candidates, nsigni.Prom.ABC,nsigni.Prom.candidates), ncol=2), alternative="two.sided")$conf.int
  
  
  signi.CRN = sum(countOverlaps(unique.Promoters.ABC, subset.SNPs.by.signi[[x]]))+sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.signi[[x]]))
  nsigni.CRN = sum(countOverlaps(unique.Promoters.ABC, subset.SNPs.by.nsigni[[x]]))+sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.signi[[x]]))
  
  
  all.signi.genome = length(subset.SNPs.by.signi[[x]])-sum(countOverlaps(regul.candidates.ABC, subset.SNPs.by.signi[[x]]))-sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.signi[[x]]))-
    sum(countOverlaps(unique.Promoters.ABC,subset.SNPs.by.signi[[x]]))-sum(countOverlaps(promoters.candidates.ABC,subset.SNPs.by.signi[[x]]))
  all.nsigni.genome = length(subset.SNPs.by.nsigni[[x]])-sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.nsigni[[x]]))-sum(countOverlaps(regul.candidates.ABC, subset.SNPs.by.nsigni[[x]]))-
    sum(countOverlaps(unique.Promoters.ABC,subset.SNPs.by.nsigni[[x]]))-sum(countOverlaps(promoters.candidates.ABC,subset.SNPs.by.nsigni[[x]]))
  
  e.CRN = fisher.test(matrix(c(signi.CRN,all.signi.genome,nsigni.CRN,all.nsigni.genome), ncol=2), alternative="two.sided")$estimate
  ci.CRN = fisher.test(matrix(c(signi.CRN,all.signi.genome,nsigni.CRN,all.nsigni.genome), ncol=2), alternative="two.sided")$conf.int
  
  
  return(list("regul" = c(e.regul, ci.regul), "prom"= c(e.prom, ci.prom), "CRN" = c(e.CRN,ci.CRN )))
})

table.enrichment.ABC = lapply(1:length(threshold), function(x)
{
  signi.Regul.ABC = sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.signi[[x]]))
  nsigni.Regul.ABC = sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.nsigni[[x]]))
  
  signi.Regul.candidates = sum(countOverlaps(regul.candidates.ABC, subset.SNPs.by.signi[[x]]))
  nsigni.Regul.candidates = sum(countOverlaps(regul.candidates.ABC, subset.SNPs.by.nsigni[[x]]))
  
  
  signi.Prom.ABC = sum(countOverlaps(unique.Promoters.ABC, subset.SNPs.by.signi[[x]]))
  nsigni.Prom.ABC = sum(countOverlaps(unique.Promoters.ABC, subset.SNPs.by.nsigni[[x]]))
  
  signi.Prom.candidates = sum(countOverlaps(promoters.candidates.ABC, subset.SNPs.by.signi[[x]]))
  nsigni.Prom.candidates = sum(countOverlaps(promoters.candidates.ABC, subset.SNPs.by.nsigni[[x]]))
  
  e.regul = matrix(c(signi.Regul.ABC, signi.Regul.candidates, nsigni.Regul.ABC,nsigni.Regul.candidates), ncol=2)
  e.prom =  matrix(c(signi.Prom.ABC, signi.Prom.candidates, nsigni.Prom.ABC,nsigni.Prom.candidates), ncol=2)
  
  
  signi.CRN = sum(countOverlaps(unique.Promoters.ABC, subset.SNPs.by.signi[[x]]))+sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.signi[[x]]))
  nsigni.CRN = sum(countOverlaps(unique.Promoters.ABC, subset.SNPs.by.nsigni[[x]]))+sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.signi[[x]]))
  
  
  all.signi.genome = length(subset.SNPs.by.signi[[x]])-sum(countOverlaps(regul.candidates.ABC, subset.SNPs.by.signi[[x]]))-sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.signi[[x]]))-
    sum(countOverlaps(unique.Promoters.ABC,subset.SNPs.by.signi[[x]]))-sum(countOverlaps(promoters.candidates.ABC,subset.SNPs.by.signi[[x]]))
  all.nsigni.genome = length(subset.SNPs.by.nsigni[[x]])-sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.nsigni[[x]]))-sum(countOverlaps(regul.candidates.ABC, subset.SNPs.by.nsigni[[x]]))-
    sum(countOverlaps(unique.Promoters.ABC,subset.SNPs.by.nsigni[[x]]))-sum(countOverlaps(promoters.candidates.ABC,subset.SNPs.by.nsigni[[x]]))
  
  e.CRN = matrix(c(signi.CRN,all.signi.genome,nsigni.CRN,all.nsigni.genome), ncol=2)
  
  
  return(list("regul" = e.regul, "prom"= e.prom, "CRN" = e.CRN))
})
saveRDS(SNPs.enrichment.ABC, "SNPs_enrichment_ABC.rds")
saveRDS(table.enrichment.ABC, "table_enrichment_ABC.rds")

SNPs.enrichment.Rao = lapply(1:length(threshold), function(x)
{
  signi.Regul.Rao = sum(countOverlaps(unique.Regul.Rao, subset.SNPs.by.signi[[x]]))
  nsigni.Regul.Rao = sum(countOverlaps(unique.Regul.Rao, subset.SNPs.by.nsigni[[x]]))
  
  signi.Regul.candidates = sum(countOverlaps(regul.candidates.Rao, subset.SNPs.by.signi[[x]]))
  nsigni.Regul.candidates = sum(countOverlaps(regul.candidates.Rao, subset.SNPs.by.nsigni[[x]]))
  
  
  signi.Prom.Rao = sum(countOverlaps(unique.Promoters.Rao, subset.SNPs.by.signi[[x]]))
  nsigni.Prom.Rao = sum(countOverlaps(unique.Promoters.Rao, subset.SNPs.by.nsigni[[x]]))
  
  signi.Prom.candidates = sum(countOverlaps(promoters.candidates.Rao, subset.SNPs.by.signi[[x]]))
  nsigni.Prom.candidates = sum(countOverlaps(promoters.candidates.Rao, subset.SNPs.by.nsigni[[x]]))
  
  e.regul = fisher.test(matrix(c(signi.Regul.Rao, signi.Regul.candidates, nsigni.Regul.Rao,nsigni.Regul.candidates), ncol=2), alternative="two.sided")$estimate
  ci.regul = fisher.test(matrix(c(signi.Regul.Rao, signi.Regul.candidates, nsigni.Regul.Rao,nsigni.Regul.candidates), ncol=2), alternative="two.sided")$conf.int
  e.prom =  fisher.test(matrix(c(signi.Prom.Rao, signi.Prom.candidates, nsigni.Prom.Rao,nsigni.Prom.candidates), ncol=2), alternative="two.sided")$estimate
  ci.prom =  fisher.test(matrix(c(signi.Prom.Rao, signi.Prom.candidates, nsigni.Prom.Rao,nsigni.Prom.candidates), ncol=2), alternative="two.sided")$conf.int
  
  
  signi.CRN = sum(countOverlaps(unique.Promoters.Rao, subset.SNPs.by.signi[[x]]))+sum(countOverlaps(unique.Regul.Rao, subset.SNPs.by.signi[[x]]))
  nsigni.CRN = sum(countOverlaps(unique.Promoters.Rao, subset.SNPs.by.nsigni[[x]]))+sum(countOverlaps(unique.Regul.Rao, subset.SNPs.by.signi[[x]]))
  
  
  all.signi.genome = length(subset.SNPs.by.signi[[x]])-sum(countOverlaps(regul.candidates.Rao, subset.SNPs.by.signi[[x]]))-sum(countOverlaps(unique.Regul.Rao, subset.SNPs.by.signi[[x]]))-
    sum(countOverlaps(unique.Promoters.Rao,subset.SNPs.by.signi[[x]]))-sum(countOverlaps(promoters.candidates.Rao,subset.SNPs.by.signi[[x]]))
  all.nsigni.genome = length(subset.SNPs.by.nsigni[[x]])-sum(countOverlaps(unique.Regul.Rao, subset.SNPs.by.nsigni[[x]]))-sum(countOverlaps(regul.candidates.Rao, subset.SNPs.by.nsigni[[x]]))-
    sum(countOverlaps(unique.Promoters.Rao,subset.SNPs.by.nsigni[[x]]))-sum(countOverlaps(promoters.candidates.Rao,subset.SNPs.by.nsigni[[x]]))
  
  e.CRN = fisher.test(matrix(c(signi.CRN,all.signi.genome,nsigni.CRN,all.nsigni.genome), ncol=2), alternative="two.sided")$estimate
  ci.CRN = fisher.test(matrix(c(signi.CRN,all.signi.genome,nsigni.CRN,all.nsigni.genome), ncol=2), alternative="two.sided")$conf.int
  
  
  return(list("regul" = c(e.regul, ci.regul), "prom"= c(e.prom, ci.prom), "CRN" = c(e.CRN,ci.CRN )))
})


table.enrichment.Rao = lapply(1:length(threshold), function(x)
{
  signi.Regul.Rao = sum(countOverlaps(unique.Regul.Rao, subset.SNPs.by.signi[[x]]))
  nsigni.Regul.Rao = sum(countOverlaps(unique.Regul.Rao, subset.SNPs.by.nsigni[[x]]))
  
  signi.Regul.candidates = sum(countOverlaps(regul.candidates.Rao, subset.SNPs.by.signi[[x]]))
  nsigni.Regul.candidates = sum(countOverlaps(regul.candidates.Rao, subset.SNPs.by.nsigni[[x]]))
  
  
  signi.Prom.Rao = sum(countOverlaps(unique.Promoters.Rao, subset.SNPs.by.signi[[x]]))
  nsigni.Prom.Rao = sum(countOverlaps(unique.Promoters.Rao, subset.SNPs.by.nsigni[[x]]))
  
  signi.Prom.candidates = sum(countOverlaps(promoters.candidates.Rao, subset.SNPs.by.signi[[x]]))
  nsigni.Prom.candidates = sum(countOverlaps(promoters.candidates.Rao, subset.SNPs.by.nsigni[[x]]))
  
  e.regul = matrix(c(signi.Regul.Rao, signi.Regul.candidates, nsigni.Regul.Rao,nsigni.Regul.candidates), ncol=2)
  e.prom =  matrix(c(signi.Prom.Rao, signi.Prom.candidates, nsigni.Prom.Rao,nsigni.Prom.candidates), ncol=2)
  
  
  signi.CRN = sum(countOverlaps(unique.Promoters.Rao, subset.SNPs.by.signi[[x]]))+sum(countOverlaps(unique.Regul.Rao, subset.SNPs.by.signi[[x]]))
  nsigni.CRN = sum(countOverlaps(unique.Promoters.Rao, subset.SNPs.by.nsigni[[x]]))+sum(countOverlaps(unique.Regul.Rao, subset.SNPs.by.signi[[x]]))
  
  
  all.signi.genome = length(subset.SNPs.by.signi[[x]])-sum(countOverlaps(regul.candidates.Rao, subset.SNPs.by.signi[[x]]))-sum(countOverlaps(unique.Regul.Rao, subset.SNPs.by.signi[[x]]))-
    sum(countOverlaps(unique.Promoters.Rao,subset.SNPs.by.signi[[x]]))-sum(countOverlaps(promoters.candidates.Rao,subset.SNPs.by.signi[[x]]))
  all.nsigni.genome = length(subset.SNPs.by.nsigni[[x]])-sum(countOverlaps(unique.Regul.Rao, subset.SNPs.by.nsigni[[x]]))-sum(countOverlaps(regul.candidates.Rao, subset.SNPs.by.nsigni[[x]]))-
    sum(countOverlaps(unique.Promoters.Rao,subset.SNPs.by.nsigni[[x]]))-sum(countOverlaps(promoters.candidates.Rao,subset.SNPs.by.nsigni[[x]]))
  
  e.CRN = matrix(c(signi.CRN,all.signi.genome,nsigni.CRN,all.nsigni.genome), ncol=2)
  
  
  return(list("regul" = e.regul, "prom"= e.prom, "CRN" = e.CRN))
})


saveRDS(SNPs.enrichment.Rao, "SNPs_enrichment_Rao.rds")
saveRDS(table.enrichment.Rao, "table_enrichment_Rao.rds")


SNPs.enrichment.DNAse = lapply(1:length(threshold), function(x)
{
  signi.Regul.DNAse = sum(countOverlaps(unique.Regul.DNAse, subset.SNPs.by.signi[[x]]))
  nsigni.Regul.DNAse = sum(countOverlaps(unique.Regul.DNAse, subset.SNPs.by.nsigni[[x]]))
  
  signi.Regul.candidates = sum(countOverlaps(regul.candidates.DNAse, subset.SNPs.by.signi[[x]]))
  nsigni.Regul.candidates = sum(countOverlaps(regul.candidates.DNAse, subset.SNPs.by.nsigni[[x]]))
  
  
  signi.Prom.DNAse = sum(countOverlaps(unique.Promoters.DNAse, subset.SNPs.by.signi[[x]]))
  nsigni.Prom.DNAse = sum(countOverlaps(unique.Promoters.DNAse, subset.SNPs.by.nsigni[[x]]))
  
  signi.Prom.candidates = sum(countOverlaps(promoters.candidates.DNAse, subset.SNPs.by.signi[[x]]))
  nsigni.Prom.candidates = sum(countOverlaps(promoters.candidates.DNAse, subset.SNPs.by.nsigni[[x]]))
  
  e.regul = fisher.test(matrix(c(signi.Regul.DNAse, signi.Regul.candidates, nsigni.Regul.DNAse,nsigni.Regul.candidates), ncol=2), alternative="two.sided")$estimate
  ci.regul = fisher.test(matrix(c(signi.Regul.DNAse, signi.Regul.candidates, nsigni.Regul.DNAse,nsigni.Regul.candidates), ncol=2), alternative="two.sided")$conf.int
  e.prom =  fisher.test(matrix(c(signi.Prom.DNAse, signi.Prom.candidates, nsigni.Prom.DNAse,nsigni.Prom.candidates), ncol=2), alternative="two.sided")$estimate
  ci.prom =  fisher.test(matrix(c(signi.Prom.DNAse, signi.Prom.candidates, nsigni.Prom.DNAse,nsigni.Prom.candidates), ncol=2), alternative="two.sided")$conf.int
  
  
  signi.CRN = sum(countOverlaps(unique.Promoters.DNAse, subset.SNPs.by.signi[[x]]))+sum(countOverlaps(unique.Regul.DNAse, subset.SNPs.by.signi[[x]]))
  nsigni.CRN = sum(countOverlaps(unique.Promoters.DNAse, subset.SNPs.by.nsigni[[x]]))+sum(countOverlaps(unique.Regul.DNAse, subset.SNPs.by.signi[[x]]))
  
  
  all.signi.genome = length(subset.SNPs.by.signi[[x]])-sum(countOverlaps(regul.candidates.DNAse, subset.SNPs.by.signi[[x]]))-sum(countOverlaps(unique.Regul.DNAse, subset.SNPs.by.signi[[x]]))-
    sum(countOverlaps(unique.Promoters.DNAse,subset.SNPs.by.signi[[x]]))-sum(countOverlaps(promoters.candidates.DNAse,subset.SNPs.by.signi[[x]]))
  all.nsigni.genome = length(subset.SNPs.by.nsigni[[x]])-sum(countOverlaps(unique.Regul.DNAse, subset.SNPs.by.nsigni[[x]]))-sum(countOverlaps(regul.candidates.DNAse, subset.SNPs.by.nsigni[[x]]))-
    sum(countOverlaps(unique.Promoters.DNAse,subset.SNPs.by.nsigni[[x]]))-sum(countOverlaps(promoters.candidates.DNAse,subset.SNPs.by.nsigni[[x]]))
  
  e.CRN = fisher.test(matrix(c(signi.CRN,all.signi.genome,nsigni.CRN,all.nsigni.genome), ncol=2), alternative="two.sided")$estimate
  ci.CRN = fisher.test(matrix(c(signi.CRN,all.signi.genome,nsigni.CRN,all.nsigni.genome), ncol=2), alternative="two.sided")$conf.int
  
  
  return(list("regul" = c(e.regul, ci.regul), "prom"= c(e.prom, ci.prom), "CRN" = c(e.CRN,ci.CRN )))
})

table.enrichment.DNAse = lapply(1:length(threshold), function(x)
{
  signi.Regul.DNAse = sum(countOverlaps(unique.Regul.DNAse, subset.SNPs.by.signi[[x]]))
  nsigni.Regul.DNAse = sum(countOverlaps(unique.Regul.DNAse, subset.SNPs.by.nsigni[[x]]))
  
  signi.Regul.candidates = sum(countOverlaps(regul.candidates.DNAse, subset.SNPs.by.signi[[x]]))
  nsigni.Regul.candidates = sum(countOverlaps(regul.candidates.DNAse, subset.SNPs.by.nsigni[[x]]))
  
  
  signi.Prom.DNAse = sum(countOverlaps(unique.Promoters.DNAse, subset.SNPs.by.signi[[x]]))
  nsigni.Prom.DNAse = sum(countOverlaps(unique.Promoters.DNAse, subset.SNPs.by.nsigni[[x]]))
  
  signi.Prom.candidates = sum(countOverlaps(promoters.candidates.DNAse, subset.SNPs.by.signi[[x]]))
  nsigni.Prom.candidates = sum(countOverlaps(promoters.candidates.DNAse, subset.SNPs.by.nsigni[[x]]))
  
  e.regul = fisher.test(matrix(c(signi.Regul.DNAse, signi.Regul.candidates, nsigni.Regul.DNAse,nsigni.Regul.candidates), ncol=2), alternative="two.sided")$estimate
  
  e.prom =  matrix(c(signi.Prom.DNAse, signi.Prom.candidates, nsigni.Prom.DNAse,nsigni.Prom.candidates), ncol=2)
  
  
  signi.CRN = sum(countOverlaps(unique.Promoters.DNAse, subset.SNPs.by.signi[[x]]))+sum(countOverlaps(unique.Regul.DNAse, subset.SNPs.by.signi[[x]]))
  nsigni.CRN = sum(countOverlaps(unique.Promoters.DNAse, subset.SNPs.by.nsigni[[x]]))+sum(countOverlaps(unique.Regul.DNAse, subset.SNPs.by.signi[[x]]))
  
  
  all.signi.genome = length(subset.SNPs.by.signi[[x]])-sum(countOverlaps(regul.candidates.DNAse, subset.SNPs.by.signi[[x]]))-sum(countOverlaps(unique.Regul.DNAse, subset.SNPs.by.signi[[x]]))-
    sum(countOverlaps(unique.Promoters.DNAse,subset.SNPs.by.signi[[x]]))-sum(countOverlaps(promoters.candidates.DNAse,subset.SNPs.by.signi[[x]]))
  all.nsigni.genome = length(subset.SNPs.by.nsigni[[x]])-sum(countOverlaps(unique.Regul.DNAse, subset.SNPs.by.nsigni[[x]]))-sum(countOverlaps(regul.candidates.DNAse, subset.SNPs.by.nsigni[[x]]))-
    sum(countOverlaps(unique.Promoters.DNAse,subset.SNPs.by.nsigni[[x]]))-sum(countOverlaps(promoters.candidates.DNAse,subset.SNPs.by.nsigni[[x]]))
  
  e.CRN = matrix(c(signi.CRN,all.signi.genome,nsigni.CRN,all.nsigni.genome), ncol=2)
  ci.CRN = matrix(c(signi.CRN,all.signi.genome,nsigni.CRN,all.nsigni.genome), ncol=2)
  
  
  return(list("regul" = e.regul, "prom"= e.prom, "CRN" = e.CRN))
})

saveRDS(SNPs.enrichment.DNAse, "SNPs_enrichment_DNAse.rds")
saveRDS(table.enrichment.DNAse, "table_enrichment_DNAse.rds")


names(SNPs.enrichment.ABC) = threshold
names(SNPs.enrichment.Rao) = threshold
names(SNPs.enrichment.DNAse) = threshold

SNPs.enrichment.ABC$`1e-05`
SNPs.enrichment.Rao$`1e-05`
SNPs.enrichment.DNAse$`1e-05`


p = data.frame(do.call("rbind", SNPs.enrichment.DNAse$`0.1`))
p$method = "DNAse"
p$pvalue= 0.1

p.1 = data.frame(do.call("rbind", SNPs.enrichment.ABC$`0.1`))
p.1$method = "ABC"
p.1$pvalue= 0.1

p.2 = data.frame(do.call("rbind", SNPs.enrichment.Rao$`0.1`))
p.2$method = "Rao"
p.2$pvalue= 0.1

p = rbind(p, p.1, p.2)
for( i in 2:length(threshold)){
  q = data.frame(do.call("rbind", SNPs.enrichment.DNAse[[i]]))
  q$method = "DNAse"
  q$pvalue= threshold[i]
  
  q.1 = data.frame(do.call("rbind", SNPs.enrichment.ABC[[i]]))
  q.1$method = "ABC"
  q.1$pvalue= threshold[i]
  
  q.2 = data.frame(do.call("rbind", SNPs.enrichment.Rao[[i]]))
  q.2$method = "Rao"
  q.2$pvalue= threshold[i]
  
  p = rbind(p,q,q.1,q.2)
}
colnames(p) = c("OR", "CI_lower", "CI_upper", "method", "pvalue")

p$element = ifelse(grepl("regul",rownames(p)), "Regul", ifelse(grepl("CRN",rownames(p)), "CRN", "prom"))

p$pvalue = factor(as.character(p$pvalue),levels = c("0.1","0.05", "0.025", "0.01","0.001", "1e-04", "1e-05", "1e-06", "1e-07"))
ggplot(p[p$element=="CRN", ], aes(x=OR, y=pvalue)) + geom_pointrange(aes(xmin=CI_lower, xmax=CI_upper, color=method),position=position_dodge(width=0.5)) +  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed")
ggplot(p[p$element=="Regul", ], aes(x=OR, y=pvalue)) + geom_pointrange(aes(xmin=CI_lower, xmax=CI_upper, color=method),position=position_dodge(width=0.5)) +  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed")
ggplot(p[p$element=="prom", ], aes(x=OR, y=pvalue)) + geom_pointrange(aes(xmin=CI_lower, xmax=CI_upper, color=method),position=position_dodge(width=0.5)) +  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed")

#Enrichment between regulatory in large vs small CRNs
small.Regul.ABC = read.table("small_ENH_ABC.bed", header=T,col.names = c("seqnames", "start", "end")) 
GR.small.Regul.ABC = makeGRangesFromDataFrame(small.Regul.ABC)
medium.Regul.ABC = read.table("medium_ENH_ABC.bed", header=T,col.names = c("seqnames", "start", "end"))
GR.medium.Regul.ABC = makeGRangesFromDataFrame(medium.Regul.ABC)
large.Regul.ABC = read.table("Large_ENH_ABC.bed", header=T,col.names = c("seqnames", "start", "end"))
GR.large.Regul.ABC =makeGRangesFromDataFrame(large.Regul.ABC)

small.Regul.Rao = read.table("small_ENH_Rao.bed", header=T,col.names = c("seqnames", "start", "end")) 
GR.small.Regul.Rao = makeGRangesFromDataFrame(small.Regul.Rao)

large.Regul.Rao = read.table("Large_ENH_Rao.bed", header=T,col.names = c("seqnames", "start", "end"))
GR.large.Regul.Rao =makeGRangesFromDataFrame(large.Regul.Rao)


small.Regul.DNAse = read.table("small_ENH_DNAse.bed", header=T,col.names = c("seqnames", "start", "end")) 
GR.small.Regul.DNAse = makeGRangesFromDataFrame(small.Regul.DNAse)

large.Regul.DNAse = read.table("Large_ENH_DNAse.bed", header=T,col.names = c("seqnames", "start", "end"))
GR.large.Regul.DNAse =makeGRangesFromDataFrame(large.Regul.DNAse)


library(epitools)

for(i in 1:length(threshold)){
print(epitab(matrix(c(sum(countOverlaps(GR.small.Regul.ABC,subset.SNPs.by.signi[[i]])),sum(countOverlaps(GR.small.Regul.ABC,subset.SNPs.by.nsigni[[i]])),sum(countOverlaps(GR.medium.Regul.ABC,subset.SNPs.by.signi[[i]])),sum(countOverlaps(GR.medium.Regul.ABC,subset.SNPs.by.nsigni[[i]])),sum(countOverlaps(GR.large.Regul.ABC,subset.SNPs.by.signi[[i]])),sum(countOverlaps(GR.large.Regul.ABC,subset.SNPs.by.nsigni[[i]]))), ncol=2,byrow = T), method="riskratio", riskratio = "small"))
}

for(i in 1:length(threshold)){
  print(epitab(matrix(c(sum(countOverlaps(GR.small.Regul.Rao,subset.SNPs.by.signi[[i]])),sum(countOverlaps(GR.small.Regul.Rao,subset.SNPs.by.nsigni[[i]])),sum(countOverlaps(GR.large.Regul.Rao,subset.SNPs.by.signi[[i]])),sum(countOverlaps(GR.large.Regul.Rao,subset.SNPs.by.nsigni[[i]]))), ncol=2,byrow = T), method="riskratio"))
}


for(i in 1:length(threshold)){
  print(epitab(matrix(c(sum(countOverlaps(GR.small.Regul.DNAse,subset.SNPs.by.signi[[i]])),sum(countOverlaps(GR.small.Regul.DNAse,subset.SNPs.by.nsigni[[i]])),sum(countOverlaps(GR.large.Regul.DNAse,subset.SNPs.by.signi[[i]])),sum(countOverlaps(GR.large.Regul.DNAse,subset.SNPs.by.nsigni[[i]]))), ncol=2,byrow = T), method="riskratio"))
}


sapply(subset.SNPs.by.signi, length)
sapply(subset.SNPs.by.nsigni, length)

subset.SNPs.by.signi[[3]]
prop =c()
for(i in 1:1000){
  j = sample(1:length(regul.candidates.ABC), 29089, replace=F)
  
  p = sum(countOverlaps(regul.candidates.ABC[j], subset.SNPs.by.signi[[3]])) / (sum(countOverlaps(regul.candidates.ABC[j], subset.SNPs.by.signi[[3]])) + sum(countOverlaps(regul.candidates.ABC[j], subset.SNPs.by.nsigni[[3]])))
  prop = c(prop,p)
}

sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.signi[[3]]))/(sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.signi[[3]]))+sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.nsigni[[3]])))



sum(countOverlaps(regul.candidates.ABC, subset.SNPs.by.signi[[3]]))
sum(countOverlaps(regul.candidates.ABC, subset.SNPs.by.nsigni[[3]]))

a = sapply(subset.SNPs.by.signi, function(x) sum(countOverlaps(unique.Regul.ABC, x)))+sapply(subset.SNPs.by.signi, function(x) sum(countOverlaps(unique.Promoters.ABC, x)))
b = sapply(subset.SNPs.by.nsigni, function(x) sum(countOverlaps(unique.Regul.ABC, x)))+sapply(subset.SNPs.by.nsigni, function(x) sum(countOverlaps(unique.Promoters.ABC, x)))

c = sapply(subset.SNPs.by.signi, length)-sapply(subset.SNPs.by.signi, function(x) sum(countOverlaps(unique.Regul.ABC, x)))-sapply(subset.SNPs.by.signi, function(x) sum(countOverlaps(unique.Promoters.ABC, x)))-
  sapply(subset.SNPs.by.signi, function(x) sum(countOverlaps(regul.candidates.ABC, x)))-sapply(subset.SNPs.by.signi, function(x) sum(countOverlaps(promoters.candidates.ABC, x)))
d = sapply(subset.SNPs.by.nsigni, length)-sapply(subset.SNPs.by.nsigni, function(x) sum(countOverlaps(unique.Regul.ABC, x)))-sapply(subset.SNPs.by.nsigni, function(x) sum(countOverlaps(unique.Promoters.ABC, x)))-
  sapply(subset.SNPs.by.nsigni, function(x) sum(countOverlaps(regul.candidates.ABC, x)))-sapply(subset.SNPs.by.nsigni, function(x) sum(countOverlaps(promoters.candidates.ABC, x)))

sapply(1:length(a), function(x) fisher.test(matrix(c(a[x],b[x],c[x],d[x]),ncol=2), alternative = "greater")$estimate)


a = sapply(subset.SNPs.by.signi, function(x) sum(countOverlaps(unique.Regul.Rao, x)))+sapply(subset.SNPs.by.signi, function(x) sum(countOverlaps(unique.Promoters.Rao, x)))
b = sapply(subset.SNPs.by.nsigni, function(x) sum(countOverlaps(unique.Regul.Rao, x)))+sapply(subset.SNPs.by.nsigni, function(x) sum(countOverlaps(unique.Promoters.Rao, x)))

c = sapply(subset.SNPs.by.signi, length)-sapply(subset.SNPs.by.signi, function(x) sum(countOverlaps(unique.Regul.Rao, x)))-sapply(subset.SNPs.by.signi, function(x) sum(countOverlaps(unique.Promoters.Rao, x)))-
  sapply(subset.SNPs.by.signi, function(x) sum(countOverlaps(regul.candidates.Rao, x)))-sapply(subset.SNPs.by.signi, function(x) sum(countOverlaps(promoters.candidates.Rao, x)))
d = sapply(subset.SNPs.by.nsigni, length)-sapply(subset.SNPs.by.nsigni, function(x) sum(countOverlaps(unique.Regul.Rao, x)))-sapply(subset.SNPs.by.nsigni, function(x) sum(countOverlaps(unique.Promoters.Rao, x)))-
  sapply(subset.SNPs.by.nsigni, function(x) sum(countOverlaps(regul.candidates.Rao, x)))-sapply(subset.SNPs.by.nsigni, function(x) sum(countOverlaps(promoters.candidates.Rao, x)))

sapply(1:length(a), function(x) fisher.test(matrix(c(a[x],b[x],c[x],d[x]),ncol=2), alternative = "greater")$estimate)


a = sapply(subset.SNPs.by.signi, function(x) sum(countOverlaps(unique.Regul.DNAse, x)))+sapply(subset.SNPs.by.signi, function(x) sum(countOverlaps(unique.Promoters.DNAse, x)))
b = sapply(subset.SNPs.by.nsigni, function(x) sum(countOverlaps(unique.Regul.DNAse, x)))+sapply(subset.SNPs.by.nsigni, function(x) sum(countOverlaps(unique.Promoters.DNAse, x)))

c = sapply(subset.SNPs.by.signi, length)-sapply(subset.SNPs.by.signi, function(x) sum(countOverlaps(unique.Regul.DNAse, x)))-sapply(subset.SNPs.by.signi, function(x) sum(countOverlaps(unique.Promoters.DNAse, x)))-
  sapply(subset.SNPs.by.signi, function(x) sum(countOverlaps(regul.candidates.DNAse, x)))-sapply(subset.SNPs.by.signi, function(x) sum(countOverlaps(promoters.candidates.DNAse, x)))
d = sapply(subset.SNPs.by.nsigni, length)-sapply(subset.SNPs.by.nsigni, function(x) sum(countOverlaps(unique.Regul.DNAse, x)))-sapply(subset.SNPs.by.nsigni, function(x) sum(countOverlaps(unique.Promoters.DNAse, x)))-
  sapply(subset.SNPs.by.nsigni, function(x) sum(countOverlaps(regul.candidates.DNAse, x)))-sapply(subset.SNPs.by.nsigni, function(x) sum(countOverlaps(promoters.candidates.DNAse, x)))

sapply(1:length(a), function(x) fisher.test(matrix(c(a[x],b[x],c[x],d[x]),ncol=2), alternative = "greater")$estimate)


for(t in seq_along(threshold)){
  print(threshold[t])
  GRanges.snps.SCZ3.signi = GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=threshold[t]]
  GRanges.snps.SCZ3.nsigni = GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>threshold[t]]
  
  sum(countOverlaps(unique.Regul.ABC, GRanges.snps.SCZ3.signi))
  sum(countOverlaps(unique.Regul.ABC, GRanges.snps.SCZ3.nsigni))
  
  print(sum(countOverlaps(unique.Regul.ABC, GRanges.snps.SCZ3.signi)) / (sum(countOverlaps(unique.Regul.ABC, GRanges.snps.SCZ3.nsigni))+sum(countOverlaps(unique.Regul.ABC, GRanges.snps.SCZ3.signi))))
  print(sum(countOverlaps(unique.Regul.Rao, GRanges.snps.SCZ3.signi)) / (sum(countOverlaps(unique.Regul.Rao, GRanges.snps.SCZ3.nsigni))+sum(countOverlaps(unique.Regul.Rao, GRanges.snps.SCZ3.signi))))
  print(sum(countOverlaps(unique.Regul.DNAse, GRanges.snps.SCZ3.signi)) / (sum(countOverlaps(unique.Regul.DNAse, GRanges.snps.SCZ3.nsigni))+sum(countOverlaps(unique.Regul.DNAse, GRanges.snps.SCZ3.signi))))
  
}

for(t in seq_along(threshold)){
  print(threshold[t])
  GRanges.snps.SCZ3.signi = GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=threshold[t]]
  GRanges.snps.SCZ3.nsigni = GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>threshold[t]]
  
  print(sum(countOverlaps(regul.candidates.ABC, GRanges.snps.SCZ3.signi)) / (sum(countOverlaps(regul.candidates.ABC, GRanges.snps.SCZ3.nsigni))+sum(countOverlaps(regul.candidates.ABC, GRanges.snps.SCZ3.signi))))
  print(sum(countOverlaps(regul.candidates.Rao, GRanges.snps.SCZ3.signi)) / (sum(countOverlaps(regul.candidates.Rao, GRanges.snps.SCZ3.nsigni))+sum(countOverlaps(regul.candidates.Rao, GRanges.snps.SCZ3.signi))))
  print(sum(countOverlaps(regul.candidates.DNAse, GRanges.snps.SCZ3.signi)) / (sum(countOverlaps(regul.candidates.DNAse, GRanges.snps.SCZ3.nsigni))+sum(countOverlaps(regul.candidates.DNAse, GRanges.snps.SCZ3.signi))))
  
}


forward = GRanges(seqnames=hic.loops$chr1, ranges=IRanges(start=hic.loops$x1,end=hic.loops$x2))
reverse = GRanges(seqnames=hic.loops$chr2, ranges=IRanges(start=hic.loops$y1,end=hic.loops$y2))

all.hic = unique(do.call(c,GRangesList(forward, reverse)))

Overlaps.Regul.ABC = findOverlaps(positive.contact, unique.Regul.ABC)
Overlaps.candidates.Regul.ABC = findOverlaps(positive.contact, regul.candidates.ABC)
Overlaps.Prom.ABC = findOverlaps(positive.contact, unique.Promoters.ABC)
Overlaps.candidates.Prom.ABC = findOverlaps(positive.contact, promoters.candidates.ABC)

subset.outside.genome.ABC = all.hic[-c(queryHits(Overlaps.Regul.ABC), queryHits(Overlaps.candidates.Regul.ABC))]

Overlaps.Regul.Rao = findOverlaps(positive.contact, unique.Regul.Rao)
Overlaps.candidates.Regul.Rao = findOverlaps(positive.contact, regul.candidates.Rao)
Overlaps.Prom.Rao = findOverlaps(positive.contact, unique.Promoters.Rao)
Overlaps.candidates.Prom.Rao = findOverlaps(positive.contact, promoters.candidates.Rao)

subset.outside.genome.Rao = all.hic[-c(queryHits(Overlaps.Regul.Rao), queryHits(Overlaps.candidates.Regul.Rao))]


Overlaps.Regul.DNAse = findOverlaps(positive.contact, unique.Regul.DNAse)
Overlaps.candidates.Regul.DNAse = findOverlaps(positive.contact, regul.candidates.DNAse)
Overlaps.Prom.DNAse = findOverlaps(positive.contact, unique.Promoters.DNAse)
Overlaps.candidates.Prom.DNAse = findOverlaps(positive.contact, promoters.candidates.DNAse)

subset.outside.genome.DNAse = all.hic[-c(queryHits(Overlaps.Regul.DNAse), queryHits(Overlaps.candidates.Regul.DNAse))]

for(t in seq_along(threshold)){
  print(threshold[t])
  GRanges.snps.SCZ3.signi = GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=threshold[t]]
  GRanges.snps.SCZ3.nsigni = GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>threshold[t]]
  
  print(sum(countOverlaps(subset.outside.genome.ABC, GRanges.snps.SCZ3.signi)) / (sum(countOverlaps(subset.outside.genome.ABC, GRanges.snps.SCZ3.nsigni))+sum(countOverlaps(subset.outside.genome.ABC, GRanges.snps.SCZ3.signi))))
  print(sum(countOverlaps(subset.outside.genome.Rao, GRanges.snps.SCZ3.signi)) / (sum(countOverlaps(subset.outside.genome.Rao, GRanges.snps.SCZ3.nsigni))+sum(countOverlaps(subset.outside.genome.Rao, GRanges.snps.SCZ3.signi))))
  print(sum(countOverlaps(subset.outside.genome.DNAse, GRanges.snps.SCZ3.signi)) / (sum(countOverlaps(subset.outside.genome.DNAse, GRanges.snps.SCZ3.nsigni))+sum(countOverlaps(subset.outside.genome.DNAse, GRanges.snps.SCZ3.signi))))
  
}

#The nearest genes with SNPs 
df.nearestPROM.ABC = data.frame(distanceToNearest(GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05],unique.Promoters.ABC))
df.nearestPROM.ABC = df.nearestPROM.ABC[df.nearestPROM.ABC$distance!=0,]
df.OverlapsRegul.ABC = data.frame(findOverlaps(GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05], unique.Regul.ABC))

df.nearestOverlaps.ABC = merge(df.nearestPROM.ABC, df.OverlapsRegul.ABC, by="queryHits")

df.nearestPROM.Rao = data.frame(distanceToNearest(GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05],unique.Promoters.Rao))
df.nearestPROM.Rao = df.nearestPROM.Rao[df.nearestPROM.Rao$distance!=0,]
df.OverlapsRegul.Rao = data.frame(findOverlaps(GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05], unique.Regul.Rao))

df.nearestOverlaps.Rao = merge(df.nearestPROM.Rao, df.OverlapsRegul.Rao, by="queryHits")

df.nearestPROM.DNAse = data.frame(distanceToNearest(GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05],unique.Promoters.DNAse))
df.nearestPROM.DNAse = df.nearestPROM.DNAse[df.nearestPROM.DNAse$distance!=0,]
df.OverlapsRegul.DNAse = data.frame(findOverlaps(GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05], unique.Regul.DNAse))

df.nearestOverlaps.DNAse = merge(df.nearestPROM.DNAse, df.OverlapsRegul.DNAse, by="queryHits")

out = lapply(unique(df.nearestOverlaps.ABC$queryHits), function(x){
  prom.id = df.nearestOverlaps.ABC[df.nearestOverlaps.ABC$queryHits==x, "subjectHits.x"]
  regul.id = df.nearestOverlaps.ABC[df.nearestOverlaps.ABC$queryHits==x, "subjectHits.y"]
  
  prom = names(unique.Promoters.ABC[prom.id])
  regul = names(unique.Regul.ABC[regul.id])
  
  
  regul%in%unique(regul.promoters.ABC[regul.promoters.ABC$TargetGene==prom, "name"])
  
})


sum(unlist(out))/length(out)
#[1] 0.9145299 --> 91% des SNPs significatifs les plus proches se situent dans un element de regulation du meme reseau que le gene

out.Rao = lapply(unique(df.nearestOverlaps.Rao$queryHits), function(x){
  prom.id = df.nearestOverlaps.Rao[df.nearestOverlaps.Rao$queryHits==x, "subjectHits.x"]
  regul.id = df.nearestOverlaps.Rao[df.nearestOverlaps.Rao$queryHits==x, "subjectHits.y"]
  
  prom = mcols(unique.Promoters.Rao)[prom.id, "geneSymbol"]
  regul = paste0(seqnames(unique.Regul.Rao[regul.id]),":", start(unique.Regul.Rao[regul.id]),":", end(unique.Regul.Rao[regul.id]))
  
  
  regul%in%unique(df.Rao[df.Rao$geneSymbol.x==prom, "name"])
  
})
sum(unlist(out.Rao))/length(out.Rao)
# [1] 0.389011 --> 38 % des SNPs significatifs les plus proches se situent dans un element de regulation du meme reseau que le gene

out.DNAse = lapply(unique(df.nearestOverlaps.DNAse$queryHits), function(x){
  prom.id = df.nearestOverlaps.DNAse[df.nearestOverlaps.DNAse$queryHits==x, "subjectHits.x"]
  regul.id = df.nearestOverlaps.DNAse[df.nearestOverlaps.DNAse$queryHits==x, "subjectHits.y"]
  
  prom = mcols(unique.Promoters.DNAse)[prom.id, "geneSymbol"]
  regul = paste0(seqnames(unique.Regul.DNAse[regul.id]),":", start(unique.Regul.DNAse[regul.id]),":", end(unique.Regul.DNAse[regul.id]))
  
  
  regul%in%unique(df.DNase[df.DNase$geneSymbol.x==prom, "name"])
  
})
sum(unlist(out.DNAse))/length(out.DNAse)
# [1] 0.2424242 --> 24% des SNPs significatifs les plus proches se situent dans un element de regulation du meme reseau que le gene

regul.promoters.ABC = add.SNPs.to.df(unique.Promoters.ABC,unique.Regul.ABC, GRanges.snps.SCZ3, regul.promoters.ABC, method="ABC",0.05)
df.Rao= add.SNPs.to.df(unique.Promoters.Rao,unique.Regul.Rao, GRanges.snps.SCZ3, df.Rao, method="Rao",0.05)
df.DNase = add.SNPs.to.df(unique.Promoters.DNAse,unique.Regul.DNAse, GRanges.snps.SCZ3, df.DNase, method="DNAse",0.05)

p = aggregate(NSigni.ENH~membership,aggregate(NSigni.ENH~membership+name, regul.promoters.ABC, unique), sum)
q = aggregate(NSigni.PROM~membership,aggregate(NSigni.PROM~membership+TargetGene, regul.promoters.ABC, unique), sum)
r = aggregate(NnSigni.ENH~membership,aggregate(NnSigni.ENH~membership+name, regul.promoters.ABC, unique), sum)
s = aggregate(NnSigni.PROM~membership,aggregate(NnSigni.PROM~membership+TargetGene, regul.promoters.ABC, unique), sum)

z = merge(p,merge(q,merge(r,s, by="membership"), by="membership"),by="membership")
z1 = z[,2]+z[,3]
z2 = z[,4]+z[,5]

asso.complexity.ABC.FULL.imput$prop.SigniSNPs = z1/(z1+z2)

nrow(asso.complexity.ABC.FULL.imput[!is.na(asso.complexity.ABC.FULL.imput$prop.SigniSNPs)&asso.complexity.ABC.FULL.imput$prop.SigniSNPs>0,])
# 353 des RRCs ont des SNPs significatifs
# 21% des RRCs

p = aggregate(NSigni.ENH~membership,aggregate(NSigni.ENH~membership+name, df.Rao, unique), sum)
q = aggregate(NSigni.PROM~membership,aggregate(NSigni.PROM~membership+geneSymbol.x, df.Rao, unique), sum)
r = aggregate(NnSigni.ENH~membership,aggregate(NnSigni.ENH~membership+name, df.Rao, unique), sum)
s = aggregate(NnSigni.PROM~membership,aggregate(NnSigni.PROM~membership+geneSymbol.x, df.Rao, unique), sum)

z = merge(p,merge(q,merge(r,s, by="membership"), by="membership"),by="membership")
z1 = z[,2]+z[,3]
z2 = z[,4]+z[,5]

asso.complexity.Rao.FULL.imput$prop.SigniSNPs = z1/(z1+z2)
nrow(asso.complexity.Rao.FULL.imput[!is.na(asso.complexity.Rao.FULL.imput$prop.SigniSNPs)&asso.complexity.Rao.FULL.imput$prop.SigniSNPs>0,])
#[1] 460
#20%

p = aggregate(NSigni.ENH~membership,aggregate(NSigni.ENH~membership+name, df.DNase, unique), sum)
q = aggregate(NSigni.PROM~membership,aggregate(NSigni.PROM~membership+geneSymbol.x, df.DNase, unique), sum)
r = aggregate(NnSigni.ENH~membership,aggregate(NnSigni.ENH~membership+name, df.DNase, unique), sum)
s = aggregate(NnSigni.PROM~membership,aggregate(NnSigni.PROM~membership+geneSymbol.x, df.DNase, unique), sum)

z = merge(p,merge(q,merge(r,s, by="membership"), by="membership"),by="membership")
z1 = z[,2]+z[,3]
z2 = z[,4]+z[,5]

asso.complexity.DNAse.FULL.imput$prop.SigniSNPs = z1/(z1+z2)
nrow(asso.complexity.DNAse.FULL.imput[!is.na(asso.complexity.DNAse.FULL.imput$prop.SigniSNPs)&asso.complexity.DNAse.FULL.imput$prop.SigniSNPs>0,])
#[1] 146
#7%

sum(countOverlaps(unique.Promoters.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=.05]))
length(which(countOverlaps(unique.Promoters.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=.05])>0))/length(unique.Promoters.ABC)
#2%

sum(countOverlaps(unique.Regul.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=.05]))
length(which(countOverlaps(unique.Regul.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=.05])>0))/length(unique.Regul.ABC)
#1%

sum(countOverlaps(unique.Promoters.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=.05]))
length(which(countOverlaps(unique.Promoters.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=.05])>0))/length(unique.Promoters.Rao)
#2%

sum(countOverlaps(unique.Regul.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=.05]))
length(which(countOverlaps(unique.Regul.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=.05])>0))/length(unique.Regul.Rao)
#10%

sum(countOverlaps(unique.Promoters.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=.05]))
length(which(countOverlaps(unique.Promoters.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=.05])>0))/length(unique.Promoters.DNAse)
#2%

sum(countOverlaps(unique.Regul.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=.05]))
length(which(countOverlaps(unique.Regul.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=.05])>0))/length(unique.Regul.DNAse)
#1%

cor.test(asso.complexity.ABC.FULL.imput[!is.na(asso.complexity.ABC.FULL.imput$prop.SigniSNPs),"prop.SigniSNPs"], asso.complexity.ABC.FULL.imput[!is.na(asso.complexity.ABC.FULL.imput$prop.SigniSNPs),"complexity"], method="spearman")
cor.test(asso.complexity.Rao.FULL.imput[!is.na(asso.complexity.Rao.FULL.imput$prop.SigniSNPs),"prop.SigniSNPs"], asso.complexity.Rao.FULL.imput[!is.na(asso.complexity.Rao.FULL.imput$prop.SigniSNPs),"complexity"], method="spearman")
cor.test(asso.complexity.DNAse.FULL.imput[!is.na(asso.complexity.DNAse.FULL.imput$prop.SigniSNPs),"prop.SigniSNPs"], asso.complexity.DNAse.FULL.imput[!is.na(asso.complexity.DNAse.FULL.imput$prop.SigniSNPs),"complexity"], method="spearman")

mean(asso.complexity.ABC.FULL.imput$prop.SigniSNPs, na.rm=T)

#Proportion of Signi SNPs for equivalent elements: 
#ABC
regul.candidates.ABC$SigniSnps = countOverlaps(regul.candidates.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=.05])
regul.candidates.ABC$NSigniSnps = countOverlaps(regul.candidates.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>.05])

regul.candidates.ABC$PropSigniSNPs = regul.candidates.ABC$SigniSnps / (regul.candidates.ABC$SigniSnps + regul.candidates.ABC$NSigniSnps)
mean(regul.candidates.ABC$PropSigniSNPs, na.rm=T)

promoters.candidates.ABC$SigniSnps = countOverlaps(promoters.candidates.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=.05])
promoters.candidates.ABC$NSigniSnps = countOverlaps(promoters.candidates.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>.05])

promoters.candidates.ABC$PropSigniSNPs = promoters.candidates.ABC$SigniSnps / (promoters.candidates.ABC$SigniSnps + promoters.candidates.ABC$NSigniSnps)
mean(promoters.candidates.ABC$PropSigniSNPs, na.rm=T)

#Rao
regul.candidates.Rao$SigniSnps = countOverlaps(regul.candidates.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=.05])
regul.candidates.Rao$NSigniSnps = countOverlaps(regul.candidates.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>.05])

regul.candidates.Rao$PropSigniSNPs = regul.candidates.Rao$SigniSnps / (regul.candidates.Rao$SigniSnps + regul.candidates.Rao$NSigniSnps)
mean(regul.candidates.Rao$PropSigniSNPs, na.rm=T)

promoters.candidates.Rao$SigniSnps = countOverlaps(promoters.candidates.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=.05])
promoters.candidates.Rao$NSigniSnps = countOverlaps(promoters.candidates.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>.05])

promoters.candidates.Rao$PropSigniSNPs = promoters.candidates.Rao$SigniSnps / (promoters.candidates.Rao$SigniSnps + promoters.candidates.Rao$NSigniSnps)
mean(promoters.candidates.Rao$PropSigniSNPs, na.rm=T)

#DNAse
regul.candidates.DNAse$SigniSnps = countOverlaps(regul.candidates.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=.05])
regul.candidates.DNAse$NSigniSnps = countOverlaps(regul.candidates.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>.05])

regul.candidates.DNAse$PropSigniSNPs = regul.candidates.DNAse$SigniSnps / (regul.candidates.DNAse$SigniSnps + regul.candidates.DNAse$NSigniSnps)
mean(regul.candidates.DNAse$PropSigniSNPs, na.rm=T)

promoters.candidates.DNAse$SigniSnps = countOverlaps(promoters.candidates.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=.05])
promoters.candidates.DNAse$NSigniSnps = countOverlaps(promoters.candidates.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>.05])

promoters.candidates.DNAse$PropSigniSNPs = promoters.candidates.DNAse$SigniSnps / (promoters.candidates.DNAse$SigniSnps + promoters.candidates.DNAse$NSigniSnps)
mean(promoters.candidates.DNAse$PropSigniSNPs, na.rm=T)


forward = GRanges(seqnames=hic.loops$chr1, ranges=IRanges(start=hic.loops$x1,end=hic.loops$x2))
reverse = GRanges(seqnames=hic.loops$chr2, ranges=IRanges(start=hic.loops$y1,end=hic.loops$y2))

all.hic = unique(do.call(c,GRangesList(forward, reverse)))

Overlaps.Regul.ABC = findOverlaps(all.hic, unique.Regul.ABC)
Overlaps.candidates.Regul.ABC = findOverlaps(all.hic, regul.candidates.ABC)
Overlaps.Prom.ABC = findOverlaps(all.hic, unique.Promoters.ABC)
Overlaps.candidates.Prom.ABC = findOverlaps(all.hic, promoters.candidates.ABC)

subset.outside.genome.ABC = all.hic[-c(queryHits(Overlaps.Regul.ABC), queryHits(Overlaps.candidates.Regul.ABC), queryHits(Overlaps.Prom.ABC), queryHits(Overlaps.candidates.Prom.ABC))]

subset.outside.genome.ABC$SigniSNPs = countOverlaps(subset.outside.genome.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=.05])
subset.outside.genome.ABC$NSigniSNPs = countOverlaps(subset.outside.genome.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>.05])

subset.outside.genome.ABC$PropSigniSNPs = subset.outside.genome.ABC$SigniSNPs/(subset.outside.genome.ABC$SigniSNPs+subset.outside.genome.ABC$NSigniSNPs)

mean(subset.outside.genome.ABC$PropSigniSNPs, na.rm=T)

prop.outside.RRCs.ABC = data.frame(c(mean(regul.candidates.ABC$PropSigniSNPs, na.rm=T), mean(promoters.candidates.ABC$PropSigniSNPs, na.rm=T),mean(subset.outside.genome.ABC$PropSigniSNPs, na.rm=T)),c("Candidate.Regul", "Candidate.Prom", "Outside"))
colnames(prop.outside.RRCs.ABC) = c("prop.SigniSNPs","Type")
prop.outside.RRCs.ABC$complexity = c(-40,-20,0)

prop.RRCs.ABC = asso.complexity.ABC.FULL.imput[,c("prop.SigniSNPs","complexity")]
prop.RRCs.ABC$Type = "ABC"

all.prop.ABC= rbind(prop.outside.RRCs.ABC,prop.RRCs.ABC)



r=ggplot(all.prop.ABC[!is.na(all.prop.ABC$prop.SigniSNPs),], aes(x=complexity, y=prop.SigniSNPs))+geom_jitter(aes(shape=Type, colour=Type))


Overlaps.Regul.Rao = findOverlaps(all.hic, unique.Regul.Rao)
Overlaps.candidates.Regul.Rao = findOverlaps(all.hic, regul.candidates.Rao)
Overlaps.Prom.Rao = findOverlaps(all.hic, unique.Promoters.Rao)
Overlaps.candidates.Prom.Rao = findOverlaps(all.hic, promoters.candidates.Rao)

subset.outside.genome.Rao = all.hic[-c(queryHits(Overlaps.Regul.Rao), queryHits(Overlaps.candidates.Regul.Rao), queryHits(Overlaps.Prom.Rao), queryHits(Overlaps.candidates.Prom.Rao))]

subset.outside.genome.Rao$SigniSNPs = countOverlaps(subset.outside.genome.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=.05])
subset.outside.genome.Rao$NSigniSNPs = countOverlaps(subset.outside.genome.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>.05])

subset.outside.genome.Rao$PropSigniSNPs = subset.outside.genome.Rao$SigniSNPs/(subset.outside.genome.Rao$SigniSNPs+subset.outside.genome.Rao$NSigniSNPs)

mean(subset.outside.genome.Rao$PropSigniSNPs, na.rm=T)

prop.outside.RRCs.Rao = data.frame(c(mean(regul.candidates.Rao$PropSigniSNPs, na.rm=T), mean(promoters.candidates.Rao$PropSigniSNPs, na.rm=T),mean(subset.outside.genome.Rao$PropSigniSNPs, na.rm=T)),c("Candidate.Regul", "Candidate.Prom", "Outside"))
colnames(prop.outside.RRCs.Rao) = c("prop.SigniSNPs","Type")
prop.outside.RRCs.Rao$complexity = c(-10,-5,0)

prop.RRCs.Rao = asso.complexity.Rao.FULL.imput[,c("prop.SigniSNPs","complexity")]
prop.RRCs.Rao$Type = "Rao"

all.prop.Rao= rbind(prop.outside.RRCs.Rao,prop.RRCs.Rao)

s=ggplot(all.prop.Rao[!is.na(all.prop.Rao$prop.SigniSNPs),], aes(x=complexity, y=prop.SigniSNPs))+geom_jitter(aes(shape=Type, colour=Type))
s


Overlaps.Regul.DNAse = findOverlaps(all.hic, unique.Regul.DNAse)
Overlaps.candidates.Regul.DNAse = findOverlaps(all.hic, regul.candidates.DNAse)
Overlaps.Prom.DNAse = findOverlaps(all.hic, unique.Promoters.DNAse)
Overlaps.candidates.Prom.DNAse = findOverlaps(all.hic, promoters.candidates.DNAse)

subset.outside.genome.DNAse = all.hic[-c(queryHits(Overlaps.Regul.DNAse), queryHits(Overlaps.candidates.Regul.DNAse), queryHits(Overlaps.Prom.DNAse), queryHits(Overlaps.candidates.Prom.DNAse))]



subset.outside.genome.DNAse$SigniSNPs = countOverlaps(subset.outside.genome.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=.05])
subset.outside.genome.DNAse$NSigniSNPs = countOverlaps(subset.outside.genome.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>.05])

subset.outside.genome.DNAse$PropSigniSNPs = subset.outside.genome.DNAse$SigniSNPs/(subset.outside.genome.DNAse$SigniSNPs+subset.outside.genome.DNAse$NSigniSNPs)

mean(subset.outside.genome.DNAse$PropSigniSNPs, na.rm=T)

prop.outside.RRCs.DNAse = data.frame(c(mean(regul.candidates.DNAse$PropSigniSNPs, na.rm=T), mean(promoters.candidates.DNAse$PropSigniSNPs, na.rm=T),mean(subset.outside.genome.DNAse$PropSigniSNPs, na.rm=T)),c("Candidate.Regul", "Candidate.Prom", "Outside"))
colnames(prop.outside.RRCs.DNAse) = c("prop.SigniSNPs","Type")
prop.outside.RRCs.DNAse$complexity = c(-10,-5,0)

prop.RRCs.DNAse = asso.complexity.DNAse.FULL.imput[,c("prop.SigniSNPs","complexity")]
prop.RRCs.DNAse$Type = "DNAse"

all.prop.DNAse= rbind(prop.outside.RRCs.DNAse,prop.RRCs.DNAse)

t=ggplot(all.prop.DNAse[!is.na(all.prop.DNAse$prop.SigniSNPs),], aes(x=complexity, y=prop.SigniSNPs))+geom_jitter(aes(shape=Type, colour=Type))
t

multiplot(r,s,t)


signi.regul.ABC = sum(countOverlaps(unique.Regul.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
nsigni.regul.ABC = sum(countOverlaps(unique.Regul.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.05]))

signi.candidates.regul.ABC = sum(countOverlaps(regul.candidates.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
nsigni.candidates.regul.ABC = sum(countOverlaps(regul.candidates.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.05]))

signi.prom.ABC = sum(countOverlaps(unique.Promoters.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
nsigni.prom.ABC = sum(countOverlaps(unique.Promoters.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.05]))

signi.candidates.prom.ABC = sum(countOverlaps(promoters.candidates.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
nsigni.candidates.prom.ABC= sum(countOverlaps(promoters.candidates.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.05]))


signi.outside.ABC = sum(countOverlaps(subset.outside.genome.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
nsigni.outside.ABC = sum(countOverlaps(subset.outside.genome.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.05]))

fisher.test(matrix(c(signi.regul.ABC, nsigni.regul.ABC,signi.candidates.regul.ABC, nsigni.candidates.regul.ABC), ncol=2), alternative="greater")
fisher.test(matrix(c(signi.prom.ABC, nsigni.prom.ABC,signi.candidates.prom.ABC, nsigni.candidates.prom.ABC), ncol=2), alternative="greater")
fisher.test(matrix(c(signi.regul.ABC+signi.prom.ABC, nsigni.regul.ABC+nsigni.prom.ABC, signi.candidates.regul.ABC+signi.candidates.prom.ABC, nsigni.candidates.regul.ABC+nsigni.candidates.prom.ABC, signi.outside.ABC, nsigni.outside.ABC), ncol=3),workspace = 2e8, alternative = "greater")


df.regul.ABC = data.frame("regul", countOverlaps(unique.Regul.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
colnames(df.regul.ABC) = c("type", "count")
df.candidates.regul.ABC = data.frame("candidate", countOverlaps(regul.candidates.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
colnames(df.candidates.regul.ABC)= colnames(df.regul.ABC)
df.outside.ABC = data.frame("outside",countOverlaps(subset.outside.genome.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
colnames(df.outside.ABC) = colnames(df.regul.ABC)

df.all.regul.ABC = rbind(df.outside.ABC,df.candidates.regul.ABC,df.regul.ABC)

df.prom.ABC = data.frame("prom", countOverlaps(unique.Promoters.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
colnames(df.prom.ABC) = c("type", "count")
df.candidates.prom.ABC = data.frame("candidate", countOverlaps(promoters.candidates.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
colnames(df.candidates.prom.ABC)= colnames(df.prom.ABC)


df.all.prom.ABC = rbind(df.outside.ABC,df.candidates.prom.ABC,df.prom.ABC)


count.model.regul.ABC = glm(count~type, df.all.regul.ABC, family=poisson)
count.model.prom.ABC = glm(count~type, df.all.prom.ABC, family=poisson)

summary(count.model.regul.ABC)
summary(count.model.prom.ABC)

signi.regul.Rao = sum(countOverlaps(unique.Regul.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
nsigni.regul.Rao = sum(countOverlaps(unique.Regul.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.0]))

signi.candidates.regul.Rao = sum(countOverlaps(regul.candidates.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
nsigni.candidates.regul.Rao = sum(countOverlaps(regul.candidates.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.05]))

signi.prom.Rao = sum(countOverlaps(unique.Promoters.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
nsigni.prom.Rao = sum(countOverlaps(unique.Promoters.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.05]))

signi.candidates.prom.Rao = sum(countOverlaps(promoters.candidates.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
nsigni.candidates.prom.Rao = sum(countOverlaps(promoters.candidates.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.05]))



signi.outside.Rao = sum(countOverlaps(subset.outside.genome.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
nsigni.outside.Rao = sum(countOverlaps(subset.outside.genome.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.05]))

fisher.test(matrix(c(signi.regul.Rao, nsigni.regul.Rao,signi.candidates.regul.Rao, nsigni.candidates.regul.Rao), ncol=2), alternative="greater")
fisher.test(matrix(c(signi.prom.Rao, nsigni.prom.Rao,signi.candidates.prom.Rao, nsigni.candidates.prom.Rao), ncol=2), alternative="greater")

fisher.test(matrix(c(signi.regul.Rao+signi.prom.Rao, nsigni.regul.Rao+nsigni.prom.Rao, signi.candidates.regul.Rao+signi.candidates.prom.Rao, nsigni.candidates.regul.Rao+nsigni.candidates.prom.Rao, signi.outside.Rao, nsigni.outside.Rao), ncol=3),workspace = 2e8, alternative = "greater")

df.regul.Rao = data.frame("regul", countOverlaps(unique.Regul.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
colnames(df.regul.Rao) = c("type", "count")
df.candidates.regul.Rao = data.frame("candidate", countOverlaps(regul.candidates.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
colnames(df.candidates.regul.Rao)= colnames(df.regul.Rao)
df.outside.Rao = data.frame("outside",countOverlaps(subset.outside.genome.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
colnames(df.outside.Rao) = colnames(df.regul.Rao)

df.all.regul.Rao = rbind(df.outside.Rao,df.candidates.regul.Rao,df.regul.Rao)

df.prom.Rao = data.frame("prom", countOverlaps(unique.Promoters.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
colnames(df.prom.Rao) = c("type", "count")
df.candidates.prom.Rao = data.frame("candidate", countOverlaps(promoters.candidates.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
colnames(df.candidates.prom.Rao)= colnames(df.prom.Rao)


df.all.prom.Rao = rbind(df.outside.Rao,df.candidates.prom.Rao,df.prom.Rao)


count.model.regul.Rao = glm(count~type, df.all.regul.Rao, family=poisson)
count.model.prom.Rao = glm(count~type, df.all.prom.Rao, family=poisson)
 

summary(count.model.regul.Rao)
summary(count.model.prom.Rao)

signi.regul.DNAse = sum(countOverlaps(unique.Regul.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
nsigni.regul.DNAse = sum(countOverlaps(unique.Regul.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.05]))

signi.candidates.regul.DNAse = sum(countOverlaps(regul.candidates.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
nsigni.candidates.regul.DNAse = sum(countOverlaps(regul.candidates.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.05]))


signi.prom.DNAse = sum(countOverlaps(unique.Promoters.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
nsigni.prom.DNAse = sum(countOverlaps(unique.Promoters.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.05]))

signi.candidates.prom.DNAse = sum(countOverlaps(promoters.candidates.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
nsigni.candidates.prom.DNAse = sum(countOverlaps(promoters.candidates.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.05]))



signi.outside.DNAse = sum(countOverlaps(subset.outside.genome.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
nsigni.outside.DNAse = sum(countOverlaps(subset.outside.genome.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.05]))

fisher.test(matrix(c(signi.regul.DNAse, nsigni.regul.DNAse,signi.candidates.regul.DNAse, nsigni.candidates.regul.DNAse), ncol=2), alternative="greater")
fisher.test(matrix(c(signi.prom.DNAse, nsigni.prom.DNAse,signi.candidates.prom.DNAse, nsigni.candidates.prom.DNAse), ncol=2), alternative="greater")

fisher.test(matrix(c(signi.regul.DNAse+signi.prom.DNAse, nsigni.regul.DNAse+nsigni.prom.DNAse, signi.candidates.regul.DNAse+signi.candidates.prom.DNAse, nsigni.candidates.regul.DNAse+nsigni.candidates.prom.DNAse, signi.outside.DNAse, nsigni.outside.DNAse), ncol=3),workspace = 2e8, alternative="greater")

df.regul.DNAse = data.frame("regul", countOverlaps(unique.Regul.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
colnames(df.regul.DNAse) = c("type", "count")
df.candidates.regul.DNAse = data.frame("candidate", countOverlaps(regul.candidates.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
colnames(df.candidates.regul.DNAse)= colnames(df.regul.DNAse)
df.outside.DNAse = data.frame("outside",countOverlaps(subset.outside.genome.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
colnames(df.outside.DNAse) = colnames(df.regul.DNAse)

df.all.regul.DNASe = rbind(df.outside.DNAse,df.candidates.regul.DNAse,df.regul.DNAse)

df.prom.DNAse = data.frame("prom", countOverlaps(unique.Promoters.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
colnames(df.prom.DNAse) = c("type", "count")
df.candidates.prom.DNAse = data.frame("candidate", countOverlaps(promoters.candidates.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05]))
colnames(df.candidates.prom.DNAse)= colnames(df.prom.DNAse)


df.all.prom.DNASe = rbind(df.outside.DNAse,df.candidates.prom.DNAse,df.prom.DNAse)


library(pscl)
count.model.regul.DNAse = glm(count~type, df.all.regul.DNASe, family=poisson)

E2 <- resid(count.model.regul.DNAse, type = "pearson")
N  <- nrow(df.all.regul.DNASe)
p  <- length(coef(count.model.regul.DNAse))  
sum(E2^2) / (N - p)

summary(count.model.regul.DNAse)


count.model.prom.DNAse = glm(count~type, df.all.prom.DNASe, family=poisson)

E2 <- resid(count.model.prom.DNAse, type = "pearson")
N  <- nrow(df.all.prom.DNASe)
p  <- length(coef(count.model.prom.DNAse))  
sum(E2^2) / (N - p)

summary(count.model.prom.DNAse)

zero.infl.regul.DNAse=  zeroinfl(count ~ type | 
                 type, 
               dist = 'poisson',
               data = df.all.regul.DNASe)

summary(zero.infl.regul.DNAse)
E2 <- resid(zero.infl.regul.DNAse, type = "pearson")
N  <- nrow(df.all.regul.DNASe)
p  <- length(coef(zero.infl.regul.DNAse))  
sum(E2^2) / (N - p)

zero.infl.prom.DNAse=  zeroinfl(count ~ type | 
                                   type, 
                                 dist = 'poisson',
                                 data = df.all.prom.DNASe)

summary(zero.infl.prom.DNAse)
E2 <- resid(zero.infl.prom.DNAse, type = "pearson")
N  <- nrow(df.all.prom.DNASe)
p  <- length(coef(zero.infl.prom.DNAse))  
sum(E2^2) / (N - p)



#Analysis of CRNs with 0 significant SNPs 
summary(asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$prop.SigniSNPs>0,"complexity"])
summary(asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$prop.SigniSNPs==0,"complexity"])

summary(asso.complexity.Rao.FULL.imput[asso.complexity.Rao.FULL.imput$prop.SigniSNPs>0,"complexity"])
summary(asso.complexity.Rao.FULL.imput[asso.complexity.Rao.FULL.imput$prop.SigniSNPs==0,"complexity"])

summary(asso.complexity.DNAse.FULL.imput[asso.complexity.DNAse.FULL.imput$prop.SigniSNPs>0,"complexity"])
summary(asso.complexity.DNAse.FULL.imput[asso.complexity.DNAse.FULL.imput$prop.SigniSNPs==0,"complexity"])

#More Complex CRNs have more Significant SNPs thant less complex networks, no matter what annotation method is considered

summary(regression.propSNPs(regul.candidates.ABC, promoters.candidates.ABC, asso.complexity.ABC.FULL.imput, model="beta"))
summary(regression.propSNPs(regul.candidates.Rao, promoters.candidates.Rao, asso.complexity.Rao.FULL.imput, model="beta"))
summary(regression.propSNPs(regul.candidates.DNAse, promoters.candidates.DNAse, asso.complexity.DNAse.FULL.imput, model="beta"))

kable(tidy(regression.propSNPs(regul.candidates.ABC, promoters.candidates.ABC, asso.complexity.ABC.FULL.imput, model="beta")),"html")%>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>% 
  save_kable(file = "./model_propSNPs_ABC.png")

kable(tidy(regression.propSNPs(regul.candidates.Rao, promoters.candidates.Rao, asso.complexity.Rao.FULL.imput, model="beta")),"html")%>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>% 
  save_kable(file = "./model_propSNPs_Rao.png")

kable(tidy(regression.propSNPs(regul.candidates.DNAse, promoters.candidates.DNAse, asso.complexity.DNAse.FULL.imput, model="beta")),"html")%>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>% 
  save_kable(file = "./model_propSNPs_DNAse.png")


#Fisher test
sum(countOverlaps(regul.candidates.ABC,GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.0001]),countOverlaps(promoters.candidates.ABC,GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.0001]))
sum(countOverlaps(regul.candidates.ABC,GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.0001]),countOverlaps(promoters.candidates.ABC,GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.0001]))

sum(countOverlaps(unique.Regul.ABC,GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.0001]),countOverlaps(unique.Promoters.ABC,GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.0001]))
sum(countOverlaps(unique.Regul.ABC,GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.0001]),countOverlaps(unique.Promoters.ABC,GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.0001]))


ty = sum(countOverlaps(others.forward.ABC,GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.0001]))+sum(countOverlaps(others.reverse.ABC,GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.0001]))
tu= sum(countOverlaps(others.forward.ABC,GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.0001]))+sum(countOverlaps(others.reverse.ABC,GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.0001]))

to = sum(countOverlaps(regul.candidates.ABC,GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.0001]))+sum(countOverlaps(promoters.candidates.ABC,GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.0001]))
tp = sum(countOverlaps(regul.candidates.ABC,GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.0001]))+sum(countOverlaps(promoters.candidates.ABC,GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.0001]))

ta= sum(countOverlaps(unique.Regul.ABC,GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.0001]))+sum(countOverlaps(unique.Promoters.ABC,GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.0001]))
ts = sum(countOverlaps(unique.Regul.ABC,GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.0001]))+sum(countOverlaps(unique.Promoters.ABC,GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.0001]))

chisq.test(matrix(c(ty,tu,ta,ts), ncol=2))
#Model for the proportion of essential genes in CRNs
summary(betareg(prop.EssentialGenes~ DI.PROM+DI.ENH +INS.PROM+INS.ENH+mean_CTCF.PROM+mean_CTCF.ENH+mean_DNAse.PROM+mean_DNAse.ENH+mean_H3K27ac.PROM+
          mean_H3K27ac.ENH +mean_H3K27me3.PROM + mean_H3K27me3.ENH +mean_H3K4me1.PROM +mean_H3K4me1.ENH+mean_H3K4me3.PROM+mean_H3K4me3.ENH+mean_p300.PROM+
          mean_p300.ENH + ScoreFire.PROM + ScoreFire.ENH + Expression.y+complexity,data=asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$prop.EssentialGenes>0&asso.complexity.ABC.FULL.imput$prop.EssentialGenes<1,], link = "logit"))

summary(betareg(prop.EssentialGenes~ mean_DI.x+mean_DI.y +mean_INS.x+mean_INS.y+mean_CTCF.x+mean_CTCF.y+mean_DNAse.x+mean_DNAse.y+mean_H3K27ac.x+
                  mean_H3K27ac.y +mean_H3K27me3.x +mean_H3K27me3.y +mean_H3K4me1.x +mean_H3K4me1.y+mean_H3K4me3.x+mean_H3K4me3.y+mean_p300.x+
                  mean_p300.y + mean_ScoreFire.x+ mean_ScoreFire.y + Expression.y+complexity,data=asso.complexity.Rao.FULL.imput[asso.complexity.Rao.FULL.imput$prop.EssentialGenes>0&asso.complexity.Rao.FULL.imput$prop.EssentialGenes<1,], link = "logit"))

summary(betareg(prop.EssentialGenes~ mean_DI.x+mean_DI.y +mean_INS.x+mean_INS.y+mean_CTCF.x+mean_CTCF.y+mean_DNAse.x+mean_DNAse.y+mean_H3K27ac.x+
                  mean_H3K27ac.y +mean_H3K27me3.x +mean_H3K27me3.y +mean_H3K4me1.x +mean_H3K4me1.y+mean_H3K4me3.x+mean_H3K4me3.y+mean_p300.x+
                  mean_p300.y + mean_ScoreFire.x+ mean_ScoreFire.y + Expression.y+complexity,data=asso.complexity.DNAse.FULL.imput[asso.complexity.DNAse.FULL.imput$prop.EssentialGenes>0&asso.complexity.DNAse.FULL.imput$prop.EssentialGenes<1,], link = "logit"))

library(clusterProfiler)
a = lapply(unique(regul.promoters.ABC$membership), function(x){
  
  u = unique(regul.promoters.ABC[regul.promoters.ABC$membership==x,"TargetGene"])
  
  return(mcols(genes.hg19)[genes.hg19$geneSymbol%in%u,"gene_id"])
})
names(a) = paste0("X", 1:length(a))
a.bis = Filter(function(x) length(x) >= 10 & length(x)<=50, a)

b = lapply(unique(df.Rao$membership), function(x){
  
  u = unique(df.Rao[df.Rao$membership==x,"geneSymbol.x"])
  
  return(mcols(genes.hg19)[genes.hg19$geneSymbol%in%u,"gene_id"])
})
names(b) = paste0("X", 1:length(b))
b.bis = Filter(function(x) length(x) >= 3 & length(x)<=50, b)

c = lapply(unique(df.DNase$membership), function(x){
  
  u = unique(df.DNase[df.DNase$membership==x,"geneSymbol.x"])
  
  return(mcols(genes.hg19)[genes.hg19$geneSymbol%in%u,"gene_id"])
})
names(c) = paste0("X", 1:length(c))
c.bis = Filter(function(x) length(x) >= 3 & length(x)<=50, c)

library(DOSE)
library(ReactomePA)

ck.Kegg = compareCluster(geneCluster = a.bis, fun = 'enrichKEGG')
ck.Do <- compareCluster(geneCluster = a.bis, fun = 'enrichDO')
ck.Pathway = compareCluster(geneCluster = a.bis, fun="enrichPathway")

dotplot(ck.Kegg)
dotplot(ck.Do)
dotplot(ck.Pathway)


ck.Kegg.Rao = compareCluster(geneCluster = b.bis, fun = 'enrichKEGG')
ck.Do.Rao <- compareCluster(geneCluster = b.bis, fun = 'enrichDO')
ck.Pathway.Rao = compareCluster(geneCluster = b.bis, fun="enrichPathway")

dotplot(ck.Kegg.Rao)
dotplot(ck.Do.Rao)
dotplot(ck.Pathway.Rao)


ck.Kegg.DNAse = compareCluster(geneCluster = c.bis, fun = 'enrichKEGG')
ck.Do.DNAse <- compareCluster(geneCluster = c.bis, fun = 'enrichDO')
ck.Pathway.DNAse = compareCluster(geneCluster = c.bis, fun="enrichPathway")

dotplot(ck.Kegg.DNAse)
dotplot(ck.Do.DNAse)
dotplot(ck.Pathway.DNAse)


p.threshold = c(0.05,0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001)

meanDist.SNPs.ABC.regul = matrix(NA, ncol=3, nrow=7)
meanDist.SNPs.ABC.prom = matrix(NA, ncol=3, nrow=7)
meanDist.SNPs.Rao.regul = matrix(NA, ncol=3, nrow=7)
meanDist.SNPs.Rao.prom = matrix(NA, ncol=3, nrow=7)
meanDist.SNPs.DNAse.regul = matrix(NA, ncol=3, nrow=7)
meanDist.SNPs.DNAse.prom = matrix(NA, ncol=3, nrow=7)

for(i in 1:length(p.threshold)){
  distance.SNPs.regul.ABC = distanceToNearest(GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=p.threshold[i]], unique.Regul.ABC)
  distance.NSNPs.regul.ABC = distanceToNearest(GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>p.threshold[i]], unique.Regul.ABC)
  
  distance.SNPs.prom.ABC = distanceToNearest(GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=p.threshold[i]], unique.Promoters.ABC)
  distance.NSNPs.prom.ABC = distanceToNearest(GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>p.threshold[i]], unique.Promoters.ABC)
  
  distance.SNPs.regul.Rao = distanceToNearest(GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=p.threshold[i]], unique.Regul.Rao)
  distance.NSNPs.regul.Rao = distanceToNearest(GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>p.threshold[i]], unique.Regul.Rao)
  
  distance.SNPs.prom.Rao = distanceToNearest(GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=p.threshold[i]], unique.Promoters.Rao)
  distance.NSNPs.prom.Rao = distanceToNearest(GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>p.threshold[i]], unique.Promoters.Rao)
  
  distance.SNPs.regul.DNAse = distanceToNearest(GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=p.threshold[i]], unique.Regul.DNAse)
  distance.NSNPs.regul.DNAse = distanceToNearest(GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>p.threshold[i]], unique.Regul.DNAse)
  
  distance.SNPs.prom.DNAse = distanceToNearest(GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=p.threshold[i]], unique.Promoters.DNAse)
  distance.NSNPs.prom.DNAse = distanceToNearest(GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>p.threshold[i]], unique.Promoters.DNAse)
  
  meanDist.SNPs.ABC.regul[i,]=c("ABC.regul",p.threshold[i],mean(data.frame(distance.SNPs.regul.ABC)$distance)-mean(data.frame(distance.NSNPs.regul.ABC)$distance))
  meanDist.SNPs.ABC.prom[i,]=c("ABC.prom",p.threshold[i],mean(data.frame(distance.SNPs.prom.ABC)$distance)-mean(data.frame(distance.NSNPs.prom.ABC)$distance))
  
  meanDist.SNPs.Rao.regul[i,]=c("Rao.regul",p.threshold[i],mean(data.frame(distance.SNPs.regul.Rao)$distance)-mean(data.frame(distance.NSNPs.regul.Rao)$distance))
  meanDist.SNPs.Rao.prom[i,]=c("Rao.prom",p.threshold[i],mean(data.frame(distance.SNPs.prom.Rao)$distance)-mean(data.frame(distance.NSNPs.prom.Rao)$distance))
  
  meanDist.SNPs.DNAse.regul[i,]=c("DNAse.regul",p.threshold[i],mean(data.frame(distance.SNPs.regul.DNAse)$distance)-mean(data.frame(distance.NSNPs.regul.DNAse)$distance))
  meanDist.SNPs.DNAse.prom[i,]=c("DNAse.prom",p.threshold[i],mean(data.frame(distance.SNPs.prom.DNAse)$distance)-mean(data.frame(distance.NSNPs.prom.DNAse)$distance))
  
}
meanDist.SNPs = data.frame(rbind(meanDist.SNPs.ABC.regul,meanDist.SNPs.ABC.prom,meanDist.SNPs.Rao.regul,meanDist.SNPs.Rao.prom,meanDist.SNPs.DNAse.regul,meanDist.SNPs.DNAse.prom))
colnames(meanDist.SNPs) = c("Annotation_method", "p.value_threshold", "diff.meanDist")
meanDist.SNPs$diff.meanDist = as.numeric(as.character(meanDist.SNPs$diff.meanDist))
ggplot(meanDist.SNPs, aes(x=p.value_threshold, y=diff.meanDist, group=Annotation_method))+geom_line(aes(color=Annotation_method))+ ylab("diff.meanDist")+xlab("p.value threshold")


OR.SNPs.enrichment.ABC.regul = matrix(NA, ncol=3, nrow=7)
OR.SNPs.enrichment.ABC.prom = matrix(NA, ncol=3, nrow=7)
OR.SNPs.enrichment.Rao.regul = matrix(NA, ncol=3, nrow=7)
OR.SNPs.enrichment.Rao.prom = matrix(NA, ncol=3, nrow=7)
OR.SNPs.enrichment.DNAse.regul = matrix(NA, ncol=3, nrow=7)
OR.SNPs.enrichment.DNAse.prom = matrix(NA, ncol=3, nrow=7)

for(i in 1:length(p.threshold)){
  OR.SNPs.enrichment.ABC.regul[i,]=c("ABC.regul",p.threshold[i],enrichments.analysis(unique.Regul.ABC, regul.candidates.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=p.threshold[i]],GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>p.threshold[i]])$estimate)
  OR.SNPs.enrichment.ABC.prom[i,]=c("ABC.prom",p.threshold[i],enrichments.analysis(unique.Promoters.ABC, promoters.candidates.ABC, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=p.threshold[i]],GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>p.threshold[i]])$estimate)
  OR.SNPs.enrichment.Rao.regul[i,]=c("Rao.regul",p.threshold[i],enrichments.analysis(unique.Regul.Rao, regul.candidates.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=p.threshold[i]],GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>p.threshold[i]])$estimate)
  OR.SNPs.enrichment.Rao.prom[i,]=c("Rao.prom",p.threshold[i],enrichments.analysis(unique.Promoters.Rao, promoters.candidates.Rao, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=p.threshold[i]],GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>p.threshold[i]])$estimate)
  OR.SNPs.enrichment.DNAse.regul[i,]=c("DNAse.regul",p.threshold[i],enrichments.analysis(unique.Regul.DNAse, regul.candidates.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=p.threshold[i]],GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>p.threshold[i]])$estimate)
  OR.SNPs.enrichment.DNAse.prom[i,]=c("DNAse.prom",p.threshold[i],enrichments.analysis(unique.Promoters.DNAse, promoters.candidates.DNAse, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=p.threshold[i]],GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>p.threshold[i]])$estimate)
  
  }

OR.SNPs.enrichment = data.frame(rbind(OR.SNPs.enrichment.ABC.regul,OR.SNPs.enrichment.ABC.prom,OR.SNPs.enrichment.Rao.regul,OR.SNPs.enrichment.Rao.prom,OR.SNPs.enrichment.DNAse.regul,OR.SNPs.enrichment.DNAse.prom))
colnames(OR.SNPs.enrichment) = c("Annotation_method", "p.value_threshold", "OR")
OR.SNPs.enrichment$OR = as.numeric(as.character(OR.SNPs.enrichment$OR))


ggplot(OR.SNPs.enrichment, aes(x=p.value_threshold, y=OR, group=Annotation_method))+geom_line(aes(color=Annotation_method))+ ylab("OR")+xlab("p.value threshold")
distance.SNPs.candidates.regul.ABC = distanceToNearest(GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.000001], regul.candidates.ABC)


#Enrichment analysis
KEGG.df.ABC = data.frame(ck.Kegg)
KEGG.df.Rao = data.frame(ck.Kegg.Rao)
KEGG.df.DNAse = data.frame(ck.Kegg.DNAse)

Do.df.ABC = data.frame(ck.Do)
Do.df.Rao = data.frame(ck.Do.Rao)
Do.df.DNAse = data.frame(ck.Do.DNAse)


pathway.df.ABC = data.frame(ck.Pathway)
pathway.df.Rao = data.frame(ck.Pathway.Rao)
pathway.df.DNAse = data.frame(ck.Pathway.DNAse)

causal.genes1.ABC = unique(strsplit(Do.df.ABC[Do.df.ABC$Description%in%c("developmental disorder of mental health"),"geneID"],"/", fixed=T))
causal.genes2.ABC = unique(strsplit(Do.df.ABC[Do.df.ABC$Description%in%c("mood disorder", "bipolar disorder"),"geneID"],"/", fixed=T))

unique(sapply(causal.genes1.ABC, function(x) {
  return(unique(regul.promoters.ABC.imput[regul.promoters.ABC.imput$TargetGene%in%mcols(genes.hg19)[genes.hg19$gene_id%in%x,"geneSymbol"],"membership"]))
}))

unique(sapply(causal.genes2.ABC, function(x) {
  return(unique(regul.promoters.ABC.imput[regul.promoters.ABC.imput$TargetGene%in%mcols(genes.hg19)[genes.hg19$gene_id%in%x,"geneSymbol"],"membership"]))
}))

pathway.functionnal.ABC = lapply(causal.genes1.ABC, function(x) sapply(x, function(y) grepl(y,pathway.df.ABC$geneID, fixed=T)))
KEGG.functionnal.ABC =  lapply(causal.genes1.ABC, function(x) sapply(x, function(y) grepl(y,KEGG.df.ABC$geneID, fixed=T)))

pathway.functionnal.ABC2 = lapply(causal.genes2.ABC, function(x) sapply(x, function(y) grepl(y,pathway.df.ABC$geneID, fixed=T)))
KEGG.functionnal.ABC2 =  lapply(causal.genes2.ABC, function(x) sapply(x, function(y) grepl(y,KEGG.df.ABC$geneID, fixed=T)))

pathway.RRCs.ABC = lapply(1:length(pathway.functionnal.ABC), function(x) lapply(1:ncol(pathway.functionnal.ABC[[x]]), function(y) pathway.df.ABC[pathway.functionnal.ABC[[x]][,y], ]))
KEGG.RRCs.ABC = lapply(1:length(KEGG.functionnal.ABC), function(x) lapply(1:ncol(KEGG.functionnal.ABC[[x]]), function(y) KEGG.df.ABC[KEGG.functionnal.ABC[[x]][,y], ]))

lapply(pathway.RRCs.ABC,function(x) unique(do.call("rbind", x)))
lapply(KEGG.RRCs.ABC,function(x) unique(do.call("rbind", x)))

pathway.RRCs.ABC2 = lapply(1:length(pathway.functionnal.ABC2), function(x) lapply(1:ncol(pathway.functionnal.ABC2[[x]]), function(y) pathway.df.ABC[pathway.functionnal.ABC2[[x]][,y], ]))
KEGG.RRCs.ABC2 = lapply(1:length(KEGG.functionnal.ABC2), function(x) lapply(1:ncol(KEGG.functionnal.ABC2[[x]]), function(y) KEGG.df.ABC[KEGG.functionnal.ABC2[[x]][,y], ]))

lapply(pathway.RRCs.ABC2,function(x) unique(do.call("rbind", x)))
lapply(KEGG.RRCs.ABC2,function(x) unique(do.call("rbind", x)))

library(rlist)
unique.Promoters.ABC.imput = imput.individual.elements(unique.Promoters.ABC,unique.Regul.ABC, regul.promoters.ABC, method="ABC")$prom
unique.Regul.ABC.imput = imput.individual.elements(unique.Promoters.ABC,unique.Regul.ABC, regul.promoters.ABC, method="ABC")$regul

permut1.ABC = compute.cluster.enrichment(regul.promoters.ABC.imput, 516, 10000, "ABC")

permut1.ABC.genes = permut1.ABC$genes
permut1.ABC.regul = permut1.ABC$regul

caract.permut1.ABC.genes = lapply(1:nrow(permut1.ABC.genes), function(x) unique.Promoters.ABC.imput[unique.Promoters.ABC.imput$TargetGene%in%permut1.ABC.genes[x,],])
caract.permut1.ABC.regul = lapply(1:nrow(permut1.ABC.regul), function(x) unique.Regul.ABC.imput[unique.Regul.ABC.imput$name%in%permut1.ABC.regul[x,],])

enrichment.functional.analysis(asso.complexity.ABC.FULL.imput, 516, caract.permut1.ABC.genes, "prom","ABC")
enrichment.functional.analysis(asso.complexity.ABC.FULL.imput, 516, caract.permut1.ABC.regul, "regul","ABC")


caract.permut1.ABC.genes = lapply(caract.permut1.ABC.genes, function(x){
  x$essentialGenes = ifelse(x$TargetGene%in%unique.essential.genes,1,0)
  return(x)
})

caract.permut1.ABC.genes = lapply(caract.permut1.ABC.genes, function(x){
  x= merge(x, unique(regul.promoters.ABC.imput[,c("TargetGene","Expression")]),by="TargetGene",all.x=T)
  return(x)
})

prop.essentialGenes.permut1.ABC.genes = sapply(caract.permut1.ABC.genes, function(x){
  sum(x$essentialGenes)/length(x$essentialGenes)
})

asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$membership==516,]

percentile.Expression.permut1.ABC.genes = sapply(caract.permut1.ABC.genes, function(x) quantile(x$Expression, prob=0.9,na.rm=T))

sum(percentile.Expression.permut1.ABC.genes>=asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$membership==516,"Expression.x"])/10000
#0.17

#Variant enrichment
caract.permut1.ABC.genes = lapply(caract.permut1.ABC.genes, function(x){
  x= merge(x, unique(regul.promoters.ABC.imput[,c("TargetGene","NSigni.PROM","NnSigni.PROM")]),by="TargetGene",all.x=T)
  return(x)
})

caract.permut1.ABC.regul = lapply(caract.permut1.ABC.regul, function(x){
  x= merge(x, unique(regul.promoters.ABC.imput[,c("name","NSigni.ENH","NnSigni.ENH")]),by="name",all.x=T)
  return(x)
})


prop.variants.permut1.ABC = sapply(1:10000, function(x){
  return(sum(caract.permut1.ABC.genes[[x]]$NSigni.PROM,caract.permut1.ABC.regul[[x]]$NSigni.ENH) /sum(caract.permut1.ABC.genes[[x]]$NSigni.PROM,caract.permut1.ABC.regul[[x]]$NSigni.ENH,
                                                                                                      caract.permut1.ABC.genes[[x]]$NnSigni.PROM,caract.permut1.ABC.regul[[x]]$NnSigni.ENH))
})
prop.variants.permut1.ABC[is.nan(prop.variants.permut1.ABC)] = 0
sum(prop.variants.permut1.ABC>=asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$membership==516,"prop.SigniSNPs"])/10000
#0.5696


permut2.ABC = compute.cluster.enrichment(regul.promoters.ABC.imput, 1432, 10000, "ABC")

permut2.ABC.genes = permut2.ABC$genes
permut2.ABC.regul = permut2.ABC$regul

caract.permut2.ABC.genes = lapply(1:nrow(permut2.ABC.genes), function(x) unique.Promoters.ABC.imput[unique.Promoters.ABC.imput$TargetGene%in%permut2.ABC.genes[x,],])
caract.permut2.ABC.regul = lapply(1:nrow(permut2.ABC.regul), function(x) unique.Regul.ABC.imput[unique.Regul.ABC.imput$name%in%permut2.ABC.regul[x,],])

enrichment.functional.analysis(asso.complexity.ABC.FULL.imput, 1432, caract.permut2.ABC.genes, "prom","ABC")
enrichment.functional.analysis(asso.complexity.ABC.FULL.imput, 1432, caract.permut2.ABC.regul, "regul","ABC")

caract.permut2.ABC.genes = lapply(caract.permut2.ABC.genes, function(x){
  x$essentialGenes = ifelse(x$TargetGene%in%unique.essential.genes,1,0)
  return(x)
})

prop.essentialGenes.permut2.ABC.genes = sapply(caract.permut2.ABC.genes, function(x){
  sum(x$essentialGenes)/length(x$essentialGenes)
})

caract.permut2.ABC.genes = lapply(caract.permut2.ABC.genes, function(x){
  x= merge(x, unique(regul.promoters.ABC.imput[,c("TargetGene","Expression")]),by="TargetGene",all.x=T)
  return(x)
})

asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$membership==1432,]

sum(prop.essentialGenes.permut2.ABC.genes>=asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$membership==1432,"prop.EssentialGenes"])/10000

percentile.Expression.permut2.ABC.genes = sapply(caract.permut2.ABC.genes, function(x) quantile(x$Expression, prob=0.9,na.rm=T))

sum(percentile.Expression.permut2.ABC.genes>=asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$membership==1432,"Expression.x"])/10000
#0.7411


#Variant enrichment

caract.permut2.ABC.genes = lapply(caract.permut2.ABC.genes, function(x){
  x= merge(x, unique(regul.promoters.ABC.imput[,c("TargetGene","NSigni.PROM","NnSigni.PROM")]),by="TargetGene",all.x=T)
  return(x)
})

caract.permut2.ABC.regul = lapply(caract.permut2.ABC.regul, function(x){
  x= merge(x, unique(regul.promoters.ABC.imput[,c("name","NSigni.ENH","NnSigni.ENH")]),by="name",all.x=T)
  return(x)
})


prop.variants.permut2.ABC = sapply(1:10000, function(x){
  return(sum(caract.permut2.ABC.genes[[x]]$NSigni.PROM,caract.permut2.ABC.regul[[x]]$NSigni.ENH) /sum(caract.permut2.ABC.genes[[x]]$NSigni.PROM,caract.permut2.ABC.regul[[x]]$NSigni.ENH,
                                                                                                      caract.permut2.ABC.genes[[x]]$NnSigni.PROM,caract.permut2.ABC.regul[[x]]$NnSigni.ENH))
})
prop.variants.permut2.ABC[is.nan(prop.variants.permut2.ABC)] = 0
sum(prop.variants.permut2.ABC>=asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$membership==1432,"prop.SigniSNPs"])/10000
#0.7783


permut3.ABC = compute.cluster.enrichment(regul.promoters.ABC.imput, 712, 10000, "ABC")

permut3.ABC.genes = permut3.ABC$genes
permut3.ABC.regul = permut3.ABC$regul

caract.permut3.ABC.genes = lapply(1:nrow(permut3.ABC.genes), function(x) unique.Promoters.ABC.imput[unique.Promoters.ABC.imput$TargetGene%in%permut3.ABC.genes[x,],])
caract.permut3.ABC.regul = lapply(1:nrow(permut3.ABC.regul), function(x) unique.Regul.ABC.imput[unique.Regul.ABC.imput$name%in%permut3.ABC.regul[x,],])

enrichment.functional.analysis(asso.complexity.ABC.FULL.imput, 712, caract.permut3.ABC.genes, "prom", "ABC")
enrichment.functional.analysis(asso.complexity.ABC.FULL.imput, 712, caract.permut3.ABC.regul, "regul","ABC")

caract.permut3.ABC.genes = lapply(caract.permut3.ABC.genes, function(x){
  x$essentialGenes = ifelse(x$TargetGene%in%unique.essential.genes,1,0)
  return(x)
})

prop.essentialGenes.permut3.ABC.genes = sapply(caract.permut3.ABC.genes, function(x){
  sum(x$essentialGenes)/length(x$essentialGenes)
})


caract.permut3.ABC.genes = lapply(caract.permut3.ABC.genes, function(x){
  x= merge(x, unique(regul.promoters.ABC.imput[,c("TargetGene","Expression")]),by="TargetGene",all.x=T)
  return(x)
})



asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$membership==712,]

sum(prop.essentialGenes.permut3.ABC.genes>=asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$membership==712,"prop.EssentialGenes"])/10000

percentile.Expression.permut3.ABC.genes = sapply(caract.permut3.ABC.genes, function(x) quantile(x$Expression, prob=0.9,na.rm=T))

sum(percentile.Expression.permut3.ABC.genes>=asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$membership==712,"Expression.x"])/10000
#0.7757



#Variant enrichment
caract.permut3.ABC.genes = lapply(caract.permut3.ABC.genes, function(x){
  x= merge(x, unique(regul.promoters.ABC.imput[,c("TargetGene","NSigni.PROM","NnSigni.PROM")]),by="TargetGene",all.x=T)
  return(x)
})

caract.permut3.ABC.regul = lapply(caract.permut3.ABC.regul, function(x){
  x= merge(x, unique(regul.promoters.ABC.imput[,c("name","NSigni.ENH","NnSigni.ENH")]),by="name",all.x=T)
  return(x)
})


prop.variants.permut3.ABC = sapply(1:10000, function(x){
  return(sum(caract.permut3.ABC.genes[[x]]$NSigni.PROM,caract.permut3.ABC.regul[[x]]$NSigni.ENH) /sum(caract.permut3.ABC.genes[[x]]$NSigni.PROM,caract.permut3.ABC.regul[[x]]$NSigni.ENH,
                                                                                                      caract.permut3.ABC.genes[[x]]$NnSigni.PROM,caract.permut3.ABC.regul[[x]]$NnSigni.ENH))
})
prop.variants.permut3.ABC[is.nan(prop.variants.permut3.ABC)] = 0
sum(prop.variants.permut3.ABC>=asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$membership==712,"prop.SigniSNPs"])/10000
#1

permut4.ABC = compute.cluster.enrichment(regul.promoters.ABC.imput, 710, 10000, "ABC")

permut4.ABC.genes = permut4.ABC$genes
permut4.ABC.regul = permut4.ABC$regul

caract.permut4.ABC.genes = lapply(1:nrow(permut4.ABC.genes), function(x) unique.Promoters.ABC.imput[unique.Promoters.ABC.imput$TargetGene%in%permut4.ABC.genes[x,],])
caract.permut4.ABC.regul = lapply(1:nrow(permut4.ABC.regul), function(x) unique.Regul.ABC.imput[unique.Regul.ABC.imput$name%in%permut4.ABC.regul[x,],])


enrichment.functional.analysis(asso.complexity.ABC.FULL.imput, 710, caract.permut4.ABC.genes, "prom", "ABC")
enrichment.functional.analysis(asso.complexity.ABC.FULL.imput, 710, caract.permut4.ABC.regul, "regul", "ABC")

#essential genes enrichmnent
caract.permut4.ABC.genes = lapply(caract.permut4.ABC.genes, function(x){
  x$essentialGenes = ifelse(x$TargetGene%in%unique.essential.genes,1,0)
  return(x)
})

prop.essentialGenes.permut4.ABC.genes = sapply(caract.permut4.ABC.genes, function(x){
  sum(x$essentialGenes)/length(x$essentialGenes)
})

#Expression enrichment
caract.permut4.ABC.genes = lapply(caract.permut4.ABC.genes, function(x){
  x= merge(x, unique(regul.promoters.ABC.imput[,c("TargetGene","Expression")]),by="TargetGene",all.x=T)
  return(x)
})
asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$membership==710,]

sum(prop.essentialGenes.permut4.ABC.genes>=asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$membership==710,"prop.EssentialGenes"])/10000

percentile.Expression.permut4.ABC.genes = sapply(caract.permut4.ABC.genes, function(x) quantile(x$Expression, prob=0.9,na.rm=T))

sum(percentile.Expression.permut4.ABC.genes>=asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$membership==710,"Expression.x"])/10000
#0.4668

#Variant enrichment

caract.permut4.ABC.genes = lapply(caract.permut4.ABC.genes, function(x){
  x= merge(x, unique(regul.promoters.ABC.imput[,c("TargetGene","NSigni.PROM","NnSigni.PROM")]),by="TargetGene",all.x=T)
  return(x)
})

caract.permut4.ABC.regul = lapply(caract.permut4.ABC.regul, function(x){
  x= merge(x, unique(regul.promoters.ABC.imput[,c("name","NSigni.ENH","NnSigni.ENH")]),by="name",all.x=T)
  return(x)
})


prop.variants.permut4.ABC = sapply(1:10000, function(x){
  return(sum(caract.permut4.ABC.genes[[x]]$NSigni.PROM,caract.permut4.ABC.regul[[x]]$NSigni.ENH) /sum(caract.permut4.ABC.genes[[x]]$NSigni.PROM,caract.permut4.ABC.regul[[x]]$NSigni.ENH,
                                                                                           caract.permut4.ABC.genes[[x]]$NnSigni.PROM,caract.permut4.ABC.regul[[x]]$NnSigni.ENH))
})
prop.variants.permut4.ABC[is.nan(prop.variants.permut4.ABC)] = 0
sum(prop.variants.permut4.ABC>=asso.complexity.ABC.FULL.imput[asso.complexity.ABC.FULL.imput$membership==710,"prop.SigniSNPs"])/10000
#1

unique.Promoters.Rao.imput = imput.individual.elements(unique.Promoters.Rao,unique.Regul.Rao, df.Rao, method="Rao")$prom
unique.Regul.Rao.imput = imput.individual.elements(unique.Promoters.Rao,unique.Regul.Rao, df.Rao, method="Rao")$regul


unique(sapply(strsplit(Do.df.Rao[Do.df.Rao$Description%in%c("developmental disorder of mental health"),"geneID"],"/", fixed=T), function(x) {
  return(unique(df.Rao.imput[df.Rao.imput$geneSymbol.x%in%mcols(genes.hg19)[genes.hg19$gene_id%in%x,"geneSymbol"],"membership"]))
}))

unique(sapply(strsplit(Do.df.Rao[Do.df.Rao$Description%in%c("autism spectrum disorder", "austistic disorder"),"geneID"],"/", fixed=T), function(x) {
  return(unique(df.Rao.imput[df.Rao.imput$geneSymbol.x%in%mcols(genes.hg19)[genes.hg19$gene_id%in%x,"geneSymbol"],"membership"]))
}))

causal.genes.Rao = unique(strsplit(Do.df.Rao[Do.df.Rao$Description%in%c("developmental disorder of mental health"),"geneID"],"/", fixed=T))

pathway.functionnal.Rao = lapply(causal.genes.Rao, function(x) sapply(x, function(y) grepl(y,pathway.df.Rao$geneID, fixed=T)))
KEGG.functionnal.Rao =  lapply(causal.genes.Rao, function(x) sapply(x, function(y) grepl(y,KEGG.df.Rao$geneID, fixed=T)))

pathway.RRCs.Rao = lapply(1:length(pathway.functionnal.Rao), function(x) lapply(1:ncol(pathway.functionnal.Rao[[x]]), function(y) pathway.df.Rao[pathway.functionnal.Rao[[x]][,y], ]))
KEGG.RRCs.Rao = lapply(1:length(KEGG.functionnal.Rao), function(x) lapply(1:ncol(KEGG.functionnal.Rao[[x]]), function(y) KEGG.df.Rao[KEGG.functionnal.Rao[[x]][,y], ]))


lapply(pathway.RRCs.Rao,function(x) unique(do.call("rbind", x)))
lapply(KEGG.RRCs.Rao,function(x) unique(do.call("rbind", x)))


permut.Rao = compute.cluster.enrichment(df.Rao.imput, 1752, 10000, "Rao")

permut.Rao.genes = permut.Rao$genes
permut.Rao.regul = permut.Rao$regul

caract.permut.Rao.genes = lapply(1:nrow(permut.Rao.genes), function(x) unique.Promoters.Rao.imput[unique.Promoters.Rao.imput$geneSymbol.x%in%permut.Rao.genes[x,],])
caract.permut.Rao.regul = lapply(1:nrow(permut.Rao.regul), function(x) unique.Regul.Rao.imput[unique.Regul.Rao.imput$name%in%permut.Rao.regul[x,],])

enrichment.functional.analysis(asso.complexity.Rao.FULL.imput, 1752, caract.permut.Rao.genes, "prom", "Rao")
enrichment.functional.analysis(asso.complexity.Rao.FULL.imput, 1752, caract.permut.Rao.regul, "regul", "Rao")

caract.permut.Rao.genes = lapply(caract.permut.Rao.genes, function(x){
  x$essentialGenes = ifelse(x$geneSymbol.x%in%unique.essential.genes,1,0)
  return(x)
})

asso.complexity.Rao.FULL.imput[asso.complexity.Rao.FULL.imput$membership==1752,]

caract.permut.Rao.genes = lapply(caract.permut.Rao.genes, function(x){
  x= merge(x, unique(df.Rao.imput[,c("geneSymbol.x","Expression")]),by="geneSymbol.x",all.x=T)
  return(x)
})

percentile.Expression.permut.Rao.genes = sapply(caract.permut.Rao.genes, function(x) quantile(x$Expression, prob=0.9,na.rm=T))

sum(percentile.Expression.permut.Rao.genes>=asso.complexity.Rao.FULL.imput[asso.complexity.Rao.FULL.imput$membership==1752,"Expression.x"])/10000
#0.9942

unique.Promoters.DNAse.imput = imput.individual.elements(unique.Promoters.DNAse,unique.Regul.DNAse, df.DNase, method="DNAse")$prom
unique.Regul.DNAse.imput = imput.individual.elements(unique.Promoters.DNAse,unique.Regul.DNAse, df.DNase, method="DNAse")$regul

unique(sapply(strsplit(Do.df.DNAse[Do.df.DNAse$Description%in%c("developmental disorder of mental health"),"geneID"],"/", fixed=T), function(x) {
  return(unique(df.DNase.imput[df.DNase.imput$geneSymbol.x%in%mcols(genes.hg19)[genes.hg19$gene_id%in%x,"geneSymbol"],"membership"]))
}))

unique(sapply(strsplit(Do.df.DNAse[Do.df.DNAse$Description%in%c("autism spectrum disorder", "autistic disorder"),"geneID"],"/", fixed=T), function(x) {
  return(unique(df.DNase.imput[df.DNase.imput$geneSymbol.x%in%mcols(genes.hg19)[genes.hg19$gene_id%in%x,"geneSymbol"],"membership"]))
}))

causal.genes.DNAse = unique(strsplit(Do.df.DNAse[Do.df.DNAse$Description%in%c("developmental disorder of mental health"),"geneID"],"/", fixed=T))

pathway.functionnal.DNAse = lapply(causal.genes.DNAse, function(x) sapply(x, function(y) grepl(y,pathway.df.DNAse$geneID, fixed=T)))
KEGG.functionnal.DNAse =  lapply(causal.genes.DNAse, function(x) sapply(x, function(y) grepl(y,KEGG.df.DNAse$geneID, fixed=T)))

pathway.RRCs.DNAse = lapply(1:length(pathway.functionnal.DNAse), function(x) lapply(1:ncol(pathway.functionnal.DNAse[[x]]), function(y) pathway.df.Rao[pathway.functionnal.DNAse[[x]][,y], ]))
KEGG.RRCs.DNAse = lapply(1:length(KEGG.functionnal.DNAse), function(x) lapply(1:ncol(KEGG.functionnal.DNAse[[x]]), function(y) KEGG.df.Rao[KEGG.functionnal.DNAse[[x]][,y], ]))


lapply(pathway.RRCs.DNAse,function(x) unique(do.call("rbind", x)))
lapply(KEGG.RRCs.DNAse,function(x) unique(do.call("rbind", x)))


asso.complexity.DNAse.FULL.imput[asso.complexity.DNAse.FULL.imput$membership==657,]

permut.DNAse = compute.cluster.enrichment(df.DNase.imput, 657, 10000, "DNAse")

permut.DNAse.genes = permut.DNAse$genes
permut.DNAse.regul = permut.DNAse$regul

caract.permut.DNAse.genes = lapply(1:nrow(permut.DNAse.genes), function(x) unique.Promoters.DNAse.imput[unique.Promoters.DNAse.imput$geneSymbol.x%in%permut.DNAse.genes[x,],])
caract.permut.DNAse.regul = lapply(1:nrow(permut.DNAse.regul), function(x) unique.Regul.DNAse.imput[unique.Regul.DNAse.imput$name%in%permut.DNAse.regul[x,],])

enrichment.functional.analysis(asso.complexity.DNAse.FULL.imput, 657, caract.permut.DNAse.genes, "prom", "DNAse")
enrichment.functional.analysis(asso.complexity.DNAse.FULL.imput, 657, caract.permut.DNAse.regul, "regul", "DNAse")

caract.permut.DNAse.genes = lapply(caract.permut.DNAse.genes, function(x){
  x= merge(x, unique(df.DNase.imput[,c("geneSymbol.x","Expression")]),by="geneSymbol.x",all.x=T)
  return(x)
})

percentile.Expression.permut.DNAse.genes = sapply(caract.permut.DNAse.genes, function(x) quantile(x$Expression, prob=0.9,na.rm=T))

sum(percentile.Expression.permut.DNAse.genes>=asso.complexity.DNAse.FULL.imput[asso.complexity.DNAse.FULL.imput$membership==657,"Expression.x"],na.rm=T)/10000
#[1] 0.9943



expression.genes.ABC = unique(regul.promoters.ABC.imput[,c("TargetGene","Expression")])
expression.genes.ABC$ABC = TRUE
expression.genes.notABC = mcols(expression.genes)[!expression.genes$geneSymbol%in%expression.genes.ABC$TargetGene,c("geneSymbol","median.Expr")]
colnames(expression.genes.notABC) = c("TargetGene", "Expression")
expression.genes.notABC$ABC = FALSE


kruskal.test(Expression~ABC, data=rbind(expression.genes.ABC,expression.genes.notABC))

aggregate(Expression~ABC, rbind(expression.genes.ABC,expression.genes.notABC),summary)

expression.genes.Rao = unique(df.Rao.imput[,c("geneSymbol.x","Expression")])
expression.genes.Rao$Rao = TRUE
expression.genes.notRao = mcols(expression.genes)[!expression.genes$geneSymbol%in%expression.genes.Rao$geneSymbol.x,c("geneSymbol","median.Expr")]
colnames(expression.genes.notRao) = c("geneSymbol.x", "Expression")
expression.genes.notRao$Rao = FALSE

kruskal.test(Expression~Rao, data=rbind(expression.genes.Rao,expression.genes.notRao))
aggregate(Expression~Rao, rbind(expression.genes.Rao,expression.genes.notRao),summary)

expression.genes.DNAse = unique(df.DNase.imput[,c("geneSymbol.x","Expression")])
expression.genes.DNAse$DNAse = TRUE
expression.genes.notDNAse = mcols(expression.genes)[!expression.genes$geneSymbol%in%expression.genes.DNAse$geneSymbol.x,c("geneSymbol","median.Expr")]
colnames(expression.genes.notDNAse) = c("geneSymbol.x", "Expression")
expression.genes.notDNAse$DNAse = FALSE

kruskal.test(Expression~DNAse, data=rbind(expression.genes.DNAse,expression.genes.notDNAse))
aggregate(Expression~DNAse, rbind(expression.genes.DNAse,expression.genes.notDNAse),summary)
