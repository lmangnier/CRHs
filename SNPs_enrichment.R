#Script pour l'analyse d'enrichissement
#Executer le script HiCrn.R avant d'executer ce code

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

library(ChIPseeker)
set.seed(1258)
setwd("/home/loic/Documents/HiC/data/4Script/NEU")
source("HiCrn.R")

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes.hg19 <- genes(txdb)
symbol <- select(org.Hs.eg.db,keys = genes.hg19$gene_id, columns = c("SYMBOL"), keytype = "ENTREZID")
genes.hg19$geneSymbol <- symbol$SYMBOL

regul.promoters.ABC = process.ABC("input_data/EnhancerPredictions.txt")
# all.putative.ABC = read.table("ABC/EnhancerPredictionsAllPutative.txt.gz", header = T)
# all.putative.ABC = all.putative.ABC[all.putative.ABC$chr!="chrX",]
# 
# all.putative.ABC.candidates = all.putative.ABC[!all.putative.ABC$name%in%regul.promoters.ABC$name&all.putative.ABC$class!="promoter",]

hic.loops = read.table("input_data/enriched_pixels_10000.bedpe", header=T)

peaks.NEU = import("input_data/NEU_DNAse.macs2_peaks.narrowPeak.sorted.candidateRegions.bed", format="bed")
peaks.NEU = peaks.NEU[seqnames(peaks.NEU)%in%paste0("chr", 1:22)]

GRanges.loops.bin1 = GRanges(seqnames = hic.loops$chr1, ranges=IRanges(start=hic.loops$x1, end=hic.loops$x2))
GRanges.loops.bin2 = GRanges(seqnames = hic.loops$chr2, ranges=IRanges(start=hic.loops$y1, end=hic.loops$y2))

Pairs.bins = Pairs(GRanges.loops.bin1, GRanges.loops.bin2,fdrBL = hic.loops$fdrBL ,fdrDonut = hic.loops$fdrDonut,fdrV=hic.loops$fdrV,fdrH = hic.loops$fdrH)
positive.contact = Pairs.bins[mcols(Pairs.bins)$fdrBL <= 0.15&mcols(Pairs.bins)$fdrDonut <=0.15&mcols(Pairs.bins)$fdrV <= 0.15&mcols(Pairs.bins)$fdrH <=0.15]


ABC.Pairs = Pairs.Enh.Prom.ABC(regul.promoters.ABC)
Rao.Pairs = process.Rao.DNAse(positive.contact, method = "Rao")
DNAse.Pairs = process.Rao.DNAse(positive.contact,peaks.NEU, method = "DNAse")

#Elements uniques inclus dans les RRCs 
unique.Promoters.ABC = unique(second(ABC.Pairs))
unique.Regul.ABC = unique(first(ABC.Pairs))
unique.Promoters.Rao = unique(first(Rao.Pairs))
unique.Regul.Rao = unique(second(Rao.Pairs))
unique.Promoters.DNAse = unique(first(DNAse.Pairs))
unique.Regul.DNAse =  unique(second(DNAse.Pairs))

#Ensembles equivalents
regul.candidates.ABC = make.comparable.set(unique.Regul.ABC,method="ABC", element="regulatory", DNAse = peaks.NEU)
promoters.candidates.ABC = make.comparable.set(unique.Promoters.ABC,method="ABC", element="promoter")
regul.candidates.Rao = make.comparable.set(unique.Regul.Rao,method="Rao", element="regulatory", contact =positive.contact )
promoters.candidates.Rao = make.comparable.set(unique.Promoters.Rao,method="Rao", element="promoter", contact = positive.contact)
regul.candidates.DNAse = make.comparable.set(unique.Promoters.DNAse,method="DNAse", element="regulatory",contact = positive.contact, DNAse = peaks.NEU)
promoters.candidates.DNAse = make.comparable.set(unique.Promoters.DNAse,method="DNAse", element="promoter", contact = positive.contact, DNAse = peaks.NEU)


#Fichier de variants clumpe
SCZ3.all.clumped = read.table("input_data/PGC3_SCZ_wave3_public.clumped.v2.tsv", header=T, fill = T)
GRanges.snps.SCZ3.clumped = GRanges(seqnames=paste0("chr",SCZ3.all.clumped$CHR), ranges=IRanges(start=SCZ3.all.clumped$BP, end=SCZ3.all.clumped$BP, names=SCZ3.all.clumped$SNP), pval=SCZ3.all.clumped$P)

#Differents seuils de significativite consideres
threshold = c(0.10,0.05,0.025,0.01,0.001,0.0001,0.00001,0.000001, 0.0000001)

subset.SNPs.by.signi = lapply(threshold, function(x) GRanges.snps.SCZ3.clumped[GRanges.snps.SCZ3.clumped$pval<=x])
subset.SNPs.by.nsigni = lapply(threshold, function(x) GRanges.snps.SCZ3.clumped[GRanges.snps.SCZ3.clumped$pval>x])

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peakAnno <- annotatePeak(peaks.NEU, tssRegion=c(-250, 250),
                         TxDb=txdb, annoDb="org.Hs.eg.db", level="gene")


unique.Regul.ABC.tmp = unique.Regul.ABC
unique.Regul.ABC.tmp$annotation = "ABC Enhancer"
unique.Regul.ABC.tmp$annotation_simplified = "ABC Enhancer"

GRange.peakAnno = as.GRanges(peakAnno)
GRange.peakAnno$annotation_simplified = ifelse(startsWith(GRange.peakAnno$annotation, "Exon"), "Exon",ifelse(startsWith(GRange.peakAnno$annotation, "Intron"), "Intron",GRange.peakAnno$annotation))

overlaps.ABC.peaks = findOverlaps(GRange.peakAnno, unique.Regul.ABC)
GRange.peakAnno.subset = GRange.peakAnno[-queryHits(overlaps.ABC.peaks)]

GRange.peakAnno.subset = do.call(c, GRangesList(GRange.peakAnno.subset,unique.Regul.ABC.tmp))
GRange.peakAnno.subset$nSNPs = countOverlaps(GRange.peakAnno.subset,GRanges.snps.SCZ3.clumped)

unique(GRange.peakAnno.subset$annotation_simplified)
aggregate(nSNPs~annotation_simplified, data.frame(GRange.peakAnno.subset), sum)$nSNPs/sum(aggregate(nSNPs~annotation_simplified, data.frame(GRange.peakAnno.subset), sum)$nSNPs)


rew = lapply(1:length(threshold), function(x) {
  GRange.peakAnno.subset$nSNPs_signi = countOverlaps(GRange.peakAnno.subset, subset.SNPs.by.signi[[x]])
  GRange.peakAnno.subset$nSNPs_nsigni = countOverlaps(GRange.peakAnno.subset, subset.SNPs.by.nsigni[[x]])
  
  signi = aggregate(nSNPs_signi~annotation_simplified, data.frame(GRange.peakAnno.subset), sum)$nSNPs_signi
  signi = signi[c(3,1,2,4,5,6,7,8,9,10)]
  nsigni = aggregate(nSNPs_nsigni~annotation_simplified, data.frame(GRange.peakAnno.subset), sum)$nSNPs_nsigni
  nsigni = nsigni[c(3,1,2,4,5,6,7,8,9,10)]
  
  names.ele = aggregate(nSNPs_signi~annotation_simplified, data.frame(GRange.peakAnno.subset), sum)$annotation_simplified
  names.ele = names.ele[c(3,1,2,4,5,6,7,8,9,10)]
  snps.tt = rbind(signi,nsigni)
  colnames(snps.tt)= names.ele
return(t(snps.tt))
})

rew.gene = lapply(1:length(threshold), function(x) {
  GRange.peakAnno$nSNPs_signi = countOverlaps(GRange.peakAnno, subset.SNPs.by.signi[[x]])
  GRange.peakAnno$nSNPs_nsigni = countOverlaps(GRange.peakAnno, subset.SNPs.by.nsigni[[x]])
  
  signi = aggregate(nSNPs_signi~annotation_simplified, data.frame(GRange.peakAnno), sum)$nSNPs_signi
  signi = signi[c(3,1,2,4,5,6,7,8,9)]
  nsigni = aggregate(nSNPs_nsigni~annotation_simplified, data.frame(GRange.peakAnno), sum)$nSNPs_nsigni
  nsigni = nsigni[c(3,1,2,4,5,6,7,8,9)]
  
  names.ele = aggregate(nSNPs_signi~annotation_simplified, data.frame(GRange.peakAnno), sum)$annotation_simplified
  names.ele = names.ele[c(3,1,2,4,5,6,7,8,9)]
  snps.tt = rbind(signi,nsigni)
  colnames(snps.tt)= names.ele
  return(t(snps.tt))
}
)


memory.limit(size=250000)
library(epitools)
lapply(1:length(threshold), function(x) epitab(rew.gene[[x]], correction=T))

#Ici pour chaque analyse d'enrichissement, on compare: 
#Les elements de regulation ou promoters avec les ensembles equivalents
#Les RRCs par rapport au reste du genome (on supprimme les elements presents dans les ensembles equivalents)
#Les sorties des differents boucles sont dans les elements RDS fournis dans le github

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
}
)

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