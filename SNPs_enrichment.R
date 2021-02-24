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

unique.Regul.ABC.tmp=unique.Regul.ABC
names(unique.Regul.ABC.tmp) = 1:length(unique.Regul.ABC.tmp)
unique.Regul.Rao.tmp=unique.Regul.Rao
names(unique.Regul.Rao.tmp) = 1:length(unique.Regul.Rao.tmp)
unique.Regul.DNAse.tmp=unique.Regul.DNAse
names(unique.Regul.DNAse.tmp) = 1:length(unique.Regul.DNAse.tmp)

sum(countOverlaps(unique.Regul.ABC,unique.Regul.Rao))
sum(countOverlaps(unique.Regul.Rao,unique.Regul.DNAse))


olaps_fromGRangesList = ssvOverlapIntervalSets(
  GenomicRanges::GRangesList(CTCF_in_10a_narrowPeak_grs))

library(ChIPpeakAnno)
getVennCounts(unique.Regul.ABC.tmp,unique.Regul.Rao.tmp,unique.Regul.DNAse.tmp, by="region", connectedPeaks = "keepAll")

res <- makeVennDiagram(Peaks=list(unique.Regul.ABC.tmp,unique.Regul.DNAse.tmp, unique.Regul.Rao.tmp),
                       NameOfPeaks=c("ABC", "DNAse", "Rao"), fill=c("#FF0000","#00FF00","#0000FF"),connectedPeaks="keepAll")

#Ensembles equivalents
regul.candidates.ABC = make.comparable.set(unique.Regul.ABC,method="ABC", element="regulatory", DNAse = peaks.NEU)
promoters.candidates.ABC = make.comparable.set(unique.Promoters.ABC,method="ABC", element="promoter")
regul.candidates.Rao = make.comparable.set(unique.Regul.Rao,method="Rao", element="regulatory", contact =positive.contact )
promoters.candidates.Rao = make.comparable.set(unique.Promoters.Rao,method="Rao", element="promoter", contact = positive.contact)
regul.candidates.DNAse = make.comparable.set(unique.Promoters.DNAse,method="DNAse", element="regulatory",contact = positive.contact, DNAse = peaks.NEU)
promoters.candidates.DNAse = make.comparable.set(unique.Promoters.DNAse,method="DNAse", element="promoter", contact = positive.contact, DNAse = peaks.NEU)


#Fichier de variants clumpe
SCZ3.all.clumped = read.table("../Genetic/PGC3_SCZ_wave3_public.clumped.v2.tsv", header=T, fill = T)
GRanges.snps.SCZ3.clumped = GRanges(seqnames=paste0("chr",SCZ3.all.clumped$CHR), ranges=IRanges(start=SCZ3.all.clumped$BP, end=SCZ3.all.clumped$BP, names=SCZ3.all.clumped$SNP), pval=SCZ3.all.clumped$P)

commom.variants = read.table("allchrs_filter.bim", header=F)
commom.variants$V1 = paste0("chr",commom.variants$V1)

GRanges.common = GRanges(seqnames=commom.variants$V1, ranges=IRanges(start=commom.variants$V4, end=commom.variants$V4, names = commom.variants$V2))

#Background proportion of common SNPs estimated by ldsc overlapping regulatory-elements
background.beh.ABC = 0.016
background.beh.cand.ABC = 0.056
background.beh.Rao = 0.026
background.beh.cand.Rao = 0.037
background.beh.DNAse = 0.0029
background.beh.cand.DNAse = 0.071


#Differents seuils de significativite consideres
threshold = c(0.95,0.5,0.10,0.05,0.025,0.01,0.001,0.0001,0.00001,0.000001, 0.0000001)

subset.SNPs.by.signi = lapply(threshold, function(x) GRanges.snps.SCZ3.clumped[GRanges.snps.SCZ3.clumped$pval<=x])
subset.SNPs.by.nsigni = lapply(threshold, function(x) GRanges.snps.SCZ3.clumped[GRanges.snps.SCZ3.clumped$pval>x])

#Are DNAse peaks enrichded in SNPs compared to the resto of genome? 
peakAnno <- annotatePeak(peaks.NEU, tssRegion=c(-250, 250),
                             TxDb=txdb, annoDb="org.Hs.eg.db", level="gene",sameStrand=F)

GRange.peakAnno = as.GRanges(peakAnno)
GRange.peakAnno$annotation_simplified = ifelse(startsWith(GRange.peakAnno$annotation, "Exon"), "Exon",ifelse(startsWith(GRange.peakAnno$annotation, "Intron"), "Intron",GRange.peakAnno$annotation))

df.enrichment.peaks = data.frame("p.value" = matrix(NA, ncol=1, nrow=11),"OR" = matrix(NA, ncol=1, nrow=11),"CI_lower" = matrix(NA, ncol=1, nrow=11),"CI_upper" = matrix(NA, ncol=1, nrow=11))

for(i in 1:11) {
  obj.f = fisher.test(matrix(c(sum(countOverlaps(GRange.peakAnno, subset.SNPs.by.signi[[i]])),length(subset.SNPs.by.signi[[i]]) - sum(countOverlaps(GRange.peakAnno, subset.SNPs.by.signi[[i]])),
                       sum(countOverlaps(GRange.peakAnno, subset.SNPs.by.nsigni[[i]])),length(subset.SNPs.by.nsigni[[i]])-sum(countOverlaps(GRange.peakAnno, subset.SNPs.by.nsigni[[i]]))),ncol=2), alternative = "two.sided")
  df.enrichment.peaks[i,] = c(threshold[i],obj.f$estimate, obj.f$conf.int[1], obj.f$conf.int[2])
}

ggplot(df.enrichment.peaks, aes(x=OR, y=factor(p.value)))+ geom_pointrange(aes(xmin=CI_lower, xmax=CI_upper),position=position_dodge(width=0.5)) +  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed")+
  ylab("P.value")

#Comparison of Proportion of fine-mapped SNPs to all common variants
#for different significance thresholds
#Update this code part with new results from ldsc 

fold.enrichment.ABC = lapply(1:11, function(x) { 
  signi.ABC = sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.signi[[x]]));nsigni.ABC=sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.nsigni[[x]]))
  signi.cand.ABC = sum(countOverlaps(regul.candidates.ABC, subset.SNPs.by.signi[[x]]));nsigni.cand.ABC=sum(countOverlaps(regul.candidates.ABC, subset.SNPs.by.nsigni[[x]]))
  
  a = signi.ABC/(signi.ABC+nsigni.ABC)
  b = signi.cand.ABC/(signi.cand.ABC+nsigni.cand.ABC)
  
  list("fold.enrichment"= (a/background.beh.ABC)/(b/background.beh.cand.ABC), "binom.test" = binom.test(c(signi.ABC, signi.ABC+nsigni.ABC), p=background.beh.ABC))
  })

fold.enrichment.Rao = lapply(1:11, function(x) { 
  
  signi.Rao = sum(countOverlaps(unique.Regul.Rao, subset.SNPs.by.signi[[x]]));nsigni.Rao=sum(countOverlaps(unique.Regul.Rao, subset.SNPs.by.nsigni[[x]]))
  signi.cand.Rao = sum(countOverlaps(regul.candidates.Rao, subset.SNPs.by.signi[[x]]));nsigni.cand.Rao=sum(countOverlaps(regul.candidates.Rao, subset.SNPs.by.nsigni[[x]]))
  
  a = signi.Rao/(signi.Rao+nsigni.Rao)
  b = signi.cand.Rao/(signi.cand.Rao+nsigni.cand.Rao)
  
  list("fold.enrichment"= (a/background.beh.Rao)/(b/background.beh.cand.Rao), "binom.test" = binom.test(c(signi.Rao, signi.Rao+nsigni.Rao), p=background.beh.Rao))
})

fold.enrichment.DNAse = lapply(1:11, function(x) { 
  signi.DNAse = sum(countOverlaps(unique.Regul.DNAse, subset.SNPs.by.signi[[x]]));nsigni.DNAse=sum(countOverlaps(unique.Regul.DNAse, subset.SNPs.by.nsigni[[x]]))
  signi.cand.DNAse = sum(countOverlaps(regul.candidates.DNAse, subset.SNPs.by.signi[[x]]));nsigni.cand.DNAse=sum(countOverlaps(regul.candidates.DNAse, subset.SNPs.by.nsigni[[x]]))
  
  a = signi.DNAse/(signi.DNAse+nsigni.DNAse)
  b = signi.cand.DNAse/(signi.cand.DNAse+nsigni.cand.DNAse)
  
  list("fold.enrichment"= (a/background.beh.DNAse)/(b/background.beh.cand.DNAse), "binom.test" = binom.test(c(signi.DNAse, signi.DNAse+nsigni.DNAse), p=background.beh.DNAse))
})

fold.enrichment.ABC[[6]]
fold.enrichment.Rao[[6]]
fold.enrichment.DNAse[[6]]

threshold[7]
enrichment.ABC = sapply(1:11, function(x) { 
  signi.ABC = sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.signi[[x]]));nsigni.ABC=sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.nsigni[[x]]))
  
  a = signi.ABC/(signi.ABC+nsigni.ABC)
  
  
  (a/background.beh.ABC)
  })


enrichment.cand.ABC = sapply(1:11, function(x) { 
  signi.cand.ABC = sum(countOverlaps(regul.candidates.ABC, subset.SNPs.by.signi[[x]]));nsigni.cand.ABC=sum(countOverlaps(regul.candidates.ABC, subset.SNPs.by.nsigni[[x]]))
  b = signi.cand.ABC/(signi.cand.ABC+nsigni.cand.ABC)
  
  b/background.beh.cand.ABC

})

enrichment.Rao = sapply(1:11, function(x) { 
  signi.Rao = sum(countOverlaps(unique.Regul.Rao, subset.SNPs.by.signi[[x]]));nsigni.Rao=sum(countOverlaps(unique.Regul.Rao, subset.SNPs.by.nsigni[[x]]))
  
  a = signi.Rao/(signi.Rao+nsigni.Rao)
  
  
  (a/background.beh.Rao)
})


enrichment.cand.Rao = sapply(1:11, function(x) { 
  signi.cand.Rao = sum(countOverlaps(regul.candidates.Rao, subset.SNPs.by.signi[[x]]));nsigni.cand.Rao=sum(countOverlaps(regul.candidates.Rao, subset.SNPs.by.nsigni[[x]]))
  b = signi.cand.Rao/(signi.cand.Rao+nsigni.cand.Rao)
  
  b/background.beh.cand.Rao
  
})

enrichment.DNAse = sapply(1:11, function(x) { 
  signi.DNAse = sum(countOverlaps(unique.Regul.DNAse, subset.SNPs.by.signi[[x]]));nsigni.DNAse=sum(countOverlaps(unique.Regul.DNAse, subset.SNPs.by.nsigni[[x]]))
  
  a = signi.DNAse/(signi.DNAse+nsigni.DNAse)
  
  
  (a/background.beh.DNAse)
})


enrichment.cand.DNAse = sapply(1:11, function(x) { 
  signi.cand.DNAse = sum(countOverlaps(regul.candidates.DNAse, subset.SNPs.by.signi[[x]]));nsigni.cand.DNAse=sum(countOverlaps(regul.candidates.DNAse, subset.SNPs.by.nsigni[[x]]))
  b = signi.cand.DNAse/(signi.cand.DNAse+nsigni.cand.DNAse)
  
  b/background.beh.cand.DNAse
  
})

df.fold.enrichment.ABC = cbind(threshold, data.frame(enrichment.cand.ABC),data.frame(enrichment.ABC))

ggplot(df.fold.enrichment.ABC, aes(x=factor(threshold), y=enrichment.ABC/enrichment.cand.ABC))+geom_point()+
  xlab("Pvalue") + ylab("Fold Enrichment") + geom_point(data=df.fold.enrichment.ABC[df.fold.enrichment.ABC$threshold=="0.001",], 
                                                        aes(x=factor(threshold),y=enrichment.ABC/enrichment.cand.ABC), 
                                                        color='red',
                                                        size=3)
df.fold.enrichment.Rao = cbind(threshold, data.frame(enrichment.cand.Rao),data.frame(enrichment.Rao))
ggplot(df.fold.enrichment.Rao, aes(x=factor(threshold), y=enrichment.Rao/enrichment.cand.Rao))+geom_point()+
  geom_point(data=df.fold.enrichment.Rao[df.fold.enrichment.Rao$threshold=="0.001",], 
             aes(x=factor(threshold),y=enrichment.Rao/enrichment.cand.Rao), 
             color='red',
             size=3)+xlab("Pvalue") + ylab("Fold Enrichment")

df.fold.enrichment.DNAse = cbind(threshold, data.frame(enrichment.cand.DNAse),data.frame(enrichment.DNAse))
ggplot(df.fold.enrichment.DNAse, aes(x=factor(threshold), y=enrichment.DNAse/enrichment.cand.DNAse))+geom_point()+
  geom_point(data=df.fold.enrichment.DNAse[df.fold.enrichment.DNAse$threshold=="0.001",], 
             aes(x=factor(threshold),y=enrichment.DNAse/enrichment.cand.DNAse), 
             color='red',
             size=3)+xlab("Pvalue") + ylab("Fold Enrichment")

signi.overlap.threshold.ABC = sapply(seq(.001,1,0.001), function(x){
  signi = GRanges.snps.SCZ3.clumped[GRanges.snps.SCZ3.clumped$pval<=x]
  nsigni = GRanges.snps.SCZ3.clumped[GRanges.snps.SCZ3.clumped$pval>x]
  
  n1 = sum(countOverlaps(unique.Regul.ABC, signi))
  n0 = sum(countOverlaps(unique.Regul.ABC, nsigni))
  
  n1/(n0+n1)
})
plot(seq(.01,1,0.01),signi.overlap.threshold.ABC)

#here I remove duplicated between regulatory and promoters for each annotation methods

library(epitools)
substract.overlap.PromRegul.ABC = findOverlaps(unique.Regul.ABC, unique.Promoters.ABC)
unique.Regul.ABC.subset = unique.Regul.ABC[-queryHits(substract.overlap.PromRegul.ABC)]

substract.overlap.PromRegul.Rao = findOverlaps(unique.Regul.Rao, unique.Promoters.Rao)
unique.Regul.Rao.subset = unique.Regul.Rao[-queryHits(substract.overlap.PromRegul.Rao)]

substract.overlap.PromRegul.DNAse = findOverlaps(unique.Regul.DNAse, unique.Promoters.DNAse)
unique.Regul.DNAse.subset = unique.Regul.DNAse[-queryHits(substract.overlap.PromRegul.DNAse)]

substract.overlap.PromRegul.candidates.ABC = findOverlaps(regul.candidates.ABC, promoters.candidates.ABC)
regul.candidates.ABC.subset = regul.candidates.ABC[-queryHits(substract.overlap.PromRegul.candidates.ABC)]

substract.overlap.PromRegul.candidates.Rao = findOverlaps(regul.candidates.Rao, promoters.candidates.Rao)
regul.candidates.Rao.subset = regul.candidates.Rao[-queryHits(substract.overlap.PromRegul.candidates.Rao)]

substract.overlap.PromRegul.candidates.DNAse = findOverlaps(regul.candidates.DNAse, promoters.candidates.DNAse)
regul.candidates.DNAse.subset = regul.candidates.DNAse[-queryHits(substract.overlap.PromRegul.candidates.DNAse)]


tab.enrichment.ABC = lapply(1:11, function(i)
  {
  signi.CRNs.ABC = sum(countOverlaps(unique.Regul.ABC.subset, subset.SNPs.by.signi[[i]])) + sum(countOverlaps(unique.Promoters.ABC, subset.SNPs.by.signi[[i]]))
  nsigni.CRNs.ABC = sum(countOverlaps(unique.Regul.ABC.subset, subset.SNPs.by.nsigni[[i]])) + sum(countOverlaps(unique.Promoters.ABC, subset.SNPs.by.nsigni[[i]]))
  
  signi.outCRNs.ABC = sum(countOverlaps(regul.candidates.ABC.subset,subset.SNPs.by.signi[[i]])) + sum(countOverlaps(promoters.candidates.ABC, subset.SNPs.by.signi[[i]]))
  nsigni.outCRNs.ABC = sum(countOverlaps(regul.candidates.ABC.subset,subset.SNPs.by.nsigni[[i]])) + sum(countOverlaps(promoters.candidates.ABC, subset.SNPs.by.nsigni[[i]]))
  
  
  signi.genome.rest = length(subset.SNPs.by.signi[[i]]) - signi.CRNs.ABC - signi.outCRNs.ABC
  nsigni.genome.rest = length(subset.SNPs.by.nsigni[[i]]) - nsigni.CRNs.ABC - nsigni.outCRNs.ABC
  
  epitab(matrix(c(signi.CRNs.ABC,signi.outCRNs.ABC,signi.genome.rest,nsigni.CRNs.ABC,nsigni.outCRNs.ABC,nsigni.genome.rest),ncol=2), method="oddsratio")$tab
  
})

prop.SNPs.CRNs.ABC = lapply(1:11, function(i)
{
  signi.CRNs.ABC = sum(countOverlaps(unique.Regul.ABC.subset, subset.SNPs.by.signi[[i]])) + sum(countOverlaps(unique.Promoters.ABC, subset.SNPs.by.signi[[i]]))
  #nsigni.CRNs.ABC = sum(countOverlaps(unique.Regul.ABC.subset, subset.SNPs.by.nsigni[[i]])) + sum(countOverlaps(unique.Promoters.ABC, subset.SNPs.by.nsigni[[i]]))
  d0 = data.frame("p.value"= threshold[i], "Element"="ABC", "NSigni" = signi.CRNs.ABC)
  signi.outCRNs.ABC = sum(countOverlaps(regul.candidates.ABC.subset,subset.SNPs.by.signi[[i]])) + sum(countOverlaps(promoters.candidates.ABC, subset.SNPs.by.signi[[i]]))
  #nsigni.outCRNs.ABC = sum(countOverlaps(regul.candidates.ABC.subset,subset.SNPs.by.nsigni[[i]])) + sum(countOverlaps(promoters.candidates.ABC, subset.SNPs.by.nsigni[[i]]))
  d1 = data.frame("p.value"= threshold[i], "Element"="Candidates", "NSigni" = signi.outCRNs.ABC)
  
  signi.genome.rest = length(subset.SNPs.by.signi[[i]]) - signi.CRNs.ABC - signi.outCRNs.ABC
  #nsigni.genome.rest = length(subset.SNPs.by.nsigni[[i]]) - nsigni.CRNs.ABC - nsigni.outCRNs.ABC
  d2 = data.frame("p.value"= threshold[i], "Element"="Genome", "NSigni" = signi.genome.rest)
  
  #epitab(matrix(c(signi.CRNs.ABC,signi.outCRNs.ABC,signi.genome.rest,nsigni.CRNs.ABC,nsigni.outCRNs.ABC,nsigni.genome.rest),ncol=2), method="oddsratio")$tab
  rbind(d0,d1,d2)
})

SNPs.prop.CRNs.ABC = do.call("rbind", prop.SNPs.CRNs.ABC)
  
ggplot(SNPs.prop.CRNs.ABC, aes(fill=Element, y=NSigni, x=as.factor(p.value))) + 
  geom_bar(position="fill", stat="identity")+scale_fill_brewer(palette="Spectral")+ xlab("P.value") + ylab("Fraction of Significant SNPs")


tab.enrichment.Rao = lapply(1:11, function(i){
  signi.CRNs.Rao = sum(countOverlaps(unique.Regul.Rao.subset, subset.SNPs.by.signi[[i]])) + sum(countOverlaps(unique.Promoters.Rao, subset.SNPs.by.signi[[i]]))
  nsigni.CRNs.Rao = sum(countOverlaps(unique.Regul.Rao.subset, subset.SNPs.by.nsigni[[i]])) + sum(countOverlaps(unique.Promoters.Rao, subset.SNPs.by.nsigni[[i]]))
  
  signi.outCRNs.Rao = sum(countOverlaps(regul.candidates.Rao.subset,subset.SNPs.by.signi[[i]])) + sum(countOverlaps(promoters.candidates.Rao, subset.SNPs.by.signi[[i]]))
  nsigni.outCRNs.Rao = sum(countOverlaps(regul.candidates.Rao.subset,subset.SNPs.by.nsigni[[i]])) + sum(countOverlaps(promoters.candidates.Rao, subset.SNPs.by.nsigni[[i]]))
  
  
  signi.genome.rest = length(subset.SNPs.by.signi[[i]]) - signi.CRNs.Rao - signi.outCRNs.Rao
  nsigni.genome.rest = length(subset.SNPs.by.nsigni[[i]]) - nsigni.CRNs.Rao - nsigni.outCRNs.Rao
  
  epitab(matrix(c(signi.CRNs.Rao,signi.outCRNs.Rao,signi.genome.rest,nsigni.CRNs.Rao,nsigni.outCRNs.Rao,nsigni.genome.rest),ncol=2), method="oddsratio")$tab
  
})

prop.SNPs.CRNs.Rao = lapply(1:11, function(i)
{
  signi.CRNs.Rao = sum(countOverlaps(unique.Regul.Rao.subset, subset.SNPs.by.signi[[i]])) + sum(countOverlaps(unique.Promoters.Rao, subset.SNPs.by.signi[[i]]))
  #nsigni.CRNs.Rao = sum(countOverlaps(unique.Regul.Rao.subset, subset.SNPs.by.nsigni[[i]])) + sum(countOverlaps(unique.Promoters.Rao, subset.SNPs.by.nsigni[[i]]))
  d0 = data.frame("p.value"= threshold[i], "Element"="Rao", "NSigni" = signi.CRNs.Rao)
  signi.outCRNs.Rao = sum(countOverlaps(regul.candidates.Rao.subset,subset.SNPs.by.signi[[i]])) + sum(countOverlaps(promoters.candidates.Rao, subset.SNPs.by.signi[[i]]))
  #nsigni.outCRNs.Rao = sum(countOverlaps(regul.candidates.Rao.subset,subset.SNPs.by.nsigni[[i]])) + sum(countOverlaps(promoters.candidates.Rao, subset.SNPs.by.nsigni[[i]]))
  d1 = data.frame("p.value"= threshold[i], "Element"="Candidates", "NSigni" = signi.outCRNs.Rao)
  
  signi.genome.rest = length(subset.SNPs.by.signi[[i]]) - signi.CRNs.Rao - signi.outCRNs.Rao
  #nsigni.genome.rest = length(subset.SNPs.by.nsigni[[i]]) - nsigni.CRNs.Rao - nsigni.outCRNs.Rao
  d2 = data.frame("p.value"= threshold[i], "Element"="Genome", "NSigni" = signi.genome.rest)
  
  #epitab(matrix(c(signi.CRNs.Rao,signi.outCRNs.Rao,signi.genome.rest,nsigni.CRNs.Rao,nsigni.outCRNs.Rao,nsigni.genome.rest),ncol=2), method="oddsratio")$tab
  rbind(d0,d1,d2)
})

SNPs.prop.CRNs.Rao = do.call("rbind", prop.SNPs.CRNs.Rao)

ggplot(SNPs.prop.CRNs.Rao, aes(fill=Element, y=NSigni, x=as.factor(p.value))) + 
  geom_bar(position="fill", stat="identity")+scale_fill_brewer(palette="Spectral")+ xlab("P.value") + ylab("Fraction of Significant SNPs")


tab.enrichment.DNAse = lapply(1:11, function(i) {
  signi.CRNs.DNAse = sum(countOverlaps(unique.Regul.DNAse.subset, subset.SNPs.by.signi[[i]])) + sum(countOverlaps(unique.Promoters.DNAse, subset.SNPs.by.signi[[i]]))
  nsigni.CRNs.DNAse = sum(countOverlaps(unique.Regul.DNAse.subset, subset.SNPs.by.nsigni[[i]])) + sum(countOverlaps(unique.Promoters.DNAse, subset.SNPs.by.nsigni[[i]]))
  
  signi.outCRNs.DNAse = sum(countOverlaps(regul.candidates.DNAse.subset,subset.SNPs.by.signi[[i]])) + sum(countOverlaps(promoters.candidates.DNAse, subset.SNPs.by.signi[[i]]))
  nsigni.outCRNs.DNAse = sum(countOverlaps(regul.candidates.DNAse.subset,subset.SNPs.by.nsigni[[i]])) + sum(countOverlaps(promoters.candidates.DNAse, subset.SNPs.by.nsigni[[i]]))
  
  
  signi.genome.rest = length(subset.SNPs.by.signi[[i]]) - signi.CRNs.DNAse - signi.outCRNs.DNAse
  nsigni.genome.rest = length(subset.SNPs.by.nsigni[[i]]) - nsigni.CRNs.DNAse - nsigni.outCRNs.DNAse
  
  epitab(matrix(c(signi.CRNs.DNAse,signi.outCRNs.DNAse,signi.genome.rest,nsigni.CRNs.DNAse,nsigni.outCRNs.DNAse,nsigni.genome.rest),ncol=2), method="oddsratio")$tab
  
})

prop.SNPs.CRNs.DNAse = lapply(1:11, function(i)
{
  signi.CRNs.DNAse = sum(countOverlaps(unique.Regul.DNAse.subset, subset.SNPs.by.signi[[i]])) + sum(countOverlaps(unique.Promoters.DNAse, subset.SNPs.by.signi[[i]]))
  #nsigni.CRNs.DNAse = sum(countOverlaps(unique.Regul.DNAse.subset, subset.SNPs.by.nsigni[[i]])) + sum(countOverlaps(unique.Promoters.DNAse, subset.SNPs.by.nsigni[[i]]))
  d0 = data.frame("p.value"= threshold[i], "Element"="DNAse", "NSigni" = signi.CRNs.DNAse)
  signi.outCRNs.DNAse = sum(countOverlaps(regul.candidates.DNAse.subset,subset.SNPs.by.signi[[i]])) + sum(countOverlaps(promoters.candidates.DNAse, subset.SNPs.by.signi[[i]]))
  #nsigni.outCRNs.DNAse = sum(countOverlaps(regul.candidates.DNAse.subset,subset.SNPs.by.nsigni[[i]])) + sum(countOverlaps(promoters.candidates.DNAse, subset.SNPs.by.nsigni[[i]]))
  d1 = data.frame("p.value"= threshold[i], "Element"="Candidates", "NSigni" = signi.outCRNs.DNAse)
  
  signi.genome.rest = length(subset.SNPs.by.signi[[i]]) - signi.CRNs.DNAse - signi.outCRNs.DNAse
  #nsigni.genome.rest = length(subset.SNPs.by.nsigni[[i]]) - nsigni.CRNs.DNAse - nsigni.outCRNs.DNAse
  d2 = data.frame("p.value"= threshold[i], "Element"="Genome", "NSigni" = signi.genome.rest)
  
  #epitab(matrix(c(signi.CRNs.DNAse,signi.outCRNs.DNAse,signi.genome.rest,nsigni.CRNs.DNAse,nsigni.outCRNs.DNAse,nsigni.genome.rest),ncol=2), method="oddsratio")$tab
  rbind(d0,d1,d2)
})

SNPs.prop.CRNs.DNAse = do.call("rbind", prop.SNPs.CRNs.DNAse)

ggplot(SNPs.prop.CRNs.DNAse, aes(fill=Element, y=NSigni, x=as.factor(p.value))) + 
  geom_bar(position="fill", stat="identity")+scale_fill_brewer(palette="Spectral") + xlab("P.value") + ylab("Fraction of Significant SNPs")


#Enrichment for all annotation methods at 1e-4 threshold level 
tab.enrichment.ABC[[8]]
tab.enrichment.Rao[[8]]
tab.enrichment.DNAse[[8]]


list.graphs.enrichment.CRNs = list()
for(i in 1:11){
  ABC.enrichment.CRNs = data.frame(tab.enrichment.ABC[[i]][2:3,c("oddsratio", "lower", "upper")])
  ABC.enrichment.CRNs$method = "ABC"
  ABC.enrichment.CRNs$element = c("Candidates", "Genome")
  Rao.enrichment.CRNs = data.frame(tab.enrichment.Rao[[i]][2:3,c("oddsratio", "lower", "upper")])
  Rao.enrichment.CRNs$method = "Rao"
  Rao.enrichment.CRNs$element = c("Candidates", "Genome")
  DNAse.enrichment.CRNs = data.frame(tab.enrichment.DNAse[[i]][2:3,c("oddsratio", "lower", "upper")])
  DNAse.enrichment.CRNs$method = "DNAse"
  DNAse.enrichment.CRNs$element = c("Candidates", "Genome")
  
  all.enrichment.CRNs = rbind(ABC.enrichment.CRNs, Rao.enrichment.CRNs, DNAse.enrichment.CRNs)

  list.graphs.enrichment.CRNs[[i]] = ggplot(all.enrichment.CRNs, aes(x=oddsratio,y=element, color = method )) + geom_pointrange(aes(xmin=lower, xmax=upper),position=position_dodge(width=0.5)) +  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed")+
    xlab("OR") + ylab("Comparison Elements")
}
list.graphs.enrichment.CRNs[[11]]



txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peakAnno.ABC <- annotatePeak(regul.candidates.ABC, tssRegion=c(-250, 250),
                         TxDb=txdb, annoDb="org.Hs.eg.db", level="gene",sameStrand=F)

GRange.peakAnno.ABC = as.GRanges(peakAnno.ABC)
GRange.peakAnno.ABC$annotation_simplified = ifelse(startsWith(GRange.peakAnno.ABC$annotation, "Exon"), "Exon",ifelse(startsWith(GRange.peakAnno.ABC$annotation, "Intron"), "Intron",GRange.peakAnno.ABC$annotation))

peakAnno.DNAse <- annotatePeak(regul.candidates.DNAse, tssRegion=c(-250, 250),
                             TxDb=txdb, annoDb="org.Hs.eg.db", level="gene",sameStrand=F)

GRange.peakAnno.DNAse = as.GRanges(peakAnno.DNAse)
GRange.peakAnno.DNAse$annotation_simplified = ifelse(startsWith(GRange.peakAnno.DNAse$annotation, "Exon"), "Exon",ifelse(startsWith(GRange.peakAnno.DNAse$annotation, "Intron"), "Intron",GRange.peakAnno.DNAse$annotation))



SNP.prop.cat.ABC = lapply(threshold, function(x) {
  
  signi = GRanges.snps.SCZ3.clumped[GRanges.snps.SCZ3.clumped$pval<=x]
  nsigni = GRanges.snps.SCZ3.clumped[GRanges.snps.SCZ3.clumped$pval>x]
  
  GRange.peakAnno.ABC$N_signiSNPs  = countOverlaps(GRange.peakAnno.ABC, signi)
  GRange.peakAnno.ABC$N_nsigniSNPs = countOverlaps(GRange.peakAnno.ABC, nsigni)
  
  
  merge.peaks.ABC = merge(aggregate(N_signiSNPs~annotation_simplified, data.frame(GRange.peakAnno.ABC), sum),aggregate(N_nsigniSNPs~annotation_simplified, data.frame(GRange.peakAnno.ABC), sum))
  merge.peaks.ABC[nrow(merge.peaks.ABC)+1,] = c("ABC", sum(countOverlaps(unique.Regul.ABC, signi)),sum(countOverlaps(unique.Regul.ABC, nsigni)))
  merge.peaks.ABC = merge.peaks.ABC[c(nrow(merge.peaks.ABC),1,2,3,4,5,6,7,8,9),]
  #rownames(merge.peaks.ABC) = merge.peaks.ABC$annotation_simplified 
  #merge.peaks.ABC$annotation_simplified = NULL
  merge.peaks.ABC$N_signiSNPs = as.numeric(merge.peaks.ABC$N_signiSNPs)
  merge.peaks.ABC$N_nsigniSNPs = as.numeric(merge.peaks.ABC$N_nsigniSNPs)
  merge.peaks.ABC$fraction_signiSNPs = merge.peaks.ABC$N_signiSNPs/sum(merge.peaks.ABC$N_signiSNPs)
  merge.peaks.ABC$p.value = x
  merge.peaks.ABC
})

SNP.prop.cat.ABC.all = do.call("rbind", SNP.prop.cat.ABC)

SNP.prop.cat.Rao = lapply(threshold, function(x) {
  
  signi = GRanges.snps.SCZ3.clumped[GRanges.snps.SCZ3.clumped$pval<=x]
  nsigni = GRanges.snps.SCZ3.clumped[GRanges.snps.SCZ3.clumped$pval>x]
  
  unique.Regul.Rao$N_signiSNPs  = countOverlaps(unique.Regul.Rao, signi)
  unique.Regul.Rao$N_nsigniSNPs = countOverlaps(unique.Regul.Rao, nsigni)
  
  other.signi = length(signi) - sum(unique.Regul.Rao$N_signiSNPs)
  other.nsigni = length(nsigni) - sum(unique.Regul.Rao$N_nsigniSNPs)
  
  rbind(data.frame("p.value"=x,"Signi" = sum(unique.Regul.Rao$N_signiSNPs), "nSigni.Rao"=sum(unique.Regul.Rao$N_nsigniSNPs), "annotation_simplified" = "Rao"),
        data.frame("p.value"=x,"Signi" = other.signi, "nSigni.Rao"=other.nsigni, "annotation_simplified" = "Other"))
})
SNP.prop.cat.Rao.all = do.call("rbind", SNP.prop.cat.Rao)

SNP.prop.cat.DNAse = lapply(threshold, function(x) {
  
  signi = GRanges.snps.SCZ3.clumped[GRanges.snps.SCZ3.clumped$pval<=x]
  nsigni = GRanges.snps.SCZ3.clumped[GRanges.snps.SCZ3.clumped$pval>x]
  
  GRange.peakAnno.DNAse$N_signiSNPs  = countOverlaps(GRange.peakAnno.DNAse, signi)
  GRange.peakAnno.DNAse$N_nsigniSNPs = countOverlaps(GRange.peakAnno.DNAse, nsigni)
  
  
  merge.peaks.DNAse = merge(aggregate(N_signiSNPs~annotation_simplified, data.frame(GRange.peakAnno.DNAse), sum),aggregate(N_nsigniSNPs~annotation_simplified, data.frame(GRange.peakAnno.DNAse), sum))
  merge.peaks.DNAse[nrow(merge.peaks.DNAse)+1,] = c("DNAse", sum(countOverlaps(unique.Regul.DNAse, signi)),sum(countOverlaps(unique.Regul.DNAse, nsigni)))
  merge.peaks.DNAse = merge.peaks.DNAse[c(nrow(merge.peaks.DNAse),1,2,3,4,5,6,7,8,9),]
  #rownames(merge.peaks.DNAse) = merge.peaks.DNAse$annotation_simplified 
  #merge.peaks.DNAse$annotation_simplified = NULL
  merge.peaks.DNAse$N_signiSNPs = as.numeric(merge.peaks.DNAse$N_signiSNPs)
  merge.peaks.DNAse$N_nsigniSNPs = as.numeric(merge.peaks.DNAse$N_nsigniSNPs)
  merge.peaks.DNAse$fraction_signiSNPs = merge.peaks.DNAse$N_signiSNPs/sum(merge.peaks.DNAse$N_signiSNPs)
  merge.peaks.DNAse$p.value = x
  merge.peaks.DNAse
})
SNP.prop.cat.DNAse.all = do.call("rbind", SNP.prop.cat.DNAse)

ggplot(SNP.prop.cat.ABC.all, aes(fill=annotation_simplified, y=N_signiSNPs, x=as.factor(p.value))) + 
  geom_bar(position="fill", stat="identity")+scale_fill_brewer(palette="Spectral")+xlab("Pvalue")+ylab("% Significant SNPs")

ggplot(SNP.prop.cat.Rao.all, aes(fill=annotation_simplified, y=Signi, x=as.factor(p.value))) + 
  geom_bar(position="fill", stat="identity")+scale_fill_brewer(palette="Spectral")

ggplot(SNP.prop.cat.DNAse.all, aes(fill=annotation_simplified, y=N_signiSNPs, x=as.factor(p.value))) + 
  geom_bar(position="fill", stat="identity")+scale_fill_brewer(palette="Spectral")+xlab("Pvalue")+ylab("% Significant SNPs")


ggplot(SNP.prop.cat.ABC.all[SNP.prop.cat.ABC.all$annotation_simplified=="ABC",], aes(x=p.value, y = fraction_signiSNPs))+geom_line(size=0.75)+
  scale_x_reverse()

ggplot(SNP.prop.cat.Rao.all[SNP.prop.cat.Rao.all$annotation_simplified=="Rao",], aes(x=p.value, y = Signi/sum(Signi)))+geom_line(size=0.75)+
  scale_x_reverse()

ggplot(SNP.prop.cat.DNAse.all[SNP.prop.cat.DNAse.all$annotation_simplified=="DNAse",], aes(x=p.value, y = fraction_signiSNPs))+geom_line(size=0.75)+
  scale_x_reverse()


enrichment.all.categories.ABC = lapply(1:11, function(x) {
GRange.peakAnno.ABC$N_signiSNPs  = countOverlaps(GRange.peakAnno.ABC, subset.SNPs.by.signi[[x]])
GRange.peakAnno.ABC$N_nsigniSNPs = countOverlaps(GRange.peakAnno.ABC, subset.SNPs.by.nsigni[[x]])


merge.peaks.ABC = merge(aggregate(N_signiSNPs~annotation_simplified, data.frame(GRange.peakAnno.ABC), sum),aggregate(N_nsigniSNPs~annotation_simplified, data.frame(GRange.peakAnno.ABC), sum))
merge.peaks.ABC[nrow(merge.peaks.ABC)+1,] = c("ABC", sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.signi[[x]])),sum(countOverlaps(unique.Regul.ABC, subset.SNPs.by.nsigni[[x]])))
merge.peaks.ABC = merge.peaks.ABC[c(nrow(merge.peaks.ABC),1,2,3,4,5,6,7,8,9),]
rownames(merge.peaks.ABC) = merge.peaks.ABC$annotation_simplified 
merge.peaks.ABC$annotation_simplified = NULL
epitab(data.matrix(merge.peaks.ABC), method="oddsratio")$tab})


peakAnno.DNAse <- annotatePeak(regul.candidates.DNAse, tssRegion=c(-250, 250),
                             TxDb=txdb, annoDb="org.Hs.eg.db", level="gene",sameStrand=F)

GRange.peakAnno.DNAse = as.GRanges(peakAnno.DNAse)
GRange.peakAnno.DNAse$annotation_simplified = ifelse(startsWith(GRange.peakAnno.DNAse$annotation, "Exon"), "Exon",ifelse(startsWith(GRange.peakAnno.DNAse$annotation, "Intron"), "Intron",GRange.peakAnno.DNAse$annotation))



enrichment.all.categories.DNAse = lapply(1:9, function(x) {
  GRange.peakAnno.DNAse$N_signiSNPs  = countOverlaps(GRange.peakAnno.DNAse, subset.SNPs.by.signi[[x]])
  GRange.peakAnno.DNAse$N_nsigniSNPs = countOverlaps(GRange.peakAnno.DNAse, subset.SNPs.by.nsigni[[x]])
  
  
  merge.peaks.DNAse = merge(aggregate(N_signiSNPs~annotation_simplified, data.frame(GRange.peakAnno.DNAse), sum),aggregate(N_nsigniSNPs~annotation_simplified, data.frame(GRange.peakAnno.DNAse), sum))
  merge.peaks.DNAse[nrow(merge.peaks.DNAse)+1,] = c("DNAse", sum(countOverlaps(unique.Regul.DNAse, subset.SNPs.by.signi[[x]])),sum(countOverlaps(unique.Regul.DNAse, subset.SNPs.by.nsigni[[x]])))
  merge.peaks.DNAse = merge.peaks.DNAse[c(nrow(merge.peaks.DNAse),1,2,3,4,5,6,7,8,9),]
  rownames(merge.peaks.DNAse) = merge.peaks.DNAse$annotation_simplified 
  merge.peaks.DNAse$annotation_simplified = NULL
  epitab(data.matrix(merge.peaks.DNAse), method="oddsratio")$tab})


#Partition of heritability by functionnal position
splitpeakAnno.ABC = split(GRange.peakAnno.ABC, GRange.peakAnno.ABC$annotation_simplified)
names(splitpeakAnno.ABC)
for(i in names(splitpeakAnno.ABC)){
  
  export.G = data.frame(seqnames(splitpeakAnno.ABC[[i]]), start(splitpeakAnno.ABC[[i]]), end(splitpeakAnno.ABC[[i]]))
  write.table(export.G, paste0(i,"ABC.bed"), col.names = F, row.names = F, sep="\t")
}


splitpeakAnno.DNAse = split(GRange.peakAnno.DNAse, GRange.peakAnno.DNAse$annotation_simplified)

for(i in names(splitpeakAnno.DNAse)){
  
  export.G = data.frame(seqnames(splitpeakAnno.DNAse[[i]]), start(splitpeakAnno.DNAse[[i]]), end(splitpeakAnno.DNAse[[i]]))
  write.table(export.G, paste0(i,"DNAse.bed"), col.names = F, row.names = F, sep="\t")
}

#Ici pour chaque analyse d'enrichissement, on compare: 
#Les elements de regulation ou promoters avec les ensembles equivalents
#Les RRCs par rapport au reste du genome (on supprimme les elements presents dans les ensembles equivalents)
#Les sorties des differents boucles sont dans les elements RDS fournis dans le github


unique.Regul.ABC.subset = unique.Regul.ABC
regul.candidates.ABC.subset = regul.candidates.ABC
unique.Regul.Rao.subset = unique.Regul.Rao
regul.candidates.Rao.subset = regul.candidates.Rao
unique.Regul.DNAse.subset = unique.Regul.DNAse
regul.candidates.DNAse.subset = regul.candidates.DNAse

findOverlaps(unique.Regul.ABC.subset, unique.Promoters.ABC)
findOverlaps(regul.candidates.ABC.subset, promoters.candidates.ABC)
findOverlaps(unique.Regul.Rao.subset, unique.Promoters.Rao)
findOverlaps(regul.candidates.Rao.subset, promoters.candidates.Rao)
findOverlaps(unique.Regul.DNAse.subset, unique.Promoters.DNAse)
findOverlaps(regul.candidates.DNAse.subset, promoters.candidates.DNAse)

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