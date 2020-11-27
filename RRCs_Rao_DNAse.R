library(rlist)
hic.loops = read.table("/home/loic/Documents/HiC/allloops_HiCCUPS/enriched_pixels_10000.bedpe", header=T)
peaks.NEU = import("/home/loic/Documents/HiC/code/Hi-C_analysis/NEU_DNAse.macs2_peaks.narrowPeak.sorted.candidateRegions.bed", format="bed")

peaks.NEU = peaks.NEU[seqnames(peaks.NEU)%in%paste0("chr", 1:22)]

GRanges.bin1 = GRanges(seqnames = hic.loops$chr1, ranges=IRanges(start=hic.loops$x1, end=hic.loops$x2))
GRanges.bin2 = GRanges(seqnames = hic.loops$chr2, ranges=IRanges(start=hic.loops$y1, end=hic.loops$y2))

Pairs.bins = Pairs(GRanges.bin1, GRanges.bin2,fdrBL = hic.loops$fdrBL ,fdrDonut = hic.loops$fdrDonut,fdrV=hic.loops$fdrV,fdrH = hic.loops$fdrH)

names(Pairs.bins) = paste0("bin", 1:length(Pairs.bins))
names(peaks.NEU) = paste0("peak", 1:length(peaks.NEU))

positive.contact = Pairs.bins[mcols(Pairs.bins)$fdrBL <= 0.15&mcols(Pairs.bins)$fdrDonut <=0.15&mcols(Pairs.bins)$fdrV <= 0.15&mcols(Pairs.bins)$fdrH <=0.15]

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes.hg19 <- genes(txdb)
symbol <- select(org.Hs.eg.db,keys = genes.hg19$gene_id, columns = c("SYMBOL"), keytype = "ENTREZID")
genes.hg19$geneSymbol <- symbol$SYMBOL

promoters = promoters(genes.hg19)
promoters = promoters[seqnames(promoters)%in%paste0("chr",1:22)]

#Genes overlaping first bin
overlaps.Prom.Bin1 = findOverlaps(promoters, first(positive.contact))
promoters.Bin1 = promoters[queryHits(overlaps.Prom.Bin1)]

#Bins for which the first bin has a gene
regulatory.elements.bin2 = second(positive.contact)[subjectHits(overlaps.Prom.Bin1)]

Pairs.prom.regulatory.1 = Pairs(promoters.Bin1, regulatory.elements.bin2)
Pairs.prom.regulatory.1 = Pairs.prom.regulatory.1[!duplicated(Pairs.prom.regulatory.1)]

#Genes overlaping second bin
overlaps.Prom.Bin2 = findOverlaps(promoters, second(positive.contact))
promoters.Bin2 = promoters[queryHits(overlaps.Prom.Bin2)]

#Bins for which the second bin has a gene
regulatory.elements.bin1 = first(positive.contact)[subjectHits(overlaps.Prom.Bin2)]

Pairs.prom.regulatory.2 = Pairs(promoters.Bin2, regulatory.elements.bin1)
Pairs.prom.regulatory.2 = Pairs.prom.regulatory.2[!duplicated(Pairs.prom.regulatory.2)]

grList.Rao = GRangesList(c(first(Pairs.prom.regulatory.1), first(Pairs.prom.regulatory.2)), c(second(Pairs.prom.regulatory.1), second(Pairs.prom.regulatory.2)))
Pairs.prom.regulatory.Rao = Pairs(grList.Rao[[1]], grList.Rao[[2]])

summary(width(reduce(first(Pairs.prom.regulatory.Rao))))
summary(width(reduce(second(Pairs.prom.regulatory.Rao))))

dist.prom.regul.Rao = start(first(Pairs.prom.regulatory.Rao)) - start(second(Pairs.prom.regulatory.Rao))
summary(dist.prom.regul.Rao)
summary(abs(dist.prom.regul.Rao))

n.aval = sum(ifelse(start(first(Pairs.prom.regulatory.Rao))<start(second(Pairs.prom.regulatory.Rao)),1,0))
n.amont= sum(ifelse(start(first(Pairs.prom.regulatory.Rao))>start(second(Pairs.prom.regulatory.Rao)),1,0))

n.aval/(n.amont+n.aval)

############################################################################################################
first(Pairs.prom.regulatory.Rao) = annotate.3D.Features(GRanges.DI, first(Pairs.prom.regulatory.Rao), kind = "DI",aggregateFunction ="mean")
first(Pairs.prom.regulatory.Rao) = annotate.3D.Features(GRanges.FIREs, first(Pairs.prom.regulatory.Rao), kind = "FIRE",aggregateFunction ="mean")
first(Pairs.prom.regulatory.Rao) = annotate.3D.Features(GRanges.INS, first(Pairs.prom.regulatory.Rao), kind = "INS",aggregateFunction ="mean")

second(Pairs.prom.regulatory.Rao) = annotate.3D.Features(GRanges.DI, second(Pairs.prom.regulatory.Rao), kind = "DI",aggregateFunction ="mean")
second(Pairs.prom.regulatory.Rao) = annotate.3D.Features(GRanges.FIREs, second(Pairs.prom.regulatory.Rao), kind = "FIRE",aggregateFunction ="mean")
second(Pairs.prom.regulatory.Rao) = annotate.3D.Features(GRanges.INS, second(Pairs.prom.regulatory.Rao), kind = "INS",aggregateFunction ="mean")

first(Pairs.prom.regulatory.Rao) = annotate.Activity(CTCFActivity.NEU, first(Pairs.prom.regulatory.Rao), kind = "CTCF")
first(Pairs.prom.regulatory.Rao) = annotate.Activity(H3K27acActivity.NEU, first(Pairs.prom.regulatory.Rao), kind = "H3K27ac")
first(Pairs.prom.regulatory.Rao) = annotate.Activity(DNAseActivity.NEU, first(Pairs.prom.regulatory.Rao), kind = "DNAse")
first(Pairs.prom.regulatory.Rao) = annotate.Activity(H3K4me1Activity.NEU, first(Pairs.prom.regulatory.Rao), kind = "H3K4me1")
first(Pairs.prom.regulatory.Rao) = annotate.Activity(H3K4me3Activity.NEU, first(Pairs.prom.regulatory.Rao), kind = "H3K4me3")
first(Pairs.prom.regulatory.Rao) = annotate.Activity(p300Activity.NEU, first(Pairs.prom.regulatory.Rao), kind = "p300")
first(Pairs.prom.regulatory.Rao) = annotate.Activity(H3K27me3Activity.NEU, first(Pairs.prom.regulatory.Rao), kind = "H3K27me3")

second(Pairs.prom.regulatory.Rao) = annotate.Activity(CTCFActivity.NEU, second(Pairs.prom.regulatory.Rao), kind = "CTCF")
second(Pairs.prom.regulatory.Rao) = annotate.Activity(H3K27acActivity.NEU, second(Pairs.prom.regulatory.Rao), kind = "H3K27ac")
second(Pairs.prom.regulatory.Rao) = annotate.Activity(DNAseActivity.NEU, second(Pairs.prom.regulatory.Rao), kind = "DNAse")
second(Pairs.prom.regulatory.Rao) = annotate.Activity(H3K4me1Activity.NEU, second(Pairs.prom.regulatory.Rao), kind = "H3K4me1")
second(Pairs.prom.regulatory.Rao) = annotate.Activity(H3K4me3Activity.NEU, second(Pairs.prom.regulatory.Rao), kind = "H3K4me3")
second(Pairs.prom.regulatory.Rao) = annotate.Activity(p300Activity.NEU, second(Pairs.prom.regulatory.Rao), kind = "p300")
second(Pairs.prom.regulatory.Rao) = annotate.Activity(H3K27me3Activity.NEU, second(Pairs.prom.regulatory.Rao), kind = "H3K27me3")

forward.Rao = unique(first(Pairs.prom.regulatory.Rao))
reverse.Rao = unique(second(Pairs.prom.regulatory.Rao))

prom.Rao = data.frame(mcols(forward.Rao))
regul.Rao = data.frame(mcols(reverse.Rao))

gl.tmp = GRangesList(forward.Rao, reverse.Rao)
loci.Rao = unique(do.call("c", as(gl.tmp, "GRangesList")))

FeaturesPlusActivity.Rao = data.frame(mcols(loci.Rao))

#DNAse peaks on first and second bins
overlaps.peaks.Bin1 = findOverlaps(peaks.NEU, first(positive.contact))
overlaps.peaks.Bin2 = findOverlaps(peaks.NEU, second(positive.contact))

peaksToProm = merge(data.frame(overlaps.peaks.Bin1), data.frame(overlaps.Prom.Bin2), by="subjectHits")
Promtopeaks = merge(data.frame(overlaps.peaks.Bin2), data.frame(overlaps.Prom.Bin1), by="subjectHits")

Pairs.peak.prom.1 = Pairs(peaks.NEU[peaksToProm$queryHits.x],promoters[peaksToProm$queryHits.y])
Pairs.peak.prom.2 = Pairs(peaks.NEU[Promtopeaks$queryHits.x],promoters[Promtopeaks$queryHits.y])

gl.Peaks = GRangesList(c(first(Pairs.peak.prom.1), first(Pairs.peak.prom.2)), c(second(Pairs.peak.prom.1), second(Pairs.peak.prom.2)))
Pairs.prom.regulatory.DNAse = Pairs(gl.Peaks[[1]], gl.Peaks[[2]])

summary(width(reduce(first(Pairs.prom.regulatory.DNAse))))
summary(width(reduce(second(Pairs.prom.regulatory.DNAse))))

#########################################################################################################
first(Pairs.prom.regulatory.DNAse) = annotate.3D.Features(GRanges.DI, first(Pairs.prom.regulatory.DNAse), kind = "DI",aggregateFunction ="mean")
first(Pairs.prom.regulatory.DNAse) = annotate.3D.Features(GRanges.FIREs, first(Pairs.prom.regulatory.DNAse), kind = "FIRE",aggregateFunction ="mean")
first(Pairs.prom.regulatory.DNAse) = annotate.3D.Features(GRanges.INS, first(Pairs.prom.regulatory.DNAse), kind = "INS",aggregateFunction ="mean")

second(Pairs.prom.regulatory.DNAse) = annotate.3D.Features(GRanges.DI, second(Pairs.prom.regulatory.DNAse), kind = "DI",aggregateFunction ="mean")
second(Pairs.prom.regulatory.DNAse) = annotate.3D.Features(GRanges.FIREs, second(Pairs.prom.regulatory.DNAse), kind = "FIRE",aggregateFunction ="mean")
second(Pairs.prom.regulatory.DNAse) = annotate.3D.Features(GRanges.INS, second(Pairs.prom.regulatory.DNAse), kind = "INS",aggregateFunction ="mean")

first(Pairs.prom.regulatory.DNAse) = annotate.Activity(CTCFActivity.NEU, first(Pairs.prom.regulatory.DNAse), kind = "CTCF")
first(Pairs.prom.regulatory.DNAse) = annotate.Activity(H3K27acActivity.NEU, first(Pairs.prom.regulatory.DNAse), kind = "H3K27ac")
first(Pairs.prom.regulatory.DNAse) = annotate.Activity(DNAseActivity.NEU, first(Pairs.prom.regulatory.DNAse), kind = "DNAse")
first(Pairs.prom.regulatory.DNAse) = annotate.Activity(H3K4me1Activity.NEU, first(Pairs.prom.regulatory.DNAse), kind = "H3K4me1")
first(Pairs.prom.regulatory.DNAse) = annotate.Activity(H3K4me3Activity.NEU, first(Pairs.prom.regulatory.DNAse), kind = "H3K4me3")
first(Pairs.prom.regulatory.DNAse) = annotate.Activity(p300Activity.NEU, first(Pairs.prom.regulatory.DNAse), kind = "p300")
first(Pairs.prom.regulatory.DNAse) = annotate.Activity(H3K27me3Activity.NEU, first(Pairs.prom.regulatory.DNAse), kind = "H3K27me3")

second(Pairs.prom.regulatory.DNAse) = annotate.Activity(CTCFActivity.NEU, second(Pairs.prom.regulatory.DNAse), kind = "CTCF")
second(Pairs.prom.regulatory.DNAse) = annotate.Activity(H3K27acActivity.NEU, second(Pairs.prom.regulatory.DNAse), kind = "H3K27ac")
second(Pairs.prom.regulatory.DNAse) = annotate.Activity(DNAseActivity.NEU, second(Pairs.prom.regulatory.DNAse), kind = "DNAse")
second(Pairs.prom.regulatory.DNAse) = annotate.Activity(H3K4me1Activity.NEU, second(Pairs.prom.regulatory.DNAse), kind = "H3K4me1")
second(Pairs.prom.regulatory.DNAse) = annotate.Activity(H3K4me3Activity.NEU, second(Pairs.prom.regulatory.DNAse), kind = "H3K4me3")
second(Pairs.prom.regulatory.DNAse) = annotate.Activity(p300Activity.NEU, second(Pairs.prom.regulatory.DNAse), kind = "p300")
second(Pairs.prom.regulatory.DNAse) = annotate.Activity(H3K27me3Activity.NEU, second(Pairs.prom.regulatory.DNAse), kind = "H3K27me3")

forward.DNAse = unique(first(Pairs.prom.regulatory.DNAse))
reverse.DNAse = unique(second(Pairs.prom.regulatory.DNAse))


regul.DNAse = data.frame(mcols(forward.DNAse))
prom.DNAse = data.frame(mcols(reverse.DNAse))

gl.tmp = GRangesList(forward.DNAse, reverse.DNAse)
loci.DNAse = unique(do.call("c", as(gl.tmp, "GRangesList")))

FeaturesPlusActivity.DNAse = data.frame(mcols(loci.DNAse))

################################################################################################################
dist.prom.regul.DNAse = start(first(Pairs.prom.peaks.DNAse)) - start(second(Pairs.prom.peaks.DNAse))
summary(dist.prom.regul.DNAse)
summary(abs(dist.prom.regul.DNAse))

n.aval = sum(ifelse( start(first(Pairs.prom.peaks.DNAse)) < start(second(Pairs.prom.peaks.DNAse)),1,0))
n.amont= sum(ifelse(start(first(Pairs.prom.peaks.DNAse)) > start(second(Pairs.prom.peaks.DNAse)),1,0))

n.aval/(n.amont+n.aval)
#Correspondance between RRCs based on Rao definition and DNAse and TADs
#Creation of Pairs
start.Pair.Rao = sapply(Pairs.prom.regulatory.Rao, function(x) min(start(first(x)), start(second(x))))
start.Pair.DNAse = sapply(Pairs.prom.regulatory.DNAse, function(x) min(start(first(x)), start(second(x))))

end.Pair.Rao =  sapply(Pairs.prom.regulatory.Rao, function(x) max(end(first(x)), end(second(x))))
end.Pair.DNAse =  sapply(Pairs.prom.regulatory.DNAse, function(x) max(end(first(x)), end(second(x))))

GRanges.Pair.Rao = GRanges(seqnames = seqnames(first(Pairs.prom.regulatory.Rao)), ranges=IRanges(start=start.Pair.Rao,end=end.Pair.Rao, names=paste0("Pair", 1:length(Pairs.prom.regulatory.Rao))))
GRanges.Pair.DNAse = GRanges(seqnames = seqnames(first(Pairs.prom.regulatory.DNAse)), ranges=IRanges(start=start.Pair.DNAse,end=end.Pair.DNAse, names=paste0("Pair",1:length(Pairs.prom.regulatory.DNAse))))

#CRNs Building for Rao and DNAse definition
prom.Rao = first(Pairs.prom.regulatory.Rao)
regul.Rao = second(Pairs.prom.regulatory.Rao)
df.Rao = cbind(data.frame(prom.Rao), data.frame(regul.Rao))

#colnames(df.Rao) = c("chr.x", "startProm", "endProm", "width.x", "strand.x", "gene_id.x", "geneSymbol.x", "chr.y", "start", "end", "width.y", "strand.y", "gene_id.y", "geneSymbol.y")
colnames(df.Rao) = c("chr.x", "start.x", "end.x", "width", "strand", "gene_id.x","geneSymbol.x","mean_DI.x","mean_ScoreFire.x","mean_INS.x","mean_CTCF.x",
                     "mean_H3K27ac.x", "mean_DNAse.x","mean_H3K4me1.x", "mean_H3K4me3.x","mean_p300.x","chr.y", "start.y", "end.y", "width", "strand", "gene_id.y", "geneSymbol.y",
                     "mean_DI.y","mean_ScoreFire.y","mean_INS.y","mean_CTCF.y","mean_H3K27ac.y", "mean_DNAse.y","mean_H3K4me1.y", "mean_H3K4me3.y","mean_p300.y")

df.Rao$name = gsub(" ", "",apply(df.Rao[,c("chr.y","start.y", "end.y")], 1, paste0, collapse=":"))

df.Rao.RRCs = df.Rao[,c("chr.x", "start.x", "end.x","geneSymbol.x", "chr.y", "start.y", "end.y","name")]


peaks.DNase1 = first(Pairs.prom.regulatory.DNAse)
peaks.DNAse2 = second(Pairs.prom.regulatory.DNAse)
df.DNAse = cbind(data.frame(peaks.DNase1), data.frame(peaks.DNAse2))

#colnames(df.DNAse) = c("chr.x", "startPeak", "endPeak", "width.x", "strand.x", "gene_id.x", "geneSymbol.x", "chr.y", "startProm", "endProm", "width.y", "strand.y", "gene_id.y", "geneSymbol.y")
colnames(df.DNAse) = c("chr.x", "startPeak", "endPeak", "width", "strand", "gene_id.x","geneSymbol.x" ,"mean_DI.x","mean_ScoreFire.x","mean_INS.x","mean_CTCF.x",
                       "mean_H3K27ac.x","mean_DNAse.x", "mean_H3K4me1.x", "mean_H3K4me3.x","mean_p300.x", "chr.y", "start.y", "end.y", "width", "strand","gene_id.y", "geneSymbol.y",
                       "mean_DI.y", "mean_ScoreFire.y", "mean_INS.y", "mean_CTCF.y", "mean_H3K27ac.y" ,"mean_DNAse.y","mean_H3K4me1.y", "mean_H3K4me3.y","mean_p300.y" )

df.DNAse$name = gsub(" ", "",apply(df.DNAse[,c("chr.x","startPeak", "endPeak")], 1, paste0, collapse=":"))

nodes.Rao = df.Rao.RRCs[,c("geneSymbol.x", "name")]
nodes.Rao = na.omit(nodes.Rao)
vertices.Rao = unique(rbind(matrix(nodes.Rao$geneSymbol.x, ncol=1), matrix(nodes.Rao$name, ncol=1)))

length(unique(nodes.Rao$geneSymbol.x))
length(unique(nodes.Rao$name))

summary(as.numeric(table(nodes.Rao$geneSymbol.x)))
summary(as.numeric(table(nodes.Rao$name)))

length(table(compo.Rao$csize))
t.test(summary(as.numeric(table(nodes.Rao$geneSymbol.x))), summary(as.numeric(table(nodes.Rao$name))))

nodes.DNAse = df.DNAse[,c("name", "geneSymbol.y")]
nodes.DNAse = na.omit(nodes.DNAse)
vertices.DNAse = unique(rbind(matrix(nodes.DNAse$name, ncol=1), matrix(nodes.DNAse$geneSymbol.y, ncol=1)))

graph_from_Rao = graph_from_data_frame(nodes.Rao, directed = F, vertices=vertices.Rao)
compo.Rao = components(graph_from_Rao)

summary(as.numeric(compo.Rao$csize))

graph_from_DNAse = graph_from_data_frame(nodes.DNAse, directed = F, vertices=vertices.DNAse)
compo.DNAse = components(graph_from_DNAse)

summary(as.numeric(compo.DNAse$csize))

summary(as.numeric(table(nodes.DNAse$name)))
summary(as.numeric(table(nodes.DNAse$geneSymbol.y)))

t.test(summary(as.numeric(table(nodes.DNAse$name))),summary(as.numeric(table(nodes.DNAse$geneSymbol.y))))

#CRNs for DNAse and Rao Methodology
df.membership.DNAse = data.frame(compo.DNAse$membership)
df.membership.DNAse$name = rownames(df.membership.DNAse)
df.membership.DNAse = df.membership.DNAse[df.membership.DNAse$name%in%nodes.DNAse$name,]

df.DNAse = merge(df.DNAse, df.membership.DNAse, by="name")

df.DNAse$minStart = apply(df.DNAse[,c("startPeak", "startProm")], 1, min)
df.DNAse$maxEnd = apply(df.DNAse[,c("endPeak", "endProm")], 1, max)

start.cluster.DNAse = aggregate(minStart~compo.DNAse.membership, df.DNAse, min)
end.cluster.DNAse = aggregate(maxEnd~compo.DNAse.membership, df.DNAse, max)

df.start.end.DNAse = unique(merge(merge(start.cluster.DNAse, end.cluster.DNAse, by="compo.DNAse.membership"), df.DNAse[,c("chr.x", "compo.DNAse.membership")], by="compo.DNAse.membership"))
GRanges.cluster.DNAse = GRanges(seqnames = df.start.end.DNAse$chr, ranges = IRanges(start=df.start.end.DNAse$minStart, end=df.start.end.DNAse$maxEnd, names=df.start.end.DNAse$compo.DNAse.membership))

df.membership.Rao = data.frame(compo.Rao$membership)
df.membership.Rao$geneSymbol.x = rownames(df.membership.Rao)
df.membership.Rao = df.membership.Rao[df.membership.Rao$geneSymbol.x%in%nodes.Rao$geneSymbol.x,]

df.Rao = merge(df.Rao, df.membership.Rao, by="geneSymbol.x")

df.Rao$minStart = apply(df.Rao[,c("startProm", "start")], 1, min)
df.Rao$maxEnd = apply(df.Rao[,c("endProm", "end")], 1, max)

start.cluster.Rao = aggregate(minStart~compo.Rao.membership, df.Rao, min)
end.cluster.Rao = aggregate(maxEnd~compo.Rao.membership, df.Rao, max)

df.start.end.Rao = unique(merge(merge(start.cluster.Rao, end.cluster.Rao, by="compo.Rao.membership"), df.Rao[,c("chr.x", "compo.Rao.membership")], by="compo.Rao.membership"))
GRanges.cluster.Rao = GRanges(seqnames = df.start.end.Rao$chr, ranges = IRanges(start=df.start.end.Rao$minStart, end=df.start.end.Rao$maxEnd, names=df.start.end.Rao$compo.Rao.membership))


#######################################################################################################
summary(width(GRanges.Pair.Rao))
summary(width(GRanges.cluster.Rao))

t.test(width(GRanges.cluster.Rao), width(GRanges.Pair.Rao))
wilcox.test(width(GRanges.cluster.Rao),  width(GRanges.Pair.Rao))


summary(width(GRanges.Pair.DNAse))
summary(width(GRanges.cluster.DNAse))

t.test(width(GRanges.cluster.DNAse), width(GRanges.Pair.DNAse))
wilcox.test(width(GRanges.cluster.DNAse),  width(GRanges.Pair.DNAse))
#######################################################################################################
#3D features and CRNs
#Rao
#INS
meanRRCsX.INS = aggregate(mean_INS.x~compo.Rao.membership, df.Rao, mean, na.rm=T)
meanRRCsY.INS = aggregate(mean_INS.y~compo.Rao.membership, df.Rao, mean, na.rm=T)

meanRRCs.Rao.INS = merge(meanRRCsX.INS,meanRRCsY.INS)
cor(meanRRCs.Rao.INS$mean_INS.x, meanRRCs.Rao.INS$mean_INS.y)
#[1] 0.4120274

meanRRCs.Rao.INS$meanINS = apply(meanRRCs.Rao.INS[,c("mean_INS.x","mean_INS.y")], 1, mean, na.rm=T)

#DI
meanRRCsX.DI = aggregate(mean_DI.x~compo.Rao.membership, df.Rao, mean, na.rm=T)
meanRRCsY.DI = aggregate(mean_DI.y~compo.Rao.membership, df.Rao, mean, na.rm=T)

meanRRCs.Rao.DI = merge(meanRRCsX.DI, meanRRCsY.DI)

cor(meanRRCs.Rao.DI$mean_DI.x, meanRRCs.Rao.DI$mean_DI.y)
#[1] 0.07281198

meanRRCs.Rao.DI$meanDI = apply(meanRRCs.Rao.DI[,c("mean_DI.x","mean_DI.y")], 1, mean, na.rm=T)

#Score-FIRE
meanRRCsX.FIRE = aggregate(mean_ScoreFire.x~compo.Rao.membership, df.Rao, mean, na.rm=T)
meanRRCsY.FIRE = aggregate(mean_ScoreFire.y~compo.Rao.membership, df.Rao, mean, na.rm=T)

meanRRCs.Rao.FIRE = merge(meanRRCsX.FIRE, meanRRCsY.FIRE)

cor(meanRRCs.Rao.FIRE$mean_ScoreFire.x, meanRRCs.Rao.FIRE$mean_ScoreFire.y)
#-0.02386936

meanRRCs.Rao.FIRE$meanFIRE = apply(meanRRCs.Rao.FIRE[,c("mean_ScoreFire.x","mean_ScoreFire.y")], 1, mean, na.rm=T)

#DNAse
meanRRCsX.DNAse = aggregate(mean_DNAse.x~compo.Rao.membership, df.Rao, mean, na.rm=T)
meanRRCsY.DNAse = aggregate(mean_DNAse.y~compo.Rao.membership, df.Rao, mean, na.rm=T)

meanRRCs.Rao.DNAse = merge(meanRRCsX.DNAse, meanRRCsY.DNAse)

cor(meanRRCs.Rao.DNAse$mean_DNAse.x, meanRRCs.Rao.DNAse$mean_DNAse.y)
#[1] 0.05788059

meanRRCs.Rao.DNAse$meanDNAse = apply(meanRRCs.Rao.DNAse[,c("mean_DNAse.x","mean_DNAse.y")], 1, mean, na.rm=T)

#CTCF
meanRRCsX.CTCF = aggregate(mean_CTCF.x~compo.Rao.membership, df.Rao, mean, na.rm=T)
meanRRCsY.CTCF = aggregate(mean_CTCF.y~compo.Rao.membership, df.Rao, mean, na.rm=T)

meanRRCs.Rao.CTCF = merge(meanRRCsX.CTCF, meanRRCsY.CTCF)

cor(meanRRCs.Rao.CTCF$mean_CTCF.x,meanRRCs.Rao.CTCF$mean_CTCF.y)
#[1] 0.05333757

meanRRCs.Rao.CTCF$meanCTCF = apply(meanRRCs.Rao.CTCF[,c("mean_CTCF.x","mean_CTCF.y")], 1, mean, na.rm=T)

#H3K27ac
meanRRCsX.H3K27ac = aggregate(mean_H3K27ac.x~compo.Rao.membership, df.Rao, mean, na.rm=T)
meanRRCsY.H3K27ac = aggregate(mean_H3K27ac.y~compo.Rao.membership, df.Rao, mean, na.rm=T)

meanRRCs.Rao.H3K27ac = merge(meanRRCsX.H3K27ac, meanRRCsY.H3K27ac)

cor(meanRRCs.Rao.H3K27ac$mean_H3K27ac.x,meanRRCs.Rao.H3K27ac$mean_H3K27ac.y)
#[1] 0.1494804

meanRRCs.Rao.H3K27ac$meanH3K27ac= apply(meanRRCs.Rao.H3K27ac[,c("mean_H3K27ac.x","mean_H3K27ac.y")], 1, mean, na.rm=T)

#H3K4me1
meanRRCsX.H3K4me1 = aggregate(mean_H3K4me1.x~compo.Rao.membership, df.Rao, mean, na.rm=T)
meanRRCsY.H3K4me1 = aggregate(mean_H3K4me1.y~compo.Rao.membership, df.Rao, mean, na.rm=T)

meanRRCs.Rao.H3K4me1 = merge(meanRRCsX.H3K4me1, meanRRCsY.H3K4me1)

cor(meanRRCs.Rao.H3K4me1$mean_H3K4me1.x,meanRRCs.Rao.H3K4me1$mean_H3K4me1.y)
#[1] 0.1688802

meanRRCs.Rao.H3K4me1$meanH3K4me1= apply(meanRRCs.Rao.H3K4me1[,c("mean_H3K4me1.x","mean_H3K4me1.y")], 1, mean, na.rm=T)

#H3K4me3
meanRRCsX.H3K4me3 = aggregate(mean_H3K4me3.x~compo.Rao.membership, df.Rao, mean, na.rm=T)
meanRRCsY.H3K4me3 = aggregate(mean_H3K4me3.y~compo.Rao.membership, df.Rao, mean, na.rm=T)

meanRRCs.Rao.H3K4me3 = merge(meanRRCsX.H3K4me3, meanRRCsY.H3K4me3)

cor(meanRRCs.Rao.H3K4me3$mean_H3K4me3.x,meanRRCs.Rao.H3K4me3$mean_H3K4me3.y)
#[1] 0.09681651

meanRRCs.Rao.H3K4me3$meanH3K4me3= apply(meanRRCs.Rao.H3K4me3[,c("mean_H3K4me3.x","mean_H3K4me3.y")], 1, mean, na.rm=T)

#p300
meanRRCsX.p300 = aggregate(mean_p300.x~compo.Rao.membership, df.Rao, mean, na.rm=T)
meanRRCsY.p300 = aggregate(mean_p300.y~compo.Rao.membership, df.Rao, mean, na.rm=T)

meanRRCs.Rao.p300 = merge(meanRRCsX.p300, meanRRCsY.p300)

cor(meanRRCs.Rao.p300$mean_p300.x,meanRRCs.Rao.p300$mean_p300.y)
#[1] 0.0540233

meanRRCs.Rao.p300$meanp300= apply(meanRRCs.Rao.p300[,c("mean_p300.x","mean_p300.y")], 1, mean, na.rm=T)


#Expression 
colnames(expression.genes.RRCs) = c("geneSymbol.x", "Expression")

df.Rao = merge(df.Rao, expression.genes.RRCs, by="geneSymbol.x")
meanRRCs.Expre.Rao = aggregate(Expression~compo.Rao.membership,df.Rao,mean,na.rm=T)

Full3D.Rao = merge(meanRRCs.Rao.FIRE[,c("compo.Rao.membership","meanFIRE")],merge(meanRRCs.Rao.DI[,c("compo.Rao.membership","meanDI")], meanRRCs.Rao.INS[,c("compo.Rao.membership","meanINS")], by="compo.Rao.membership"), by="compo.Rao.membership")
FullAct.Rao = merge(meanRRCs.Rao.DNAse[,c("compo.Rao.membership","meanDNAse")],merge(meanRRCs.Rao.H3K27ac[,c("compo.Rao.membership","meanH3K27ac")], meanRRCs.Rao.CTCF[,c("compo.Rao.membership","meanCTCF")], by="compo.Rao.membership"), by="compo.Rao.membership")
FullAct.Rao = merge(FullAct.Rao, merge(meanRRCs.Expre.Rao,merge(meanRRCs.Rao.H3K4me3[,c("compo.Rao.membership","meanH3K4me3")],merge (meanRRCs.Rao.H3K4me1[,c("compo.Rao.membership","meanH3K4me1")],meanRRCs.Rao.p300[,c("compo.Rao.membership", "meanp300")], by="compo.Rao.membership"), by = "compo.Rao.membership"), by="compo.Rao.membership"), by="compo.Rao.membership")

Full3DAct.Rao = merge(Full3D.Rao, FullAct.Rao, by="compo.Rao.membership")

rownames(Full3DAct.Rao) = Full3DAct.Rao$compo.Rao.membership
Full3DAct.Rao$compo.Rao.membership = NULL


#DNAse
#INS
meanRRCspeak1.INS = aggregate(mean_INS.x~compo.DNAse.membership, df.DNAse, mean, na.rm=T)
meanRRCspeak2.INS = aggregate(mean_INS.y~compo.DNAse.membership, df.DNAse, mean, na.rm=T)

meanRRCs.DNAse.INS = merge(meanRRCspeak1.INS,meanRRCspeak2.INS)
cor(meanRRCs.DNAse.INS$mean_INS.x, meanRRCs.DNAse.INS$mean_INS.y)
#[1] 0.4232474

meanRRCs.DNAse.INS$meanINS = apply(meanRRCs.DNAse.INS[,c("mean_INS.x","mean_INS.y")], 1, mean, na.rm=T)

#DI
meanRRCspeak1.DI = aggregate(mean_DI.x~compo.DNAse.membership, df.DNAse, mean, na.rm=T)
meanRRCspeak2.DI = aggregate(mean_DI.y~compo.DNAse.membership, df.DNAse, mean, na.rm=T)

meanRRCs.DNAse.DI = merge(meanRRCspeak1.DI, meanRRCspeak2.DI)

cor(meanRRCs.DNAse.DI$mean_DI.x, meanRRCs.DNAse.DI$mean_DI.y)
#[1] 0.04278193

meanRRCs.DNAse.DI$meanDI = apply(meanRRCs.DNAse.DI[,c("mean_DI.x","mean_DI.y")], 1, mean, na.rm=T)

#Score-FIRE
meanRRCspeak1.FIRE = aggregate(mean_ScoreFire.x~compo.DNAse.membership, df.DNAse, mean, na.rm=T)
meanRRCspeak2.FIRE = aggregate(mean_ScoreFire.y~compo.DNAse.membership, df.DNAse, mean, na.rm=T)

meanRRCs.DNAse.FIRE = merge(meanRRCspeak1.FIRE, meanRRCspeak2.FIRE)

cor(meanRRCs.DNAse.FIRE$mean_ScoreFire.x, meanRRCs.DNAse.FIRE$mean_ScoreFire.y)
#[1] 0.03335056

meanRRCs.DNAse.FIRE$meanFIRE = apply(meanRRCs.DNAse.FIRE[,c("mean_ScoreFire.x","mean_ScoreFire.y")], 1, mean, na.rm=T)

#DNAse
meanRRCspeak1.DNAse = aggregate(mean_DNAse.x~compo.DNAse.membership, df.DNAse, mean, na.rm=T)
meanRRCspeak2.DNAse = aggregate(mean_DNAse.y~compo.DNAse.membership, df.DNAse, mean, na.rm=T)

meanRRCs.DNAse.DNAse = merge(meanRRCspeak1.DNAse, meanRRCspeak2.DNAse)

cor(meanRRCs.DNAse.DNAse$mean_DNAse.x, meanRRCs.DNAse.DNAse$mean_DNAse.y)
#[1] 0.07382295

meanRRCs.DNAse.DNAse$meanDNAse = apply(meanRRCs.DNAse.DNAse[,c("mean_DNAse.x","mean_DNAse.y")], 1, mean, na.rm=T)

#CTCF

meanRRCspeak1.CTCF = aggregate(mean_CTCF.x~compo.DNAse.membership, df.DNAse, mean, na.rm=T)
meanRRCspeak2.CTCF = aggregate(mean_CTCF.y~compo.DNAse.membership, df.DNAse, mean, na.rm=T)

meanRRCs.DNAse.CTCF = merge(meanRRCspeak1.CTCF, meanRRCspeak2.CTCF)

cor(meanRRCs.DNAse.CTCF$mean_CTCF.x,meanRRCs.DNAse.CTCF$mean_CTCF.y)
#[1] 0.1142575

meanRRCs.DNAse.CTCF$meanCTCF = apply(meanRRCs.DNAse.CTCF[,c("mean_CTCF.x","mean_CTCF.y")], 1, mean, na.rm=T)

#H3K27ac
meanRRCspeak1.H3K27ac = aggregate(mean_H3K27ac.x~compo.DNAse.membership, df.DNAse, mean, na.rm=T)
meanRRCspeak2.H3K27ac = aggregate(mean_H3K27ac.y~compo.DNAse.membership, df.DNAse, mean, na.rm=T)

meanRRCs.DNAse.H3K27ac = merge(meanRRCspeak1.H3K27ac, meanRRCspeak2.H3K27ac)

cor(meanRRCs.DNAse.H3K27ac$mean_H3K27ac.x,meanRRCs.DNAse.H3K27ac$mean_H3K27ac.y)
#[1] 0.1230133

meanRRCs.DNAse.H3K27ac$meanH3K27ac= apply(meanRRCs.DNAse.H3K27ac[,c("mean_H3K27ac.x","mean_H3K27ac.y")], 1, mean, na.rm=T)

#H3K4me1
meanRRCspeak1.H3K4me1 = aggregate(mean_H3K4me1.x~compo.DNAse.membership, df.DNAse, mean, na.rm=T)
meanRRCspeak2.H3K4me1 = aggregate(mean_H3K4me1.y~compo.DNAse.membership, df.DNAse, mean, na.rm=T)

meanRRCs.DNAse.H3K4me1 = merge(meanRRCspeak1.H3K4me1, meanRRCspeak2.H3K4me1)

cor(meanRRCs.DNAse.H3K4me1$mean_H3K4me1.x,meanRRCs.DNAse.H3K4me1$mean_H3K4me1.y)
#[1] 0.07653141

meanRRCs.DNAse.H3K4me1$meanH3K4me1= apply(meanRRCs.DNAse.H3K4me1[,c("mean_H3K4me1.x","mean_H3K4me1.y")], 1, mean, na.rm=T)

#H3K4me3
meanRRCspeak1.H3K4me3 = aggregate(mean_H3K4me3.x~compo.DNAse.membership, df.DNAse, mean, na.rm=T)
meanRRCspeak2.H3K4me3 = aggregate(mean_H3K4me3.y~compo.DNAse.membership, df.DNAse, mean, na.rm=T)

meanRRCs.DNAse.H3K4me3 = merge(meanRRCspeak1.H3K4me3, meanRRCspeak2.H3K4me3)

cor(meanRRCs.DNAse.H3K4me3$mean_H3K4me3.x,meanRRCs.DNAse.H3K4me3$mean_H3K4me3.y)
#[1] 0.08338789

meanRRCs.DNAse.H3K4me3$meanH3K4me3= apply(meanRRCs.DNAse.H3K4me3[,c("mean_H3K4me3.x","mean_H3K4me3.y")], 1, mean, na.rm=T)


#p300
meanRRCspeak1.p300 = aggregate(mean_p300.x~compo.DNAse.membership, df.DNAse, mean, na.rm=T)
meanRRCspeak2.p300 = aggregate(mean_p300.y~compo.DNAse.membership, df.DNAse, mean, na.rm=T)

meanRRCs.DNAse.p300 = merge(meanRRCspeak1.p300, meanRRCspeak2.p300)

cor(meanRRCs.DNAse.p300$mean_p300.x,meanRRCs.DNAse.p300$mean_p300.y)
#[1] 0.0491476

meanRRCs.DNAse.p300$meanp300= apply(meanRRCs.DNAse.p300[,c("mean_p300.x","mean_p300.y")], 1, mean, na.rm=T)

#H3K27me3
meanRRCspeak1.H3K27me3 = aggregate(mean_H3K27me3.x~compo.DNAse.membership, df.DNAse, mean, na.rm=T)
meanRRCspeak2.H3K27me3 = aggregate(mean_H3K27me3.y~compo.DNAse.membership, df.DNAse, mean, na.rm=T)

meanRRCs.DNAse.H3K27me3 = merge(meanRRCspeak1.H3K27me3, meanRRCspeak2.H3K27me3)

cor(meanRRCs.DNAse.H3K27me3$mean_H3K27me3.x,meanRRCs.DNAse.H3K27me3$mean_H3K27me3.y)
#[1] [1] 0.05603322

meanRRCs.DNAse.H3K27me3$meanH3K27me3= apply(meanRRCs.DNAse.H3K27me3[,c("mean_H3K27me3.x","mean_H3K27me3.y")], 1, mean, na.rm=T)

colnames(expression.genes.RRCs) = c("geneSymbol.y", "Expression")

df.DNAse = merge(df.DNAse, expression.genes.RRCs, by="geneSymbol.y")
meanRRCs.Expre.DNAse = aggregate(Expression~compo.DNAse.membership,df.DNAse,mean,na.rm=T)


Full3D.DNAse = merge(meanRRCs.DNAse.FIRE[,c("compo.DNAse.membership","meanFIRE")],merge(meanRRCs.DNAse.DI[,c("compo.DNAse.membership","meanDI")], meanRRCs.DNAse.INS[,c("compo.DNAse.membership","meanINS")], by="compo.DNAse.membership"), by="compo.DNAse.membership")
FullAct.DNAse = merge(meanRRCs.DNAse.DNAse[,c("compo.DNAse.membership","meanDNAse")],merge(meanRRCs.DNAse.H3K27ac[,c("compo.DNAse.membership","meanH3K27ac")], meanRRCs.DNAse.CTCF[,c("compo.DNAse.membership","meanCTCF")], by="compo.DNAse.membership"), by="compo.DNAse.membership")
FullAct.DNAse = merge(FullAct.DNAse, merge(meanRRCs.Expre.DNAse, merge(meanRRCs.DNAse.H3K4me3[,c("compo.DNAse.membership","meanH3K4me3")],merge (meanRRCs.DNAse.H3K4me1[,c("compo.DNAse.membership","meanH3K4me1")],meanRRCs.DNAse.p300[,c("compo.DNAse.membership", "meanp300")], by="compo.DNAse.membership"), by = "compo.DNAse.membership"), by="compo.DNAse.membership"), by="compo.DNAse.membership")

Full3DAct.DNAse = merge(Full3D.DNAse, FullAct.DNAse, by="compo.DNAse.membership")

rownames(Full3DAct.DNAse) = Full3DAct.DNAse$compo.DNAse.membership
Full3DAct.DNAse$compo.DNAse.membership = NULL

###################################################################################################

df.Rao$mean.Pair.FIREs = apply(df.Rao[,c("mean_ScoreFire.x","mean_ScoreFire.y")], 1, mean,na.rm=T)
df.Rao$mean.Pair.INS = apply(df.Rao[,c("mean_INS.x","mean_INS.y")], 1, mean,na.rm=T)
df.Rao$mean.Pair.DI = apply(df.Rao[,c("mean_DI.y","mean_DI.y")], 1, mean, na.rm=T)
df.Rao$mean.Pair.H3K27ac = apply(df.Rao[,c("mean_H3K27ac.x","mean_H3K27ac.y")], 1, mean,na.rm=T)
df.Rao$mean.Pair.H3K27me3 = apply(df.Rao[,c("mean_H3K27me3.x","mean_H3K27me3.y")], 1, mean,na.rm=T)
df.Rao$mean.Pair.CTCF = apply(df.Rao[,c("mean_CTCF.x","mean_CTCF.y")], 1, mean, na.rm=T)
df.Rao$mean.Pair.DNAse = apply(df.Rao[,c("mean_DNAse.x","mean_DNAse.y")], 1, mean, na.rm=T)
df.Rao$mean.Pair.H3K4me1 = apply(df.Rao[,c("mean_H3K4me1.x","mean_H3K4me1.y")], 1, mean, na.rm=T)
df.Rao$mean.Pair.H3K4me3 = apply(df.Rao[,c("mean_H3K4me3.x","mean_H3K4me3.y")], 1, mean, na.rm=T)
df.Rao$mean.Pair.p300 = apply(df.Rao[,c("mean_p300.x","mean_p300.y")], 1, mean, na.rm=T)

df.Rao.subset  = df.Rao[ave(1:nrow(df.Rao), df.Rao$compo.Rao.membership, FUN = length)!=1,]

ICC.Rao = lapply(col2investiguate, 
                 function(x) ICC(random.model(x,"compo.Rao.membership",df.Rao.subset)))

Rao.ICC=data.frame(cbind(col2investiguate, do.call(rbind, ICC.Rao)))
colnames(Rao.ICC) = c("Metric", "ICC")
Rao.ICC$Method = "Rao"


df.DNAse$mean.Pair.FIREs = apply(df.DNAse[,c("mean_ScoreFire.x","mean_ScoreFire.y")], 1, mean,na.rm=T)
df.DNAse$mean.Pair.INS = apply(df.DNAse[,c("mean_INS.x","mean_INS.y")], 1, mean,na.rm=T)
df.DNAse$mean.Pair.DI = apply(df.DNAse[,c("mean_DI.y","mean_DI.y")], 1, mean,na.rm=T)
df.DNAse$mean.Pair.H3K27ac = apply(df.DNAse[,c("mean_H3K27ac.x","mean_H3K27ac.y")], 1, mean,na.rm=T)
df.DNAse$mean.Pair.H3K27me3 = apply(df.DNAse[,c("mean_H3K27me3.x","mean_H3K27me3.y")], 1, mean,na.rm=T)
df.DNAse$mean.Pair.CTCF = apply(df.DNAse[,c("mean_CTCF.x","mean_CTCF.y")], 1, mean,na.rm=T)
df.DNAse$mean.Pair.DNAse = apply(df.DNAse[,c("mean_DNAse.x","mean_DNAse.y")], 1, mean,na.rm=T)
df.DNAse$mean.Pair.H3K4me1 = apply(df.DNAse[,c("mean_H3K4me1.x","mean_H3K4me1.y")], 1, mean,na.rm=T)
df.DNAse$mean.Pair.H3K4me3 = apply(df.DNAse[,c("mean_H3K4me3.x","mean_H3K4me3.y")], 1, mean,na.rm=T)
df.DNAse$mean.Pair.p300 = apply(df.DNAse[,c("mean_p300.x","mean_p300.y")], 1, mean,na.rm=T)

df.DNAse.subset  = df.DNAse[ave(1:nrow(df.DNAse), df.DNAse$compo.DNAse.membership, FUN = length)!=1,]


ICC.DNAse = lapply(col2investiguate, 
                 function(x) ICC(random.model(x,"compo.DNAse.membership",df.DNAse.subset)))

DNAse.ICC = data.frame(cbind(col2investiguate, do.call(rbind, ICC.DNAse)))
colnames(DNAse.ICC) = c("Metric", "ICC")
DNAse.ICC$Method = "DNAse"

par(mfrow=c(3,1))

##################################################################################
colnames(expression.genes.RRCs) = c("geneSymbol.x", "Expression")
df.Rao = merge(df.Rao, expression.genes.RRCs, by="geneSymbol.x")

expr.Rao.model = lmer(Expression~(1|geneSymbol.x), df.Rao, REML=T)
summary(expr.Rao.model)
ICC(expr.Rao.model)

colnames(expression.genes.RRCs) = c("geneSymbol.y", "Expression")
df.DNAse = merge(df.DNAse, expression.genes.RRCs, by="geneSymbol.y")

expr.DNAse.model = lmer(Expression~(1|geneSymbol.y), df.DNAse, REML=T)
summary(expr.DNAse.model)

ICC(expr.DNAse.model, nested = T)

ggplot(unique(df.DNAse[df.DNAse$compo.DNAse.membership%in%1:100,c("Expression", "compo.DNAse.membership")]), aes(y=Expression, x=compo.DNAse.membership)) + geom_point(position = position_jitter(width=0.1, height=0.0))

fd = aggregate(geneSymbol.y~compo.DNAse.membership, df.DNAse, function(x) length(unique(x)))
fd.subset = df.DNAse[df.DNAse$compo.DNAse.membership%in% fd[fd$geneSymbol.y>1,"compo.DNAse.membership"],]

fd.1 = aggregate(geneSymbol.x~compo.Rao.membership, df.Rao, function(x) length(unique(x)))
fd.1.subset = df.Rao[df.Rao$compo.Rao.membership%in% fd.1[fd.1$geneSymbol.x>1,"compo.Rao.membership"],]

fd.2 = aggregate(TargetGene~membership, enhancers.promoters, function(x) length(unique(x)))
fd.2.subset = enhancers.promoters[enhancers.promoters$membership%in% fd.2[fd.2$TargetGene>1,"membership"],]

ICC(lmer(Expression~(1|membership),fd.2.subset, REML=T), nested=F)
ICC(lmer(Expression~(1|compo.DNAse.membership),fd.subset, REML=T), nested=F)
ICC(lmer(Expression~(1|compo.Rao.membership),fd.1.subset, REML=T), nested=F)


##########################################################################################
#Expression
par(mfrow=c(3,1))
plot(fd.2.subset[fd.2.subset$membership%in%1:100, "membership"], fd.2.subset[fd.2.subset$membership%in%1:100, "Expression"],
     ylab="Expression", xlab="Component", main="ABC")
lines(aggregate(Expression~membership,fd.2.subset[fd.2.subset$membership%in%1:100,], mean)$membership,aggregate(Expression~membership,fd.2.subset[fd.2.subset$membership%in%1:100,], mean)$Expression, col="red", lwd=2)

plot(fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100, "compo.Rao.membership"], fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100, "Expression"],
     ylab="Expression", xlab="Component", main="Rao")
lines(aggregate(Expression~compo.Rao.membership,fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100,], mean)$compo.Rao.membership,aggregate(Expression~compo.Rao.membership,fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100,], mean)$Expression, col="red")

plot(fd.subset[fd.subset$compo.DNAse.membership%in%1:100, "compo.DNAse.membership"], fd.subset[fd.subset$compo.DNAse.membership%in%1:100, "Expression"],
     ylab="Expression", xlab="Component", main="Rao+DNAse")
lines(aggregate(Expression~compo.DNAse.membership,fd.subset[fd.subset$compo.DNAse.membership%in%1:100,], mean)$compo.DNAse.membership,aggregate(Expression~compo.DNAse.membership,fd.subset[fd.subset$compo.DNAse.membership%in%1:100,], mean)$Expression, col="red")

###########################################################################################
#H3K27ac
plot(fd.2.subset[fd.2.subset$membership%in%1:100, "membership"], fd.2.subset[fd.2.subset$membership%in%1:100, "mean.Pair.H3K27ac"], main=paste("ICC:", round(as.numeric(as.character(ABC.ICC[ABC.ICC$Metric=="mean.Pair.H3K27ac", "ICC"])),2)),
     ylab="Mean H3K27ac By Pair", xlab="Component")
lines(aggregate(mean.Pair.H3K27ac~membership,fd.2.subset[fd.2.subset$membership%in%1:100,], mean)$membership,aggregate(mean.Pair.H3K27ac~membership,fd.2.subset[fd.2.subset$membership%in%1:100,], mean)$mean.Pair.H3K27ac, col="red", lwd=2)

plot(fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100, "compo.Rao.membership"], fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100, "mean.Pair.H3K27ac"],main=paste("ICC:", round(as.numeric(as.character(Rao.ICC[Rao.ICC$Metric=="mean.Pair.H3K27ac", "ICC"])),2)),
     ylab="Mean H3K27ac By Pair", xlab="Component")
lines(aggregate(mean.Pair.H3K27ac~compo.Rao.membership,fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100,], mean)$compo.Rao.membership,aggregate(mean.Pair.H3K27ac~compo.Rao.membership,fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100,], mean)$mean.Pair.H3K27ac, col="red")

plot(fd.subset[fd.subset$compo.DNAse.membership%in%1:100, "compo.DNAse.membership"], fd.subset[fd.subset$compo.DNAse.membership%in%1:100, "mean.Pair.H3K27ac"],main=paste("ICC:", round(as.numeric(as.character(DNAse.ICC[DNAse.ICC$Metric=="mean.Pair.H3K27ac", "ICC"])),2)),
     ylab="Mean H3K27ac By Pair", xlab="Component")
lines(aggregate(mean.Pair.H3K27ac~compo.DNAse.membership,fd.subset[fd.subset$compo.DNAse.membership%in%1:100,], mean)$compo.DNAse.membership,aggregate(mean.Pair.H3K27ac~compo.DNAse.membership,fd.subset[fd.subset$compo.DNAse.membership%in%1:100,], mean)$mean.Pair.H3K27ac, col="red")
########################################################################################
#DNAse
plot(fd.2.subset[fd.2.subset$membership%in%1:100, "membership"], fd.2.subset[fd.2.subset$membership%in%1:100, "mean.Pair.DNAse"], main=paste("ICC:", round(as.numeric(as.character(ABC.ICC[ABC.ICC$Metric=="mean.Pair.DNAse", "ICC"])),2)),
     ylab="Mean DNAse By Pair", xlab="Component")
lines(aggregate(mean.Pair.DNAse~membership,fd.2.subset[fd.2.subset$membership%in%1:100,], mean)$membership,aggregate(mean.Pair.DNAse~membership,fd.2.subset[fd.2.subset$membership%in%1:100,], mean)$mean.Pair.DNAse, col="red", lwd=2)

plot(fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100, "compo.Rao.membership"], fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100, "mean.Pair.DNAse"],main=paste("ICC:", round(as.numeric(as.character(Rao.ICC[Rao.ICC$Metric=="mean.Pair.H3K27ac", "ICC"])),2)),
     ylab="Mean DNAse By Pair", xlab="Component")
lines(aggregate(mean.Pair.DNAse~compo.Rao.membership,fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100,], mean)$compo.Rao.membership,aggregate(mean.Pair.DNAse~compo.Rao.membership,fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100,], mean)$mean.Pair.DNAse, col="red")

plot(fd.subset[fd.subset$compo.DNAse.membership%in%1:100, "compo.DNAse.membership"], fd.subset[fd.subset$compo.DNAse.membership%in%1:100, "mean.Pair.DNAse"],main=paste("ICC:", round(as.numeric(as.character(DNAse.ICC[DNAse.ICC$Metric=="mean.Pair.DNAse", "ICC"])),2)),
     ylab="Mean DNAse By Pair", xlab="Component")
lines(aggregate(mean.Pair.DNAse~compo.DNAse.membership,fd.subset[fd.subset$compo.DNAse.membership%in%1:100,], mean)$compo.DNAse.membership,aggregate(mean.Pair.DNAse~compo.DNAse.membership,fd.subset[fd.subset$compo.DNAse.membership%in%1:100,], mean)$mean.Pair.DNAse, col="red")
########################################################################################
#p300
plot(fd.2.subset[fd.2.subset$membership%in%1:100, "membership"], fd.2.subset[fd.2.subset$membership%in%1:100, "mean.Pair.p300"], main=paste("ICC:", round(as.numeric(as.character(ABC.ICC[ABC.ICC$Metric=="mean.Pair.p300", "ICC"])),2)),
     ylab="Mean p300 By Pair", xlab="Component")
lines(aggregate(mean.Pair.p300~membership,fd.2.subset[fd.2.subset$membership%in%1:100,], mean)$membership,aggregate(mean.Pair.p300~membership,fd.2.subset[fd.2.subset$membership%in%1:100,], mean)$mean.Pair.p300, col="red", lwd=2)

plot(fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100, "compo.Rao.membership"], fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100, "mean.Pair.p300"],main=paste("ICC:", round(as.numeric(as.character(Rao.ICC[Rao.ICC$Metric=="mean.Pair.p300", "ICC"])),2)),
     ylab="Mean p300 By Pair", xlab="Component")
lines(aggregate(mean.Pair.p300~compo.Rao.membership,fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100,], mean)$compo.Rao.membership,aggregate(mean.Pair.p300~compo.Rao.membership,fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100,], mean)$mean.Pair.p300, col="red")

plot(fd.subset[fd.subset$compo.DNAse.membership%in%1:100, "compo.DNAse.membership"], fd.subset[fd.subset$compo.DNAse.membership%in%1:100, "mean.Pair.p300"],main=paste("ICC:", round(as.numeric(as.character(DNAse.ICC[DNAse.ICC$Metric=="mean.Pair.p300", "ICC"])),2)),
     ylab="Mean p300 By Pair", xlab="Component")
lines(aggregate(mean.Pair.p300~compo.DNAse.membership,fd.subset[fd.subset$compo.DNAse.membership%in%1:100,], mean)$compo.DNAse.membership,aggregate(mean.Pair.p300~compo.DNAse.membership,fd.subset[fd.subset$compo.DNAse.membership%in%1:100,], mean)$mean.Pair.p300, col="red")
########################################################################################
#H3K4me1
plot(fd.2.subset[fd.2.subset$membership%in%1:100, "membership"], fd.2.subset[fd.2.subset$membership%in%1:100, "mean.Pair.H3K4me1"], main=paste("ICC:", round(as.numeric(as.character(ABC.ICC[ABC.ICC$Metric=="mean.Pair.H3K4me1", "ICC"])),2)),
     ylab="Mean H3K4me1 By Pair", xlab="Component")
lines(aggregate(mean.Pair.H3K4me1~membership,fd.2.subset[fd.2.subset$membership%in%1:100,], mean)$membership,aggregate(mean.Pair.H3K4me1~membership,fd.2.subset[fd.2.subset$membership%in%1:100,], mean)$mean.Pair.H3K4me1, col="red", lwd=2)

plot(fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100, "compo.Rao.membership"], fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100, "mean.Pair.H3K4me1"],main=paste("ICC:", round(as.numeric(as.character(Rao.ICC[Rao.ICC$Metric=="mean.Pair.H3K4me1", "ICC"])),2)),
     ylab="Mean H3K4me1 By Pair", xlab="Component")
lines(aggregate(mean.Pair.H3K4me1~compo.Rao.membership,fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100,], mean)$compo.Rao.membership,aggregate(mean.Pair.H3K4me1~compo.Rao.membership,fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100,], mean)$mean.Pair.H3K4me1, col="red")

plot(fd.subset[fd.subset$compo.DNAse.membership%in%1:100, "compo.DNAse.membership"], fd.subset[fd.subset$compo.DNAse.membership%in%1:100, "mean.Pair.H3K4me1"],main=paste("ICC:", round(as.numeric(as.character(DNAse.ICC[DNAse.ICC$Metric=="mean.Pair.H3K4me1", "ICC"])),2)),
     ylab="Mean H3K4me1 By Pair", xlab="Component")
lines(aggregate(mean.Pair.H3K4me1~compo.DNAse.membership,fd.subset[fd.subset$compo.DNAse.membership%in%1:100,], mean)$compo.DNAse.membership,aggregate(mean.Pair.H3K4me1~compo.DNAse.membership,fd.subset[fd.subset$compo.DNAse.membership%in%1:100,], mean)$mean.Pair.H3K4me1, col="red")

#########################################################################################
#H3K4me3
plot(fd.2.subset[fd.2.subset$membership%in%1:100, "membership"], fd.2.subset[fd.2.subset$membership%in%1:100, "mean.Pair.H3K4me3"], main=paste("ICC:", round(as.numeric(as.character(ABC.ICC[ABC.ICC$Metric=="mean.Pair.H3K4me3", "ICC"])),2)),
     ylab="Mean H3K4me3 By Pair", xlab="Component")
lines(aggregate(mean.Pair.H3K4me3~membership,fd.2.subset[fd.2.subset$membership%in%1:100,], mean)$membership,aggregate(mean.Pair.H3K4me3~membership,fd.2.subset[fd.2.subset$membership%in%1:100,], mean)$mean.Pair.H3K4me3, col="red", lwd=2)

plot(fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100, "compo.Rao.membership"], fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100, "mean.Pair.H3K4me3"],main=paste("ICC:", round(as.numeric(as.character(Rao.ICC[Rao.ICC$Metric=="mean.Pair.H3K4me3", "ICC"])),2)),
     ylab="Mean H3K4me3 By Pair", xlab="Component")
lines(aggregate(mean.Pair.H3K4me3~compo.Rao.membership,fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100,], mean)$compo.Rao.membership,aggregate(mean.Pair.H3K4me3~compo.Rao.membership,fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100,], mean)$mean.Pair.H3K4me3, col="red")

plot(fd.subset[fd.subset$compo.DNAse.membership%in%1:100, "compo.DNAse.membership"], fd.subset[fd.subset$compo.DNAse.membership%in%1:100, "mean.Pair.H3K4me3"],main=paste("ICC:", round(as.numeric(as.character(DNAse.ICC[DNAse.ICC$Metric=="mean.Pair.H3K4me3", "ICC"])),2)),
     ylab="Mean H3K4me3 By Pair", xlab="Component")
lines(aggregate(mean.Pair.H3K4me3~compo.DNAse.membership,fd.subset[fd.subset$compo.DNAse.membership%in%1:100,], mean)$compo.DNAse.membership,aggregate(mean.Pair.H3K4me3~compo.DNAse.membership,fd.subset[fd.subset$compo.DNAse.membership%in%1:100,], mean)$mean.Pair.H3K4me3, col="red")
###########################################################################################
#CTCF
plot(fd.2.subset[fd.2.subset$membership%in%1:100, "membership"], fd.2.subset[fd.2.subset$membership%in%1:100, "mean.Pair.CTCF"], main=paste("ICC:", round(as.numeric(as.character(ABC.ICC[ABC.ICC$Metric=="mean.Pair.CTCF", "ICC"])),2)),
     ylab="Mean CTCF By Pair", xlab="Component")
lines(aggregate(mean.Pair.CTCF~membership,fd.2.subset[fd.2.subset$membership%in%1:100,], mean)$membership,aggregate(mean.Pair.CTCF~membership,fd.2.subset[fd.2.subset$membership%in%1:100,], mean)$mean.Pair.CTCF, col="red", lwd=2)

plot(fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100, "compo.Rao.membership"], fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100, "mean.Pair.CTCF"],main=paste("ICC:", round(as.numeric(as.character(Rao.ICC[Rao.ICC$Metric=="mean.Pair.CTCF", "ICC"])),2)),
     ylab="Mean CTCF By Pair", xlab="Component")
lines(aggregate(mean.Pair.CTCF~compo.Rao.membership,fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100,], mean)$compo.Rao.membership,aggregate(mean.Pair.CTCF~compo.Rao.membership,fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100,], mean)$mean.Pair.CTCF, col="red")

plot(fd.subset[fd.subset$compo.DNAse.membership%in%1:100, "compo.DNAse.membership"], fd.subset[fd.subset$compo.DNAse.membership%in%1:100, "mean.Pair.CTCF"],main=paste("ICC:", round(as.numeric(as.character(DNAse.ICC[DNAse.ICC$Metric=="mean.Pair.CTCF", "ICC"])),2)),
     ylab="Mean CTCF By Pair", xlab="Component")
lines(aggregate(mean.Pair.CTCF~compo.DNAse.membership,fd.subset[fd.subset$compo.DNAse.membership%in%1:100,], mean)$compo.DNAse.membership,aggregate(mean.Pair.CTCF~compo.DNAse.membership,fd.subset[fd.subset$compo.DNAse.membership%in%1:100,], mean)$mean.Pair.CTCF, col="red")
############################################################################################
#FIREs
plot(fd.2.subset[fd.2.subset$membership%in%1:100, "membership"], fd.2.subset[fd.2.subset$membership%in%1:100, "mean.Pair.FIREs"], main=paste("ICC:", round(as.numeric(as.character(ABC.ICC[ABC.ICC$Metric=="mean.Pair.FIREs", "ICC"])),2)),
     ylab="Mean FIREs By Pair", xlab="Component")
lines(aggregate(mean.Pair.FIREs~membership,fd.2.subset[fd.2.subset$membership%in%1:100,], mean)$membership,aggregate(mean.Pair.FIREs~membership,fd.2.subset[fd.2.subset$membership%in%1:100,], mean)$mean.Pair.FIREs, col="red", lwd=2)

plot(fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100, "compo.Rao.membership"], fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100, "mean.Pair.FIREs"],main=paste("ICC:", round(as.numeric(as.character(Rao.ICC[Rao.ICC$Metric=="mean.Pair.FIREs", "ICC"])),2)),
     ylab="Mean FIREs By Pair", xlab="Component")
lines(aggregate(mean.Pair.FIREs~compo.Rao.membership,fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100,], mean)$compo.Rao.membership,aggregate(mean.Pair.FIREs~compo.Rao.membership,fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100,], mean)$mean.Pair.FIREs, col="red")

plot(fd.subset[fd.subset$compo.DNAse.membership%in%1:100, "compo.DNAse.membership"], fd.subset[fd.subset$compo.DNAse.membership%in%1:100, "mean.Pair.FIREs"],main=paste("ICC:", round(as.numeric(as.character(DNAse.ICC[DNAse.ICC$Metric=="mean.Pair.FIREs", "ICC"])),2)),
     ylab="Mean FIREs By Pair", xlab="Component")
lines(aggregate(mean.Pair.FIREs~compo.DNAse.membership,fd.subset[fd.subset$compo.DNAse.membership%in%1:100,], mean)$compo.DNAse.membership,aggregate(mean.Pair.FIREs~compo.DNAse.membership,fd.subset[fd.subset$compo.DNAse.membership%in%1:100,], mean)$mean.Pair.FIREs, col="red")
############################################################################################
#INS
plot(fd.2.subset[fd.2.subset$membership%in%1:100, "membership"], fd.2.subset[fd.2.subset$membership%in%1:100, "mean.Pair.INS"], main=paste("ICC:", round(as.numeric(as.character(ABC.ICC[ABC.ICC$Metric=="mean.Pair.INS", "ICC"])),2)),
     ylab="Mean INS By Pair", xlab="Component")
lines(aggregate(mean.Pair.INS~membership,fd.2.subset[fd.2.subset$membership%in%1:100,], mean)$membership,aggregate(mean.Pair.INS~membership,fd.2.subset[fd.2.subset$membership%in%1:100,], mean)$mean.Pair.INS, col="red", lwd=2)

plot(fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100, "compo.Rao.membership"], fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100, "mean.Pair.INS"],main=paste("ICC:", round(as.numeric(as.character(Rao.ICC[Rao.ICC$Metric=="mean.Pair.INS", "ICC"])),2)),
     ylab="Mean INS By Pair", xlab="Component")
lines(aggregate(mean.Pair.INS~compo.Rao.membership,fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100,], mean)$compo.Rao.membership,aggregate(mean.Pair.INS~compo.Rao.membership,fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100,], mean)$mean.Pair.INS, col="red")

plot(fd.subset[fd.subset$compo.DNAse.membership%in%1:100, "compo.DNAse.membership"], fd.subset[fd.subset$compo.DNAse.membership%in%1:100, "mean.Pair.INS"],main=paste("ICC:", round(as.numeric(as.character(DNAse.ICC[DNAse.ICC$Metric=="mean.Pair.INS", "ICC"])),2)),
     ylab="Mean INS By Pair", xlab="Component")
lines(aggregate(mean.Pair.INS~compo.DNAse.membership,fd.subset[fd.subset$compo.DNAse.membership%in%1:100,], mean)$compo.DNAse.membership,aggregate(mean.Pair.INS~compo.DNAse.membership,fd.subset[fd.subset$compo.DNAse.membership%in%1:100,], mean)$mean.Pair.INS, col="red")
############################################################################################
#DI
plot(fd.2.subset[fd.2.subset$membership%in%1:100, "membership"], fd.2.subset[fd.2.subset$membership%in%1:100, "mean.Pair.DI"], main=paste("ICC:", round(as.numeric(as.character(ABC.ICC[ABC.ICC$Metric=="mean.Pair.DI", "ICC"])),2)),
     ylab="Mean DI By Pair", xlab="Component")
lines(aggregate(mean.Pair.DI~membership,fd.2.subset[fd.2.subset$membership%in%1:100,], mean)$membership,aggregate(mean.Pair.DI~membership,fd.2.subset[fd.2.subset$membership%in%1:100,], mean)$mean.Pair.DI, col="red", lwd=2)

plot(fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100, "compo.Rao.membership"], fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100, "mean.Pair.DI"],main=paste("ICC:", round(as.numeric(as.character(Rao.ICC[Rao.ICC$Metric=="mean.Pair.DI", "ICC"])),2)),
     ylab="Mean DI By Pair", xlab="Component")
lines(aggregate(mean.Pair.DI~compo.Rao.membership,fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100,], mean)$compo.Rao.membership,aggregate(mean.Pair.DI~compo.Rao.membership,fd.1.subset[fd.1.subset$compo.Rao.membership%in%1:100,], mean)$mean.Pair.DI, col="red")

plot(fd.subset[fd.subset$compo.DNAse.membership%in%1:100, "compo.DNAse.membership"], fd.subset[fd.subset$compo.DNAse.membership%in%1:100, "mean.Pair.DI"],main=paste("ICC:", round(as.numeric(as.character(DNAse.ICC[DNAse.ICC$Metric=="mean.Pair.DI", "ICC"])),2)),
     ylab="Mean DI By Pair", xlab="Component")
lines(aggregate(mean.Pair.DI~compo.DNAse.membership,fd.subset[fd.subset$compo.DNAse.membership%in%1:100,], mean)$compo.DNAse.membership,aggregate(mean.Pair.DI~compo.DNAse.membership,fd.subset[fd.subset$compo.DNAse.membership%in%1:100,], mean)$mean.Pair.DI, col="red")
