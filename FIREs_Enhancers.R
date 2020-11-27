FIREs = read.table("/home/loic/Documents/HiC/code/FIREs/FIREs_NEU.txt", header=T)
superFIREs = read.table("/home/loic/Documents/HiC/code/FIREs/NEU_SuperFIREs.txt", header=T)

GRanges.FIREs = GRanges(seqnames=FIREs$chr, ranges=IRanges(start=FIREs$start, end=FIREs$end), ScoreFire = FIREs$NEU_neg_ln_pval)
GRanges.FIREs.signi = GRanges.FIREs[GRanges.FIREs$ScoreFire>=3]

GRanges.SuperFIREs = GRanges(seqnames=superFIREs$chr, ranges=IRanges(start=superFIREs$start, end=superFIREs$end), ScoreFire =superFIREs$cum_FIRE_score)

sample.chr15 = FIREs[FIREs$chr=="chr15",]
plot(sample.chr15[2000:4000, "NEU_neg_ln_pval"], type="l")
abline(h=3)

distance.EnhancersToFires = mcols(distanceToNearest(unique.Enhancers.RRCs, GRanges.FIREs.signi))$distance
distance.RegulToFires = mcols(distanceToNearest(first(Pairs.prom.regulatory.Rao), GRanges.FIREs.signi))$distance
distance.peakToFires = mcols(distanceToNearest(first(Pairs.prom.peaks.DNAse), GRanges.FIREs.signi))$distance

quantile(distance.EnhancersToFires, probs=seq(0,1,0.01))
quantile(distance.RegulToFires, probs=seq(0,1,0.01))
quantile(distance.peakToFires, probs=seq(0,1,0.01))

start(GRanges.FIREs.signi) = start(GRanges.FIREs.signi) - 200000
end(GRanges.FIREs.signi) = end(GRanges.FIREs.signi) + 200000

Overlaps.FIRE.Enhancers = findOverlaps(GRanges.FIREs, unique.Enhancers.RRCs)
subset.FIRE.Enhancers = GRanges.FIREs[subjectHits(Overlaps.FIRE.Enhancers)]

mcols(unique.Enhancers.RRCs)[unique(subjectHits(Overlaps.FIRE.Enhancers)), "FIRE"] = aggregate(subset.FIRE.Enhancers$ScoreFire, list(subjectHits(Overlaps.FIRE.Enhancers)), mean, na.rm=T)$x


Overlaps.INS.Enhancers = findOverlaps(GRanges.INS, unique.Enhancers.RRCs)
subset.INS.Enhancers = GRanges.INS[queryHits(Overlaps.INS.Enhancers)]

mcols(unique.Enhancers.RRCs)[unique(subjectHits(Overlaps.INS.Enhancers)), "INS"] = aggregate(subset.INS.Enhancers$normalizedScore, list(subjectHits(Overlaps.INS.Enhancers)), mean, na.rm=T)$x


v = data.frame(mcols(unique.Enhancers.RRCs))
v= na.omit(v)
cor(v, method="spearman")

Overlaps.INS.Enhancers = findOverlaps(GRanges.INS, unique.Enhancers.RRCs)
subset.INS.Enhancers = GRanges.INS[queryHits(Overlaps.INS.Enhancers)]

mcols(unique.Enhancers.RRCs)[unique(subjectHits(Overlaps.INS.Enhancers)), "INS"] = aggregate(subset.INS.Enhancers$normalizedScore, list(subjectHits(Overlaps.INS.Enhancers)), mean, na.rm=T)$x


start(GRanges.SuperFIREs) = start(GRanges.SuperFIREs) -200000
end(GRanges.SuperFIREs) = end(GRanges.SuperFIREs) +200000

table(countOverlaps(unique.Enhancers.RRCs, GRanges.FIREs.signi))/sum(table(countOverlaps(unique.Enhancers.RRCs, GRanges.FIREs.signi)))
table(countOverlaps(unique.Promoters.RRCs, GRanges.FIREs.signi))/sum(table(countOverlaps(unique.Promoters.RRCs, GRanges.FIREs.signi)))

table(countOverlaps(unique.Enhancers.RRCs, GRanges.SuperFIREs))/sum(table(countOverlaps(unique.Enhancers.RRCs, GRanges.SuperFIREs)))
table(countOverlaps(unique.Promoters.RRCs, GRanges.SuperFIREs))/sum(table(countOverlaps(unique.Promoters.RRCs, GRanges.SuperFIREs)))

table(countOverlaps(first(Pairs.prom.regulatory.Rao), GRanges.FIREs.signi))/sum(table(countOverlaps(first(Pairs.prom.regulatory.Rao), GRanges.FIREs.signi)))
table(countOverlaps(second(Pairs.prom.regulatory.Rao), GRanges.FIREs.signi))/sum(table(countOverlaps(second(Pairs.prom.regulatory.Rao), GRanges.FIREs.signi)))

table(countOverlaps(first(Pairs.prom.regulatory.Rao), GRanges.SuperFIREs))/sum(table(countOverlaps(first(Pairs.prom.regulatory.Rao), GRanges.SuperFIREs)))
table(countOverlaps(second(Pairs.prom.regulatory.Rao), GRanges.SuperFIREs))/sum(table(countOverlaps(second(Pairs.prom.regulatory.Rao), GRanges.SuperFIREs)))

table(countOverlaps(first(Pairs.prom.peaks.DNAse), GRanges.FIREs.signi))/sum(table(countOverlaps(first(Pairs.prom.peaks.DNAse), GRanges.FIREs.signi)))
table(countOverlaps(second(Pairs.prom.peaks.DNAse), GRanges.FIREs.signi))/sum(table(countOverlaps(second(Pairs.prom.peaks.DNAse), GRanges.FIREs.signi)))

table(countOverlaps(first(Pairs.prom.peaks.DNAse), GRanges.SuperFIREs))/sum(table(countOverlaps(first(Pairs.prom.peaks.DNAse), GRanges.SuperFIREs)))
table(countOverlaps(second(Pairs.prom.peaks.DNAse), GRanges.SuperFIREs))/sum(table(countOverlaps(second(Pairs.prom.peaks.DNAse), GRanges.SuperFIREs)))

