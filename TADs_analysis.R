#TADs COnstruction and Descriptive Analysis

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
#Because the output of software developed by Crane et al., 2015 gives boundary TADs, we have to process the file to have whole domains
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
TADs.Rao = import("AnnotationsForJuicer/TADs/TADs_10Kb_Rao.bedpe", format="bedpe")

length(first(TADs.Rao))
#4252

sum(width(first(TADs.Rao)))
#1 701 990 000

summary(width(first(TADs.Rao)))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#120000  200000  280000  400280  410000 5470000



#Descriptive analysis for comparing TADs 
sd(width(GRanges.TADs.DI));sd(width(TADs.INS));sd(width(first(TADs.Rao)))

hist(width(TADs.INS))
hist(width(GRanges.TADs.DI))
hist(width(first(TADs.Rao)))

#Do the TAD coverage distributions come from the same distribution?
ks.test(width(TADs.INS), width(GRanges.TADs.DI))
ks.test(width(TADs.INS), width(first(TADs.Rao)))
ks.test(width(first(TADs.Rao)), width(GRanges.TADs.DI))

OverlappingINSRao = subsetByOverlaps(TADs.INS,TADs.Rao)
OverlappingINSDI = subsetByOverlaps(TADs.INS,GRanges.TADs.DI)
OverlappingRaoDI = subsetByOverlaps(first(TADs.Rao),GRanges.TADs.DI)
#3866 overlaping TADs between Arrowhead and Insulation Score
#5109 overlaping TADs between DI and Insulation Score
#4235 overlaping TADs between DI and Arrowhead

OverlappingINSRaoDI = subsetByOverlaps(subsetByOverlaps(TADs.INS,TADs.Rao),GRanges.TADs.DI)
#3841 Consensus TADs: first attempt (have to be improved)

#Link between pairs and TADs

##############################################################################################
################################Pairs & TADs#################################################
##############################################################################################
#Overlapping between pairs of elements and TADs
#The analysis is done for each method of TAD construction

chisq.PairesTADS= matrix(0, ncol=2, nrow=3)

#DI
overlapsPairsABCDI = table(countOverlaps(GRanges.Pair.ABC , GRanges.TADs.DI))
#0     1     2     3     4     5     6     7     8     9    10    11    12   103   104 
#2684 50471  7048  1663   495   148    81    39    16     3     3     2     1     1     3 
overlapsPairsABCDI[-1]/sum(overlapsPairsABCDI[-1])

#Percentage of pairs overlaping at least one TAD
#           1            2            3            4            5            6            7            8            9           10 
#8.415480e-01 1.175176e-01 2.772868e-02 8.253577e-03 2.467736e-03 1.350585e-03 6.502818e-04 2.667823e-04 5.002168e-05 5.002168e-05 
#11           12          103          104 
#3.334778e-05 1.667389e-05 1.667389e-05 5.002168e-05 

overlapsPairsRaoDI = table(countOverlaps(GRanges.Pair.Rao , GRanges.TADs.DI))
overlapsPairsRaoDI[-1]/sum(overlapsPairsRaoDI[-1])

overlapsPairsDNAseDI = table(countOverlaps(GRanges.Pair.DNAse , GRanges.TADs.DI))
overlapsPairsDNAseDI[-1]/sum(overlapsPairsDNAseDI[-1])



DI.1=overlapsPairsDI[-1][1]
sup1DI = sum(overlapsPairsDI[-1][2:length(overlapsPairsDI[-1])])

chisq.PairesTADS[1,] = c(DI.1, sup1DI)

#INS
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

#Arrowhead
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

chisq.test(chisq.PairesTADS, correct = F)
#Dependance between the building TAD methodology and the number of overlaping pairs 
chisq.RRCsTADS= matrix(0, ncol=2, nrow=3)

#############################################################################################
#################################CRNs and TADs###############################################
#############################################################################################

#Same kind of analysis as previously with pairs
overlapsABCDI = table(countOverlaps(GRanges.cluster.ABC , GRanges.TADs.DI))
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   21   22   24   30   42  109 
#56 1044  167   70   45   45   37   30   31   22   20    9   13    9    8    6    3    4    2    4    1    3    1    1    1    1 

overlapsRaoDI = table(countOverlaps(GRanges.cluster.Rao , GRanges.TADs.DI))
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   16   24   26   27   29   31   33   44   51   52   60   67   71   96 
# 4 1777  326   44   17   11   17   12    2    2    3    1    1    6    1    1    2    1    1    1    3    1    1    1    2    3    1    1    2 
# 99  112  153  217  277  322 
# 1    1    1    1    1    1 

overlapsDNAseDI = table(countOverlaps(GRanges.cluster.DNAse , GRanges.TADs.DI))
# 0    1    2    3    4    5    6    7    8    9   10   11   13   14   16   17   24   26   27   28   30   31   45   49   52   55   60   67   71 
# 5 2328  356   53   15    2   10    6    1    5    1    2    1    1    1    1    1    1    1    1    1    2    1    5    2    1    1    1    1 
# 86   96   98  153  172  277 
# 1    1    1    1    1    1 


overlapsABCDI[-1]/sum(overlapsABCDI[-1])
#1            2            3            4            5            6            7            8            9           10 
#0.6620164870 0.1058972733 0.0443880786 0.0285351934 0.0285351934 0.0234622701 0.0190234623 0.0196575777 0.0139505390 0.0126823082 
#11           12           13           14           15           16           17           18           19           21 
#0.0057070387 0.0082435003 0.0057070387 0.0050729233 0.0038046925 0.0019023462 0.0025364616 0.0012682308 0.0025364616 0.0006341154 
#22           24           30           42          109 
#0.0019023462 0.0006341154 0.0006341154 0.0006341154 0.0006341154 


overlapsRaoDI[-1]/sum(overlapsRaoDI[-1])
# 1            2            3            4            5            6            7            8            9           10           11 
# 0.7908322207 0.1450823320 0.0195816644 0.0075656431 0.0048954161 0.0075656431 0.0053404539 0.0008900757 0.0008900757 0.0013351135 0.0004450378 
# 12           13           14           16           24           26           27           29           31           33           44 
# 0.0004450378 0.0026702270 0.0004450378 0.0004450378 0.0008900757 0.0004450378 0.0004450378 0.0004450378 0.0013351135 0.0004450378 0.0004450378 
# 51           52           60           67           71           96           99          112          153          217          277 
# 0.0004450378 0.0008900757 0.0013351135 0.0004450378 0.0004450378 0.0008900757 0.0004450378 0.0004450378 0.0004450378 0.0004450378 0.0004450378 
# 322 
# 0.0004450378 

overlapsDNAseDI[-1]/sum(overlapsDNAseDI[-1])
# 1            2            3            4            5            6            7            8            9           10           11 
# 0.8290598291 0.1267806268 0.0188746439 0.0053418803 0.0007122507 0.0035612536 0.0021367521 0.0003561254 0.0017806268 0.0003561254 0.0007122507 
# 13           14           16           17           24           26           27           28           30           31           45 
# 0.0003561254 0.0003561254 0.0003561254 0.0003561254 0.0003561254 0.0003561254 0.0003561254 0.0003561254 0.0003561254 0.0007122507 0.0003561254 
# 49           52           55           60           67           71           86           96           98          153          172 
# 0.0017806268 0.0007122507 0.0003561254 0.0003561254 0.0003561254 0.0003561254 0.0003561254 0.0003561254 0.0003561254 0.0003561254 0.0003561254 
# 277 
# 0.0003561254 

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

#Proportion of CRNs fully included in TADs (analysis done by method)
OverlapsRRCINS = findOverlaps(GRanges.cluster, TADs.INS)

#INS
uniqueOverlapsRRCINS = names(table(queryHits(OverlapsRRCINS))[table(queryHits(OverlapsRRCINS)) == 1])
subsetOverlapsRRCINS = OverlapsRRCINS[queryHits(OverlapsRRCINS)%in%uniqueOverlapsRRCINS]

RRCsPleinementInclusINS = sum(ifelse(start(GRanges.cluster[queryHits(subsetOverlapsRRCINS)]) >= start(TADs.INS[subjectHits(subsetOverlapsRRCINS)]) & end(GRanges.cluster[queryHits(subsetOverlapsRRCINS)])<=end(TADs.INS[subjectHits(subsetOverlapsRRCINS)]), 1,0))
#Proportion of one-overlaping CRNs which are fully included in TAD (denominator: TADs which overlap one TAD)
RRCsPleinementInclusINS/overlapsRRCsINS[2]
#Global proportion of one-overlaping CRNs which are fully included in TAD (denominator: all TADs)
RRCsPleinementInclusINS/sum(overlapsRRCsINS[-1])

#DI
OverlapsRRCDI = findOverlaps(GRanges.cluster, GRanges.TADs.DI)

uniqueOverlapsRRCDI = names(table(queryHits(OverlapsRRCDI))[table(queryHits(OverlapsRRCDI)) == 1])
subsetOverlapsRRCDI = OverlapsRRCDI[queryHits(OverlapsRRCDI)%in%uniqueOverlapsRRCDI]

RRCsPleinementInclusDI = sum(ifelse(start(GRanges.cluster[queryHits(subsetOverlapsRRCDI)]) >= start(GRanges.TADs.DI[subjectHits(subsetOverlapsRRCDI)]) & end(GRanges.cluster[queryHits(subsetOverlapsRRCDI)])<=end(GRanges.TADs.DI[subjectHits(subsetOverlapsRRCDI)]), 1,0))

RRCsPleinementInclusDI/overlapsRRCsDI[2]
RRCsPleinementInclusDI/sum(overlapsRRCsDI[-1])

#Rao
OverlapsRRCRao = findOverlaps(GRanges.cluster, TADs.Rao)

uniqueOverlapsRRCRao = names(table(queryHits(OverlapsRRCRao))[table(queryHits(OverlapsRRCRao)) == 1])
subsetOverlapsRRCRao = OverlapsRRCRao[queryHits(OverlapsRRCRao)%in%uniqueOverlapsRRCRao]

RRCsPleinementInclusRao = sum(ifelse(start(GRanges.cluster[queryHits(subsetOverlapsRRCRao)]) >= start(TADs.Rao[subjectHits(subsetOverlapsRRCRao)]) & end(GRanges.cluster[queryHits(subsetOverlapsRRCRao)])<=end(TADs.Rao[subjectHits(subsetOverlapsRRCRao)]), 1,0))

RRCsPleinementInclusRao/overlapsRRCsRao[2]
RRCsPleinementInclusRao/sum(overlapsRRCsRao[-1])


#Relationship between TAD coverage and number of overlaping CRNs
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



distance.EnhancersToBorders = mcols(distanceToNearest(unique.Enhancers.RRCs, GRanges.boundaries.noverlap))$distance
distance.RegulToBorders = mcols(distanceToNearest(first(Pairs.prom.regulatory.Rao), GRanges.boundaries.noverlap))$distance
distance.peakToBorders = mcols(distanceToNearest(first(Pairs.prom.peaks.DNAse), GRanges.boundaries.noverlap))$distance

summary(distance.EnhancersToBorders)
summary(distance.RegulToBorders)
summary(distance.peakToBorders)
