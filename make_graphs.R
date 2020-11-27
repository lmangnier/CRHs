library(ggplot2)
library(ggpubr)
detach(dplyr)

hist.ABC = hist(abs(dist.prom.enh)[abs(dist.prom.enh)<1000000], plot = F, breaks = 50)
hist.ABC$density = hist.ABC$counts/sum(hist.ABC$counts)

hist.DNAse = hist(abs(dist.prom.regul.DNAse)[abs(dist.prom.regul.DNAse)<=1000000], plot = F, breaks = 50)
hist.DNAse$density = hist.DNAse$counts/sum(hist.DNAse$counts)

hist.Rao = hist(abs(dist.prom.regul.Rao)[abs(dist.prom.regul.Rao)<1000000], plot=F, breaks=50)
hist.Rao$density = hist.Rao$counts/sum(hist.Rao$counts)

plot(hist.ABC, col=rgb(0,0,1,1/4),freq=F,xlab="Distance",ylab="%", main="Distance between pair of promoter and regulatory element")
plot(hist.Rao, col=rgb(1,0,0,1/4), add=T, freq=F) 
plot(hist.DNAse, col=rgb(0,1,0,1/4), add=T,freq=F)
legend("topright", c("ABC.Score", "Rao","DNAse"), col=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4),rgb(0,1,0,1/4)), lwd=10)


#Overlaping between Pairs  and TADs
#DI
par(mfrow=c(3,1))
countOverlaps.Pairs.RaoDI = as.numeric(countOverlaps(GRanges.Pair.Rao,GRanges.TADs.DI))
table.Overlaps.Pairs.Rao.DI = table(countOverlaps.Pairs.RaoDI)
aggregate.table.Overlaps.Pairs.Rao.DI = cbind(data.frame(matrix(table.Overlaps.Pairs.Rao.DI[0:5], ncol=5)), sum(table.Overlaps.Pairs.Rao.DI[6:length(table.Overlaps.Pairs.Rao.DI)]))
colnames(aggregate.table.Overlaps.Pairs.Rao.DI) = c("0", "1","2","3","4",">=5")
barplot(data.matrix(aggregate.table.Overlaps.Pairs.Rao.DI)/sum(data.matrix(aggregate.table.Overlaps.Pairs.Rao.DI)),col=rgb(1,0,0,1/4),main="Pairs Based on Rao", ylab="%", xlab="#overlapings")

countOverlaps.Pairs.DNAseDI = as.numeric(countOverlaps(GRanges.Pair.DNAse,GRanges.TADs.DI))
table.Overlaps.Pairs.DNAse.DI = table(countOverlaps.Pairs.DNAseDI)
aggregate.table.Overlaps.Pairs.DNAse.DI = cbind(data.frame(matrix(table.Overlaps.Pairs.DNAse.DI[0:5], ncol=5)), sum(table.Overlaps.Pairs.DNAse.DI[6:length(table.Overlaps.Pairs.DNAse.DI)]))
colnames(aggregate.table.Overlaps.Pairs.DNAse.DI) = c("0", "1","2","3","4",">=5")
barplot(data.matrix(aggregate.table.Overlaps.Pairs.DNAse.DI)/sum(data.matrix(aggregate.table.Overlaps.Pairs.DNAse.DI)),col = rgb(0,1,0,1/4),main="Pairs Based on DNAse", ylab="%", xlab="#overlapings")

countOverlaps.Pairs.ABCDI = as.numeric(countOverlaps(GRanges.Pair.ABC, GRanges.TADs.DI))
table.Overlaps.Pairs.ABC.DI = table(countOverlaps.Pairs.ABCDI)
aggregate.table.Overlaps.Pairs.ABC.DI = cbind(data.frame(matrix(table.Overlaps.Pairs.ABC.DI[0:5], ncol=5)), sum(table.Overlaps.Pairs.ABC.DI[6:length(table.Overlaps.Pairs.ABC.DI)]))
colnames(aggregate.table.Overlaps.Pairs.ABC.DI) = c("0", "1","2","3","4",">=5")
barplot(data.matrix(aggregate.table.Overlaps.Pairs.ABC.DI)/sum(data.matrix(aggregate.table.Overlaps.Pairs.ABC.DI)),col=rgb(0,0,1,1/4),main="Pairs Based on ABC-Score", ylab="%", xlab="#overlapings")

#INS
countOverlaps.Pairs.RaoINS = as.numeric(countOverlaps(GRanges.Pair.Rao,TADs.INS))
table.Overlaps.Pairs.Rao.INS = table(countOverlaps.Pairs.RaoINS)
aggregate.table.Overlaps.Pairs.Rao.INS = cbind(data.frame(matrix(table.Overlaps.Pairs.Rao.INS[0:5], ncol=5)), sum(table.Overlaps.Pairs.Rao.INS[6:length(table.Overlaps.Pairs.Rao.INS)]))
colnames(aggregate.table.Overlaps.Pairs.Rao.INS) = c("0", "1","2","3","4",">=5")
barplot(data.matrix(aggregate.table.Overlaps.Pairs.Rao.INS)/sum(data.matrix(aggregate.table.Overlaps.Pairs.Rao.INS)),col=rgb(1,0,0,1/4),main="Pairs Based on Rao", ylab="%", xlab="#overlapings")

countOverlaps.Pairs.DNAseINS = as.numeric(countOverlaps(GRanges.Pair.DNAse,TADs.INS))
table.Overlaps.Pairs.DNAse.INS = table(countOverlaps.Pairs.DNAseINS)
aggregate.table.Overlaps.Pairs.DNAse.INS = cbind(data.frame(matrix(table.Overlaps.Pairs.DNAse.INS[0:5], ncol=5)), sum(table.Overlaps.Pairs.DNAse.INS[6:length(table.Overlaps.Pairs.DNAse.INS)]))
colnames(aggregate.table.Overlaps.Pairs.DNAse.INS) = c("0", "1","2","3","4",">=5")
barplot(data.matrix(aggregate.table.Overlaps.Pairs.DNAse.INS)/sum(data.matrix(aggregate.table.Overlaps.Pairs.DNAse.INS)),col = rgb(0,1,0,1/4),main="Pairs Based on DNAse", ylab="%", xlab="#overlapings")


countOverlaps.Pairs.ABCINS = as.numeric(countOverlaps(GRanges.Pair.ABC, TADs.INS))
table.Overlaps.Pairs.ABC.INS = table(countOverlaps.Pairs.ABCINS)
aggregate.table.Overlaps.Pairs.ABC.INS = cbind(data.frame(matrix(table.Overlaps.Pairs.ABC.INS[0:5], ncol=5)), sum(table.Overlaps.Pairs.ABC.INS[6:length(table.Overlaps.Pairs.ABC.INS)]))
colnames(aggregate.table.Overlaps.Pairs.ABC.INS) = c("0", "1","2","3","4",">=5")
barplot(data.matrix(aggregate.table.Overlaps.Pairs.ABC.INS)/sum(data.matrix(aggregate.table.Overlaps.Pairs.ABC.INS)),col=rgb(0,0,1,1/4),main="Pairs Based on ABC-Score", ylab="%", xlab="#overlapings")

#Rao
#Pairs
countOverlaps.Pairs.RaoRao = as.numeric(countOverlaps(GRanges.Pair.Rao,TADs.Rao))
table.Overlaps.Pairs.Rao.Rao = table(countOverlaps.Pairs.RaoRao)
aggregate.table.Overlaps.Pairs.Rao.Rao = cbind(data.frame(matrix(table.Overlaps.Pairs.Rao.Rao[0:5], ncol=5)), sum(table.Overlaps.Pairs.Rao.Rao[6:length(table.Overlaps.Pairs.Rao.Rao)]))
colnames(aggregate.table.Overlaps.Pairs.Rao.Rao) = c("0", "1","2","3","4",">=5")
barplot(data.matrix(aggregate.table.Overlaps.Pairs.Rao.Rao)/sum(data.matrix(aggregate.table.Overlaps.Pairs.Rao.Rao)),col=rgb(1,0,0,1/4),main="Pairs Based on Rao", ylab="%", xlab="#overlapings")

countOverlaps.Pairs.DNAseRao = as.numeric(countOverlaps(GRanges.Pair.DNAse,TADs.Rao))
table.Overlaps.Pairs.DNAse.Rao = table(countOverlaps.Pairs.DNAseRao)
aggregate.table.Overlaps.Pairs.DNAse.Rao = cbind(data.frame(matrix(table.Overlaps.Pairs.DNAse.Rao[0:5], ncol=5)), sum(table.Overlaps.Pairs.DNAse.Rao[6:length(table.Overlaps.Pairs.DNAse.Rao)]))
colnames(aggregate.table.Overlaps.Pairs.DNAse.Rao) = c("0", "1","2","3","4",">=5")
barplot(data.matrix(aggregate.table.Overlaps.Pairs.DNAse.Rao)/sum(data.matrix(aggregate.table.Overlaps.Pairs.DNAse.Rao)),col = rgb(0,1,0,1/4),main="Pairs Based on DNAse", ylab="%", xlab="#overlapings")

countOverlaps.Pairs.ABCRao = as.numeric(countOverlaps(GRanges.Pair.ABC, TADs.Rao))
table.Overlaps.Pairs.ABC.Rao = table(countOverlaps.Pairs.ABCRao)
aggregate.table.Overlaps.Pairs.ABC.Rao = cbind(data.frame(matrix(table.Overlaps.Pairs.ABC.Rao[0:5], ncol=5)), sum(table.Overlaps.Pairs.ABC.Rao[6:length(table.Overlaps.Pairs.ABC.Rao)]))
colnames(aggregate.table.Overlaps.Pairs.ABC.Rao) = c("0", "1","2","3","4",">=5")
barplot(data.matrix(aggregate.table.Overlaps.Pairs.ABC.Rao)/sum(data.matrix(aggregate.table.Overlaps.Pairs.ABC.Rao)),col=rgb(0,0,1,1/4),main="Pairs Based on ABC-Score", ylab="%", xlab="#overlapings")

#CRNs and TADs
countOverlaps.CRNs.RaoINS = as.numeric(countOverlaps(GRanges.cluster.Rao,TADs.INS))
table.Overlaps.CRNs.RaoIns = table(countOverlaps.CRNs.RaoINS)
aggregate.table.Overlaps.CRNs.Rao.INS = cbind(data.frame(matrix(table.Overlaps.CRNs.RaoIns[0:5], ncol=5)), sum(table.Overlaps.CRNs.RaoIns[6:length(table.Overlaps.CRNs.RaoIns)]))
colnames(aggregate.table.Overlaps.CRNs.Rao.INS) = c("0", "1","2","3","4",">=5")


countOverlaps.CRNs.RaoRao  = as.numeric(countOverlaps(GRanges.cluster.Rao,TADs.Rao))
table.Overlaps.CRNs.RaoRao = table(countOverlaps.CRNs.RaoRao)
aggregate.table.Overlaps.CRNs.Rao.Rao = cbind(data.frame(matrix(table.Overlaps.CRNs.RaoRao[0:5], ncol=5)), sum(table.Overlaps.CRNs.RaoRao[6:length(table.Overlaps.CRNs.RaoRao)]))
colnames(aggregate.table.Overlaps.CRNs.Rao.Rao) = c("0", "1","2","3","4",">=5")

countOverlaps.CRNs.RaoDI = as.numeric(countOverlaps(GRanges.cluster.Rao,GRanges.TADs.DI))
table.Overlaps.CRNs.RaoDI = table(countOverlaps.CRNs.RaoDI)
aggregate.table.Overlaps.CRNs.Rao.DI = cbind(data.frame(matrix(table.Overlaps.CRNs.RaoDI[0:5], ncol=5)), sum(table.Overlaps.CRNs.RaoDI[6:length(table.Overlaps.CRNs.RaoDI)]))
colnames(aggregate.table.Overlaps.CRNs.Rao.DI) = c("0", "1","2","3","4",">=5")

countOverlaps.CRNs.DNAseINS = as.numeric(countOverlaps(GRanges.cluster.DNAse,TADs.INS))
table.Overlaps.CRNs.DNAseINS = table(countOverlaps.CRNs.DNAseINS)
aggregate.table.Overlaps.CRNs.DNAse.INS = cbind(data.frame(matrix(table.Overlaps.CRNs.DNAseINS[0:5], ncol=5)), sum(table.Overlaps.CRNs.DNAseINS[6:length(table.Overlaps.CRNs.DNAseINS)]))
colnames(aggregate.table.Overlaps.CRNs.DNAse.INS) = c("0", "1","2","3","4",">=5")


countOverlaps.CRNs.DNAseRao = as.numeric(countOverlaps(GRanges.cluster.DNAse,TADs.Rao))
table.Overlaps.CRNs.DNAseRao = table(countOverlaps.CRNs.DNAseRao)
aggregate.table.Overlaps.CRNs.DNAse.Rao = cbind(data.frame(matrix(table.Overlaps.CRNs.DNAseRao[0:5], ncol=5)), sum(table.Overlaps.CRNs.DNAseRao[6:length(table.Overlaps.CRNs.DNAseRao)]))
colnames(aggregate.table.Overlaps.CRNs.DNAse.Rao) = c("0", "1","2","3","4",">=5")

countOverlaps.CRNs.DNAseDI = as.numeric(countOverlaps(GRanges.cluster.DNAse,GRanges.TADs.DI))
table.Overlaps.CRNs.DNAseDI = table(countOverlaps.CRNs.DNAseDI)
aggregate.table.Overlaps.CRNs.DNAse.DI = cbind(data.frame(matrix(table.Overlaps.CRNs.DNAseDI[0:5], ncol=5)), sum(table.Overlaps.CRNs.DNAseDI[6:length(table.Overlaps.CRNs.DNAseDI)]))
colnames(aggregate.table.Overlaps.CRNs.DNAse.DI) = c("0", "1","2","3","4",">=5")



countOverlaps.CRNs.ABCINS = as.numeric(countOverlaps(GRanges.cluster.ABC,TADs.INS))
table.Overlaps.CRNs.ABCINS = table(countOverlaps.CRNs.ABCINS)
aggregate.table.Overlaps.CRNs.ABC.INS = cbind(data.frame(matrix(table.Overlaps.CRNs.ABCINS[0:5], ncol=5)), sum(table.Overlaps.CRNs.ABCINS[6:length(table.Overlaps.CRNs.ABCINS)]))
colnames(aggregate.table.Overlaps.CRNs.ABC.INS) = c("0", "1","2","3","4",">=5")

countOverlaps.CRNs.ABCRao = as.numeric(countOverlaps(GRanges.cluster.ABC,TADs.Rao))
table.Overlaps.CRNs.ABCRao = table(countOverlaps.CRNs.ABCRao)
aggregate.table.Overlaps.CRNs.ABC.Rao = cbind(data.frame(matrix(table.Overlaps.CRNs.ABCRao[0:5], ncol=5)), sum(table.Overlaps.CRNs.ABCRao[6:length(table.Overlaps.CRNs.ABCRao)]))
colnames(aggregate.table.Overlaps.CRNs.ABC.Rao) = c("0", "1","2","3","4",">=5")

countOverlaps.CRNs.ABCDI = as.numeric(countOverlaps(GRanges.cluster.ABC,GRanges.TADs.DI))
table.Overlaps.CRNs.ABCDI = table(countOverlaps.CRNs.ABCDI)
aggregate.table.Overlaps.CRNs.ABC.DI = cbind(data.frame(matrix(table.Overlaps.CRNs.ABCDI[0:5], ncol=5)), sum(table.Overlaps.CRNs.ABCDI[6:length(table.Overlaps.CRNs.ABCDI)]))
colnames(aggregate.table.Overlaps.CRNs.ABC.DI) = c("0", "1","2","3","4",">=5")



#Vizualization through piechart
par(mfrow=c(3,3))
paste0(round((data.matrix(aggregate.table.Overlaps.Pairs.Rao.Rao)/sum(data.matrix(aggregate.table.Overlaps.Pairs.Rao.Rao)))*100),"%")

pie(data.matrix(aggregate.table.Overlaps.Pairs.Rao.Rao), labels = names(aggregate.table.Overlaps.Pairs.Rao.Rao), main="Rao")
legend(x=-2.75,y=0.5, "AH")
pie(data.matrix(aggregate.table.Overlaps.Pairs.DNAse.Rao), labels = names(aggregate.table.Overlaps.Pairs.Rao.Rao), main="DNAse")
pie(data.matrix(aggregate.table.Overlaps.Pairs.ABC.Rao), labels = names(aggregate.table.Overlaps.Pairs.Rao.Rao), main="ABC-Score")


pie(data.matrix(aggregate.table.Overlaps.Pairs.Rao.DI), labels = names(aggregate.table.Overlaps.Pairs.Rao.Rao), main="Rao")
legend(x=-2.75,y=0.5, "DI")
pie(data.matrix(aggregate.table.Overlaps.Pairs.DNAse.DI), labels = names(aggregate.table.Overlaps.Pairs.Rao.Rao), main="DNAse")
pie(data.matrix(aggregate.table.Overlaps.Pairs.ABC.DI), labels = names(aggregate.table.Overlaps.Pairs.Rao.Rao), main="ABC-Score")

pie(data.matrix(aggregate.table.Overlaps.Pairs.Rao.INS), labels = names(aggregate.table.Overlaps.Pairs.Rao.Rao), main="Rao")
legend(x=-2.75,y=0.5, "INS")
pie(data.matrix(aggregate.table.Overlaps.Pairs.DNAse.INS), labels = names(aggregate.table.Overlaps.Pairs.Rao.Rao), main="DNAse")
pie(data.matrix(aggregate.table.Overlaps.Pairs.ABC.INS), labels = names(aggregate.table.Overlaps.Pairs.Rao.Rao), main="ABC-Score")

pie(data.matrix(aggregate.table.Overlaps.CRNs.Rao.Rao), labels = names(aggregate.table.Overlaps.CRNs.Rao.Rao), main="Rao")
legend(x=-2.75,y=0.5, "AH")
pie(data.matrix(aggregate.table.Overlaps.CRNs.DNAse.Rao), labels = names(aggregate.table.Overlaps.CRNs.DNAse.Rao), main = "DNAse")
pie(data.matrix(aggregate.table.Overlaps.CRNs.ABC.Rao), labels = names(aggregate.table.Overlaps.CRNs.ABC.Rao), main="ABC-Score")

pie(data.matrix(aggregate.table.Overlaps.CRNs.Rao.DI), labels = names(aggregate.table.Overlaps.CRNs.Rao.DI), main="Rao")
legend(x=-2.75,y=0.5, "DI")
pie(data.matrix(aggregate.table.Overlaps.CRNs.DNAse.DI), labels = names(aggregate.table.Overlaps.CRNs.DNAse.DI), main = "DNAse")
pie(data.matrix(aggregate.table.Overlaps.CRNs.ABC.DI), labels = names(aggregate.table.Overlaps.CRNs.ABC.DI),main="ABC-Score")

pie(data.matrix(aggregate.table.Overlaps.CRNs.Rao.INS), labels = names(aggregate.table.Overlaps.CRNs.Rao.INS), main="Rao")
legend(x=-2.75,y=0.5, "INS")
pie(data.matrix(aggregate.table.Overlaps.CRNs.DNAse.INS), labels = names(aggregate.table.Overlaps.CRNs.DNAse.INS), main = "DNAse")
pie(data.matrix(aggregate.table.Overlaps.CRNs.ABC.INS), labels = names(aggregate.table.Overlaps.CRNs.ABC.INS), main="ABC-Score")


#Overlaping between Pairs and genes
boxplot(as.numeric(countOverlaps(GRanges.Pair.Rao, genes.hg19)))
barplot(table(as.numeric(countOverlaps(GRanges.Pair.DNAse, genes.hg19))))
barplot(table(as.numeric(countOverlaps(GRanges.Pair.ABC, genes.hg19))))

#Numbers of connecting elements per method
enhancers.promoters$TargetGene = droplevels(enhancers.promoters$TargetGene)

barplot(table(as.numeric(table(first(Pairs.prom.regulatory.DNAse)))), col = rgb(0,1,0,1/4), xlab="#promoter/regulatory-elements" ,ylab="#elements", main="Number of connected elements by promoter/regulatory-element")
lines(x=1:length(table(as.numeric(table(first(Pairs.prom.regulatory.DNAse))))),y=table(as.numeric(table(first(Pairs.prom.regulatory.DNAse)))),col="green", lwd=3)
points(x=1:length(table(as.numeric(table(first(Pairs.prom.regulatory.DNAse))))),y=table(as.numeric(table(first(Pairs.prom.regulatory.DNAse)))), col="green", lwd=3)

barplot(table(as.numeric(table(enhancers.promoters$TargetGene))),col=rgb(0,0,1,1/4),add=T, xaxt="n")
lines(x=1:length(table(as.numeric(table(enhancers.promoters$TargetGene)))),y=table(as.numeric(table(enhancers.promoters$TargetGene))),col="blue", lwd=3)
points(x=1:length(table(as.numeric(table(enhancers.promoters$TargetGene)))),y=table(as.numeric(table(enhancers.promoters$TargetGene))), col="blue", lwd=3)


barplot(table(as.numeric(table(first(Pairs.prom.regulatory.Rao)))), col=rgb(1,0,0,1/4), add = T, xaxt="n")
lines(x=1:length(table(as.numeric(table(first(Pairs.prom.regulatory.Rao))))),y=table(as.numeric(table(first(Pairs.prom.regulatory.Rao)))),col="red", lwd=3)
points(x=1:length(table(as.numeric(table(first(Pairs.prom.regulatory.Rao))))),y=table(as.numeric(table(first(Pairs.prom.regulatory.Rao)))), col="red", lwd=3)
legend("topright", c("ABC.Score", "Rao","DNAse"), col=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4),rgb(0,1,0,1/4)), lwd=10)

#Number of CRN elements 
par(mfrow=c(1,1))
csizeABC = data.frame(matrix(table(compo.graph.ABC$csize), ncol=length(table(compo.graph.ABC$csize))))
colnames(csizeABC) = names(table(compo.graph.ABC$csize))
csizeDNAse = data.frame(matrix(table(compo.DNAse$csize), ncol=length(table(compo.DNAse$csize))))
colnames(csizeDNAse) = names(table(compo.DNAse$csize))
csizeRao = data.frame(matrix(table(compo.Rao$csize), ncol=length(table(compo.Rao$csize))))
colnames(csizeRao) = names(table(compo.Rao$csize))

library(dplyr)
full.elements.RRC = data.matrix(bind_rows(csizeABC,csizeDNAse,csizeRao))

barplot(full.elements.RRC[,1:19],col = c(rgb(0,0,1,1/4),rgb(0,1,0,1/4),rgb(1,0,0,1/4)), beside = T,xlab = "#Elements per CRN",ylab = "Count" ,main="Number of CRN elements by annotation types")
legend("topright", c("ABC.Score", "Rao","DNAse"), col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4),rgb(0,1,0,1/4)), lwd=10)

wilcox.test(compo.graph.ABC$csize, compo.Rao$csize)
wilcox.test(compo.graph.ABC$csize, compo.DNAse$csize)
wilcox.test(compo.Rao$csize, compo.DNAse$csize)

boxplot(compo.graph.ABC$csize[compo.graph.ABC$csize<=80], compo.Rao$csize[compo.Rao$csize<=80], compo.DNAse$csize[compo.DNAse$csize<=80],
        names=c("ABC.Score", "Rao", "DNAse"), col=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4),rgb(0,1,0,1/4)), ylab="#elements", main="Distribution of CRN size by annotation types")
points(c(26.31,3.232,6.036), bg=24,pch=24,col="red")

#Complexity
Complexity.ABC = Complexity.Crn(graph.ABC, enhancers.promoters$name, enhancers.promoters$TargetGene, extract.central.genes = T)
Complexity.Rao = Complexity.Crn(graph_from_Rao, df.Rao$name, df.Rao$geneSymbol.x, extract.central.genes = T)
Complexity.DNAse = Complexity.Crn(graph_from_DNAse, df.DNAse$name, df.DNAse$geneSymbol.y, extract.central.genes = T)

Complexity.ABC$monogamous;Complexity.Rao$monogamous;Complexity.DNAse$monogamous
Complexity.ABC$`11N`;Complexity.Rao$`11N`;Complexity.DNAse$`11N`

Complexity.ABC$Centrality$n.genes/length(compo.graph.ABC$csize)
Complexity.Rao$Centrality$n.genes/length(compo.Rao$csize)
Complexity.DNAse$Centrality$n.genes/length(compo.DNAse$csize)

Complexity.ABC$monogamous;Complexity.Rao$monogamous;Complexity.DNAse$monogamous
matrix.complexity = matrix(c(0.018,0.18,0.32, 0.028,0.28,0.50, 0.008,0.09,0.24,0.067,0.62,0.64)*100, ncol=4)
colnames(matrix.complexity) = c("Promoters 1-1", "Promoters 1-1-N", "Regulatory 1-1", "Regulatory 1-1-N")
barplot(matrix.complexity, col = c(rgb(0,0,1,1/4),rgb(0,1,0,1/4),rgb(1,0,0,1/4)), beside = T, ylab = "%", main = "Promoter/Regulatory-Element Complexity Relationships for AL1C Definition")
legend("topleft", c("ABC.Score", "Rao","DNAse"), col=c(rgb(0,0,1,1/4), rgb(1,0,0,1/4), rgb(0,1,0,1/4)), lwd=10)

par(mfrow=c(1,1))
#Central Genes
barplot(c(0.2823,0.2847,0.1603)*100,col = c(rgb(0,0,1,1/4), rgb(0,1,0,1/4), rgb(1,0,0,1/4)), ylab="%", names=c("ABC.Score", "DNAse", "Rao"),
        main="% of CRNs with central genes")


#Justification of using CRNs in relation to naive-based pair approach
density.Pair.Cluster.Rao = data.frame(matrix(c(width(GRanges.cluster.Rao)[width(GRanges.cluster.Rao)<=3000000], width(GRanges.Pair.Rao)[width(GRanges.Pair.Rao)<=3000000]), ncol=2))
colnames(density.Pair.Cluster.Rao) = c("CRN", "Pair")

dens.Rao <- apply(density.Pair.Cluster.Rao, 2, density, na.rm=T)

plot(NA,ylab="Density" ,xlab="Distance",xlim=range(sapply(dens.Rao, "[", "x")), ylim=range(sapply(dens.Rao, "[", "y")), main="Density of distance CRNs versus Pairs")
mapply(lines, dens.Rao, col=1:length(dens.Rao))

legend("topright", legend=names(dens.Rao), fill=1:length(dens.Rao))
ks.test(width(GRanges.Pair.Rao), width(GRanges.cluster.Rao))


density.Pair.Cluster.DNAse = data.frame(matrix(c(width(GRanges.cluster.DNAse)[width(GRanges.cluster.DNAse)<=3000000], width(GRanges.Pair.DNAse)[width(GRanges.Pair.DNAse)<=3000000]), ncol=2))
colnames(density.Pair.Cluster.DNAse) = c("CRN", "Pair")

dens.DNAse <- apply(density.Pair.Cluster.DNAse, 2, density, na.rm=T)

plot(NA,ylab="Density" ,xlab="Distance",xlim=range(sapply(dens.DNAse, "[", "x")), ylim=range(sapply(dens.DNAse, "[", "y")),main="Density of distance CRNs versus Pairs")
mapply(lines, dens.DNAse, col=1:length(dens.DNAse))

legend("topright", legend=names(dens.DNAse), fill=1:length(dens.DNAse))
ks.test(width(GRanges.Pair.DNAse), width(GRanges.cluster.DNAse))

density.Pair.Cluster.ABC = data.frame(matrix(c(width(GRanges.cluster.ABC)[width(GRanges.cluster.ABC)<=3000000], width(GRanges.Pair.ABC)[width(GRanges.Pair.ABC)<=3000000]), ncol=2))
colnames(density.Pair.Cluster.ABC) = c("CRN", "Pair")

dens.ABC <- apply(density.Pair.Cluster.ABC, 2, density, na.rm=T)

plot(NA,ylab="Density" ,xlab="Distance",xlim=range(sapply(dens.ABC, "[", "x")), ylim=range(sapply(dens.ABC, "[", "y")),main="Density of distance CRNs versus Pairs")
mapply(lines, dens.ABC, col=1:length(dens.ABC))
legend("topright", legend=names(dens.ABC), fill=1:length(dens.ABC))

ks.test(width(GRanges.Pair.ABC), width(GRanges.cluster))

################################################################################
library(karyoploteR)
library(regioneR)

mcols(GRanges.cluster)[names(GRanges.cluster)%in%Full.3D.RRCs.mean$compo.graph.membership,c("meanFIRE", "meanDI", "meanINS")] = Full.3D.RRCs.mean[,c("meanFIRE", "meanDI", "meanINS")]

GRanges.cluster.INS = GRanges.cluster[,"meanINS"]
GRanges.cluster.INS$y = GRanges.cluster.INS$meanINS
GRanges.cluster.INS$meanINS = NULL

GRanges.cluster.DI = GRanges.cluster[,"meanDI"]
GRanges.cluster.DI$y = GRanges.cluster.DI$meanDI
GRanges.cluster.DI$meanDI = NULL


GRanges.cluster.FIRE = GRanges.cluster[,"meanFIRE"]
GRanges.cluster.FIRE$y = GRanges.cluster.FIRE$meanFIRE
GRanges.cluster.FIRE$meanFIRE = NULL

kp <- plotKaryotype(plot.type = 3, chromosomes = chroms[chroms!="chrX"])
kpDataBackground(kp, data.panel = 1)
#kpRect(kp, data = GRanges.cluster.INS[!is.na(mcols(GRanges.cluster.INS)$y)], y0=0, y1=1, col="#FFDDDD", border=NA, r0=0, r1=1)
kpAxis(kp, ymin = 0, ymax = 1, r0=0, r1=1, numticks = 5, col="#666666", cex=2)
kpPoints(kp, data=GRanges.cluster.INS[!is.na(mcols(GRanges.cluster.INS)$y)], pch=16, cex=2, col="orange", r0=0, r1=1, ymin = 0, ymax = 1)
kpAbline(kp, h=mean(GRanges.INS$normalizedScore), col="blue", ymin=0, ymax=1, r0=0, r1=1, lwd=2)
kpRect(kp, chr="chr1", x0=100000000, x1=450000000, y0=0.25, y1=0.5, r0=0, r1=1, col="#EEEEEE", border="#666666")
kpText(kp, chr="chr1",x=300000000, y=0.375, col="red", r0=0, r1=1, labels="INS", cex=2)

kp <- plotKaryotype(plot.type = 3, chromosomes = chroms[chroms!="chrX"])
#kpRect(kp, data = GRanges.cluster.DI[!is.na(mcols(GRanges.cluster.DI)$y)], y0=0, y1=1, col="#FFDDDD", border=NA, r0=0, r1=1)
kpDataBackground(kp, data.panel = 1)
kpAxis(kp, ymin = 0, ymax = 600, r0=0, r1=1, numticks = 5, col="#666666", cex=2)
kpPoints(kp, data=GRanges.cluster.DI[!is.na(mcols(GRanges.cluster.DI)$y)], pch=16, cex=2, ymin = 0, ymax = 600,col="orange", r0=0, r1=1)
kpAbline(kp, h=mean(GRanges.DI$DI), col="blue", ymin=0, ymax=600, r0=0, r1=1, lwd=2)
kpRect(kp, chr="chr1", x0=100000000, x1=450000000, y0=0.5, y1=0.65, r0=0, r1=1, col="#EEEEEE", border="#666666")
kpText(kp, chr="chr1",x=300000000, y=0.57, col="red", r0=0, r1=1, labels="DI", cex=2)

kp <- plotKaryotype(plot.type = 3, chromosomes = chroms[chroms!="chrX"])
kpDataBackground(kp, data.panel = 1)
#kpRect(kp, data = GRanges.cluster.FIRE[!is.na(mcols(GRanges.cluster.FIRE)$y)], y0=0, y1=1, col="#FFDDDD", border=NA, r0=0, r1=1)
kpAxis(kp, ymin = 0, ymax = 4, r0=0, r1=1, numticks = 5, col="#666666", cex=2)
kpPoints(kp, data=GRanges.cluster.FIRE[!is.na(mcols(GRanges.cluster.FIRE)$y)], pch=16, cex=2, ymin = 0, ymax = 4,col="orange", r0=0, r1=1)
kpAbline(kp, h=mean(GRanges.FIREs$ScoreFire), col="blue", ymin=0, ymax=4, r0=0, r1=1, lwd=2)
kpRect(kp, chr="chr1", x0=100000000, x1=450000000, y0=0.75, y1=0.90, r0=0, r1=1, col="#EEEEEE", border="#666666")
kpText(kp, chr="chr1",x=275000000, y=0.82, col="red", r0=0, r1=1, labels="FIRE", cex=2)
