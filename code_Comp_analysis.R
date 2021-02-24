load("taloggcc.RData")


GRanges.Pair.ABC = coverage.By.Pair(ABC.Pairs)
GRanges.Pair.Rao = coverage.By.Pair(Rao.Pairs)
GRanges.Pair.DNAse = coverage.By.Pair(DNAse.Pairs)

#Merging consecutive GRanges

compartments <- unlist(reduce(split(taloggcc, ~Compartment)))
compartments$Compartment <- names(compartments)

names(compartments) = 1:length(compartments)

overlaps.ABC.Compartments = findOverlaps(GRanges.Pair.ABC, compartments)
overlaps.Rao.Compartments = findOverlaps(GRanges.Pair.Rao, compartments)
overlaps.DNAse.Compartments = findOverlaps(GRanges.Pair.DNAse, compartments)

for(q in unique(queryHits(overlaps.Rao.Compartments))){
  subj = subjectHits(overlaps.Rao.Compartments[queryHits(overlaps.Rao.Compartments) == q])
  mcols(GRanges.Pair.Rao)[q,"Compartment"] = paste(mcols(compartments)[subj,"Compartment"][c(1,length(mcols(compartments)[subj,"Compartment"]))], collapse="")
  
}

for(q in unique(queryHits(overlaps.ABC.Compartments))){
  subj = subjectHits(overlaps.ABC.Compartments[queryHits(overlaps.ABC.Compartments) == q])
  mcols(GRanges.Pair.ABC)[q,"Compartment"] = paste(mcols(compartments)[subj,"Compartment"][c(1,length(mcols(compartments)[subj,"Compartment"]))], collapse="")
  
}

for(q in unique(queryHits(overlaps.DNAse.Compartments))){
  subj = subjectHits(overlaps.DNAse.Compartments[queryHits(overlaps.DNAse.Compartments) == q])
  mcols(GRanges.Pair.DNAse)[q,"Compartment"] = paste(mcols(compartments)[subj,"Compartment"][c(1,length(mcols(compartments)[subj,"Compartment"]))], collapse="")
  
}

GRanges.Pair.ABC$Compartment = ifelse(GRanges.Pair.ABC$Compartment=="BA","AB", GRanges.Pair.ABC$Compartment)
GRanges.Pair.Rao$Compartment = ifelse(GRanges.Pair.Rao$Compartment=="BA","AB", GRanges.Pair.Rao$Compartment)
GRanges.Pair.DNAse$Compartment = ifelse(GRanges.Pair.DNAse$Compartment=="BA","AB", GRanges.Pair.DNAse$Compartment)

df.Rao.Comp = data.frame(cbind("Rao",matrix(table(GRanges.Pair.Rao$Compartment), ncol=1), c("AA","AB","BB")))
df.ABC.Comp = data.frame(cbind("ABC",matrix(table(GRanges.Pair.ABC$Compartment), ncol=1), c("AA","AB","BB")))
df.DNAse.Comp = data.frame(cbind("DNAse",matrix(table(GRanges.Pair.DNAse$Compartment), ncol=1), c("AA","AB","BB")))

df.all.Comp = rbind(df.Rao.Comp,df.ABC.Comp,df.DNAse.Comp)
df.all.Comp$X2 = as.numeric(as.character(df.all.Comp$X2))
ggplot(df.all.Comp, aes(x=X1,y=X2,fill=X3))+ geom_bar(position="fill", stat="identity")+scale_y_continuous(labels = scales::percent_format())+
  ylab("% Overlap") + xlab("Annotation Method") + labs(fill = "Compartment")

GRanges.cluster.ABC = coverage.By.Crn(graph_from_ABC, regul.promoters.ABC, method = "ABC")$clusters
GRanges.cluster.Rao = coverage.By.Crn(graph_from_Rao, df.Rao, method = "Rao")$clusters
GRanges.cluster.DNAse = coverage.By.Crn(graph_from_DNAse, df.DNase, method="DNAse")$clusters

boxplot(width(GRanges.cluster.ABC)[-width(GRanges.cluster.ABC))])
boxplot(width(GRanges.cluster.Rao)[width(GRanges.cluster.Rao)<=1.0e+08])
summary(width(GRanges.cluster.DNAse))

overlaps.ABC.CRNs.Compartments = findOverlaps(GRanges.cluster.ABC, compartments)
overlaps.Rao.CRNs.Compartments = findOverlaps(GRanges.cluster.Rao, compartments)
overlaps.DNAse.CRNs.Compartments = findOverlaps(GRanges.cluster.DNAse, compartments)

for(q in unique(queryHits(overlaps.ABC.CRNs.Compartments))){
  subj = subjectHits(overlaps.ABC.CRNs.Compartments[queryHits(overlaps.ABC.CRNs.Compartments) == q])
  mcols(GRanges.cluster.ABC)[q,"Compartment"] = paste(mcols(compartments)[subj,"Compartment"][c(1,length(mcols(compartments)[subj,"Compartment"]))], collapse="")
  
}

for(q in unique(queryHits(overlaps.Rao.CRNs.Compartments))){
  subj = subjectHits(overlaps.Rao.CRNs.Compartments[queryHits(overlaps.Rao.CRNs.Compartments) == q])
  mcols(GRanges.cluster.Rao)[q,"Compartment"] = paste(mcols(compartments)[subj,"Compartment"][c(1,length(mcols(compartments)[subj,"Compartment"]))], collapse="")
  
}

for(q in unique(queryHits(overlaps.DNAse.CRNs.Compartments))){
  subj = subjectHits(overlaps.DNAse.CRNs.Compartments[queryHits(overlaps.DNAse.CRNs.Compartments) == q])
  mcols(GRanges.cluster.DNAse)[q,"Compartment"] = paste(mcols(compartments)[subj,"Compartment"][c(1,length(mcols(compartments)[subj,"Compartment"]))], collapse="")

}

GRanges.cluster.ABC$Compartment = ifelse(GRanges.cluster.ABC$Compartment=="BA","AB", GRanges.cluster.ABC$Compartment)
GRanges.cluster.Rao$Compartment = ifelse(GRanges.cluster.Rao$Compartment=="BA","AB", GRanges.cluster.Rao$Compartment)
GRanges.cluster.DNAse$Compartment = ifelse(GRanges.cluster.DNAse$Compartment=="BA","AB", GRanges.cluster.DNAse$Compartment)


df.Rao.Comp.CRN = data.frame(cbind("Rao",matrix(table(GRanges.cluster.Rao$Compartment), ncol=1), c("AA","AB","BB")))
df.ABC.Comp.CRN = data.frame(cbind("ABC",matrix(table(GRanges.cluster.ABC$Compartment), ncol=1), c("AA","AB","BB")))
df.DNAse.Comp.CRN = data.frame(cbind("DNAse",matrix(table(GRanges.cluster.DNAse$Compartment), ncol=1), c("AA","AB","BB")))

df.all.Comp.CRNs = rbind(df.Rao.Comp.CRN,df.ABC.Comp.CRN,df.DNAse.Comp.CRN)
df.all.Comp.CRNs$X2 = as.numeric(as.character(df.all.Comp.CRNs$X2))
ggplot(df.all.Comp.CRNs, aes(x=X1,y=X2,fill=X3))+ geom_bar(position="fill", stat="identity")+scale_y_continuous(labels = scales::percent_format())+
  ylab("% Overlap") + xlab("Annotation Method") + labs(fill = "Compartment")


agg.elements.byCRN.ABC = merge(aggregate(TargetGene~membership, regul.promoters.ABC, function(x) length(unique(x))),
                               aggregate(name~membership, regul.promoters.ABC, function(x) length(unique(x))), by="membership")


agg.elements.byCRN.ABC$n = agg.elements.byCRN.ABC$TargetGene + agg.elements.byCRN.ABC$name
GRanges.cluster.ABC$n_CRN = agg.elements.byCRN.ABC$n

df.cluster.ABC = data.frame(GRanges.cluster.ABC)
overlap.one.ABC = queryHits(overlaps.ABC.CRNs.Compartments)[!(duplicated(queryHits(overlaps.ABC.CRNs.Compartments))| duplicated(queryHits(overlaps.ABC.CRNs.Compartments), fromLast = TRUE))]

df.cluster.ABC.subset = df.cluster.ABC[overlap.one.ABC,]

summary(aov(n_CRN~Compartment,))
wilcox.test(df.cluster.ABC.subset[df.cluster.ABC.subset$Compartment=="AA", "n_CRN"], df.cluster.ABC.subset[df.cluster.ABC.subset$Compartment=="BB", "n_CRN"])


agg.elements.byCRN.Rao = merge(aggregate(geneSymbol.x~membership, df.Rao, function(x) length(unique(x))),
                               aggregate(name~membership, df.Rao, function(x) length(unique(x))), by="membership")


agg.elements.byCRN.Rao$n = agg.elements.byCRN.Rao$geneSymbol.x + agg.elements.byCRN.Rao$name
GRanges.cluster.Rao$n_CRN = agg.elements.byCRN.Rao$n

df.cluster.Rao = data.frame(GRanges.cluster.Rao)
overlap.one.Rao = queryHits(overlaps.Rao.CRNs.Compartments)[!(duplicated(queryHits(overlaps.Rao.CRNs.Compartments))| duplicated(queryHits(overlaps.Rao.CRNs.Compartments), fromLast = TRUE))]

df.cluster.Rao.subset = df.cluster.Rao[overlap.one.Rao,]

wilcox.test(df.cluster.Rao.subset[df.cluster.Rao.subset$Compartment=="AA", "n_CRN"], df.cluster.Rao.subset[df.cluster.Rao.subset$Compartment=="BB", "n_CRN"])


agg.elements.byCRN.DNAse = merge(aggregate(geneSymbol.x~membership, df.DNase, function(x) length(unique(x))),
                               aggregate(name~membership, df.DNase, function(x) length(unique(x))), by="membership")


agg.elements.byCRN.DNAse$n = agg.elements.byCRN.DNAse$geneSymbol.x + agg.elements.byCRN.DNAse$name
GRanges.cluster.DNAse$n_CRN = agg.elements.byCRN.DNAse$n

df.cluster.DNAse = data.frame(GRanges.cluster.DNAse)
overlap.one.DNAse = queryHits(overlaps.DNAse.CRNs.Compartments)[!(duplicated(queryHits(overlaps.DNAse.CRNs.Compartments))| duplicated(queryHits(overlaps.DNAse.CRNs.Compartments), fromLast = TRUE))]

df.cluster.DNAse.subset = df.cluster.DNAse[overlap.one.DNAse,]

wilcox.test(df.cluster.DNAse.subset[df.cluster.DNAse.subset$Compartment=="AA", "n_CRN"], df.cluster.DNAse.subset[df.cluster.DNAse.subset$Compartment=="BB", "n_CRN"])

df.cluster.ABC.subset$method = "ABC"
df.cluster.Rao.subset$method = "Rao"
df.cluster.DNAse.subset$method = "DNAse"

df.cluster.all = rbind(df.cluster.ABC.subset,df.cluster.Rao.subset,df.cluster.DNAse.subset)

library(ggpubr)

p = ggboxplot(df.cluster.all, x = "Compartment", y = "n_CRN",
          color = "Compartment", palette = "jco",add = "jitter",
          facet.by = "method", short.panel.labs = FALSE)+
  stat_compare_means()
ggpar(p,ylab="N Elements")


ggplot(df.cluster.all, aes(x=method,y=n_CRN, colour=Compartment))+geom_boxplot()
