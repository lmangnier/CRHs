#Additional code non used in the main core code 

#Whole file of ABC-Score: 
#Fichier des RRCs sur la base du score-ABC contenant l'ensemble des contacts genes-enhancers: will be used for enrichment questions 
all.contacts = read.table("EnhancerPredictionsAllPutative.txt", header=T)
all.contacts.WX = all.contacts[all.contacts$chr!="chrX", ]
all.contacts.WX$chr = droplevels(all.contacts.WX$chr)

all.contacts.nRRCs = all.contacts.WX[all.contacts.WX$ABC.Score<0.012,c("chr", "start", "end", "name", "class", "TargetGene")]

rm(all.contacts)

#Recuperation des enhancers non inclus dans les RRCs
#Fulco et al., 2019 ont demontre que le Score-ABC ne performe pas bien pour les contacts promoters-promoters
enhancers.nRRCs = all.contacts.nRRCs[!(all.contacts.nRRCs$name%in%enhancers)&(all.contacts.nRRCs$class!="promoter"),]
enhancers.nRRCs = na.omit(enhancers.nRRCs)

enhancers.promoters.nRRCs = merge(enhancers.nRRCs, df.promoters, by="TargetGene")

GRanges.Promoters.nRRCs = GRanges(seqnames=enhancers.promoters.nRRCs$chr.y, ranges=IRanges(start=enhancers.promoters.nRRCs$startProm, end=enhancers.promoters.nRRCs$endProm))
GRanges.Enhancers.nRRCs = GRanges(seqnames=enhancers.promoters.nRRCs$chr.x, ranges=IRanges(start=enhancers.promoters.nRRCs$start, end=enhancers.promoters.nRRCs$end))


GRanges.Enhancers.Prom.nRRCs = Pairs(GRanges.Enhancers.nRRCs,GRanges.Promoters.nRRCs)

#Promoters and Enhancers are different based on association behaviour ? 
#Les enhancers presentent-ils un profil d'association similaire aux genes ?
asso.enhancers = as.numeric(table(EnhancersPred.WX$name))
asso.genes = as.numeric(table(EnhancersPred.WX$TargetGene))

summary(asso.enhancers)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   1.000   1.000   2.198   2.000  71.000

summary(asso.genes)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   2.000   4.000   4.386   6.000  26.000 

sd(asso.enhancers); sd(asso.genes)

e = data.frame(cbind("enhancers", asso.enhancers))
g = data.frame(cbind("genes", asso.genes))
colnames(e) = c("type", "n")
colnames(g) = colnames(e)
eg = rbind(e,g) 
eg$n = as.numeric(as.character(eg$n))

#Non-parametric counterpart to know if we find the same conclusions
anova(lm(n~type, eg)) #Average difference between the two groups
kruskal.test(n~type,eg) #Difference betwenn mean rank 

#Are promoters more associated than enhancers ?  
#Effect direction:
t.test(asso.genes,asso.enhancers, alternative = "greater")
#Promoters are more connected than enhancers when we consider ABC-Score


#Are Enhancers included in CRNs enriched in FIREs?
GRanges.enhancers$FIREs.nsigni = countOverlaps(GRanges.enhancers, GRanges.allFIREs[GRanges.allFIREs$pvalues>0.05])
GRanges.enhancers$FIREs.signi = countOverlaps(GRanges.enhancers, GRanges.allFIREs[GRanges.allFIREs$pvalues<=0.05])

GRanges.enhancers.nRRCs$FIREs.nsigni = countOverlaps(GRanges.enhancers.nRRCs, GRanges.allFIREs[GRanges.allFIREs$pvalues>0.05])
GRanges.enhancers.nRRCs$FIREs.signi = countOverlaps(GRanges.enhancers.nRRCs, GRanges.allFIREs[GRanges.allFIREs$pvalues<=0.05])

table(GRanges.enhancers.nRRCs$FIREs.signi)
enrichissement.FIREs = matrix(c(sum(GRanges.enhancers$FIREs.signi), sum(GRanges.enhancers.nRRCs$FIREs.signi), sum(GRanges.enhancers$FIREs.nsigni), sum(GRanges.enhancers.nRRCs$FIREs.nsigni)), 
                              ncol=2, nrow=2, byrow=T)

fisher.test(enrichissement.FIREs, alternative = "greater")

#Disctinction between genic and intergenic enhancers
enrichissement.geniqueVSintergenique = matrix(c(sum(GRanges.enhancers[GRanges.enhancers$typeOf=="intergenic"]$FIREs.signi), sum(GRanges.enhancers[GRanges.enhancers$typeOf=="genic"]$FIREs.signi), sum(GRanges.enhancers[GRanges.enhancers$typeOf=="intergenic"]$FIREs.nsigni), sum(GRanges.enhancers[GRanges.enhancers$typeOf=="genic"]$FIREs.nsigni)), 
                                              ncol=2, nrow=2, byrow=T)

(enrichissement.geniqueVSintergenique[1,1] * enrichissement.geniqueVSintergenique[2,2]) / (enrichissement.geniqueVSintergenique[1,2] * enrichissement.geniqueVSintergenique[2,1])
#1.14471

fisher.test(enrichissement.geniqueVSintergenique, alternative = "greater")
#Profile of enhancers seems to not impact FIRE enrichment


