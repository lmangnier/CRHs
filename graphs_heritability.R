
library(tidyverse)

heri.Enh.ABC = read.table("results_weights_EUR/NEU_ENH_ABC_allchrs_EUR.baseline.results", header=T)
heri.Cand.Regul.ABC = read.table("results_weights_EUR/NEU_CandidateRegul_ABC_allchrs_EUR.baseline.results", header=T)

heri.Enh.Rao = read.table("results_weights_EUR/NEU_ENH_Rao_allchrs_EUR.baseline.results", header=T)
heri.Cand.Regul.Rao = read.table("results_weights_EUR/NEU_CandidateRegul_Rao_allchrs_EUR.baseline.results", header=T)

heri.Enh.DNAse = read.table("results_weights_EUR/NEU_ENH_DNAse_allchrs_EUR.baseline.results", header=T)
heri.Cand.Regul.DNAse = read.table("results_weights_EUR/NEU_CandidateRegul_DNAse_allchrs_EUR.baseline.results", header=T)


heri.G.ABC = read.table("results_weights_EUR/NEU_Genes_ABC_allchrs_EUR.baseline.results", header=T)
heri.Cand.G.ABC = read.table("results_weights_EUR/NEU_CandidatesGenes_ABC_allchrs_EUR.baseline.results", header=T)

heri.G.Rao = read.table("results_weights_EUR/NEU_Genes_Rao_allchrs_EUR.baseline.results", header=T)
heri.Cand.G.Rao = read.table("results_weights_EUR/NEU_CandidatesGenes_Rao_allchrs_EUR.baseline.results", header=T)

heri.G.DNAse = read.table("results_weights_EUR/NEU_Genes_DNAse_allchrs_EUR.baseline.results", header=T)
heri.Cand.G.DNAse = read.table("results_weights_EUR/NEU_CandidatesGenes_DNAse_allchrs_EUR.baseline.results", header=T)


heri.Enh.Cand.all = rbind(heri.Enh.ABC,heri.Cand.Regul.ABC,heri.Enh.Rao,heri.Cand.Regul.Rao,heri.Enh.DNAse,heri.Cand.Regul.DNAse,
                          heri.G.ABC,heri.Cand.G.ABC,heri.G.Rao,heri.Cand.G.Rao,heri.G.DNAse,heri.Cand.G.DNAse)

heri.Enh.Cand.all$method = c("ABC", "ABC","Rao", "Rao","DNAse", "DNAse","ABC", "ABC","Rao", "Rao","DNAse", "DNAse")
heri.Enh.Cand.all$type = c("Regul", "Candidate Regul","Regul", "Candidate Regul","Regul", "Candidate Regul","Prom", "Candidate Prom","Prom", "Candidate Prom","Prom", "Candidate Prom")

write.table(heri.Enh.Cand.all, "heritability_Enh_G_Cand.txt", col.names = T, row.names = F, sep="\t")



heri.SML.Enh.ABC = read.table("results_weights_EUR/SML_Enh_ABC_EUR.baseline.results", header=T)
heri.SML.Enh.ABC$method = "ABC"
heri.SML.Enh.ABC$size = c("small", "medium", "large")
heri.SL.Enh.Rao = read.table("results_weights_EUR/small_Large_Enh_Rao_EUR.baseline.results", header=T)
heri.SL.Enh.Rao$method = "Rao"
heri.SL.Enh.Rao$size = c("small", "large")
heri.SL.Enh.DNAse = read.table("results_weights_EUR/small_Large_Enh_DNAse_EUR.baseline.results", header=T)
heri.SL.Enh.DNAse$method = "DNAse"
heri.SL.Enh.DNAse$size = c("small", "large")

heri.Enh.size.all = rbind(heri.SML.Enh.ABC, heri.SL.Enh.Rao,heri.SL.Enh.DNAse)
write.table(heri.Enh.size.all, "heritability_Enh_size.txt", col.names = T, row.names = F, sep="\t")

heri.SML.G.ABC = read.table("results_weights_EUR/SML_G_ABC_EUR.baseline.results", header=T)
heri.SML.G.ABC$method = "ABC"
heri.SML.G.ABC$size = c("small", "medium", "large")
heri.SL.G.Rao = read.table("results_weights_EUR/small_Large_G_Rao_EUR.baseline.results", header=T)
heri.SL.G.Rao$method = "Rao"
heri.SL.G.Rao$size = c("small", "large")
heri.SL.G.DNAse = read.table("results_weights_EUR/small_Large_G_DNAse_EUR.baseline.results", header=T)
heri.SL.G.DNAse$method = "DNAse"
heri.SL.G.DNAse$size = c("small", "large")

heri.G.size.all = rbind(heri.SML.G.ABC, heri.SL.G.Rao,heri.SL.G.DNAse)
write.table(heri.G.size.all, "heritability_G_size.txt", col.names = T, row.names = F, sep="\t")

heri.Enh.allAnnot = read.table("results_weights_EUR/NEU_ENH_allAnnot_allchrs_EUR.baseline.results", header=T)
heri.Enh.allAnnot$method = c("ABC", "Rao", "DNAse")
write.table(heri.Enh.allAnnot,"heritability_Enh_Adjust_allAnnot.txt", col.names = T, row.names = F)

heri.Enh.G.Adjust_allAnnot = read.table("results_weights_EUR/NEU_ENH_CandidateRegul_allAnnot_allchrs_EUR.baseline.results", header=T)
heri.Enh.G.Adjust_allAnnot$method = c("ABC", "ABC", "Rao", "Rao", "DNAse", "DNAse")
heri.Enh.G.Adjust_allAnnot$type = c("Regul", "Prom","Regul", "Prom","Regul", "Prom")

write.table(heri.Enh.G.Adjust_allAnnot, "heri_Enh_G_Adjust_allAnnot.txt", col.names = T, row.names = F)

ggplot(heri.Enh.Cand.all, aes(x = type, y = log(Enrichment), fill = method)) +
  geom_col(position = "dodge",alpha=0.5) +
  coord_polar()+ theme(text = element_text(size=10))+ylab("Log(Enrichment)") + xlab("Element type") +
  annotate("text", x = 1.7, y = 2.7, label = "***") +annotate("text", x = 2, y = 2.5, label = "***") +annotate("text", x = 2.3, y = 3, label = "***")+
  annotate("text", x = 1, y = 1, label = "***") + annotate("text", x = 1.2, y = 1, label = "***") +annotate("text", x = 0.8, y = 1, label = "***")+
  annotate("text", x = 2.7, y = 1, label = "***") + annotate("text", x = 3, y = 1.3, label = "***") +annotate("text", x = 3.3, y = 1.3, label = "***")+
  annotate("text", x = 4.3, y = 4, label = "***") + annotate("text", x = 4, y = 5.5, label = "***") +annotate("text", x = 3.7, y = 4, label = "***")

ggplot(heri.Enh.size.all, aes(x = size, y = log(Enrichment), fill = method)) +
  geom_col(position = "dodge",alpha=0.5) +
  coord_polar()+ theme(text = element_text(size=10))+ annotate("text", x = 2.7, y = 5, label = "**") +
  annotate("text", x = 3, y = 5.7, label = "****") +
  annotate("text", x = 3.3, y = 5, label = "***") +annotate("text", x = 2, y = 5, label = "***") +annotate("text", x = 0.7, y = 5, label = "***")+
  annotate("text", x = 1, y = 5.7, label = "***") + ylab("Log(Enrichment)") + xlab("CRN Size")

ggplot(heri.G.size.all, aes(x = size, y = log(Enrichment), fill = method)) +
  geom_col(position = "dodge",alpha=0.5) +
  coord_polar()+
  theme(text = element_text(size=10))+ annotate("text", x = 3, y = 1, label = "**") + annotate("text", x = 3.3, y = 1, label = "***")+
  annotate("text", x = 0.7, y = 1, label = "*") +annotate("text", x = 2, y = 1, label = "***") +
  annotate("text", x = 1, y = 1, label = "***") + ylab("Log(Enrichment)") + xlab("CRN Size")
