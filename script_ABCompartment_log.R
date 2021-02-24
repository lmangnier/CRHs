#args = commandArgs(trailingOnly = TRUE)

#Fichiers dump produits par Juicer sur une resolution de 1Mb et 500Kb
#Observe/Attendu + normalisation KR
#oe.dense.500Kb = read.table(args[1], header=F)
#oe.dense.1Mb = read.table(args[2], header=F)
#chr = args[3]

for (i in 1:22)
{
cat(i,"\n")
oe.dense.500Kb = read.table(paste0("observed_expected_matrix/chr",i,"/dense_500Kb.matrix"), header=F)
oe.dense.1Mb = read.table(paste0("observed_expected_matrix/chr",i,"/dense_1Mb.matrix"), header=F)
chr = paste0("chr",i)
dim(oe.dense.500Kb)
dim(oe.dense.1Mb)

# Suppression des valeurs manquantes
if (any(is.na(oe.dense.500Kb)))
{
	cat("Valeurs manquantes 500kb\n")
	oe.dense.500Kb[is.na(oe.dense.500Kb)] = 0
}
if (any(is.na(oe.dense.1Mb)))
{
	cat("Valeurs manquantes 1Mb\n")
	oe.dense.1Mb[is.na(oe.dense.1Mb)] = 0
}

#Suppression des colonnes pour lesquels l'ecart-type est nul
pos.500Kb = (1:nrow(oe.dense.500Kb))[-which(apply(oe.dense.500Kb, 2,function(x) sd(x)==0))]
pos.1Mb = (1:nrow(oe.dense.1Mb))[-which(apply(oe.dense.1Mb, 2,function(x) sd(x)==0))]
oe.dense.500Kb = oe.dense.500Kb[,pos.500Kb]
oe.dense.1Mb = oe.dense.1Mb[,pos.1Mb]
dim(oe.dense.500Kb)
dim(oe.dense.1Mb)
#Suppression des lignes pour lesquels l'ecart-type est nul
oe.dense.500Kb = oe.dense.500Kb[pos.500Kb,]
oe.dense.1Mb = oe.dense.1Mb[pos.1Mb,]

#Visualisation de la matrice observee sur attendue
# jpeg(paste0("obs_expect/",chr,"_obs_expect_1Mb.jpg"))
# heatmap(data.matrix(oe.dense.1Mb),Rowv=NA,Colv=NA)
# dev.off()

# jpeg(paste0("obs_expect/",chr,"_obs_expect_500Kb.jpg"))
# heatmap(data.matrix(oe.dense.500Kb),Rowv=NA,Colv=NA)
# dev.off()

#Calcul de la matrice de correlation--> correlation pour chaque paire de cellules
#Intuitivement la correlation va nous donner l'association d'une paire de cellules en lien avec l'ensemble des autres cellules du chromosome
#Presentent-elles des valeurs proches pour les memes cellules pour lesquelles elles sont associees?
# corMap.500Kb=cor(oe.dense.500Kb)
# corMap.500Kb[is.na(corMap.500Kb)] = 0
# dim(corMap.500Kb)

# corMap.1Mb = cor(oe.dense.1Mb)
# dim(corMap.1Mb)

#Visualisation de la correlation intra-chromosomale
# jpeg(paste0(chr,"_corrMap_500Kb.jpg"))
# heatmap(corMap.500Kb)
# dev.off()

# jpeg(paste0(chr,"_corrMap_1Mb.jpg"))
# heatmap(corMap.1Mb)
# dev.off()

#PC.500Kb = prcomp(corMap.500Kb)

# On plafonne les valeurs aux 99e percentile et on applique un plancher au 1e percentile
# suivant ce que Luis Chumpitaz a fait
#q99 = quantile(unlist(oe.dense.1Mb),0.99)
#q01 = quantile(unlist(oe.dense.1Mb),0.01)
#oe.dense.clipped.1Mb = pmax(pmin(as.matrix(oe.dense.1Mb),q99),q01)


#q99 = quantile(unlist(oe.dense.500Kb),0.99)
#q01 = quantile(unlist(oe.dense.500Kb),0.01)
#oe.dense.clipped.500Kb = pmax(pmin(as.matrix(oe.dense.500Kb),q99),q01)

# Analyse en composantes principales
# PC.oe.500Kb = prcomp(oe.dense.clipped.500Kb-1)
# PC.oe.1Mb = prcomp(oe.dense.clipped.1Mb-1)

# jpeg(paste0("PCplots/",chr,"_E1_500Kb.jpg"))
# plot(PC.oe.500Kb$rotation[,1],type="l",ylab="E1")
# if (!any(is.na(corMap.500Kb))) lines(PC.500Kb$rotation[,1],col="red")
# legend(10,-0.1,c("O/E - 1","correlation"),lty=1,col=c("black","red"),bty="n")
# abline(h=0)
# dev.off()

# #Exportation des 3 premieres composantes pour les analyses en lien avec la densite de genes pour choisir la 1CP (Gorkin et al., 2019)
# write.table(cbind(chr,1:nrow(PC.500Kb$rotation),PC.500Kb$rotation[,1:3]), paste0("PC_oe/",chr,"_PC_oe.500Kb.txt"), sep="\t",col.names=F,row.names=F,quote=F)
# write.table(cbind(chr,1:nrow(PC.1Mb$rotation),PC.1Mb$rotation[,1:3]), paste0("PC_oe/",chr,"_PC_oe.1Mb.txt"), sep="\t",col.names=F,row.names=F,quote=F)

PC.logoe.500Kb = prcomp(log(oe.dense.500Kb+0.0001))
PC.logoe.1Mb = prcomp(log(oe.dense.1Mb+0.0001))

# jpeg(paste0("PCplots/",chr,"_log_E1_500Kb.jpg"))
# plot(PC.logoe.500Kb$rotation[,1],type="l",ylab="E1")
# lines(PC.500Kb$rotation[,1],col="red")
# legend(10,-0.1,c("log(O/E)","correlation"),lty=1,col=c("black","red"),bty="n")
# abline(h=0)
# dev.off()

#Exportation des 3 premieres composantes pour les analyses en lien avec la densite de genes pour choisir la 1CP (Gorkin et al., 2019)
write.table(cbind(chr,pos.500Kb,PC.logoe.500Kb$rotation[,1:3]), paste0("PC_logoe/",chr,"_PC_logoe.500Kb.txt"), sep="\t",col.names=F,row.names=F,quote=F)
write.table(cbind(chr,pos.1Mb,PC.logoe.1Mb$rotation[,1:3]), paste0("PC_logoe/",chr,"_PC_logoe.1Mb.txt"), sep="\t",col.names=F,row.names=F,quote=F)

}

# Ensuite: cat chr*_PC_logoe.500Kb.txt > allchrs_PC_logoe.500Kb.txt
