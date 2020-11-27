args = commandArgs(trailingOnly = TRUE)

#Fichiers dump produits par Juicer sur une resolution de 1Mb et 500Kb
#Observe/Attendu + normalisation KR
oe.dense.500Kb = read.table(args[1], header=F)
oe.dense.1Mb = read.table(args[2], header=F)
chr = args[3]

#Suppression des colonnes pour lesquels l'ecart-type est nul
oe.dense.500Kb = oe.dense.500Kb[,-which(apply(oe.dense.500Kb, 2,function(x) sd(x)==0))]
dim(oe.dense.500Kb)
oe.dense.1Mb = oe.dense.1Mb[,-which(apply(oe.dense.1Mb ,2, function(x) sd(x)==0))]
dim(oe.dense.1Mb)

#Visualisation de la matrice observee sur attendue
jpeg(paste0(chr,"_obs_expect_1Mb.jpg"))
heatmap(data.matrix(oe.dense.1Mb))
dev.off()

jpeg(paste0(chr,"_obs_expect_500Kb.jpg"))
heatmap(data.matrix(oe.dense.500Kb))
dev.off()

#Calcul de la matrice de correlation--> correlation pour chaque paire de cellules
#Intuitivement la correlation va nous donner l'association d'une paire de cellules en lien avec l'ensemble des autres cellules du chromosome
#Presentent-elles des valeurs proches pour les memes cellules pour lesquelles elles sont associees?
corMap.500Kb=cor(oe.dense.500Kb)
corMap.500Kb[is.na(corMap.500Kb)] = 0
dim(corMap.500Kb)

corMap.1Mb = cor(oe.dense.1Mb)
corMap.1Mb[is.na(corMap.1Mb)] = 0
dim(corMap.1Mb)

#Visualisation de la correlation intra-chromosomale
jpeg(paste0(chr,"_corrMap_500Kb.jpg"))
heatmap(corMap.500Kb)
dev.off()

jpeg(paste0(chr,"_corrMap_1Mb.jpg"))
heatmap(corMap.1Mb)
dev.off()

PC.500Kb = prcomp(corMap.500Kb)
PC.1Mb = prcomp(corMap.1Mb)

#Exportation des 3 premieres composantes pour les analyses en lien avec la densite de genes pour choisir la 1CP (Gorkin et al., 2019)
write.table(PC.500Kb$rotation[,1:3], paste0(chr,"_PC.500Kb.txt"), sep="\t")
write.table(PC.1Mb$rotation[,1:3], paste0(chr,"_PC.1Mb.txt"), sep="\t")
