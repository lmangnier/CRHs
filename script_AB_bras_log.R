library(rCGH)
# ou load("hg19.RData")
p = paste(paste0("chr",hg19$chrom),1,hg19$centromerStart,sep="_")
q = paste(paste0("chr",hg19$chrom),hg19$centromerEnd+1,hg19$length,sep="_")

resol = c("100Kb","500Kb","1Mb")
resoln = c(1e5,5e5,1e6)
for (i in 1:22)
{
cat(i,"\n")
for (j in 1:length(resol))
{
oe.dense.p = read.table(paste0("output_HiC_dense_KR/",resol[j],"/",p[i],"_",resol[j],".txt"), header=F)
oe.dense.q = read.table(paste0("output_HiC_dense_KR/",resol[j],"/",q[i],"_",resol[j],".txt"), header=F)
chr = paste0("chr",i)
dim(oe.dense.p)
dim(oe.dense.q)

# Suppression des valeurs manquantes
if (any(is.na(oe.dense.p)))
{
	cat("Valeurs manquantes p\n")
	oe.dense.p[is.na(oe.dense.p)] = 0
}
if (any(is.na(oe.dense.q)))
{
	cat("Valeurs manquantes q\n")
	oe.dense.q[is.na(oe.dense.q)] = 0
}

#Suppression des colonnes pour lesquels l'ecart-type est nul
pos.p = (1:nrow(oe.dense.p))[-which(apply(oe.dense.p, 2,function(x) sd(x)==0))]
pos.q = (1:nrow(oe.dense.q))[-which(apply(oe.dense.q, 2,function(x) sd(x)==0))]
oe.dense.p = oe.dense.p[,pos.p]
oe.dense.q = oe.dense.q[,pos.q]
dim(oe.dense.p)
dim(oe.dense.q)
#Suppression des lignes pour lesquels l'ecart-type est nul
oe.dense.p = oe.dense.p[pos.p,]
oe.dense.q = oe.dense.q[pos.q,]

# Bras p
if (all(dim(oe.dense.p))>0)
{
#Visualisation de la matrice observee sur attendue
# jpeg(paste0("obs_expect/",chr,"p_obs_expect",resol[j],".jpg"))
# heatmap(data.matrix(oe.dense.p),Rowv=NA,Colv=NA)
# dev.off()


# Analyse en composantes principales
PC.logoe.p = prcomp(log(oe.dense.p+0.0001))

jpeg(paste0("PCplots/",chr,"p_log_E1",resol[j],".jpg"))
plot(PC.logoe.p$rotation[,1],type="l",ylab="E1")
abline(h=0)
dev.off()
#Exportation des 3 premieres composantes pour les analyses en lien avec la densite de genes pour choisir la 1CP (Gorkin et al., 2019)
write.table(cbind(chr,pos.p,PC.logoe.p$rotation[,1:3]), paste0("PC_logoe/",chr,"p_",resol[j],"_PC_logoe.txt"), sep="\t",col.names=F,row.names=F,quote=F)
}

# Bras q
if(all(dim(oe.dense.q)>0))
{
#Visualisation de la matrice observee sur attendue
# jpeg(paste0("obs_expect/",chr,"q_obs_expect",resol[j],".jpg"))
# heatmap(data.matrix(oe.dense.q),Rowv=NA,Colv=NA)
# dev.off()

# Analyse en composantes principales
PC.logoe.q = prcomp(log(oe.dense.q+0.0001))

jpeg(paste0("PCplots/",chr,"q_log_E1",resol[j],".jpg"))
plot(PC.logoe.q$rotation[,1],type="l",ylab="E1")
abline(h=0)
dev.off()

#Exportation des 3 premieres composantes pour les analyses en lien avec la densite de genes pour choisir la 1CP (Gorkin et al., 2019)
write.table(cbind(chr,pos.q,PC.logoe.q$rotation[,1:3]), paste0("PC_logoe/",chr,"q_",resol[j],"_PC_logoe.txt"), sep="\t",col.names=F,row.names=F,quote=F)
}
}
}
# Ensuite: cat chr*_100Kb_PC_logoe.txt > allchrarms_PC_logoe.100Kb.txt
# cat chr*_500Kb_PC_logoe.txt > allchrarms_PC_logoe.500Kb.txt
