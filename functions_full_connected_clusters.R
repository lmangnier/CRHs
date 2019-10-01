#Functions of "Full-linked" cluster based definition

test.overlapping <- function(TSS.gene,TES.gene, enhS, enhE) {
  #Overlapping analysis function
  #This function tests if an enhancer is dowstream, upstream or overlapping his target gene
  if(TSS.gene > enhE){
    
    return("downstream")
  }
  else if (TES.gene < enhS){
    
    return("upstream")
  }
  else{
    
    return("overlapping")
  }
}

contact.func <- function(data, col.to.agg){
  #This function takes as inputs a DataFrame and a column to aggregate
  #The function returns a qualitative and quantitative results based on "full-linked" cluster
  #This definition allows that each element of a cluster has to be linked each other
  library(data.table)
  tmp <- copy(data)
  
  #If column which wants to aggregate has a type list, one more conversion step is needed
  tmp$col1.bis <- as.character(lapply(tmp[,col.to.agg], FUN = function(x) paste(x, collapse = "-")))
  
  if(col.to.agg == "geneSymbol") {
    tmp.agg.qual <- aggregate(enhancer ~ col1.bis, data = tmp, length)
    tmp.agg.qual$N.genes <- as.numeric(lapply(strsplit(tmp.agg.qual$col1.bis, "-"), length))
    tmp.agg.quant <- as.data.frame(table(tmp.agg.qual$N.genes, tmp.agg.qual$enhancer))
    colnames(tmp.agg.quant) <- c("N.Genes", "N.Enhancers", "N")
    
    tmp.agg.quant <- tmp.agg.quant[!tmp.agg.quant$N == 0,]
    
  }
  
  else{
    tmp.agg.qual <- aggregate(geneSymbol ~ col1.bis, data = tmp, length)
    tmp.agg.qual$N.genes <- as.numeric(lapply(strsplit(tmp.agg.qual$col1.bis, "-"), length))
    tmp.agg.quant <- as.data.frame(table(tmp.agg.qual$N.genes, tmp.agg.qual$geneSymbol))
    colnames(tmp.agg.quant) <- c("N.Enhancers", "N.Genes", "N")
    
    tmp.agg.quant <- tmp.agg.quant[!tmp.agg.quant$N == 0,]
  }
  
  
  list("qualitative_asso"= tmp.agg.qual, "quantitative_asso"= tmp.agg.quant)
}

#"Full-connected cluster" definition on FANTOM5 enhancer methodology: 
#The code is not executed
  
#agg.table.enhancers <- aggregate(geneSymbol~enhancer, data=genes.enhancers.contact, unique, na.rm=TRUE)
#agg.table.genes <- aggregate(enhancer~geneSymbol, data=genes.enhancers.contact, unique, na.rm=TRUE)
  
#table(sapply(agg.table.enhancers$geneSymbol, length))
#table(sapply(agg.table.genes$enhancer, length))
  
#barplot(table(sapply(agg.table.enhancers$geneSymbol, length)), main = "Distribution of number of genes associated by enhancer", xlab = "Number of genes associated", ylab = "Frequency")
#barplot(table(sapply(agg.table.genes$enhancer, length)), main = "Distribution of number of enhancers associated by genes", xlab = "Number of enhancers associated", ylab = "Frequency")
  
#genes.enhancers.contact$overlapping <- mapply(test.overlapping, genes.enhancers.contact$TSS, genes.enhancers.contact$TES, genes.enhancers.contact$enhancerstart, genes.enhancers.contact$enhancerstop)

#barplot(table(genes.enhancers.contact$overlapping), main = "Barplot of enhancer position in relation to his gene")

#Quantile Coverage analysis by overlapping status 
#genes.enhancers.contact$width.enh <- genes.enhancers.contact$enhancerstop - genes.enhancers.contact$enhancerstart
#quant.overlapping <- aggregate(width.enh ~ overlapping, data=genes.enhancers.contact, FUN = "quantile" ,probs=c(0.25, 0.50,0.75))
#quant.overlapping

#analysis.g_e <- contact.func("geneSymbol",agg.table.enhancers)
#head(analysis.g_e$quantitative_asso)

#analysis.e_g <- contact.func("enhancer", agg.table.genes)
#head(analysis.e_g$quantitative_asso)


#"Full-connected cluster" definition on Epigenomic Roadmap enhancer methodology: 
#agg.table.enhancers.epg <- aggregate(geneSymbol~enhancer, data=genes.enhancers.EGRM, unique, na.rm=TRUE)
#agg.table.genes.epg <- aggregate(enhancer~geneSymbol, data=genes.enhancers.EGRM, unique, na.rm=TRUE)

#table(sapply(agg.table.enhancers.epg$geneSymbol, length))
#table(sapply(agg.table.genes.epg$enhancer, length))

#barplot(table(sapply(agg.table.enhancers.epg$geneSymbol, length)), main = "Distribution of number of genes associated by enhancer", xlab = "Number of genes associated", ylab = "Frequency")
#barplot(table(sapply(agg.table.genes.epg$enhancer, length)), main = "Distribution of number of enhancers associated by genes", xlab = "Number of enhancers associated", ylab = "Frequency")

#genes.enhancers.EGRM$overlapping <- mapply(test.overlapping, genes.enhancers.EGRM$TSS, genes.enhancers.EGRM$TES, genes.enhancers.EGRM$enhancerstart, genes.enhancers.EGRM$enhancerstop)

#barplot(table(genes.enhancers.EGRM$overlapping), main = "Barplot of enhancer position in relation to his gene")

#Quantile Coverage analysis by overlapping status 
#genes.enhancers.EGRM$width.enh <- genes.enhancers.EGRM$enhancerstop - genes.enhancers.EGRM$enhancerstart
#quant.overlapping <- aggregate(width.enh ~ overlapping, data=dgenes.enhancers.EGRM FUN = "quantile" ,probs=c(0.25, 0.50,0.75))
#quant.overlapping

#analysis.g_e.epg <- contact.func("geneSymbol",agg.table.enhancers.epg)
#head(analysis.g_e.epg$quantitative_asso)

#analysis.e_g.epg <- contact.func("enhancer", agg.table.genes.epg)
#head(analysis.e_g.epg$quantitative_asso)
