library(rlist)

process.ABC = function(ABC.inputFile) {
  #Function taking in input the EnhancerPredictions.txt output file from Score-ABC software
  #X chromsome is removed from analysis and gene promoters was added
  
  ABC.input = read.table(ABC.inputFile, header=T)
  #Removing of X chromsome due to specific organization of 3D architecture from Dosage Compensation
  ABC.input.WX = ABC.input[ABC.input$chr != "chrX", ]
  ABC.input.WX$chr = droplevels(ABC.input.WX$chr)
  #Distinction between genic enhancers and intergenic enhancers
  ABC.input.WX$typeOf = sub("\\|.*", "", ABC.input.WX$name)
  
  genes = unique(ABC.input.WX$TargetGene)
  
  #Genes for promoter annotation step
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  genes.hg19 <- genes(txdb)
  symbol <- select(org.Hs.eg.db,keys = genes.hg19$gene_id, columns = c("SYMBOL"), keytype = "ENTREZID")
  genes.hg19$geneSymbol <- symbol$SYMBOL
  #Retriv
  genesInRRCs = genes.hg19[genes.hg19$geneSymbol%in%genes]
  #gene promoters
  promoters = promoters(genesInRRCs)
  
  df.promoters = data.frame(promoters)
  df.promoters = df.promoters[,c("seqnames", "start", "end", "geneSymbol")]
  colnames(df.promoters) = c("chr", "startProm", "endProm", "TargetGene")
  
  enhancers.promoters = merge(ABC.input.WX, df.promoters, by="TargetGene")
  return(enhancers.promoters)
}

Pairs.Enh.Prom.ABC = function(input) {
  #Function takes as input the process.ABC output 
  #The function returns a Pairs of GRanges with enhancers and promoters
  GRanges.enhancers = GRanges(seqnames=input$chr.x, ranges=IRanges(start=input$start, end=input$end, names=input$name), ABC.Score=input$ABC.Score)
  Granges.promoters = GRanges(seqnames=input$chr.y, ranges=IRanges(start=input$startProm, end=input$endProm, names=input$TargetGene), ABC.Score=input$ABC.Score)
  
  return(Pairs(GRanges.enhancers, Granges.promoters))
                              
}

distance.between.Pairs = function(PairObject, absolute=T){
  
  if(absolute){
    dist.prom.enh = abs(start(second(PairObject)) - start(first(PairObject)))
  }
  else{dist.prom.enh = start(second(PairObject)) - start(first(PairObject))}
  return(dist.prom.enh)
  
}

coverage.By.Pair = function(PairObject) {
  start.Pair = sapply(PairObject, function(x) min(start(first(x)), start(second(x))))
  end.Pair = sapply(PairObject, function(x) max(end(first(x)), end(second(x))))
  return(GRanges.Pair.ABC = GRanges(seqnames = seqnames(first(GRanges.Enhancers.Prom.RRCs)), ranges=IRanges(start.Pair, end.Pair, names=paste0("Pair",1:length(GRanges.Enhancers.Prom.RRCs)))))
}

AL1C.Crn.ABC = function(input) { 
  #The function creates CRN based on AL1C definition and ABC-Score annotation
  #The function returns a grpah from igraph package 
  nodes = input[, c("name", "TargetGene")]
  enhancers = unique(nodes$name)
  vertices =  data.frame(rbind(matrix(unique(nodes$name), ncol=1), matrix(unique(nodes$TargetGene), ncol=1)))
  
  colnames(vertices) = c("name")
  vertices$type = vertices[,"name"]%in%enhancers
  vertices$col = ifelse(vertices[,"name"]%in%enhancers, "blue", "red")
  
  graph = graph_from_data_frame(nodes, directed=F, vertices=vertices)
  return(graph)
    
}

Complexity.Crn = function(graph,enhancers,genes ,extract.central.genes=T) {
  
  compo = components(graph)
  subgraphs = decompose(graph)
  
  n.unique.genes = length(unique(genes))
  n.unique.enhancers = length(unique(enhancers))
  
  monogamous.genes = table(compo$csize)[1][[1]] / n.unique.genes
  monogamous.enhancers = table(compo$csize)[1][[1]] / n.unique.enhancers
  
  prom.enhancer.clusters = tapply(names(compo$membership),compo$membership,function(vec,genes) table(vec%in%genes),genes=genes)
  prom.enhancer.clusters.mat = matrix(unlist(prom.enhancer.clusters),length(unlist(prom.enhancer.clusters))/2,2,byrow = T)
  
  OneOneN.genes = sum(prom.enhancer.clusters.mat[prom.enhancer.clusters.mat[,1]==1,2])/n.unique.genes
  OneOneN.enhancers = sum(prom.enhancer.clusters.mat[prom.enhancer.clusters.mat[,2]==1,1])/n.unique.enhancers
  
  if(extract.central.genes){
    #If extract.central.genes the function returns the genes connected with several enhancers and the percentage of genes which are central
    #in network
    
    n.central.genes = lapply(subgraphs, function(x) {
      if(sum(names(V(x)) %in% unique(genes))==1 & length(V(x))>2){
        return(1)
      }
    })
    
    
    central.genes = lapply(subgraphs, function(x) {
      if(sum(names(V(x)) %in% unique(enhancers))==1 & length(V(x))>2){
        return(x)
      }
    })
    
    central.genes = list.clean(central.genes,fun = is.null, recursive = F)
    n.central.genes = sum(unlist(n.central.genes))
    
    return(list("monogamous" = list("genes"= monogamous.genes, "enhancers"=monogamous.enhancers), "11N" = list("genes" = OneOneN.genes, "enhancers" = OneOneN.enhancers), "Centrality" = list("list.genes" = central.genes , "n.genes"= n.central.genes )))
  }
  
  else {
    return(list("monogamous" = list("genes"= monogamous.genes, "enhancers"=monogamous.enhancers), "11N" = list("genes" = OneOneN.genes, "enhancers" = OneOneN.enhancers)))
  }
  
}

.add.membership = function(graph, df, method="ABC"){
  membership = components(graph)$membership
  df.membership = data.frame(membership)
  df.membership$name = rownames(df.membership)
  
  if(method=="ABC"){
    df.membership = df.membership[df.membership$name%in%df$name,]
    df = merge(df, df.membership, by="name")
  }
  
  return(df)
  
}

coverage.By.Crn = function(graph, df, method="ABC") { 
    
    df = .add.membership(graph, df, method)
    
    df$minStart = apply(df[,c("start", "startProm")], 1, min)
    df$maxEnd = apply(df[,c("end", "endProm")], 1, max)

    start.cluster = aggregate(minStart~membership, df, min)
    end.cluster = aggregate(maxEnd~membership, df, max)


    df.start.end = unique(merge(merge(start.cluster, end.cluster, by="membership"), df[,c("chr.x", "membership")], by="membership"))

    GRanges.cluster = GRanges(seqnames = df.start.end$chr, ranges = IRanges(start=df.start.end$minStart, end=df.start.end$maxEnd, names=df.start.end$compo.graph.membership ) )
  
    return(list("clusters"=GRanges.cluster, "data.frame"=df))
}


.make.GRanges.compartments = function(PC.by.locus, resolution=1000000){
  #Function which creates GRanges for Compartments analysis
  #by default the 3 first PC are include
  chroms = paste0("chr",1:22)
  PC.compartments = read.table(PC.by.locus, header=F)
  colnames(PC.compartments) = c("chr", "bin", "PC1", "PC2", "PC3")
  
  PC.compartments.WX = PC.compartments[PC.compartments$chr!="chrX", ]
  
  list.PC.compartments = list()
  
  for(chr in chroms){
    PC.compartments.tmp = PC.compartments.WX[PC.compartments.WX$chr ==chr,]
    
    windows = seq(0, resolution*(nrow(PC.compartments.tmp)), resolution)
    start.stop = data.frame(cbind(windows[-length(windows)], windows[-1]))
    colnames(start.stop) = c("start", "stop")
    
    list.PC.compartments[[chr]] = cbind(PC.compartments.tmp, start.stop)
  }
  
  PC.compartments.analyse = do.call(rbind, list.PC.compartments)
  PC.compartments.analyse[,3:5] = PC.compartments.analyse[,3:5]*-1
  
  GRanges.PC= GRanges(seqnames = PC.compartments.analyse$chr, ranges=IRanges(start = PC.compartments.analyse$start,end =PC.compartments.analyse$stop,names=paste0("bin", 1:nrow(PC.compartments.analyse))),
                           PC1=PC.compartments.analyse$PC1, PC2=PC.compartments.analyse$PC2, PC3=PC.compartments.analyse$PC3)
  
  return(GRanges.PC)
}

define.active.compartments = function(PC.by.locus,resolution=1000000 ,genes, corByChr = TRUE){
  GRanges.PC = .make.GRanges.compartments(PC.by.locus)
  
  GRanges.PC$ngenes = countOverlaps(GRanges.PC, genes)
  
  df.PC = data.frame(GRanges.PC)
  
  corPC1Genes = cor(df.PC$PC1, df.PC$ngenes, method="spearman")
  corPC2Genes = cor(df.PC$PC2, df.PC$ngenes, method="spearman")
  corPC3Genes = cor(df.PC$PC3, df.PC$ngenes, method="spearman")
  
  print(c(corPC1Genes, corPC2Genes,corPC3Genes))
  indexMaxCor = which.max(c(corPC1Genes,corPC2Genes,corPC3Genes))
  
  if(corByChr){
    cor.PC1.genes = ddply(df.PC, .(seqnames), summarise, "corr" = cor(PC1, ngenes, method = "spearman"))
    cor.PC2.genes = ddply(df.PC, .(seqnames), summarise, "corr" = cor(PC2, ngenes, method = "spearman"))
    cor.PC3.genes = ddply(df.PC, .(seqnames), summarise, "corr" = cor(PC3, ngenes, method = "spearman"))
    
    cor.densite.genes.PC = merge(merge(cor.PC1.genes, cor.PC2.genes, by="seqnames"),cor.PC3.genes, by="seqnames")
    colnames(cor.densite.genes.PC) = c("chr", "COR.PC1.GDENSITE", "COR.PC2.GDENSITE", "COR.PC3.GDENSITE")
    
    print(cor.densite.genes.PC)
  }
  
  #density.mediane.genes = median(df.PC$ngenes)
  #GRanges.PC$Compartment = ifelse(GRanges.PC$ngenes>density.mediane.genes, "A", "B")

  mneg = mean(mcols(GRanges.PC)[mcols(GRanges.PC)[,indexMaxCor] < 0, "ngenes"])
  mpos = mean(mcols(GRanges.PC)[mcols(GRanges.PC)[,indexMaxCor] > 0, "ngenes"])
  
  if(mneg > mpos){
    GRanges.PC[,indexMaxCor] =  GRanges.PC[,indexMaxCor] *-1
    GRanges.PC$Compartment = ifelse(mcols(GRanges.PC)[,indexMaxCor]>0, "A", "B")
    return(GRanges.PC[,c(indexMaxCor,5)])
  }
  
  else{
    GRanges.PC$Compartment = ifelse(mcols(GRanges.PC)[,indexMaxCor]>0, "A", "B")
    return(GRanges.PC[,c(indexMaxCor,5)])
  }
}

CRNs.in.Acompartments = function(GRanges.CRNs, GRanges.PC, df){
  Overlap.RRCs.Comp = findOverlaps(GRanges.CRNs, GRanges.PC)
  GRanges.PC$nRRCs = countOverlaps(GRanges.PC, GRanges.CRNs)
  
  query.Overlap.RRCs.Comp = queryHits(Overlap.RRCs.Comp)
  subject.Overlap.RRCs.Comp = subjectHits(Overlap.RRCs.Comp)
  
  list.RRCs.Comp = list()
  
  for(q in unique(query.Overlap.RRCs.Comp)){
  
     genes.enhancers.tmp = df[df$membership==q,c("chr.x","name","start","end","TargetGene", "startProm", "endProm")]
  
     GRanges.prom = GRanges(seqnames=genes.enhancers.tmp$chr.x, ranges=IRanges(start=genes.enhancers.tmp$startProm,end=genes.enhancers.tmp$endProm,names=genes.enhancers.tmp$TargetGene))
     GRanges.enh = GRanges(seqnames=genes.enhancers.tmp$chr.x, ranges=IRanges(start=genes.enhancers.tmp$start,end=genes.enhancers.tmp$end,names=genes.enhancers.tmp$name))
  
     subset.subjects = subjectHits(Overlap.RRCs.Comp[queryHits(Overlap.RRCs.Comp) == q])
  
     GRanges.subset = GRanges.PC[subset.subjects]
  
     GRanges.subset$npromRRCsOverlap = countOverlaps(GRanges.subset,unique(GRanges.prom))
     GRanges.subset$nenhancersRRCsOverlap = countOverlaps(GRanges.subset,unique(GRanges.enh))
  
     list.RRCs.Comp[[q]] = GRanges.subset
   }
  
   proportion.prom.A = sapply(list.RRCs.Comp, function(x){
     if(length(x[x$Compartment=="A"])>0) {
       sum(x[x$Compartment=="A"]$npromRRCsOverlap)/ sum(x$npromRRCsOverlap)}
   }
   )
   proportion.prom.A = proportion.prom.A[!is.nan(proportion.prom.A)]
   
   
   proportion.enhancers.A = sapply(list.RRCs.Comp, function(x){
     if(length(x[x$Compartment=="A"])>0) {
       sum(x[x$Compartment=="A"]$nenhancersRRCsOverlap)/ sum(x$nenhancersRRCsOverlap)}
   }
   )
   proportion.prom.A = proportion.enhancers.A[!is.nan(proportion.enhancers.A)]
   
  return(list("CRNsinA"=(sum(GRanges.PC$nRRCs[GRanges.PC$Compartment=="A"]) / sum(GRanges.PC$nRRCs))*100, "EnhinA" = mean(unlist(proportion.enhancers.A)), "nProminA" = mean(unlist(proportion.prom.A))))

}

annotate.3D.Features = function(GRanges.3D, GRanges.to.annotate, kind=c("DI", "INS", "FIRE"), aggregateFunction = c("mean", "sum")){
  
  Overlaps.3D = findOverlaps(GRanges.to.annotate, GRanges.3D)
  subset.3D = GRanges.3D[subjectHits(Overlaps.3D)]
  
  if(aggregateFunction=="sum"){
    if(kind=="FIRE"){
      mcols(GRanges.to.annotate)[unique(queryHits(Overlaps.3D)), "sum_ScoreFire"] = aggregate(subset.3D$ScoreFire, list(queryHits(Overlaps.3D)), sum, na.rm=T)$x
    }
    else if(kind =="DI"){
      
      mcols(GRanges.to.annotate)[unique(queryHits(Overlaps.3D)), "sum_DI"] = abs(aggregate(subset.3D$DI, list(queryHits(Overlaps.3D)), sum, na.rm=T)$x)
      
    }
    else{
      mcols(GRanges.to.annotate)[unique(queryHits(Overlaps.3D)), "sum_INS"] = aggregate(subset.3D$normalizedScore, list(queryHits(Overlaps.3D)), sum, na.rm=T)$x
      mcols(GRanges.to.annotate)$sum_INS = ifelse(mcols(GRanges.to.annotate)$sum_INS>1, 1, mcols(GRanges.to.annotate)$sum_INS)
    }
  }
  
  else {
    if(kind=="FIRE"){
    mcols(GRanges.to.annotate)[unique(queryHits(Overlaps.3D)), "mean_ScoreFire"] = aggregate(subset.3D$ScoreFire, list(queryHits(Overlaps.3D)), mean, na.rm=T)$x
    }
  else if(kind =="DI"){
    
    mcols(GRanges.to.annotate)[unique(queryHits(Overlaps.3D)), "mean_DI"] = abs(aggregate(subset.3D$DI, list(queryHits(Overlaps.3D)), mean, na.rm=T)$x)
    
  }
  else{
    mcols(GRanges.to.annotate)[unique(queryHits(Overlaps.3D)), "mean_INS"] = aggregate(subset.3D$normalizedScore, list(queryHits(Overlaps.3D)), mean, na.rm=T)$x
    
  }
    
}
  
  
  return(GRanges.to.annotate)
}

annotate.Activity = function(GRanges.Activity, GRanges.to.annotate, kind=c("CTCF", "DNAse", "H3K27ac", "H3K4me1", "p300")){
  
  Overlaps.Activity = findOverlaps(GRanges.to.annotate, GRanges.Activity)
  subset.Activity = GRanges.Activity[subjectHits(Overlaps.Activity)]
  
  mcols(GRanges.to.annotate)[unique(queryHits(Overlaps.Activity)), paste0("mean_", kind)] = aggregate(subset.Activity$signalValue, list(queryHits(Overlaps.Activity)), mean, na.rm=T)$x
  
  return(GRanges.to.annotate)
}

correlationByPairs = function(df, method=c("pearson", "spearman")) {
  combinationOfColumns = combn(colnames(df), 2)
  
  
  PcorByPairs = lapply(1:ncol(combinationOfColumns), function(x) {
    tmp = na.omit(df[,combinationOfColumns[,x]])
    
    return(list("cor"=cor(na.omit(df[,combinationOfColumns[,x]]), method = method),"p.value" = cor.test(tmp[,combinationOfColumns[1,x]], tmp[,combinationOfColumns[2,x]],method=method)$p.value))
  }
  )
  return(PcorByPairs)
  
}

permuteRegionsMetadata <- function(GRangeToPermute, index.col) {
  GRangeToPermute.tmp = GRangeToPermute
  mcols(GRangeToPermute.tmp)[,index.col] <- mcols(GRangeToPermute.tmp)[sample(length(GRangeToPermute.tmp)),index.col]
  return(GRangeToPermute.tmp)
}

PermutationCorrelation = function(df, ActivityCol="", FeaturesCol="") {
  f = formula(paste(ActivityCol,FeaturesCol,sep="~"))
  return(spearman_test(f, df, distribution="approximate"))
}

featuresPerCluster = function(Grange, mcolsToExtract, type=c("enhancer", "promoter")){
  df = data.frame(Grange)
  df = df[,mcolsToExtract]
  colnames(df) = paste0(type,mcolsToExtract)
  return(df)
}  

Pairs.with.Compartments = function(GRanges.Pairs, GRanges.Compartments) {
  
  overlapsPairsCompart = findOverlaps(GRanges.Pairs,GRanges.Compartments)
  d = unique(queryHits(overlapsPairsCompart)[duplicated(queryHits(overlapsPairsCompart))])
  
  non.duplicated.Pairs.Compart = overlapsPairsCompart[!queryHits(overlapsPairsCompart)%in%d]
  duplicated.Pairs.Compart = overlapsPairsCompart[queryHits(overlapsPairsCompart) %in% d]
  
  e = table(mcols(GRanges.Compartments)[subjectHits(non.duplicated.Pairs.Compart),"Compartment"])
  e[3] = e[2]
  e[2] = 0 
  names(e) = c("A", "AB", "B")
  
  multiple.overlaps = lapply(unique(queryHits(duplicated.Pairs.Compart)), function(x) mcols(GRanges.Compartments)[subjectHits(overlapsPairsCompart[queryHits(overlapsPairsCompart)==x]),"Compartment"])
  r = lapply(multiple.overlaps, function(x) x[c(1,length(x))])
  
  g=do.call(rbind,r)
  e1 = table(ifelse(g[,1]=="A" & g[,2] =="A","A", ifelse(g[,1]=="B"&g[,2]=="B","B","AB")))
  
  j = data.frame(rbind(e,e1))
  j = colSums(j)
  return(j)
}


