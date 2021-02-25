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
  #Retrieve only genes with known promoters
  genesInRRCs = genes.hg19[genes.hg19$geneSymbol%in%genes]
  #gene promoters
  promoters = promoters(genesInRRCs)
  
  df.promoters = data.frame(promoters)
  df.promoters = df.promoters[,c("seqnames", "start", "end", "geneSymbol")]
  colnames(df.promoters) = c("chr", "startProm", "endProm", "TargetGene")
  
  enhancers.promoters = merge(ABC.input.WX, df.promoters, by="TargetGene")
  return(enhancers.promoters)
}

process.Rao.DNAse = function(positive.contact, DNAse.file=NULL, method=c("Rao", "DNAse")){
  #This function create dataframe for Rao or DNAse annotation method
  
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  genes.hg19 <- genes(txdb)
  symbol <- select(org.Hs.eg.db,keys = genes.hg19$gene_id, columns = c("SYMBOL"), keytype = "ENTREZID")
  genes.hg19$geneSymbol <- symbol$SYMBOL
  
  promoters = promoters(genes.hg19)
  promoters = promoters[seqnames(promoters)%in%paste0("chr",1:22)]
  
  #Genes overlaping first bin
  overlaps.Prom.Bin1 = findOverlaps(promoters, first(positive.contact))
  #Genes overlaping second bin
  overlaps.Prom.Bin2 = findOverlaps(promoters, second(positive.contact))
  
  if(method=="Rao"){
    
    promoters.Bin1 = promoters[queryHits(overlaps.Prom.Bin1)]
    
    #Bins for which the first bin has a gene
    regulatory.elements.bin2 = second(positive.contact)[subjectHits(overlaps.Prom.Bin1)]
    
    Pairs.prom.regulatory.1 = Pairs(promoters.Bin1, regulatory.elements.bin2)
    Pairs.prom.regulatory.1 = Pairs.prom.regulatory.1[!duplicated(Pairs.prom.regulatory.1)]
    
    promoters.Bin2 = promoters[queryHits(overlaps.Prom.Bin2)]
    
    #Bins for which the second bin has a gene
    regulatory.elements.bin1 = first(positive.contact)[subjectHits(overlaps.Prom.Bin2)]
    
    Pairs.prom.regulatory.2 = Pairs(promoters.Bin2, regulatory.elements.bin1)
    Pairs.prom.regulatory.2 = Pairs.prom.regulatory.2[!duplicated(Pairs.prom.regulatory.2)]
    
    grList.Rao = GRangesList(c(first(Pairs.prom.regulatory.1), first(Pairs.prom.regulatory.2)), c(second(Pairs.prom.regulatory.1), second(Pairs.prom.regulatory.2)))
    Pairs.prom.regulatory.Rao = unique(Pairs(grList.Rao[[1]], grList.Rao[[2]]))
    return(Pairs.prom.regulatory.Rao)
  }
  
  else{
    overlaps.peaks.Bin1 = findOverlaps(DNAse.file, first(positive.contact))
    overlaps.peaks.Bin2 = findOverlaps(DNAse.file, second(positive.contact))
    
    peaksToProm = merge(data.frame(overlaps.peaks.Bin1), data.frame(overlaps.Prom.Bin2), by="subjectHits")
    Promtopeaks = merge(data.frame(overlaps.peaks.Bin2), data.frame(overlaps.Prom.Bin1), by="subjectHits")
    
    Pairs.peak.prom.1 = Pairs(DNAse.file[peaksToProm$queryHits.x],promoters[peaksToProm$queryHits.y])
    Pairs.peak.prom.2 = Pairs(DNAse.file[Promtopeaks$queryHits.x],promoters[Promtopeaks$queryHits.y])
    
    gl.Peaks = GRangesList(c(second(Pairs.peak.prom.1), second(Pairs.peak.prom.2)),c(first(Pairs.peak.prom.1), first(Pairs.peak.prom.2)))
    Pairs.prom.regulatory.DNAse = Pairs(gl.Peaks[[1]], gl.Peaks[[2]])
    
    # mcols(first(Pairs.prom.regulatory.DNAse))[,c("name", "score", "signalValue", "pValue", "qValue", "peak")] = NULL
    # mcols(second(Pairs.prom.regulatory.DNAse))[,c("name", "score", "signalValue", "pValue", "qValue", "peak")] = NULL
    
    return(Pairs.prom.regulatory.DNAse)
  }
}

annotate.dataframe = function(Pair.object, method=c("ABC", "Rao", "DNAse"), df.ABC = NULL){
  #Function takes as input a GRanges Pair annotated for 3D features and epigenetic activity
  
  first.element = first(Pair.object)
  second.element = second(Pair.object)
  
  
  if(method=="ABC"){
    df.Promoters.RRCs = data.frame(unique(second.element))
    df.Promoters.RRCs = df.Promoters.RRCs[,c("seqnames", "start", "end","ABC.Score" ,"mean_DI", "mean_ScoreFire","mean_INS", "mean_CTCF", "mean_H3K27ac", "mean_DNAse","mean_H3K4me1","mean_H3K4me3", "mean_p300", "mean_H3K27me3")]
    colnames(df.Promoters.RRCs) = c("chr.x", "startProm", "endProm","ABC.Score.PROM","DI.PROM", "ScoreFire.PROM", "INS.PROM","mean_CTCF.PROM", "mean_H3K27ac.PROM","mean_DNAse.PROM", "mean_H3K4me1.PROM", "mean_H3K4me3.PROM","mean_p300.PROM", "mean_H3K27me3.PROM")
    
    
    df.Enhancers.RRCs = data.frame(unique(first.element))
    df.Enhancers.RRCs = df.Enhancers.RRCs[,c("seqnames", "start", "end","ABC.Score" ,"mean_DI", "mean_ScoreFire","mean_INS", "mean_CTCF", "mean_H3K27ac", "mean_DNAse","mean_H3K4me1","mean_H3K4me3", "mean_p300", "mean_H3K27me3")]
    colnames(df.Enhancers.RRCs) = c("chr.y", "start", "end","ABC.Score.ENH","DI.ENH", "ScoreFire.ENH", "INS.ENH","mean_CTCF.ENH", "mean_H3K27ac.ENH","mean_DNAse.ENH", "mean_H3K4me1.ENH", "mean_H3K4me3.ENH","mean_p300.ENH", "mean_H3K27me3.ENH")
  
    df.ABC = merge(df.ABC, df.Promoters.RRCs, by=c("chr.x", "startProm", "endProm"))
    df.ABC = merge(df.ABC, df.Enhancers.RRCs, by=c("chr.y", "start", "end"))
    return(df.ABC)
  }
  
  
  else{
    df.concat = cbind(data.frame(first.element), data.frame(second.element))
    
    colnames(df.concat) = c("chr.x", "start.x", "end.x", "width", "strand", "gene_id.x","geneSymbol.x","mean_DI.x","mean_ScoreFire.x","mean_INS.x","mean_CTCF.x",
                            "mean_H3K27ac.x", "mean_DNAse.x","mean_H3K4me1.x", "mean_H3K4me3.x","mean_p300.x","mean_H3K27me3.x","chr.y", "start.y", "end.y", "width", "strand", "gene_id.y", "geneSymbol.y",
                            "mean_DI.y","mean_ScoreFire.y","mean_INS.y","mean_CTCF.y","mean_H3K27ac.y", "mean_DNAse.y","mean_H3K4me1.y", "mean_H3K4me3.y","mean_p300.y","mean_H3K27me3.y")
    
    df.concat$name = gsub(" ", "",apply(df.concat[,c("chr.y","start.y", "end.y")], 1, paste0, collapse=":"))
    
    return(df.concat)
  }
  
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
  return(GRanges.Pair.ABC = GRanges(seqnames = seqnames(first(PairObject)), ranges=IRanges(start.Pair, end.Pair, names=paste0("Pair",1:length(PairObject)))))
}

create.Crn = function(input, method=c("ABC", "Rao", "DNAse")) { 
  
  if(method =="ABC"){
    #The function creates CRN based on AL1C definition and ABC-Score annotation
    #The function returns a graph from igraph package 
    nodes = input[, c("name", "TargetGene")]
    enhancers = unique(nodes$name)
    vertices =  data.frame(rbind(matrix(unique(nodes$name), ncol=1), matrix(unique(nodes$TargetGene), ncol=1)))
    
    colnames(vertices) = c("name")
    vertices$type = vertices[,"name"]%in%enhancers
    vertices$col = ifelse(vertices[,"name"]%in%enhancers, "blue", "red")
    
    graph = graph_from_data_frame(nodes, directed=F, vertices=vertices)
    return(graph)
  }
  
  else {
    nodes = input[,c("geneSymbol.x","name")]
    nodes = na.omit(nodes)
    vertices = unique(rbind(matrix(nodes$geneSymbol.x, ncol=1), matrix(nodes$name, ncol=1)))
    graph = graph_from_data_frame(nodes, directed = F, vertices=vertices)
    return(graph)
  }
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

add.membership = function(graph, df){
  membership = components(graph)$membership
  df.membership = data.frame(membership)
  df.membership$name = rownames(df.membership)
  
  df.membership = df.membership[df.membership$name%in%df$name,]
  df = merge(df, df.membership, by="name")
  
  return(df)
  
}

coverage.By.Crn = function(graph, df, method=c("ABC", "Rao", "DNAse")) { 
    
    #df = add.membership(graph, df)
    
  if(method=="ABC"){
    df$minStart = apply(df[,c("start", "startProm")], 1, min)
    df$maxEnd = apply(df[,c("end", "endProm")], 1, max)
  }
  else{
    df$minStart = apply(df[,c("start.x", "start.y")], 1, min)
    df$maxEnd = apply(df[,c("end.x", "end.y")], 1, max)
  }

    start.cluster = aggregate(minStart~membership, df, min)
    end.cluster = aggregate(maxEnd~membership, df, max)


    df.start.end = unique(merge(merge(start.cluster, end.cluster, by="membership"), df[,c("chr.x", "membership")], by="membership"), all.y=T)

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

individual.correlation =function(RegulGRangesOrPromoterGRanges,annotation=c("ABC", "Rao", "DNAse"),method=c("pearson", "spearman"), check.3D.circularity = T){
  #The function computes the correlation for all pairs of metacolumns of GRanges
  #The GRanges has to have only unique elements as rows
  
  metacolumns = data.frame(mcols(RegulGRangesOrPromoterGRanges))
  
  if(annotation=="ABC"){
    colnames(metacolumns) = c("ABC.Score","DI", "ScoreFire","Insulation", "CTCF", "H3K27ac", "DNAse", "H3K4me1", "H3K4me3", "p300", "H3K27me3")
    metacolumns.final = metacolumns[,-1]
    PerformanceAnalytics::chart.Correlation(metacolumns.final, method = method, histogram = T)
    
    if(check.3D.circularity){
      PerformanceAnalytics::chart.Correlation(metacolumns[,c("ABC.Score", "DI", "ScoreFire","Insulation")], method=method)
    }
  }
  else{
    metacolumns.final = metacolumns[,-c(1,2)]
    colnames(metacolumns.final) = colnames(metacolumns) = c("DI", "ScoreFire","Insulation", "CTCF", "H3K27ac", "DNAse", "H3K4me1", "H3K4me3", "p300", "H3K27me3")
    PerformanceAnalytics::chart.Correlation(metacolumns.final, method = method, histogram = T)
  }
  
}
make.comparable.set = function(GRange,method=c("ABC", "Rao", "DNAse"), element = c("promoter", "regulatory"), DNAse = NULL, contact=NULL){
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  genes.hg19 <- genes(txdb)
  symbol <- select(org.Hs.eg.db,keys = genes.hg19$gene_id, columns = c("SYMBOL"), keytype = "ENTREZID")
  genes.hg19$geneSymbol <- symbol$SYMBOL
  
  
  
  
  if(method=="ABC"){
    if(element=="promoter"){
      candidates.genes.ABC = genes.hg19[!genes.hg19$geneSymbol%in%names(GRange)]
      return(candidates.genes.ABC)
    }
    else{
      DNAseActivity.NEU.WX = DNAse[seqnames(DNAse)%in%paste0("chr",1:22)]
      overlaps.peaks.ABC = findOverlaps(GRange, DNAseActivity.NEU.WX)
      candidates.regions.ABC = DNAseActivity.NEU.WX[-subjectHits(overlaps.peaks.ABC)]
      return(candidates.regions.ABC)
    }
    
  }
  
  else if(method=="Rao"){
    if(element=="promoter"){
      candidates.genes.Rao = genes.hg19[!genes.hg19$geneSymbol%in%GRange$geneSymbol]
      return(candidates.genes.Rao)
    }
    
    else{
      overlaps.regions.Rao = findOverlaps(second(contact), GRange)
      canditates.regions.Rao = second(contact)[-queryHits(overlaps.regions.Rao)]
      return(canditates.regions.Rao)
      
    }
  }
  
  else{
    if(element=="promoter"){
      candidates.genes.DNAse = genes.hg19[!genes.hg19$geneSymbol%in%GRange$geneSymbol]
      return(candidates.genes.DNAse)
    }
    else{
      DNAseActivity.NEU.WX = DNAse[seqnames(DNAse)%in%paste0("chr",1:22)]
      overlaps.regions.DNAse = findOverlaps(GRange,DNAseActivity.NEU.WX)
      canditates.regions.DNAse = DNAseActivity.NEU.WX[-subjectHits(overlaps.regions.DNAse)]
      return(canditates.regions.DNAse)
    }
  }
  
  
}
enrichments.analysis = function(GRangeTOTest, GRangesOFCandidatesRegions, GRangeOFSigniRegions, GRangeOFNSigniRegions){
  GRangeTOTest.tmp = GRangeTOTest
  GRangesOFCandidatesRegions.tmp = GRangesOFCandidatesRegions
 
  
  GRangeTOTest.tmp$signi = countOverlaps(GRangeTOTest.tmp,GRangeOFSigniRegions)
  GRangeTOTest.tmp$nsigni = countOverlaps(GRangeTOTest.tmp,GRangeOFNSigniRegions)
  
  GRangesOFCandidatesRegions.tmp$signi = countOverlaps(GRangesOFCandidatesRegions, GRangeOFSigniRegions)
  GRangesOFCandidatesRegions.tmp$nsigni = countOverlaps(GRangesOFCandidatesRegions, GRangeOFNSigniRegions)
  
  matrix.enrichment = matrix(c(sum(GRangeTOTest.tmp$signi), sum(GRangeTOTest.tmp$nsigni), sum(GRangesOFCandidatesRegions.tmp$signi), sum(GRangesOFCandidatesRegions.tmp$nsigni)), ncol=2, nrow=2)
  return(fisher.test(matrix.enrichment))
}

plot.OR = function(fisher.output.ABC, fisher.output.Rao, fisher.output.DNAse){
  OR.ABC = data.frame(fisher.output.ABC$estimate[[1]], fisher.output.ABC$conf.int[1],fisher.output.ABC$conf.int[2])
  OR.Rao = data.frame(fisher.output.Rao$estimate[[1]],fisher.output.Rao$conf.int[1], fisher.output.Rao$conf.int[2])
  OR.DNAse = data.frame(fisher.output.DNAse$estimate[[1]], fisher.output.DNAse$conf.int[1],fisher.output.DNAse$conf.int[2])
  
  colnames(OR.ABC) = c("OR", "CI_lower", "CI_upper")
  colnames(OR.Rao) = colnames(OR.ABC)
  colnames(OR.DNAse) = colnames(OR.ABC)
  OR.ABC$Method = "ABC"
  OR.Rao$Method = "Rao"
  OR.DNAse$Method = "DNAse"
  
  OR.ALL = rbind(OR.ABC,OR.Rao,OR.DNAse)
  ggplot(OR.ALL, aes(x=OR, y=Method)) + geom_pointrange(aes(xmin=CI_lower, xmax=CI_upper, color=Method)) +  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed")+
    ylab("Annotation Method")
  
}
cor.by.CRNSize = function(df,col1, col2){
  q = aggregate(name~membership, df, function(x) length(unique(x)))
  q.subset = q[q$name>1,]
  
  r=by(df[df$membership%in%q.subset$membership,], df[df$membership%in%q.subset$membership,"membership"], FUN = function(X) cor(X[,col1], X[,col2], method = "pearson"))
  
  r.df = data.frame(membership = dimnames(r)[[1]], corr = as.vector(r))
  s = na.omit(merge(r.df,q, by="membership"))
  s = data.frame(s)
  s$t = (s$corr*sqrt(s$name-2))/sqrt(1-s$corr^2)
  s$p.value = 2*pt(abs(s$t), s$name-2, lower.tail = F)
  return(s)
}
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

add.Complexity.RRCs = function(df, method=c("ABC", "Rao", "DNAse"), decompose.graph){
  if(method=="ABC"){
    
   
    unique.EssentialGenes.CRN = unique(df[,c("membership", "TargetGene", "essentialGenes", "Expression")])
    
    asso.complexity = aggregate(essentialGenes~membership, unique.EssentialGenes.CRN, sum)
    asso.complexity$complexity = sapply(seq_along(decompose.graph), function(x) sum(degree(decompose.graph[[x]])))
    
    length.CRN = aggregate(TargetGene~membership, unique.EssentialGenes.CRN, function(x) length(unique(x)))
    length.regul.CRN = aggregate(name~membership, df, function(x) length(unique(x)))
    asso.complexity = merge(asso.complexity, length.CRN)
    asso.complexity$prop.EssentialGenes = asso.complexity$essentialGenes/asso.complexity$TargetGene
    return(asso.complexity)
  }
  
  else{
  
    
    unique.EssentialGenes.CRN = unique(df[,c("membership", "geneSymbol.x", "essentialGenes", "Expression")])
    
    asso.complexity = aggregate(essentialGenes~membership, unique.EssentialGenes.CRN, sum)
    asso.complexity$complexity = sapply(seq_along(decompose.graph), function(x) sum(degree(decompose.graph[[x]])))
    
    length.CRN = aggregate(geneSymbol.x~membership, unique.EssentialGenes.CRN, function(x) length(unique(x)))
    length.regul.CRN = aggregate(name~membership, df, function(x) length(unique(x)))
    asso.complexity = merge(asso.complexity, length.CRN)
    asso.complexity$prop.EssentialGenes = asso.complexity$essentialGenes/asso.complexity$geneSymbol.x
    
    return(asso.complexity)
  }
}

Signal.By.RRCs = function(df, method=c("ABC", "Rao", "DNAse"), aggfun = c("mean", "90th")){
  
  format = c("mean_DNAse","mean_DNAse", "mean_H3K27ac","mean_H3K27ac", "mean_p300","mean_p300",
             "mean_H3K4me3","mean_H3K4me3", "mean_H3K4me1","mean_H3K4me1", "mean_H3K27me3","mean_H3K27me3",
             "mean_CTCF","mean_CTCF","ScoreFire","ScoreFire", "INS","INS", "DI","DI")
  
  if(method=="ABC"){
    
    cols = c("mean_DNAse.PROM","mean_DNAse.ENH", "mean_H3K27ac.PROM","mean_H3K27ac.ENH", "mean_p300.PROM","mean_p300.ENH",
             "mean_H3K4me3.PROM","mean_H3K4me3.ENH", "mean_H3K4me1.PROM","mean_H3K4me1.ENH", "mean_H3K27me3.PROM","mean_H3K27me3.ENH",
             "mean_CTCF.PROM","mean_CTCF.ENH","ScoreFire.PROM","ScoreFire.ENH", "INS.PROM","INS.ENH", "DI.PROM", "DI.ENH")
    
    pairwise.cols = split(cols,format)
    
      if(aggfun=="mean"){
        
          list.df = lapply(pairwise.cols, function(x){
          formula.PROM = as.formula(paste0(x[1],"~membership"))
          formula.ENH= as.formula(paste0(x[2],"~membership"))
          
          meanRRCsPROM = aggregate(formula.PROM, df, mean, na.rm=T)
          meanRRCsEnh = aggregate(formula.ENH, df, mean, na.rm=T)
          
          meanRRCs = merge(meanRRCsPROM,meanRRCsEnh)
        } )
        
        return(Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "membership", all.x = TRUE),
               list.df))
       
      }
      
      else{
        list.dfs = lapply(pairwise.cols, function(x){
          formula.PROM = as.formula(paste0(x[1],"~membership"))
          formula.ENH= as.formula(paste0(x[2],"~membership"))
          
          NinetyPercRRCsPROM = aggregate(formula.PROM, df, function(x) quantile(x, probs=0.90, na.rm=T))
          NinetyPercRRCsENH = aggregate(formula.ENH, df, function(x) quantile(x, probs=0.90, na.rm=T))
          
          NinetyPercRRCs = merge(NinetyPercRRCsPROM,NinetyPercRRCsENH)
          return(NinetyPercRRCs)
        } )
        return(Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "membership", all.x = TRUE),
                      list.dfs)) 
      }
    }
    
   else{
     cols = c("mean_DNAse.x","mean_DNAse.y", "mean_H3K27ac.x","mean_H3K27ac.y", "mean_p300.x","mean_p300.y",
              "mean_H3K4me3.x","mean_H3K4me3.y", "mean_H3K4me1.x","mean_H3K4me1.y", "mean_H3K27me3.x","mean_H3K27me3.y",
              "mean_CTCF.x","mean_CTCF.y","mean_ScoreFire.x","mean_ScoreFire.y", "mean_INS.x","mean_INS.y", "mean_DI.x", "mean_DI.y")
     
     pairwise.cols = split(cols,format)
     
     if(aggfun=="mean"){
       
       list.df = lapply(pairwise.cols, function(x){
         formula.PROM = as.formula(paste0(x[1],"~membership"))
         formula.ENH= as.formula(paste0(x[2],"~membership"))
         
         meanRRCsPROM = aggregate(formula.PROM, df, mean, na.rm=T)
         meanRRCsEnh = aggregate(formula.ENH, df, mean, na.rm=T)
         
         meanRRCs = merge(meanRRCsPROM,meanRRCsEnh)
         return(meanRRCs)
       } )
       
       return(Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "membership", all.x = TRUE),
                     list.df))
       
     }
     
     else{
       list.dfs = lapply(pairwise.cols, function(x){
         formula.PROM = as.formula(paste0(x[1],"~membership"))
         formula.ENH= as.formula(paste0(x[2],"~membership"))
         
         NinetyPercRRCsPROM = aggregate(formula.PROM, df, function(x) quantile(x, probs=0.90, na.rm=T))
         NinetyPercRRCsENH = aggregate(formula.ENH, df, function(x) quantile(x, probs=0.90, na.rm=T))
         
         NinetyPercRRCs = merge(NinetyPercRRCsPROM,NinetyPercRRCsENH)
         return(NinetyPercRRCs)
       } )
       return(Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "membership", all.x = TRUE),
                     list.dfs)) 
     }
   }
    
   
}

add.SNPs.to.df = function(GRange1,GRange2, GRange.SNPs, df, method=c("ABC", "Rao", "DNAse"),threshold=0.05){
  
  
  if(method=="ABC"){
    
    NSigniSNPs.prom = data.frame(countOverlaps(GRange1, GRange.SNPs[GRange.SNPs$pval<=threshold]))
    colnames(NSigniSNPs.prom) = "NSigni.PROM"
    NSigniSNPs.prom$TargetGene = rownames(NSigniSNPs.prom)
    
    NnSigniSNPs.prom = data.frame(countOverlaps(GRange1, GRange.SNPs[GRange.SNPs$pval>threshold]))
    colnames(NnSigniSNPs.prom) = "NnSigni.PROM"
    NnSigniSNPs.prom$TargetGene = rownames(NnSigniSNPs.prom)
    
    NSigniSNPs.regul = data.frame(countOverlaps(GRange2, GRange.SNPs[GRange.SNPs$pval<=threshold]))
    colnames(NSigniSNPs.regul) = "NSigni.ENH"
    NSigniSNPs.regul$name = rownames(NSigniSNPs.regul)
    
    NnSigniSNPs.regul = data.frame(countOverlaps(GRange2, GRange.SNPs[GRange.SNPs$pval>threshold]))
    colnames(NnSigniSNPs.regul) = "NnSigni.ENH"
    NnSigniSNPs.regul$name = rownames(NnSigniSNPs.regul)
    
    df = merge(merge(df, NSigniSNPs.prom, by="TargetGene", all.x=T),NnSigniSNPs.prom,by="TargetGene", all.x=T )
    df = merge(merge(df,NSigniSNPs.regul,by="name", all.x=T), NnSigniSNPs.regul, by="name", all.x=T)
    
  }
  
  else{
    NSigniSNPs.prom = data.frame(countOverlaps(GRange1, GRange.SNPs[GRange.SNPs$pval<=0.05]))
    colnames(NSigniSNPs.prom) = "NSigni.PROM"
    NSigniSNPs.prom$geneSymbol.x = GRange1$geneSymbol
    
    NnSigniSNPs.prom = data.frame(countOverlaps(GRange1, GRange.SNPs[GRange.SNPs$pval>0.05]))
    colnames(NnSigniSNPs.prom) = "NnSigni.PROM"
    NnSigniSNPs.prom$geneSymbol.x = GRange1$geneSymbol
    
    name = paste0(seqnames(GRange2),":",start(GRange2),":", end(GRange2))
    
    NSigniSNPs.regul = data.frame(countOverlaps(GRange2, GRange.SNPs[GRange.SNPs$pval<=0.05]))
    colnames(NSigniSNPs.regul) = "NSigni.ENH"
    NSigniSNPs.regul$name =name
    
    NnSigniSNPs.regul = data.frame(countOverlaps(GRange2, GRange.SNPs[GRange.SNPs$pval>0.05]))
    colnames(NnSigniSNPs.regul) = "NnSigni.ENH"
    NnSigniSNPs.regul$name = name
    
    df = merge(merge(df, NSigniSNPs.prom, by="geneSymbol.x", all.x=T),NnSigniSNPs.prom,by="geneSymbol.x", all.x=T )
    df = merge(merge(df,NSigniSNPs.regul,by="name", all.x=T), NnSigniSNPs.regul, by="name", all.x=T)
    
  }
  
  return(df)
}

imput.individual.elements = function(GRangeProm2Imput,GRangeRegul2Imput, original.df, method=c("ABC", "Rao", "DNAse")){
  prom.tmp = GRangeProm2Imput
  regul.tmp = GRangeRegul2Imput
  
 
  
  if(method=="ABC"){
    
    impl=mice(data.frame(mcols(prom.tmp))[,-1], m = 5)
    prom.tmp = complete(impl)
    
    impl.Regul=mice(data.frame(mcols(regul.tmp))[,-1], m = 5)
    regul.tmp = complete(impl.Regul)
    
    colnames(prom.tmp) = c("DI.PROM", "ScoreFire.PROM", "INS.PROM","mean_CTCF.PROM", "mean_H3K27ac.PROM","mean_DNAse.PROM", "mean_H3K4me1.PROM", "mean_H3K4me3.PROM","mean_p300.PROM", "mean_H3K27me3.PROM")
    prom.tmp$TargetGene = rownames(prom.tmp)
    
    
    colnames(regul.tmp) = c("DI.ENH", "ScoreFire.ENH", "INS.ENH","mean_CTCF.ENH", "mean_H3K27ac.ENH","mean_DNAse.ENH", "mean_H3K4me1.ENH", "mean_H3K4me3.ENH","mean_p300.ENH", "mean_H3K27me3.ENH")
    regul.tmp$name = rownames(regul.tmp)
    
    df.imput = original.df
    
    df.imput[,c("DI.PROM",  "ScoreFire.PROM","INS.PROM", "mean_CTCF.PROM","mean_H3K27ac.PROM" ,
                "mean_DNAse.PROM",    "mean_H3K4me1.PROM",  "mean_H3K4me3.PROM",  "mean_p300.PROM","mean_H3K27me3.PROM", "ABC.Score.ENH",
                "DI.ENH",   "ScoreFire.ENH", "INS.ENH",  "mean_CTCF.ENH", "mean_H3K27ac.ENH","mean_DNAse.ENH",   
                "mean_H3K4me1.ENH",   "mean_H3K4me3.ENH",   "mean_p300.ENH", "mean_H3K27me3.ENH")] = NULL
    
    df.imput = merge(merge(df.imput, prom.tmp, by="TargetGene", all.x=T), regul.tmp,by="name", all.x=T)
  }
  
  else{
    impl=mice(data.frame(mcols(prom.tmp))[,-c(1,2)], m = 5)
    prom.tmp = complete(impl)
    
    impl.Regul=mice(data.frame(mcols(regul.tmp))[,-c(1,2)], m = 5)
    regul.tmp = complete(impl.Regul)
    
    colnames(prom.tmp) = c("mean_DI.x", "mean_ScoreFire.x", "mean_INS.x","mean_CTCF.x", "mean_H3K27ac.x","mean_DNAse.x", "mean_H3K4me1.x", "mean_H3K4me3.x","mean_p300.x", "mean_H3K27me3.x")
    prom.tmp$geneSymbol.x = GRangeProm2Imput$geneSymbol
    
    
    colnames(regul.tmp) = c("mean_DI.y", "mean_ScoreFire.y", "mean_INS.y","mean_CTCF.y", "mean_H3K27ac.y","mean_DNAse.y", "mean_H3K4me1.y", "mean_H3K4me3.y","mean_p300.y", "mean_H3K27me3.y")
    regul.tmp$name = paste0(seqnames(GRangeRegul2Imput),":",start(GRangeRegul2Imput),":", end(GRangeRegul2Imput))
    
    df.imput = original.df
    
    df.imput[,c("mean_DI.x", "mean_ScoreFire.x", "mean_INS.x","mean_CTCF.x", "mean_H3K27ac.x","mean_DNAse.x", 
                "mean_H3K4me1.x", "mean_H3K4me3.x","mean_p300.x", "mean_H3K27me3.x","mean_DI.y", "mean_ScoreFire.y", "mean_INS.y","mean_CTCF.y", "mean_H3K27ac.y",
                "mean_DNAse.y", "mean_H3K4me1.y", "mean_H3K4me3.y","mean_p300.y", "mean_H3K27me3.y")] = NULL
    
    df.imput = merge(merge(df.imput, prom.tmp, by="geneSymbol.x", all.x=T), regul.tmp,by="name", all.x=T)
    
  }
  
  
  return(list("prom"=prom.tmp, "regul"=regul.tmp, "all"=df.imput))
  
}
CS.by.CRN = function(GRange.Regul,GRange.CS,df, asso.complex.df, method=c("ABC", "DNAse","Rao")){
  
  tmp.Regul = GRange.Regul
  tmp.CS = GRange.CS
  tmp.asso = asso.complex.df
  
  overlaps.CS = findOverlaps(tmp.Regul, tmp.CS)
  subset.CS = tmp.CS[queryHits(overlaps.CS)]
  
  if(method=="ABC"){
    mcols(tmp.CS)[subjectHits(overlaps.CS),"Regul.name"] = names(tmp.Regul[queryHits(overlaps.CS)])
  }
  else{
    names(tmp.Regul) = paste0(seqnames(tmp.Regul),":",start(tmp.Regul),":",end(tmp.Regul))
    mcols(tmp.CS)[subjectHits(overlaps.CS),"Regul.name"] = names(tmp.Regul[queryHits(overlaps.CS)])
  }
  print(tmp.CS)
  states = data.frame(tmp.CS[!is.na(tmp.CS$Regul.name)])
  crossed.table = table(states$Regul.name, states$name)
  
  df.crossed.table = data.frame(matrix(crossed.table,ncol=18))
  colnames(df.crossed.table) = colnames(crossed.table)
  df.crossed.table$name = rownames(crossed.table)
  
  
  df.add.CS = merge(df, df.crossed.table, by="name", all.x=T)
  
  l = lapply(colnames(crossed.table), function(x) {
    f = as.formula(paste0("`",x,"`","~", "membership"))
    aggregate(f, df.add.CS, sum)
  }
  )
  df.l = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "membership", all.x = TRUE),l)
  
  tmp.asso = merge(tmp.asso,df.l, by="membership", all.x=T)
  
  prop.CS = tmp.asso[,-c(1,2,3,4,5)]/rowSums(tmp.asso[,-c(1,2,3,4,5)])
  
  #max.prop.CS = apply(prop.CS, 1, max)
  all.prop.CS = apply(prop.CS, 1, function(x) cumsum(sort(x, decreasing=T)))
  
  eighty.perc.by.RRC = lapply(1:length(all.prop.CS), function(x){
    if(any(all.prop.CS[[x]]<=0.8)){
      all.prop.CS[[x]][all.prop.CS[[x]]<=0.8]
    }
    else{
      all.prop.CS[[x]][1]
    }}
    
    )
  
  
  return(eighty.perc.by.RRC)
}

regression.propSNPs = function(candidates.regul, candidates.prom, df, model=c("beta", "gam")){
  candidates.regul$SNPs.signi = countOverlaps(candidates.regul, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05])
  candidates.regul$SNPs.nsigni = countOverlaps(candidates.regul, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.05])
  candidates.regul$prop.SigniSNPs = candidates.regul$SNPs.signi/(candidates.regul$SNPs.signi+candidates.regul$SNPs.nsigni)
  
  candidates.prom$SNPs.signi = countOverlaps(candidates.prom, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval<=0.05])
  candidates.prom$SNPs.nsigni = countOverlaps(candidates.prom, GRanges.snps.SCZ3[GRanges.snps.SCZ3$pval>0.05])
  candidates.prom$prop.SigniSNPs = candidates.prom$SNPs.signi/(candidates.prom$SNPs.signi+candidates.prom$SNPs.nsigni)
  
  candidates.regul$complexity = 1
  candidates.regul$inCRNs = FALSE
  regul.candidates.ABC.WX = candidates.regul[seqnames(candidates.regul)!="chrX"]
  
  candidates.prom$complexity = 1
  candidates.prom$inCRNs = FALSE
  promoters.candidates.ABC.WX = candidates.prom[seqnames(candidates.prom)!="chrX"]
  
  
  df.promoters.candidates.ABC.WX = data.frame(promoters.candidates.ABC.WX)
  df.promoters.candidates.ABC.WX$prop.SigniSNPs[is.na(df.promoters.candidates.ABC.WX$prop.SigniSNPs)] = 0
  df.regul.candidates.ABC.WX = data.frame(regul.candidates.ABC.WX)
  df.regul.candidates.ABC.WX$prop.SigniSNPs[is.na(df.regul.candidates.ABC.WX$prop.SigniSNPs)] = 0
  
  df.regul.candidates.ABC.WX = df.regul.candidates.ABC.WX[,c("complexity", "prop.SigniSNPs", "inCRNs")]
  df.promoters.candidates.ABC.WX = df.promoters.candidates.ABC.WX[,c("complexity", "prop.SigniSNPs", "inCRNs")]
  
  subset.df = df[,c("complexity", "prop.SigniSNPs")]
  subset.df$inCRNs = TRUE
  
  df.combined = rbind(subset.df,df.promoters.candidates.ABC.WX,df.regul.candidates.ABC.WX)
  
  if(model=="beta"){
    return(betareg(prop.SigniSNPs~complexity+inCRNs, data=df.combined[df.combined$prop.SigniSNPs>0&df.combined$prop.SigniSNPs<1,]))
  }
  else{
    return(gam(log(prop.SigniSNPs+1)~s(complexity)+inCRNs, data=df.combined, method="REML"))
  }
}

compute.cluster.enrichment = function(df,membership,npermut, method=c("ABC", "Rao", "DNAse")){
  if(method=="ABC"){
    l.genes = length(unique(df[df$membership==membership,"TargetGene"]))
    l.name = length(unique(df[df$membership==membership,"name"]))
    
    random.genes.ABC = t(replicate(npermut,sample(names(unique.Promoters.ABC), l.genes, replace = FALSE )))
    random.regul.ABC = t(replicate(npermut,sample(names(unique.Regul.ABC), l.name, replace = FALSE )))
  
    
    return(list("genes"=random.genes.ABC,"regul"=random.regul.ABC))
  }
  
  else if(method=="Rao"){
    l.genes = length(unique(df[df$membership==membership,"geneSymbol.x"]))
    l.name = length(unique(df[df$membership==membership,"name"]))
    
    random.genes.Rao = t(replicate(npermut,sample(unique.Promoters.Rao$geneSymbol, l.genes, replace = FALSE )))
    unique.Regul.Rao$name = paste0(seqnames(unique.Regul.Rao),":", start(unique.Regul.Rao),":", end(unique.Regul.Rao))
    random.regul.Rao = t(replicate(npermut,sample(unique.Regul.Rao$name, l.name, replace = FALSE )))
  
    
    return(list("genes"=random.genes.Rao,"regul"=random.regul.Rao))
    
  }
  
  else{
    l.genes = length(unique(df[df$membership==membership,"geneSymbol.x"]))
    l.name = length(unique(df[df$membership==membership,"name"]))
    
    random.genes.DNAse = t(replicate(npermut,sample(unique.Promoters.DNAse$geneSymbol, l.genes, replace = FALSE )))
    unique.Regul.DNAse$name = paste0(seqnames(unique.Regul.DNAse),":", start(unique.Regul.DNAse),":", end(unique.Regul.DNAse))
    random.regul.DNAse = t(replicate(npermut,sample(unique.Regul.DNAse$name, l.name, replace = FALSE )))
    
    return(list("genes"=random.genes.DNAse,"regul"=random.regul.DNAse))
  }
  
  
  
}
enrichment.functional.analysis = function(df, membership, list.perm, kind=c("prom", "regul"), method=c("ABC", "Rao", "DNAse")){
  
  if(method=="ABC"){
    if(kind=="prom"){
      init.values.prom = df[df$membership==membership,c("DI.PROM","ScoreFire.PROM","INS.PROM", "mean_CTCF.PROM","mean_H3K27ac.PROM","mean_DNAse.PROM",
                                                        "mean_H3K4me1.PROM","mean_H3K4me3.PROM","mean_p300.PROM", "mean_H3K27me3.PROM"
      )]
      return(lapply(1:ncol(init.values.prom), function(x){
        sum(sapply(list.perm, function(y) quantile(y[,x],prob=0.9))>=init.values.prom[,x])/10000
      }))
    }
    
    else{
      init.values.enh = df[df$membership==membership,c("DI.ENH","ScoreFire.ENH","INS.ENH", "mean_CTCF.ENH","mean_H3K27ac.ENH","mean_DNAse.ENH",
                                                       "mean_H3K4me1.ENH","mean_H3K4me3.ENH","mean_p300.ENH", "mean_H3K27me3.ENH")]
      return(lapply(1:ncol(init.values.enh), function(x){
        sum(sapply(list.perm, function(y) quantile(y[,x],prob=0.9))>=init.values.enh[,x])/10000
      }))
      
    }
  }
  
  else{
    if(kind=="prom"){
      init.values.prom = df[df$membership==membership,c("mean_DI.x", "mean_ScoreFire.x", "mean_INS.x", "mean_CTCF.x", "mean_H3K27ac.x", "mean_DNAse.x", "mean_H3K4me1.x", "mean_H3K4me3.x", "mean_p300.x",
                                                        "mean_H3K27me3.x"
      )]
      return(lapply(1:ncol(init.values.prom), function(x){
        sum(sapply(list.perm, function(y) quantile(y[,x],prob=0.9))>=init.values.prom[,x])/10000
      }))
    }
    
    else{
      init.values.enh = df[df$membership==membership,c("mean_DI.y", "mean_ScoreFire.y", "mean_INS.y", "mean_CTCF.y", "mean_H3K27ac.y", "mean_DNAse.y", "mean_H3K4me1.y", "mean_H3K4me3.y", "mean_p300.y",
                                                       "mean_H3K27me3.y")]
      return(lapply(1:ncol(init.values.enh), function(x){
        sum(sapply(list.perm, function(y) quantile(y[,x],prob=0.9))>=init.values.enh[,x])/10000
      }))
      
    }
  }
  
  
  
  
}
CI.regression = function(model.object, method=c("ABC", "Rao", "DNAse"), alpha=0.05){
  coefficients = coef(model.object)
  summary.model = summary(model.object)
  df = summary.model$df[2]
  
  firstpart = qt(c(alpha/2,1-alpha/2),df[2])
  se = sqrt(diag(vcovHC(model.object, type="HC0")))
  
  robust.CI = coef(model.object) + se %o% tt
  
  coefficients_name = rownames(robust.CI)
  CI.lower = robust.CI[,1]
  CI.upper = robust.CI[,2]
  
  df.coeff = data.frame(cbind(coefficients_name, coefficients,CI.lower,CI.upper ))
  colnames(df.coeff) = c("name", "Estimate", "CI.l", "CI.u")
  
  df.coeff$Estimate = as.numeric(as.character(df.coeff$Estimate))
  df.coeff$CI.l = as.numeric(as.character(df.coeff$CI.l))
  df.coeff$CI.u = as.numeric(as.character(df.coeff$CI.u))
  df.coeff$Method=method
  return(df.coeff)
}

export.GRanges.for.Genesannot = function(GRanges, type=c("elements", "candidates"), method=c("ABC", "Rao", "DNAse")){
  
  if(type=="elements"){
    if(method=="ABC"){
      GeneSet = data.frame(matrix(names(GRanges), ncol=1))
      GeneCoord = cbind(GeneSet, data.frame(GRanges)[,1:3])
      colnames(GeneCoord) = c("GENE", "CHR", "START", "END")
      return(list("GSet" = na.omit(GeneSet), "GCoor" = na.omit(GeneCoord)))
    }
    else{
      names(GRanges) = GRanges$geneSymbol
      GeneSet = data.frame(matrix(names(GRanges), ncol=1))
      GeneCoord = cbind(GeneSet, data.frame(GRanges)[,1:3])
      colnames(GeneCoord) = c("GENE", "CHR", "START", "END")
      return(list("GSet" = na.omit(GeneSet), "GCoor" = na.omit(GeneCoord)))
    }
  }
  else{
    names(GRanges) = GRanges$geneSymbol
    GeneSet = data.frame(matrix(names(GRanges), ncol=1))
    GeneCoord = cbind(GeneSet, data.frame(GRanges)[,1:3])
    colnames(GeneCoord) = c("GENE", "CHR", "START", "END")
    return(list("GSet" = na.omit(GeneSet), "GCoor" = na.omit(GeneCoord)))
  }
  
}
export.GRanges.for.Regulannot = function(GRanges){
  df = data.frame(GRanges)
  df[,1:3]
}
#' @description Function which creates GRanges for Compartments analysis by chromosome
#' by default the 3 first PC are included
#' @param PC.by.locus Name of the file containing the chromosome, bin and values of the 3 first principal components
#' @param resolution Width of a bin in nucleotides
#' @param chrl list of the chromosome lengths in nucleotides 
.make.GRanges.compartments = function(PC.by.locus, resolution=1000000,chrl){
  chroms = paste0("chr",1:22)
  PC.compartments = read.table(PC.by.locus, header=F)
  colnames(PC.compartments) = c("chr", "bin", "PC1", "PC2", "PC3")
  
  PC.compartments.WX = PC.compartments[PC.compartments$chr!="chrX", ]
  
  list.PC.compartments = list()
  
  for(chr in chroms){
    PC.compartments.tmp = PC.compartments.WX[PC.compartments.WX$chr ==chr,]
    
    #    windows = seq(0, resolution*(nrow(PC.compartments.tmp)), resolution)
    windows = c(0,PC.compartments.tmp$bin*resolution)
    # On s'assure que le dernier intervalle ne dépasse pas la fin du chromosome
    if (!missing(chrl)) windows[length(windows)] = min(windows[length(windows)],chrl[[chr]])
    start.stop = data.frame(cbind(windows[-length(windows)]+1, windows[-1]))
    colnames(start.stop) = c("start", "stop")
    
    list.PC.compartments[[chr]] = cbind(PC.compartments.tmp, start.stop)
  }
  
  PC.compartments.analyse = do.call(rbind, list.PC.compartments)
  PC.compartments.analyse[,3:5] = PC.compartments.analyse[,3:5]*-1
  
  GRanges.PC= GRanges(seqnames = PC.compartments.analyse$chr, ranges=IRanges(start = PC.compartments.analyse$start,end =PC.compartments.analyse$stop,names=paste0("bin", 1:nrow(PC.compartments.analyse))),
                      PC1=PC.compartments.analyse$PC1, PC2=PC.compartments.analyse$PC2, PC3=PC.compartments.analyse$PC3)
  
  return(GRanges.PC)
}

#' @description Function which creates GRanges for Compartments analysis by chromosome arm
#' by default the 3 first PC are included
#' @param PC.by.locus Name of the file containing the chromosome, bin and values of the 3 first principal components
#' @param resolution Width of a bin in nucleotides  
.make.GRanges.compartments.arms = function(PC.by.locus, resolution=1000000){
  chroms = paste0("chr",1:22)
  PC.compartments = read.table(PC.by.locus, header=F)
  colnames(PC.compartments) = c("chr", "bin", "PC1", "PC2", "PC3")
  
  list.PC.compartments = list()
  
  for(chr in 1:22){
    PC.compartments.tmp = PC.compartments[PC.compartments$chr ==chroms[chr],]
    windows = c((PC.compartments.tmp$bin[1]-1),PC.compartments.tmp$bin)*resolution
    
    # Bras p (s'il y en a un)
    if (PC.compartments.tmp$bin[1]*resolution < hg19$centromerStart[chr])
    {
      bras = paste0(chroms[chr],"p")    
      # On s'assure que le dernier intervalle ne dépasse pas le début du centromère
      windows.p = c(windows[windows<hg19$centromerStart[chr]],hg19$centromerStart[chr])
      start.stop = data.frame(cbind(windows.p[-length(windows.p)]+1, windows.p[-1]))
      colnames(start.stop) = c("start", "stop")
      # On s'assure que le premier intervalle du bras q démarre après la fin du centromère et que le dernier intervalle ne dépasse pas la fin du chromosome
      windows.q = pmin(c(max(hg19$centromerEnd[chr],min(windows[windows>hg19$centromerEnd[chr]])-resolution)+1,windows[windows>hg19$centromerEnd[chr]]),hg19$length[chr])
      
      list.PC.compartments[[bras]] = cbind(bras,PC.compartments.tmp[1:nrow(start.stop),], start.stop)
    }
    # On s'assure que le premier intervalle démarre après la fin du centromère et que le dernier intervalle ne dépasse pas la fin du chromosome
    else windows.q = pmin(windows[windows>hg19$centromerEnd[chr]],hg19$length[chr])
    # Bras q
    bras = paste0(chroms[chr],"q")    
    start.stop = data.frame(cbind(windows.q[-length(windows.q)]+1, windows.q[-1]))
    colnames(start.stop) = c("start", "stop")
    list.PC.compartments[[bras]] = cbind(bras,PC.compartments.tmp[(PC.compartments.tmp$bin-1)*resolution>=hg19$centromerEnd[chr],], start.stop)
  }
  
  PC.compartments.analyse = do.call(rbind, list.PC.compartments)
  PC.compartments.analyse[,3:5] = PC.compartments.analyse[,3:5]*-1
  
  GRanges.PC= GRanges(seqnames = PC.compartments.analyse$chr, ranges=IRanges(start = PC.compartments.analyse$start,end =PC.compartments.analyse$stop,names=paste0("bin", 1:nrow(PC.compartments.analyse))),bras=PC.compartments.analyse$bras,
                      PC1=PC.compartments.analyse$PC1, PC2=PC.compartments.analyse$PC2, PC3=PC.compartments.analyse$PC3)
  
  return(GRanges.PC)
}

# The compute_GCcontent function is taken from the Bioconductor package gcapc
compute_GCcontent = function(region,genome){
  seqs <- getSeq(genome,region) # slow
  gcpos <- startIndex(vmatchPattern("S", seqs, fixed=F))
  npos <- startIndex(vmatchPattern("N", seqs, fixed=T))
  round((sapply(gcpos,length) - sapply(npos,length))/sapply(seqs,length),3)
}

#' @param GRanges.PC GRanges object with the 3 first principal components and the GC content for each bin of the genome
#' @param corByChr boolean: if TRUE, the Pearson correlation between GC content and each of the 3 first PCs is computed for each chromosome
define.active.compartments.GC = function(GRanges.PC, corByChr = TRUE){
  
  df.PC = data.frame(GRanges.PC)
  
  if(corByChr){
    cor.PC1.genes = ddply(df.PC, .(seqnames), summarise, "corr" = cor(PC1, gc, method = "pearson"))
    cor.PC2.genes = ddply(df.PC, .(seqnames), summarise, "corr" = cor(PC2, gc, method = "pearson"))
    cor.PC3.genes = ddply(df.PC, .(seqnames), summarise, "corr" = cor(PC3, gc, method = "pearson"))
    
    cor.gc.PC = merge(merge(cor.PC1.genes, cor.PC2.genes, by="seqnames"),cor.PC3.genes, by="seqnames")
    colnames(cor.gc.PC) = c("chr", "COR.PC1.GC", "COR.PC2.GC", "COR.PC3.GC")
    
    print(cor.gc.PC)
  }
  return(cor.gc.PC)
}

#' @param GRanges.PC GRanges object with the 3 first principal components and the GC content for each bin of the genome
#' @param corByChr boolean: if TRUE, the Pearson correlation between GC content and each of the 3 first PCs is returned for each chromosome arm, if FALSE a GRanges object with a column "Compartment" indicating the compartment assigned to each bin is returned.
define.active.compartments.arms.GC = function(GRanges.PC, corByChr = TRUE){
  
  df.PC = data.frame(GRanges.PC)
  
  cor.PC1.gc = ddply(df.PC, .(bras), summarise, "corr" = cor(PC1, gc, method = "pearson"))
  cor.PC2.gc = ddply(df.PC, .(bras), summarise, "corr" = cor(PC2, gc, method = "pearson"))
  cor.PC3.gc = ddply(df.PC, .(bras), summarise, "corr" = cor(PC3, gc, method = "pearson"))
  
  if(corByChr){
    cor.gc.PC = merge(merge(cor.PC1.gc, cor.PC2.gc, by="bras"),cor.PC3.gc, by="bras")
    colnames(cor.gc.PC) = c("chr", "COR.PC1.GC", "COR.PC2.GC", "COR.PC3.GC")
    
    print(cor.gc.PC)
    return(cor.gc.PC)
  }
  else
  {
    # moyenne du contenu GC dans les régions où la 1re composante principale est positive ou négative par bras de chromosome
    mneg = tapply(mcols(GRanges.PC)[mcols(GRanges.PC)[,"PC1"] < 0, "gc"],mcols(GRanges.PC)[mcols(GRanges.PC)[,"PC1"] < 0,"bras"],mean)
    mpos = tapply(mcols(GRanges.PC)[mcols(GRanges.PC)[,"PC1"] > 0, "gc"],mcols(GRanges.PC)[mcols(GRanges.PC)[,"PC1"] > 0,"bras"],mean)
    # On répète les statistiques par bras de chromosome pour toutes les bins de chaque bras
    bin.par.chr = table(mcols(GRanges.PC)[,"bras"])
    mneg.rep = rep(mneg,bin.par.chr)
    mpos.rep = rep(mpos,bin.par.chr)
    cor.PC1.gc.rep = rep(cor.PC1.gc$corr,bin.par.chr)
    cor.PC2.gc.rep = rep(cor.PC2.gc$corr,bin.par.chr)
    cor.PC3.gc.rep = rep(cor.PC3.gc$corr,bin.par.chr)
    # Si la 1re PC est la plus grande, on change le signe de la 1re composante principale s'il y a plus de gènes dans les régions où la 1re PC est négative, sinon on met à NA
    GRanges.PC$PC1bonsigne =  ifelse(abs(cor.PC1.gc.rep)>abs(cor.PC2.gc.rep) & abs(cor.PC1.gc.rep)>abs(cor.PC3.gc.rep), GRanges.PC$PC1 *ifelse(mneg.rep > mpos.rep,-1,1),NA)
    GRanges.PC$Compartment = ifelse(GRanges.PC$PC1bonsigne>0, "A", "B")
    return(GRanges.PC[,c("bras","Compartment")])
    
  }
}

#' @param PC.by.locus Name of the file containing the chromosome, bin and values of the 3 first principal components
#' @param resolution Width of a bin in nucleotides
#' @param genome Name of the human genome build
#' @param genes GRanges object with the coordinates of genes in the selected human genome build
#' @param corByChr boolean: if TRUE, the Pearson correlation between GC content and each of the 3 first PCs is returned for each chromosome arm, if FALSE a GRanges object with a column "Compartment" indicating the compartment assigned to each bin is returned.
define.active.compartments.arms = function(PC.by.locus,resolution=1000000, genome="hg19", genes, corByChr = TRUE){
  GRanges.PC = .make.GRanges.compartments.arms(PC.by.locus,resolution=resolution)
  
  GRanges.PC$ngenes = countOverlaps(GRanges.PC, genes)
  #  GRanges.PC$meth=method
  
  df.PC = data.frame(GRanges.PC)
  
  # corPC1Genes = cor(df.PC$PC1, df.PC$ngenes, method="spearman")
  # corPC2Genes = cor(df.PC$PC2, df.PC$ngenes, method="spearman")
  # corPC3Genes = cor(df.PC$PC3, df.PC$ngenes, method="spearman")
  
  # print(c(corPC1Genes, corPC2Genes,corPC3Genes))
  # indexMaxCor = which.max(c(corPC1Genes,corPC2Genes,corPC3Genes))
  
  cor.PC1.genes = ddply(df.PC, .(bras), summarise, "corr" = cor(PC1, ngenes, method = "pearson"))
  cor.PC2.genes = ddply(df.PC, .(bras), summarise, "corr" = cor(PC2, ngenes, method = "pearson"))
  cor.PC3.genes = ddply(df.PC, .(bras), summarise, "corr" = cor(PC3, ngenes, method = "pearson"))
  
  if(corByChr){
    cor.densite.genes.PC = merge(merge(cor.PC1.genes, cor.PC2.genes, by="bras"),cor.PC3.genes, by="bras")
    colnames(cor.densite.genes.PC) = c("chr", "COR.PC1.GDENSITE", "COR.PC2.GDENSITE", "COR.PC3.GDENSITE")
    
    print(cor.densite.genes.PC)
    return(cor.densite.genes.PC)
  }
  else {
    #density.mediane.genes = median(df.PC$ngenes)
    #GRanges.PC$Compartment = ifelse(GRanges.PC$ngenes>density.mediane.genes, "A", "B")
    
    # moyenne du nombre de gènes dans les régions où la 1re composante principale est positive ou négative par bras de chromosome
    mneg = tapply(mcols(GRanges.PC)[mcols(GRanges.PC)[,"PC1"] < 0, "ngenes"],mcols(GRanges.PC)[mcols(GRanges.PC)[,"PC1"] < 0,"bras"],mean)
    mpos = tapply(mcols(GRanges.PC)[mcols(GRanges.PC)[,"PC1"] > 0, "ngenes"],mcols(GRanges.PC)[mcols(GRanges.PC)[,"PC1"] > 0,"bras"],mean)
    # On répète les statistiques par bras de chromosome pour toutes les bins de chaque bras
    bin.par.chr = table(mcols(GRanges.PC)[,"bras"])
    mneg.rep = rep(mneg,bin.par.chr)
    mpos.rep = rep(mpos,bin.par.chr)
    cor.PC1.genes.rep = rep(cor.PC1.genes$corr,bin.par.chr)
    cor.PC2.genes.rep = rep(cor.PC2.genes$corr,bin.par.chr)
    cor.PC3.genes.rep = rep(cor.PC3.genes$corr,bin.par.chr)
    # Si la 1re PC est la plus grande, on change le signe de la 1re composante principale s'il y a plus de gènes dans les régions où la 1re PC est négative, sinon on met à NA
    GRanges.PC$PC1bonsigne =  ifelse(abs(cor.PC1.genes.rep)>abs(cor.PC2.genes.rep) & abs(cor.PC1.genes.rep)>abs(cor.PC3.genes.rep), GRanges.PC$PC1 *ifelse(mneg.rep > mpos.rep,-1,1),NA)
    GRanges.PC$Compartment = ifelse(GRanges.PC$PC1bonsigne>0, "A", "B")
    return(GRanges.PC[,c("bras","Compartment")])
  }
}


#' @param PC.by.locus Name of the file containing the chromosome, bin and values of the 3 first principal components
#' @param resolution Width of a bin in nucleotides
#' @param genome Name of the human genome build
#' @param genes GRanges object with the coordinates of genes in the selected human genome build
#' @param corByChr boolean: if TRUE, the Pearson correlation between gene content and each of the 3 first PCs is computed for each chromosome
define.active.compartments = function(PC.by.locus,resolution=1000000, genome="hg19", genes, corByChr = TRUE){
  genome <- getBSgenome(genome)
  chrl = list()
  for (chr in paste0("chr",1:22)) chrl[chr] = length(genome[[chr]])
  GRanges.PC = .make.GRanges.compartments(PC.by.locus,resolution=resolution,chrl=chrl)
  
  GRanges.PC$ngenes = countOverlaps(GRanges.PC, genes)
  #  GRanges.PC$meth=method
  
  df.PC = data.frame(GRanges.PC)
  
  # corPC1Genes = cor(df.PC$PC1, df.PC$ngenes, method="spearman")
  # corPC2Genes = cor(df.PC$PC2, df.PC$ngenes, method="spearman")
  # corPC3Genes = cor(df.PC$PC3, df.PC$ngenes, method="spearman")
  
  # print(c(corPC1Genes, corPC2Genes,corPC3Genes))
  # indexMaxCor = which.max(c(corPC1Genes,corPC2Genes,corPC3Genes))
  
  if(corByChr){
    cor.PC1.genes = ddply(df.PC, .(seqnames), summarise, "corr" = cor(PC1, ngenes, method = "pearson"))
    cor.PC2.genes = ddply(df.PC, .(seqnames), summarise, "corr" = cor(PC2, ngenes, method = "pearson"))
    cor.PC3.genes = ddply(df.PC, .(seqnames), summarise, "corr" = cor(PC3, ngenes, method = "pearson"))
    
    cor.densite.genes.PC = merge(merge(cor.PC1.genes, cor.PC2.genes, by="seqnames"),cor.PC3.genes, by="seqnames")
    colnames(cor.densite.genes.PC) = c("chr", "COR.PC1.GDENSITE", "COR.PC2.GDENSITE", "COR.PC3.GDENSITE")
    
    print(cor.densite.genes.PC)
  }
  
  #density.mediane.genes = median(df.PC$ngenes)
  #GRanges.PC$Compartment = ifelse(GRanges.PC$ngenes>density.mediane.genes, "A", "B")
  
  # mneg = mean(mcols(GRanges.PC)[mcols(GRanges.PC)[,indexMaxCor] < 0, "ngenes"])
  # mpos = mean(mcols(GRanges.PC)[mcols(GRanges.PC)[,indexMaxCor] > 0, "ngenes"])
  
  # if(mneg > mpos){
  # GRanges.PC[,indexMaxCor] =  GRanges.PC[,indexMaxCor] *-1
  # GRanges.PC$Compartment = ifelse(mcols(GRanges.PC)[,indexMaxCor]>0, "A", "B")
  # return(GRanges.PC[,c(indexMaxCor,5)])
  # }
  
  # else{
  # GRanges.PC$Compartment = ifelse(mcols(GRanges.PC)[,indexMaxCor]>0, "A", "B")
  # return(GRanges.PC[,c(indexMaxCor,5)])
  # }
  return(cor.densite.genes.PC)
}