process.ABC = function(ABC.inputFile) {
  #Function taking in input the EnhancerPredictions.txt output file from Score-ABC software
  #X chromsome is removed from analysis and gene promoters was added
  
  ABC.input = read.table(ABC.inputFile, header=T)
  
  #Removing of X chromsome due to specific organization of 3D architecture from Dosage Compensation
  ABC.input$chr = as.factor(ABC.input$chr)
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
  #This function creates dataframe for the 2 other annotation methods 
  #positive.contact parameter is the output from juicer software
  #DNAse.file parameter is a narrowpeak file 
  
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
    DNAse.file = import(DNAse.file, format="narrowpeak")
    
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

Pairs.from.ABC = function(input) {
  #Function takes as input the process.ABC output 
  #The function returns a Pairs of GRanges with enhancers and promoters
  GRanges.enhancers = GRanges(seqnames=input$chr.x, ranges=IRanges(start=input$start, end=input$end, names=input$name), ABC.Score=input$ABC.Score)
  Granges.promoters = GRanges(seqnames=input$chr.y, ranges=IRanges(start=input$startProm, end=input$endProm, names=input$TargetGene), ABC.Score=input$ABC.Score)
  
  return(Pairs(GRanges.enhancers, Granges.promoters))
  
}

coverage.by.Pair = function(Pair.Object) {
  #This function takes as input the object returned by Pairs.from.ABC function
  start.Pair = sapply(Pair.Object, function(x) min(start(first(x)), start(second(x))))
  end.Pair = sapply(Pair.Object, function(x) max(end(first(x)), end(second(x))))
  
  return(GRanges.Pair = GRanges(seqnames = seqnames(first(Pair.Object)), ranges=IRanges(start.Pair, end.Pair, names=paste0("Pair",1:length(Pair.Object)))))
}

create.CRHs.ABC = function(process.ABC.output) { 
  #The function takes as input the data.frame returned by process.ABC 
  #The function returns a graph from igraph package 
  nodes = process.ABC.output[, c("name", "TargetGene")]
  enhancers = unique(nodes$name)
  vertices =  data.frame(rbind(matrix(unique(nodes$name), ncol=1), matrix(unique(nodes$TargetGene), ncol=1)))
    
  colnames(vertices) = c("name")
  vertices$type = vertices[,"name"]%in%enhancers
  vertices$col = ifelse(vertices[,"name"]%in%enhancers, "blue", "red")
    
  graph = graph_from_data_frame(nodes, directed=F, vertices=vertices)
  return(graph)
}

add.membership = function(graph, df){
  #This function takes as parameters the graph returned by create.CRHs.ABC
  #and the annotated data.frame  
  
  membership = components(graph)$membership
  df.membership = data.frame(membership)
  df.membership$name = rownames(df.membership)
  
  df.membership = df.membership[df.membership$name%in%df$name,]
  df = merge(df, df.membership, by="name")
  
  return(df)
  
}

coverage.by.CRH = function(graph, df) { 
  #This function takes as input 
  
  df$minStart = apply(df[,c("start", "startProm")], 1, min)
  df$maxEnd = apply(df[,c("end", "endProm")], 1, max)
  
  start.cluster = aggregate(minStart~membership, df, min)
  end.cluster = aggregate(maxEnd~membership, df, max)
  
  
  df.start.end = unique(merge(merge(start.cluster, end.cluster, by="membership"), df[,c("chr.x", "membership")], by="membership"), all.y=T)
  
  GRanges.cluster = GRanges(seqnames = df.start.end$chr, ranges = IRanges(start=df.start.end$minStart, end=df.start.end$maxEnd, names=df.start.end$compo.graph.membership ) )
  
  return(list("clusters"=GRanges.cluster, "data.frame"=df))
}

Complexity.CRHs = function(graph,enhancers,genes ,extract.central.genes=T) {
  
  compo = components(graph)
  
  
  n.unique.genes = length(unique(genes))
  n.unique.enhancers = length(unique(enhancers))
  
  polygamous.genes = 1-(table(compo$csize)[1][[1]] / n.unique.genes)
  polygamous.enhancers = 1-(table(compo$csize)[1][[1]] / n.unique.enhancers)
  
  prom.enhancer.clusters = tapply(names(compo$membership),compo$membership,function(vec,genes) table(vec%in%genes),genes=genes)
  prom.enhancer.clusters.mat = matrix(unlist(prom.enhancer.clusters),length(unlist(prom.enhancer.clusters))/2,2,byrow = T)
  
  OneN.genes = sum(prom.enhancer.clusters.mat[prom.enhancer.clusters.mat[,1]==1&prom.enhancer.clusters.mat[,2]>1,2])/n.unique.genes
  OneN.enhancers = sum(prom.enhancer.clusters.mat[prom.enhancer.clusters.mat[,2]==1&prom.enhancer.clusters.mat[,1]>1,1])/n.unique.enhancers
  
  if(extract.central.genes){
    
    subgraphs = decompose(graph)
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
    
    return(list("polygamous" = list("genes"= polygamous.genes, "enhancers"=polygamous.enhancers), "1N" = list("genes" = OneN.genes, "enhancers" = OneN.enhancers), "Centrality" = list("list.genes" = central.genes , "n.genes"= n.central.genes )))
  }
  
  else {
    return(list("polygamous" = list("genes"= polygamous.genes, "enhancers"=polygamous.enhancers), "1N" = list("genes" = OneN.genes, "enhancers" = OneN.enhancers)))
  }
  
}

precompute_AB = function(AB_file) {
  AB_file_split <- unlist(reduce(split(taloggcc.A006, ~Compartment)))
  AB_file_split$Compartment <- names(AB_file)
  
  names(AB_file) = 1:length(AB_file)
  AB_file
}

overlaps_CRHs_to_AB = function(GRanges_CRHs, AB_Compartments){
  overlaps_CRHs_AB = findOverlaps(GRanges_CRHs, AB_Compartments)
  
  unique_Overlaps = unique(queryHits(overlaps_CRHs_AB)) 
  
  for(q in unique_Overlaps){
    subj = subjectHits(overlaps_CRHs_AB[queryHits(overlaps_CRHs_AB) == q])
    mcols(GRanges_CRHs)[q,"Compartment"] = paste(mcols(AB_Compartments)[subj,"Compartment"][c(1,length(mcols(AB_Compartments)[subj,"Compartment"]))], collapse="")
    
  }
  
  GRanges_CRHs$Compartment = ifelse(GRanges_CRHs$Compartment=="BA","AB",GRanges_CRHs$Compartment)
}