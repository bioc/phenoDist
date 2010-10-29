clusterDist <- function(x, distMatrix, clusterFun='hclust', ...){
  distMatrix <- cleanDistMatrix(x, distMatrix)
  unames <- rownames(distMatrix)
  genes <- getWellFeatures(x, unames, 'LocusID')

  ## average replicates
  unique.genes <- unique(genes)
  distMatrix <- lapply(unique.genes, function(g){
    z <- which(genes == g)
    colMeans(distMatrix[z,])
  })
  distMatrix <- do.call(rbind, distMatrix)
  distMatrix <- lapply(unique.genes, function(g){
    z <- which(genes == g)
    rowMeans(distMatrix[,z])
  })
  distMatrix <- do.call(cbind, distMatrix)
  rownames(distMatrix) <- colnames(distMatrix) <- unique.genes
  
  tree <- do.call(clusterFun, list(d=as.dist(distMatrix), ...))
  tree
}

enrichAnalysis <- function(x, cl, terms=c('GO','KEGG'), verbose=FALSE, ...){
  match.arg(terms)
  universe <- names(cl)
  paraName <- paste(terms, 'HyperGParams', sep='')
  if(verbose) cat('current cluster  -  total clusters\n')
  res <- lapply(1:max(cl), function(c){
    if(verbose) cat(c,'-',max(cl),'\n')
    selected <- names(cl[which(cl == c)])
    params <- new(paraName, geneIds=selected, universeGeneIds=universe, ...)
    score <- try(hyperGTest(params))
    score
  })
  res
}