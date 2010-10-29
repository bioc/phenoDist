ctlSeparatn <- function(x, pheno, neg='rluc', pos='ubc', ...){
  unames <- names(pheno)
  controlStatus <- getWellFeatures(x, unames, 'controlStatus')
  negUnames <- unames[which(controlStatus==neg)]
  posUnames <- unames[which(controlStatus==pos)]

  negScores <- pheno[negUnames]
  negScores <- negScores[which(!is.na(negScores))]
  posScores <- pheno[posUnames]
  posScores <- posScores[which(!is.na(posScores))]

  zprime(negScores, posScores, ...)
}

repDistRank <- function(x, distMatrix){
  distMatrix <- cleanDistMatrix(x, distMatrix)
  unames <- rownames(distMatrix)
  sapply(unames,function(u){
    repUnames <- getReplicate(x,u)
    rl <- rank(distMatrix[u,])
    res <- median(rl[repUnames])
    res
  })
}

repCorr <- function(x, pheno, ...){
  unames <- names(pheno)
  rep1 <- pheno[which(unames %in% getUnames(x,replicate=1))]
  rep2 <- pheno[which(unames %in% getUnames(x,replicate=2))]
  cor(rep1, rep2, use='complete.obs', ...)
}
