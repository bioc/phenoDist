distToNeg <- function(x, distMatrix, neg='rluc'){

  ## make sure distMatrix is a symmetric matrix
  if(class(distMatrix) == 'dist') distMatrix <- as.matrix(distMatrix)
  else distMatrix <- (distMatrix+t(distMatrix))/2

  ## mark bad wells listed in screenLog with NAs
  badWells <- getBadWells(x)
  if(length(badWells)>0){
    z <- which(rownames(distMatrix) %in% badWells)
    if(length(z) > 0){
      for(i in z){
        distMatrix[i,] <- rep(NA, times=ncol(distMatrix))
        distMatrix[,i] <- rep(NA, times=nrow(distMatrix))
      }
    }
  }

  unames <- rownames(distMatrix)
  res <- sapply(unames, function(u){
  ## to compare with the negative control wells on the same plate
    prw <- uname2prw(u)
    negUnames <- getUnames(x, plate=prw$plate, replicate=prw$replicate, content=neg)
    negUnames <- setdiff(negUnames, u)
    negUnames <- intersect(unames, negUnames)
    if(length(negUnames)>0) d <- mean(distMatrix[u, negUnames], na.rm=T) else d <- NA
    d
  })
  res
}

