getReplicate <- function(x, uname){
  prw <- uname2prw(uname)
  rep.unames <- getUnames(x, plate=prw$plate, row=prw$row, col=prw$col)
  setdiff(rep.unames, uname)
}

getBadWells <- function(x){
  prw2uname(list(plate=x@screenLog$Plate, replicate=x@screenLog$Sample, well=x@screenLog$Well))
}

cleanDistMatrix <- function(x,distMatrix){

  ## convert the dist or matrix object into a symatric matrix
  if(class(distMatrix) == 'dist') distMatrix <- as.matrix(distMatrix)
  distMatrix <- (distMatrix+t(distMatrix))/2
  for(row in 1:nrow(distMatrix)) distMatrix[row,row] <- NA

  ## filtering out non-sample wells and bad wells listed in the screenlog file
  controlStatus <- getWellFeatures(x, rownames(distMatrix), 'controlStatus')
  z <- which(controlStatus == 'sample')
  distMatrix <- distMatrix[z,z]
  unames <- setdiff(rownames(distMatrix), getBadWells(x))
  distMatrix <- distMatrix[unames, unames]

  ## filtering for genes with unique entrez ID and all replicates
  genes <- getWellFeatures(x, unames, 'LocusID')
  count <- sapply(unique(genes), function(g) length(which(genes == g)))
  z <- which(count != length(x$replicate))
  if(length(z) > 0){
    unames <- unames[which(genes %in% names(count)[-z])]
    distMatrix <- distMatrix[unames, unames]
  }
  distMatrix
}
