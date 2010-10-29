PDMByWellAvg <- function(profiles, selectedWellFtrs, transformMethod=c('none', 'scale', 'PCA'), distMethod=c('manhattan', 'euclidean', 'correlation', 'mahalanobis'), nPCA){

  ## only keep selectedFtrs
  if('uname' %in% colnames(profiles)){
    rownames(profiles) <- profiles$uname
    profiles$uname <- NULL
  }
  if(!missing(selectedWellFtrs)) profiles <- profiles[,selectedWellFtrs]
  
  ## transform well features
  transformMethod <- match.arg(transformMethod)
  if(transformMethod=='scale') profiles <- scale(profiles)
  else if(transformMethod=='PCA'){
    profiles <- prcomp(profiles, scale=TRUE)$x
    if(ncol(profiles) > nPCA) profiles <- profiles[,1:nPCA]
    else warning(sprintf('PCA transformation only returns %d dimensions', ncol(profiles)))
  }
  
  ## calculate distance matrix
  distMethod <- match.arg(distMethod)
  if(distMethod %in% c('manhattan', 'euclidean')) res <- as.matrix(dist(profiles, method=distMethod))
  else if(distMethod=='correlation') res <- 1-abs(cor(t(profiles)))
  else if(distMethod=='mahalanobis'){
    require('MASS')
    Sx <- cov(profiles)
    res <- lapply(1:nrow(profiles),function(row){
      mahalanobis(profiles, profiles[row,], cov=Sx)
    })
    res <- do.call(rbind, res)
  }

  rownames(res) <- colnames(res) <- rownames(profiles)
  res
}


PDMByFactorAnalysis <- function(x, unames, selectedCellFtrs, distMethod=c('manhattan','euclidean', 'correlation','mahalanobis'), nFactors, scores=c('regression','Bartlett'), ...){

## load features for individual cells
  cellFtrs <- collectCellFeatures(x, unames)

cellUnames <- cellFtrs$uname
## remove undesired features
  if(!missing(selectedCellFtrs)) cellFtrs <- cellFtrs[,selectedCellFtrs]
  else{
    cellFtrs$uname <- NULL
    cellFtrs$spot <- NULL
    cellFtrs$id <- NULL
  }

## factor analysis for cell features
  cellFtrsFA <- factanal(x=cellFtrs, factors=nFactors, scores=scores, ...)

## average FA scores by well
  profilesFA <- lapply(unique(cellUnames), function(uname){
    z <- which(cellUnames==uname)
    res <- colSums(cellFtrsFA$scores[z,])/nFactors
  })
  profilesFA <- do.call(rbind, profilesFA)
  rownames(profilesFA) <- unames
  
  PDMByWellAvg(profiles=profilesFA, transformMethod='none', distMethod=distMethod)
}


PDMByKS <- function(x, unames, neg='rluc', selectedCellFtrs, distMethod=c('manhattan', 'euclidean', 'correlation', 'mahalanobis')){

  if(missing(selectedCellFtrs)){
    ftrs <- collectCellFeatures(x, unames[1])
    selectedCellFtrs <- setdiff(colnames(ftrs),c('uname', 'spot', 'id', 'c.g.x', 'c.g.y', 'n.g.x', 'n.g.y'))
  }
  
## load negative control features

  negPlates <- uname2prw(unames)[,1:2]
  negPlates <- negPlates[!duplicated(negPlates),]
  negFtrsPool <- lapply(1:nrow(negPlates), function(row){
    negUnames <- getUnames(x, content=neg, plate=negPlates$plate[row], replicate=negPlates$replicate[row])
    negUnames <- setdiff(negUnames, getBadWells(x))
    negFtrs <- lapply(negUnames, function(uname){ 
      ftrs <- collectCellFeatures(x,uname)
      ftrs <- ftrs[, selectedCellFtrs]
    })
    names(negFtrs) <- negUnames
    negFtrs
  })

  profilesKS <- lapply(unames, function(uname){
    ## load sample features
    sampleFtrs <- collectCellFeatures(x, uname)
    sampleFtrs <- sampleFtrs[, selectedCellFtrs]
  
    ## negative controls on the same plate
    prw <- uname2prw(uname)
    negFtrs <- negFtrsPool[[intersect(which(x@plateList$Plate==prw$plate), which(x@plateList$Replicate==prw$replicate))]]
    negUnames=names(negFtrs)

    ## KS statistics for each feature (averaging multiple negative controls)
    sapply(1:ncol(sampleFtrs), function(col){
      KSstatistics <- sapply(negUnames, function(n){
        ks.test(sampleFtrs[,col], negFtrs[[n]][,col])$statistic
      })
      median(KSstatistics)
    })
  })
  profilesKS <- do.call(rbind, profilesKS)
  rownames(profilesKS) <- unames

  PDMByWellAvg(profiles=profilesKS, transformMethod='none', distMethod=distMethod)
}


PDMBySvmWeightVector <- function(x, unames, neg='rluc', selectedCellFtrs, distMethod=c('manhattan','euclidean', 'correlation','mahalanobis'), verbose=FALSE, kernel='linear', ...){
  
  if(missing(selectedCellFtrs)){
    ftrs <- collectCellFeatures(x, unames[1])
    selectedCellFtrs <- setdiff(colnames(ftrs),c('uname', 'spot', 'id', 'c.g.x', 'c.g.y', 'n.g.x', 'n.g.y'))
  }
  
  ## load negative control features
  negPlates <- uname2prw(unames)[,1:2]
  negPlates <- negPlates[!duplicated(negPlates),]
  negFtrsPool <- lapply(1:nrow(negPlates), function(row){
    negUnames <- getUnames(x, content=neg, plate=negPlates$plate[row], replicate=negPlates$replicate[row])
    negUnames <- setdiff(negUnames, getBadWells(x))
    negFtrs <- lapply(negUnames, function(uname){ 
      ftrs <- collectCellFeatures(x,uname)
      ftrs <- ftrs[, selectedCellFtrs]
    })
    names(negFtrs) <- negUnames
    negFtrs
  })
  
  profilesSvmWV <- lapply(unames, function(uname){
    ## load sample features
    sampleFtrs <- collectCellFeatures(x, uname)
    sampleFtrs <- sampleFtrs[, selectedCellFtrs]

    ## negative controls on the same plate
    prw <- uname2prw(uname)
    negFtrs <- negFtrsPool[[intersect(which(x@plateList$Plate==prw$plate), which(x@plateList$Replicate==prw$replicate))]]
    negUnames=names(negFtrs)

    ## for each neg well, calculate svm weight vector
    if(verbose) cat('calculating svm weight vector for', uname,'... ')
    svmWV <- lapply(negUnames, function(n){
      negFtrs <- negFtrs[[n]]
      negftrs <- negFtrs[, selectedCellFtrs]

      xts <- rbind(sampleFtrs, negFtrs)
      yts <- factor(c(rep(1, times=nrow(sampleFtrs)), rep(2, times=nrow(negFtrs))))
      model <- svm(x=xts, y=yts, kernel=kernel, ...)
      w <- t(model$coefs) %*% model$SV
      w
    })
    svmWV <- do.call(rbind, svmWV)
    rownames(svmWV) <- negUnames

    ## report median of all neg wells
    svmWV <- apply(svmWV, 2, median)
    if(verbose) cat('done.\n')
    svmWV
  })

  profilesSvmWV <- do.call(rbind, profilesSvmWV)
  rownames(profilesSvmWV) <- unames

  PDMByWellAvg(profiles=profilesSvmWV, transformMethod='none', distMethod=distMethod)
}


PDMBySvmAccuracy <- function(x, unames, selectedCellFtrs, cross=5, verbose=FALSE, ...){

  if(missing(selectedCellFtrs)){
    ftrs <- collectCellFeatures(x, unames[1])
    selectedCellFtrs <- setdiff(colnames(ftrs),c('uname', 'spot', 'id', 'c.g.x', 'c.g.y', 'n.g.x', 'n.g.y'))
  }

  ## load features for individual cells
  ftrsPool <- lapply(unames, function(uname){
    collectCellFeatures(x, uname)[,selectedCellFtrs]
  })
  names(ftrsPool) <- unames
  
  res <- lapply(1:length(unames), function(row){
    ftrs1 <- ftrsPool[[unames[row]]]

    sapply(1:length(unames), function(col){
      if(col == row) return(NA)
      ftrs2 <- ftrsPool[[unames[col]]]
      x <- rbind(ftrs1,ftrs2)
      y <- factor(c(rep(1, times=nrow(ftrs1)), rep(2,times=nrow(ftrs2))))
      if(verbose) cat('calculating svm dist. between', unames[row], 'and', unames[col], '... ')
      svmDist <- svm(x, y, cross=cross, ...)$tot.accuracy/100
      if(verbose) cat('done.\n')
      svmDist
    })
  })
  res <- do.call(rbind, res)
  rownames(res) <- colnames(res) <- unames
  res
}

