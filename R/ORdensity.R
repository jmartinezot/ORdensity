#' @title
#' Class for representing ORdensity
#'
#' @description
#' Class for representing ORdensity
#'
#' @name ORdensity
#' @rdname ORdensity
#' @exportClass ORdensity
#'
#' @author Itziar Irigoien, Concepcion Arenas, Jose Maria Martinez-Otzeta
#'
#' @export plotFPvsOR
#' @export silhouetteAnalysis
#' @export clusplotk
#' @export compute.ORdensity
#' @export findDEgenes
#' @export summary
#' @export preclusteredData

ORdensity <- setClass(
	"ORdensity",
	slots = c(positive="matrix", negative="matrix", labels="character", 
	          B="numeric", scale="logical", alpha="numeric", 
	          fold="numeric", weights="numeric", K="numeric", 
	          out="list", OR="numeric", FP="numeric", dFP="numeric", char="data.frame", bestK = "numeric", 
	          verbose="logical", parallel="logical", replicable="logical", seed="numeric"),
	prototype = list(positive=matrix(), negative=matrix(), labels=character(), 
	                 B=numeric(), scale=logical(), alpha=numeric(), 
	                 fold=numeric(), weights=numeric(), K=numeric(), 
	                 out=list(), OR=numeric(), FP=numeric(), dFP=numeric(), char=data.frame(), bestK = numeric(), 
	                 verbose=logical(), parallel=logical(), replicable=logical(), seed=numeric())
)

# setGeneric("summary")

# setGeneric("summary.ORdensity", function(object, ...) standardGeneric("summary.ORdensity"))

#' @title summary
#' @param 
#' @return 
#' @examples
#' 
#' @rdname summary
#' @export
setMethod("summary",
          signature = "ORdensity",
          definition = function(object){
              clustering <- pam(dist(scale(object@char)), object@bestK)$clustering
              result_prov <- list()
              meanOR <- list()
              for (k in 1:object@bestK)
              {
                result_prov[[k]] <- object@out$summary[clustering==k,]
                result_prov[[k]] <- result_prov[[k]][,!colnames(result_prov[[k]]) %in% c("DifExp",  "minFP", "maxFP", "radius")]
                meanOR[[k]] <- mean(result_prov[[k]][,'OR'])
              }
              clusters_ordering <- order(as.numeric(meanOR), decreasing = TRUE)
              clusters <- list()
              for (k in 1:object@bestK)
              {
                clusters[[k]] <- result_prov[[clusters_ordering[k]]]
              }
              # cat("The ORdensity method has found that the optimal clustering of the data consists of",object@bestK,"clusters\n")
              prop <- object@out$prop
              neighbours <- prop[3]
              p0 <- prop[2]
              cat("This is the proposed clustering made by the ORdensity method\n")
              cat("For the computation of FP and dFP a total of", neighbours, "neighbours have been taken into account\n")
              cat("The expected number of false positives neighbours is", p0*neighbours, "\n")
              cat("The ORdensity method has found that the optimal clustering of the data consists of",object@bestK,"clusters\n\n")
              return(list("neighbours"=neighbours, "expectedFalsePositiveNeighbours"=p0*neighbours, "clusters"=clusters))
        }
)

setGeneric("preclusteredData", function(object, ...) standardGeneric("preclusteredData"))

#' @title preclusteredData
#' @param 
#' @return 
#' @examples
#' 
#' @rdname preclusteredData
#' @export
setMethod("preclusteredData",
          signature = "ORdensity",
          definition = function(object){
              prop <- object@out$prop
              neighbours <- prop[3]
              p0 <- prop[2]
              preclustered_data <- as.data.frame(object@out$summary)
              preclustered_data$DifExp <- NULL
              preclustered_data$minFP <- NULL
              preclustered_data$maxFP <- NULL
              preclustered_data$radius <- NULL
              preclustered_data$Strong <- ifelse(preclustered_data$FP == 0, "S", "")
              preclustered_data$Flexible <- ifelse(preclustered_data$FP < p0 * neighbours, "F", "")
              cat("Columns \"Strong\" and \"Flexible\" show the genes identified as DE genes\n")
              cat("They denote the strong selection (FP=0) with S and the flexible selection (FP < expectedFalsePositives) with F\n")
              preclustered_data
          }
)

#setMethod("show",
#          signature = "ORdensity",
#          definition = function(object) {
#            cat("Positive data matrix: ", object@positive, "\n", sep = " ")
#            cat("Negative data matrix: ", object@negative, "\n", sep = " ")
#          }
#)

setGeneric("plotFPvsOR", function(object, ...) standardGeneric("plotFPvsOR"))

#' @title plotFPvsOR
#' @param 
#' @return 
#' @examples
#' 
#' @rdname plotFPvsOR
#' @export
setMethod("plotFPvsOR",
  signature = "ORdensity",
  definition = function(object, k = object@bestK){
    clustering <- pam(dist(scale(object@char)), k)$clustering
    legend_text <- sprintf("cluster %s",seq(1:k))
    plot(object@FP, object@OR, type="n",main="Potential genes",xlab="FP",ylab="OR")
    points(object@FP, object@OR, cex=1/(0.5+object@dFP),  col = clustering)
    legend("topright", legend=legend_text, pch=16, col=unique(clustering))
    }
  )

setGeneric("silhouetteAnalysis", function(object, ...) standardGeneric("silhouetteAnalysis"))

#' @title silhouetteAnalysis
#' @param 
#' @return 
#' @examples
#' 
#' @rdname silhouetteAnalysis
#' @export
setMethod("silhouetteAnalysis",
  signature = "ORdensity",
  definition = function(object){
    library(cluster)
    s <- rep(NA, 10)
    for (k in 2:10)
    {
      aux <- pam(dist(scale(object@char)), k)
      s[k] <- mean(silhouette(aux)[, "sil_width"])
    }
    plot(s, type="b", ylim=c(0,1), main="Clustering goodness", xlab = "K value", ylab = "silhouette")
  }
  )

setGeneric("clusplotk", function(object, ...) standardGeneric("clusplotk"))

#' @title clusplotk
#' @param 
#' @return 
#' @examples
#' 
#' @rdname clusplotk
#' @export
setMethod("clusplotk",
  signature = "ORdensity",
  definition = function(object, k = object@bestK){
    aa <- pam(dist(scale(object@char)), k)
    clusplot(aa, main = paste("Clustering with k = ", k))
  }
  )

getQuantilesDifferencesWeighted <- function(positiveCases, negativeCases, scale, weights, probs){
  numGenes <- dim(positiveCases)[1]
  quantilesPositiveCases <- t(apply(positiveCases, 1, quantile, probs=probs))
  quantilesNegativeCases <- t(apply(negativeCases, 1, quantile, probs=probs))
  quantilesDifferences <- cbind(quantilesPositiveCases - quantilesNegativeCases)
  numProbs <- length(probs)
  if(scale){
    interquartileRangePositiveCases <- quantilesPositiveCases[,3]-quantilesPositiveCases[,1]
    interquartileRangeNegativeCases  <- quantilesNegativeCases[,3]-quantilesNegativeCases[,1]
    maxInterquartileRange <- apply(cbind(interquartileRangePositiveCases, interquartileRangeNegativeCases), 1, max)
    if(any(maxInterquartileRange==0))
    {stop('Can\'t scale the data')}
    quantilesDifferences <- quantilesDifferences/maxInterquartileRange
  }
  quantilesDifferencesWeighted <- quantilesDifferences*matrix(rep(weights, numGenes), byrow=TRUE, ncol=numProbs)
  return (quantilesDifferencesWeighted)
}

getBootstrapSample <- function(allCases, numPositiveCases)
{
  numCases <- dim(allCases)[2]
  s1 <- sample(1:numCases, numPositiveCases, replace=FALSE)
  s2 <- (1:numCases)[-s1]
  aux1 <- allCases[, s1]
  aux2 <- allCases[, s2]
  return(list("positives"=aux1, "negatives"=aux2))
}

getOR <- function(distObject)
{	
  distMatrix <- as.matrix(distObject)
  vgR <- median(distMatrix^2)/2
  I <- apply(distMatrix, 1,  IindexRobust, vg=vgR)
  OR <- 1/I
  return(OR)
}

setGeneric("compute.ORdensity", function(object, ...) standardGeneric("compute.ORdensity"))

#' @title compute.ORdensity
#' @param 
#' @return 
#' @examples
#' 
#' @rdname compute.ORdensity
#' @export
setMethod("compute.ORdensity",
	signature = "ORdensity",
	definition =  function(object, B=100, scale=FALSE, alpha=0.05, fold=floor(B/10), weights=c(1/4,1/2,1/4), K = 10, verbose=FALSE, parallel = FALSE, replicable = TRUE, seed = 0) {
	    a <- system.time ({
	    positiveCases <- as.matrix(object@positive)
		  negativeCases <- as.matrix(object@negative)
		  numGenes <- dim(positiveCases)[1]
		  numPositiveCases <- dim(positiveCases)[2]
		  numNegativeCases <- dim(negativeCases)[2]
		  numCases <- numPositiveCases + numNegativeCases
		  probs=c(0.25, 0.5, 0.75)
		  numProbs <- length(probs)
		  numFolds <- fold})

		  if (verbose) {print('Time after first chunk'); print(a)}

		  b <- system.time ({
		  quantilesDifferencesWeighted <- getQuantilesDifferencesWeighted(positiveCases, negativeCases, scale, weights, probs)
      })

		  if (verbose) {print('Time after second chunk'); print(b)}

		  c <- system.time ({Dxy <- dist(quantilesDifferencesWeighted)})

		  if (verbose) {print('Time after third chunk'); print(c)}

		  d <- system.time ({
		  allCases <- cbind(positiveCases, negativeCases)
		  ORbootstrap <- matrix(0, nrow=numGenes, ncol=B)
		  quantilesDifferencesWeighted.null <- array(0, dim=c(numGenes, numProbs, B))

		  if (parallel)
		  {
		    require(foreach)
		    if (replicable){
		      require(doRNG)
		      set.seed(seed)
		    }

		  # require(bigstatsr)
	    nproc <- parallel::detectCores()
		  cl <- parallel::makeCluster(nproc)
		  doParallel::registerDoParallel(cl)
      res_par <- foreach(b = 1:B, .combine = 'c', .options.RNG=seed) %dorng% {
        bootstrapSample <- getBootstrapSample(allCases, numPositiveCases)
  			res_one <- list()
  			res_one[[1]] <- getQuantilesDifferencesWeighted(bootstrapSample$positives, bootstrapSample$negatives, scale, weights, probs)
  			res_one[[2]] <- getOR(dist(res_one[[1]]))
  			res_one
		  }
      parallel::stopCluster(cl)

      for (b in 1:B) {
        quantilesDifferencesWeighted.null[ , ,b] <- res_par[[b*2-1]]
        ORbootstrap[, b] <- res_par[[b*2]]
      }
		} # end if
		else
		{
		  if (replicable){
        set.seed(seed)
		  }
		  for (b in 1:B)
		  {
		    bootstrapSample <- getBootstrapSample(allCases, numPositiveCases)
		    quantilesDifferencesWeighted.null[ , ,b] <- getQuantilesDifferencesWeighted(bootstrapSample$positives, bootstrapSample$negatives, scale, weights, probs)
		    ORbootstrap[, b] <- getOR(dist(quantilesDifferencesWeighted.null[,,b]))
		  }
		}
		  })

		  glquantilesDifferencesWeighted.null <<- quantilesDifferencesWeighted.null
		  if (verbose) {print('Time after fourth chunk'); print(d)}

		# OR values for original data
		  e <- system.time ({
		  ORoriginal <- getOR(Dxy)

		# Find cut point
		   cutPoint <- (sort(c(ORbootstrap)))[floor((1-alpha)*numGenes*B)]

		# Find individuals beyond threshold
		   suspicious <- ORoriginal > cutPoint
		   numSuspicious <- sum(suspicious)
		   
		   glORbootstrap <<- ORbootstrap
		   glORoriginal <<- ORoriginal
		   glcutPoint <<- cutPoint
		   glsuspicious <<- suspicious
		   glnumSuspicious <<- numSuspicious

		   # the indices are in the form (case, bootstrap_sample)
		  indicesBiDimORbootstrapBeyondCutPoint <- which(ORbootstrap > cutPoint, arr.ind=TRUE)
		  
		  numORbootstrapBeyondCutPoint <- dim(indicesBiDimORbootstrapBeyondCutPoint)[1]
		  
		  glindicesBiDimORbootstrapBeyondCutPoint <<- indicesBiDimORbootstrapBeyondCutPoint
		  glnumORbootstrapBeyondCutPoint <<- numORbootstrapBeyondCutPoint
		#
		  # vector of integers (assigned fold) of size numGenes * B * alpha
		  assignFoldToBootstrapBeyondCutPoint <- sample(1:numFolds, numORbootstrapBeyondCutPoint, replace=TRUE) #
		  glassignFoldToBootstrapBeyondCutPoint <<- assignFoldToBootstrapBeyondCutPoint
		  
		  # create a zero-filled 3D matrix with a 2D matrix of dim (numSuspicious, 3) for each fold
		  # created to store FPneighbourghood, densityFP and radius
		  originalDataFPStatistics <- array(0, dim=c(numSuspicious, 3, numFolds)) })

		  if (verbose) {print('Time after fifth chunk'); print(e)}

		  f <- system.time ({
		  for (j in 1:numFolds)
		  {
		    # for every fold, we see how is the distribution of the OR statistic along the boostrap samples and the original ones
		      currentFold <- assignFoldToBootstrapBeyondCutPoint == j
		      numInCurrentFold <- sum(currentFold)
		      quantilesBootstrapFold <- matrix(0, nrow=numInCurrentFold, ncol=numProbs)
		      glcurrentFold <<- currentFold
		      glnumInCurrentFold <<- numInCurrentFold
		      glquantilesBootstrapFoldOriginal <<- quantilesBootstrapFold
		      cont <- 1
		      indicesBeyondCutPointCurrentFold <- (1:numORbootstrapBeyondCutPoint)[currentFold]
		      for (i in indicesBeyondCutPointCurrentFold)
		      {
  			   numGene <- indicesBiDimORbootstrapBeyondCutPoint[i, 1]
  			   numBootstrap <- indicesBiDimORbootstrapBeyondCutPoint[i, 2]
  			   quantilesBootstrapFold[cont, ] <- quantilesDifferencesWeighted.null[numGene, , numBootstrap]
  			   cont <- cont + 1
		      }
		      glquantilesBootstrapFoldAfterLoop <<- quantilesBootstrapFold
		      quantilesOriginalPlusBootstrapFold <- rbind(quantilesDifferencesWeighted[suspicious, ], quantilesBootstrapFold)
		      glquantilesOriginalPlusBootstrapFold <<- quantilesOriginalPlusBootstrapFold
		      # after joining the original data with the bootstrap, we need the labels to find which is which
		      label <- c(rep(1, numSuspicious), rep(0, numInCurrentFold))
		      gllabel <<- label

		      Dmix <- dist(quantilesOriginalPlusBootstrapFold)
		      Dmix <- as.matrix(Dmix)
		      originalDataFPStatisticsByFold <- matrix(0, nrow=numSuspicious, ncol=3)
		      colnames(originalDataFPStatisticsByFold) <- c( "FPneighbourghood", "densityFP", "radius")

		      DOriginal <- Dmix[1:numSuspicious, ]
		      for(i in 1:numSuspicious)
		      {
		        originalDataFPStatisticsByFold[i, ] <- c(density(DOriginal[i, -i], label=label[-i], K))
		      }
		      glDOriginal <<- DOriginal
		      gloriginalDataFPStatisticsByFold <<- originalDataFPStatisticsByFold
		      originalDataFPStatistics[ , , j] <- originalDataFPStatisticsByFold
		  } # end for (j in 1:numFolds)
		    })

		  gloriginalDataFPStatistics <<- originalDataFPStatistics
		  if (verbose) {print('Time after sixth chunk'); print(f)}

		  g <- system.time ({
		    # means with respect to the folds
		    originalDataFPStatisticsMeans <- t(plyr::aaply(originalDataFPStatistics, c(2,1), mean))
		    originalDataFPNeighboursStats <- t(apply(originalDataFPStatistics[, 1, ], 1, function(x){c(min(x), mean(x), max(x))}))
		    percentageSuspiciousOverPositives <- numSuspicious/(numSuspicious+numORbootstrapBeyondCutPoint/numFolds)
		    percentageBoostrapOverPositives <- (numORbootstrapBeyondCutPoint/numFolds)/(numSuspicious+numORbootstrapBeyondCutPoint/numFolds)
		    
		    gloriginalDataFPStatisticsMeans <<- originalDataFPStatisticsMeans
		    gloriginalDataFPNeighboursStats <<- originalDataFPNeighboursStats
		    glpercentageSuspiciousOverPositives <<- percentageSuspiciousOverPositives
		    glpercentageBoostrapOverPositives <<- percentageBoostrapOverPositives
		    
		    diffOverExpectedFPNeighbours <- originalDataFPStatisticsMeans[, 1] - percentageBoostrapOverPositives * K
		    ####
		    glnumGenes <<- numGenes
		    glsuspicious <<- suspicious
		    labelGenes <- object@labels
		    # genes <- (1:numGenes)[suspicious]
		    genes <- labelGenes[suspicious]
		    print(genes)
		    #finalResult <- cbind(genes, ORoriginal[suspicious], diffOverExpectedFPNeighbours, originalDataFPNeighboursStats, originalDataFPStatisticsMeans[, -1])
		    #row.names(finalResult) <- NULL
		    #colnames(finalResult) <- c("id", "OR", "DifExp",  "minFP", "FP", "maxFP", "dFP", "radius")
		    #glfinalResult <<- finalResult
		    #finalOrdering <- order(finalResult[, 3], -finalResult[, 2])
		    finalResult <- data.frame("id"=genes, "OR"=ORoriginal[suspicious], "DifExp"=diffOverExpectedFPNeighbours,
		                              "minFP"=originalDataFPNeighboursStats[,1], "FP"= originalDataFPNeighboursStats[,2],
		                              "maxFP"=originalDataFPNeighboursStats[,3], "dFP"=originalDataFPStatisticsMeans[, 2],
		                              "radius"=originalDataFPStatisticsMeans[, 3])
		    finalResult$id <- as.character(finalResult$id)
		    row.names(finalResult) <- NULL
		    finalResult <<- finalResult
		    finalOrdering <- order(finalResult[, 3], -finalResult[, 2])
		   # print(res[oo,])
		   # print(numSuspicious)
		   # print(c(ps, p0))
		   # print(list("summary"=res[oo, ], "ns"=numSuspicious, "prop"=c(ps, p0)))
		   })

		   if (verbose) {print('Time after seventh chunk'); print(g)}

		   object@out <- list("summary"=finalResult[finalOrdering, ], "ns"=numSuspicious, "prop"=c(percentageSuspiciousOverPositives, percentageBoostrapOverPositives , K))
		   # object
		}
	)

# a <- ORdensity(positive = matrix(rnorm(500), nrow=50, ncol=10), negative = matrix(rnorm(500), nrow=50, ncol=10))
# a <- compute.ORdensity(a)

setValidity("ORdensity", function(object) {
  valid <- TRUE
  msg <- NULL
  if (length(object@positive) == 0) {
    valid <- FALSE
    msg <- c(msg, "There is no positive data")
  }
  if (length(object@negative) == 0) {
    valid <- FALSE
    msg <- c(msg, "There is no negative data")
  }
  if (nrow(object@positive) != nrow(object@negative)) {
    valid <- FALSE
    msg <- c(msg, "Positive and negative number of rows do not match")
  }
  if (valid) TRUE else msg
}
)

setMethod("initialize", "ORdensity", function(.Object, Exp_cond_1, Exp_cond_2, labels = NULL, B=100, scale=FALSE, alpha=0.05, 
                                              fold=floor(B/10), weights=c(1/4,1/2,1/4), K = 10, 
                                              out, OR, FP, dFP, char, bestK, verbose = FALSE, 
                                              parallel = FALSE, replicable = TRUE, seed = 0) {
  .Object@positive <- Exp_cond_1
  .Object@negative <- Exp_cond_2
  # validObject(.Object)
  if (is.null(labels))
  { 
    .Object@labels<- paste("Gene", 1:nrow(positive), sep="")
  }
  else
  {
    .Object@labels <- labels
  }
  .Object@B <- B
  .Object@scale <- scale
  .Object@alpha <- alpha
  .Object@fold <- fold
  .Object@weights <- weights
  .Object@K <- K
  .Object@verbose <- verbose
  .Object@parallel <- parallel
  .Object@replicable <- replicable
  .Object@seed <- seed
  .Object@out <- compute.ORdensity(.Object, B = .Object@B, scale = .Object@scale, alpha = .Object@alpha, fold = .Object@fold, 
                                   weights = .Object@weights, K = .Object@K, 
                                   verbose = .Object@verbose, parallel = .Object@parallel, replicable = .Object@replicable, 
                                   seed = .Object@seed)
  .Object@OR <- .Object@out$summary[, "OR"]
  .Object@FP <- .Object@out$summary[, "FP"]
  .Object@dFP <- .Object@out$summary[, "dFP"]
  .Object@char <- data.frame(.Object@OR, .Object@FP, .Object@dFP)
  require(cluster)
  .Object@bestK <- findbestK(.Object)
  .Object
})

setGeneric("findbestK", function(object, ...) standardGeneric("findbestK"))

#' @title findbestK
#' @param 
#' @return 
#' @examples
#' 
#' @rdname findbestK
#' @export
setMethod("findbestK",
          signature = "ORdensity",
          definition = function(object){
            s <- rep(NA, 10)
            for (k in 2:10)
            {
              aux <- pam(dist(scale(object@char)), k)
              s[k] <- mean(silhouette(aux)[, "sil_width"])
            }
            best_k <- which(s == max(s, na.rm = TRUE))
            return (best_k)
          }
)

setGeneric("findDEgenes", function(object, ...) standardGeneric("findDEgenes"))

#' @title findDEgenes
#' @param 
#' @return 
#' @examples
#' 
#' @rdname findDEgenes
#' @export
setMethod("findDEgenes",
          signature = "ORdensity",
          definition = function(object){
            clustering <- pam(dist(scale(object@char)), object@bestK)$clustering
            result_prov <- list()
            meanOR <- list()
            for (k in 1:object@bestK)
            {
              result_prov[[k]] <- object@out$summary[clustering==k,]
              meanOR[[k]] <- mean(result_prov[[k]][,'OR'])
            }
            clusters_ordering <- order(as.numeric(meanOR), decreasing = TRUE)
            clusters <- list()
            for (k in 1:object@bestK)
            {
              clusters[[k]] <- result_prov[[clusters_ordering[k]]]
            }
            DFgenes <- list()
            for (k in 1:object@bestK)
            {
              # DFgenes[[k]] <- list("cluster_number"=k, "genes"=sort(clusters[[k]][,'id']), "meanOR"=mean(clusters[[k]][,'OR']))
              DFgenes[[k]] <- list("cluster_number"=k, "numberOfGenes"=length(clusters[[k]][,'id']), "meanOR"=mean(clusters[[k]][,'OR']))
            }
            cat("The ORdensity method has found that the optimal clustering of the data consists of",object@bestK,"clusters\n")
            # return(list("DFgenes"=DFgenes,"clusters"=clusters))
            return(DFgenes)
          }
)
