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
#' 

ORdensity <- setClass(
	"ORdensity",
	slots = c(positive="matrix", negative="matrix", labels="character", 
	          B="numeric", scale="logical", alpha="numeric", 
	          fold="numeric", probs="numeric", weights="numeric", K="numeric", 
	          out="list", OR="numeric", FP="numeric", dFP="numeric", char="data.frame", bestK = "numeric", 
	          verbose="logical", parallel="logical", replicable="logical", seed="numeric"),
	prototype = list(positive=matrix(), negative=matrix(), labels=character(), 
	                 B=numeric(), scale=logical(), alpha=numeric(), 
	                 fold=numeric(), probs=numeric(), weights=numeric(), K=numeric(), 
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
              clustering <- cluster::pam(distances::distances(scale(object@char)), object@bestK)$clustering
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

#' preclusteredData
#' 
#' This function returns the data after applying the ORdensity procedure, but before doing any clustering
#' 
#' @param object Object of type ORdensity
#' @return The data before doing any clustering
#' @examples
#' randomA <- matrix(rnorm(30*100), nrow = 100)
#' randomB <- matrix(rnorm(30*100), nrow = 100)
#' myORdensity <- new("ORdensity", Exp_cond_1 = randomA, Exp_cond_2 = randomB, parallel = TRUE)
#' preclusteredData(myORdensity)
#' @rdname preclusteredData
#' @docType methods
#' @export
setGeneric("preclusteredData", function(object, ...) standardGeneric("preclusteredData"))

setMethod("preclusteredData",
          signature = "ORdensity",
          definition = function(object, verbose=TRUE){
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
              if (verbose) {
                cat("Columns \"Strong\" and \"Flexible\" show the genes identified as DE genes\n")
                cat("They denote the strong selection (FP=0) with S and the flexible selection (FP < expectedFalsePositives) with F\n")
              }
              preclustered_data
          }
)

#' @title print
#' @param 
#' @return 
#' @examples
#' 
#' @rdname print
#' @export
setMethod("show",
           signature = "ORdensity",
           definition = function(object) {
             preClustering <- preclusteredData(object, verbose=FALSE)
             numGenes <- nrow(preClustering)
             cat("The ORdensity method has detected", numGenes, "suspicious genes\n", sep = " ")
             cat("The data before clustering is: \n")
             cat("Columns \"Strong\" and \"Flexible\" show the genes identified as DE genes\n")
             cat("They denote the strong selection (FP=0) with S and the flexible selection (FP < expectedFalsePositives) with F\n")
             print(preClustering)
           }
)

#' @title plotFPvsOR
#' @param 
#' @return 
#' @examples
#' 
#' @rdname plotFPvsOR
#' @export
setGeneric("plotFPvsOR", function(object, ...) standardGeneric("plotFPvsOR"))

setMethod("plotFPvsOR",
  signature = "ORdensity",
  definition = function(object, k = object@bestK){
    clustering <- cluster::pam(distances::distances(scale(object@char)), k)$clustering
    legend_text <- sprintf("cluster %s",seq(1:k))
    plot(object@FP, object@OR, type="n",main="Potential genes",xlab="FP",ylab="OR")
    points(object@FP, object@OR, cex=1/(0.5+object@dFP),  col = clustering)
    legend("topright", legend=legend_text, pch=16, col=unique(clustering))
    }
  )

#' @title silhouetteAnalysis
#' @param 
#' @return 
#' @examples
#' 
#' @rdname silhouetteAnalysis
#' @export
setGeneric("silhouetteAnalysis", function(object, ...) standardGeneric("silhouetteAnalysis"))

setMethod("silhouetteAnalysis",
  signature = "ORdensity",
  definition = function(object){
    s <- rep(NA, 10)
    for (k in 2:10)
    {
      aux <- cluster::pam(distances::distances(scale(object@char)), k)
      s[k] <- mean(cluster::silhouette(aux)[, "sil_width"])
    }
    plot(s, type="b", ylim=c(0,1), main="Clustering goodness", xlab = "K value", ylab = "silhouette")
  }
  )

#' @title clusplotk
#' @param 
#' @return 
#' @examples
#' 
#' @rdname clusplotk
#' @export
setGeneric("clusplotk", function(object, ...) standardGeneric("clusplotk"))

setMethod("clusplotk",
  signature = "ORdensity",
  definition = function(object, k = object@bestK){
    aa <- cluster::pam(distances::distances(scale(object@char)), k)
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
  vgR <- Rfast::med(distMatrix^2)/2
  I <- apply(distMatrix, 1,  IindexRobust, vg=vgR)
  OR <- 1/I
  return(OR)
}

#' @title compute.ORdensity
#' @param 
#' @return 
#' @examples
#' 
#' @rdname compute.ORdensity
#' @export
setGeneric("compute.ORdensity", function(object, ...) standardGeneric("compute.ORdensity"))

setMethod("compute.ORdensity",
	signature = "ORdensity",
	definition =  function(object, B=100, scale=FALSE, alpha=0.05, fold=floor(B/10), probs=c(0.25, 0.5, 0.75), weights=c(1/4,1/2,1/4), K = 10, verbose=FALSE, 
	                       parallel = FALSE, replicable = TRUE, seed = 0) {
	    a <- system.time ({
	    bootstrap_time_estimated <- FALSE
	      
	    positiveCases <- as.matrix(object@positive)
		  negativeCases <- as.matrix(object@negative)
		  numGenes <- dim(positiveCases)[1]
		  
		  cat("An object of size", format(object.size(1.0) * numGenes * numGenes / 7, unit="auto"), "is going to be created in memory. ")
		  cat("If the parallel option is enabled, as many objects of that size as the number of processors in your computer, ")
		  cat("are going to be created at the same time. Please consider that when running this code.\n")
		  
		  #cat("An object of size", format(object.size(1.0) * numGenes * length(probs) * B, unit="auto"), 
		  #    "is going to be created in memory. Please consider that when running this code. ")
		  
		  numPositiveCases <- dim(positiveCases)[2]
		  numNegativeCases <- dim(negativeCases)[2]
		  numCases <- numPositiveCases + numNegativeCases
		  numProbs <- length(probs)
		  numFolds <- fold})
	    
		  if (verbose) {print('Time after first chunk'); print(a)}

		  b <- system.time ({
		  quantilesDifferencesWeighted <- getQuantilesDifferencesWeighted(positiveCases, negativeCases, scale, weights, probs)
      })

		  if (verbose) {print('Time after second chunk'); print(b)}

		  c <- system.time ({Dxy <- distances::distances(quantilesDifferencesWeighted)})

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
		  cl <- parallel::makeForkCluster(nproc)
		  doParallel::registerDoParallel(cl)
      res_par <- foreach(b = 1:B, .combine = 'c', .options.RNG=seed) %dorng% {
        bootstrapSample <- getBootstrapSample(allCases, numPositiveCases)
  			res_one <- list()
  			res_one[[1]] <- getQuantilesDifferencesWeighted(bootstrapSample$positives, bootstrapSample$negatives, scale, weights, probs)
  			res_one[[2]] <- getOR(distances::distances(res_one[[1]]))
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
		    w <- system.time({
		      w1 <- system.time({
		    bootstrapSample <- getBootstrapSample(allCases, numPositiveCases)})
		      if (verbose) {print('Time after a non-parallel bootstrap replication (step 1)'); print(w1)}
		      w2 <- system.time({
		    quantilesDifferencesWeighted.null[ , ,b] <- getQuantilesDifferencesWeighted(bootstrapSample$positives, bootstrapSample$negatives, scale, weights, probs)})
		      if (verbose) {print('Time after a non-parallel bootstrap replication (step 2)'); print(w2)}
		      w3 <- system.time({
		       ORbootstrap[, b] <- getOR(distances::distances(quantilesDifferencesWeighted.null[,,b]))
		      if (verbose) {print('Time after a non-parallel bootstrap replication (step 3)'); print(w3)}
		    })
		    })
		    if (!bootstrap_time_estimated)
		    {
		      bootstrap_time_estimated <- TRUE
		      bootstrap_time <- w['elapsed']
		      cat("A bootstrap replication takes", bootstrap_time, "seconds, and you have requested", B, "bootstrap replications.\n")
		    }
		    if (verbose) {print('Time after a non-parallel bootstrap replication'); print(w)}
		  }
		}
   })
		  if (verbose) {print('Time after fourth chunk'); print(d)}

		# OR values for original data
		  e <- system.time ({
		  ORoriginal <- getOR(Dxy)

		# Find cut point
		   cutPoint <- (sort(c(ORbootstrap)))[floor((1-alpha)*numGenes*B)]

		# Find individuals beyond threshold
		   suspicious <- ORoriginal > cutPoint
		   numSuspicious <- sum(suspicious)

		   # the indices are in the form (case, bootstrap_sample)
		  indicesBiDimORbootstrapBeyondCutPoint <- which(ORbootstrap > cutPoint, arr.ind=TRUE)
		  
		  numORbootstrapBeyondCutPoint <- dim(indicesBiDimORbootstrapBeyondCutPoint)[1]
		  
		#
		  # vector of integers (assigned fold) of size numGenes * B * alpha
		  assignFoldToBootstrapBeyondCutPoint <- sample(1:numFolds, numORbootstrapBeyondCutPoint, replace=TRUE) #
		  
		  # create a zero-filled 3D matrix with a 2D matrix of dim (numSuspicious, 3) for each fold
		  # created to store FPneighbourghood, densityFP and radius
		  originalDataFPStatistics <- array(0, dim=c(numSuspicious, 3, numFolds)) })

		  if (verbose) {print('Time after fifth chunk'); print(e)}

		  globalquantilesDifferencesWeighted.null <<- quantilesDifferencesWeighted.null
		  
		  f <- system.time ({
		  for (j in 1:numFolds)
		  {
		    # for every fold, we see how is the distribution of the OR statistic along the boostrap samples and the original ones
		      currentFold <- assignFoldToBootstrapBeyondCutPoint == j
		      numInCurrentFold <- sum(currentFold)
		      quantilesBootstrapFold <- matrix(0, nrow=numInCurrentFold, ncol=numProbs)
		      cont <- 1
		      indicesBeyondCutPointCurrentFold <- (1:numORbootstrapBeyondCutPoint)[currentFold]
		      for (i in indicesBeyondCutPointCurrentFold)
		      {
  			   numGene <- indicesBiDimORbootstrapBeyondCutPoint[i, 1]
  			   numBootstrap <- indicesBiDimORbootstrapBeyondCutPoint[i, 2]
  			   quantilesBootstrapFold[cont, ] <- quantilesDifferencesWeighted.null[numGene, , numBootstrap]
  			   cont <- cont + 1
		      }
		      quantilesOriginalPlusBootstrapFold <- rbind(quantilesDifferencesWeighted[suspicious, ], quantilesBootstrapFold)
		      # after joining the original data with the bootstrap, we need the labels to find which is which
		      label <- c(rep(1, numSuspicious), rep(0, numInCurrentFold))

		      Dmix <- distances::distances(quantilesOriginalPlusBootstrapFold)
		      Dmix <- as.matrix(Dmix)
		      originalDataFPStatisticsByFold <- matrix(0, nrow=numSuspicious, ncol=3)
		      colnames(originalDataFPStatisticsByFold) <- c( "FPneighbourghood", "densityFP", "radius")

		      DOriginal <- Dmix[1:numSuspicious, ]
		      for(i in 1:numSuspicious)
		      {
		        originalDataFPStatisticsByFold[i, ] <- c(density(DOriginal[i, -i], label=label[-i], K))
		      }
		      originalDataFPStatistics[ , , j] <- originalDataFPStatisticsByFold
		  } # end for (j in 1:numFolds)
		    })

		  if (verbose) {print('Time after sixth chunk'); print(f)}

		  g <- system.time ({
		    # means with respect to the folds
		    originalDataFPStatisticsMeans <- t(plyr::aaply(originalDataFPStatistics, c(2,1), mean))
		    originalDataFPNeighboursStats <- t(apply(originalDataFPStatistics[, 1, ], 1, function(x){c(min(x), mean(x), max(x))}))
		    percentageSuspiciousOverPositives <- numSuspicious/(numSuspicious+numORbootstrapBeyondCutPoint/numFolds)
		    percentageBoostrapOverPositives <- (numORbootstrapBeyondCutPoint/numFolds)/(numSuspicious+numORbootstrapBeyondCutPoint/numFolds)
		    
		    diffOverExpectedFPNeighbours <- originalDataFPStatisticsMeans[, 1] - percentageBoostrapOverPositives * K
		    ####
		    labelGenes <- object@labels
		    # genes <- (1:numGenes)[suspicious]
		    genes <- labelGenes[suspicious]
		    # print(genes)
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
                                              fold=floor(B/10), probs=c(0.25, 0.5, 0.75), weights=c(1/4,1/2,1/4), K = 10, 
                                              out, OR, FP, dFP, char, bestK, verbose = FALSE, 
                                              parallel = FALSE, replicable = TRUE, seed = 0) {
  .Object@positive <- Exp_cond_1
  .Object@negative <- Exp_cond_2
  # validObject(.Object)
  if (is.null(labels))
  { 
    .Object@labels<- paste("Gene", 1:nrow(Exp_cond_1), sep="")
  }
  else
  {
    .Object@labels <- labels
  }
  .Object@B <- B
  .Object@scale <- scale
  .Object@alpha <- alpha
  .Object@fold <- fold
  .Object@probs <- probs
  .Object@weights <- weights
  .Object@K <- K
  .Object@verbose <- verbose
  .Object@parallel <- parallel
  .Object@replicable <- replicable
  .Object@seed <- seed
  .Object@out <- compute.ORdensity(.Object, B = .Object@B, scale = .Object@scale, alpha = .Object@alpha, fold = .Object@fold, 
                                   probs = .Object@probs, weights = .Object@weights, K = .Object@K, 
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

#' @title findbestK
#' @param 
#' @return 
#' @examples
#' 
#' @rdname findbestK
#' @export
setGeneric("findbestK", function(object, ...) standardGeneric("findbestK"))

setMethod("findbestK",
          signature = "ORdensity",
          definition = function(object){
            s <- rep(NA, 10)
            # len(object@char) could be less than 10
            for (k in 2:10)
            {
              shit <<- object@char
              aux <- cluster::pam(distances::distances(scale(object@char)), k)
              s[k] <- mean(cluster::silhouette(aux)[, "sil_width"])
            }
            best_k <- which(s == max(s, na.rm = TRUE))
            return (best_k)
          }
)

#' @title findDEgenes
#' @param 
#' @return 
#' @examples
#' 
#' @rdname findDEgenes
#' @export
setGeneric("findDEgenes", function(object, ...) standardGeneric("findDEgenes"))

setMethod("findDEgenes",
          signature = "ORdensity",
          definition = function(object, numclusters=NULL, verbose=FALSE){
            KForClustering <- object@bestK
            if (!is.null(numclusters))
            {
              KForClustering <- numclusters
            }
            clustering <- cluster::pam(distances::distances(scale(object@char)), KForClustering)$clustering
            result_prov <- list()
            meanOR <- list()
            for (k in 1:KForClustering)
            {
              result_prov[[k]] <- object@out$summary[clustering==k,]
              meanOR[[k]] <- mean(result_prov[[k]][,'OR'])
            }
            clusters_ordering <- order(as.numeric(meanOR), decreasing = TRUE)
            clusters <- list()
            for (k in 1:KForClustering)
            {
              clusters[[k]] <- result_prov[[clusters_ordering[k]]]
            }
            DFgenes <- list()
            for (k in 1:KForClustering)
            {
              # DFgenes[[k]] <- list("cluster_number"=k, "genes"=sort(clusters[[k]][,'id']), "meanOR"=mean(clusters[[k]][,'OR']))
              if (verbose)
              {
                DFgenes[[k]] <- list("cluster_number"=k, "numberOfGenes"=length(clusters[[k]][,'id']), 
                                     "meanOR"=mean(clusters[[k]][,'OR']), "genes"=sort(clusters[[k]][,'id']))
              }
              else
              {
                DFgenes[[k]] <- list("cluster_number"=k, "numberOfGenes"=length(clusters[[k]][,'id']), "meanOR"=mean(clusters[[k]][,'OR']))
              }
            }
            cat("The ORdensity method has found that the optimal clustering of the data consists of",object@bestK,"clusters\n")
            if (!is.null(numclusters))
            {
              cat("The user has chosen a clustering of",numclusters,"clusters\n")
            }
            # return(list("DFgenes"=DFgenes,"clusters"=clusters))
            return(DFgenes)
          }
)
