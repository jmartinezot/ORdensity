#' @title
#' Class for representing ORdensity
#'
#' @description
#' Class for representing ORdensity
#'
#' @name ORdensity
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
	slots = c(positive="matrix", negative="matrix", out="list", OR="numeric", FP="numeric", dFP="numeric", char="data.frame", bestK = "numeric", verbose="logical", parallel="logical", replicable="logical", seed="numeric"),
	prototype = list(positive=matrix(), negative=matrix(), out=list(), OR=numeric(), FP=numeric(), dFP=numeric(), char=data.frame(), bestK = numeric(), verbose=logical(), parallel=logical(), replicable=logical(), seed=numeric())
)

# setGeneric("summary")

# setGeneric("summary.ORdensity", function(object, ...) standardGeneric("summary.ORdensity"))

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
              cat("Column Strong denotes the cases when FP=0\n")
              cat("Column Flexible denotes the cases when FP < expectedFalsePositives\n")
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
		  numProbs <- length(probs)})

		  if (verbose) {print('Time after first chunk'); print(a)}

		  b <- system.time ({
		  quantilesDifferencesWeighted <- getQuantilesDifferencesWeighted(positiveCases, negativeCases, scale, weights, probs)
      })

		  if (verbose) {print('Time after second chunk'); print(b)}

		  c <- system.time ({Dxy <- dist(quantilesDifferencesWeighted)})

		  if (verbose) {print('Time after third chunk'); print(c)}

		  d <- system.time ({
		  allCases <- cbind(positiveCases, negativeCases)
		  ORnull <- matrix(0, nrow=numGenes, ncol=B)
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
  			res <- list()
  			res[[1]] <- getQuantilesDifferencesWeighted(bootstrapSample$positives, bootstrapSample$negatives, scale, weights, probs)
  			res[[2]] <- getOR(dist(res[[1]]))
  			res
		  }
      parallel::stopCluster(cl)

      for (b in 1:B) {
        quantilesDifferencesWeighted.null[ , ,b] <- res_par[[b*2-1]]
        ORnull[, b] <- res_par[[b*2]]
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
		    ORnull[, b] <- getOR(dist(quantilesDifferencesWeighted.null[,,b]))
		  }
		}
		  })

		  if (verbose) {print('Time after fourth chunk'); print(d)}

		# OR values for original data
		  e <- system.time ({
		  OR <- getOR(Dxy)

		# Find cut point
		   corte <- (sort(c(ORnull)))[floor((1-alpha)*numGenes*B)]

		# Find individuals beyond threshold
		   suspicious <- OR > corte
		   ns <- sum(suspicious)

		  selec.null <- which(ORnull > corte, arr.ind=TRUE)
		  nnull <- dim(selec.null)[1]
		#
		  mirar <- sample(1:fold, nnull, replace=TRUE) #
		  # print(ns)
		  # print(fold)
		  apilar <- array(0, dim=c(ns, 3, fold)) })

		  if (verbose) {print('Time after fifth chunk'); print(e)}

		# Que pasa en la capa 1
		  f <- system.time ({
		  for (j in 1:fold)
		  {
		      capa <- mirar == j
		      nnull.red <- sum(capa)
		      dat <- matrix(0, nrow=nnull.red, ncol=numProbs)
		      cont <- 1
		      for (i in (1:nnull)[capa])
		      {
			   f <- selec.null[i, 1]
			   b <- selec.null[i, 2]
			   dat[cont, ] <- quantilesDifferencesWeighted.null[f, , b]
			   cont <- cont + 1
		       }
		      dat <- rbind(quantilesDifferencesWeighted[suspicious, ], dat)
		      label <- c(rep(1, ns), rep(0, nnull.red))

		    #  Dmix <- ICGE::dgower(dat, type=list(cuant=1:p))
		    # dat <- dat*matrix(rep(weights, ns+nnull.red), byrow=TRUE, ncol=p) NO HAY QUE VOLVER A PONER LOS PESOS!!!
		      Dmix <- dist(dat)
		      Dmix <- as.matrix(Dmix)
		      res <- matrix(0, nrow=ns, ncol=3)
		      colnames(res) <- c( "FPneighbourghood", "densityFP", "radius")


		      Dred <- Dmix[1:ns, ]
		      for(i in 1:ns)
		      {
			  res[i, ] <- c(density(Dred[i, -i], label=label[-i], K))
		      }
		      apilar[ , , j] <- res
		  } })

		  if (verbose) {print('Time after sixth chunk'); print(f)}

		  g <- system.time ({
		   aux <- t(plyr::aaply(apilar, c(2,1), mean))
		   auxx <- t(apply(apilar[, 1, ], 1, function(x){c(min(x), mean(x), max(x))}))
		   ps <- ns/(ns+nnull/fold)
		   p0 <- (nnull/fold)/(ns+nnull/fold)
		   segununiforme <- aux[, 1] - p0
		   genes <- (1:numGenes)[suspicious]
		   print(genes)
		   res <- cbind(genes, OR[suspicious], segununiforme, auxx, aux[, -1])
		   row.names(res) <- NULL
		   colnames(res) <- c("id", "OR", "DifExp",  "minFP", "FP", "maxFP", "dFP", "radius")
		   oo <- order(res[, 3], -res[, 2])
		   # print(res[oo,])
		   # print(ns)
		   # print(c(ps, p0))
		   # print(list("summary"=res[oo, ], "ns"=ns, "prop"=c(ps, p0)))
		   })

		   if (verbose) {print('Time after seventh chunk'); print(g)}

		   object@out <- list("summary"=res[oo, ], "ns"=ns, "prop"=c(ps, p0, K))
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

setMethod("initialize", "ORdensity", function(.Object, positive, negative, out, OR, FP, dFP, char, bestK, verbose = FALSE, parallel = FALSE, replicable = TRUE, seed = 0) {
  .Object@positive <- positive
  .Object@negative <- negative
  # validObject(.Object)
  .Object@verbose <- verbose
  .Object@parallel <- parallel
  .Object@replicable <- replicable
  .Object@seed <- seed
  .Object@out <- compute.ORdensity(.Object, verbose = .Object@verbose, parallel = .Object@parallel, replicable = .Object@replicable, seed = .Object@seed)
  .Object@OR <- .Object@out$summary[, "OR"]
  .Object@FP <- .Object@out$summary[, "FP"]
  .Object@dFP <- .Object@out$summary[, "dFP"]
  .Object@char <- data.frame(.Object@OR, .Object@FP, .Object@dFP)
  require(cluster)
  .Object@bestK <- findbestK(.Object)
  .Object
})

setGeneric("findbestK", function(object, ...) standardGeneric("findbestK"))

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
              DFgenes[[k]] <- list("cluster_number"=k, "genes"=sort(clusters[[k]][,'id']), "meanOR"=mean(clusters[[k]][,'OR']))
            }
            cat("The ORdensity method has found that the optimal clustering of the data consists of",object@bestK,"clusters\n")
            # return(list("DFgenes"=DFgenes,"clusters"=clusters))
            return(DFgenes)
          }
)
