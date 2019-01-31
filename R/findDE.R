#' @title
#' Class for representing findDE
#'
#' @description
#' Class for representing findDE
#'
#' @name findDE
#' @exportClass findDE
#'
#' @author Itziar Irigoien, Concepcion Arenas, Jose Maria Martinez-Otzeta
#'
#' @export plotFPvsOR
#' @export plotclusters
#' @export clusplotk
#' @export compute.findDE
#' @export findDEgenes
#' @export summary
#' @export preclusteredData

findDE <- setClass(
	"findDE",
	slots = c(positive="matrix", negative="matrix", out="list", OR="numeric", FP="numeric", dFP="numeric", char="data.frame", bestK = "numeric", verbose="logical", parallel="logical", replicable="logical", seed="numeric"),
	prototype = list(positive=matrix(), negative=matrix(), out=list(), OR=numeric(), FP=numeric(), dFP=numeric(), char=data.frame(), bestK = numeric(), verbose=logical(), parallel=logical(), replicable=logical(), seed=numeric())
)

# setGeneric("summary")

# setGeneric("summary.findDE", function(object, ...) standardGeneric("summary.findDE"))

setMethod("summary",
          signature = "findDE",
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
              cat("The expected number of false positives neighbours is", expectedFalsePositiveNeighbours, "\n")
              cat("The ORdensity method has found that the optimal clustering of the data consists of",object@bestK,"clusters\n\n")
              return(list("neighbours"=neighbours, "expectedFalsePositiveNeighbours"=p0*neighbours, "clusters"=clusters))
        }
)

setGeneric("preclusteredData", function(object, ...) standardGeneric("preclusteredData"))

setMethod("preclusteredData",
          signature = "findDE",
          definition = function(object){
              prop <- object@out$prop
              neighbours <- prop[3]
              p0 <- prop[2]
              preclustered_data <- as.data.frame(object@out$summary)
              preclustered_data$DifExp <- NULL
              preclustered_data$minFP <- NULL
              preclustered_data$maxFP <- NULL
              preclustered_data$radius <- NULL
              preclustered_data$S <- ifelse(preclustered_data$FP == 0, "S", "")
              preclustered_data$F <- ifelse(preclustered_data$FP < p0 * neighbours, "F", "")
              cat("Column S denotes the cases when FP=0\n")
              cat("Column F denotes the cases when FP < expectedFalsePositives")
              preclustered_data
          }
)

#setMethod("show",
#          signature = "findDE",
#          definition = function(object) {
#            cat("Positive data matrix: ", object@positive, "\n", sep = " ")
#            cat("Negative data matrix: ", object@negative, "\n", sep = " ")
#          }
#)

setGeneric("plotFPvsOR", function(object, ...) standardGeneric("plotFPvsOR"))

setMethod("plotFPvsOR",
  signature = "findDE",
  definition = function(object, k = object@bestK){
    clustering <- pam(dist(scale(object@char)), k)$clustering
    legend_text <- sprintf("cluster %s",seq(1:k))
    plot(object@FP, object@OR, type="n",main="Potential genes",xlab="FP",ylab="OR")
    points(object@FP, object@OR, cex=1/(0.5+object@dFP),  col = clustering)
    legend("topright", legend=legend_text, pch=16, col=unique(clustering))
    }
  )

setGeneric("plotclusters", function(object, ...) standardGeneric("plotclusters"))

setMethod("plotclusters",
  signature = "findDE",
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

setGeneric("clusplotk", function(object, k, ...) standardGeneric("clusplotk"))

setMethod("clusplotk",
  signature = "findDE",
  definition = function(object, k){
    aa <- pam(dist(scale(object@char)), k)
    clusplot(aa, main = paste("Clustering with k = ", k))
  }
  )


setGeneric("bootstrap", function(object, N, N1, scale, weights, G, p, ...) standardGeneric("bootstrap"))

setMethod("bootstrap",
  signature = "findDE",
  definition = function(object, N, N1, scale, weights, G, p){
	s1 <- sample(1:N, N1, replace=FALSE)
	s2 <- (1:N)[-s1]
	aux1 <- z[, s1]
	aux2 <- z[, s2]


	qx.b <- t(apply(aux1, 1, quantile, probs=c(0.25, 0.5, 0.75)))
	#    mx.b <- apply(aux1, 1, mean)
	qy.b <- t(apply(aux2, 1, quantile, probs=c(0.25, 0.5, 0.75)))
	#   my.b <- apply(aux2, 1, mean)
	mz <- cbind(qx.b - qy.b) #, mx.b-my.b)

	if(scale)
	{
	RIx.b <- qx.b[,3]-qx.b[,1]
	RIy.b <- qy.b[,3]-qy.b[,1]
	maxRI.b <- apply(cbind(RIx.b, RIy.b), 1, max)
	mz <- mz/maxRI.b
	}
	mv.null[ , ,b] <- mz*matrix(rep(weights, G), byrow=TRUE, ncol=p)

	#D0 <- as.matrix(ICGE::dgower(mz, type=list(cuant=1:p)))
	D0 <- as.matrix(dist(mv.null[,,b]))
	vgR0 <- median(D0^2)/2
	I <- apply(D0, 1,  IindexRobust, vg=vgR0)
	ORnull[, b] <- 1/I
  }
  )

setGeneric("compute.findDE", function(object, ...) standardGeneric("compute.findDE"))

setMethod("compute.findDE",
	signature = "findDE",
	definition =  function(object, B=100, scale=FALSE, alpha=0.05, fold=floor(B/10), weights=c(1/4,1/2,1/4), K = 10, verbose=FALSE, parallel = FALSE, replicable = TRUE, seed = 0) {
	    a <- system.time ({
	    x <- as.matrix(object@positive)
		  y <- as.matrix(object@negative)
		  G <- dim(x)[1]
		  N1 <- dim(x)[2]
		  N2 <- dim(y)[2]
		  N <- N1 + N2 })

		  gloN <<- N
		  gloN1 <<- N1

		  if (verbose) {print('a'); print(a)}

		  b <- system.time ({
		  qx <- t(apply(x, 1, quantile, probs=c(0.25, 0.5, 0.75)))
		 # mx <- apply(x, 1, mean)
		  qy <- t(apply(y, 1, quantile, probs=c(0.25, 0.5, 0.75)))
		 # my <- apply(y, 1, mean)
		  mv <- cbind(qx - qy) #, mx-my)

		  p <- dim(mv)[2]

		  if(scale){
		    RIx <- qx[,3]-qx[,1]
		    RIy <- qy[,3]-qy[,1]
		    maxRI <- apply(cbind(RIx, RIy), 1, max)
		    if(any(maxRI==0))
		    {stop('Can\'t scale the data')}
		    mv <- mv/maxRI
		  } })

		  if (verbose) {print('b'); print(b)}

		  c <- system.time ({
		  mv <- mv*matrix(rep(weights, G), byrow=TRUE, ncol=p)
		  Dxy <- dist(mv) #ICGE::dgower(mv, type=list(cuant=1:p))
		  })

		  if (verbose) {print('c'); print(c)}

		  d <- system.time ({
		  z <- cbind(x, y)
		  ORnull <- matrix(0, nrow=G, ncol=B)
		  mv.null <- array(0, dim=c(G, p, B))

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
  			s1 <- sample(1:N, N1, replace=FALSE)
  			s2 <- (1:N)[-s1]
  			aux1 <- z[, s1]
  			aux2 <- z[, s2]

  			qx.b <- t(apply(aux1, 1, quantile, probs=c(0.25, 0.5, 0.75))) # se puede tener antes
  			#    mx.b <- apply(aux1, 1, mean)
  			qy.b <- t(apply(aux2, 1, quantile, probs=c(0.25, 0.5, 0.75))) # se puede tener antes
  			#   my.b <- apply(aux2, 1, mean)
  			mz <- cbind(qx.b - qy.b) #, mx.b-my.b)

  			if(scale)
  			{
  			RIx.b <- qx.b[,3]-qx.b[,1] # se puede tener antes
  			RIy.b <- qy.b[,3]-qy.b[,1] # se puede tener antes
  			maxRI.b <- apply(cbind(RIx.b, RIy.b), 1, max)
  			mz <- mz/maxRI.b
  			}
  			res <- list()
  			# mv.null[ , ,b] <- mz*matrix(rep(weights, G), byrow=TRUE, ncol=p)
  			res[[1]] <- mz*matrix(rep(weights, G), byrow=TRUE, ncol=p)

  			#D0 <- as.matrix(ICGE::dgower(mz, type=list(cuant=1:p)))
  			# D0 <- as.matrix(dist(mv.null[,,b]))
  			D0 <- as.matrix(dist(res[[1]]))
  			vgR0 <- median(D0^2)/2
  			I <- apply(D0, 1,  IindexRobust, vg=vgR0)
  			# ORnull[, b] <- 1/I
  			res[[2]] <- 1/I
  			res
		  }
      parallel::stopCluster(cl)

      glo <<- res_par

      for (b in 1:B) {
        mv.null[ , ,b] <- res_par[[b*2-1]]
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
		    s1 <- sample(1:N, N1, replace=FALSE)
		    s2 <- (1:N)[-s1]
		    aux1 <- z[, s1]
		    aux2 <- z[, s2]

		    qx.b <- t(apply(aux1, 1, quantile, probs=c(0.25, 0.5, 0.75)))
		    #    mx.b <- apply(aux1, 1, mean)
		    qy.b <- t(apply(aux2, 1, quantile, probs=c(0.25, 0.5, 0.75)))
		    #   my.b <- apply(aux2, 1, mean)
		    mz <- cbind(qx.b - qy.b) #, mx.b-my.b)

		    if(scale)
		    {
		      RIx.b <- qx.b[,3]-qx.b[,1]
		      RIy.b <- qy.b[,3]-qy.b[,1]
		      maxRI.b <- apply(cbind(RIx.b, RIy.b), 1, max)
		      mz <- mz/maxRI.b
		    }

		    mv.null[ , ,b] <- mz*matrix(rep(weights, G), byrow=TRUE, ncol=p)


		    #D0 <- as.matrix(ICGE::dgower(mz, type=list(cuant=1:p)))
		    D0 <- as.matrix(dist(mv.null[,,b]))
		    vgR0 <- median(D0^2)/2
		    I <- apply(D0, 1,  IindexRobust, vg=vgR0)
		    ORnull[, b] <- 1/I
		  }
		}
		  })

		  if (verbose) {print('d'); print(d)}

		# OR values for original data
		  e <- system.time ({
		  Dxy <- as.matrix(Dxy)
		  vgR <- median(Dxy^2)/2
		  I <- apply(Dxy, 1,  IindexRobust, vg=vgR)
		  OR <- 1/I

		# Find cut point
		#  M <- apply(ORnull, 2, sort)
		#  ORnullord <- apply(M, 1, median)
		#  corte <- ORnullord[floor((1-alpha)*G)]
		   corte <- (sort(c(ORnull)))[floor((1-alpha)*G*B)]

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

		  if (verbose) {print('e'); print(e)}

		# Que pasa en la capa 1
		  f <- system.time ({
		  for (j in 1:fold)
		  {
		      capa <- mirar == j
		      nnull.red <- sum(capa)
		      dat <- matrix(0, nrow=nnull.red, ncol=p)
		      cont <- 1
		      for (i in (1:nnull)[capa])
		      {
			   f <- selec.null[i, 1]
			   b <- selec.null[i, 2]
			   dat[cont, ] <- mv.null[f, , b]
			   cont <- cont + 1
		       }
		      dat <- rbind(mv[suspicious, ], dat)
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

		  if (verbose) {print(f)}

		  g <- system.time ({
		   aux <- t(plyr::aaply(apilar, c(2,1), mean))
		   auxx <- t(apply(apilar[, 1, ], 1, function(x){c(min(x), mean(x), max(x))}))
		   ps <- ns/(ns+nnull/fold)
		   p0 <- (nnull/fold)/(ns+nnull/fold)
		   segununiforme <- aux[, 1] - p0
		   genes <- (1:G)[suspicious]
		   res <- cbind(genes, OR[suspicious], segununiforme, auxx, aux[, -1])
		   row.names(res) <- NULL
		   colnames(res) <- c("id", "OR", "DifExp",  "minFP", "FP", "maxFP", "dFP", "radius")
		   oo <- order(res[, 3], -res[, 2])
		   # print(res[oo,])
		   # print(ns)
		   # print(c(ps, p0))
		   # print(list("summary"=res[oo, ], "ns"=ns, "prop"=c(ps, p0)))
		   })

		   if (verbose) {print('g'); print(g)}

		   object@out <- list("summary"=res[oo, ], "ns"=ns, "prop"=c(ps, p0, K))
		   # object
		}
	)

# a <- findDE(positive = matrix(rnorm(500), nrow=50, ncol=10), negative = matrix(rnorm(500), nrow=50, ncol=10))
# a <- compute.findDE(a)

setValidity("findDE", function(object) {
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

setMethod("initialize", "findDE", function(.Object, positive, negative, out, OR, FP, dFP, char, bestK, verbose = FALSE, parallel = FALSE, replicable = TRUE, seed = 0) {
  .Object@positive <- positive
  .Object@negative <- negative
  # validObject(.Object)
  .Object@verbose <- verbose
  .Object@parallel <- parallel
  .Object@replicable <- replicable
  .Object@seed <- seed
  .Object@out <- compute.findDE(.Object, verbose = .Object@verbose, parallel = .Object@parallel, replicable = .Object@replicable, seed = .Object@seed)
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
          signature = "findDE",
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
          signature = "findDE",
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
