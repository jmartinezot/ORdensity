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
#' @export
findDE <- setClass(
	"findDE",
	slots = c(positive="matrix", negative="matrix", out="list", OR="numeric", FP="numeric", dFP="numeric", char="data.frame", clustering2="integer", clustering4="integer")
)

setMethod("show",
	signature = "findDE",
	definition = function(object) {
		cat("Positive data matrix: ", object@positive, "\n", sep = "")
		cat("Negative data matrix: ", object@negative, "\n", sep = "")
		}
	)

setGeneric("summaries", function(object, ...) standardGeneric("summaries"))

setMethod("summaries",
  signature = "findDE",
  definition = function(object){
    list("summaryOR"=object@OR,
      "summarymeanFP"=object@FP,
      "summarydFP"=object@dFP)
    }
  )

setGeneric("plotFPvsOR", function(object, ...) standardGeneric("plotFPvsOR"))

setMethod("plotFPvsOR",
  signature = "findDE",
  definition = function(object){
    plot(object@FP, object@OR, type="n")
    points(object@FP, object@OR, cex=1/(0.5+object@dFP))
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
    plot(s, type="b", ylim=c(0,1))
  }
  )

setGeneric("clusplotk", function(object, k, ...) standardGeneric("clusplotk"))

setMethod("clusplotk",
  signature = "findDE",
  definition = function(object, k){
    aa <- pam(dist(scale(object@char)), k)
    clusplot(aa)
  }
  )
          
setGeneric("compute.findDE", function(object, ...) standardGeneric("compute.findDE"))

setMethod("compute.findDE", 
	signature = "findDE",
	definition =  function(object, B=100, scale=FALSE, alpha=0.05, fold=floor(B/10), weights=c(1/4,1/2,1/4), show_time=FALSE) {
		  a <- system.time ({
	    x <- as.matrix(object@positive)
		  y <- as.matrix(object@negative)
		  G <- dim(x)[1]
		  N1 <- dim(x)[2] 
		  N2 <- dim(y)[2]
		  N <- N1 + N2 })
		  
		  if (show_time) {print('a'); print(a)}
		  
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
		  
		  if (show_time) {print('b'); print(b)}
		  
		  c <- system.time ({
		  mv <- mv*matrix(rep(weights, G), byrow=TRUE, ncol=p)
		  Dxy <- dist(mv) #ICGE::dgower(mv, type=list(cuant=1:p))
		  })
		  
		  if (show_time) {print('c'); print(c)}
		  
		  d <- system.time ({
		  z <- cbind(x, y)
		  ORnull <- matrix(0, nrow=G, ncol=B)
		  mv.null <- array(0, dim=c(G, p, B))
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
		  } })
		  
		  if (show_time) {print('d'); print(d)}
		  
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
		  apilar <- array(0, dim=c(ns, 3, fold)) })
		  
		  if (show_time) {print('e'); print(e)}
		  
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
			  res[i, ] <- c(density(Dred[i, -i], label=label[-i], K=10))
		      }
		      apilar[ , , j] <- res
		  } })
		  
		  if (show_time) {print(f)}
		  
		  g <- system.time ({
		   aux <- t(plyr::aaply(apilar, c(2,1), mean))
		   auxx <- t(apply(apilar[, 1, ], 1, function(x){c(min(x), mean(x), max(x))}))
		   ps <- ns/(ns+nnull/fold)
		   p0 <- (nnull/fold)/(ns+nnull/fold)
		   segununiforme <- aux[, 1] - p0
		   genes <- (1:G)[suspicious]
		   res <- cbind(genes, OR[suspicious], segununiforme, auxx, aux[, -1])
		   row.names(res) <- NULL
		   colnames(res) <- c("id", "OR", "DifExp",  "minFP", "meanFP", "maxFP", "density", "radio")
		   oo <- order(res[, 3], -res[, 2])
		   object@out <- list("summary"=res[oo, ], "ns"=ns, "prop"=c(ps, p0))
		   # object
		  })
		  
		  if (show_time) {print('g'); print(g)}
		  
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

setMethod("initialize", "findDE", function(.Object, positive, negative, out, OR, FP, dFP, char, clustering2, clustering4) {
  .Object@positive <- positive
  .Object@negative <- negative
  # validObject(.Object)
  .Object@out <- compute.findDE(.Object)
  .Object@OR <- .Object@out$summary[, "OR"]
  .Object@FP <- .Object@out$summary[, "meanFP"]
  .Object@dFP <- .Object@out$summary[, "density"]
  .Object@char <- data.frame(.Object@OR, .Object@FP, .Object@dFP)
  require(cluster)
  .Object@clustering2 <- pam(dist(scale(.Object@char)), 2)$clustering
  .Object@clustering4 <- pam(dist(scale(.Object@char)), 4)$clustering
  .Object
})
