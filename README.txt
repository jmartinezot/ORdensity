Process:

Remove the package from the current installation using

sudo R

remove.packages('findDE')

install again doing

sudo R

library('devtools')
install_github('jmartinezot/findDE')
library('findDE')

installing from local source

install.packages("/home/otzeta/Github/findDE", repos = NULL, type = "source")
library('findDE')

remove all variables in environment

rm(list=ls())

initialize an findDE object
x <- example[, 3:32]
y <- example[, 33:62]
positive <- as.matrix(x)
negative <- as.matrix(y)
randomA <- matrix(rnorm(30*100), nrow = 100)
randomB <- matrix(rnorm(30*100), nrow = 100)
myfindDE <- new("findDE", positive = positive, negative = negative)
myfindDE <- new("findDE", positive = randomA, negative = randomB)

let us see what it contains
myfindDE@positive
myfindDE@negative
myfindDE@.out # this has to be empty

initialize the .out field
myfindDE <- compute.findDE(myfindDE)
out <- myfindDE@.out

out$summary
OR <- out$summary[, "OR"]
FP <- out$summary[, "mediaFP"] # Aldatu izena
dFP <- out$summary[, "density"]
plot(FP, OR, type="n")
points(FP, OR, cex=1/(0.5+dFP))

library(cluster)
char <- data.frame(OR, FP, dFP)
s <- rep(NA, 10)
for (k in 2:10)
{
  aux <- pam(dist(scale(char)), k)
  s[k] <- mean(silhouette(aux)[, "sil_width"])
}
plot(s, type="b", ylim=c(0,1))
# Try 2, 3, and 4 clusters
k <- 4
clustering <- pam(dist(scale(char)), k)$clustering
by(char, clustering, summary)

veryimportant <- out$summary[clustering==1, "id"]
table(dat$DEgen[veryimportant])
table( dat$gap[veryimportant])
important <- out$summary[clustering==2, "id"]
table(dat$DEgen[important])
table( dat$gap[important])
doubtfull <- out$summary[clustering==3, "id"]
table(dat$DEgen[doubtfull])
table( dat$gap[doubtfull])
nonimportant <- out$summary[clustering==4, "id"]
table(dat$DEgen[nonimportant])

clustering <- pam(dist(scale(char)), 2)$clustering
veryimportant <- out$summary[clustering==1, "id"]
table(dat$DEgen[veryimportant])
table( dat$gap[veryimportant])
nonimportant <- out$summary[clustering==2, "id"]
table(dat$DEgen[nonimportant])


aa <- pam(dist(scale(char)), 2)
clusplot(aa)
aa <- pam(dist(scale(char)), 4)
ama <- clusplot(aa)
aux <- prcomp(char, scale=TRUE)
x <- aux$x[, 1:2]
getAnywhere(clusplot.default)
