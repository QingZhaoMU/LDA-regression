#===============
# simulate data
#===============
# create independent x-values 
x1 <- rnorm(nl, 0, 1)
x2 <- rnorm(nl, 0, 1)
x <- cbind(x1, x2)
nx <- dim(x)[2]

eta <- alpha + x %*% beta
vmat <- matrix(, nrow=nl, ncol=nc)
for (i in 1:nl) {
  vmat[i,] <- inv.logit(rnorm(nc, eta[i,], exp(sigma)))
} # i
vmat[,nc] <- 1

theta <- vmat2theta(vmat)
#colMeans(theta)
#table(rowSums(theta))

par(mfrow=c(1,2))
par(mar=c(0,0,0,0))
plot(theta[,1] ~ x1, type='n', ylim=c(0,1), axes=F)
for (i in 1:nc) {
  points(theta[,i] ~ x1, col=i, pch=4)
  box()
}

plot(theta[,1] ~ x2, type='n', ylim=c(0,1), axes=F)
for (i in 1:nc) {
  points(theta[,i] ~ x2, col=i, pch=4)
  box()
}

phi <- matrix(rbeta(nc*ns, .01, .99), nc, ns)
for (i in 1:nc) {
  phi[i,] <- phi[i,] / sum(phi[i,])
}

pi <- theta %*% phi

y <- matrix(rbinom(nl*ns, nmat, pi), nl, ns)
skeep <- which(colMeans(y) > nl/100)
y <- y[,skeep]
ns <- dim(y)[2]
phi <- phi[,skeep]
phi.jump <- phi.jump[,skeep]



