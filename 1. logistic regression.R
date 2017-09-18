rm(list=ls())
library(boot)
library(MASS)

path <- 'c:/A. UFL/1. simulation/'

#==============
# basic values
#==============
nl <- 200
nc <- 3
nc.max <- 5
ns <- 1000
nmat <- 100

niter <- 6000
nburn <- round(niter * .6)
adapt <- round(nburn * .6)

alpha <- logit(.85)
sigma <- log(1)
beta1 <- c(-1.4, .3, .2, -.2, -.1, .1)
beta2 <- c(.3, -.3, -.2, .2, .1, -.1)
beta <- rbind(beta1, beta2)[,1:nc]

theta.jump <- matrix(.1, nrow=nl, ncol=nc.max)
phi.jump <- matrix(.006, nrow=nc.max, ncol=ns)
sigma.jump <- .02
beta.jump <- .02

#==================
# define functions
#==================
source(paste(c(path, 'code/Functions.R'), collapse=''))

#===============
# simulate data
#===============
source(paste(c(path, 'code/Data Sims.R'), collapse=''))

#======
# test
#======
param <- list()

param$vmat <- cbind(vmat[,-nc], matrix(.5, nl, nc.max-nc+1))
param$vmat[,nc.max] <- 1
param$theta <- vmat2theta(param$vmat)
param$phi <- rbind(phi, matrix(1/ns, nc.max-nc, ns))

param$beta <- cbind(beta, matrix(0, nrow=nx, ncol=nc.max-nc))
param$sigma <- sigma

vmat.post <- theta.post <- theta.jump.post <- theta.accept <- array(, dim=c(nl, nc.max, niter))
vmat.post[,,1] <- param$vmat
theta.post[,,1] <- param$theta
theta.jump.post[,,1] <- theta.jump
theta.accept[,,1] <- FALSE

phi.post <- phi.jump.post <- phi.accept <- array(, dim=c(nc.max, ns, niter))
phi.post[,,1] <- param$phi
phi.jump.post[,,1] <- phi.jump
phi.accept[,,1] <- FALSE

beta.post <- beta.accept <- array(, dim=c(nx, nc.max, niter))
beta.post[,,1] <- param$beta
beta.jump.post <- numeric(niter)
beta.jump.post[1] <- beta.jump
beta.accept <- logical(niter)
beta.post.mat <- matrix(, nrow=nx*nc.max, ncol=niter)
beta.post.mat[,1] <- param$beta
beta.cor <- matrix(0, nx*nc.max, nx*nc.max)
diag(beta.cor) <- 1
beta.cor.post <- array(, dim=c(nx*nc.max, nx*nc.max, niter))
beta.cor.post[,,1] <- beta.cor

sigma.post <- sigma.jump.post <- numeric(niter)
sigma.post[1] <- param$sigma
sigma.jump.post[1] <- sigma.jump
sigma.accept <- logical(niter)

for (i in 2:niter) {
  theta.up <- update.theta(param, jump=theta.jump, nl, nc.max, y, x, nmat)
  vmat.post[,,i] <- param$vmat <- theta.up$vmat
  theta.post[,,i] <- param$theta <- theta.up$theta
  theta.accept[,,i] <- theta.up$accept

  phi.up <- update.phi(param, jump=phi.jump, nc.max, ns, y, nmat, a.phi=.01, b.phi=.99)
  phi.post[,,i] <- param$phi <- phi.up$phi
  phi.accept[,,i] <- phi.up$accept

  theta.jump.up <- thetaphi.jumpTune(accept=theta.accept, jump=theta.jump, ni=i, adapt=adapt, low=.3, high=.8)
  theta.jump.post[,,i] <- theta.jump <- theta.jump.up$jump

  phi.jump.up <- thetaphi.jumpTune(accept=phi.accept, jump=phi.jump, ni=i, adapt=adapt, low=.3, high=.8)
  phi.jump.post[,,i] <- phi.jump <- phi.jump.up$jump

  beta.up <- update.beta(param, jump=beta.jump, nl, nc.max, x)
  beta.post[,,i] <- param$beta <- beta.up$beta
  beta.accept[i] <- beta.up$accept
  beta.post.mat[,i] <- beta.post[,,i]
  if (i %% 50 == 0 & i > 50 & i < 2000) {
    beta.cor <- cor(t(beta.post.mat[,(i-49):i]))
  } else {
    beta.cor <- beta.cor
  }
  beta.cor[which(is.na(beta.cor))] <- 0
  beta.cor.post[,,i] <- beta.cor

  beta.jump.up <- betasigma.jumpTune(accept=beta.accept, jump=beta.jump, ni=i, adapt=adapt, low=.3, high=.8)
  beta.jump.post[i] <- beta.jump <- beta.jump.up$jump

  sigma.up <- update.sigma(param, jump=sigma.jump, nl, nc.max, x)
  sigma.post[i] <- param$sigma <- sigma.up$sigma
  sigma.accept[i] <- sigma.up$accept

  sigma.jump.up <- betasigma.jumpTune(accept=sigma.accept, jump=sigma.jump, ni=i, adapt=adapt, low=.3, high=.8)
  sigma.jump.post[i] <- sigma.jump <- sigma.jump.up$jump
}

theta.est <- apply(theta.post[,,nburn:niter], 1:2, median)
phi.est <- apply(phi.post[,,nburn:niter], 1:2, median)
beta.est <- apply(beta.post[,,nburn:niter], 1:2, median)

#===============
# graph results
#===============
par(mfrow=c(5,nc*2))
par(mar=c(4,4,1,1))

for (i in 1:nc) {
  plot(theta.est[,i] ~ theta[,i], xlim=c(0,1), ylim=c(0,1))
  abline(0, 1, col=2)
}

for (i in 1:nc) {
  plot(phi.est[i,] ~ phi[i,], xlim=c(0,1), ylim=c(0,1))
  abline(0, 1, col=2)
}

for (i in 1:nc) {
  theta.plot <- sample(1:nl, 10, replace=F)
  plot(theta.post[theta.plot[1],i,], type='l', ylim=range(theta.post[theta.plot,i,]))
  for (j in 1:length(theta.plot)) {
    lines(theta.post[theta.plot[j],i,], col=j)
  }
}

for (i in 1:nc) {
  phi.plot <- sample(1:ns, 10, replace=F)
  plot(phi.post[i,phi.plot[1],], type='l', ylim=range(phi.post[i,phi.plot,]))
  for (j in 1:length(phi.plot)) {
    lines(phi.post[i,phi.plot[j],], col=j)
  }
}

for (i in 1:nx) {
  for (j in 1:(nc-1)) {
    plot(beta.post[i,j,], type='l')
    abline(h=beta[i,j], col=2)
  }
}

plot(exp(sigma.post), type='l')
abline(h=exp(sigma), col=2)

boxplot(theta.est, 
        border=c(rep('red',nc),rep('blue',nc.max-nc)))

for (i in 1:nc) {
  plot(theta.jump.post[theta.plot[1],i,], type='l', ylim=range(theta.jump.post[theta.plot,i,]))
  for (j in 1:length(theta.plot)) {
    lines(theta.jump.post[theta.plot[j],i,], col=j)
  }
}

for (i in 1:nc) {
  plot(phi.jump.post[i,phi.plot[1],], type='l', ylim=range(phi.jump.post[i,phi.plot,]))
  for (j in 1:length(phi.plot)) {
    lines(phi.jump.post[i,phi.plot[j],], col=j)
  }
}

plot(beta.jump.post, type='l')
legend('topright', legend='', bty='n', 
       title=round(sum(beta.accept[nburn:niter]) / length(beta.accept[nburn:niter]), digits=3))

plot(sigma.jump.post, type='l')
legend('topright', legend='', bty='n', 
       title=round(sum(sigma.accept[nburn:niter]) / length(sigma.accept[nburn:niter]), digits=3))

theta.accept.rate <- apply(theta.accept, 2:3, sum) / nl
plot(theta.accept.rate[1,], type='n', ylim=range(theta.accept.rate[1:nc,]))
for (i in 1:nc) {
  lines(theta.accept.rate[i,], col=i)
}

phi.accept.rate <- apply(phi.accept, c(1,3), sum) / nl
plot(phi.accept.rate[1,], type='n', ylim=range(phi.accept.rate[1:nc,]))
for (i in 1:nc) {
  lines(phi.accept.rate[i,], col=i)
}


