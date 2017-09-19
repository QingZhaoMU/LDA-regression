#==================
# define functions
#==================
vmat2theta <- function(vmat) {
  nl <- dim(vmat)[1]
  nc <- dim(vmat)[2]
  theta <- vmat
  for (i in 1:nc) {
    if (i==1) {
      prod <- rep(1,nl)
    } else {
      prod <- prod * (1-vmat[,i-1])
    }
    theta[,i] <- vmat[,i] * prod
  } # i
  return(theta)
} # vmat2theta

tnorm <- function(n, lo, hi, mu, sig) {
  # generates truncated normal variates based on cumulative normal distribution normal truncated lo and hi

  if (length(lo) == 1 & length(mu) > 1) 
    lo <- rep(lo, length(mu))
  if (length(hi) == 1 & length(mu) > 1) 
    hi <- rep(hi, length(mu))
    
  q1 <- pnorm(lo, mu, sig)  #cumulative distribution
  q2 <- pnorm(hi, mu, sig)  #cumulative distribution
    
  z <- runif(n, q1, q2)
  z <- qnorm(z, mu, sig)
  z[z == -Inf] <- lo[z == -Inf]
  z[z == Inf] <- hi[z == Inf]
  z
} # tnorm

fix.MH <- function(lo, hi, old1, new1, jump) {
  jold <- pnorm(hi, mean = old1, sd = jump) - pnorm(lo, mean = old1, sd = jump)
  jnew <- pnorm(hi, mean = new1, sd = jump) - pnorm(lo, mean = new1, sd = jump)
  log(jold) - log(jnew)  #add this to pnew
} # fix.MH

get.logl <- function(theta, phi, y, nmat) {
  prob <- theta %*% phi
  cond <- prob < 1e-05
  prob[cond] <- 1e-05
  cond <- prob > 0.99999
  prob[cond] <- 0.99999
  dbinom(y, size = nmat, prob = prob, log = T)
} # get.logl

acceptMH <- function(p0, p1, x0, x1, BLOCK) {
  # accept for M, M-H if BLOCK, then accept as a block, otherwise, accept individually
    
  nz <- length(x0)  #no. to accept
  if (BLOCK) 
    nz <- 1
    
  a <- exp(p1 - p0)  #acceptance PR
  z <- runif(nz, 0, 1)
  keep <- which(z < a)
    
  if (BLOCK & length(keep) > 0) 
    x0 <- x1
  if (!BLOCK) 
    x0[keep] <- x1[keep]
  accept <- length(keep)
    
  list(x = x0, accept = accept)
} # acceptMH

update.theta <- function(param, jump, nl, nc, y, x, nmat) {
  phi <- param$phi
  beta <- param$beta
  sigma <- param$sigma

  eta <- alpha + x %*% beta

  vmat.ori <- vmat.old <- param$vmat
  vmat.tmp <- tnorm(nl*(nc-1), lo=0, hi=1, mu=vmat.old[,-nc], sig=jump[,-nc])
  proposed <- cbind(matrix(vmat.tmp, nl, nc-1), 1)
  adj <- fix.MH(lo=0, hi=1, old1=vmat.old, new1=proposed, jump=jump)

  for (j in 1:(nc-1)) {
    # last column has to be 1
    vmat.new <- vmat.old
    vmat.new[,j] <- proposed[,j]

    theta.old <- vmat2theta(vmat=vmat.old)
    theta.new <- vmat2theta(vmat=vmat.new)

    prob.old <- get.logl(theta=theta.old, phi=phi, y=y, nmat=nmat)
    prob.new <- get.logl(theta=theta.new, phi=phi, y=y, nmat=nmat)

    pold <- rowSums(prob.old) + dnorm(logit(vmat.old[,j]), eta[,j], exp(sigma), log = T)
    pnew <- rowSums(prob.new) + dnorm(logit(vmat.new[,j]), eta[,j], exp(sigma), log = T)

    k <- acceptMH(p0=pold, p1=pnew+adj[,j], x0=vmat.old[,j], x1=vmat.new[,j], BLOCK=F)
    vmat.old[,j] <- k$x
  }
  vmat <- vmat.old
  theta <- vmat2theta(vmat=vmat)
  list(theta=theta, vmat=vmat, accept=vmat.ori!=vmat.old)
} # update.theta

update.phi <- function(param, jump, nc, ns, y, nmat, a.phi, b.phi) {
  theta <- param$theta

  phi.ori <- phi.old <- param$phi
  proposed <- matrix(tnorm(nc*ns, lo=0, hi=1, mu=phi.old, sig=jump), nc, ns)
  adj <- fix.MH(lo=0, hi=1, old1=phi.old, new1=proposed, jump=jump)

  for (j in 1:nc) {
    phi.new <- phi.old
    phi.new[j,] <- proposed[j,]

    prob.old <- get.logl(theta=theta, phi=phi.old, y=y, nmat=nmat)
    prob.new <- get.logl(theta=theta, phi=phi.new, y=y, nmat=nmat)

    pold <- colSums(prob.old) + dbeta(phi.old[j,], a.phi, b.phi, log=T)
    pnew <- colSums(prob.new) + dbeta(phi.new[j,], a.phi, b.phi, log=T)

    k <- acceptMH(p0=pold, p1=pnew + adj[j,], x0=phi.old[j,], x1=phi.new[j,], BLOCK=F)
    phi.old[j,] <- k$x
  }
  phi <- phi.old
  list(phi=phi, accept=phi.ori!=phi.old)
} # update.phi

update.beta <- function(param, jump, nl, nc, x) {
  vmat <- param$vmat
  sigma <- param$sigma

  beta.old <- param$beta
  beta.new <- matrix(mvrnorm(n=1, mu=as.vector(beta.old), Sigma=jump*beta.cor), 
                     nrow=nx, ncol=nc.max, byrow=F)

  eta.old <- alpha + x %*% beta.old
  eta.new <- alpha + x %*% beta.new

  prob.old <- dnorm(logit(vmat[,-nc]), eta.old[,-nc], exp(sigma), log=T)
  prob.new <- dnorm(logit(vmat[,-nc]), eta.new[,-nc], exp(sigma), log=T)

  pold <- sum(prob.old) + sum(dnorm(beta.old[,-nc], 0, 5, log=T))
  pnew <- sum(prob.new) + sum(dnorm(beta.new[,-nc], 0, 5, log=T))

  k <- acceptMH(p0=pold, p1=pnew, x0=beta.old, x1=beta.new, BLOCK=F)
  beta <- k$x
  list(beta=beta, accept=beta[1,1]!=beta.old[1,1])
} # update.beta

update.sigma <- function(param, jump, nl, nc, x) {
  vmat <- param$vmat
  beta <- param$beta
  eta <- alpha + x %*% beta

  sigma.old <- param$sigma
  sigma.new <- rnorm(1, sigma.old, jump)

  prob.old <- dnorm(logit(vmat[,-nc]), eta[,-nc], exp(sigma.old), log=T)
  prob.new <- dnorm(logit(vmat[,-nc]), eta[,-nc], exp(sigma.new), log=T)

  pold <- sum(prob.old) + dnorm(sigma.old, 0, 5, log=T)
  pnew <- sum(prob.new) + dnorm(sigma.new, 0, 5, log=T)

  k <- acceptMH(p0=pold, p1=pnew, x0=sigma.old, x1=sigma.new, BLOCK=F)
  sigma <- k$x
  list(sigma=sigma, accept=sigma!=sigma.old)
} # update.sigma

thetaphi.jumpTune <- function(accept, jump, ni, adapt=2000, low=.3, high=.8) {
  nstart <- ifelse(ni>=100, ni-99, 1)
  if (ni > adapt) {
    jump.new <- jump
  } else if (ni %% 10 > 0) {
    jump.new <- jump
  } else {
    accept.rate <- apply(accept[,,nstart:ni], 1:2, mean)
    jump.new <- jump
    jump.new[which(accept.rate < low)] <- 
      jump[which(accept.rate < low)] * rnorm(length(which(accept.rate < low)),.5,.01)
    jump.new[which(accept.rate > high)] <- 
      jump[which(accept.rate > high)] * rnorm(length(which(accept.rate > high)),2,.01)
  }
  jump <- jump.new
  list(jump=jump)
} # jumpTune

betasigma.jumpTune <- function(accept, jump, ni, adapt=2000, low=.3, high=.8) {
  nstart <- ifelse(ni>=100, ni-99, 1)
  if (ni > adapt) {
    jump.new <- jump
  } else if (ni %% 10 > 0) {
    jump.new <- jump
  } else {
    accept.rate <- mean(accept[nstart:ni], na.rm=T)
    if (accept.rate < low) {
      jump.new <- jump * rnorm(1,.5,.05)
    } else if (accept.rate > high) {
      jump.new <- jump * rnorm(1,2,.05)
    } else {
      jump.new <- jump
    }
  }
  jump <- jump.new
  list(jump=jump)
} # sigma.jumpTune
