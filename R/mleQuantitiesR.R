#' Log of likelihood qexp
#' @param param vector of paramters \eqn{(\mu,\theta)}
#' @param n sample size
#' @export
logQexp <- function(param, x, n) {
  mu <- param[1]
  theta <- param[2]
  return(-(n * log((theta + 1) / (theta * mu)) - (theta + 2) * sum(log(1 + x / (mu * theta)))))
}

#' MLE estimation
#' @param sample Sample of size n
#' @export
estMle <- function(sample) {
  mu.est <- median(sample)
  theta.est <- mean(sample^0.1) / (mean(sample)^0.1)
  est.mu.theta <- optim(par = c(mu.est, theta.est), logQexp, n = length(sample), x = sample)
  return(est.mu.theta$par)
}

#' Expected Fisher Information
#' @param mu True value of mu
#' @param theta True value of theta
#' @export
KMuTheta <- function(mu, theta) {
  KF <- matrix(0, 2, 2)
  KF[1, 1] <- (theta+1)/(mu^2*(theta+3)) #((theta + 2)^2) / (mu^2 * (theta^2 + 5 * theta + 6))
  #KF[1, 1] <- ((theta + 2)^2) / (mu^2 * (theta^2 + 5 * theta + 6))
  KF[2, 2] <- (theta^2 + theta + 2) / (theta^2 * (theta + 1)^2 * (theta^2 + 5 * theta + 6))
  KF[1, 2] <- 2 / (mu * theta * (theta^2 + 5 * theta + 6))
  KF[2, 1] <- KF[1, 2]
  return(KF)
}


#' Theoretical variances
#' @param mu True value of mu
#' @param theta True value of theta
#' @export
varTeoMle <- function(mu, theta) {
  var.mu <- solve(KMuTheta(mu, theta))[1, 1]
  var.theta <- solve(KMuTheta(mu, theta))[2, 2]
  return(list(var.mu = var.mu, var.theta = var.theta))
}


#' Sample quantities MLE
#'
#' Theoretical quantities for the mle
#' @param samples matrix of samples Nxn
#' @param mu true value of mu
#' @param theta true value of theta
#' @export
mleSamples <- function(samples, mu, theta) {
  n <- ncol(samples)
  estMle <- t(apply(samples, 1, estMle))
  estMle.mu <- estMle[, 1]
  estMle.theta <- estMle[, 2]
  estMle.q <- (3 + estMle.theta) / (2 + estMle.theta)
  q.real <- (theta + 3) / (theta + 2)
  var.mu <- solve(KMuTheta(mu, theta))[1, 1] ## Of mu
  var.theta <- solve(KMuTheta(mu, theta))[2, 2] ## Of theta
  mu.pad.mle <- sqrt(n) * (estMle.mu - mu) / sqrt(var.mu)
  theta.pad.mle <- sqrt(n) * (estMle.theta - theta) / sqrt(var.theta)
  q.pad.mle <- sqrt(n) * (estMle.q - q.real) / sqrt(var.theta * ((theta + 2)^(-4)))

  return(list(
    estMle.mu = estMle.mu, estMle.theta = estMle.theta,
    estMle.q = estMle.q, mu.pad.mle = mu.pad.mle,
    theta.pad.mle = theta.pad.mle, q.pad.mle = q.pad.mle
  ))
}


