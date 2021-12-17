# Plot g(theta) and the values of the statistic Tn = T_v^u/T_u^v


# theta is vector


reject.mfmm(log(3), c(1 / 9, 1, 9))



# -----


# plots mfmm
# organize the data


library(tidyverse)
library(reshape2)


mu <- log(3)
N <- 10000
x <- data.frame(
  theta = rep(c(1 / 9, 1, 9), 2), uv = rep(c(1, 2, 3), 2), v = c(rep(c(0.1, 0.2, 0.3), 2)),
  u = rep(c(0.2, 0.3, 0.4), 2), n = c(52, 79, 1101, rep(5000, 3))
)
x <- x %>% mutate(
  mu = mu,
  N = N,
  id = 1:n()
)
#create samples
samples <- mapply(qexp.samples, n = x$n, theta = x$theta, mu = x$mu, N = x$N)

get.mfmm.est <- list()
for (i in 1:nrow(x)) {
  get.mfmm.est[[i]] <- mfmm.samples(
    n = x$n[i], mu = x$mu[i], theta = x$theta[i],
    v = x$v[i], u = x$u[i], samples = samples[[i]]
  )
}


theta.hat.st <- lapply(get.mfmm.est, function(y) {
  y$theta.hat.pad
}) %>% melt()

colnames(theta.hat.st) <- c("value", "id")

q.hat.st <- lapply(get.mfmm.est, function(y) {
  y$q.hat.pad
}) %>% melt()
colnames(q.hat.st) <- c("value", "id")

theta.hat.st <- merge(x, theta.hat.st)
q.hat.st <- merge(x, q.hat.st)


theta.hat.st <- full_join(x, theta.hat.st, by = "id")
q.hat.st <- full_join(x, q.hat.st, by = "id")

colNames <- c("theta", "uv", "v", "u" ,
              "n",     "mu", "N", "id"  ,
              "theta.y", "uv.y", "v.y", "u.y" ,
              "n.y",     "mu.y", "N.y", "value")


colnames(theta.hat.st) <- colNames
colnames(q.hat.st) <- colNames

theta = c("1/9", "1", "9")
n <- c("5000", "1101", "79", "52")
plotCompDens.mfmm(data = theta.hat.st, n = n, theta = theta, parameter = "theta", standard = TRUE)
plotCompDens.mfmm(data = q.hat.st, n = n, theta = theta, parameter = "q", standard = TRUE)




