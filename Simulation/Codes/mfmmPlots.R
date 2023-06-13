# Plot g(theta) and the values of the statistic Tn = T_v^u/T_u^v

library(qExp)
# theta is vector
library(tidyr)
library(dplyr)

reject.mfmm(mu = log(3), theta = c(1 / 9, 1, 9))

# -----


# plots mfmm
# organize the data


library(tidyverse)
library(reshape2)


mu <- log(3)
N <- 10000
x <- data.frame(
  theta = rep(c(1/9, 1, 9), 2), uv = rep(c(1, 2, 3), 2), v = c(rep(c(0.07, 0.22, 0.49), 2)),
  u = rep(c(0.15, 0.33, 0.5), 2), n = c(48, 78, 817, rep(5000, 3))
)
x <- x  |>  mutate(
  mu = mu,
  N = N,
  id = 1:n()
)
#create samples
set.seed(2022)
samples <- mapply(qexp.samples, n = x$n, theta = x$theta, mu = x$mu, N = x$N)
get.mfmm.est <- list()
for (i in 1:nrow(x)) {
  get.mfmm.est[[i]] <- mfmmSamples(
    n = x$n[i], mu = x$mu[i], theta = x$theta[i],
    v = x$v[i], u = x$u[i], samples = samples[[i]]
  )
}


theta.hat.st <- lapply(get.mfmm.est, function(y) {
  y$theta.hat.pad
}) |> melt()

colnames(theta.hat.st) <- c("value", "id")

q.hat.st <- lapply(get.mfmm.est, function(y) {
  y$q.hat.pad
}) |> melt()
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
n <- c("5000", "817", "78", "48")
plotCompDens.mfmm(data = theta.hat.st, n = n, theta = theta, parameter = "theta", standard = TRUE)
plotCompDens.mfmm(data = q.hat.st, n = n, theta = theta, parameter = "q", standard = TRUE)



# get values no std

theta.hat <- lapply(get.mfmm.est, function(y) {
  y$theta.hat
}) |> melt()
colnames(theta.hat) <- c("value", "id")

q.hat <- lapply(get.mfmm.est, function(y) {
  y$q.hat
}) |> melt()
colnames(q.hat) <- c("value", "id")

theta.hat <- merge(x, theta.hat)
q.hat <- merge(x, q.hat)


theta.hat <- full_join(x, theta.hat, by = "id")
q.hat <- full_join(x, q.hat, by = "id")

colNames <- c("theta", "uv", "v", "u" ,
              "n",     "mu", "N", "id"  ,
              "theta.y", "uv.y", "v.y", "u.y" ,
              "n.y",     "mu.y", "N.y", "value")


colnames(theta.hat) <- colNames
colnames(q.hat) <- colNames

theta = c("1/9", "1", "9")
n <- c("5000", "817", "78", "48")
plotCompDens.mfmm(data = theta.hat, n = n, theta = theta, parameter = "theta", standard = FALSE)
plotCompDens.mfmm(data = q.hat, n = n, theta = theta, parameter = "q", standard = FALSE)






