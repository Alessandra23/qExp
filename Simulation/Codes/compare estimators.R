## Simulation study
library(qExp)
library(ggplot2)

N <- 10000
n = c(50, 500, 1000, 5000)
theta <- c(1/9, 1, 9)
mu <- c(1/10, log(3), 10)
q = round((3 + theta) / (2 + theta), 2)

allcomb <- expand.grid(N = N,
                       n = n,
                       mu = mu,
                       theta = theta)
ncomb <- nrow(allcomb)

# data <- list()
# for (i in 1:ncomb) {
#   set.seed(02)
#   comb <- allcomb[i,]
#   N <- comb$N
#   n <- comb$n
#   mu <- comb$mu
#   theta <- comb$theta
#
#   data[[i]] <- qexp.samples(N = N,n = n, theta = theta, mu = mu)
# }

# theta = 1/9, mu = log(3)

data1 <- list()
data1$N50 <- qexp.samples(N = N,n = n[1], theta = theta[1], mu = mu[2])
data1$N500 <- qexp.samples(N = N,n = n[2], theta = theta[1], mu = mu[2])
data1$N1000 <- qexp.samples(N = N,n = n[3], theta = theta[1], mu = mu[2])
data1$N5000 <- qexp.samples(N = N,n = n[4], theta = theta[1], mu = mu[2])

# filename <- paste('mu', round(mu[2],2),
#                   'theta', round(theta[1],2),
#                   sep='')
#
# data_filename <- paste('Simulation/Data/list', filename, '_data.RData', sep='')
# save(data, file = data_filename)


# MFMM
estMFMM1 <- list()
estMFMM1$N50 <- mfmmSamples(n = 50, mu = mu[2], theta = theta[1],
                           v = 0.07, u = 0.15, samples = data1$N50,
                           lower = 0, upper = 10000)$q.hat |> quantile(c(0.25, 0.5, 0.75))

estMFMM1$N500 <- mfmmSamples(n = 500, mu = mu[2], theta = theta[1],
                            v = 0.07, u = 0.15, samples = data1$N500,
                           lower = 0, upper = 10000)$q.hat |> quantile(c(0.25, 0.5, 0.75))

estMFMM1$N1000 <- mfmmSamples(n = 1000, mu = mu[2], theta = theta[1],
                            v = 0.07, u = 0.15, samples = data1$N1000,
                            lower = 0, upper = 10000)$q.hat |> quantile(c(0.25, 0.5, 0.75))

estMFMM1$N5000 <- mfmmSamples(n = 5000, mu = mu[2], theta = theta[1],
                            v = 0.07, u = 0.15, samples = data1$N5000,
                            lower = 0, upper = 10000)$q.hat |> quantile(c(0.25, 0.5, 0.75))
dfestMFMM1 <- plyr::ldply(estMFMM1)
dfestMFMM1
dfestMFMM1$Method <- rep('MFMM', 4)

# ML

estML1 <- list()

estML1$N50 <- mleSamples(data1$N50, mu[2], theta[1])
estML1$N50 <- estML1$N50$estMle.q |> quantile(c(0.25, 0.5, 0.75))

estML1$N500 <- mleSamples(data1$N500, mu[2], theta[1])
estML1$N500 <- estML1$N500$estMle.q |> quantile(c(0.25, 0.5, 0.75))

estML1$N1000 <- mleSamples(data1$N1000, mu[2], theta[1])
estML1$N1000 <- estML1$N1000$estMle.q |> quantile(c(0.25, 0.5, 0.75))

estML1$N5000 <- mleSamples(data1$N5000, mu[2], theta[1])
estML1$N5000 <- estML1$N5000$estMle.q |> quantile(c(0.25, 0.5, 0.75))

dfestML1 <- plyr::ldply(estML1)
dfestML1$Method <- rep('ML', 4)


# join the dfs

df1 <- rbind(dfestMFMM1, dfestML1)
df1

df1 |> ggplot(aes(x = factor(.id, level = c('N50', 'N500', 'N1000', 'N5000')), y = `50%`, linetype = Method)) +
  geom_point(position=position_dodge(0.2)) +
  geom_errorbar(aes(ymin=`25%`, ymax=`75%`), width=.2,
                position=position_dodge(0.2)) +
  theme_classic(base_size = 16) +
  labs(x = 'Sample size', y = expression(hat(q))) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  geom_hline(yintercept = q[1], color = 'gray', size = 0.4) +
  scale_x_discrete(labels = c('50', '500', '1000', '5000'))


# theta = 1, mu = log(3)

data2 <- list()
data2$N50 <- qexp.samples(N = N,n = n[1], theta = theta[2], mu = mu[2])
data2$N500 <- qexp.samples(N = N,n = n[2], theta = theta[2], mu = mu[2])
data2$N1000 <- qexp.samples(N = N,n = n[3], theta = theta[2], mu = mu[2])
data2$N5000 <- qexp.samples(N = N,n = n[4], theta = theta[2], mu = mu[2])

# filename <- paste('mu', round(mu[2],2),
#                   'theta', round(theta[1],2),
#                   sep='')
#
# data_filename <- paste('Simulation/Data/list', filename, '_data.RData', sep='')
# save(data, file = data_filename)


# MFMM
estMFMM2 <- list()
estMFMM2$N50 <- mfmmSamples(n = 50, mu = mu[2], theta = theta[2],
                           v = 0.22, u = 0.33, samples = data2$N50,
                           lower = 0, upper = 10000)$q.hat |> quantile(c(0.25, 0.5, 0.75))

estMFMM2$N500 <- mfmmSamples(n = 500, mu = mu[2], theta = theta[2],
                            v = 0.22, u = 0.33, samples = data2$N500,
                            lower = 0, upper = 10000)$q.hat |> quantile(c(0.25, 0.5, 0.75))

estMFMM2$N1000 <- mfmmSamples(n = 1000, mu = mu[2], theta = theta[2],
                             v = 0.22, u = 0.33, samples = data2$N1000,
                             lower = 0, upper = 10000)$q.hat |> quantile(c(0.25, 0.5, 0.75))

estMFMM2$N5000 <- mfmmSamples(n = 5000, mu = mu[2], theta = theta[2],
                             v = 0.22, u = 0.33, samples = data2$N5000,
                             lower = 0, upper = 10000)$q.hat |> quantile(c(0.25, 0.5, 0.75))
dfestMFMM2 <- plyr::ldply(estMFMM2)
dfestMFMM2
dfestMFMM2$Method <- rep('MFMM', 4)

# ML

estML2 <- list()

estML2$N50 <- mleSamples(data2$N50, mu[2], theta[2])
estML2$N50 <- estML2$N50$estMle.q |> quantile(c(0.25, 0.5, 0.75))

estML2$N500 <- mleSamples(data2$N500, mu[2], theta[2])
estML2$N500 <- estML2$N500$estMle.q |> quantile(c(0.25, 0.5, 0.75))

estML2$N1000 <- mleSamples(data2$N1000, mu[2], theta[2])
estML2$N1000 <- estML2$N1000$estMle.q |> quantile(c(0.25, 0.5, 0.75))

estML2$N5000 <- mleSamples(data2$N5000, mu[2], theta[2])
estML2$N5000 <- estML2$N5000$estMle.q |> quantile(c(0.25, 0.5, 0.75))

dfestML2 <- plyr::ldply(estML2)
dfestML2$Method <- rep('ML', 4)


# join the dfs

df2 <- rbind(dfestMFMM2, dfestML2)
df2


df2 |> ggplot(aes(x = factor(.id, level = c('N50', 'N500', 'N1000', 'N5000')), y = `50%`, linetype = Method)) +
  geom_point(position=position_dodge(0.2)) +
  geom_errorbar(aes(ymin=`25%`, ymax=`75%`), width=.2,
                position=position_dodge(0.2)) +
  theme_classic(base_size = 16) +
  labs(x = 'Sample size', y = expression(hat(q))) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  geom_hline(yintercept = q[2], color = 'gray', size = 0.4) +
  scale_x_discrete(labels = c('50', '500', '1000', '5000'))




# theta = 9, mu = log(3)

data3 <- list()
data3$N50 <- qexp.samples(N = N,n = n[1], theta = theta[3], mu = mu[2])
data3$N500 <- qexp.samples(N = N,n = n[2], theta = theta[3], mu = mu[2])
data3$N1000 <- qexp.samples(N = N,n = n[3], theta = theta[3], mu = mu[2])
data3$N5000 <- qexp.samples(N = N,n = n[4], theta = theta[3], mu = mu[2])

# filename <- paste('mu', round(mu[2],2),
#                   'theta', round(theta[1],2),
#                   sep='')
#
# data_filename <- paste('Simulation/Data/list', filename, '_data.RData', sep='')
# save(data, file = data_filename)


# MFMM
estMFMM3 <- list()
estMFMM3$N50 <- mfmmSamples(n = 50, mu = mu[2], theta = theta[3],
                           v = 0.49, u = 0.5, samples = data3$N50,
                           lower = 0, upper = 10000)$q.hat |> quantile(c(0.25, 0.5, 0.75))

estMFMM3$N500 <- mfmmSamples(n = 500, mu = mu[2], theta = theta[3],
                            v = 0.49, u = 0.5, samples = data3$N500,
                            lower = 0, upper = 10000)$q.hat |> quantile(c(0.25, 0.5, 0.75))

estMFMM3$N1000 <- mfmmSamples(n = 1000, mu = mu[2], theta = theta[3],
                             v = 0.49, u = 0.5, samples = data3$N1000,
                             lower = 0, upper = 10000)$q.hat |> quantile(c(0.25, 0.5, 0.75))

estMFMM3$N5000 <- mfmmSamples(n = 5000, mu = mu[2], theta = theta[3],
                             v = 0.49, u = 0.5, samples = data3$N5000,
                             lower = 0, upper = 10000)$q.hat |> quantile(c(0.25, 0.5, 0.75))
dfestMFMM3 <- plyr::ldply(estMFMM3)
dfestMFMM3
dfestMFMM3$Method <- rep('MFMM', 4)

# ML

estML3 <- list()

estML3$N50 <- mleSamples(data3$N50, mu[2], theta[3])
estML3$N50 <- estML3$N50$estMle.q |> quantile(c(0.25, 0.5, 0.75))

estML3$N500 <- mleSamples(data3$N500, mu[2], theta[3])
estML3$N500 <- estML3$N500$estMle.q |> quantile(c(0.25, 0.5, 0.75))

estML3$N1000 <- mleSamples(data3$N1000, mu[2], theta[3])
estML3$N1000 <- estML3$N1000$estMle.q |> quantile(c(0.25, 0.5, 0.75))

estML3$N5000 <- mleSamples(data3$N5000, mu[2], theta[3])
estML3$N5000 <- estML3$N5000$estMle.q |> quantile(c(0.25, 0.5, 0.75))

dfestML3 <- plyr::ldply(estML3)
dfestML3$Method <- rep('ML', 4)


# join the dfs

df3 <- rbind(dfestMFMM3, dfestML3)
df3


df3 |> ggplot(aes(x = factor(.id, level = c('N50', 'N500', 'N1000', 'N5000')), y = `50%`, linetype = Method)) +
  geom_point(position=position_dodge(0.2)) +
  geom_errorbar(aes(ymin=`25%`, ymax=`75%`), width=.2,
                position=position_dodge(0.2)) +
  theme_classic(base_size = 16) +
  labs(x = 'Sample size', y = expression(hat(q))) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  geom_hline(yintercept = q[3], color = 'gray', size = 0.4) +
  scale_x_discrete(labels = c('50', '500', '1000', '5000'))


## Make just one plot

df1$theta <- rep('q = 1.47', 8)
df2$theta <- rep('q = 1.33', 8)
df3$theta <- rep('q = 1.09', 8)

dfFinal <- rbind(df1, df2, df3)

hlineData <- data.frame(q = q, theta = c(df1$theta[1], df2$theta[1], df3$theta[1]))
geom_hline(aes(yintercept = z), hline.data)

dfFinal |> ggplot(aes(x = factor(.id, level = c('N50', 'N500', 'N1000', 'N5000')), y = `50%`, linetype = Method)) +
  geom_point(position=position_dodge(0.2)) +
  geom_errorbar(aes(ymin=`25%`, ymax=`75%`), width=.2,
                position=position_dodge(0.2)) +
  theme_classic(base_size = 20) +
  labs(x = 'Sample size', y = expression(hat(q))) +
  scale_linetype_manual(values = c('solid','dotted')) +
  geom_hline(aes(yintercept = q), hlineData, color = 'gray', size = 0.4) +
  scale_x_discrete(labels = c('50', '500', '1000', '5000')) +
  facet_wrap(~theta, scales = 'free_y')




# Compare the theoretical variances ---------------------------------------


vars <- function(n, mu){
  varMfmm <- vector()
  varMl <- vector()
  nn <- seq(1, 10, 0.1)
  for (i in 1:length(nn)) {
    varMfmm[i] <- mfmmTheo(mu = mu, theta = nn[i], v = 0.49, u = 0.5, n = n)$kappa2/n
    varMl[i] <- varTeoMle(mu = mu, theta = nn[i])$var.theta/n
  }

  dfVars <- data.frame(MFMM = varMfmm,
                       ML = varMl,
                       id = nn)

  return(dfVars)
}

varsdf <- vars(n = 50, mu = log(3))
varsdf <- varsdf |>  tidyr::pivot_longer(-id, names_to = "Method", values_to = "y")

varsdf |> ggplot(aes(x = id, y = y, linetype = Method))+
  geom_line() +
  theme_classic(base_size = 20) +
  labs(x = expression(hat(theta)), y = expression(paste('Var(', hat(theta), ')')))



varsmu <- function(n, theta){
  varMfmm <- vector()
  varMl <- vector()
  nnmu <- seq(1, 10, 0.1)

  for (i in 1:length(nnmu)) {
    varMfmm[i] <- (nnmu[i]^2 * (theta + 1)/(theta - 1))/n
    varMl[i] <- varTeoMle(mu = nnmu[i], theta = theta)$var.mu/n
  }

  dfVars <- data.frame(MFMM = varMfmm,
                       ML = varMl,
                       id = nnmu)

  return(dfVars)
}


varsmudf <- varsmu(n = 50, theta = 2)
varsmudf <- varsmudf |>  tidyr::pivot_longer(-id, names_to = "Method", values_to = "y")

varsmudf |> ggplot(aes(x = id, y = y, linetype = Method))+
  geom_line() +
  theme_classic(base_size = 20) +
  labs(x = expression(hat(mu)), y = expression(paste('Var(', hat(mu), ')')))



vardflist <- list(theta = varsdf,
                  mu = varsmudf) |>
            plyr::ldply()

parnames <- c(
  expression(hat(mu)),
  expression(hat(theta))
)

vardflist |> ggplot(aes(x = id, y = y, linetype = Method))+
  geom_line() +
  facet_wrap(~.id, scales = 'free_y',
             labeller = labeller(.id = as_labeller(parnames, label_parsed)))+
  theme_classic(base_size = 20) +
  labs(x = 'Parameter values', y = 'Variance')



# Wald test


testWald <- function(N, n, mu, theta, v, u, seed = 02){

  set.seed(seed = seed)
  data <- qexp.samples(N = N, n = n, theta = theta, mu = mu)

  estMfmm <- mfmmSamples(n = n, mu = mu, theta = theta,
                         v = v, u = u, samples = data,
                         lower = 0, upper = 10000)
  wMfmm <- (estMfmm$theta.hat - theta)/sqrt(estMfmm$kappa2/n)
  pvaluesWmfmm <- pchisq(wMfmm^2, df = 1, lower.tail = FALSE)
  #hist(pvaluesWmfmm, prob = TRUE)


  estMl <- mleSamples(data, mu, theta)
  varml <- varTeoMle(mu = mu, theta = theta)$var.theta/n
  wMl <- (estMl$estMle.theta - theta)/sqrt(varml)
  pvaluesWml <- pchisq(wMl^2, df = 1, lower.tail = FALSE)
  #hist(pvaluesWml, prob = TRUE)

  return(list(pvaluesWmfmm = pvaluesWmfmm,
              pvaluesWml = pvaluesWml))
}

# n = 5000
pvTheta1 <- testWald(N = N, n = n[4], mu = mu[2], theta = theta[1], v = 0.07, u = 0.15, seed = 02)
pvTheta2 <- testWald(N = N, n = n[4], mu = mu[2], theta = theta[2], v = 0.22, u = 0.33, seed = 02)
pvTheta3 <- testWald(N = N, n = n[4], mu = mu[2], theta = theta[3], v = 0.49, u = 0.5, seed = 02)

listdfs <- list(theta1 = pvTheta1,
                theta2 = pvTheta2,
                theta3 = pvTheta3)


dfpvalues <- reshape2::melt(listdfs)
dfpvalues$L2 <- as.factor(dfpvalues$L2)
dfpvalues$L1 <- as.factor(dfpvalues$L1)

levels(dfpvalues$L1) <- c(`0.1` = "theta[0]~`=`~1/9",
                          `1` = "theta[0]~`=`~1",
                          `9` = "theta[0]~`=`~9")

levels(dfpvalues$L2) <- c('MFMM', 'MLE')

dfpvalues |> ggplot(aes(x = value)) +
  geom_histogram(col = 'darkgray', fill = 'white', aes(y = after_stat(density)))+
  theme_classic(base_size = 20)  +
  facet_grid(L2 ~ L1,
             labeller = labeller(L1 = as_labeller(label_parsed)))+
  labs(x = 'p-values', y = 'Density')


# n = 50

pvTheta1 <- testWald(N = N, n = n[1], mu = mu[2], theta = theta[1], v = 0.07, u = 0.15, seed = 02)
pvTheta2 <- testWald(N = N, n = n[1], mu = mu[2], theta = theta[2], v = 0.22, u = 0.33, seed = 02)
pvTheta3 <- testWald(N = N, n = n[1], mu = mu[2], theta = theta[3], v = 0.49, u = 0.5, seed = 01)

listdfs <- list(theta1 = pvTheta1,
                theta2 = pvTheta2,
                theta3 = pvTheta3)


dfpvalues <- reshape2::melt(listdfs)
dfpvalues$L2 <- as.factor(dfpvalues$L2)
dfpvalues$L1 <- as.factor(dfpvalues$L1)

levels(dfpvalues$L1) <- c(`0.1` = "theta[0]~`=`~1/9",
                          `1` = "theta[0]~`=`~1",
                          `9` = "theta[0]~`=`~9")

levels(dfpvalues$L2) <- c('MFMM', 'MLE')

dfpvalues |> ggplot(aes(x = value)) +
  geom_histogram(col = 'darkgray', fill = 'white', aes(y = after_stat(density)))+
  theme_classic(base_size = 20)  +
  facet_grid(L2 ~ L1,
             labeller = labeller(L1 = as_labeller(label_parsed)))+
  labs(x = 'p-values', y = 'Density')



## Compare time

data = qexp.samples(N = 1, n = 1000, theta = theta[2], mu = mu[2]) |> as.vector()

mfmm <- function(dat){
  v <-  0.22
  u <- 0.33
  mu_uv <- (mean(dat^v))^u/(mean(dat^u))^v # T_v^u/T_u^v
  theta_hat_mf <- g.theta.inv(mu_uv, u = u, v = v, upper = 10000)
  mu_hat_mf <- mean(dat)
  return(theta_hat_mf)
}


mle <- function(dat){
  mu_in <- mean(dat)
  theta_in <- mean(dat^0.1)/(mean(dat)^0.1)
  tt <- optim(par = c(mu_in,theta_in), logQexp, x = dat, n = length(dat))
  return(tt$par)
}

times <- microbenchmark::microbenchmark("MFMM" = {mfmm(dat = data)},
                               "MLE" = {mle(dat = data)},
                               times = 1000,
                               unit = 'seconds')

times |> autoplot() + theme_classic()


#print(xtable::xtable(times, type = "latex"), file = "Simulation/Plots/times.tex")


