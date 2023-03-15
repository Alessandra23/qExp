load("~/Documents/GitHub/qExp/Simulation/Data/n50mu1.1theta0.11_data.RData")
data

N = 10000
n <- c(20, 30, 50, 100, 500,1000)
mu <- c(1/10, log(3), 10)
theta <- c(1/9, 1, 9)

allcomb <- expand.grid(N = N,
                       n = n,
                       mu = mu,
                       theta = theta)


ncomb <- nrow(allcomb)
model <- list()

for (i in 1:ncomb) {
  comb <- allcomb[i,]
  N <- comb$N
  n <- comb$n
  mu <- comb$mu
  theta <- comb$theta

  filename <- paste('n', n,
                    'mu', round(mu,2),
                    'theta', round(theta,2),
                    sep='')

  data_filename <- paste('~/Documents/GitHub/qExp/Simulation/Data/', filename, '_data.RData', sep='')
  load(data_filename)

  model[[i]] <- mle.samples(data, mu, theta)
}

# saveRDS(model, '~/Documents/GitHub/qExp/Simulation/Data/model.rds')

model <- readRDS("~/Documents/GitHub/qExp/Simulation/Data/model.rds")


rows <- which(allcomb$n %in% c(100,500,1000))
dfrows <- allcomb[rows,]
mudf <- list()
thetadf <- list()
for (i in 1:length(rows)) {
  n = dfrows[i, "n"]
  mu = dfrows[i, "mu"]
  theta = dfrows[i, "theta"]
  estmu = model[[rows[i]]]$mu.pad.mle
  esttheta = model[[rows[i]]]$theta.pad.mle

  mudf[[i]] <- data.frame(n = rep(n, 10000), mu = rep(mu, 10000),
                          theta = rep(theta, 10000), est = estmu)
  thetadf[[i]] <- data.frame(n = rep(n, 10000), mu = rep(mu, 10000),
                             theta = rep(theta, 10000), est = esttheta)
}

mudfN <- do.call(rbind,mudf)
colnames(mudfN) <- c("n", "mu" ,  "theta",  "muest")
thetadfN <- do.call(rbind,thetadf)
colnames(thetadfN) <- c("n", "mu" ,  "theta",  "thetaest")


dfF <- cbind(mudfN, thetadfN$thetaest)
colnames(dfF) <- c("n", "mu" ,  "theta", "muest" ,"thetaest")
dfF[dfF$mu == 0.1,"mu"] <- paste0('\u03BC'," = 1/10")
dfF[dfF$mu == log(3),"mu"] <- paste0('\u03BC'," = log(3)")
dfF[dfF$mu == 10,"mu"] <- paste0('\u03BC'," = 10")
dfF$mu <- factor(dfF$mu, levels = unique(dfF$mu))

dfF[dfF$theta == 1/9,"theta"] <- paste0('\u03B8'," = 1/9")
dfF[dfF$theta == 1,"theta"] <- paste0('\u03B8'," = 1")
dfF[dfF$theta == 9,"theta"] <- paste0('\u03B8'," = 9")
dfF$theta <- factor(dfF$theta, levels = unique(dfF$theta))

yy <- colorspace::sequential_hcl(palette = "Light Grays", n = (3 + 1))
nn <- c("Normal", c(100, 500,1000))

nameLabel <- expression(hat(mu)[SML])
dfF %>% ggplot() +
stat_function(fun = dnorm, n = 1000, args = list(mean = 0, sd = 1), aes(colour = "Normal", linetype = "Normal"), size = 0.8) +
  geom_density(aes(x = muest, colour = factor(n), linetype = factor(n)), size = 0.8) +
  xlim(c(-5, 5)) +
  theme_classic(base_size = 20) +
  facet_grid(vars(theta), vars(mu)) +
  scale_linetype_manual(breaks = rev(nn), values = rev(1:length(nn))) +
  scale_color_manual(breaks = rev(nn), values = rev(yy)) +
  labs(y = "Density", x = nameLabel, colour = "n", linetype = "n")

names <- paste0('\u03B8'," = ", theta)
nameLabel <- expression(hat(theta)[SML])
dfF %>% ggplot() +
  stat_function(fun = dnorm, n = 1000, args = list(mean = 0, sd = 1), aes(colour = "Normal", linetype = "Normal"), size = 0.8) +
  geom_density(aes(x = thetaest, colour = factor(n), linetype = factor(n)), size = 0.8) +
  xlim(c(-5, 5)) +
  theme_classic(base_size = 20) +
  facet_grid(mu ~ theta) +
  scale_linetype_manual(breaks = rev(nn), values = rev(1:length(nn))) +
  scale_color_manual(breaks = rev(nn), values = rev(yy)) +
  labs(y = "Density", x = nameLabel, colour = "n", linetype = "n")


dfFT <- dfF[dfF$mu == paste0('\u03BC'," = log(3)"),]

nameLabel <- expression(hat(mu)[SML])
dfFT %>% ggplot() +
  stat_function(fun = dnorm, n = 1000, args = list(mean = 0, sd = 1), aes(colour = "Normal", linetype = "Normal"), size = 0.8) +
  geom_density(aes(x = muest, colour = factor(n), linetype = factor(n)), size = 0.8) +
  xlim(c(-5, 5)) +
  theme_classic(base_size = 20) +
  facet_wrap(~theta) +
  scale_linetype_manual(breaks = rev(nn), values = rev(1:length(nn))) +
  scale_color_manual(breaks = rev(nn), values = rev(yy)) +
  labs(y = "Density", x = nameLabel, colour = "n", linetype = "n")

names <- paste0('\u03B8'," = ", theta)
nameLabel <- expression(hat(theta)[SML])
dfF %>% ggplot() +
  stat_function(fun = dnorm, n = 1000, args = list(mean = 0, sd = 1), aes(colour = "Normal", linetype = "Normal"), size = 0.8) +
  geom_density(aes(x = thetaest, colour = factor(n), linetype = factor(n)), size = 0.8) +
  xlim(c(-5, 5)) +
  theme_classic(base_size = 20) +
  facet_wrap(~theta) +
  scale_linetype_manual(breaks = rev(nn), values = rev(1:length(nn))) +
  scale_color_manual(breaks = rev(nn), values = rev(yy)) +
  labs(y = "Density", x = nameLabel, colour = "n", linetype = "n")


# plots q -----------------------------------------------------------------

#small sample size

rows <- which(allcomb$mu == log(3) & allcomb$n %in% c(20,30,50))
dfrows <- allcomb[rows,]
qdf <- list()
for (i in 1:length(rows)) {
  n = dfrows[i, "n"]
  mu = dfrows[i, "mu"]
  theta = dfrows[i, "theta"]
  estq = model[[rows[i]]]$q.pad.mle
  q = round((3 + theta) / (2 + theta), 2)

  qdf[[i]] <- data.frame(n = rep(n, 10000), q = rep(q), est = estq)
}
unique(qdfN$q)
qdfN <- do.call(rbind, qdf)
qdfN[qdfN$q == 1.47,"q"] <- "q = 1.47"
qdfN[qdfN$q == 1.33,"q"] <- "q = 1.33"
qdfN[qdfN$q == 1.09,"q"] <- "q = 1.09"
qdfN$q <- factor(qdfN$q, levels = unique(qdfN$q))

yy <- colorspace::sequential_hcl(palette = "Light Grays", n = (3 + 1))
nn <- c("Normal", c(20, 30,50))

nameLabel <- expression(hat(q)[SML])
qdfN %>% ggplot() +
  stat_function(fun = dnorm, n = 1000, args = list(mean = 0, sd = 1), aes(colour = "Normal", linetype = "Normal"), size = 0.8) +
  geom_density(aes(x = est, colour = factor(n), linetype = factor(n)), size = 0.8) +
  xlim(c(-5, 5)) +
  theme_classic(base_size = 20) +
  facet_wrap(~q) +
  scale_linetype_manual(breaks = rev(nn), values = rev(1:length(nn))) +
  scale_color_manual(breaks = rev(nn), values = rev(yy)) +
  labs(y = "Density", x = nameLabel, colour = "n", linetype = "n")



# large sample size

rows <- which(allcomb$mu == log(3) & allcomb$n %in% c(100,500,1000))
dfrows <- allcomb[rows,]
qdf <- list()
for (i in 1:length(rows)) {
  n = dfrows[i, "n"]
  mu = dfrows[i, "mu"]
  theta = dfrows[i, "theta"]
  estq = model[[rows[i]]]$q.pad.mle
  q = round((3 + theta) / (2 + theta), 2)

  qdf[[i]] <- data.frame(n = rep(n, 10000), q = rep(q), est = estq)
}
unique(qdfN$q)
qdfN <- do.call(rbind, qdf)
qdfN[qdfN$q == 1.47,"q"] <- "q = 1.47"
qdfN[qdfN$q == 1.33,"q"] <- "q = 1.33"
qdfN[qdfN$q == 1.09,"q"] <- "q = 1.09"
qdfN$q <- factor(qdfN$q, levels = unique(qdfN$q))

yy <- colorspace::sequential_hcl(palette = "Light Grays", n = (3 + 1))
nn <- c("Normal", c(100,500,1000))

nameLabel <- expression(hat(q)[SML])
qdfN %>% ggplot() +
  stat_function(fun = dnorm, n = 1000, args = list(mean = 0, sd = 1), aes(colour = "Normal", linetype = "Normal"), size = 0.8) +
  geom_density(aes(x = est, colour = factor(n), linetype = factor(n)), size = 0.8) +
  xlim(c(-5, 5)) +
  theme_classic(base_size = 20) +
  facet_wrap(~q) +
  scale_linetype_manual(breaks = rev(nn), values = rev(1:length(nn))) +
  scale_color_manual(breaks = rev(nn), values = rev(yy)) +
  labs(y = "Density", x = nameLabel, colour = "n", linetype = "n")
