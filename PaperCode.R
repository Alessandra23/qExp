# Code to generate the paper results

#devtools::install_github("Alessandra23/qExp")
library(qExp)
library(dplyr)
library(ggplot2)
library(reshape2)

# Figure 1 ----------------------------------------------------------------

theta.1 = 0.5
theta.2 = 9

# Fig 1(a)

# setting mu = 1 and mu = 10

x <- seq(0,20,length=1000)
v1.mt.1 <- dtsal(x,shape = theta.1+1, scale = 1*theta.1)
v1.mt.2 <- dtsal(x,shape = theta.2+1, scale = 1*theta.2)
v1.mt.3 <- dtsal(x,shape = theta.1+10, scale = 10*theta.1)
v1.mt.4 <- dtsal(x,shape = theta.2+10, scale = 10*theta.2)

df.mt.paper <- data.frame(x = x, v1 = v1.mt.1, v2 = v1.mt.2, v3 = v1.mt.3, v4 = v1.mt.4)


legenda.paper <- c(expression(paste(mu, " = 1.0,  ", ~theta," = 0.5 "),
                              paste(~mu, " = 1.0,  ", ~theta," = 9.0 "),
                              paste(~mu, " = 10,   ", ~theta," = 0.5 "),
                              paste(~mu, " = 10,   ", ~theta," = 9.0 ")))

p.mt.paper <- ggplot(melt(df.mt.paper,id.vars='x'), aes(x,value,group=variable))+
  geom_line(aes(linetype = variable)) + ylab("f(x)")+
  theme_classic(base_size = 22)+
  scale_linetype_discrete(labels=legenda.paper, name=" ")+
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    legend.title = element_text(size = 8),
    legend.text.align = 1
    #legend.box.background = element_rect(colour = "black")
  )


p.mt.paper

# Fig 1(b)

theta = 1 ; mu = 1
x = seq(0,15,.10)
d.qexp = dtsal(x,shape = theta+1,scale = mu*theta)
d.exp = dexp(x,1)

df.comp <- data.frame(x,d.qexp,d.exp)
p.comp.1 <- ggplot(melt(df.comp,id.vars='x'), aes(x,value,group=variable)) +
  geom_line(aes(linetype = variable)) +
  scale_linetype_discrete("", labels = c("q-Exponential", "Exponential"))+
  xlab("x") + ylab("f(x)") +
  theme_classic(base_size = 22)+
  theme(legend.position=c(0.8,0.8))
#legend.box.background = element_rect(colour = "black"))

p.comp.1
p.comp.2 <- ggplot(melt(df.comp,id.vars='x'), aes(x,value,group=variable)) +
  geom_line(aes(linetype = variable), show.legend = FALSE) +
  coord_cartesian(xlim=c(7,15), ylim = c(0,0.005)) +
  #scale_linetype_discrete("", labels = c("q-Exponential", "Exponential")) +
  # xlab("x") + ylab("f(x)") +
  xlab("") + ylab("") +
  theme(legend.position=c(0.8,0.8)) +
  theme_classic(base_size = 22)

p.comp.2

p.comp.1 + annotation_custom(ggplotGrob(p.comp.2), xmin = 6, xmax = 15,
                             ymin = 0.01, ymax = 1)


# Figure 2 ----------------------------------------------------------------

v <- c(seq(-0.4, -0.1, 0.1), seq(0.1, 0.5, 0.1))
u <- c(seq(-0.4, -0.1, 0.1), seq(0.1, 0.5, 0.1))
x <- expand_grid(v = v, u = u) %>% filter(v != u)
theta <- c(0:100)
df <- expand_grid(x, theta)
df <- df %>% mutate(g = g.theta(theta, v, u),
                    vu = paste0(v,", ", u))
# Fig 2(a)
df %>% filter(v>0&u>0&v<u) %>%
  ggplot() + geom_line(aes(x = theta, y = g, linetype = vu)) + theme_classic(base_size = 22)+
  labs(x=expression(theta), y = bquote(paste("g(",theta,")")), linetype= "v,u")

# Fig 2(b)
df %>% filter(v>0&u>0&v>u) %>%
  ggplot() + geom_line(aes(x = theta, y = g, linetype = vu)) + theme_classic(base_size = 22) +
  labs(x=expression(theta), y = " ",  linetype= "v,u")


# Figure 3 ----------------------------------------------------------------

# The function f.vu is in the file qExp/R/gThetaFunction and the plotContour function is in the file
# qExp/R/plots

# theta = 1/9
cc_1 <- curve3d(f.vu(c(x,y), theta = 1/9), xlim = c(-0.5,0.5), ylim = c(-0.5, 0.5),
                sys3d = "image", xlab = "v", ylab = "u")
# Fig 3(a)
plotContour(contour = cc_1, xlim = 0, ylim = 2500, type = 2)

cc_1_zoom <- curve3d(f.vu(c(x,y), theta = 1/9), xlim = c(0,0.5), ylim = c(0, 0.5),
                     sys3d = "image", xlab = "v", ylab = "u")
# Fig 3(d)
plotContour(cc_1_zoom, xlim = seq(0,0.5,0.1), ylim = seq(0,0.5,0.1))


# theta = 1
cc_2 <- curve3d(f.vu(c(x,y), theta = 1), xlim = c(-0.5,0.5), ylim = c(-0.5, 0.5),
                sys3d = "image", xlab = "v", ylab = "u")
# Fig 3(b)
plotContour(contour = cc_2, xlim = 0, ylim = 2500, type = 2)

cc_2_zoom <- curve3d(f.vu(c(x,y), theta = 1), xlim = c(0,0.5), ylim = c(0, 0.5),
                     sys3d = "image", xlab = "v", ylab = "u")
# Fig 3(e)
plotContour(cc_2_zoom, xlim = seq(0,0.5,0.1), ylim = seq(0,0.5,0.1))


# theta = 9
cc_3 <- curve3d(f.vu(c(x,y), theta = 9), xlim = c(-0.5,0.5), ylim = c(-0.5, 0.5),
                sys3d = "image", xlab = "v", ylab = "u")
# Fig 3(c)
plotContour(contour = cc_3, xlim = 0, ylim = 2500, type = 2)

cc_3_zoom <- curve3d(f.vu(c(x,y), theta = 9), xlim = c(0,0.5), ylim = c(0, 0.5),
                     sys3d = "image", xlab = "v", ylab = "u")
# Fig 3(f)
plotContour(cc_3_zoom, xlim = seq(0,0.5,0.1), ylim = seq(0,0.5,0.1))


# Figures 4, 5 and 16 ----------------------------------------------------------------

# Generate the data
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
#
#create samples
# set.seed(2022)
# samples.mfmm <- mapply(qexp.samples, n = x$n, theta = x$theta, mu = x$mu, N = x$N)
# save(samples.mfmm, file="data/samples.mfmm.RData")
# mfmm.est <- list()
# for (i in 1:nrow(x)) {
#   mfmm.est[[i]] <- mfmmSamples(
#     n = x$n[i], mu = x$mu[i], theta = x$theta[i],
#     v = x$v[i], u = x$u[i], samples = samples.mfmm[[i]]
#   )
# }
#save(mfmm.est, file="data/mfmm.est.RData")

# load data
data('mfmm.est')
# Fig 4 (theta hat MFMM)
theta.hat <- lapply(mfmm.est, function(y) {
  y$theta.hat
}) |> reshape2::melt()
colnames(theta.hat) <- c("value", "id")
theta.hat <- merge(x, theta.hat)

theta = c("1/9", "1", "9")
n <- c("5000", "817", "78", "48")
plotCompDens.mfmm(data = theta.hat, n = n, theta = theta, parameter = "theta", standard = FALSE)


# Fig 5 (q hat MFMM standardized)
q.hat.st <- lapply(mfmm.est, function(y) {
  y$q.hat.pad
}) |> melt()
colnames(q.hat.st) <- c("value", "id")
q.hat.st <- merge(x, q.hat.st)

theta = c("1/9", "1", "9")
n <- c("5000", "817", "78", "48")
plotCompDens.mfmm(data = q.hat.st, n = n, theta = theta, parameter = "q", standard = TRUE)


#Fig 16 (theta hat MFMM standardized)
theta.hat.st <- lapply(mfmm.est, function(y) {
  y$theta.hat.pad
}) |> melt()
colnames(theta.hat.st) <- c("value", "id")
theta.hat.st <- merge(x, theta.hat.st)

theta = c("1/9", "1", "9")
n <- c("5000", "817", "78", "48")
plotCompDens.mfmm(data = theta.hat.st, n = n, theta = theta, parameter = "theta", standard = TRUE)

#Fig 17 (q hat MFMM)
# q.hat <- lapply(mfmm.est, function(y) {
#   y$q.hat
# }) |> melt()
# colnames(q.hat) <- c("value", "id")
# q.hat <- merge(x, q.hat)
#
# theta = c("1/9", "1", "9")
# n <- c("5000", "817", "78", "48")
# plotCompDens.mfmm(data = q.hat, n = n, theta = theta, parameter = "q", standard = FALSE)


# Figure 6, 7, 8 and 9 ----------------------------------------------------------------

N = 10000
n <- c(20, 30, 50, 100, 500,1000)
mu <- c(1/10, log(3), 10)
theta <- c(1/9, 1, 9)
x.mle <- expand.grid(N = N,
                     n = n,
                     mu = mu,
                     theta = theta)

# #create samples
# set.seed(2022)
# samples.mle <- mapply(qexp.samples, n = x.mle$n, theta = x.mle$theta, mu = x.mle$mu, N = x.mle$N)
#
# mle.est <- list()
# for (i in 1:length(samples.mle)) {
#   mle.est[[i]] <- mleSamples(samples.mle[[i]], mu = x.mle$mu[i], theta = x.mle$theta[i])
# }
# save(mle.est, file="data/mle.est.RData")

# load data
data('mle.est')

set.seed(2022)
mu = log(3)
theta = 9
samp <- qexp.samples(n = 100, theta = theta, mu = mu, N = 10000)
ests <- mleSamples(samp, mu = mu, theta = theta)

ests$theta.pad.mle |> as.data.frame() |>
  ggplot() +
  stat_function(fun = dnorm, n = 1000, args = list(mean = 0, sd = 1), aes(colour = "Normal"), linewidth = 0.8) +
  geom_density(aes(x = ests$theta.pad.mle), linewidth = 0.8) +
  xlim(c(-5, 5))


# Fig 6

rows <- which(x.mle$n %in% c(100,500,1000))
dfrows <- x.mle[rows,]
mudf <- list()
thetadf <- list()
for (i in 1:length(rows)) {
  n = dfrows[i, "n"]
  mu = dfrows[i, "mu"]
  theta = dfrows[i, "theta"]
  estmu = mle.est[[rows[i]]]$mu.pad.mle
  esttheta = mle.est[[rows[i]]]$theta.pad.mle

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

#yy <- colorspace::sequential_hcl(palette = "Light Grays", n = (3 + 1))
yy <- colorRampPalette(c("black","gray80"))
yy <- yy(3 + 1)
nn <- c("Normal", c(100, 500,1000))

# Fig 17
nameLabel <- expression(hat(mu)[SML])
dfF %>% ggplot() +
  stat_function(fun = dnorm, n = 1000, args = list(mean = 0, sd = 1), aes(colour = "Normal", linetype = "Normal"), size = 0.8) +
  geom_density(aes(x = muest, colour = factor(n), linetype = factor(n)), size = 0.8) +
  xlim(c(-5, 5)) +
  theme_classic(base_size = 22) +
  facet_wrap(~mu) +
  scale_linetype_manual(breaks = rev(nn), values = rev(1:length(nn))) +
  scale_color_manual(breaks = rev(nn), values = rev(yy)) +
  labs(y = "Density", x = nameLabel, colour = "n", linetype = "n")


# Fig 6
dfFT <- dfF[dfF$mu == paste0('\u03BC'," = log(3)"),]
nameLabel <- expression(hat(mu)[SML])
dfFT %>% ggplot() +
  stat_function(fun = dnorm, n = 1000, args = list(mean = 0, sd = 1), aes(colour = "Normal", linetype = "Normal"), size = 0.8) +
  geom_density(aes(x = muest, colour = factor(n), linetype = factor(n)), size = 0.8) +
  xlim(c(-5, 5)) +
  theme_classic(base_size = 22) +
  facet_wrap(~theta) +
  scale_linetype_manual(breaks = rev(nn), values = rev(1:length(nn))) +
  scale_color_manual(breaks = rev(nn), values = rev(yy)) +
  labs(y = "Density", x = nameLabel, colour = "n", linetype = "n")


# Fig 7
names <- paste0('\u03B8'," = ", theta)
nameLabel <- expression(hat(theta)[SML])
dfF %>% ggplot() +
  stat_function(fun = dnorm, n = 1000, args = list(mean = 0, sd = 1), aes(colour = "Normal", linetype = "Normal"), size = 0.8) +
  geom_density(aes(x = thetaest, colour = factor(n), linetype = factor(n)), size = 0.8) +
  xlim(c(-5, 5)) +
  theme_classic(base_size = 22) +
  facet_wrap(~theta) +
  scale_linetype_manual(breaks = rev(nn), values = rev(1:length(nn))) +
  scale_color_manual(breaks = rev(nn), values = rev(yy)) +
  labs(y = "Density", x = nameLabel, colour = "n", linetype = "n")


# Fig 8
rows <- which(x.mle$mu == log(3) & x.mle$n %in% c(20,30,50))
dfrows <- x.mle[rows,]
qdf <- list()
for (i in 1:length(rows)) {
  n = dfrows[i, "n"]
  mu = dfrows[i, "mu"]
  theta = dfrows[i, "theta"]
  estq = mle.est[[rows[i]]]$q.pad.mle
  q = round((3 + theta) / (2 + theta), 2)

  qdf[[i]] <- data.frame(n = rep(n, 10000), q = rep(q), est = estq)
}
unique(qdfN$q)
qdfN <- do.call(rbind, qdf)
qdfN[qdfN$q == 1.47,"q"] <- "q = 1.47"
qdfN[qdfN$q == 1.33,"q"] <- "q = 1.33"
qdfN[qdfN$q == 1.09,"q"] <- "q = 1.09"
qdfN$q <- factor(qdfN$q, levels = unique(qdfN$q))

#yy <- colorspace::sequential_hcl(palette = "Light Grays", n = (3 + 1))
yy <- colorRampPalette(c("black","gray80"))
yy <- yy(3 + 1)
nn <- c("Normal", c(20, 30,50))

nameLabel <- expression(hat(q)[SML])
qdfN %>% ggplot() +
  stat_function(fun = dnorm, n = 1000, args = list(mean = 0, sd = 1), aes(colour = "Normal", linetype = "Normal"), size = 0.8) +
  geom_density(aes(x = est, colour = factor(n), linetype = factor(n)), size = 0.8) +
  xlim(c(-5, 5)) +
  theme_classic(base_size = 22) +
  facet_wrap(~q) +
  scale_linetype_manual(breaks = rev(nn), values = rev(1:length(nn))) +
  scale_color_manual(breaks = rev(nn), values = rev(yy)) +
  labs(y = "Density", x = nameLabel, colour = "n", linetype = "n")


# Fig 9
rows <- which(x.mle$mu == log(3) & x.mle$n %in% c(100,500,1000))
dfrows <- x.mle[rows,]
qdf <- list()
for (i in 1:length(rows)) {
  n = dfrows[i, "n"]
  mu = dfrows[i, "mu"]
  theta = dfrows[i, "theta"]
  estq = mle.est[[rows[i]]]$q.pad.mle
  q = round((3 + theta) / (2 + theta), 2)

  qdf[[i]] <- data.frame(n = rep(n, 10000), q = rep(q), est = estq)
}
unique(qdfN$q)
qdfN <- do.call(rbind, qdf)
qdfN[qdfN$q == 1.47,"q"] <- "q = 1.47"
qdfN[qdfN$q == 1.33,"q"] <- "q = 1.33"
qdfN[qdfN$q == 1.09,"q"] <- "q = 1.09"
qdfN$q <- factor(qdfN$q, levels = unique(qdfN$q))

#yy <- colorspace::sequential_hcl(palette = "Light Grays", n = (3 + 1))
yy <- colorRampPalette(c("black","gray80"))
yy <- yy(3 + 1)
nn <- c("Normal", c(100,500,1000))

nameLabel <- expression(hat(q)[SML])
qdfN %>% ggplot() +
  stat_function(fun = dnorm, n = 1000, args = list(mean = 0, sd = 1), aes(colour = "Normal", linetype = "Normal"), size = 0.8) +
  geom_density(aes(x = est, colour = factor(n), linetype = factor(n)), size = 0.8) +
  xlim(c(-5, 5)) +
  theme_classic(base_size = 22) +
  facet_wrap(~q) +
  scale_linetype_manual(breaks = rev(nn), values = rev(1:length(nn))) +
  scale_color_manual(breaks = rev(nn), values = rev(yy)) +
  labs(y = "Density", x = nameLabel, colour = "n", linetype = "n")


# Figure 10 ---------------------------------------------------------------


# Figure 11 ---------------------------------------------------------------


# Figure 12 ---------------------------------------------------------------


# Figure 13 ---------------------------------------------------------------


# Figure 14 ---------------------------------------------------------------


# Figure 15 ---------------------------------------------------------------

# Figure 16 ---------------------------------------------------------------

# Figure 17 ---------------------------------------------------------------


