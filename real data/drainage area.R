# drainage area (mi^2) - McCuen (2016), example 6.3, pag 122
library(qExp)
library(ggplot2)
library(GoFKernel)

dat <- data.frame(x = c(0.30 , 1.19 , 1.70 , 2.18 , 2.30 ,2.82, 3.85, 6.73 , 6.82 , 6.92 , 6.98,8.50,10.4,18.5 , 24.0,24.3 ,26.8 ,28.0 , 30.2 ,39.5 ,54.8 ,98.4))
n <- length(dat$x)
dat |> ggplot(aes(x=x)) +
  geom_histogram(bins = 10, aes(y = after_stat(density)),fill="white",color="black")+
  ylab("Density") + xlab(expression("Drainage area (mi"^2~")")) +
  coord_cartesian(ylim = c(0,0.042))+
  theme_classic(base_size = 16)


# Obtain the estimates

### Estimando via momentos fracionários

v <-  0.1
u <- -v
mu_uv <- (mean(dat$x^v))^u/(mean(dat$x^u))^v # T_v^u/T_u^v
theta_hat_mf <- g.theta.inv(mu_uv, u = u, v = v, upper = 10000)
mu_hat_mf <- mean(dat$x)
c(mu_hat_mf,theta_hat_mf)

hist(dat$x, prob = T, breaks = 10)
curve(dtsal(x, shape = theta_hat_mf+1, scale = mu_hat_mf*theta_hat_mf), add = T, col = "blue")

## Estimando via máxima verossimilhança


mu_in <- mean(dat$x)
theta_in <- mean(dat$x^0.1)/(mean(dat$x)^0.1)
tt <- optim(par = c(mu_in,theta_in), logQexp, x = dat$x, n = length(dat$x))
tt$par


hist(dat$x, prob = T,  breaks = 10, main = " ", xlab = "LOSS", ylab = "Densidade")
curve(dtsal(x,shape = tt$par[2]+1, scale = tt$par[1]*tt$par[2]),add = T, col = "blue")
curve(dtsal(x, shape = theta_hat_mf+1, scale = mu_hat_mf*theta_hat_mf), add = T, col = "red")
legend("right", legend=c("MMFM", "MMV"),
       lty=1, col=c("red","blue"), bty="n")


curve(dexp(x,1/mean(dat$x)), add = T, col = "green")

# tests

prEstML <- tt$par[2]+1
secEsML <-tt$par[1]*tt$par[2]

goftest::ad.test(dat$x,ptsal, shape=prEstML, scale = secEsML)
ks.test(dat$x,ptsal, shape=prEstML, scale = secEsML)

prEstMF <- theta_hat_mf+1
secEsMF <- mu_hat_mf*theta_hat_mf

goftest::ad.test(dat$x,ptsal, shape=prEstMF, scale = secEsMF)
ks.test(dat$x,ptsal, shape=prEstMF, scale = secEsMF)


# plots

ggplot(dat, aes(x=x)) +
  geom_histogram(bins = 8, aes(y = ..density..),fill="white",color="black", boundary = 0)+
  stat_function(fun = dtsal, args = list(shape = tt$par[2]+1, scale = tt$par[1]*tt$par[2]),show.legend=T, aes(linetype="solid"),xlim = c(0,120))+ ## ajuste por MV
  stat_function(fun = dtsal, args = list(shape = theta_hat_mf+1, scale = mu_hat_mf*theta_hat_mf),show.legend=T,aes(linetype="dotted"),xlim = c(0,120))+ ## ajuste por MMFM
  scale_linetype_manual("Method",values = c("solid", "dotted"),
                        labels = c("MFMM","MLE"))+
  ylab("Density") + xlab(expression("Drainage area (mi"^2~")")) +
  coord_cartesian(ylim = c(0,0.042))+
  theme_classic(base_size = 16)


# Using ML estimates
paramsDatML  = list(shape = prEstML, scale = secEsML)

qqDatML <- ggplot(dat, aes(sample = x))+
  stat_qq(distribution = qtsal, dparams = paramsDatML)+
  stat_qq_line(distribution = qtsal, dparams = paramsDatML)+
  ylab("Quantis amostrais") + xlab("Quantis teóricos") +
  theme_bw()
qqDatML


# Using MF estimates
paramsDatMF  = list(shape = prEstMF, scale = secEsMF)

di <- 'tsal'
qqDatMF <- ggplot(dat, aes(sample = x))+
  ggplot2::stat_qq(distribution = qtsal, dparams = paramsDatMF)+
  ggplot2::stat_qq_line(distribution = qtsal, dparams = paramsDatMF)+
  ylab("Sample Quantiles") + xlab("Theoretical Quantiles") +
  theme_classic(base_size = 16) +
  qqplotr::stat_qq_band(distribution = di, dparams = paramsDatMF, alpha = 0.3)
qqDatMF


# DescTools::PlotQQ(dat$x,
#                   qdist = function(p) qtsal(p, shape = prEstMF, scale = secEsMF),
#                   main = '', pch=1, cex=1.4, cex.lab = 1.4, cex.axis = 1.4)



## Confidence intervals

mfmmValues <- mfmm.theo(mu = mu_hat_mf, theta = theta_hat_mf, v = v, u = u,  prob = 0.95, n = length(dat$x))
kappa2 <- mfmmValues$gamma2*(mfmmValues$der.inv)^2 # variance
n <- length(dat$x)

kappa2 <- mfmmTheo(mu = mu_hat_mf, theta = theta_hat_mf, v = 0.49, u = 0.5, n = n)$kappa2

mu.sahie.MF.inf <- mu_hat_mf - qnorm(0.975)*(1/sqrt(n))*sqrt((theta_hat_mf*theta_hat_mf)*(beta(3,theta_hat_mf-1)-(theta_hat_mf+1)*beta(2,theta_hat_mf)^2))
mu.sahie.MF.sup <- mu_hat_mf - qnorm(0.025)*(1/sqrt(n))*sqrt((mu_hat_mf*theta_hat_mf)*(beta(3,theta_hat_mf-1)-(theta_hat_mf+1)*beta(2,theta_hat_mf)^2))
IC.mu.sahie.MF <- c(mu.sahie.MF.inf,mu.sahie.MF.sup)
round(IC.mu.sahie.MF,2)

theta.sahie.MF.inf <- theta_hat_mf - qnorm(0.975)*(1/sqrt(n))*sqrt(kappa2)
theta.sahie.MF.sup <- theta_hat_mf - qnorm(0.025)*(1/sqrt(n))*sqrt(kappa2)
IC.theta.sahie.MF <- c(theta.sahie.MF.inf,theta.sahie.MF.sup)
round(IC.theta.sahie.MF,2)

## MMV

varmu <- varTeoMle(tt$par[1], tt$par[2])$var.mu ## Of mu
vartheta <- varTeoMle(tt$par[1], tt$par[2])$var.theta #(theta+1)^2*(theta+2)^2

mu.sahie.MV.inf <- tt$par[1] - qnorm(0.975)*sqrt(varmu/n)
mu.sahie.MV.sup <- tt$par[1] -  qnorm(0.025)*sqrt(varmu/n)
IC.mu.sahie.MV <- c(mu.sahie.MV.inf,mu.sahie.MV.sup)
round(IC.mu.sahie.MV,2)

theta.sahie.MV.inf <- tt$par[2] - qnorm(0.975)*sqrt(vartheta/n)
theta.sahie.MV.sup <- tt$par[2]- qnorm(0.025)*sqrt(vartheta/n)
IC.theta.sahie.MV <- c(theta.sahie.MV.inf,theta.sahie.MV.sup)
round(IC.theta.sahie.MV,2)










