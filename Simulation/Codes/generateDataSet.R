library(devtools)
devtools::install_github("Alessandra23/qExp")
library(qExponential)

#save_file = "~/Documents/GitHub/qExp/Simulation/Data"

#N <- 10000
#n <- c(20,100,500,1000, 5000, 10000)
#mu <- c(1/10, log(3), 10)
#theta <- c(1/9, 1, 9, 100)


N <- 10000
n <- c(4179, 9987, 25000)
mu <- log(3)
theta <- 9


N <- 50
n <- 10
mu <- log(3)
theta <- 9

allcomb <- expand.grid(N = N,
                       n = n,
                       mu = mu,
                       theta = theta)


ncomb <- nrow(allcomb)
seed  <- 4

for (i in 1:ncomb) {
  set.seed(seed)
  comb <- allcomb[i,]
  N <- comb$N
  n <- comb$n
  mu <- comb$mu
  theta <- comb$theta

  data <- qexp.samples(N = N,n = n, theta = theta, mu = mu)

  filename <- paste('N', N,
                    'n', n,
                    'mu', round(mu,2),
                    'theta', round(theta,2),
                    sep='')

  data_filename <- paste(filename, '_data.RData', sep='')
  save(data, file = paste(save_file, data_filename, sep=''))

}

