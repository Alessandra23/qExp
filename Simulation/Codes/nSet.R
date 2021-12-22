library(devtools)
library(dplyr)
library(tidyverse)
install_github("Alessandra23/qExp")

f.n.uv <- function(mu, theta, v, u, prob = 0.97) {
  (n <- mfmm.theo(mu, theta, v, u, prob = prob))
  return(n$n.min)
}

mu <- log(3)
theta <- c(1 / 9, 1, 9)
v <- c(seq(-0.4, -0.1, 0.1), seq(0.1, 0.5, 0.1))
u <- c(seq(-0.4, -0.1, 0.1), seq(0.1, 0.5, 0.1))
x <- expand_grid(v = v, u = u) %>% filter(v != u)
x <- merge(x, theta) %>%
  mutate(mu = mu) %>%
  rename("theta" = y)
df.n <- x %>% mutate(n = f.n.uv(mu, theta, v, u))

df.n %>%
  group_by(theta) %>%
  slice(which.min(n))

