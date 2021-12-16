# plots mfmm
# Plot g(theta) and the values of the statistic Tn = T_v^u/T_u^v

library(colorspace)

points.mfmm <- function(mfmm.est, v, u, xlimit = 100, colour = "black", size = 0.25) {
  lim.sup <- mfmm.est$lim.sup
  mu.uv <- mfmm.est$mu.uv
  prop.rejec <- mfmm.est$prop.rejec
  n <- length(mfmm.est$mu.uv)

  data <- data.frame(x = 0, y = mu.uv)

  p1 <- data %>% ggplot(aes(x = x, y = y)) +
    stat_function(fun = g.theta, args = list(v = v, u = u)) +
    geom_hline(yintercept = lim.sup, linetype = 5, size = 0.5) +
    geom_point(colour = colour, size = size) +
    labs(x = expression(theta), y = bquote(paste("g(", theta, ")"))) +
    xlim(0, xlimit) +
    theme_bw()

  p2 <- p1 + labs(title = bquote(paste("n = ", .(n), ", Rejected proportion = ", .(prop.rejec)))) +
    theme(plot.title = element_text(hjust = 0.55, size = rel(1)))

  return(list(p1, p2))
}


# theta is vector
reject.mfmm <- function(mu, theta, nrow = 1) {
  v <- c(0.1, 0.2, 0.3, 0.4)
  u <- c(0.1, 0.2, 0.3, 0.4)
  x <- expand_grid(v = v, u = u)
  x <- x %>% filter(v < u)
  # y <- -v
  # x <- rbind(x,cbind(v,u = y)) %>% mutate(uv = 1:n())
  x <- x %>% mutate(uv = 1:n())


  # create the limits of the sequences
  lim <- sapply(theta, mfmm.theo, mu = mu, v = v, u = -v)[3, ]
  lim <- lapply(lim, function(x) {
    ll <- ifelse(x < 300, 300, x)
    return(max(ll))
  })

  # create the sequences for each theta
  seqs <- lapply(lim, function(x) round(seq(0, x, length = 1000), 0))
  seqs <- melt(seqs)
  seqs$L1 <- factor(seqs$L1, labels = theta)
  seqs$L1 <- as.numeric(levels(seqs$L1))[seqs$L1]
  colnames(seqs) <- c("n", "theta")

  v.n.theta <- expand_grid(x, seqs) %>% mutate(mu = mu)

  probs <- mapply(mfmm.theo,
    theta = v.n.theta$theta, mu = v.n.theta$mu,
    v = v.n.theta$v, u = v.n.theta$u, n = v.n.theta$n
  )[4, ]


  df <- cbind(v.n.theta, probs = unlist(probs))
  df$theta <- round(df$theta, 1)

  labelNames <- paste0(x$v, ", ", x$u)
  theta.names <- c(
    `0.1` = "theta~`=`~1/9",
    `1` = "theta~`=`~1",
    `9` = "theta~`=`~9"
  )

  yy <- colorspace::sequential_hcl(palette = "Light Grays", n = 6)


  p <- df %>%
    group_by(theta) %>%
    ggplot() +
    geom_line(aes(x = n, y = probs, colour = factor(uv), linetype = factor(uv)), size = 0.5) +
    geom_hline(yintercept = 0.97, linetype = 5, size = 0.5) +
    theme_bw() +
    facet_wrap(~theta,
      scales = "free_x", nrow = nrow,
      labeller = labeller(theta = as_labeller(theta.names, label_parsed))
    ) +
    scale_linetype_manual(values = c(1, 2, 3, 4, 6, 7), labels = labelNames) +
    scale_color_manual(values = yy, labels = labelNames) +
    # scale_color_discrete_sequential(palette = "Light Grays", n = 6, labels = labelNames)+
    labs(
      x = expression(n), y = expression(Phi(sqrt(n) * (l[uv] - mu[uv]) / kappa)),
      colour = "v, u", linetype = "v, u"
    )


  return(p)
}

reject.mfmm(log(3), c(1 / 9, 1, 9))



# -----

plotCompDens <- function(data, theta = TRUE, mfmm = TRUE) {
  if (theta) {
    names <- c(
      `0.1` = "theta~`=`~1/9",
      `1` = "theta~`=`~1",
      `9` = "theta~`=`~9"
    )
  } else {
    names <- c(
      `0.1` = "q~`=`~1.47",
      `1` = "q~`=`~1",
      `9` = "q~`=`~1.09"
    )
  }

  if (theta & mfmm) {
    nameLabel <- expression(hat(theta)[SMF])
  }

  if (theta & mfmm == FALSE) {
    nameLabel <- expression(hat(theta)[SML])
  }

  if (theta == FALSE & mfmm) {
    nameLabel <- expression(hat(q)[SMF])
  }

  if (theta == FALSE & mfmm == FALSE) {
    nameLabel <- expression(hat(q)[SML])
  }

  yy <- colorspace::sequential_hcl(palette = "Light Grays", n = (length(n) + 1))

  p <- data %>% ggplot() +
    stat_function(fun = dnorm, n = 1000, args = list(mean = 0, sd = 1), aes(colour = "Normal", linetype = "Normal"), size = 0.6) +
    geom_density(aes(x = value, colour = factor(n), linetype = factor(n)), size = 0.6) +
    xlim(c(-5, 5)) +
    theme_bw() +
    facet_wrap(~theta.x, labeller = labeller(theta.x = as_labeller(names, label_parsed))) +
    scale_linetype_manual(breaks = rev(n), values = rev(1:(length(n) + 1))) +
    scale_color_manual(breaks = rev(n), values = rev(yy)) +
    labs(y = "Density", x = nameLabel, colour = "n", linetype = "n")

  return(p)
}





# do manunally

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
x$theta <- round(x$theta, 1)
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

yy <- colorspace::sequential_hcl(palette = "Light Grays", n = 5)

theta.names <- c(
  `0.1` = "theta~`=`~1/9",
  `1` = "theta~`=`~1",
  `9` = "theta~`=`~9"
)

q.names <- c(
  `0.1` = "q~`=`~1.47",
  `1` = "q~`=`~1.33",
  `9` = "q~`=`~1.09"
)


p.theta <- theta.hat.st %>% ggplot() +
  stat_function(fun = dnorm, n = 1000, args = list(mean = 0, sd = 1), aes(colour = "Normal", linetype = "Normal"), size = 0.6) +
  geom_density(aes(x = value, colour = factor(n.x), linetype = factor(n.x)), size = 0.6) +
  xlim(c(-5, 5)) +
  theme_bw() +
  facet_wrap(~theta.x, labeller = labeller(theta.x = as_labeller(theta.names, label_parsed))) +
  scale_linetype_manual(breaks = rev(c("Normal", "5000", "1101", "79", "52")), values = rev(c(1, 2, 3, 4, 5))) +
  scale_color_manual(breaks = rev(c("Normal", "5000", "1101", "79", "52")), values = rev(yy)) +
  labs(y = "Density", x = expression(hat(theta)[SMF]), colour = "n", linetype = "n")
p.theta

p.q <- q.hat.st %>% ggplot() +
  geom_density(aes(x = value, colour = factor(n.x), linetype = factor(n.x)), size = 0.6) +
  stat_function(fun = dnorm, n = 1000, args = list(mean = 0, sd = 1), aes(colour = "Normal", linetype = "Normal"), size = 0.6) +
  xlim(c(-5, 5)) +
  theme_bw() +
  facet_wrap(~theta.x, labeller = labeller(theta.x = as_labeller(q.names, label_parsed))) +
  scale_linetype_manual(breaks = rev(c("Normal", "5000", "1101", "79", "52")), values = rev(c(1, 2, 3, 4, 5))) +
  scale_color_manual(breaks = rev(c("Normal", "5000", "1101", "79", "52")), values = rev(yy)) +
  labs(y = "Density", x = expression(hat(q)[SMF]), colour = "n", linetype = "n")
p.q



