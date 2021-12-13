# plots mfmm
# Plot g(theta) and the values of the statistic Tn = T_v^u/T_u^v

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
reject.mfmm <- function(mu, theta, allV = FALSE, colour = TRUE, nrow = 1) {
  v <- c(0.1,0.2,0.3,0.4)
  u <- c(0.1, 0.2, 0.3, 0.4)
  x <- expand_grid(v = v, u = u)
  x <- x %>% filter(v<u)
  #y <- -v
  #x <- rbind(x,cbind(v,u = y)) %>% mutate(uv = 1:n())
  x <- x %>% mutate(uv = 1:n())
  x

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
                  v = v.n.theta$v, u = v.n.theta$u, n = v.n.theta$n)[4, ]


  df <- cbind(v.n.theta, probs = unlist(probs))

  labelNames <- paste0(x$v,", ",x$u)



  if (colour) {

    p <- df %>%
      group_by(theta) %>%
      ggplot() +
      geom_line(aes(x = n, y = probs, colour = factor(uv)), size = 0.5) +
      labs(x = expression(n), y = expression(Phi(sqrt(n) * (l[uv] - mu[uv]) / kappa))) +
      geom_hline(yintercept = 0.97, linetype = 5, size = 0.5) +
      theme_bw() +
      facet_wrap(~theta, scales = "free_x", nrow = nrow, labeller = label_bquote(paste(theta, " = ", .(round(theta, 1)))))+
      scale_color_discrete_sequential(name = "v, u", palette = "Plasma", n = 6, labels = labelNames) #+

    } else {
    p <- df %>%
      group_by(theta) %>%
      ggplot() +
      geom_line(aes(x = n, y = probs, linetype = factor(uv))) + # change colour for linetype
      labs(x = expression(n), y = expression(Phi(sqrt(n) * (l[uv] - mu[uv]) / kappa)), linetype = "v, u") +
      geom_hline(yintercept = 0.97, linetype = 5, size = 0.5) +
      theme_bw() +
      facet_wrap(~theta, scales = "free_x", nrow = nrow, labeller = label_bquote(paste(theta, " = ", .(round(theta, 1)))))
  }
p

  return(p)
}

reject.mfmm(log(3), c(1 / 9, 1, 9), colour = TRUE, nrow = 1, allV = TRUE)



# -----


x <- data.frame(theta = rep(c(1/9, 1, 9), 2) , uv= rep(c(1,2,3), 2), v = c(rep(c(0.1,0.2,0.3),2)) ,
                u = rep(c(0.2, 0.3, 0.4), 2), n = c(52, 79, 2000, rep(5000, 3)))
x <- x %>% mutate(mu = mu,
                  N = N)

x
samples <- mapply(qexp.samples, n = x$n, theta = x$theta, mu = x$mu, N = x$N)

get.mfmm.est <- list()

for(i in 1: nrow(x)){
  get.mfmm.est[[i]] <-  mfmm.samples(n = size_sample$n[i], mu = size_sample$mu[i], theta = size_sample$theta[i],
                                     v = size_sample$v[i], u = size_sample$u[i], samples = samples[[i]])

}


theta.hat.st <- lapply(get.mfmm.est, function(y){
  y$theta.hat.pad
})


## repeat each row of data frame 10000 (tamanho de cada elemento da lista)



