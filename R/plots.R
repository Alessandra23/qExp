#' points mfmm
#' Plot points
#' @import colorspace
#' @import tidyr
#' @param  v tuning parameter
#' @param u tuning parameter
#' @param mfmm.est object obtained from mfmm.samples function
#' @export
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


#' Reject mfmm
#'
#' Plot asymptotic approximation for different values of n
#'
#' @param mu paramter mu
#' @param theta vector of values of theta
#' @param nrow number of rows of facet grid
reject.mfmm <- function(mu, theta, nrow = 1) {
  v <- c(0.1, 0.2, 0.3, 0.4)
  u <- c(0.1, 0.2, 0.3, 0.4)
  x <- expand_grid(v = v, u = u)
  x <- x %>% filter(v < u)
  # y <- -v
  # x <- rbind(x,cbind(v,u = y)) %>% mutate(uv = 1:n())
  x <- x %>% mutate(uv = 1:n())


  # create the limits of the sequences
  lim <- sapply(theta, mfmmTheo, mu = mu, v = v, u = -v)[3, ]
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

  probs <- mapply(mfmmTheo,
                  theta = v.n.theta$theta, mu = v.n.theta$mu,
                  v = v.n.theta$v, u = v.n.theta$u, n = v.n.theta$n
  )[5,]


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
    geom_line(aes(x = n, y = probs, colour = factor(uv), linetype = factor(uv)), linewidth = 0.8) +
    geom_hline(yintercept = 0.97, linetype = 5, size = 0.5) +
    theme_classic(base_size = 16) +
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

#' Plot densities
#' Comparison of densities of mfmm estimators and Standart Normal
#' @param data data frame
#' @param n vector of samples sizes
#' @param theta vector of values of theta (as character)
#' @param parameter if is "theta" or "q"
#' @param standard if TRUE plot the standard stimates
#' @export
plotCompDens.mfmm <- function(data, n, theta, parameter = "theta", standard = TRUE) {


  if (parameter == "theta") {
    names <- paste0('\u03B8'," = ", theta)
    #names <- paste0(expression(theta)," = ", theta)
  } else {
    vTheta <-  unique(data$theta)
    q <- sapply(vTheta, function(x) {
      round((3 + x) / (2 + x), 2)
    })
    names <- paste0("q = ", q)
  }

  custLab <- function (x){
    names
  }

  if (parameter == "theta" & standard) {
    nameLabel <- expression(hat(theta)[SMF])
  }

  if (parameter == "theta" & standard == FALSE) {
    nameLabel <- expression(hat(theta)[MF])
  }

  if (parameter == "q" & standard) {
    nameLabel <- expression(hat(q)[SMF])
  }

  if (parameter == "q" & standard == FALSE) {
    nameLabel <- expression(hat(q)[MF])
  }


  yy <- colorspace::sequential_hcl(palette = "Light Grays", n = (length(n) + 1))
  nn <- c("Normal", n)

  if(standard){
    p <- data %>% ggplot() +
      stat_function(fun = dnorm, n = 1000, args = list(mean = 0, sd = 1), aes(colour = "Normal", linetype = "Normal"), size = 0.8) +
      geom_density(aes(x = value, colour = factor(n), linetype = factor(n)), size = 0.8) +
      xlim(c(-5, 5)) +
      theme_classic(base_size = 20) +
      facet_wrap(~theta, labeller =  as_labeller(custLab)) +
      scale_linetype_manual(breaks = rev(nn), values = rev(1:length(nn))) +
      scale_color_manual(breaks = rev(nn), values = rev(yy)) +
      labs(y = "Density", x = nameLabel, colour = "n", linetype = "n")
  } else{
    p <- data %>% ggplot() +
      geom_density(aes(x = value, colour = factor(n), linetype = factor(n)), size = 0.8) +
      xlim(c(0, 10)) +
      theme_classic(base_size = 20) +
      facet_wrap(~theta, labeller =  as_labeller(custLab)) +
      scale_linetype_manual(breaks = rev(n), values = rev(1:length(n))) +
      scale_color_manual(breaks = rev(n), values = rev(yy)) +
      labs(y = "Density", x = nameLabel, colour = "n", linetype = "n")
  }

  return(p)
}
