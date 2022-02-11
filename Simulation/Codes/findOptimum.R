library(emdbook)
library(ggplot2)

f.vu <- function(par, mu = 1, theta = 1/9, prob = 0.97){
  v <- par[1]
  u <- par[2]

  ev <- (mu * theta)^v * (theta + 1) * beta(v + 1, theta + 1 - v)
  eu <- (mu * theta)^u * (theta + 1) * beta(u + 1, theta + 1 - u)
  sig.v <- (mu * theta)^(2 * v) * (theta + 1) * (beta(2 * v + 1, theta + 1 - 2 * v) - (theta + 1) * (beta(v + 1, theta + 1 - v)^2))
  sig.u <- (mu * theta)^(2 * u) * (theta + 1) * (beta(2 * u + 1, theta + 1 - 2 * u) - (theta + 1) * (beta(u + 1, theta + 1 - u)^2))
  sig.uv <- ((mu * theta)^(u + v)) * (theta + 1) * (beta(u + v + 1, theta + 1 - u - v) - (theta + 1) * beta(u + 1, theta + 1 - u) * beta(v + 1, theta + 1 - v))
  p1 <- (u^2) * sig.v * ev^(2 * u - 2) * eu^(-2 * v)
  p2 <- u * v * ev^(2 * u - 1) * eu^(-2 * v - 1) * sig.uv
  p3 <- (v^2) * sig.u * ev^(2 * u) * eu^(-2 * v - 2)
  gamma2 <- p1 - 2 * p2 + p3
  Esp.Tn <- (ev^u) / (eu^v)
  lim.sup <- (gamma(v+1)^u/gamma(u+1)^v)

  return(ceiling((qnorm(prob)*sqrt(gamma2)/(lim.sup-Esp.Tn))^2))
}


plotContour <- function(contour, xlim, ylim){
  x <- contour$x
  y <- contour$y
  z <- as.vector(contour$z)
  df <- tibble(expand_grid(x,y), z)
  #yy <- colorspace::sequential_hcl(palette = "Light Grays", n = n)

  # p <- ggplot(df, aes(x, y, z = z))+ geom_contour_filled(bins = n) +
  #   scale_fill_manual(values = rev(yy)) +theme_classic(base_size = 15) +
  #   labs(x = "v", y = "u", fill = "n") +
  #   scale_x_continuous(breaks = xlim)+
  #   scale_y_continuous(breaks = ylim)

  rng = range(z[!is.na(z)& z!=Inf])

  p <- ggplot(df, aes(x, y, fill= z)) +
    geom_tile() +
    scale_fill_gradient(low = "grey92",
                        high = "black",
                        na.value = "white") +
    labs(x = "v", y = "u", fill = "n") +
    theme_classic(base_size = 20) +
    #theme(legend.position = "none")+
    scale_x_continuous(breaks = xlim) +
    scale_y_continuous(breaks = ylim)


  print(p)
}


cc <- curve3d(f.vu(c(x,y), theta = 1/9), xlim = c(-0.5,0.5), ylim = c(-0.5, 0.5),
              sys3d = "image", xlab = "v", ylab = "u")
plotContour(cc, xlim = seq(-0.5,0.5,0.1), ylim = seq(-0.5,0.5,0.1))

cc <- curve3d(f.vu(c(x,y), theta = 1/9), xlim = c(0,0.5), ylim = c(0, 0.5),
              sys3d = "image", xlab = "v", ylab = "u")
plotContour(cc, xlim = seq(0,0.5,0.1), ylim = seq(0,0.5,0.1))


# contour <- cc
#
# x <- contour$x
# y <- contour$y
# z <- as.vector(contour$z)
# df <- tibble(expand_grid(x,y), z)
# yy <- colorspace::sequential_hcl(palette = "Light Grays", n = 6)
#
# p <- ggplot(df, aes(x, y, z = z))+ geom_contour(colour = "black", bins = 100)  +
#  # geom_contour_filled()+
#   scale_fill_manual(values = rev(yy)) +theme_bw(base_size = 15) +
#   labs(x = "v", y = "u", fill = "n") +
#   scale_x_continuous(breaks = xlim)+
#   scale_y_continuous(breaks = ylim)
# p


cc <- curve3d(f.vu(c(x,y), theta = 1), xlim = c(0,0.5), ylim = c(0, 0.5),
        sys3d = "image", xlab = "v", ylab = "u")

plotContour(cc, xlim = seq(0,0.5,0.1), ylim = seq(0,0.5,0.1))


# contour <- cc
#
# x <- contour$x
# y <- contour$y
# z <- as.vector(contour$z)
# df <- tibble(expand_grid(x,y), z)
# yy <- colorspace::sequential_hcl(palette = "Light Grays", n = 6)
#
# p <- ggplot(df, aes(x, y, z = z))+ geom_contour(colour = "black", bins = 12)  +
#   # geom_contour_filled()+
#   scale_fill_manual(values = rev(yy)) +theme_bw(base_size = 15) +
#   labs(x = "v", y = "u", fill = "n") +
#   scale_x_continuous(breaks = xlim)+
#   scale_y_continuous(breaks = ylim)
# p


cc <- curve3d(f.vu(c(x,y), theta = 9), xlim = c(0,0.5), ylim = c(0, 0.5),
        sys3d = "image", xlab = "v", ylab = "u")

plotContour(cc, xlim = seq(0,0.5,0.1), ylim = seq(0,0.5,0.1))


# contour <- cc
#
# x <- contour$x
# y <- contour$y
# z <- as.vector(contour$z)
# df <- tibble(expand_grid(x,y), z)
# yy <- colorspace::sequential_hcl(palette = "Light Grays", n = 6)
#
# p <- ggplot(df, aes(x, y, z = z))+ geom_contour(colour = "black", bins = 6)  +
#   # geom_contour_filled()+
#   scale_fill_manual(values = rev(yy)) +theme_bw(base_size = 15) +
#   labs(x = "v", y = "u", fill = "n") +
#   scale_x_continuous(breaks = xlim)+
#   scale_y_continuous(breaks = ylim)
# p



# Código antigo e optim ---------------------------------------------------


cc <- curve3d(f.vu(c(x,y), theta = 1/9), xlim = c(-0.1,0.5), ylim = c(-0.1, 0.5),
              sys3d = "image", xlab = "v", ylab = "u")
cc



# Quando faço o gráfico com intervalos menores, vejo que quando theta é pequeno,
# valores pequenos de 'v' e 'u' me dão os menores tamanhos amostrais.
# Já quando theta theta cresce, o oposto ocorre, 'v' e 'u' grabdes me dão os menores tamanhos amostrais.
# Fazer o gráfico para isso:

# Todo o intervalo
par(mfrow = c(1,3))
cc <- curve3d(f.vu(c(x,y), theta = 1/9), xlim = c(-0.5,0.5), ylim = c(-0.5, 0.5),
              sys3d = "image", xlab = "v", ylab = "u")
plotContour(cc, xlim = c(-0.5,0.5), ylim = c(-0.5, 0.5))
curve3d(f.vu(c(x,y), theta = 1), xlim = c(-0.5,0.5), ylim = c(-0.5, 0.5),
        sys3d = "image", xlab = "v", ylab = "u")
curve3d(f.vu(c(x,y), theta = 9), xlim = c(-0.5,0.5), ylim = c(-0.5, 0.5),
        sys3d = "image", xlab = "v", ylab = "u")


# Intervalo menor
# v>0
curve3d(f.vu(c(x,y), theta = 1/9), xlim = c(0,0.5), ylim = c(-0.5, 0.5),
        sys3d = "image", xlab = "v", ylab = "u")
curve3d(f.vu(c(x,y), theta = 1), xlim = c(0,0.5), ylim = c(-0.5, 0.5),
        sys3d = "image", xlab = "v", ylab = "u")
curve3d(f.vu(c(x,y), theta = 9), xlim = c(0,0.5), ylim = c(-0.5, 0.5),
        sys3d = "image", xlab = "v", ylab = "u")

# u>0
curve3d(f.vu(c(x,y), theta = 1/9), xlim = c(-0.5,0.5), ylim = c(0, 0.5),
        sys3d = "image", xlab = "v", ylab = "u")
curve3d(f.vu(c(x,y), theta = 1), xlim = c(-0.5,0.5), ylim = c(0, 0.5),
        sys3d = "image", xlab = "v", ylab = "u")
curve3d(f.vu(c(x,y), theta = 9), xlim = c(-0.5,0.5), ylim = c(0, 0.5),
        sys3d = "image", xlab = "v", ylab = "u")

# v,u>0
curve3d(f.vu(c(x,y), theta = 1/9), xlim = c(0,0.5), ylim = c(0, 0.5),
        sys3d = "image", xlab = "v", ylab = "u")
curve3d(f.vu(c(x,y), theta = 1), xlim = c(0,0.5), ylim = c(0, 0.5),
        sys3d = "image", xlab = "v", ylab = "u")
curve3d(f.vu(c(x,y), theta = 9), xlim = c(0,0.5), ylim = c(0, 0.5),
        sys3d = "image", xlab = "v", ylab = "u")


par(mfrow = c(2,5))
cc <- curve3d(f.vu(c(x,y), theta = 1), xlim = c(0,0.5), ylim = c(0, 0.5),
              sys3d = "image", xlab = "v", ylab = "u")
cc$x

df <- data.frame(v = cc$x, u = cc$y, n = cc$z )

df

v <- ggplot(df, aes(v, u, z = n))
v + geom_contour() + theme_bw()

# find the optimum --------------------------------------------------------

optimx(par = c(0.25, 0.4), fn = f.vu,  mu = 1, theta = 9)
optimx(par = c(0.1, 0.5), fn = f.vu,  lower = c(0.0001, 0.0002), upper = c(0.49, 0.5),  mu = 1, theta = 1/9)


oo = optim(par = c(0.1, 0.2), f.vu, theta = 1/9)
oo



oo = optim(par = c(0.3, 0.4), f.vu, theta = 1/9,
           lower = c(0.19, 0.2),
           upper = c(0.49, 0.5), method="L-BFGS-B")

oo



oo = optim(par = c(0.3, 0.4), f.vu, theta = 1/9,
           lower = c(0,0.001),
           upper = c(0.49, 0.5), method="L-BFGS-B")

oo


oo <- optimx(par = c(0.1, 0.5), fn = f.vu,  lower = c(0.0001, 0.0002), upper = c(0.49, 0.5),  mu = 1, theta = 1/9)
oo
c(round(oo$p1,2), round(oo$p2,2))
f.vu(c(round(oo$p1,2), round(oo$p2,2)), theta = 1/9)

oo <- optimx(par = c(0.17, 0.30), fn = f.vu,  lower = c(0.0001, 0.0002), upper = c(0.49, 0.5),  mu = 1, theta = 1)
oo
c(round(oo$p1,2), round(oo$p2,2))
f.vu(c(oo$p1, oo$p2), theta = 1)

oo <- optimx(par = c(0.4, 0.45), fn = f.vu,  lower = c(0.0001, 0.0002), upper = c(0.49, 0.5),  mu = 1, theta = 9)
oo
f.vu(c(oo$p1, oo$p2), theta = 9)

f.vu(c(0.499, 0.5), theta = 9)




