library(emdbook)
library(ggplot2)
library(optimx)

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



plotContour <- function(contour, xlim, ylim, type = 2){
  x <- contour$x
  y <- contour$y
  z <- as.vector(contour$z)
  df <- tibble::tibble(tidyr::expand_grid(x,y), z)
  #yy <- colorspace::sequential_hcl(palette = "Light Grays", n = n)

  # p <- ggplot(df, aes(x, y, z = z))+ geom_contour_filled(bins = n) +
  #   scale_fill_manual(values = rev(yy)) +theme_classic(base_size = 15) +
  #   labs(x = "v", y = "u", fill = "n") +
  #   scale_x_continuous(breaks = xlim)+
  #   scale_y_continuous(breaks = ylim)

  #rng = range(z[!is.na(z)& z!=Inf])

  pal <- colorRampPalette(c("grey90","black"))
  pal <- pal(100)

  if(type == 1){
    p <- ggplot(df, aes(x, y, fill= z)) +
      geom_tile() +
      scale_fill_gradient(low = "grey90",
                          high = "black",
                          na.value = "white",
                          limits = c(1,7500000)) +
      # scale_fill_gradientn(na.value = "white",
      #                      limits = c(1,7500000),
      #                      colors = rev(colorspace::sequential_hcl(palette = 'light grays', n = 100)),
      #                      breaks = c(1, 1875001, 3750000, 5625000, 7500000)) +
      labs(x = "v", y = "u", fill = "n") +
      theme_classic(base_size = 16) +
      #theme(legend.position = "none")+
      scale_x_continuous(breaks = xlim) +
      scale_y_continuous(breaks = ylim)
  }

  if(type == 2){
    p <- ggplot(df, aes(x, y, fill= z)) +
      geom_tile() +
      scale_fill_gradient(low = "grey90",
                          high = "black",
                          na.value = "white") +
      labs(x = "v", y = "u", fill = "n") +
      theme_classic(base_size = 16) +
      #theme(legend.position = "none")+
      scale_x_continuous(breaks = xlim) +
      scale_y_continuous(breaks = ylim)
  }

  if(type == 3){
    lims <- c(xlim, ylim)
    p <-   ggplot(df, aes(x = x, y = y, fill = z)) +
      geom_tile() +
      scale_fill_gradientn(colors = pal,
                           limits = lims,
                           name = 'n',
                           oob = scales::squish,
                           na.value = "white",
                           guide = guide_colorbar(
                             order = 2,
                             frame.colour = "white",
                             ticks.colour = "white"
                           )
      ) +
      labs(x = "v", y = "u") +
      theme_classic(base_size = 16) +
      #theme(legend.position = "none")+
      scale_x_continuous(breaks = seq(-0.5,0.5, 0.1)) +
      scale_y_continuous(breaks = seq(-0.5,0.5, 0.1))
  }


  print(p)
}


# theta = 1/9
cc <- curve3d(f.vu(c(x,y), theta = 1/9), xlim = c(-0.5,0.5), ylim = c(-0.5, 0.5),
              sys3d = "image", xlab = "v", ylab = "u")
plotContour(cc, xlim = seq(-0.5,0.5,0.1), ylim = seq(-0.5,0.5,0.1))

cc_1_zoom <- curve3d(f.vu(c(x,y), theta = 1/9), xlim = c(0,0.5), ylim = c(0, 0.5),
                     sys3d = "image", xlab = "v", ylab = "u")
plotContour(cc_1_zoom, xlim = seq(0,0.5,0.1), ylim = seq(0,0.5,0.1))


# theta = 1
cc <- curve3d(f.vu(c(x,y), theta = 1), xlim = c(-0.5,0.5), ylim = c(-0.5, 0.5),
              sys3d = "image", xlab = "v", ylab = "u")
plotContour(cc, xlim = seq(-0.5,0.5,0.1), ylim = seq(-0.5,0.5,0.1))

cc_2_zoom <- curve3d(f.vu(c(x,y), theta = 1), xlim = c(0,0.5), ylim = c(0, 0.5),
                     sys3d = "image", xlab = "v", ylab = "u")

plotContour(cc_2_zoom, xlim = seq(0,0.5,0.1), ylim = seq(0,0.5,0.1))


# theta = 9
cc <- curve3d(f.vu(c(x,y), theta = 9), xlim = c(-0.5,0.5), ylim = c(-0.5, 0.5),
              sys3d = "image", xlab = "v", ylab = "u")
plotContour(cc, xlim = seq(-0.5,0.5,0.1), ylim = seq(-0.5,0.5,0.1))

cc_3_zoom <- curve3d(f.vu(c(x,y), theta = 9), xlim = c(0,0.5), ylim = c(0, 0.5),
                     sys3d = "image", xlab = "v", ylab = "u")
plotContour(cc_3_zoom, xlim = seq(0,0.5,0.1), ylim = seq(0,0.5,0.1))



## just the entire range

# theta = 1/9
cc_1 <- curve3d(f.vu(c(x,y), theta = 1/9), xlim = c(-0.5,0.5), ylim = c(-0.5, 0.5),
                sys3d = "image", xlab = "v", ylab = "u")
plotContour(cc_1, xlim = seq(-0.5,0.5,0.1), ylim = seq(-0.5,0.5,0.1))

# theta = 1
cc_2 <- curve3d(f.vu(c(x,y), theta = 1), xlim = c(-0.5,0.5), ylim = c(-0.5, 0.5),
                sys3d = "image", xlab = "v", ylab = "u")
plotContour(cc_2, xlim = seq(-0.5,0.5,0.1), ylim = seq(-0.5,0.5,0.1))


# theta = 9
cc_3 <- curve3d(f.vu(c(x,y), theta = 9), xlim = c(-0.5,0.5), ylim = c(-0.5, 0.5),
                sys3d = "image", xlab = "v", ylab = "u")
plotContour(cc_3, xlim = seq(-0.5,0.5,0.1), ylim = seq(-0.5,0.5,0.1))

library(patchwork)

plotContour(cc_1, xlim = seq(-0.5,0.5,0.1), ylim = seq(-0.5,0.5,0.1), type = 1) +
  plotContour(cc_2, xlim = seq(-0.5,0.5,0.1), ylim = seq(-0.5,0.5,0.1), type = 1) +
  plotContour(cc_3, xlim = seq(-0.5,0.5,0.1), ylim = seq(-0.5,0.5,0.1), type = 1) +

plotContour(cc_1, xlim = seq(-0.5,0.5,0.1), ylim = seq(-0.5,0.5,0.1), type = 2) +
  plotContour(cc_2, xlim = seq(-0.5,0.5,0.1), ylim = seq(-0.5,0.5,0.1), type = 2) +
  plotContour(cc_3, xlim = seq(-0.5,0.5,0.1), ylim = seq(-0.5,0.5,0.1), type = 2) +

plotContour(cc_1_zoom, xlim = seq(0,0.5,0.1), ylim = seq(0,0.5,0.1), type = 2) +
  plotContour(cc_2_zoom, xlim = seq(0,0.5,0.1), ylim = seq(0,0.5,0.1), type = 2) +
  plotContour(cc_3_zoom, xlim = seq(0,0.5,0.1), ylim = seq(0,0.5,0.1), type = 2)



# Código antigo e optim ---------------------------------------------------


cc <- curve3d(f.vu(c(x,y), theta = 1/9), xlim = c(-0.1,0.5), ylim = c(-0.1, 0.5),
              sys3d = "image", xlab = "v", ylab = "u")
cc



# Quando faço o gráfico com intervalos menores, vejo que quando theta é pequeno,
# valores pequenos de 'v' e 'u' me dão os menores tamanhos amostrais.
# Já quando theta theta cresce, o oposto ocorre, 'v' e 'u' grandes me dão os menores tamanhos amostrais.
# Fazer o gráfico para isso:

# Todo o intervalo
#par(mfrow = c(1,3))
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


## optimx
oo <- optimx(par = c(0.1, 0.5), fn = f.vu,  lower = c(0.0001, 0.0002), upper = c(0.49, 0.5),  mu = 1, theta = 1/9)
oo
c(round(oo$p1,2), round(oo$p2,2))
f.vu(c(round(oo$p1,2), round(oo$p2,2)), theta = 1/9)

oo <- optimx(par = c(0.17, 0.30), fn = f.vu,  lower = c(0.0001, 0.0002), upper = c(0.49, 0.5),  mu = 1, theta = 1)
oo
c(round(oo$p1,2), round(oo$p2,2))
f.vu(c(oo$p1, oo$p2), theta = 1)

oo <- optimx(par = c(0.45, 0.49), fn = f.vu,  lower = c(0.0001, 0.0002), upper = c(0.49, 0.5),  mu = 1, theta = 9)
oo
f.vu(c(oo$p1, oo$p2), theta = 9)

f.vu(c(0.499, 0.5), theta = 9)

oo = optim(par = c(0.4, 0.45), f.vu, theta = 9,
           lower = c(0.1,0.2),
           upper = c(0.49, 0.5), method="L-BFGS-B")
oo





# New function for the plots


pf <- function(df, zval, pal, lims){
  ggplot(df, aes(x = x, y = y, fill = zval)) +
    geom_tile() +
    scale_fill_gradientn(colors = pal,
                         limits = lims,
                         name = 'n',
                         oob = scales::squish,
                         na.value = "white",
                         guide = guide_colorbar(
                           order = 2,
                           frame.colour = "white",
                           ticks.colour = "white"
                         )
    ) +
    labs(x = "v", y = "u") +
    theme_classic(base_size = 16) +
    #theme(legend.position = "none")+
    scale_x_continuous(breaks = seq(-0.5,0.5, 0.1)) +
    scale_y_continuous(breaks = seq(-0.5,0.5, 0.1))
}



cc_1 <- curve3d(f.vu(c(x,y), theta = 1/9), xlim = c(-0.5,0.5), ylim = c(-0.5, 0.5),
                sys3d = "image", xlab = "v", ylab = "u")
x <- cc_1$x
y <- cc_1$y
z <- as.vector(cc_1$z)
df <- tibble::tibble(tidyr::expand_grid(x,y), z)
range(z[!is.na(df$z)& z!=Inf])

plotContour(contour = cc_1, xlim = 0, ylim = 2500, type = 3)


# set colour palette
pal <- colorRampPalette(c("grey90","black"))
pal <- pal(100)

p1 <- pf(df = df, zval = df$z, pal = pal, lims = c(0,2500))
p1

cc_1 <- curve3d(f.vu(c(x,y), theta = 1), xlim = c(-0.5,0.5), ylim = c(-0.5, 0.5),
                sys3d = "image", xlab = "v", ylab = "u")
x <- cc_1$x
y <- cc_1$y
z <- as.vector(cc_1$z)
df <- tibble::tibble(tidyr::expand_grid(x,y), z)
range(z[!is.na(df$z)& z!=Inf])

# set colour palette
pal <- colorRampPalette(c("grey90","black"))
pal <- pal(100)

p2 <- pf(df = df, zval = df$z, pal = pal, lims = c(0,2500))
p2


cc_1 <- curve3d(f.vu(c(x,y), theta = 9), xlim = c(-0.5,0.5), ylim = c(-0.5, 0.5),
                sys3d = "image", xlab = "v", ylab = "u")
x <- cc_1$x
y <- cc_1$y
z <- as.vector(cc_1$z)
df <- tibble::tibble(tidyr::expand_grid(x,y), z)
range(z[!is.na(df$z)& z!=Inf])

# set colour palette
pal <- colorRampPalette(c("grey90","black"))
pal <- pal(100)

p3 <- pf(df = df, zval = df$z, pal = pal, lims = c(0,2500))
p3

p1+p2+p3


