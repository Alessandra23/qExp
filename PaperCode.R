# Code to generate the paper results

library(qExp)

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
  theme_classic(base_size = 20)+
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
  theme_classic(base_size = 20)+
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
  theme_classic(base_size = 20)

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
  ggplot() + geom_line(aes(x = theta, y = g, linetype = vu)) + theme_classic(base_size = 20)+
  labs(x=expression(theta), y = bquote(paste("g(",theta,")")), linetype= "v,u")

# Fig 2(b)
df %>% filter(v>0&u>0&v>u) %>%
  ggplot() + geom_line(aes(x = theta, y = g, linetype = vu)) + theme_classic(base_size = 20) +
  labs(x=expression(theta), y = " ",  linetype= "v,u")


# Figure 3 ----------------------------------------------------------------

# The function f.vu is in the file qExp/R/gThetaFunction and the plotContour function is in the file
# qExp/R/plots

# theta = 1/9
cc_1 <- curve3d(f.vu(c(x,y), theta = 1/9), xlim = c(-0.5,0.5), ylim = c(-0.5, 0.5),
                sys3d = "image", xlab = "v", ylab = "u")
# Fig 3(a)
plotContour(contour = cc_1, xlim = 0, ylim = 2500, type = 3)

cc_1_zoom <- curve3d(f.vu(c(x,y), theta = 1/9), xlim = c(0,0.5), ylim = c(0, 0.5),
                     sys3d = "image", xlab = "v", ylab = "u")
# Fig 3(d)
plotContour(cc_1_zoom, xlim = seq(0,0.5,0.1), ylim = seq(0,0.5,0.1))


# theta = 1
cc_2 <- curve3d(f.vu(c(x,y), theta = 1), xlim = c(-0.5,0.5), ylim = c(-0.5, 0.5),
                sys3d = "image", xlab = "v", ylab = "u")
# Fig 3(b)
plotContour(contour = cc_2, xlim = 0, ylim = 2500, type = 3)

cc_2_zoom <- curve3d(f.vu(c(x,y), theta = 1), xlim = c(0,0.5), ylim = c(0, 0.5),
                     sys3d = "image", xlab = "v", ylab = "u")
# Fig 3(e)
plotContour(cc_2_zoom, xlim = seq(0,0.5,0.1), ylim = seq(0,0.5,0.1))


# theta = 9
cc_3 <- curve3d(f.vu(c(x,y), theta = 9), xlim = c(-0.5,0.5), ylim = c(-0.5, 0.5),
                sys3d = "image", xlab = "v", ylab = "u")
# Fig 3(c)
plotContour(contour = cc_3, xlim = 0, ylim = 2500, type = 3)

cc_3_zoom <- curve3d(f.vu(c(x,y), theta = 9), xlim = c(0,0.5), ylim = c(0, 0.5),
                     sys3d = "image", xlab = "v", ylab = "u")
# Fig 3(f)
plotContour(cc_3_zoom, xlim = seq(0,0.5,0.1), ylim = seq(0,0.5,0.1))


# Figure 4 ----------------------------------------------------------------



# Figure 5 ----------------------------------------------------------------


# Figure 6 ----------------------------------------------------------------


# Figure 7 ----------------------------------------------------------------


# Figure 8 ----------------------------------------------------------------


# Figure 9 ----------------------------------------------------------------


# Figure 10 ---------------------------------------------------------------


# Figure 11 ---------------------------------------------------------------


# Figure 12 ---------------------------------------------------------------


# Figure 13 ---------------------------------------------------------------


# Figure 14 ---------------------------------------------------------------


# Figure 15 ---------------------------------------------------------------

# Figure 16 ---------------------------------------------------------------

# Figure 17 ---------------------------------------------------------------


