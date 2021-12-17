library(gridExtra)
library(reshape2)
library(ggplot2)
theme_set(theme_bw())
lty <- c("11", "18", "1f", "81", "88", "8f", "f1", "f8", "ff")

# values of v>0
v1 <- seq(0.1,0.5,0.1)
theta <- c(0:100)
df.v1 <- outer(theta, v1,g.theta)
df.v1 <- data.frame(df.v1, x = rep(0:100,length(v1)))

p.v1 <- ggplot(melt(df.v1,id.vars='x'), aes(x,value,group=variable))+
  geom_line(aes(linetype = variable), size = 0.4) +
  scale_linetype_manual(values = rev(lty[1:5]),labels=v1, name="v")+
  labs(x=expression(theta), y = bquote(paste("g(",theta,")")))+
  #theme(legend.position="top")
  theme(legend.position=c(0.85,0.2))


p.v1


# values of v<0
v2 <- seq(-0.5,-0.1,0.1)
theta <- c(0:100)
df.v2 <- outer(theta, v2,g.theta)
df.v2 <- data.frame(df.v2, x = rep(0:100,length(v2)))

p.v2 <- ggplot(melt(df.v2,id.vars='x'), aes(x,value,group=variable))+
  geom_line(aes(linetype = variable), size = .5) +
  scale_linetype_manual(values = lty[1:5],labels=v2, name="v")+
  labs(x=expression(theta), y = "")+
  theme(legend.position=c(0.85,0.8))


p.v2


grid.arrange(p.v1,p.v2, nrow = 1 )


# v = -0.5:0.5

v3 <- c(v1,v2)
theta <- c(0:100)
df.v3 <- outer(theta, v3,g.theta)
df.v3 <- data.frame(df.v3, x = rep(0:100,length(v3)))

p.v3 <- ggplot(melt(df.v3,id.vars='x'), aes(x,value,group=variable))+
  geom_line(aes(linetype = variable), size = 1) +
  scale_linetype_discrete(labels=v3, name="v")+
  labs(x=expression(theta), y = bquote(paste("g(",theta,")")))


p.v3


# Colour
v4 <- c(v1, v2)
theta <- c(0:100)
df.v4 <- outer(theta, v4,g.theta)
df.v4 <- data.frame(df.v4, x = rep(0:100,length(v4)))
df.new <- reshape2::melt(df.v4, )
library(RColorBrewer)
library(colorspace)
# display.brewer.all(colorblindFriendly = TRUE)

p.v4 <- ggplot(melt(df.v4,id.vars='x'), aes(x,value,colour=variable))+
  geom_line(aes(color = variable), size = 0.25) +
  labs(x=expression(theta), y = bquote(paste("g(",theta,")")), color = "v") +
  scale_color_discrete_sequential(palette = "Plasma", n = 10, labels = v4)

p.v4



#

v <- c(-0.1,-0.2, -0.3, -0.4, 0.1, 0.2, 0.3, 0.4)
u <- c(-0.1,-0.2, -0.3, -0.4, 0.1, 0.2, 0.3, 0.4)
x <- expand_grid(v, u) %>% filter(v!=u)
theta <- c(0:100)
df <- expand_grid(x, theta)

df1 <- df %>% filter(v ==0.4&u==0.1)

mapply(g.theta, v = df1$v, u = df1$u, theta = df1$theta)  %>% as.data.frame() %>%
  `colnames<-`("y") %>% mutate(x = 1:n()) %>%
  ggplot() + geom_line(aes(x = x, y = y)) +theme_bw()


