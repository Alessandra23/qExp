
# Plots qexp kappa and q

# setting kappa

x <- seq(0,4,length=1000)
v1.qk.k <- dtsal(x,shape = 0.1,scale = 1)
v2.qk.k <- dtsal(x,shape = 1,scale = 1)
v3.qk.k <- dtsal(x,shape = 2,scale = 1)
v4.qk.k <- dtsal(x,shape = 5,scale = 1)

df.qk.k <- data.frame(x = x, v1 = v1.qk.k, v2 = v2.qk.k, v3 = v3.qk.k, v4 = v4.qk.k)


legenda.k <- c(expression(paste(q, " = 1, ", kappa," = 0.1" ),
                          paste(q, " = 1, ", kappa," = 1 "),
                          paste(q, " = 1, ", kappa," = 2 "),
                          paste(q, " = 1, ", kappa," = 5 ")))

p.qk.k <- ggplot(melt(df.qk.k,id.vars='x'), aes(x,value,group=variable))+
  geom_line(aes(linetype = variable)) + ylab("f(x)")+
  scale_linetype_discrete(labels=legenda.k, name=" ")+
  theme(legend.position=c(0.8,0.8))

p.qk.k

####



# setting q

v1.qk.q <- dtsal(x,shape = 1,scale = 1.001)
v2.qk.q <- dtsal(x,shape = 1,scale = 1.2)
v3.qk.q <- dtsal(x,shape = 1,scale = 1.5)

df.qk.q <- data.frame(x = x, v1 = v1.qk.q, v2 = v2.qk.q, v3 = v3.qk.q)


legenda.q <- c(expression(paste(q, " = 1, ", kappa," = 1 " ),
                          paste(q, " = 1.2, ", kappa," = 1 "),
                          paste(q, " = 1.5, ", kappa," = 1 ")))

p.qk.q <- ggplot(melt(df.qk.q,id.vars='x'), aes(x,value,group=variable))+
  geom_line(aes(linetype = variable)) + ylab("f(x)")+
  scale_linetype_discrete(labels=legenda.q, name=" ")+
  theme(legend.position=c(0.8,0.8))

p.qk.q




grid.arrange(p.qk.k,p.qk.q, nrow = 1)

#####

# Plots qexp mu and theta

theta.1 = 0.5
theta.2 = 1
theta.3 = 6
theta.4 = 100


# setting mu = 1

x <- seq(0,5,length=1000)
v1.mt.1 <- dtsal(x,shape = theta.1+1, scale = 1*theta.1)
v1.mt.2 <- dtsal(x,shape = theta.2+1, scale = 1*theta.2)
v1.mt.3 <- dtsal(x,shape = theta.3+1, scale = 1*theta.3)
v1.mt.4 <- dtsal(x,shape = theta.4+1, scale = 1*theta.4)

df.mt.1 <- data.frame(x = x, v1 = v1.mt.1, v2 = v1.mt.2, v3 = v1.mt.3, v4 = v1.mt.4)


legenda.1 <- c(bquote(paste(mu, " = 1.0, ", theta," = 0.5 " ),
                          paste(mu, " = 1.0, ", theta," = 1 "),
                          paste(mu, " = 1.0, ", theta," = 6 "),
                          paste(mu, " = 1.0, ", theta," = 100 ")))


p.mt.1 <- ggplot(melt(df.mt.1,id.vars='x'), aes(x,value,group=variable))+
  geom_line(aes(linetype = variable)) + ylab("f(x)")+
  scale_linetype_discrete(labels=c(bquote(paste(mu, " = 1.0, ", theta," = 0.5 ")),
                                   bquote(paste(mu, " = 1.0, ", theta," = 1 ")),
                                   bquote(paste(mu, " = 1.0, ", theta," = 6 ")),
                                   bquote(paste(mu, " = 1.0, ", theta," = 100 "))),
                          name=" ")+
  theme_classic(base_size = 20)+
  theme(legend.position=c(0.8,0.8), legend.box.just = "left")

p.mt.1


# setting mu = 10

x <- seq(0,50,length=1000)
v2.mt.1 <- dtsal(x,shape = theta.1+1, scale = 10*theta.1)
v2.mt.2 <- dtsal(x,shape = theta.2+1, scale = 10*theta.2)
v2.mt.3 <- dtsal(x,shape = theta.3+1, scale = 10*theta.3)
v2.mt.4 <- dtsal(x,shape = theta.4+1, scale = 10*theta.4)

df.mt.2 <- data.frame(x = x, v1 = v2.mt.1, v2 = v2.mt.2, v3 = v2.mt.3, v4 = v2.mt.4)


legenda.2 <- c(expression(paste(mu, " = 10, ", theta," = 0.5 " ),
                          paste(mu, " = 10, ", theta," = 1 "),
                          paste(mu, " = 10, ", theta," = 6 "),
                          paste(mu, " = 10, ", theta," = 100 ")))

p.mt.2 <- ggplot(melt(df.mt.2,id.vars='x'), aes(x,value,group=variable))+
  geom_line(aes(linetype = variable)) + ylab("f(x)")+
  scale_linetype_discrete(labels=legenda.2, name=" ")+
  theme(legend.position=c(0.8,0.8), legend.box.just = "left")

p.mt.2


# setting mu = 100

x <- seq(0,500,length=1000)
v3.mt.1 <- dtsal(x,shape = theta.1+1, scale = 100*theta.1)
v3.mt.2 <- dtsal(x,shape = theta.2+1, scale = 100*theta.2)
v3.mt.3 <- dtsal(x,shape = theta.3+1, scale = 100*theta.3)
v3.mt.4 <- dtsal(x,shape = theta.4+1, scale = 100*theta.4)

df.mt.3 <- data.frame(x = x, v1 = v3.mt.1, v2 = v3.mt.2, v3 = v3.mt.3, v4 = v3.mt.4)


legenda.3 <- c(expression(paste(mu, " = 100, ", theta," = 0.5 " ),
                          paste(mu, " = 100, ", theta," = 1 "),
                          paste(mu, " = 100, ", theta," = 6 "),
                          paste(mu, " = 100, ", theta," = 100 ")))

p.mt.3 <- ggplot(melt(df.mt.3,id.vars='x'), aes(x,value,group=variable))+
  geom_line(aes(linetype = variable)) + ylab("f(x)")+
  scale_linetype_discrete(labels=legenda.3, name=" ")+
  theme(legend.position=c(0.8,0.8), legend.box.just = "left")

p.mt.3


grid.arrange(p.mt.1,p.mt.2, p.mt.3, nrow = 1)



## Paper - 4 cases



theta.1 = 0.5
theta.2 = 9


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




# Comparison of the tails of the Exponential and qExponential distributions (theta = 1 and mu = 1)

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

grid.arrange(p.comp.1,p.comp.2, nrow = 1)

p.comp.1 + annotation_custom(ggplotGrob(p.comp.2), xmin = 6, xmax = 15,
                             ymin = 0.01, ymax = 1)
