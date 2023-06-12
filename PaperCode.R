# Code to generate the paper results

library(qExp)

# Figure 1 ----------------------------------------------------------------

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


# Figure 3 ----------------------------------------------------------------


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


