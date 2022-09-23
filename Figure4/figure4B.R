
########################
# R script for figure 4 B: Density and CDF curves of a random sample
########################


######################
# load libraries
library(plyr)

######################
# Two random samples
Sample1 = c(-0.4944729, 0.4747895, 0.5487083, -0.6990714, 1.1185330, 1.0494700, -0.6771520, 0.8839922, 0.9312416, 1.6370360)
Sample2 = c(1.65800000, 1.65267300, -0.06678641, 1.72421500, -0.37933150, -1.09030800,  1.16830900, -0.39002530,  2.01759800,  1.44725000)

######################
# Melt the data
df1 = data.frame(Sample1,Sample2)
df2 = melt(df1)

######################
# Calculate median of each sample
mu <- ddply(df2, "variable", summarise, grp.mean=median(value))
head(mu)

######################
# Setup ggplot parameters
Y.text <- element_text( color = "black", size = 13, angle = 0)
X.text <- element_text( color = "black",  size = 13, angle = 0) 
a.text = element_text(size=16, color = "black", angle = 0)
t.text = element_text(size=16, color = "black")  



######################
# Fig 1
# (a) plot density curves, (b) plot CDF curves

a = ggplot(df2, aes(x=value, color=variable)) +
    geom_density(size=.7) + 
    geom_vline(data=mu, aes(xintercept=grp.mean, color=variable), linetype="dashed") +
    scale_color_manual(values=c("black", "#E69F00")) +

    # set up the y-axis labels
    scale_y_continuous(limits=c(0,.6), breaks=c(0, 0.1, 0.2,
                                                   0.3, 0.4, 0.5, 0.6), expand = c(0,0)) +
    theme_bw() + 
    theme(legend.position = 'none') +
    theme(axis.text.y = Y.text, axis.text.x = X.text, axis.title = a.text, title = t.text, plot.title = element_text(size=17) ) +
    labs(x = NULL, y = "Density curves", title = NULL)


b = ggplot(df2, aes(x=value, color=variable)) + stat_ecdf()+
    theme_bw() + 
    scale_color_manual(values=c("black", "#E69F00")) +
    theme(legend.position = 'none') +
    scale_y_continuous(limits=c(0,1), breaks=c(0, 0.2, 0.4,
                                              0.6, 0.8, 1), expand = c(0,0)) +
    labs(title=NULL, y = "Cumulative distribution", x = "Subsample feature data") +
    theme(axis.text.y = Y.text, axis.text.x = X.text, axis.title = a.text, title = t.text, plot.title = element_text(size=15) ) 

  
######################
# Plot the two figures together
  grid.arrange(a,b,nrow = 2)


  
  
  
  ######################
  # Fig 2
  # Show vertical KS distance on the CDF

  sample1 = Sample1
  sample2 = Sample2

  group <- c(rep("sample1", length(sample1)), rep("sample2", length(sample2)))
  dat <- data.frame(KSD = c(sample1,sample2), group = group)

  ######################
  # create ECDF of data
  cdf1 <- ecdf(sample1) 
  cdf2 <- ecdf(sample2) 

  ######################
  # find min and max statistics to draw line between points of greatest distance
  minMax <- seq(min(sample1, sample2), max(sample1, sample2), length.out=length(sample1)) 
  x0 <- minMax[which( abs(cdf1(minMax) - cdf2(minMax)) == max(abs(cdf1(minMax) - cdf2(minMax))) )] 
  y0 <- cdf1(x0) 
  y1 <- cdf2(x0) 

  ######################
  # Show vertical KS distance on the CDF
  
  ggplot(dat, aes(x = KSD, group = group, color = group))+
  stat_ecdf(size=1) +
  theme_bw(base_size = 28) +
  theme(legend.position ="right") +
  # xlab("Sample") +
  # ylab("ECDF") +

  geom_segment(aes(x = x0[1], y = y0[1], xend = x0[1], yend = y1[1]),
               linetype = "dashed", color = "red") +
  geom_point(aes(x = x0[1] , y= y0[1]), color="red", size=4) +
  geom_point(aes(x = x0[1] , y= y1[1]), color="red", size=4) +
  # ggtitle("KS Test: \nMax discrepency between ECDFs") +
  scale_color_manual(values=c("black", "#E69F00")) +
  theme(legend.title=element_blank()) +
  theme(legend.position = 'none') +
  labs(y = NULL, x=NULL) +
  theme(axis.text.y = Y.text, axis.text.x = X.text, axis.title = a.text, title = t.text, plot.title = element_text(size=20) ) 






  
  
  ######################
  # Fig 3
  ###################### 
  # non ggplot example
  
  # CDF 1
  plot(cdf1, xlim=range(sample1, sample2), verticals=TRUE, do.points=FALSE, col="black", lwd = 3, main = NULL,
     col.main="blue", col.lab="red", col.sub="red") 
  # CDF 2
  plot(cdf2, verticals=TRUE, do.points=FALSE, col="#E69F00", lwd = 3, add=TRUE) 
  # Show KS distance with vertical red line
  points(c(x0, x0), c(y0, y1), pch=16, col="red") 
  segments(x0, y0, x0, y1, col="red", lty="dotted") 

  
  
  
  
  ######################
  # Fig 4
  ###################### 
  # Using line segments to shade in the area between two ECDFs

  ecdf0 <- cdf1 
  ecdf1 <- cdf2

  plot(ecdf0, xlim=range(sample1, sample2), verticals=TRUE, do.points=FALSE, col="black", lwd = 3, main = NULL,
     col.main="blue", col.lab="red", col.sub="red") 
  plot(ecdf1, verticals=TRUE, do.points=FALSE, col="#E69F00", lwd = 3, add=TRUE) 

  # ----- this will plot lots of segments (2001 of them?)
  tmpx <- seq(min(sample1, sample2), max(sample1, sample2), len=2001) 
  p0 <- ecdf0(tmpx) # ----- y upper and y lower
  p1 <- ecdf1(tmpx) # ----- y upper and y lower
  segments(x0=tmpx, y0=p0, y1=p1, col=adjustcolor(ifelse(p0<p1, "#E69F00", "grey"), alpha.f=0.1)) 

  # Show KS distance with vertical red line
  points(c(x0, x0), c(y0, y1), pch=16, col="red") 
  segments(x0, y0, x0, y1, col="red", lty="dotted") 







