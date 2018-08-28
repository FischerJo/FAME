
require(ggplot2)
require(ggthemes)
require(scales)

names <- c("Bismark", "BSMAP", "FAME")
runtimes <- c(230620 + 50013, 98108 + 44324, 25881)
runtimes <- runtimes / 60

data <- data.frame(Method=names, Runtimes=runtimes)

data$Method <- factor(data$Method, levels = c("Bismark", "BSMAP", "FAME"))

gg <- ggplot(data, aes(x=Method, y=Runtimes))
gg <- gg + geom_bar(width=0.4, stat='identity') 
gg <- gg + theme_minimal()
gg <- gg + geom_hline(yintercept=(60*24),color="red")
gg <- gg + annotate("text", "Bismark", 60*26, label = "24h                                     ", color="red", size=3)
gg <- gg + scale_y_continuous(labels = comma)
gg <- gg + theme(plot.margin = unit(c(2,2,2,2),"cm"))
gg <- gg + theme(axis.title.x = element_text(margin=margin(20,0,0,0)),
                 axis.title.y = element_text(margin=margin(0,20,0,0)),
                 plot.title = element_text(margin=margin(0,0,40,0)),
		 axis.text.x = element_text(angle = 45, hjust = 1)
                )
gg <- gg + labs(title="Runtimes on ENCODE data set",y="Runtime in minutes")

names2 <- c("Bismark", "BSMAP", "FAME STL", "FAME sparsepp", "FAME dense")
runtimes2 <- c(230620 + 50013, 98108 + 44324, 25881, 22993, 14836)
runtimes2 <- runtimes2 / 60

data2 <- data.frame(Method=names2, Runtimes=runtimes2)

data2$Method <- factor(data2$Method, levels = c("Bismark", "BSMAP", "FAME STL", "FAME sparsepp", "FAME dense"))

gg2 <- ggplot(data2, aes(x=Method, y=Runtimes))
gg2 <- gg2 + geom_bar(width=0.4, stat='identity') 
gg2 <- gg2 + theme_minimal()
gg2 <- gg2 + geom_hline(yintercept=(60*24),color="red")
gg2 <- gg2 + annotate("text", "Bismark", 60*26, label = "24h                    ", color="red", size=3)
gg2 <- gg2 + scale_y_continuous(labels = comma)
gg2 <- gg2 + theme(plot.margin = unit(c(2,2,2,2),"cm"))
gg2 <- gg2 + theme(axis.title.x = element_text(margin=margin(20,0,0,0)),
                 axis.title.y = element_text(margin=margin(0,20,0,0)),
                 plot.title = element_text(margin=margin(0,0,40,0)),
		 axis.text.x = element_text(angle = 45, hjust = 1)
                )
gg2 <- gg2 + labs(title="Runtimes on ENCODE data set",y="Runtime in minutes")

names3 <- c("Bismark", "BSMAP", "FAME")
runtimes3 <- c(230620 + 50013, 98108 + 44324, 14836)
runtimes3 <- runtimes3 / (60 * 60)

data3 <- data.frame(Method=names3, Runtimes=runtimes3)

data3$Method <- factor(data3$Method, levels = c("Bismark", "BSMAP", "FAME"))
breaks <- 10^(-10:10)
minor_breaks <- rep(c(1:9), 21)*(10^rep(-10:10, each=9))
gg3 <- ggplot(data3, aes(x=Method, y=Runtimes))
gg3 <- gg3 + geom_bar(width=0.4, stat='identity') 
gg3 <- gg3 + theme_minimal()
gg3 <- gg3 + geom_hline(yintercept=(24),color="red")
gg3 <- gg3 + annotate("text", "Bismark", 27, label = "24h                                 ", color="red", size=3)
gg3 <- gg3 + scale_y_log10(breaks = breaks,
   labels = sapply(log10(breaks), scales::math_format(10^.x)), minor_breaks = minor_breaks)
#gg3 <- gg3 + coord_trans(y="log10")
#gg3 <- gg3 + scale_y_continuous(trans = log10_trans(),
  #                   breaks = trans_breaks("log10", function(x) 10^x),
 #                    labels = trans_format("log10", math_format(10^.x)))
gg3 <- gg3 + theme(plot.margin = unit(c(2,2,2,2),"cm"))
gg3 <- gg3 + theme(axis.title.x = element_text(margin=margin(20,0,0,0)),
                 axis.title.y = element_text(margin=margin(0,20,0,0)),
                 plot.title = element_text(margin=margin(0,0,40,0)),
		 axis.text.x = element_text(angle = 45, hjust = 1)
                )
gg3 <- gg3 + labs(title="Runtimes on ENCODE data set*",y="Runtime in hours (log scale)")

pdf("runtime_real.pdf")

gg
gg2
gg3

dev.off()
