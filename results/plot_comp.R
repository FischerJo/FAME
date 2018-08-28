
require(ggplot2)
require(ggthemes)
require(scales)

names <- c("Bismark", "BSMAP", "BSSeeker2", "FAME 500", "FAME 2000", "FAME 5000")
rmse <- c(0.1214, 0.1418, 0.1265, 0.1236, 0.1119, 0.107)
runtimes <- c(110908 + 34524, 8023 + 33849, 755520 + (33*60 + 39)*60 + 33, (2*60 + 51)*60 + 28, (3*60 + 23)*60 + 19, (5*60 + 58)*60 + 57)
runtimes <- runtimes / (60*60)

data <- data.frame(Method=names, RMSE=rmse, Runtimes=runtimes)

data$Method <- factor(data$Method, levels = c("Bismark", "BSMAP", "BSSeeker2", "FAME 500", "FAME 2000", "FAME 5000"))

pdf("runtime_and_rmse.pdf")

gg <- ggplot(data, aes(x=Method, y=RMSE))
gg <- gg + geom_bar(width=0.6, stat='identity') 
gg <- gg + theme_minimal()
gg <- gg + scale_y_continuous(limits = c(0, 0.2))
gg <- gg + theme(plot.margin = unit(c(2,2,2,2),"cm"))
gg <- gg + theme(axis.title.x = element_text(margin=margin(20,0,0,0)),
                 axis.title.y = element_text(margin=margin(0,20,0,0)),
                 plot.title = element_text(margin=margin(0,0,40,0)),
		 axis.text.x = element_text(angle = 45, hjust = 1)
                )
gg <- gg + labs(title="RMSE of methylation rates on chromosome 22 data set")

gg

breaks <- 10^(-10:10)
minor_breaks <- rep(c(1:9), 21)*(10^rep(-10:10, each=9))
gg <- ggplot(data, aes(x=Method, y=Runtimes))
gg <- gg + geom_bar(width=0.4, stat='identity') 
gg <- gg + theme_minimal()
gg <- gg + scale_y_log10(breaks = breaks,
   labels = sapply(log10(breaks), scales::math_format(10^.x)), minor_breaks = minor_breaks)
gg <- gg + geom_hline(yintercept=24,color="red")
gg <- gg + annotate("text", "Bismark", 27, label = "24h                   ", color="red", size=3)
gg <- gg + theme(plot.margin = unit(c(2,2,2,2),"cm"))
gg <- gg + theme(axis.title.x = element_text(margin=margin(20,0,0,0)),
                 axis.title.y = element_text(margin=margin(0,20,0,0)),
                 plot.title = element_text(margin=margin(0,0,40,0)),
		 axis.text.x = element_text(angle = 45, hjust = 1)
                )
#gg <- gg + labs(title="Runtimes on whole genome data set",y="Runtime in hours (log scale)")
gg <- gg + labs(title="Runtimes on synthetic data set",y="Runtime in hours (log scale)")

gg

dev.off()
