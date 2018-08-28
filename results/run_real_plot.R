
require(ggplot2)
require(ggthemes)
require(scales)

names <- c("Bismark", "BSMAP", "Metal STL", "Metal sparsepp", "Metal dense")
runtimes <- c(230620 + 0, 98108 + 44324, (7*60+11)*60+21, (6*60+23)*60+13, (4*60+7)*60+16)
runtimes <- runtimes / 60

data <- data.frame(Method=names, Runtimes=runtimes)

data$Method <- factor(data$Method, levels = c("Bismark", "BSMAP", "Metal STL", "Metal sparsepp", "Metal dense"))

pdf("runtime_real.pdf")


gg <- ggplot(data, aes(x=Method, y=Runtimes))
gg <- gg + geom_bar(width=0.6, stat='identity') 
gg <- gg + theme_minimal()
#gg <- gg + geom_text(aes(label=Method), size = 3, nudge_y=0.01, hjust=0, parse="True")
#gg <- gg + scale_color_manual(values=colos)
gg <- gg + geom_hline(yintercept=(60*24),color="red")
gg <- gg + annotate("text", "Bismark", 60*26, label = "24h                       ", color="red", size=3)
gg <- gg + scale_y_continuous(labels = comma)
gg <- gg + theme(plot.margin = unit(c(2,2,2,2),"cm"))
gg <- gg + theme(axis.title.x = element_text(margin=margin(20,0,0,0)),
                 axis.title.y = element_text(margin=margin(0,20,0,0)),
                 plot.title = element_text(margin=margin(0,0,40,0)),
		 axis.text.x = element_text(angle = 45, hjust = 1)
                )
gg <- gg + labs(title="Runtimes on ENCODE data set ENCFF414AAO",y="Runtime in minutes")

gg

dev.off()
