require(ggplot2)
require(ggthemes)
require(scales)

get_secs <- function(h,m,s)
{
	return((h*60 + m)*60 + s)
}

get_hours <- function(h,m,s)
{
	return(h + m/60 + s/(60*60))
}

data <- data.frame(
names = c("BSMAP", "BratNova", "Bismark", "FAME[def]", "FAME[acc]"),
#runtime = c(get_secs(7+10,22+49,30+55), get_secs(25+3,42+3,24+56), get_secs(35+13,3+9,25+53), get_secs(4,28,32)),
runtime = c(get_hours(1,59+54,24+58), get_hours(1+1,58+26,30+18), get_hours(1,30+14,34+47), get_hours(0,15,34), get_hours(0,27,12)),
unmapped = c(583528,44800,12371,13949,9713),
rmse = c(0.31,0.07,0.06,0.06,0.048),
corr = c(0.3,0.97,0.96,0.97,0.98)
)

# Todo log scale?

gg <- ggplot(data, aes(x=runtime, y=rmse))
gg <- gg + geom_point(aes(colour=names,size=unmapped), alpha=0.8,stroke=0.1)
gg <- gg + geom_text(aes(label=names),size=3,nudge_y=0.01,nudge_x=0.2,parse=T)
gg <- gg + ylim(0,0.35) + xlim(0,4)
gg <- gg + geom_vline(xintercept=(24),color="red", lwd=0.2)
#gg <- gg + annotate("text", 22, 0, label = "24h", color="red", size=3)
#gg <- gg + coord_trans(x="log2")
gg <- gg + theme_tufte()
gg <- gg + theme(plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
gg <- gg + theme(plot.background = element_blank(),
                 panel.border = element_blank(),
                 axis.title.x = element_text(margin=margin(20,0,0,0)),
                 axis.title.y = element_text(margin=margin(0,20,0,0)),
                 plot.title = element_text(margin=margin(0,0,40,0))
                 )
gg <- gg + scale_color_manual(values=c("FAME[def]"="firebrick","FAME[acc]"="firebrick", "Bismark"="midnightblue", "BratNova"="tan4","BSMAP"="goldenrod1"), 
                       		name="Method",
                       		breaks=c("FAME[def]", "FAME[acc]", "Bismark", "BratNova", "BSMAP"),
                       		labels=c("FAME[def]", "FAME[acc]", "Bismark", "BratNova", "BSMAP")
)
gg <- gg + guides(color=FALSE)
gg <- gg + scale_size_continuous(name="#Unmapped\nCpGs", breaks=c(15000,50000,100000,250000,500000), labels = comma)

gg <- gg + labs(x="Runtime (hours)",
	        y="RMSE",
		title="Performance of methylation calling on synthetic data")

gg2 <- ggplot(data, aes(x=runtime, y=corr))
gg2 <- gg2 + geom_point(aes(colour=names,size=unmapped), alpha=0.8,stroke=0.1)
gg2 <- gg2 + geom_text(aes(label=names),size=3,nudge_y=-0.02,nudge_x=0.2,parse=T)
gg2 <- gg2 + ylim(0,1) + xlim(0,4)
gg2 <- gg2 + geom_vline(xintercept=(24),color="red", lwd=0.2)
#gg2 <- gg2 + annotate("text", 22, 0, label = "24h", color="red", size=3)
#gg2 <- gg2 + coord_trans(x="log2")
gg2 <- gg2 + theme_tufte()
gg2 <- gg2 + theme(plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
gg2 <- gg2 + theme(plot.background = element_blank(),
                 panel.border = element_blank(),
                 axis.title.x = element_text(margin=margin(20,0,0,0)),
                 axis.title.y = element_text(margin=margin(0,20,0,0)),
                 plot.title = element_text(margin=margin(0,0,40,0))
                 )
gg2 <- gg2 + scale_color_manual(values=c("FAME[def]"="firebrick","FAME[acc]"="firebrick", "Bismark"="midnightblue", "BratNova"="tan4","BSMAP"="goldenrod1"), 
                       		name="Method",
                       		breaks=c("FAME[def]", "FAME[acc]", "Bismark", "BratNova", "BSMAP"),
                       		labels=c("FAME[def]", "FAME[acc]", "Bismark", "BratNova", "BSMAP")
)
gg2 <- gg2 + guides(color=FALSE)
gg2 <- gg2 + scale_size_continuous(name="#Unmapped\nCpGs", breaks=c(15000,50000,100000,250000,500000), labels = comma)

gg2 <- gg2 + labs(x="Runtime (hours)",
	        y="Spearman Correlation Coefficient",
		title="Performance of methylation calling on synthetic data")

pdf("summary_plot.pdf")

print(gg)
print(gg2)

dev.off()

