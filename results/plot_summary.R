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
get_minutes <- function(h,m,s)
{
	return(h*60 + m + s/60)
}

data <- data.frame(
names = c("BSMAP", "BratNova", "Bismark", "FAME", "Segemehl"),
#runtime = c(get_secs(7+10,22+49,30+55), get_secs(25+3,42+3,24+56), get_secs(35+13,3+9,25+53), get_secs(4,28,32)),
runtime = c(get_minutes(1,59+54,24+58), get_minutes(1+1,58+26,30+18), get_minutes(1+2,46+1,27+26), get_minutes(0,13,35), get_minutes(2,5+8+31,45+50+28)),
unmapped = c(583528,44800,12371,12408,57),
rmse = c(0.61,0.16,0.1,0.099,0.04),
corr = c(-0.34,0.85,0.93,0.94,0.99)
)

# Todo log scale?

gg <- ggplot(data, aes(x=runtime, y=rmse))
gg <- gg + geom_point(aes(fill=names,size=unmapped), shape=21, color="black", alpha=0.8,stroke=0.3)
gg <- gg + geom_text(aes(label=names),size=4,nudge_y=0.02,parse=T)
gg <- gg + ylim(0,0.65)
gg <- gg + theme_tufte()
gg <- gg + scale_x_log10(breaks=c(1,10,100), minor_breaks=c(seq(10,100,20), seq(100,1000,200)), labels=c("1","10","100"), limits = c(10,300))
gg <- gg + theme(plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
gg <- gg + theme(plot.background = element_blank(),
                 panel.border = element_blank(),
		 axis.text = element_text(size=14),
                 axis.title.x = element_text(size=16,margin=margin(20,0,0,0)),
                 axis.title.y = element_text(size=16,margin=margin(0,20,0,0)),
                 plot.title = element_text(size=20),
		 plot.subtitle = element_text(size=16,margin=margin(0,0,40,0)),
		 panel.grid.major = element_line(colour="gray", size=0.5),
		 panel.grid.minor = element_line(colour="gray", size=0.2),
		 legend.title=element_text(size=16),
		 legend.text=element_text(size=14),
                 )
gg <- gg + scale_fill_manual(values=c("FAME"="firebrick", "Bismark"="midnightblue", "BratNova"="tan4","BSMAP"="goldenrod1", "Segemehl"="darkgreen"), 
                       		name="Method",
                       		breaks=c("FAME", "Bismark", "BratNova", "BSMAP", "Segemehl"),
                       		labels=c("FAME", "Bismark", "BratNova", "BSMAP", "Segemehl")
)
gg <- gg + guides(fill=FALSE)
gg <- gg + scale_size_continuous(name="#Unmapped\nCpGs", breaks=c(15000,50000,100000,250000,500000), labels = comma)

gg <- gg + labs(x="Runtime in minutes (log scale)",
	        y="RMSE",
		title="Performance of methylation calling",
		subtitle="Synthetic data")

gg2 <- ggplot(data, aes(x=runtime, y=corr))
gg2 <- gg2 + geom_point(aes(fill=names,size=unmapped), shape=21, color="black", alpha=0.8,stroke=0.3)
gg2 <- gg2 + geom_text(aes(label=names),size=4,nudge_y=-0.03,parse=T)
gg2 <- gg2 + ylim(-0.4,1)
gg2 <- gg2 + theme_tufte()
gg2 <- gg2 + scale_x_log10(breaks=c(1,10,100), minor_breaks=c(seq(10,100,20), seq(100,1000,200)), labels=c("1","10","100"), limits = c(10,300))
gg2 <- gg2 + theme(plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
gg2 <- gg2 + theme(plot.background = element_blank(),
                 panel.border = element_blank(),
		 axis.text = element_text(size=14),
                 axis.title.x = element_text(size=16,margin=margin(20,0,0,0)),
                 axis.title.y = element_text(size=16,margin=margin(0,20,0,0)),
                 plot.title = element_text(size=20),
		 plot.subtitle = element_text(size=16,margin=margin(0,0,40,0)),
		 panel.grid.major = element_line(colour="gray", size=0.5),
		 panel.grid.minor = element_line(colour="gray", size=0.2),
		 legend.title=element_text(size=16),
		 legend.text=element_text(size=14),
                 )
gg2 <- gg2 + scale_fill_manual(values=c("FAME"="firebrick", "Bismark"="midnightblue", "BratNova"="tan4","BSMAP"="goldenrod1", "Segemehl"="darkgreen"), 
                       		name="Method",
                       		breaks=c("FAME", "Bismark", "BratNova", "BSMAP", "Segemehl"),
                       		labels=c("FAME", "Bismark", "BratNova", "BSMAP", "Segemehl")
)
gg2 <- gg2 + guides(fill=FALSE)
gg2 <- gg2 + scale_size_continuous(name="#Unmapped\nCpGs", breaks=c(15000,50000,100000,250000,500000), labels = comma)

gg2 <- gg2 + labs(x="Runtime in minutes (log scale)",
	        y="Spearman Correlation Coefficient",
		title="Performance of methylation calling",
		subtitle="Synthetic data")

pdf("summary_plot.pdf")

print(gg)
print(gg2)

dev.off()

