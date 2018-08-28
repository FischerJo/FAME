require(ggplot2)
require(reshape2)
require(ggthemes)
require(energy)

# Arguments:
# 1)  file with real values (how they were generated)
# 2)  file with estimated values
# 3)  output file name
# 4)  plot title

args <- commandArgs(trailingOnly = TRUE)

data_real <- read.table(args[1],header=F,sep="\t",nrows=29401360,colClasses=c("integer","integer","integer","integer","double"))[,-5]
colnames(data_real) <- c("Chromosome","Position","Unmethylated","Methylated")
data_real <- transform(data_real, RealRate = Methylated / (Unmethylated + Methylated))
data_real$RealRate <- replace(data_real$RealRate, is.na(data_real$RealRate), mean(data_real$RealRate, na.rm=T))

data_real$Position <- data_real$Position + 1
#
data_bss <- read.table(args[2],nrows=29401360,colClasses=c("character","character","integer","character","character","character","integer","integer"))[,c(-2,-4,-5,-6)]
colnames(data_bss) <- c("Chrom","Position","Meth","Unmeth")
data_bss$Unmeth <- data_bss$Unmeth - data_bss$Meth


data <- merge(data_real, data_bss, by=c("Position"), all.x=TRUE)

data$Meth[is.na(data$Meth)] <- 0
data$Unmeth[is.na(data$Unmeth)] <- 0
data <- transform(data, PredRate = Meth / (Unmeth + Meth))
mean_meth <- mean(data$PredRate, na.rm=T)

data$PredRate <- replace(data$PredRate, is.na(data$PredRate), mean_meth)



rmse_bss <- round(sqrt(sum((data$Methylated - data$Meth)^2)/nrow(data)),digits=4)
pdf(args[3])
# plot the correlation between aligned reads
gg <- ggplot(data,aes(x=Methylated, y=Meth))
gg <- gg + geom_hex(bins = 100) + scale_fill_gradientn("", colours = c("gray88","lightgoldenrod1","orange2","firebrick","black"), values=c(0,0.01,0.1,0.5,1))
gg <- gg + theme_tufte()
gg <- gg + theme(plot.margin = unit(c(2,2,2,2),"cm"))
gg <- gg + theme(plot.background = element_blank(),
                              panel.border = element_blank(),
                              axis.title.x = element_text(margin=margin(20,0,0,0)),
                              axis.title.y = element_text(margin=margin(0,20,0,0)),
                              plot.title = element_text(margin=margin(0,0,40,0))
                              )
gg <- gg + labs(x="Number of generated methylated reads for this Position", y="Aligned methylated reads")
gg <- gg + labs(title="Aligned methylated reads (BSseeker2) vs\n actual number of methylated reads",subtitle=args[4])

gg <- gg +  annotate("text", color="black", size=3, x = 10, y = max(data_bss$Meth)-3, label = paste("RMSE:", rmse_bss, sep=" "))
gg
rmse_bss <- round(sqrt(sum((data$Unmethylated - data$Unmeth)^2)/nrow(data)),digits=4)
gg <- ggplot(data,aes(x=Unmethylated, y=Unmeth))
gg <- gg + geom_hex(bins = 100) + scale_fill_gradientn("", colours = c("gray88","lightgoldenrod1","orange2","firebrick","black"), values=c(0,0.01,0.1,0.5,1))
gg <- gg + theme_tufte()
gg <- gg + theme(plot.margin = unit(c(2,2,2,2),"cm"))
gg <- gg + theme(plot.background = element_blank(),
                              panel.border = element_blank(),
                              axis.title.x = element_text(margin=margin(20,0,0,0)),
                              axis.title.y = element_text(margin=margin(0,20,0,0)),
                              plot.title = element_text(margin=margin(0,0,40,0))
                              )
gg <- gg + labs(x="Number of generated unmethylated reads for this Position", y="Aligned unmethylated reads")
gg <- gg + labs(title="Aligned unmethylated reads (BSseeker2) vs\n actual number of unmethylated reads",subtitle=args[4])

gg <- gg +  annotate("text", color="black", size=3, x = 10, y = max(data_bss$Unmeth)-10, label = paste("RMSE:", rmse_bss, sep=" "))
gg


cor_pe_bss <- round(cor(data$RealRate, data$PredRate, method="pearson",use='complete.obs'), digits=2)
cor_sp_bss <- round(cor(data$RealRate, data$PredRate, method="spearman",use='complete.obs'), digits=2)
rmse_bss <- round(sqrt(sum((data$RealRate - data$PredRate)^2)/nrow(data)),digits=4)


# plot bssmark data
gg <- ggplot(data,aes(x=PredRate, y=RealRate))
#    geom_point(alpha=0.2,size=0.3)
gg <- gg + geom_hex(bins = 100) + scale_fill_gradientn("", colours = c("gray88","lightgoldenrod1","orange2","firebrick","black"), values=c(0,0.01,0.1,0.5,1))
gg <- gg + theme_tufte()
gg <- gg + theme(plot.margin = unit(c(2,2,2,2),"cm"))
gg <- gg + theme(plot.background = element_blank(),
                              panel.border = element_blank(),
                              axis.title.x = element_text(margin=margin(20,0,0,0)),
                              axis.title.y = element_text(margin=margin(0,20,0,0)),
                              plot.title = element_text(margin=margin(0,0,40,0))
                              )
gg <- gg + labs(x="Predicted methylation rate", y="Real methylation rate")

gg <- gg +  annotate("text", color="darkred", size=3, x = 0.75, y = 0.15, label = paste("Spearman:", cor_sp_bss, sep=" "))

gg <- gg +  annotate("text", color="darkred", size=3, x = 0.75, y = 0.1, label = paste("Pearson:", cor_pe_bss, sep=" "))
gg <- gg +  annotate("text", color="darkred", size=3, x = 0.75, y = 0.2, label = paste("RMSE:", rmse_bss, sep=" "))
#gg <- gg +  annotate("text", color="darkred", size=3, x = 0.75, y = 0.9, label = paste("Proportion of NAs:",na_perc, sep=" "))

gg <- gg + annotate("rect", xmin = 0.6, xmax = 0.9, ymin = 0.08, ymax = 0.22, alpha = .1)
gg <- gg + labs(title="Methylation rates of BSseeker2",subtitle=args[4])

gg

dev.off()
