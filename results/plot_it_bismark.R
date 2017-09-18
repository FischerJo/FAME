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


#
data_bis <- read.table(args[2],nrows=29401360,colClasses=c("character","integer","character","integer","integer","character","character"))[,c(-3,-6,-7)]
colnames(data_bis) <- c("Chrom","Pos","Meth","Unmeth")
data_bis <- transform(data_bis, PredRate = Meth / (Unmeth + Meth))

mean_meth <- mean(data_bis$PredRate, na.rm=T)

print(sum(data_bis$Unmeth + data_bis$Meth))
rmse_bis <- round(sqrt(sum((data_real$Methylated - data_bis$Meth)^2)/nrow(data_real)),digits=4)
pdf(args[3])
# plot the correlation between aligned reads
gg <- ggplot(cbind(data_real,data_bis),aes(x=Methylated, y=Meth))
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
gg <- gg + labs(title="Aligned methylated reads (Bismark) vs\n actual number of methylated reads",subtitle=args[4])

gg <- gg +  annotate("text", color="black", size=3, x = 10, y = max(data_bis$Meth)-3, label = paste("RMSE:", rmse_bis, sep=" "))
gg
rmse_bis <- round(sqrt(sum((data_real$Unmethylated - data_bis$Unmeth)^2)/nrow(data_real)),digits=4)
gg <- ggplot(cbind(data_real,data_bis),aes(x=Unmethylated, y=Unmeth))
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
gg <- gg + labs(title="Aligned unmethylated reads (Bismark) vs\n actual number of unmethylated reads",subtitle=args[4])

gg <- gg +  annotate("text", color="black", size=3, x = 10, y = max(data_bis$Unmeth)-10, label = paste("RMSE:", rmse_bis, sep=" "))
gg

na_perc <- round(sum(is.na(data_bis$PredRate))/nrow(data_bis),digits=2)

data_bis$PredRate <- replace(data_bis$PredRate, is.na(data_bis$PredRate), mean_meth)


cor_pe_bis <- round(cor(data_real[,5], data_bis[,5], method="pearson",use='complete.obs'), digits=2)
cor_sp_bis <- round(cor(data_real[,5], data_bis[,5], method="spearman",use='complete.obs'), digits=2)
rmse_bis <- round(sqrt(sum((data_real$RealRate - data_bis$PredRate)^2)/nrow(data_real)),digits=4)

#sample_ids <- sample(nrow(data_real),20000)
#round(dcor(data_real[sample_ids,5], data_bis[sample_ids,5]),digits=2)

# plot bismark data
gg <- ggplot(cbind(data_real,data_bis),aes(x=PredRate, y=RealRate))
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

gg <- gg +  annotate("text", color="darkred", size=3, x = 0.75, y = 0.15, label = paste("Spearmen:", cor_sp_bis, sep=" "))

gg <- gg +  annotate("text", color="darkred", size=3, x = 0.75, y = 0.1, label = paste("Pearson:", cor_pe_bis, sep=" "))
gg <- gg +  annotate("text", color="darkred", size=3, x = 0.75, y = 0.2, label = paste("RMSE:", rmse_bis, sep=" "))
#gg <- gg +  annotate("text", color="darkred", size=3, x = 0.75, y = 0.9, label = paste("Proportion of NAs:",na_perc, sep=" "))

gg <- gg + annotate("rect", xmin = 0.6, xmax = 0.9, ymin = 0.08, ymax = 0.22, alpha = .1)
gg <- gg + labs(title="Methylation rates of Bismark",subtitle=args[4])

gg

dev.off()
