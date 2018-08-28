require(ggplot2)
require(reshape2)
require(ggthemes)
require(energy)
require(dplyr)
# Arguments:
# 1)  file with real values of forward strand CpGs (how they were generated)
# 2)  file with real values of reverse strand CpGs (how they were generated)
# 3)  file with estimated values
# 4)  output file name
# 5)  plot title

#args <- commandArgs(trailingOnly = TRUE)
args <- c("/MMCI/MS/MethylationAlgos/work/BisulfitePipeline/Metal/Synth/paired_hg19_CHR22_poisson0_5_cpginfo_fwd.tsv",
"/MMCI/MS/MethylationAlgos/work/BisulfitePipeline/Metal/Synth/paired_hg19_CHR22_poisson0_5_cpginfo_rev.tsv",
"/MMCI/MS/EpiregDeep2/nobackup/Jonas_Directory/paired_hg19_CHR22_poisson0_5_p1_bismark_bt2_pe.CpG_report.txt",
"bismark_chr22_pe_poisson_performance_union_5r.pdf",
"25 Million reads on hg19 Chromosome 22")

covered_cpgs <- do.call(cbind,sapply(c("bratnova","bismark","fame","bsmap"), function(x) {return(as.vector(read.table(paste("/MMCI/MS/MethylationAlgos/work/BisulfitePipeline/Metal/results/",x,"_mapped_regions.tsv",sep=""), head=F)))}))
covered_cpgs <- covered_cpgs[,1] | covered_cpgs[,2] | covered_cpgs[,3] | covered_cpgs[,4]


data_real_fwd <- read.table(args[1],header=T,sep="\t",colClasses=c("character","integer","integer","integer","double"))[,-5]
data_real_rev <- read.table(args[2],header=T,sep="\t",colClasses=c("character","integer","integer","integer","double"))[,-5]
data_real_fwd <- cbind(data_real_fwd, rep('+',nrow(data_real_fwd)) )
data_real_rev <- cbind(data_real_rev, rep('-',nrow(data_real_rev)) )
colnames(data_real_fwd) <- c("Chromosome","Position","Unmethylated","Methylated","Strand")
colnames(data_real_rev) <- c("Chromosome","Position","Unmethylated","Methylated","Strand")

data_real_fwd$Position <- data_real_fwd$Position + 1
data_real_rev$Position <- data_real_rev$Position + 2

data_real <- rbind(data_real_fwd, data_real_rev)
data_real <- data_real[grep("chr22",data_real$Chromosome),]


rmse <- function(x1, x2)
{
	return(sqrt(mean((x1-x2)^2,na.rm=T)))
}

colnames(data_real) <- c("Chromosome","Position","Unmethylated","Methylated","Strand")
data_real <- transform(data_real, RealRate = Methylated / (Unmethylated + Methylated))
#data_real <- data_real[!is.na(data_real$RealRate) & !is.nan(data_real$RealRate),]
#data_real$Chromosome <- rep('chr22',nrow(data_real))


data_bis <- read.table(args[3],sep="\t",stringsAsFactors=F,colClasses=c("character","integer","character","integer","integer","character","character"))[,c(-6,-7)]
colnames(data_bis) <- c("Chrom","Pos","Strand","Meth","Unmeth")
data_bis <- transform(data_bis, PredRate = Meth / (Unmeth + Meth))


data_org <- left_join(data_real, data_bis, by=c("Chromosome" = "Chrom","Position" = "Pos","Strand"="Strand"))

data_org <- data_org[!is.na(data_org$RealRate) & !is.nan(data_org$RealRate),]

#print(paste("#Unmapped CpGs: ", sum((is.na(data_org$PredRate) | is.nan(data_org$PredRate)) ), sep=""))
#write(!is.na(data_org$PredRate) & !is.nan(data_org$PredRate) & ((data_org$Meth + data_org$Unmeth) > 5),"bismark_mapped_regions.tsv",ncolumns=1)


data <- data_org[!is.na(data_org$PredRate) & !is.nan(data_org$PredRate),]

# replace NAs with 0.5
data_org <- data_org[covered_cpgs,]
data_org[is.na(data_org$PredRate) | is.nan(data_org$PredRate),]$PredRate <- 0.5


plot_all <- function()
{


cor_pe_my <- round(cor(data$RealRate, data$PredRate, method="pearson",use='complete.obs'), digits=2)
cor_sp_my <- round(cor(data$RealRate, data$PredRate, method="spearman",use='complete.obs'), digits=2)
rmse_my <- round(rmse(data$RealRate, data$PredRate), digits=2)

rmse_meth <- round(rmse(data$Meth, data$Methylated), digits=3)
rmse_unmeth <- round(rmse(data$Unmeth, data$Unmethylated), digits=3)
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
gg <- gg + labs(title="Aligned methylated reads (Bismark) vs\n actual number of methylated reads",subtitle=args[5])

gg <- gg +  annotate("text", color="black", size=3, x = 10, y = max(data$Meth)-1/5*max(data$Meth), label = paste("RMSE:", rmse_meth, sep=" "))

print(gg)

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
gg <- gg + labs(title="Aligned unmethylated reads (Bismark) vs\n actual number of unmethylated reads",subtitle=args[5])
gg <- gg +  annotate("text", color="black", size=3, x = 10, y = max(data$Unmeth)-1/5*max(data$Unmeth), label = paste("RMSE:", rmse_unmeth, sep=" "))
print(gg)


gg <- ggplot(data,aes(x=PredRate, y=RealRate))
gg <- gg + geom_hex(bins = 100) + scale_fill_gradientn("", colours = c("gray88","lightgoldenrod1","orange2","firebrick","black"), values=c(0,0.01,0.1,0.5,1), limits=c(0,7000))
gg <- gg + theme_tufte()
gg <- gg + theme(plot.margin = unit(c(2,2,2,2),"cm"))
gg <- gg + theme(plot.background = element_blank(),
                              panel.border = element_blank(),
                              axis.title.x = element_text(margin=margin(20,0,0,0)),
                              axis.title.y = element_text(margin=margin(0,20,0,0)),
                              plot.title = element_text(margin=margin(0,0,40,0))
                              )
gg <- gg + labs(x="Predicted methylation rate", y="Real methylation rate")

gg <- gg +  annotate("text", color="darkred", size=3, x = 0.75, y = 0.2, label = paste("RMSE:", rmse_my, sep=" "))

gg <- gg +  annotate("text", color="darkred", size=3, x = 0.75, y = 0.15, label = paste("Spearman:", cor_sp_my, sep=" "))

gg <- gg +  annotate("text", color="darkred", size=3, x = 0.75, y = 0.1, label = paste("Pearson:", cor_pe_my, sep=" "))

gg <- gg + annotate("rect", xmin = 0.6, xmax = 0.9, ymin = 0.08, ymax = 0.22, alpha = .1)
gg <- gg + labs(title="Methylation rates of Bismark",subtitle=args[5])

print(gg)
}


pdf(args[4])

plot_all()
data <- data_org
plot_all()

dev.off()
