require(VennDiagram)

# Arguments:
# 1)  file with real values (how they were generated)
# 2)  file with estimated values Metal
# 3)  file with estimated values Bismark
# 3)  output file name


args <- commandArgs(trailingOnly = TRUE)

data_real <- read.table(args[1],header=F,sep="\t",nrows=29401360,colClasses=c("integer","integer","integer","integer","double"))[,-5]
colnames(data_real) <- c("Chromosome","Position","Unmethylated","Methylated")
data_real <- transform(data_real, RealRate = Methylated / (Unmethylated + Methylated))
data_real$RealRate <- replace(data_real$RealRate, is.na(data_real$RealRate), mean(data_real$RealRate, na.rm=T))



data_my <- read.table(args[2],sep="\t",nrows=29401360,colClasses=c("integer","integer","integer","integer","integer","integer"))[,c(-5,-6)]
colnames(data_my) <- c("Chrom","Pos","Meth","Unmeth")
data_my <- transform(data_my, PredRate = Meth / (Unmeth + Meth))

my_notNA <- which(!is.na(data_my$PredRate))
#
data_bis <- read.table(args[3],nrows=29401360,colClasses=c("character","integer","character","integer","integer","character","character"))[,c(-3,-6,-7)]
colnames(data_bis) <- c("Chrom","Pos","Meth","Unmeth")
data_bis <- transform(data_bis, PredRate = Meth / (Unmeth + Meth))
bis_notNA <- which(!is.na(data_bis$PredRate))


pdf(args[4])


draw.triple.venn(length(my_notNA), length(bis_notNA), nrow(data_real), length(intersect(my_notNA, bis_notNA)), length(bis_notNA), length(my_notNA), length(intersect(my_notNA, bis_notNA)), euler.d=T, scaled=T,fill=c("red","green","black"),alpha=c(0.2,0.2,0.2), lwd=c(0.5,0.5,0.5),category=c("Metal","Bismark","All"), overrideTriple=1)




dev.off()
