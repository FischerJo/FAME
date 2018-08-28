require(ggplot2)
require(ggthemes)
require(reshape2)
library(plyr)

args <- commandArgs(trailingOnly = TRUE)

# Arguments:
# 1) Metal file with predicted values
# 2) Bismark file -"-
# 3) BSMAP file -"-
# 4) BSSeeker file -"-

data_my <- read.table(args[1],sep="\t",nrows=29401360,colClasses=c("integer","integer","integer","integer","integer","integer"))[,c(-5,-6)]
colnames(data_my) <- c("Chrom","Position","Methy","Unmethy")
data_my <- transform(data_my, PredRate = Methy / (Unmethy + Methy))
data_my$PredRate <- replace(data_my$PredRate, is.na(data_my$PredRate), -0.01)
data_my$Position <- data_my$Position + 1

data_bis <- read.table(args[2],nrows=29401360,colClasses=c("character","integer","character","integer","integer","character","character"))[,c(-3,-6,-7)]
colnames(data_bis) <- c("Chrom","Pos","Meth","Unmeth")
data_bis <- transform(data_bis, PredRate = Meth / (Unmeth + Meth))
data_bis$PredRate <- replace(data_bis$PredRate, is.na(data_bis$PredRate), -0.01)

data <- data.frame(Bismark=data_bis$PredRate, Metal=data_my$PredRate)
data <- data[(data$Bismark != -0.01 | data$Metal != -0.01),]

gg1 <- ggplot(data, aes(Bismark,Metal))
gg1 <- gg1 + theme_tufte()
gg1 <- gg1 + geom_hex(bins = 100)
gg1 <- gg1 + scale_fill_gradientn("", colours = c("gray88","lightgoldenrod1","orange2","firebrick","black"), values=c(0,0.1,0.2,0.6,1))
gg1 <- gg1 + theme(plot.margin = unit(c(2,2,2,2),"cm"))
gg1 <- gg1 + theme(plot.background = element_blank(),
                              panel.border = element_blank(),
                              axis.title.x = element_text(margin=margin(20,0,0,0)),
                              axis.title.y = element_text(margin=margin(0,20,0,0)),
                              plot.title = element_text(margin=margin(0,0,40,0))
                              )
gg1 <- gg1 + labs(x="Predicted methylation rate Bismark", y="Predicted methylation rate Metal")



data_bsm <- read.table(args[3])[,c(2,5,7,8,9,10)]
data_bsm[,3] <- data_bsm[,3] + data_bsm[,5]
data_bsm[,4] <- data_bsm[,4] + data_bsm[,6]
data_bsm <- data_bsm[,c(-5,-6)]
colnames(data_bsm) <- c("Position","PredRate2","Meth","Unmeth")
data_bsm$Unmeth <- data_bsm$Unmeth - data_bsm$Meth

data2 <- merge(data_my, data_bsm, by=c("Position"), all.x=TRUE)
data2$Meth[is.na(data2$Meth)] <- 0
data2$Unmeth[is.na(data2$Unmeth)] <- 0

data2$PredRate2 <- data2$Meth / (data2$Unmeth + data2$Meth)
data2$PredRate2 <- replace(data2$PredRate2, is.na(data2$PredRate2), -0.01)

data2 <- data.frame(BSMAP=data2$PredRate2, Metal=data2$PredRate)
data2 <- data2[(data2$BSMAP != -0.01 | data2$Metal != -0.01),]


gg2 <- ggplot(data2, aes(BSMAP,Metal))
gg2 <- gg2 + theme_tufte()
gg2 <- gg2 + geom_hex(bins = 100)
gg2 <- gg2 + scale_fill_gradientn("", colours = c("gray88","lightgoldenrod1","orange2","firebrick","black"), values=c(0,0.1,0.2,0.6,1))
gg2 <- gg2 + theme(plot.margin = unit(c(2,2,2,2),"cm"))
gg2 <- gg2 + theme(plot.background = element_blank(),
                              panel.border = element_blank(),
                              axis.title.x = element_text(margin=margin(20,0,0,0)),
                              axis.title.y = element_text(margin=margin(0,20,0,0)),
                              plot.title = element_text(margin=margin(0,0,40,0))
                              )
gg2 <- gg2 + labs(x="Predicted methylation rate BSMAP", y="Predicted methylation rate Metal")

data_bss <- read.table(args[4],nrows=29401360,colClasses=c("character","character","integer","character","character","character","integer","integer"))[,c(-2,-4,-5,-6)]
colnames(data_bss) <- c("Chrom","Position","Meth","Unmeth")
data_bss$Unmeth <- data_bss$Unmeth - data_bss$Meth
data3 <- merge(data_my, data_bss, by=c("Position"), all.x=TRUE)

data3$Meth[is.na(data3$Meth)] <- 0
data3$Unmeth[is.na(data3$Unmeth)] <- 0
data3 <- transform(data3, PredRate3 = Meth / (Unmeth + Meth))

data3$PredRate3 <- replace(data3$PredRate3, is.na(data3$PredRate3), -0.01)
data3 <- data.frame(BSSeeker2=data3$PredRate3, Metal=data3$PredRate)
data3 <- data3[(data3$BSSeeker2 != -0.01 | data3$Metal != -0.01),]

gg3 <- ggplot(data3, aes(BSSeeker2,Metal))
gg3 <- gg3 + theme_tufte()
gg3 <- gg3 + geom_hex(bins = 100)
gg3 <- gg3 + scale_fill_gradientn("", colours = c("gray88","lightgoldenrod1","orange2","firebrick","black"), values=c(0,0.1,0.2,0.6,1))
gg3 <- gg3 + theme(plot.margin = unit(c(2,2,2,2),"cm"))
gg3 <- gg3 + theme(plot.background = element_blank(),
                              panel.border = element_blank(),
                              axis.title.x = element_text(margin=margin(20,0,0,0)),
                              axis.title.y = element_text(margin=margin(0,20,0,0)),
                              plot.title = element_text(margin=margin(0,0,40,0))
                              )
gg3 <- gg3 + labs(x="Predicted methylation rate BSSeeker2", y="Predicted methylation rate Metal")



pdf("ratio_corrs.pdf")

gg1
gg2
gg3

dev.off()
