
require(ggplot2)
require(reshape2)
require(ggthemes)

randData <- read.table("hashedRand_1bill.tsv", sep='\t', col.names=c('Bin', 'Random'))
distData <- read.table("hashedDist_1bill.tsv", sep='\t', col.names=c('Bin', 'Distinct'))



wholeData <- as.data.frame(cbind(randData, distData[,2]))
colnames(wholeData) <- c('Bin', 'Random', 'Distinct')



meltedData <- melt(wholeData, id.vars='Bin')


# generate counts
maxi <- max(wholeData[,-1])
randStats <- cbind(1:maxi, tabulate(randData[,2], maxi))
distStats <- cbind(1:maxi, tabulate(distData[,2], maxi))

colnames(randStats) <- c('collisions', 'Random')
colnames(distStats) <- c('collisions', 'Distinct')

wholeStats <- as.data.frame(cbind(randStats, distStats[,2]))
colnames(wholeStats) <- c('collisions', 'Random', 'Distinct')
meltedStats <- melt(wholeStats, id.vars='collisions')
colnames(meltedStats) <- c('collisions', 'variable', 'value')

pdf('ntHash_stats_1Bill.pdf')

gg1 <- ggplot(randData, aes(x=Bin, y=Random))
gg1 <- gg1 + theme_tufte() + geom_hex(binwidth=c(10000000,1.5))
# gg1 <- gg1 + theme(axis.ticks=element_blank())
gg1 <- gg1 + scale_y_continuous("# hashed values", limits=c(0,100))
gg1 <- gg1 + scale_x_continuous("Bin", breaks=c(0, 500000000, 1000000000), labels=c('0', '0.5bill', '1bill.'))
gg1 <- gg1 + ggtitle("3^20 random 20-mers hashed into 1 Billion bins")
gg1 <- gg1 + theme(plot.margin = unit(c(0.5,0,1,1),"cm"))
gg1

gg1_2 <- ggplot(distData, aes(x=Bin, y=Distinct))
gg1_2 <- gg1_2 + theme_tufte() + geom_hex(binwidth=c(10000000,1.3))
# gg1_2 <- gg1_2 + theme(axis.ticks=element_blank())
gg1_2 <- gg1_2 + scale_y_continuous("# hashed values", limits=c(0,100))
gg1_2 <- gg1_2 + scale_x_continuous("Bin", breaks=c(0, 500000000, 1000000000), labels=c('0', '0.5bill', '1bill'))
gg1_2 <- gg1_2 + ggtitle("All distinct 20-mers hashed into 1 Billion bins")
gg1_2 <- gg1_2 + theme(plot.margin = unit(c(0.5,0,1,1),"cm"))
gg1_2


gg2 <- ggplot(meltedStats, aes(x=collisions, y=value, color=variable, fill= variable, group=variable))
gg2 <- gg2 + theme_tufte() + geom_point(size=1)
gg2 <- gg2 + scale_y_log10("# of encounters of so many collisions (logscale)", breaks=c(10,1000,100000,10000000),labels=c('10','1thsd','100thsd','10mill'))
gg2 <- gg2 + scale_fill_manual(values=c("darkred", "goldenrod1")) 
gg2 <- gg2 + scale_color_manual(values=c("darkred", "goldenrod1")) 
gg2 <- gg2 + ggtitle("Hashing statistic about collisions")
gg2 <- gg2 + theme(plot.margin = unit(c(0.5,0.5,1,1),"cm"))
gg2

dev.off()
