

k_par <- c(rep(c(23,25,28),times=c(20,20,20)))
t_par <- rep(rep(c("500","1500","3500", "5000", "10000"),times=c(4,4,4,4,4)),times=3)
m_par <- rep(c(1024,2048,4096,8192),15)

runtimes <- c(
# k=23
 # t=500
 285,
 264,
 252,
 245,

 # t=1500
 455,
 432,
 402,
 374,

 # t=3500
 783,
 707,
 619,
 577,

 # t=5000
 868,
 833,
 855,
 684,

 # t=10000
 1163,
 1340,
 1088,
 1298,

# k=25
 # t=500
 232,
 175,
 206,
 261,

 # t=1500
 271,
 292,
 282,
 306,

 # t=3500
 551,
 535,
 525,
 523,

 # t=5000
 663,
 629,
 662,
 643,


 # t=10000
 1039,
 1068,
 1015,
 1121,


# k=28
 # t=500
 180,
 189,
 183,
 233,

 # t=1500
 291,
 274,
 365,
 414,

 # t=3500
 523,
 568,
 618,
 897,

 # t=5000
 600,
 685,
 848,
 1035,

 # t=10000
 993,
 1003,
 1148,
 1823


# k=25
 # t=no
 #61478,
 #57354,
 #62171,
# k=28
 # t=no
 #61672,
 #49563,
 #44102

)

rmse <- c(
# k=23
 # t=500
 0.1365,
 0.1339,
 0.1346,
 0.135,

 # t=1500
 0.1246,
 0.1229,
 0.1229,
 0.1217,

 # t=3500
 0.1083,
 0.1154,
 0.1153,
 0.1136,

 # t=5000
 0.1045,
 0.1125,
 0.1125,
 0.111,

 # t=10000
 0.1002,
 0.1087,
 0.1085,
 0.1073,

# k=25
 # t=500
 0.1229,
 0.1201,
 0.1207,
 0.1213,

 # t=1500
 0.1132,
 0.1115,
 0.1116,
 0.1101,

 # t=3500
 0.0995,
 0.1073,
 0.1072,
 0.1053,

 # t=5000
 0.0973,
 0.1056,
 0.1057,
 0.1043,

 # t=10000
 0.096,
 0.1038,
 0.1037,
 0.1028,


# k=28
 # t=500
 0.108,
 0.1045,
 0.1054,
 0.106,

 # t=1500
 0.1044,
 0.1027,
 0.1029,
 0.101,

 # t=3500
 0.0942,
 0.1021,
 0.1021,
 0.0996,

 # t=5000
 0.0941,
 0.1018,
 0.1017,
 0.1005,

 # t=10000
 0.0954,
 0.1017,
 0.1016,
 0.1007


# k=25
 # t=no
 #0.1054,
 #0.1037,
 #0.1026,

# k=28
 # t=no
# 0.1055,
 #0.1037,
 #0.1027

)


labels <- rep(c(c("m:1024","m:2048","m:4096","m:8192"),rep(NA,4*4)),3)

# convert to minutes
runtimes <- runtimes / 60

data <- data.frame(k=k_par,t=t_par,m=m_par,labels=labels,runtimes=runtimes,RMSE=rmse)
data$m <- as.factor(data$m)

require(ggplot2)
require(ggthemes)
require(viridis)

breaks <- 10^(-10:10)
minor_breaks <- rep(c(1:9), 21)*(10^rep(-10:10, each=9))


pdf("grid_res.pdf")

gg <- ggplot(as.data.frame(data),aes(x=t,y=runtimes,group=m))
gg <- gg + geom_line(size=0.3,aes(color=m)) + theme_bw()
gg <- gg + geom_point(size=0.3,aes(fill=m)) + facet_wrap(~ k,labeller = label_both)
gg <- gg + theme(strip.text = element_text(size=7,lineheight=5.0),
strip.background = element_rect(fill="white", colour="black",
size=0.8))
gg <- gg + scale_x_discrete(limits=c("500","1500","3500","5000","10000"))
gg <- gg + scale_y_log10(breaks = breaks,
   labels = sapply(log10(breaks), scales::math_format(10^.x)), minor_breaks = minor_breaks)
gg <- gg + theme(plot.margin = unit(c(2,2,2,2),"cm"))
gg <- gg + theme(panel.border = element_blank(),
                 axis.title.x = element_text(margin=margin(20,0,0,0)),
                 axis.title.y = element_text(margin=margin(0,20,0,0)),
                 axis.text.x = element_text(angle = 45, hjust = 1)
                )
gg <- gg + labs(x="Filter threshold", y="Runtime (logscale minutes)")

gg

gg2 <- ggplot(as.data.frame(data),aes(x=t,y=RMSE,fill=m))
gg2 <- gg2 + geom_bar(stat='identity',position='dodge',width=0.5) + theme_bw()
gg2 <- gg2 + theme(plot.margin = unit(c(2,2,2,2),"cm"))
gg2 <- gg2 + facet_wrap(~ k,labeller = label_both)
gg2 <- gg2 + theme(strip.text = element_text(size=7,lineheight=5.0),
strip.background = element_rect(fill="white", colour="black",
size=0.8))
gg2 <- gg2 + theme(panel.border = element_blank(),
		 axis.title.x = element_text(margin=margin(20,0,0,0)),
                 axis.title.y = element_text(margin=margin(0,20,0,0)),
                 axis.text.x = element_text(angle = 45, hjust = 1)
                )
gg2 <- gg2 + scale_x_discrete(limits=c("500","1500","3500","5000","10000"))
gg2 <- gg2 + labs(x="Filter threshold")
gg2

dev.off()
