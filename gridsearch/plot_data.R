

k_par <- c(rep(c(25,28),times=c(9,9)), rep(c(25,28),times=c(3,3)))
#t_type <- rep(c("On","Off"),times=c(18,6))
#t_par <- c(rep(rep(c(500,2000,5000),times=c(3,3,3)),times=2),rep(-1,times=6))
t_par <- c(rep(rep(c("500","2000","5000"),times=c(3,3,3)),times=2),rep("Off",times=6))
m_par <- rep(c(512,1024,2048),8)

runtimes <- c(
# k=25
 # t=500
 2424,
 315,
 314,

 # t=2000
 2959,
 666,
 650,

 # t=5000
 3580,
 1284,
 1348,

# k=28
 # t=500
 285,
 315,
 315,

 # t=2000
 568,
 612,
 697,

 # t=5000
 1241,
 1419,
 1524,

# k=25
 # t=no
 61478,
 57354,
 62171,
# k=28
 # t=no
 61672,
 49563,
 44102

)

rmse <- c(
# k=25
 # t=500
 0.1264,
 0.1251,
 0.1236,

 # t=2000
 0.1148,
 0.1131,
 0.1119,

 # t=5000
 0.1101,
 0.1083,
 0.107,


# k=28
 # t=500
 0.1128,
 0.1104,
 0.1089,

 # t=2000
 0.1082,
 0.106,
 0.1046,

 # t=5000
 0.1068,
 0.1047,
 0.1033,

# k=25
 # t=no
 0.1054,
 0.1037,
 0.1026,

# k=28
 # t=no
 0.1055,
 0.1037,
 0.1027

)

#labels <- c(rep(c("m:512","m:1024","m:2048"),2),rep(NA,6),rep(c("m:512","m:1024","m:2048"),2),rep(NA,6))
labels <- c(rep(c("m:512","m:1024","m:2048")),rep(NA,9),rep(c("m:512","m:1024","m:2048")),rep(NA,9))


#data <- data.frame(k=k_par,t=t_par,Filter=t_type,m=m_par,labels=labels,runtimes=runtimes,RMSE=rmse)
data <- data.frame(k=k_par,t=t_par,m=m_par,labels=labels,runtimes=runtimes,RMSE=rmse)
data$m <- as.factor(data$m)

require(ggplot2)
require(ggthemes)
require(viridis)



pdf("grid_res.pdf")

gg <- ggplot(as.data.frame(data),aes(x=t,y=runtimes,group=m))
gg <- gg + geom_line(size=0.3,aes(color=m)) + theme_bw()
#gg <- gg + geom_point(size=1,aes(color=RMSE,fill=RMSE)) + facet_wrap(~ k + Filter,labeller = label_both)
#gg <- gg + geom_point(size=0.3) + facet_wrap(~ k + Filter,labeller = label_both)
gg <- gg + geom_point(size=0.3,aes(fill=m)) + facet_wrap(~ k,labeller = label_both)
#gg <- gg + scale_fill_viridis(option="B")
#gg <- gg + scale_color_viridis(option="B")
gg <- gg + theme(strip.text = element_text(size=7,lineheight=5.0),
strip.background = element_rect(fill="white", colour="black",
size=0.8))
gg <- gg + scale_x_discrete(limits=c("500","2000","5000","Off"))
gg <- gg + theme(plot.margin = unit(c(2,2,2,2),"cm"))
gg <- gg + theme(panel.border = element_blank(),
                 axis.title.x = element_text(margin=margin(20,0,0,0)),
                 axis.title.y = element_text(margin=margin(0,20,0,0))
                )
gg <- gg + labs(x="Filter threshold", y="Runtime (in seconds)")
#gg <- gg + labs(title="Gridsearch results on chromosome 22 performance", subtitle="10 Million reads drawn from chromosome 22\nUp to two errors, bimodal methylation rates")
#gg <- gg + geom_text(label=labels,na.rm=T,nudge_x="500",nudge_y=2000,size=2)

gg

gg2 <- ggplot(as.data.frame(data),aes(x=t,y=RMSE,fill=m))
gg2 <- gg2 + geom_bar(stat='identity',position='dodge',width=0.5) + theme_bw()
#gg2 <- gg2 + coord_cartesian(ylim=c(0.19,0.205))
gg2 <- gg2 + theme(plot.margin = unit(c(2,2,2,2),"cm"))
gg2 <- gg2 + facet_wrap(~ k,labeller = label_both)
gg2 <- gg2 + theme(strip.text = element_text(size=7,lineheight=5.0),
strip.background = element_rect(fill="white", colour="black",
size=0.8))
gg2 <- gg2 + theme(panel.border = element_blank(),
		 axis.title.x = element_text(margin=margin(20,0,0,0)),
                 axis.title.y = element_text(margin=margin(0,20,0,0))
                )
gg2 <- gg2 + scale_x_discrete(limits=c("500","2000","5000","Off"))
gg2 <- gg2 + labs(x="Filter threshold")
gg2

dev.off()
