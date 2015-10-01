library(ggplot2)

h3k4 = read.table("case_control_LGD.txt")
str(h3k4)
head(h3k4)
colnames(h3k4) = c("gene","EncodeNormalPercent","EncodeNormalLength","EncodeCancerPercent","EncodeCancerLength","RoadmapNormalPercent","RoadmapNormalLength","case")

sapply(h3k4[,2:7],max) #check for potential outlier

ggplot(h3k4,aes(x=case,y=RoadmapNormalLength,fill=case)) + geom_boxplot() + geom_violin(fill=NA,col='blue') + guides(fill=FALSE) +xlab('')+ylab("H3K4me3 peak length (bp)")+theme_grey(base_size = 20) + scale_y_log10() 
wilcox.test(h3k4[,"EncodeNormalLength"]~h3k4[,"case"],alternative="greater")
#W = 20724, p-value = 0.9847
