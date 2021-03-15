library(ggplot2)
library(vroom)
library(plyr)
library(dplyr)
library(tidyverse)
library(scales)


options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
name <- args[2]

#setwd('~/Work/TempAnalysis/denbicloud/RIssMed')
#file <- 'InfoSubset_20_7.tsv.gz'
#name <- 'GC20_Cons7'

data <- vroom(file, delim="\t", col_names=c("Delta_acc", "Distance", "Acc_raw", "Zscore","Type"), col_types=c("Delta_acc"="n","Distance"="i","Acc_raw"="n","Zscore"="n","Type"="c"), num_threads=6)

data <- data %>% mutate(bin = cut_width(Distance, width = 10, center = 0, closed="left")) %>% group_by(bin)

#data$Distbin <- findInterval(data$Distance, c(seq(min(data$Distance),max(data$Distance))), all.inside=T)
#p <- ggplot(data, aes(x = factor(round_any(Distbin,10)), y = Delta_acc)) + geom_violin() + stat_summary(fun=mean, colour="black", geom="point", shape=18, size=3, show.legend = FALSE) + scale_fill_manual(values = c("firebrick"), alpha(.7)) + guides(fill=FALSE)

p <- ggplot(data, aes(x=Distance, y=Delta_acc)) + geom_jitter(aes(color='blue'),alpha=0.2) + geom_violin(fill="bisque",color="black",alpha=0.3) + stat_summary(fun=mean, colour="black", geom="point", shape=18, size=3, show.legend = FALSE)# + scale_fill_manual(values = c("firebrick"), alpha(.4)) + guides(fill=FALSE) + theme_minimal()
p <- p + theme(aspect.ratio=0.4)
p <- p + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6))
p <- p + theme(axis.title.y = element_text(angle=90))
p <- p + ggtitle(name)
p <- p + xlab("Distance to constraint")
p <- p + ylab("Delta Accessibility")
p
out <- paste("Accessibility_", file,".png",sep="")
ggsave(filename=out, path="./", width=7.2, height=7.2)


#PLot Correlation
field1 <- "Delta_acc"
field2 <- "Distance"
x <- pull(data,field1)
y <- pull(data,field2)

a <- as.character(round(cor(x=x, y=y, use="everything", method="pearson"), digits = 4))
b <- cor.test(x,y, use="everything", method="pearson")
d <- format.pval(round(b$p.value,digits=4))
c <- as.character(round(cor(x=x,y=y, use="everything", method="spearman"), digits = 4))
t <- cor.test(x=x,y=y, use="everything", method="spearman")
p <- format.pval(round(t$p.value,digits=4))
l <- bquote(paste("Spearman-R = ", .(c),"\n p-value = ",.(p),"\n Pearson-R = ",.(a),"\n p-value = ",.(d)))

p <- ggplot(data, aes(x=x, y=y, guides=FALSE)) + geom_point() + scale_colour_discrete(drop=TRUE,guide=FALSE) + theme_classic()
p <- p + theme(axis.text.x=element_text(angle=0, hjust=0.5, vjust=0, size = 18))
p <- p + theme(axis.text.y=element_text(angle=0, hjust=0, vjust=0.5, size = 18))
p <- p + theme(axis.title.y = element_text(angle=90,size=10))
p <- p + theme(axis.title.x = element_text(angle=0,size=10))
p <- p + xlab(field1)
p <- p + ylab(field2)
p <- p + theme(title = element_text(angle=0,size=10))
p <- p + labs(title=paste("Correlation",file,field1,field2,sep='_'))
p <- p + geom_text(data = data.frame(), aes(x=max(data[,field1])-(max(data[,field1])/5),y=min(data[,field2])+(max(data[,field2])/1.5), label=paste("Spearman-R = ",(c),"\n p-value = ",(p),"\n Pearson-R = ",(a),"\n p-value = ",(d)), colour = "black"))
out <- paste("Correlation",file,field1,field2,".png",sep="_")
ggsave(filename=out, path=".")
