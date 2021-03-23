library(ggplot2)
library(vroom)
library(plyr)
library(dplyr)
library(tidyverse)
library(scales)
library(hexbin)

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
name <- args[2]

#setwd('~/Work/TempAnalysis/denbicloud/RIssMed')
#file <- 'InfoSubset_20_7.tsv.gz'
#name <- 'GC20_Cons7'

data <- vroom(file, delim="\t", col_names=c("Delta_acc", "Distance", "Acc_raw", "Zscore","Type"), col_types=c("Delta_acc"="n","Distance"="i","Acc_raw"="n","Zscore"="n","Type"="c"), num_threads=6)
#data <- data %>% sample_frac(size=.01) %>% mutate(bin = cut_width(Distance, width = 10, center = 0, closed="left")) %>% group_by(bin)
data <- data %>% mutate(bin = cut_width(Distance, width = 10, center = 0, closed="left")) %>% group_by(bin)


p <- ggplot(data, aes(x=bin, y=Delta_acc)) + geom_violin(fill="bisque",color="black",alpha=0.3) + stat_summary(fun=mean, colour="black", geom="point", shape=18, size=3, show.legend = FALSE) + guides(fill=FALSE) + theme_minimal()# + scale_fill_manual(values = c("firebrick"), alpha(.4)) + geom_jitter(aes(color='blue'),alpha=0.2)
p <- p + theme(aspect.ratio=0.4)
p <- p + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6))
p <- p + theme(axis.title.y = element_text(angle=90))
p <- p + ggtitle(name)
p <- p + xlab("Distance to constraint")
p <- p + ylab("Delta Accessibility")
p
out <- paste("Accessibility_", file,".svg",sep="")
ggsave(filename=out, path="./", width=7.2, height=7.2)


#PLot Correlation
field1 <- "Distance"
field2 <- "Delta_acc"
x <- pull(data,field1) %>% abs()
y <- pull(data,field2) %>% abs()

a <- as.character(round(cor(x=x, y=y, use="everything", method="pearson"), digits = 4))
b <- cor.test(x,y, use="everything", method="pearson")
d <- format.pval(round(b$p.value,digits=4))
c <- as.character(round(cor(x=x,y=y, use="everything", method="spearman"), digits = 4))
t <- cor.test(x=x,y=y, use="everything", method="spearman")
v <- format.pval(round(t$p.value,digits=4))
l <- bquote(paste("Spearman-R = ", .(c),"\n p-value = ",.(v),"\n Pearson-R = ",.(a),"\n p-value = ",.(d)))

p <- ggplot(data, aes(x=x, y=y)) + geom_point(color='blue', alpha=.1) + geom_rug(color='brown', alpha=.1) + theme_classic()
p <- p + theme(axis.text.x=element_text(angle=0, hjust=0.5, vjust=0, size = 18))
p <- p + theme(axis.text.y=element_text(angle=0, hjust=0, vjust=0.5, size = 18))
p <- p + theme(axis.title.y = element_text(angle=90,size=10))
p <- p + theme(axis.title.x = element_text(angle=0,size=10))
p <- p + xlab(field1)
p <- p + ylab(field2)
p <- p + theme(title = element_text(angle=0,size=10))
p <- p + labs(title=paste("Correlation",file,field1,field2,sep='_'))
p <- p + geom_text(data = data.frame(), aes(x=max(x)-(max(x)/5),y=min(y)+(max(y)/1.5), label=paste("Spearman-R = ",(c),"\n p-value = ",(v),"\n Pearson-R = ",(a),"\n p-value = ",(d))), colour = "black")
p
out <- paste(paste("Correlation", file, field1, field2, sep="_"), ".svg",sep="")
ggsave(filename=out, path=".")
