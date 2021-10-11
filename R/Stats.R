library(ggplot2)
library(vroom)
library(plyr)
library(dplyr)
library(tidyverse)
#library(scales)
#library(hexbin)

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
name <- args[2]
cores <- args[3]

#setwd('~/Work/TempAnalysis/denbicloud/RIssMed')
#file <- 'InfoSubset_20_7.tsv.gz'
#name <- 'GC20_Cons7'

### Total output of CollectConstResults.py
#data <- vroom(file, delim="\t", col_names=c("Delta_acc", "Distance", "Acc_raw", "Zscore","Type"), col_types=c("Delta_acc"="n","Distance"="i","Acc_raw"="n","Zscore"="n","Type"="c"), num_threads=cores)
###Filtered output
data <- vroom(file, delim="\t", col_names=c("Constraint", "Distance", "Scores"), col_types=c("Constraint"="c", "Distance"="i","Scores"="c"), num_threads=cores)

#https://stackoverflow.com/questions/4350440/split-data-frame-string-column-into-multiple-columns

split_into_multiple <- function(column, pattern = ", ", into_prefix){
    cols <- str_split_fixed(column, pattern, n = Inf)
    # Sub out the ""'s returned by filling the matrix to the right, with NAs which are useful
    cols[which(cols == "")] <- NA
    cols <- as.tibble(cols)
    # name the 'cols' tibble as 'into_prefix_1', 'into_prefix_2', ..., 'into_prefix_m' 
    # where m = # columns of 'cols'
    m <- dim(cols)[2]
    
    names(cols) <- paste(into_prefix, 1:m, sep = "_")
    return(cols)
}

data <- data %>% 
    bind_cols(split_into_multiple(.$Scores, "\\|", "Scores")) %>% 
    # selecting those that start with 'type_' will remove the original 'type' column
    select(Constraint, Distance, starts_with("Scores_")) %>% rename(PreAcc = Scores_1, delta_Acc = Scores_2, zScore = Scores_3)


# Density plot
p <- ggplot(data, aes(x=`delta_Acc`)) + geom_density() + theme_minimal()
p <- p + theme(aspect.ratio=0.4)
p <- p + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6))
p <- p + theme(axis.title.y = element_text(angle=90))
p <- p + ggtitle(name)
p <- p + xlab("Delta Accessibility")
p <- p + ylab("Density")
#p <- p + scale_x_log10()
p
out <- paste("DeltaAccessibilityDensity_", file,".svg",sep="")
ggsave(filename=out, path="./", width=7.2, height=7.2)


p <- ggplot(data, aes(x=`Acc_raw`)) + geom_density() + theme_minimal()
p <- p + theme(aspect.ratio=0.4)
p <- p + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6))
p <- p + theme(axis.title.y = element_text(angle=90))
p <- p + ggtitle(name)
p <- p + xlab("Accessibility")
p <- p + ylab("Density")
p <- p + scale_x_log10()
p
out <- paste("AccessibilityDensity_", file,".svg",sep="")
ggsave(filename=out, path="./", width=7.2, height=7.2)

p <- ggplot(data, aes(x=log10(`Acc_raw`))) + stat_ecdf(geom = "step", pad = FALSE) + theme_minimal()
p <- p + theme(aspect.ratio=0.4)
p <- p + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6))
p <- p + theme(axis.title.y = element_text(angle=90))
p <- p + ggtitle(name)
p <- p + xlab("Accessibility")
p <- p + ylab("Density")
p
out <- paste("AccessibilityECDF_", file,".svg",sep="")
ggsave(filename=out, path="./", width=7.2, height=7.2)



# Binning
#data <- data %>% mutate(bin = cut_width(Distance, width = 10, center = 0, closed="left")) %>% group_by(bin)
#
#p <- ggplot(data, aes(x=`bin`, y=`Delta_acc`)) + geom_violin(fill="bisque",color="black",alpha=0.3) + stat_summary(fun=mean, colour="black", geom="point", shape=18, size=3, show.legend = FALSE) + guides(fill=FALSE) + theme_minimal()# + scale_fill_manual(values = c("firebrick"), alpha(.4)) + geom_jitter(aes(color='blue'),alpha=0.2)
#p <- p + theme(aspect.ratio=0.4)
#p <- p + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6))
#p <- p + theme(axis.title.y = element_text(angle=90))
#p <- p + ggtitle(name)
#p <- p + xlab("Distance to constraint")
#p <- p + ylab("Delta Accessibility")
#p
#out <- paste("Accessibility_", file,".svg",sep="")
#ggsave(filename=out, path="./", width=7.2, height=7.2)
#
#
##PLot Correlation
#field1 <- "Distance"
#field2 <- "Delta_acc"
#x <- pull(data,field1) %>% abs()
#y <- pull(data,field2) %>% abs()
#
#a <- as.character(round(cor(x=x, y=y, use="everything", method="pearson"), digits = 4))
#b <- cor.test(x,y, use="everything", method="pearson")
#d <- format.pval(round(b$p.value,digits=4))
#c <- as.character(round(cor(x=x,y=y, use="everything", method="spearman"), digits = 4))
#t <- cor.test(x=x,y=y, use="everything", method="spearman")
#v <- format.pval(round(t$p.value,digits=4))
#l <- bquote(paste("Spearman-R = ", .(c),"\n p-value = ",.(v),"\n Pearson-R = ",.(a),"\n p-value = ",.(d)))
#
#p <- ggplot(data, aes(x=x, y=y)) + geom_point(color='blue', alpha=.1) + geom_rug(color='brown', alpha=.1) + theme_classic()
#p <- p + theme(axis.text.x=element_text(angle=0, hjust=0.5, vjust=0, size = 18))
#p <- p + theme(axis.text.y=element_text(angle=0, hjust=0, vjust=0.5, size = 18))
#p <- p + theme(axis.title.y = element_text(angle=90,size=10))
#p <- p + theme(axis.title.x = element_text(angle=0,size=10))
#p <- p + xlab(field1)
#p <- p + ylab(field2)
#p <- p + theme(title = element_text(angle=0,size=10))
#p <- p + labs(title=paste("Correlation",file,field1,field2,sep='_'))
#p <- p + geom_text(data = data.frame(), aes(x=max(x)-(max(x)/5),y=min(y)+(max(y)/1.5), label=paste("Spearman-R = ",(c),"\n p-value = ",(v),"\n Pearson-R = ",(a),"\n p-value = ",(d))), colour = "black")
#p
#out <- paste(paste("Correlation", file, field1, field2, sep="_"), ".png",sep="")
#ggsave(filename=out, path=".")
