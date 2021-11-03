library(ggplot2)
library(vroom)
library(plyr)
library(dplyr)
library(tidyverse)
library(stringr)
#library(scales)
#library(hexbin)

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
name <- args[2]
cores <- as.integer(args[3])

fn <- name
#fn <- paste(fn[3],fn[4],sep="_")

### Total output of CollectConstResults.py
#data <- vroom(file, delim="\t", col_names=c("Delta_acc", "Distance", "Acc_raw", "Zscore","Type"), col_types=c("Delta_acc"="n","Distance"="i","Acc_raw"="n","Zscore"="n","Type"="c"), num_threads=cores)
###Filtered output
data <- vroom(file, delim="\t", col_names=c("Constraint", "Distance", "Raw", "Scores"), col_types=c("Constraint"="c", "Distance"="i", "Raw"="n", "Scores"="c"), num_threads=cores)

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
    select(Constraint, Distance, Raw, starts_with("Scores_")) %>% rename(PreAcc = Scores_1, delta_Acc = Scores_2, zScore = Scores_3) %>% mutate(across(c(Distance, Raw, PreAcc, delta_Acc, zScore), as.numeric))


###Bin by Raw
toplot <- data %>% mutate(bin = cut(Raw, seq(0,1,0.1), right = FALSE)) %>% select(delta_Acc, bin)

# Density plot
p <- ggplot(toplot, aes(x=`delta_Acc`)) + geom_density() + theme_classic()
p <- p + theme(aspect.ratio=0.4)
p <- p + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6))
p <- p + theme(axis.title.y = element_text(angle=90))
p <- p + ggtitle(name)
p <- p + xlab("log Delta Accessibility")
p <- p + ylab("Density")
#p <- p + scale_x_log10()
p <- p + scale_x_continuous(trans = scales::pseudo_log_trans())
p <- p + scale_y_log10()
p <- p + facet_wrap(~ bin)
p
out <- paste("DeltaAccessibilityDensity_bin_Raw", fn,".png",sep="")
ggsave(filename=out, path="./", width=7.2, height=7.2)

###Bin by Relative Dist
toplot <- data %>% mutate(Distance = Distance/(-1*min(Distance))) %>%  mutate(bin = cut(Distance, seq(-1,1,0.05), right = FALSE)) %>% select(delta_Acc, bin)
    #mutate(bin = cut(Distance, 25, right = FALSE)) %>% select(delta_Acc, bin)

# Density plot
p <- ggplot(toplot, aes(x=`delta_Acc`)) + geom_density() + theme_classic()
p <- p + theme(aspect.ratio=0.4)
p <- p + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6))
p <- p + theme(axis.title.y = element_text(angle=90))
p <- p + ggtitle(name)
p <- p + xlab("log Delta Accessibility")
p <- p + ylab("Density")
#p <- p + scale_x_log10()
p <- p + scale_x_continuous(trans = scales::pseudo_log_trans())
p <- p + scale_y_log10()
p <- p + facet_wrap(~ bin)
p
out <- paste("DeltaAccessibilityDensity_bin_Dist", fn,".png",sep="")
ggsave(filename=out, path="./", width=7.2, height=7.2)

###Bin by Relative Constraint Position
toplot <- data %>% mutate(Constraint = as.numeric(str_split(.$Constraint, "\\-")[[1]][2])/1000) %>%  mutate(bin = cut(Constraint, 10, right = FALSE)) %>% select(delta_Acc, bin)

# Density plot
p <- ggplot(toplot, aes(x=`delta_Acc`)) + geom_density() + theme_classic()
p <- p + theme(aspect.ratio=0.4)
p <- p + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6))
p <- p + theme(axis.title.y = element_text(angle=90))
p <- p + ggtitle(name)
p <- p + xlab("log Delta Accessibility")
p <- p + ylab("Density")
#p <- p + scale_x_log10()
p <- p + scale_x_continuous(trans = scales::pseudo_log_trans())
p <- p + scale_y_log10()
p <- p + facet_wrap(~ bin)
p
out <- paste("DeltaAccessibilityDensity_bin_Constraint", fn,".png",sep="")
ggsave(filename=out, path="./", width=7.2, height=7.2)


###Bin by Zscore
toplot <- data %>% mutate(bin = cut(zScore, 10, right = FALSE)) %>% select(delta_Acc, bin)

# Density plot
p <- ggplot(toplot, aes(x=`delta_Acc`)) + geom_density() + theme_classic()
p <- p + theme(aspect.ratio=0.4)
p <- p + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6))
p <- p + theme(axis.title.y = element_text(angle=90))
p <- p + ggtitle(name)
p <- p + xlab("log Delta Accessibility")
p <- p + ylab("Density")
#p <- p + scale_x_log10()
p <- p + scale_x_continuous(trans = scales::pseudo_log_trans())
p <- p + scale_y_log10()
p <- p + facet_wrap(~ bin)
p
out <- paste("DeltaAccessibilityDensity_bin_zScore", fn,".png",sep="")
ggsave(filename=out, path="./", width=7.2, height=7.2)


### Plot Distribution
### 

###Bin by Relative Dist
toplot <- data %>% mutate(Distance = Distance/(-1*min(Distance))) %>%  mutate(bin = cut(Distance, seq(-1,1,0.1), right = FALSE)) %>% select(delta_Acc, bin)
#mutate(bin = cut(Distance, 25, right = FALSE)) %>% select(delta_Acc, bin)

p <- ggplot(toplot, aes(x=`bin`, y=`delta_Acc`)) + geom_violin(fill="bisque",color="black",alpha=0.3) + stat_summary(fun=mean, colour="black", geom="point", shape=18, size=3, show.legend = FALSE) + guides(fill=FALSE) + theme_classic()# + scale_fill_manual(values = c("firebrick"), alpha(.4)) + geom_jitter(aes(color='blue'),alpha=0.2)
p <- p + theme(aspect.ratio=0.4)
p <- p + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6))
p <- p + theme(axis.title.y = element_text(angle=90))
p <- p + ggtitle(name)
p <- p + xlab("Distance to constraint")
p <- p + ylab("Delta Accessibility")
p
out <- paste("Accessibility_", fn,".png",sep="")
ggsave(filename=out, path="./", width=7.2, height=7.2)


### Plot Correlation Distance-deltaAcc
## Only pos
field1 <- "Distance"
field2 <- "delta_Acc"
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
#p <- p + xlim(c(min(data$Distance),max(data$Distance)))
p <- p + xlab(field1)
p <- p + ylab(field2)
p <- p + theme(title = element_text(angle=0,size=10))
p <- p + labs(title=paste("Correlation",fn,field1,field2,sep='_'))
p <- p + geom_text(data = data.frame(), aes(x=max(x)-(max(x)/4),y=min(y)+(max(y)/1.5), label=paste("Spearman-R = ",(c),"\n p-value = ",(v),"\n Pearson-R = ",(a),"\n p-value = ",(d))), colour = "black")
p
out <- paste(paste("Correlation", fn, field1, field2, sep="_"), ".png",sep="")
ggsave(filename=out, path=".")

## Twosides
field1 <- "Distance"
field2 <- "delta_Acc"
x <- pull(data,field1) #%>% abs()
y <- pull(data,field2) #%>% abs()

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
p <- p + labs(title=paste("Correlation",fn,field1,field2,sep='_'))
p <- p + geom_text(data = data.frame(), aes(x=max(x)-(max(x)/4),y=min(y)+(max(y)/1.5), label=paste("Spearman-R = ",(c),"\n p-value = ",(v),"\n Pearson-R = ",(a),"\n p-value = ",(d))), colour = "black")
p
out <- paste(paste("Correlation_plusmin", fn, field1, field2, sep="_"), ".png",sep="")
ggsave(filename=out, path=".")


### Plot Correlation zScore-deltaAcc
### Only pos
field1 <- "zScore"
field2 <- "delta_Acc"
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
p <- p + labs(title=paste("Correlation",fn,field1,field2,sep='_'))
p <- p + geom_text(data = data.frame(), aes(x=max(x)-(max(x)/5),y=min(y)+(max(y)/1.5), label=paste("Spearman-R = ",(c),"\n p-value = ",(v),"\n Pearson-R = ",(a),"\n p-value = ",(d))), colour = "black")
p
out <- paste(paste("Correlation", fn, field1, field2, sep="_"), ".png",sep="")
ggsave(filename=out, path=".")

### Twosides
field1 <- "zScore"
field2 <- "delta_Acc"
x <- pull(data,field1) #%>% abs()
y <- pull(data,field2) #%>% abs()

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
p <- p + labs(title=paste("Correlation",fn,field1,field2,sep='_'))
p <- p + geom_text(data = data.frame(), aes(x=max(x)-(max(x)/5),y=min(y)+(max(y)/1.5), label=paste("Spearman-R = ",(c),"\n p-value = ",(v),"\n Pearson-R = ",(a),"\n p-value = ",(d))), colour = "black")
p
out <- paste(paste("Correlation_plusmin", fn, field1, field2, sep="_"), ".png",sep="")
ggsave(filename=out, path=".")


### Plot Correlation Raw-deltaAcc in relation to distance
### 
### 

###Bin by Relative Dist
toplot <- data %>% mutate(Distance = Distance/(-1*min(Distance))) %>%  mutate(bin = cut(Distance, seq(-1,1,0.1), right = FALSE)) 

cor_map <- toplot %>%
    split(., toplot$bin) %>%
    vapply(.,
           function(di) {
               cor(di$Raw, di$delta_Acc)
           }, numeric(1))

toplot$cor <- cor_map[toplot$bin]
toplot$content_and_cor <- paste(toplot$bin,
                                      toplot$cor)

field1 <- "Raw"
field2 <- "delta_Acc"

p <- ggplot(toplot, aes(x=Raw, y=delta_Acc)) + geom_point(color='blue', alpha=.1) + geom_rug(color='brown', alpha=.1) + theme_classic()
p <- p + theme(axis.text.x=element_text(angle=0, hjust=0.5, vjust=0, size = 5))
p <- p + theme(axis.text.y=element_text(angle=0, hjust=0, vjust=0.5, size = 5))
p <- p + theme(axis.title.y = element_text(angle=90,size=12))
p <- p + theme(axis.title.x = element_text(angle=0,size=12))
p <- p + xlab(field1)
p <- p + ylab(field2)
p <- p + theme(title = element_text(angle=0,size=10))
p <- p + facet_wrap(~content_and_cor, labeller = "label_value", ncol=3)
#p <- p + labs(title=paste("Correlation",fn,field1,field2,sep='_'))
#p <- p + geom_text(data = data.frame(), aes(x=max(x)-(max(x)/5),y=min(y)+(max(y)/1.5), label=paste("Spearman-R = ",(c),"\n p-value = ",(v),"\n Pearson-R = ",(a),"\n p-value = ",(d))), colour = "black")
p
out <- paste(paste("Correlation_binned_by_Dist", fn, field1, field2, sep="_"), ".png",sep="")
ggsave(filename=out, path=".")

### Plot Correlation Raw-deltaAcc in relation to zScore
### 
### 

###Bin by Relative Dist

toplot <- data %>% mutate(bin = cut(zScore, 10, right = FALSE))
#toplot <- data %>% mutate(zScore = zScore/max(zScore)) %>%  mutate(bin = cut(zScore, seq(-1,1,0.1), right = FALSE)) 

cor_map <- toplot %>%
    split(., toplot$bin) %>%
    vapply(.,
           function(di) {
               cor(di$Raw, di$delta_Acc)
           }, numeric(1))

toplot$cor <- cor_map[toplot$bin]
toplot$content_and_cor <- paste(toplot$bin,
                                toplot$cor)

field1 <- "Raw"
field2 <- "delta_Acc"

p <- ggplot(toplot, aes(x=Raw, y=delta_Acc)) + geom_point(color='blue', alpha=.1) + geom_rug(color='brown', alpha=.1) + theme_classic()
p <- p + theme(axis.text.x=element_text(angle=0, hjust=0.5, vjust=0, size = 5))
p <- p + theme(axis.text.y=element_text(angle=0, hjust=0, vjust=0.5, size = 5))
p <- p + theme(axis.title.y = element_text(angle=90,size=12))
p <- p + theme(axis.title.x = element_text(angle=0,size=12))
p <- p + xlab(field1)
p <- p + ylab(field2)
p <- p + theme(title = element_text(angle=0,size=10))
p <- p + facet_wrap(~content_and_cor, labeller = "label_value", ncol=3)
#p <- p + labs(title=paste("Correlation",fn,field1,field2,sep='_'))
#p <- p + geom_text(data = data.frame(), aes(x=max(x)-(max(x)/5),y=min(y)+(max(y)/1.5), label=paste("Spearman-R = ",(c),"\n p-value = ",(v),"\n Pearson-R = ",(a),"\n p-value = ",(d))), colour = "black")
p
out <- paste(paste("Correlation_binned_by_zScore", fn, field1, field2, sep="_"), ".png",sep="")
ggsave(filename=out, path=".")



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
