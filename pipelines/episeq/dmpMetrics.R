library(ggplot2)
library(gplots)
library(plyr)
library(knitr)

# usage:
# Rscript dmpMetrics.R dmp_file beta_file [cases] [controls] output_dir contrast_name

# Helper functions
parse.input.list <- function(x) {
    x <- (strsplit(gsub('\\[', '', gsub('\\]', '', x)), ', ')[[1]])
    sapply(x, function(x) { gsub("'", "", gsub('"', "", x)) })
}
saveimg <- function(plot.name, width=size, height=size){
    filename <- paste(output.dir, "/", contrast.name, ".", plot.name, ".png", sep="")
    ggsave(filename, units="in", width=width, height=height)
} 

# Set font size
theme_set(theme_grey(base_size=4) + theme(title=element_text(size=rel(1.1))))

# Parse arguments
args = commandArgs(trailingOnly=TRUE)
dmp.file <- args[1]
beta.file <- args[2]
cases <- parse.input.list(args[3])
controls <- parse.input.list(args[4])
output.dir <- args[5]
contrast.name <- args[6]

# Read in data
dmps <- read.csv(dmp.file)
betas <- read.csv(beta.file)
rownames(betas) <- betas$X

# image parameters
dpi <- 300
pixels <- 600
size <- pixels / dpi

# dmps by average delta beta
ggplot(dmps, aes(x=Avg.Delta.Beta)) + 
    geom_histogram(binwidth=0.05) + 
    scale_x_continuous(name="Average Delta Beta", breaks=seq(from=-1, to=1, by=0.2)) + 
    ggtitle("DMPs by Average Delta Beta")
saveimg("dmps_by_avg_delta_beta")
    
# Beta value pca - all points
global.pcadata <- as.data.frame(prcomp(t(na.omit(betas[c(controls, cases)])))$x)
global.pcadata$state <- ifelse(rownames(global.pcadata) %in% cases, "case", "control")
ggplot(global.pcadata, aes(x=PC1, y=PC2, col=state)) + 
    geom_point() + 
    theme(legend.key.size=unit(0.2, "cm")) +
    ggtitle("PCA of Methylation Values - All Positions")
saveimg("global_beta_pca", height=size*0.75)

# Beta value pca - dmps
dmp.betas <- betas[as.character(dmps$Row.names), c(controls, cases)]
dmp.pcadata <- as.data.frame(prcomp(t(dmp.betas))$x)
dmp.pcadata$state <- ifelse(rownames(dmp.pcadata) %in% cases, "case", "control")
ggplot(dmp.pcadata, aes(x=PC1, y=PC2, col=state)) + 
    geom_point() + 
    theme(legend.key.size=unit(0.2, "cm")) +
    ggtitle("PCA of Methylation Values - DMPs Only")
saveimg("dmp_beta_pca", height=size*0.75)

# Beta value heatmap
num.points <- 75
sort <- order(dmps$Avg.Delta.Beta, decreasing=TRUE)[1:num.points]
ind <- dmps[sort, 'Row.names']
most.variation <- betas[as.character(ind), c(cases,controls)]
# row.labels <- dmps[sort, ]$Row.names
row.labels <- paste(dmps[sort, 'seqnames'], ':', dmps[sort, 'start'], sep="")
colors <- ifelse(colnames(most.variation) %in% cases, "red", "green")
png(paste(output.dir, "/", contrast.name, ".", "dmp_heatmap.png", sep=""), width=700, height=1000)
heatmap.2(as.matrix(most.variation), key=TRUE, col=colorRampPalette(c('red', 'yellow', 'blue')), 
          colCol=colors, trace='none', breaks=21, labRow=row.labels, margins=c(10,5), 
          main=paste("Methylation Values of", num.points, "most Differential Positions", sep=" "))
dev.off()

# metrics
metrics <- list(total.pos=nrow(betas))
metrics$num.dmp <- nrow(dmps)
metrics$percent.differential <- (metrics$num.dmp / metrics$total.pos) * 100
metrics <- as.data.frame(metrics)
rownames(metrics) <- c(contrast.name)

# append to metrics file
dmp.metrics.file <- paste(output.dir, "dmp.metrics.csv", sep="/")
dmp.metrics.table <- paste(output.dir, "dmp.metrics.table", sep="/")
if (file.exists(dmp.metrics.file)) {
    metrics.file <- read.csv(dmp.metrics.file, row.names=1)
    metrics <- rbind(metrics, metrics.file)
}
write.csv(metrics, file=dmp.metrics.file)
cat(kable(metrics, row.names=TRUE), file=dmp.metrics.table, sep="\n")
