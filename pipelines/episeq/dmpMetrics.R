library(ggplot2)
library(gplots)
library(plyr)
library(knitr)

# usage:
# Rscript dmpMetrics.R dmp_file beta_file [cases] [controls] output_dir contrast_name

# Helper functions
parse.input.list <- function(x) {
    x <- (strsplit(gsub('\\[', '', gsub('\\]', '', x)), ', ')[[1]])
    unlist(lapply(x, function(x) { gsub("'", "", gsub('"', "", x)) }))
}
saveimg <- function(plot.name){
    filename <- paste(output.dir, "/", contrast.name, ".", plot.name, ".png", sep="")
    ggsave(filename, units="in", width=size, height=size)
} 

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

# image parameters
dpi <- 300
pixels <- 800
size <- pixels / dpi

# dmps by average delta beta
try({
    ggplot(dmps, aes(x=Avg.Delta.Beta)) + 
    geom_density() +
    scale_x_continuous(name="Average Delta Beta", breaks=seq(from=-1, to=1, by=0.25))
    saveimg("dmps_by_avg_delta_beta")
})

# Beta value pca - all points
try({
    global.pcadata <- as.data.frame(prcomp(t(na.omit(betas)))$x)
    ggplot(global.pcadata, aes(x=PC1, y=PC2)) + geom_point()
    saveimg("global_beta_pca")
})

# Beta value pca - dmps
try({
    dmp.betas <- dmps[c(cases, controls)]
    dmp.pcadata <- as.data.frame(prcomp(t(dmp.betas))$x)
    ggplot(dmp.pcadata, aes(x=PC1, y=PC2)) + geom_point()
    saveimg("dmp_beta_pca")
})

# Beta value heatmap
try({
    num.points <- 300
    sort <- order(dmps$Avg.Delta.Beta, decreasing=TRUE)[1:num.points]
    most.variation <- dmps[sort, ][c(cases, controls)]
    row.labels <- dmps[sort, ]$Row.names
    colors <- ifelse(colnames(most.variation) %in% cases, "red", "green")
    png(paste(output.dir, "/", contrast.name, ".", "dmp_heatmap.png", sep=""), width=700, height=1200, margins=(10,5))
    heatmap.2(as.matrix(most.variation), key=TRUE, col=colorRampPalette(c('red', 'yellow', 'blue')), colCol=colors, trace='none', breaks=21, labRow=row.labels)
    dev.off()
})

# metrics
try({
    metrics <- list(total.pos=length(betas))
    metrics$filtered.pos <- length(dmps)
    metrics$percent.kept <- (metrics$filtered.pos / metrics$total.pos) * 100
    metrics <- as.data.frame(metrics)

    # append to metrics file
    dmp.metrics.file <- paste(output.dir, "dmp.metrics.csv", sep="/")
    dmp.metrics.table <- paste(output.dir, "dmp.metrics.table", sep="/")
    if (file.exists(dmp.metrics.file)) {
        metrics.file <- read.csv(dmp.metrics.file)
        metrics <- rbind(metrics, metrics.file)
    }
    cat(kable(metrics), file=dmp.metrics.table, sep="\n")
})
