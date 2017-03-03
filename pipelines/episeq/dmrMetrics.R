library(reshape2)
library(ggplot2)
library(plyr)
library(knitr)

# usage:
# Rscript dmpMetrics.R dmr_file beta_file output_dir contrast_name

# Helper functions
parse.input.list <- function(x) {
    x <- (strsplit(gsub('\\[', '', gsyb('\\]', '', input.files)), ', ')[[1]])
    unlist(lapply(x, function(x) { gsub("'", "", gsub('"', "", x)) }))
}
saveimg <- function(plot.name){
    filename <- paste(output.dir, "/", contrast.name, ".", plot.name, ".png", sep="")
    ggsave(filename, units="in", width=size, height=size)
}

# Set font size
theme_set(theme_grey(base_size=5))

# Parse arguments
args = commandArgs(trailingOnly=TRUE)
dmr.file <- args[1]
beta.file <- args[2]
output.dir <- args[3]
contrast.name <- args[4]

# Read in data
dmrs <- read.csv(dmr.file)
betas <- read.csv(beta.file)

# image parameters
dpi <- 300
pixels <- 600
size <- pixels / dpi

# dmrs by value
ggplot(dmrs, aes(x=value)) + 
    geom_density() + 
    scale_x_continuous(name="Average difference in methylation", breaks=seq(from=-1, to=1, by=0.25)) +
    ggtitle("Differentially Methylated Regions by Average Difference in Methylation")
saveimg("dmrs_by_value")

# dmrs by area
# Group largest 1% as one bin, and divide the rest of the data into 20 bins
area <- dmrs$area
len <- length(area)
percentile99 <- area[order(area, decreasing=TRUE)][len / 100]
binwidth <- round(percentile99 / 20, digits=1)
breaks <- seq(from = 0, to = binwidth * 20, by = binwidth)
breaks <- c(breaks, ceiling(max(area)))
ggplot(dmrs, aes(x=area)) + 
    geom_histogram(binwidth=binwidth) +
    ggplot("Differentially Methylated Regions by Area of the Bump")
saveimg("dmrs_by_area")

###### DMRS BY AREA - VARIABLE BIN WIDTH !?!?!?!?
# area <- dmrs$area
# len <- length(area)
# percentile99 <- area[order(area, decreasing=TRUE)][len / 100]
# binwidth <- round(percentile99 / 10, digits=1)
# breaks <- c(seq(from=0, to=binwidth*10, by=binwidth), ceiling(max(area)))
# dmrs$bin <- cut(dmrs$area, breaks)
# ggplot(dmrs, aes(x=bin)) +
#     geom_histogram(stat='count') +
#     ggtitle("Differentially Methylated Regions by Area of the Bump")
# saveimg("dmrs_by_area")
