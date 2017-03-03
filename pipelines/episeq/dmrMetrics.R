library(reshape2)
library(ggplot2)
library(plyr)
library(knitr)

# usage:
# Rscript dmpMetrics.R dmr_file beta_file output_dir contrast_name

# Helper functions
get.freq <- function(by.column) {
    function(x) {
        data.frame(table(x[by.column], dnn=by.column))
    }
}
merge.by <- function(by.column){
    function(x, y){
        merge(x, y, by=by.column, all=TRUE)
    }
}
parse.input.list <- function(x) {
    x <- (strsplit(gsub('\\[', '', gsyb('\\]', '', input.files)), ', ')[[1]])
    unlist(lapply(x, function(x) { gsub("'", "", gsub('"', "", x)) }))
}
saveimg <- function(plot.name){
    filename <- paste(output.dir, "/", contrast.name, ".", plot.name, ".png", sep="")
    ggsave(filename, units="in", width=size, height=size)
}

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
pixels <- 800
size <- pixels / dpi

# dmrs by value
try({
    ggplot(dmrs, aes(x=value)) + 
    geom_density() + 
    scale_x_continuous(name="Average difference in methylation", breaks=seq(from=-1, to=1, by=0.25))
    saveimg("dmrs_by_value")
})

# dmrs by area
try({
    ggplot(dmrs, aes(x=area)) +
    geom_density() +
    scale_x_continuous(name="Area of differentially methylated region")
    saveimg("dmrs_by_area")
})

# dmrs by area - barplot
try({
    # Group largest 1% as one bin, and divide the rest of the data into 20 bins
    area <- dmrs$area
    len <- length(area)
    percentile99 <- area[order(area, decreasing=TRUE)][len / 100]
    binwidth <- round(percentile99 / 20)
    breaks <- seq(from = 1, to = binwidth * 20, by = binwidth)

    ggplot(dmrs, aes(x=area)) + geom_bar(binwidth = breaks) + scale_x_continuous(labels=c(0, breaks))
    saveimg("area_barplot1")

    ggplot(dmrs, aes(x=area)) + geom_bar() + scale_x_continuous(labels=c(0, breaks), breaks=breaks)
    saveimg("area_barplot2")
})
