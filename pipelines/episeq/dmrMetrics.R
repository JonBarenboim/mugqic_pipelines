library(reshape2)
library(ggplot2)
library(plyr)
library(knitr)

# usage:
# Rscript dmpMetrics.R dmr_file output_dir contrast_name

# Helper functions
saveimg <- function(plot.name, width=size, height=size){
    filename <- paste(output.dir, "/", contrast.name, ".", plot.name, ".png", sep="")
    ggsave(filename, units="in", width=width, height=height)
}

# Set font size
theme_set(theme_grey(base_size=4) + theme(title=element_text(size=rel(1.1))))

# Read arguments
args = commandArgs(trailingOnly=TRUE)
dmr.file <- args[1]
output.dir <- args[2]
contrast.name <- args[3]

# Read in data
dmrs <- read.csv(dmr.file)

# image parameters
dpi <- 300
pixels <- 600
size <- pixels / dpi

# dmrs by value
ggplot(dmrs, aes(x=value)) + 
    geom_histogram(binwidth=0.02) + 
    scale_x_continuous(name="Average difference in methylation", breaks=seq(from=-1, to=1, by=0.25)) +
    ggtitle("DMRs by Average Difference in Methylation")
saveimg("dmrs_by_value")

# dmrs by area
break_func <- function(limits) {
    limits <- c(floor(limits[1]), ceiling(limits[2]))
    min <- if(limits[1] %% 2 == 0) limits[1] else limits[1] - 1
    max <- if(limits[2] %% 2 == 0) limits[2] else limits[2] + 1
    seq(from=min, to=max, by=2)
}
ggplot(dmrs, aes(x=area)) + 
    geom_histogram(binwidth=0.2) +
    scale_x_continuous(breaks=break_func) + 
    ggtitle("DMRs by Area of the Bump")
saveimg("dmrs_by_area")

