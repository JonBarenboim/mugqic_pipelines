library(ggplot2)
library(minfi)


# usage:
# Rscript betaMetrics.R beta_file output_dir

# image parameters
pixels <- 800

# read arguments
args = commandArgs(trailingOnly=TRUE)
beta.file <- args[1]
output.dir <- args[2]

# read in data
beta <- read.csv(beta.file)

png(paste(output.dir, 'beta_beanplot.png', sep='/'), units="px", width=pixels, height=pixels)
densityBeanPlot(as.matrix(beta)) 
dev.off()
