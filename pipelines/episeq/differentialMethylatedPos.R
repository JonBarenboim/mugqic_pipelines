library(reshape2)
library(ggplot2)
library(plyr)

# usage:
# -i [input.files] or -i input_file
# -o output.basename


# Helper functions
remove.quotes <- function(x) {
    gsub("'", "", gsub('"', "", x))
}

get.contrast <- function(x) {
    gsub("_RRBS_differential_methylated_pos.csv", "", basename(x))
}

get.freq <- function(by.column){
    function(x){
        data.frame(table(seqnames=x$seqnames))
    }
}

merge.by <- function(by.column) {
    function(x, y) {
        merge(x, y, by=by.column, all=TRUE)
    }
}

#parse args
args = commandArgs()
for (i in seq_along(args)) {
    if (args[i] == "-i") { 
        input.files = args[i+1]
    } else if (args[i] == "-o") { 
        output.basename = args[i+1]
    }
}
if (substr(input.files, 1, 1) == "[") {
   input.files <- strsplit(gsub("\\[", "", gsub("\\]", "", input.files)), ", " )[[1]]
   input.files <- lapply(input.files, remove.quotes)
}

# parameters

# Read in data
csvs <- lapply(input.files, read.csv)
names(csvs) <- lapply(input.files, get.contrast)

#dmps by chromosome
freqs <- lapply(csvs, get.freq("seqnames"))
freqs <- Reduce(merge.by("seqnames"), freqs)
freqs[is.na(freqs)] <- 0
freqs.melt <- melt(freqs, id.var="seqnames")
ggplot(freqs.melt, aes(x=variable, y=value, xlab="contrast", ylab="number of dmps", fill=seqnames)) +
    geom_bar(stat="identity", position="stack")
ggsave(paste(output.basename, "dmps_by_sequence.png", sep="/"), dpi=100)


# dmps by delta beta
deltas <- ldply(csvs, rbind)
deltas <- deltas[c('.id', 'Row.names', 'Avg.Delta.Beta')]
ggplot(deltas, aes(x=Avg.Delta.Beta, fill=.id)) +
    geom_histogram(position="dodge", binwidth=0.25, center=0.125) +
    scale_x_continuous(name="Average Delta Beta", breaks=seq(from=-1, to=1, by=0.25))
ggsave(paste(output.basename,"dmps_by_avg_delta_beta.png", sep="/"), dpi=100)


#deltas <- lapply(csvs, function(x){data.frame(table(bin=cut(x$Avg.Delta.Beta, seq(from=-1, to=1, by=0.25))))})
#deltas <- Reduce (function(x,y){merge(x, y, by="bin", all=True)}, deltas)
#deltas.melt <- melt(deltas, id.var="bin")
#ggplot(deltas.melt, aes(x=bin, y=value, xlab="Avg Delta Beta", ylab="Frequency", fill="variable")) + geom_bar(stat="identity", position="dodge")
#    geom_histogram(position="dodge", stat="identity", binwidth=0.25, center=0.125)
#ggsave(paste(output.basename, "dmps_by_avg_delta_beta.pdf", sep="/"))




#deltas <- Reduce(function(x,y){cbind(names(x)=x$Avg.Delta.Beta, names(y)=y$Avg.Delta.Beta)})
#deltas.melt <- melt(deltas, id.var="Row.names")
#ggplot(deltas) + geom_histogram(binwidth=0.25, center=0.125
#ggplot(freqs, aes(x="variable", y="value"))
#ggsave(paste(output.basename, "dmps_by_avg_delta_beta.pdf", sep="/"))
                                      
