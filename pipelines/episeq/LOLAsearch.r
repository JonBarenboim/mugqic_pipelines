# Search LOLA core collections for region DBs
#   collections: collections to search. 
#   filename: a list or vector of strings. Will match any file containing the string
#   description: a list or vector of strings. Will match any file where the description contains the string
#   any: a list or vector of strings. Matches any column containing the string
#   ...: any other column names you want to search against. Same format as above arguments
#   genome: default "hg19"
#   basedir: rott directory of genomes. Must be specified
#   All matching is "or".
#
# Returns a list of filenames with additional metadata
#
# Example:
#   LOLAsearch(collection=NA, 
#              description=c("histone", "ChIP"), 
#              additional=list(tissue=c("breast", "liver"), antibody=c("AR"), 
#              genome="hg19"))
#   

# Creates a regex pattern that matches any of the strings in lst
makeRegex <- function(lst) {
    paste(sapply(lst, function(x) paste("(", x, ")", sep="")), collapse="|")
}


#returns a regionDB, a list with names:
# $collectionAnnos:     data.table. metadata for the collections
# $regionAnnos:         data.table. metadata for the regions
# $regionGRL:           GRangesList. Each GRanges object the regions from one bed file

buildRegionDB <- function(collection=list(), filename=list(), description=list(), any=list(),
                          genome="hg19", rootdir=stop("directory must be specified"), ...) {
    require(LOLA)
    require(data.table)

    additional <- list(...)
    collection <- collection[!is.na(collection)]
    filename <- filename[!is.na(filename)]
    description <- description[!is.na(description)]
    any <- any[!is.na(any)]

    if (basename(rootdir) != genome) rootdir <- file.path(rootdir, genome)
    if (!dir.exists(rootdir)) stop(paste(rootdir, "is not a valid directory"))

    regionDB <- list()

    # find collections to be searched
    dirs <- list.dirs(rootdir, recursive=FALSE)
    if (length(collection) > 0) {
        basenames <- intersect(basename(dirs), collection)
        if(length(setdiff(collection, basenames)) > 0) {
            warning(paste("The following collections could not be found:",
                          paste(setdiff(collection, basenames), collapse=", ")))
        }
        if (length (basenames) == 0){
            stop("no valid collections found")
        }
        dirs <- dirs[which(basename(dirs) %in% basenames)]
    }
    collection <- basename(dirs)

    # get metadata for collections
    regionDB$collectionAnno <- readCollectionAnnotation(rootdir, collections=collection)

    regionSetAnnos <- data.table()

    # Perform search if any terms specified. Otherwise, return entire collection(s) 
    if ((length(filename) + length(description) + length(any) + length(additional)) > 0) {
        for (collect in collection) {
            regionAnno <- readRegionSetAnnotation(rootdir, collections=collect)
            indeces <- c()
            
            # search filenames
            if (length(filename) > 0) {
                ind <- which(grepl(makeRegex(filename), regionAnno$filename, ignore.case=TRUE))
                indeces <- union(indeces, ind)
            }

            # search description
            if (length(description) > 0) {
                ind <- which(grepl(makeRegex(description), regionAnno$description, ignore.case=TRUE))
                indeces <- union(indeces, ind)
            }

            # search all other columns
            if (length(any) > 0) {
                for (name in colnames(regionAnno)){
                    ind <- which(grepl(makeRegex(any), regionAnno[[name]], ignore.case=TRUE))
                    indeces <- union(indeces, ind)
                }
            }

            # search additional columns by name
            for (name in names(additional)) {
                if (length(additional[[name]]) == 0) next 
                ind <- which(grepl(makeRegex(additional[[name]]), regionAnno[[name]], ignore.case=TRUE))
                indeces <- union(indeces, ind)
            }

            regionSetAnnos <- rbind(regionSetAnnos, regionAnno[indeces, ], fill=TRUE)
        }
    } else {
        regionSetAnnos <- readRegionSetAnnotation(rootdir, collections=collection)
    }

    # load region data
    regionDB$regionGRL <- readRegionGRL(rootdir, regionSetAnnos)
    regionDB$regionAnno <- regionSetAnnos

    return (regionDB)

}

LOLAsearch <- function(collections=list(), filename=list(), description=list(), 
                       additional=list(), any= list(),  genome="hg19", 
                       basedir=stop("directory must be specified")) {

    collections <- collections[!is.na(collections)]
    filename <- filename[!is.na(filename)]
    description <- description[!is.na(description)]
    additional <- additional[!is.na(additional)]

    root <- file.path(basedir, genome)
    if (!dir.exists(root)) stop (paste(root, "is not a valid directory"))

    # find collections to be searched
    dirs <- list.dirs(root, recursive=FALSE)
    if (length(collections) > 0){
        basenames <- intersect(basename(dirs), collections)
        if(length(setdiff(collections, basenames)) > 0) {
            warning(paste("The following collections could not be found:", 
                          paste(setdiff(collections, basenames), collapse=", ")))
        }
        if (length(basenames) == 0) {
            stop("no valid collections found")
        }
        dirs <- dirs[which(basename(dirs) %in% basenames)]
    }

    results=list()
    for (d in dirs) {
        d.results <- c() 
        # FIXME    delimiter is not standard! Usually tab, but sometimes commas!
        load(paste(file.path(d, paste(basename(d), "_files.RData"))))
        indedx <- ret
        load(paste(file.path(d, paste(basename(d), ".RData"))))
        DBlist <- ret

        # search filenames
        if(!(length(filename) == 0 || is.na(filename))){
            ind <- which(grepl(makeRegex(filename), index[['filename', exact=FALSE]]))
            d.results <- union(d.results, ind)
        }

        # search description
        if(!(length(description) == 0 || is.na(description))){
            ind <- which(grepl(makeRegex(description), index[['description', exact=FALSE]]))
            d.results <- union(d.results, ind)
        }

        # search additional columns
        for (name in names(additional)) {
            strs <- additional[[name]]
            ind <- which(grepl(makeRegex(strs), index[[name, exact=FALSE]]))
            d.results <- union(d.results, ind)
        }


        # Append to results
        results[[d]] <- index[d.results, ]

    }

    # remove collections with no results
    results <- Filter(function(x) { nrow(x) > 0 }, results)

    # Add a 'collection' column and reshape into one data frame 
    results <- Map(function(dat, name) {dat$collection <- basename(name); dat}, results, names(results))
    results <- rbind.fill(results)

    # order columns
    ordered <- c("collection", "filename", "description")
    unordered <- setdiff(colnames(results), ordered)
    if (! "description" %in% colnames(results)) results$description <- NA
    return (results[c(ordered, unordered)])

}
