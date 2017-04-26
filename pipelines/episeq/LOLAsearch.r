# Search LOLA core collections for region DBs
#
# ARGUMENTS
#   basedir:       Root directory of region set database. Must be specified
#   genome:        Assembly to search. Default "hg19"
#   collection:    Collections to search. If not given, use all collections for the genome
# search terms. All are character vectors 
#   filename:      Will match any file containing the string in its name
#   description:   Will match any file where the description contains the string
#   any:           Matches any column containing the string
#   ...:           Any other column names you want to search against.
#   
#
# A region set will be used if it matches ANY of the patterns in ANY of the arguments 
#   filename, description, any, or `...`. If NONE of these arguments are given, ALL region
#   sets (in the correct genome and collection) will be used 
#
# Regular expressions are supported. Search is case insensitive
#
# returns a regionDB, a list with the following names: 
#   $collectionAnnos:     data.table. metadata for the collections
#   $regionAnnos:         data.table. metadata for the regions
#   $regionGRL:           GRangesList. Each GRanges object contains the regions from one bed file
#
# Example:
#   res <- LOLAsearch(colleciton=NA, description=c("histone", "ChIP"), filename="^199[1-9]?\\.bed$"
#                           antibody="AR", tissue=c("breast", "liver"), genome="hg19")


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
        dirs <- dirs[basename(dirs) %in% basenames]
    }
    collection <- basename(dirs)

    # get metadata for collections
    regionDB$collectionAnno <- readCollectionAnnotation(rootdir, collections=collection)

    # Search for any of the strings in lst in regionAnno[anno.col]]
    # Add rows containing matching the search pattern to the list of indeces
    # Return the updated list of indeces, maintaining uniqueness
    search_anno <- function(anno.col, lst) {
        if (length(lst) == 0) return(indeces)
        regex <- paste("(", lst, ")", sep="", collapse="|")
        ind <- which(grepl(regex, regionAnno[[anno.col]], ignore.case=TRUE))
        union(indeces, ind)
    }

    # Perform search if any terms specified. Otherwise, return entire collection(s)
    if ((length(filename) + length(description) + length(any) + length(additional)) == 0) {    
       regionSetAnnos <- readRegionSetAnnotation(rootdir, collections=collection)
    } else {
        regionSetAnnos <- data.table()
        for (collect in collection) {
            regionAnno <- readRegionSetAnnotation(rootdir, collection=collect)
            indeces <- c()
            indeces <- search_anno('filename', filename)
            indeces <- search_anno('description', description)
            for (name in colnames(regionAnno)) indeces <- search_anno(name, any)
            for (name in names(additional)) indeces <- search_anno(name, additional[[name]])
            regionSetAnnos <- rbind(regionSetAnnos, regionAnno[indeces, ], fill=TRUE)
        }
    } 

    # load region data
    regionDB$regionGRL <- readRegionGRL(rootdir, regionSetAnnos)
    regionDB$regionAnno <- regionSetAnnos

    return (regionDB)

}

