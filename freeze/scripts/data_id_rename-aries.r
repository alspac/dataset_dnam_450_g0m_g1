
args = commandArgs(trailingOnly=TRUE)


#  timecode = NA as we want all time points.

tc = NA

datdir = args[1]
outdir <- args[2]
mapping <- read.csv(args[3])
woc_file <- args[4]

names(mapping)[1] <- "aln"
names(mapping)[16] <- "new" ### Change to id col we want

##remove consent
consent <- read.table(woc_file,h=T,sep=",")
withdrawnConsent <- unique(paste(consent$ALN))

## identify and remove the overlap between the mapping file and withdrawnConsent
a <- which(mapping$aln %in% withdrawnConsent)
randomID <- paste("Sample",sample(c(1:100000),length(a)),sep="_")
mapping$new[a] <- paste(randomID)


#' 3. create a copy of ARIES with the new identifiers (about 1 hour)
source("./rename-aries-function.r")
rename.aries(mapping, datdir, outdir,
             time.codes = tc ) # if time.codes not set, include all samples;
                                        # otherwise, set to a subset of
                                        # antenatal,cord,F17,F7,FOM,TF1-3,TF3,15up
					# children = cord, F7, TF3, F17 (time_codes)
					# 15up is a "time_point" that is equal to TF3 and F17
					# mothers = antenatal, FOM, TF1-3

#' 4 (optional). check that the renaming was performed correctly
## NOTE: this script does not work when subsetting !
source("./compare-releases-function.r")
compare.releases(mapping, datdir , outdir)

