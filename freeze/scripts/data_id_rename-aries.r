args = commandArgs(trailingOnly=TRUE)

datdir = args[1]
outdir <- args[2]
mapping <- read.csv(args[3])
woc_file <- args[4]

names(mapping)[which(names(mapping)=="aln")] <- "aln"
names(mapping)[which(names(mapping)=="dnam_450_g0m_g1")] <- "new"

##remove non-consenting
consent <- read.table(woc_file,h=T,sep=",")
mapping <- mapping[!(as.character(mapping$aln) %in% as.character(consent$ALN)),]

## create a copy of ARIES with the new identifiers (about 1 hour)
source("./rename-aries-function.r")
rename.aries(mapping, datdir, outdir)
  ## takes about 20 minutes to run

## check that the renaming was performed correctly
source("./compare-releases-function.r")
compare.releases(mapping, datdir, outdir)
  ## takes about 5 minutes to run
