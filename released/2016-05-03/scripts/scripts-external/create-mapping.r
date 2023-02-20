#' generate a set of new ALNs for renaming ARIES
#' as well as a mapping between original ALNs and new 

load("third_release/data/samplesheet/data.Robj")
mapping <- data.frame(aln=unique(as.character(samplesheet$ALN)), stringsAsFactors=F)
set.seed(20160317)
mapping$new <- paste("ALN",sample(1:nrow(mapping)), sep="")
write.csv(mapping, file="mapping.csv",row.names=F)
