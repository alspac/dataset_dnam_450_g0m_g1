#' 1. Copy the dataset to scratch space or some place accessible to a compute node.
#' !!! Do not attempt to run this on a bluecrystal login node!  Large data files need to be loaded. !!!
#' e.g.
#' cp -rv /projects/ALSPAC/studies/originals/alspac/epigenetic/methylation/450k/aries/dev/third_release .

#' 2. load ALN mapping 
mapping <- read.csv("mapping.csv")
#mapping$aln <- mapping$...
#mapping$new <- mapping$...

#' 3. create a copy of ARIES with the new identifiers (about 1 hour)
source("rename-aries-function.r")
rename.aries(mapping,
             "third_release/data",
             "external_release")
## ... copy documentation as well (remove lines with email addresses or "Josine")
## grep -v

#' 4 (optional). check that the renaming was performed correctly
source("compare-releases-function.r")
compare.releases(mapping,
                 "third_release/data",
                 "external_release")

#' 5. move the resulting dataset to its final destination
#' e.g.
#' cp -rv external_release /projects/ALSPAC/studies/originals/alspac/epigenetic/methylation/450k/aries/dev/
