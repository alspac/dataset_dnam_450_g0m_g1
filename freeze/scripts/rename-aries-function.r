#' Remove identifying information from ARIES
#'
#' Make a copy of the ARIES dataset with reassigned ALNs.
#'
#' @param mapping Data frame mapping ALSPAC ALNs (column 'aln') to
#' new identifiers (column 'new').
#' @param data.dir The directory containing the dataset.
#' @param output.dir The directory in which to place the modified dataset.
#' @return The samplesheet data frame for the new dataset.
#'
#' Upon completion, the directory denoted by 'output.dir' will
#' contain a set of files identical to that of 'data.dir'
#' but with all ALNs replaced with new identifiers
#' and any other identifying information
#' such as SNPs removed.
#'
#' @examples \dontrun{
#' ## create a sample mapping
#' load("/projects/ALSPAC/studies/originals/alspac/epigenetic/methylation/450k/aries/dev/third_release/data/samplesheet/data.Robj")
#' mapping <- data.frame(aln=samplesheet$ALN, new=paste("s",1:nrow(samplesheet), sep=""))
#' write.csv(mapping, file="mapping.csv",row.names=F)
#'
#' ## load the sample identifier mapping
#' mapping <- read.csv("mapping.csv")
#' ## create renamed dataset in directory 'data-renamed/' based on the dataset files in 'data/'.
#' rename.aries(mapping, "data", "data-renamed")
#' }
#'
#' @author Matthew Suderman (matthew.suderman@bristol.ac.uk)
rename.aries <- function(mapping, data.dir, output.dir) {
    stopifnot(is.data.frame(mapping) && all(c("aln","new") %in% colnames(mapping)))
    stopifnot(file.exists(data.dir))
    for (filename in file.path(data.dir,
                               c("samplesheet/data.Robj",
                                 "control_matrix/data.txt",
                                 "qc.objects_all/data.Robj",
                                 "qc.objects_clean/data.Robj",
                                 "detection_p_values/data.Robj",
                                 "betas/data.Robj")))
        stopifnot(file.exists(filename))

    mapping$aln <- as.character(mapping$aln)
    mapping$new <- as.character(mapping$new)

    my.write.table <- function(x, file) {
        write.table(x,
                    sep="\t",
                    col.names=T,
                    row.names=F,
                    quote=F,
                    file=file)
    }
    my.file.path <- function(...) {
        filename <- file.path(...)
        dir <- dirname(filename)
        if (!file.exists(dir))
            dir.create(dir, recursive=T)
        filename
    }
    my.load <- function(filename) {
        envir <- sys.frame(sys.nframe())
        object_names <- load(filename, envir=envir)
        L <- lapply(object_names, get, envir=envir)
        names(L) <- object_names
        L
    }

    if (!file.exists(output.dir))
        dir.create(output.dir, recursive=T)
    output.dir <- normalizePath(output.dir)

    cat(date(), "change to data directory\n")
    current.dir <- getwd()
    on.exit(setwd(current.dir))
    setwd(data.dir)

    cat(date(), "load samplesheet\n")
    filename <- "samplesheet/data.Robj"
    samplesheet <- my.load(filename)[[1]]
    stopifnot("ALN" %in% colnames(samplesheet))
    stopifnot("Slide" %in% colnames(samplesheet))
    stopifnot("BCD_plate" %in% colnames(samplesheet))
    stopifnot("Sample_Name" %in% colnames(samplesheet))

    samplesheet$ALN <- as.character(samplesheet$ALN)
    samplesheet$Slide <- as.character(samplesheet$Slide)
    samplesheet$BCD_plate <- as.character(samplesheet$BCD_plate)
    samplesheet$time_code <- as.character(samplesheet$time_code)
    samplesheet$time_point <- as.character(samplesheet$time_point)

    cat(date(), "restrict to samples in both the dataset and the mapping\n")
    samplesheet <- samplesheet[which(samplesheet$ALN %in% mapping$aln),,drop=F]
    stopifnot(nrow(samplesheet) > 0)

    cat(date(), "convert the samplesheet\n")
    samplesheet <- samplesheet[sample(1:nrow(samplesheet)),]
    samplesheet$ALN <- mapping$new[match(samplesheet$ALN, mapping$aln)]

    samplesheet <- samplesheet[,c("Sample_Name",
                                  "ALN", ##ALN
                                  "QLET",
                                  "Slide",
                                  "sentrix_row",
                                  "sentrix_col",
                                  "time_code",
                                  "time_point",
                                  "Sex",
                                  "BCD_plate",
                                  "sample_type",
                                  "additive",
                                  "age",
                                  "duplicate.rm",
                                  "genotypeQCkids",
                                  "genotypeQCmums")]

    colnames(samplesheet) <-
                               c("Sample_Name",
                                  "dnam_450_g0m_g1", ##ALN
                                  "QLET",
                                  "Slide",
                                  "sentrix_row",
                                  "sentrix_col",
                                  "time_code",
                                  "time_point",
                                  "Sex",
                                  "BCD_plate",
                                  "sample_type",
                                  "additive",
                                  "age",
                                  "duplicate.rm",
                                  "genotypeQCkids",
                                  "genotypeQCmums")

    save(samplesheet, file=my.file.path(output.dir, filename))

    cat(date(), "saving cell counts files\n")
    filenames <- list.files(path="derived/cellcounts",
                            recursive=T,
                            full.names=T)
    for (filename in filenames) {
        counts <- read.table(filename, sep="\t", header=T)
        counts <- counts[counts[[1]] %in% samplesheet$Sample_Name,]
        my.write.table(counts, file=my.file.path(output.dir, filename))
    }

    cat(date(), "saving control matrix\n")
    filename <- "control_matrix/data.txt"
    control.matrix <- read.table(filename, sep="\t", header=T)
    control.matrix <- control.matrix[match(samplesheet$Sample_Name, control.matrix[["Sample_Name"]]),]
    my.write.table(control.matrix, file=my.file.path(output.dir, filename))

    cat(date(), "saving qc.objects\n")
    for (filename in c("qc.objects_all/data.Robj", "qc.objects_clean/data.Robj")) {
        qc.objects <- my.load(filename)[[1]]
        qc.objects <- qc.objects[match(samplesheet$Sample_Name, names(qc.objects))]
        for (i in 1:length(qc.objects)) {
            qc.objects[[i]]$basename <- NULL
            qc.objects[[i]]$snp.betas <- NULL
            qc.objects[[i]]$samplesheet <- samplesheet[i,]
        }
        save(qc.objects, file=my.file.path(output.dir, filename))
        rm(qc.objects)
        gc()
    }

    cat(date(), "saving betas\n")
    filename <- "betas/data.Robj"
    betas <- my.load(filename)[[1]]
    betas <- betas[,match(samplesheet$Sample_Name, colnames(betas))]
    idx <- which(substring(rownames(betas),1,2) == "rs")
    if (length(idx) > 0)
        betas <- betas[-idx,,drop=F]
    save(betas, file=my.file.path(output.dir, filename))
    rm(betas)
    gc()

    cat(date(), "saving detection p-values\n")
    filename <- "detection_p_values/data.Robj"
    detection.p <- my.load(filename)[[1]]
    detection.p <- detection.p[,match(samplesheet$Sample_Name, colnames(detection.p))]
    idx <- which(substring(rownames(detection.p),1,2) == "rs")
    if (length(idx) > 0)
        detection.p <- detection.p[-idx,,drop=F]
    save(detection.p, file=my.file.path(output.dir, filename))

    invisible(samplesheet)
}
