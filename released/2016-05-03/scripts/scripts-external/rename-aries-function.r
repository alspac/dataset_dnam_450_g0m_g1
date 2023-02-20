#' Remove identifying information from ARIES
#'
#' Make a copy of the ARIES dataset with identifying information
#' removed and sample names replaced with new identifiers.
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
#' such as SNPs and microarray identifiers removed.
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
    require(digest)

    mapping$aln <- as.character(mapping$aln)
    mapping$new <- as.character(mapping$new)

    cat(date(), "generating random seed based on mapping\n")
    seed <- digest(mapping, algo = "md5", serialize = TRUE)
    digits <- floor(log(.Machine$integer.max, 16))
    seed <- strtoi(substring(seed,nchar(seed)-digits+1, nchar(seed)), base=16)
    set.seed(seed)
       
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

    cat(date(), "restrict to samples in both the dataset and the mapping\n")
    samplesheet <- samplesheet[which(samplesheet$ALN %in% mapping$aln),,drop=F]
    stopifnot(nrow(samplesheet) > 0)

    cat(date(), "generating random slide identifiers\n")
    slides <- data.frame(original=unique(samplesheet$Slide), stringsAsFactors=F)
    slides$new <- sample(paste("SLIDE", 1:nrow(slides), sep=""))

    cat(date(), "generating random plate identifiers\n")
    plates <- data.frame(original=unique(samplesheet$BCD_plate), stringsAsFactors=F)
    plates$new <- sample(paste("PLATE", 1:nrow(plates), sep=""))

    cat(date(), "convert the samplesheet\n")
    samplesheet <- samplesheet[sample(1:nrow(samplesheet)),]
    samplesheet$ALN <- mapping$new[match(samplesheet$ALN, mapping$aln)]
    samplesheet$Slide <- slides$new[match(samplesheet$Slide, slides$original)]
    samplesheet$BCD_plate <- plates$new[match(samplesheet$BCD_plate, plates$original)]
    
    sample.names <- with(samplesheet, data.frame(original=Sample_Name,
                                                 new=paste(Slide,
                                                     "_R", sentrix_row,
                                                     "C", sentrix_col, sep=""),
                                                 stringsAsFactors=F))
    samplesheet$Sample_Name <- sample.names$new
    
    samplesheet <- samplesheet[,c("Sample_Name",
                                  "ALN",
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
    save(samplesheet, file=my.file.path(output.dir, filename))

    convert.data.frame <- function(x, column) {
        sample.names <- sample.names[which(sample.names$original %in% as.character(x[[column]])),]
        stopifnot(nrow(sample.names) > 0)
        idx <- match(sample.names$original, as.character(x[[column]]))
        x <- x[idx,,drop=F]
        x[[column]] <- sample.names$new
        x
    }
    convert.matrix <- function(x) {
        sample.names <- sample.names[which(sample.names$original %in% colnames(x)),]
        stopifnot(nrow(sample.names) > 0)
        idx <- match(sample.names$original, colnames(x))
        x <- x[,idx,drop=F]
        colnames(x) <- sample.names$new
        x
    }
    convert.list <- function(x) {
        sample.names <- sample.names[which(sample.names$original %in% names(x)),]
        idx <- match(sample.names$original, names(x))
        stopifnot(nrow(sample.names) > 0)
        x <- x[idx]
        names(x) <- sample.names$new
        x        
    }

    cat(date(), "convert cell counts files\n")
    filenames <- list.files(path="derived/cellcounts",
                            recursive=T,
                            full.names=T)
    for (filename in filenames) {
        counts <- read.table(filename, sep="\t", header=T)
        counts <- convert.data.frame(counts, 1)
        my.write.table(counts, file=my.file.path(output.dir, filename))
    }

    cat(date(), "convert control matrix\n")
    filename <- "control_matrix/data.txt"
    control.matrix <- read.table(filename, sep="\t", header=T)
    control.matrix <- convert.data.frame(control.matrix, "Sample_Name")
    my.write.table(control.matrix, file=my.file.path(output.dir, filename))

    cat(date(), "convert qc.objects\n")
    for (filename in c("qc.objects_all/data.Robj", "qc.objects_clean/data.Robj")) {
        qc.objects <- my.load(filename)[[1]]
        qc.objects <- convert.list(qc.objects)
        for (i in 1:length(qc.objects)) {
            qc.objects[[i]]$sample.name <- names(qc.objects)[i]
            qc.objects[[i]]$basename <- NULL
            qc.objects[[i]]$snp.betas <- NULL
            qc.objects[[i]]$samplesheet <- samplesheet[i,]
        }
        save(qc.objects, file=my.file.path(output.dir, filename))
        rm(qc.objects)
        gc()
    } 
    
    cat(date(), "convert betas\n")
    filename <- "betas/data.Robj"
    betas <- my.load(filename)[[1]]
    betas <- convert.matrix(betas)
    idx <- which(substring(rownames(betas),1,2) == "rs")
    if (length(idx) > 0)
        betas <- betas[-idx,,drop=F]
    save(betas, file=my.file.path(output.dir, filename))
    rm(betas)
    gc()

    cat(date(), "convert detection p-values\n")
    filename <- "detection_p_values/data.Robj"
    detection.p <- my.load(filename)[[1]]
    detection.p <- convert.matrix(detection.p)
    idx <- which(substring(rownames(detection.p),1,2) == "rs")
    if (length(idx) > 0)
        detection.p <- detection.p[-idx,,drop=F]
    save(detection.p, file=my.file.path(output.dir, filename))

    invisible(samplesheet)
}
