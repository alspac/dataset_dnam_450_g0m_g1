compare.releases <- function(mapping, x.dir, y.dir) {
    stopifnot(is.data.frame(mapping) && all(c("aln","new") %in% colnames(mapping)))
    stopifnot(file.exists(x.dir))
    stopifnot(file.exists(y.dir))

    mapping$aln <- as.character(mapping$aln)
    mapping$new <- as.character(mapping$new)

    my.load <- function(filename) {
        envir <- sys.frame(sys.nframe())
        object_names <- load(filename, envir=envir)
        L <- lapply(object_names, get, envir=envir)
        names(L) <- object_names
        L
    }
        
    cat(date(), "comparing samplesheets ...\n")
    x.samplesheet <- my.load(file.path(x.dir, "samplesheet/data.Robj"))[[1]]
    y.samplesheet <- my.load(file.path(y.dir, "samplesheet/data.Robj"))[[1]]

    ## synchronize samplesheets and sample names
    y.samplesheet$ALN <- mapping$aln[match(y.samplesheet$ALN, mapping$new)]
    x.id <- with(x.samplesheet, paste(ALN, QLET, time_code, duplicate.rm))
    y.id <- with(y.samplesheet, paste(ALN, QLET, time_code, duplicate.rm))
    x.samplesheet <- x.samplesheet[match(y.id, x.id),]
    sample.names <- data.frame(original=x.samplesheet$Sample_Name,
                               new=y.samplesheet$Sample_Name)

    cat(date(), "mapping slide identifiers ...\n")
    slides <- data.frame(original=as.character(x.samplesheet$Slide),
                         new=as.character(y.samplesheet$Slide),
                         stringsAsFactors=F)
    slides <- unique(slides)
    if (any(table(slides$original) > 1) || any(table(slides$new) > 1)) {
        stop("not a 1-1 mapping between slide identifiers")
    }

    cat(date(), "mapping plate identifiers ...\n")
    plates <- data.frame(original=as.character(x.samplesheet$BCD_plate),
                         new=as.character(y.samplesheet$BCD_plate),
                         stringsAsFactors=F)
    plates <- unique(plates)
    if (any(table(plates$original) > 1) || any(table(plates$new) > 1)) {
        stop("not a 1-1 mapping between plate identifiers")
    }

    compare.data.frames <- function(x,y,column) {
        if (!(all(colnames(x) %in% colnames(y))
              && all(colnames(y) %in% colnames(x))))
            return(FALSE)
        rownames(x) <- x[,1]
        rownames(y) <- sample.names$original[match(y[,1], sample.names$new)]
        x <- x[match(rownames(y), rownames(x)),]
        x <- x[,-1]
        y <- y[,-1]
        x <- x[,match(colnames(y), colnames(x))]
        diff <- max(abs(as.matrix(x) - as.matrix(y)), na.rm=T)
        return(diff <= 2e-16)
    }
    
    compare.matrices <- function(x,y) {
        if (!(all(rownames(x) %in% rownames(y))
              && all(rownames(y) %in% rownames(x))))
            return(FALSE)
        colnames(y) <- sample.names$original[match(colnames(y), sample.names$new)]
        x <- x[match(rownames(y), rownames(x)),
               match(colnames(y), colnames(x)), drop=F]
        diff <- max(abs(x-y))
        return(diff <= 2e-16)
    }

    compare.qc.objects <- function(x,y) {
        names(y) <- sample.names$original[match(names(y), sample.names$new)]
        idx <- match(names(y), names(x))
        if (any(is.na(idx)))
            return(FALSE)
        x <- x[idx]
        x.controls <- sapply(x,function(x) x$controls)
        y.controls <- sapply(y,function(y) y$controls)
        diff <- max(abs(x.controls - y.controls))
        return(diff <= 2e-16)
    }
    

    cat(date(), "checking cell counts estimates ...\n")
    x.filenames <- list.files(path=file.path(x.dir, "derived/cellcounts"), recursive=T, full.names=T)
    names(x.filenames) <- basename(x.filenames)
    y.filenames <- list.files(path=file.path(y.dir, "derived/cellcounts"), recursive=T, full.names=T)
    names(y.filenames) <- basename(y.filenames)
    
    for (name in union(names(x.filenames), names(y.filenames))) {
        x.counts <- read.table(x.filenames[name], sep="\t", header=T)
        y.counts <- read.table(y.filenames[name], sep="\t", header=T)
        if (!compare.data.frames(x.counts, y.counts))
            stop(paste("Cell count estimates", name, "do not match"))
    }

    cat(date(), "checking control matrices ...\n")
    x.control.matrix <- read.table(file.path(x.dir, "control_matrix/data.txt"), sep="\t", header=T)
    y.control.matrix <- read.table(file.path(y.dir, "control_matrix/data.txt"), sep="\t", header=T)
    if (!compare.data.frames(x.control.matrix, y.control.matrix))
        stop(paste("control matrices", name, "do not match"))

    cat(date(), "check qc objects ...\n")
    for (filename in c("qc.objects_all/data.Robj", "qc.objects_clean/data.Robj")) {
        x.qc.objects <- my.load(file.path(x.dir, filename))[[1]]
        y.qc.objects <- my.load(file.path(y.dir, filename))[[1]]
        if (!compare.qc.objects(x.qc.objects, y.qc.objects))
            stop(paste("qc objects in", filename, "do not match"))
    }

    cat(date(), "checking beta matrices ...\n")
    x.betas <- my.load(file.path(x.dir, "betas/data.Robj"))[[1]]
    sites <- sample(rownames(x.betas), 10000, replace=F)
    x.betas <- x.betas[sites,]
    y.betas <- my.load(file.path(y.dir, "betas/data.Robj"))[[1]]
    sites <- intersect(sites, rownames(x.betas))
    x.betas <- x.betas[sites,]
    y.betas <- y.betas[sites,]
    if (!compare.matrices(x.betas, y.betas)) 
        stop("beta matrices do not match")
    
    cat(date(), "checking detection p-value matrices ...\n")
    x.detect <- my.load(file.path(x.dir, "detection_p_values/data.Robj"))[[1]]
    x.detect <- x.detect[sites,]
    y.detect <- my.load(file.path(y.dir, "detection_p_values/data.Robj"))[[1]]
    y.detect <- y.detect[sites,]
    if (!compare.matrices(x.detect, y.detect)) 
        stop("detection p-value matrices do not match")
}
