compare.releases <- function(mapping, old.dir, new.dir) {
    stopifnot(is.data.frame(mapping) && all(c("aln","new") %in% colnames(mapping)))
    stopifnot(file.exists(old.dir))
    stopifnot(file.exists(new.dir))

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
    old.samplesheet <- my.load(file.path(old.dir, "samplesheet/data.Robj"))[[1]]
    new.samplesheet <- my.load(file.path(new.dir, "samplesheet/data.Robj"))[[1]]

    cat("There are",
        nrow(old.samplesheet), "samples in the old and",
        nrow(new.samplesheet), "in the new.\n")

    stopifnot(all(new.samplesheet$ALN %in% mapping$new))
    
    ## synchronize samplesheets and sample names
    new.samplesheet <- new.samplesheet[which(new.samplesheet$ALN %in% mapping$new),,drop=F]
    if (nrow(new.samplesheet) == 0)
        stop("mapping$new does not match ALN's in the second dataset")
    
    new.aln <- mapping$aln[match(new.samplesheet$ALN, mapping$new)]
    old.id <- with(old.samplesheet, paste(ALN, QLET, time_code, duplicate.rm))
    new.id <- with(new.samplesheet, paste(new.aln, QLET, time_code, duplicate.rm))

    if (!all(new.id %in% old.id))
        stop("new dataset contains data not in the old dataset")

    ids <- intersect(old.id, new.id)
    if (length(ids) == 0)
        stop("there is something wrong with the ALN mapping")

    old.samplesheet <- old.samplesheet[match(ids, old.id),,drop=F]
    new.samplesheet <- new.samplesheet[match(ids, new.id),,drop=F]

    sample.name.mapping <- data.frame(
        old=old.samplesheet$Sample_Name,
        new=new.samplesheet$Sample_Name)
    
    compare.data.frames <- function(old,new,mapping) {
        rownames(old) <- old[,1]
        rownames(new) <- mapping$old[match(new[,1],mapping$new)]
        if (!all(rownames(new) %in% rownames(old)))
            return(FALSE)
        old <- old[match(rownames(new), rownames(old)),]
        old <- old[,-1]
        new <- new[,-1]
        old <- old[,match(colnames(new), colnames(old))]
        diff <- max(abs(as.matrix(old) - as.matrix(new)), na.rm=T)
        return(diff <= 2e-16)
    }
    
    compare.matrices <- function(old,new,mapping) {
        if (!identical(sort(rownames(old)),sort(rownames(new))))
            return(FALSE)
        colnames(new) <- mapping$old[match(colnames(new),mapping$new)]
        if (!all(colnames(new) %in% colnames(old)))
            return(FALSE)        
        old <- old[match(rownames(new), rownames(old)),
                   match(colnames(new), colnames(old)), drop=F]
        diff <- max(abs(old-new))
        return(diff <= 2e-16)
    }

    compare.qc.objects <- function(old,new,mapping) {
        compare.matrices(
            sapply(old,function(old) old$controls),
            sapply(new,function(new) new$controls),
            mapping)
    }
    

    cat(date(), "checking cell counts estimates ...\n")
    old.filenames <- list.files(path=file.path(old.dir, "derived/cellcounts"), recursive=T, full.names=T)
    names(old.filenames) <- basename(dirname(old.filenames))
    new.filenames <- list.files(path=file.path(new.dir, "derived/cellcounts"), recursive=T, full.names=T)
    names(new.filenames) <- basename(dirname(new.filenames))
   
    for (name in union(names(old.filenames), names(new.filenames))) {
        old.counts <- read.table(old.filenames[name], sep="\t", header=T)
        new.counts <- read.table(new.filenames[name], sep="\t", header=T)
        if (!compare.data.frames(old.counts, new.counts, sample.name.mapping))
            stop(paste("Cell count estimates", name, "do not match"))		
    }

    cat(date(), "checking control matrices ...\n")
    old.control.matrix <- read.table(file.path(old.dir, "control_matrix/data.txt"), sep="\t", header=T)
    new.control.matrix <- read.table(file.path(new.dir, "control_matrix/data.txt"), sep="\t", header=T)
    if (!compare.data.frames(old.control.matrix, new.control.matrix, sample.name.mapping))
        stop(paste("control matrices", name, "do not match"))

    cat(date(), "check qc objects ...\n")
    for (filename in c("qc.objects_all/data.Robj", "qc.objects_clean/data.Robj")) {
        old.qc.objects <- my.load(file.path(old.dir, filename))[[1]]
        new.qc.objects <- my.load(file.path(new.dir, filename))[[1]]
        if (!compare.qc.objects(old.qc.objects, new.qc.objects, sample.name.mapping))
            stop(paste("qc objects in", filename, "do not match"))
    }

    cat(date(), "checking beta matrices ...\n")
    old.betas <- my.load(file.path(old.dir, "betas/data.Robj"))[[1]]
    sites <- sample(rownames(old.betas), 10000, replace=F)
    old.betas <- old.betas[sites,]
    new.betas <- my.load(file.path(new.dir, "betas/data.Robj"))[[1]]
    sites <- intersect(sites, rownames(old.betas))
    old.betas <- old.betas[sites,]
    new.betas <- new.betas[sites,]
    if (!compare.matrices(old.betas, new.betas,sample.name.mapping)) 
        stop("beta matrices do not match")
    
    cat(date(), "checking detection p-value matrices ...\n")
    old.detect <- my.load(file.path(old.dir, "detection_p_values/data.Robj"))[[1]]
    old.detect <- old.detect[sites,]
    new.detect <- my.load(file.path(new.dir, "detection_p_values/data.Robj"))[[1]]
    new.detect <- new.detect[sites,]
    if (!compare.matrices(old.detect, new.detect,sample.name.mapping)) 
        stop("detection p-value matrices do not match")
}
