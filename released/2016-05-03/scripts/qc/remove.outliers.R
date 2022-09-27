library(meffil)

load("/projects/MRC-IEU/groups/ARIES/test/QC/ARIESall_qc.Robj")
length(qc.objects)

outliers<-read.table("~/ARIES/methylation/normalisation/outliers_ARIES160915.txt",header=T,sep="\t",stringsAsFactors=F)
outliers<-unique(outliers$sample.name)

outliers2<-read.table("~/ARIES/methylation/normalisation/pca_outliers_pc15.txt",header=T,sep="\t",stringsAsFactors=F)
outliers2<-unique(outliers2$Sample_Name)

outliers<-unique(c(outliers,outliers2))

qc.objects <- meffil.remove.samples(qc.objects,outliers)
save(qc.objects,file="/projects/MRC-IEU/groups/ARIES/test/QC/qc.objects.clean160915.Robj")
length(qc.objects)

qc.parameters <-meffil.qc.parameters(beadnum.samples.threshold=0.1,detectionp.samples.threshold=0.1,snp.concordance.threshold=0.95,detectionp.cpgs.threshold = 0.1, beadnum.cpgs.threshold = 0.1,sex.outlier.sd=5,sample.genotype.concordance.threshold=0.8)

load("/projects/MRC-IEU/groups/ARIES/test/QC/samplesheet_aries.Robj")
plinkfiles <- list.files("~/ARIES/methylation/normalisation/genotypes/",pattern="*.raw")
plinkfiles<-paste("~/ARIES/methylation/normalisation/genotypes/",plinkfiles,sep="")
genotypes <- meffil.extract.genotypes(plinkfiles)

genotypes <- genotypes[,match(samplesheet$id, colnames(genotypes))] #4904
colnames(genotypes) <- samplesheet$Sample_Name

qc.summary <- meffil.qc.summary(qc.objects,parameters=qc.parameters,genotypes=genotypes)
traceback()
meffil.qc.report(qc.summary, output.file="ARIESall_qc.clean-report160915.html")
traceback()
save(qc.summary,file="~/ARIES/methylation/normalisation/qcsummary.clean160915.Robj")
