unlink(".RData")
gc()
library(meffil)
options(mc.cores=8)
load("/projects/MRC-IEU/groups/ARIES/normalisation/qc.objects.clean160915.Robj")

for (i in 1:length(qc.objects)){
basename<-paste("/projects/MRC-IEU/groups/ARIES/ariesidats_all",qc.objects[[i]]$samplesheet$Sample_Name,sep="/")
qc.objects[[i]]$samplesheet$Basename<-basename
qc.objects[[i]]$basename<-basename
}

#cc<-mclapply(qc.objects,function (qc.object) meffil.estimate.cell.counts(qc.object,cell.type.reference="blood gse35069 complete")$counts)
#cc2<-data.frame(IID=names(cc),t(do.call(cbind,cc)))

#save(cc,file="cell.counts.Robj")
#write.table(cc2,"cellcounts.5469.txt",sep="\t",quote=F,row.names=F,col.names=T)

#meffil.list.cell.type.references()
#[1] "blood gse35069"             "blood gse35069 complete"   
#[3] "cord blood gse68456"        "gervin and lyle cord blood"

load("samplesheet_aries.160915.Robj")
samplesheet<-samplesheet[samplesheet$time_point=="cord",]

m<-match(samplesheet$Sample_Name,names(qc.objects))
qc.objects<-qc.objects[m]
###
cc<-mclapply(qc.objects,function (qc.object) meffil.estimate.cell.counts(qc.object,cell.type.reference="cord blood gse68456")$counts)
cc2<-data.frame(IID=names(cc),t(do.call(cbind,cc)))

save(cc,file="cell.counts.cord.gse68456.Robj")
write.table(cc2,"cellcounts.cord.gse68456.txt",sep="\t",quote=F,row.names=F,col.names=T)
####
cc<-mclapply(qc.objects,function (qc.object) meffil.estimate.cell.counts(qc.object,cell.type.reference="gervin and lyle cord blood")$counts)
cc2<-data.frame(IID=names(cc),t(do.call(cbind,cc)))

save(cc,file="cell.counts.cord.gervinandlyle.Robj")
write.table(cc2,"cellcounts.cord.gervinandlyle.txt",sep="\t",quote=F,row.names=F,col.names=T)





