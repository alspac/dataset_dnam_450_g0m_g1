library(meffil)
pc=10
options(mc.cores=16)

load("/projects/MRC-IEU/groups/ARIES/normalisation/qc.objects.clean160915.Robj")
gc()

for (i in 1:length(qc.objects)){
basename<-paste("/projects/MRC-IEU/groups/ARIES/ariesidats_all",qc.objects[[i]]$samplesheet$Sample_Name,sep="/")
qc.objects[[i]]$samplesheet$Basename<-basename
qc.objects[[i]]$basename<-basename
}

detp<-meffil.load.detection.pvalues(qc.objects)
save(detp,file="detp.Robj")

#load("./release/samplesheet.release3.4854indiv.inclreplicates_popstrat.Robj")



