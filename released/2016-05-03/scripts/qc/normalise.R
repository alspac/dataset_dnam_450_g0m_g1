args = (commandArgs(TRUE));
pc= as.numeric(args[1]);
cat(pc,"\n")

library(meffil)
#library(devtools)
#devtools::install_github('kaneplusplus/bigmemory')

#packageVersion("bigmemory")

options(mc.cores=16)
samplesheet <- meffil.create.samplesheet("/projects/MRC-IEU/groups/ARIES/ariesidats_all")
#samplesheet <- meffil.create.samplesheet("/projects/MRC-IEU/groups/ARIES/test/samplesheet_aries.Robj")

#meffil.set.current.cell.type.reference("blood gse35069")

#load("/projects/MRC-IEU/groups/ARIES/test/qc.objects.clean.Robj")
load("/projects/MRC-IEU/groups/ARIES/test/QC/qc.objects.clean160915.Robj")
length(qc.objects)

for (i in 1:length(qc.objects)){
basename<-paste("/projects/MRC-IEU/groups/ARIES/ariesidats_all",qc.objects[[i]]$samplesheet$Sample_Name,sep="/")
qc.objects[[i]]$samplesheet$Basename<-basename
qc.objects[[i]]$basename<-basename
}
#[1] 4930

#load(file="/projects/MRC-IEU/groups/ARIES/test/qcsummary.clean.Robj")
load(file="~/ARIES/methylation/normalisation/qcsummary.clean160915.Robj")
norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=pc)
save(norm.objects,file=paste("/data/epzjlm/ARIES/methylation/normalisation/ARIESall_norm.obj.pc",pc,".160915.Robj",sep="")) 

#load(paste("/data/epzjlm/methylation/normalisation/ARIESall_norm.obj.pc",pc,".Robj",sep=""))
#length(norm.objects)
#[1] 4930
 

#norm.beta <- meffil.normalize.samples(norm.objects,tempdir=".")

norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove=qc.summary$bad.cpgs$name)
save(norm.beta,file=paste("/data/epzjlm/ARIES/methylation/normalisation/aries_funnorm.norandomeffect.pc",pc,".160915.Robj",sep=""))
norm.parameters<-meffil.normalization.parameters(norm.objects,variables=c("Slide","BCD_plate","time_point","time_code","sample_type2"),control.pcs=1:10,probe.pcs=1:10,probe.range=20000)

#norm.summary <- meffil.normalization.summary(norm.beta, norm.objects,parameters=norm.parameters)
#meffil.normalization.report(norm.summary, output.file=paste("/data/epzjlm/ARIES/methylation/normalisation/normalization-report160915.pc",pc,".html",sep=""))

rm(norm.beta,norm.summary)
gc()
#load(file=paste("/data/epzjlm/ARIES/methylation/normalisation/ARIESall_norm.obj.pc",pc,".Robj",sep="")) 
#load(file=paste("/data/epzjlm/ARIES/methylation/normalisation/ARIESall_funnorm.pc",pc,".Robj",sep=""))

norm.objects.random <- meffil.normalize.quantiles(qc.objects, random.effects="Slide", number.pcs=pc)
save(norm.objects.random,file=paste("/data/epzjlm/ARIES/methylation/normalisation/ARIESall_norm.obj.randomeffect.pc",pc,".160915.Robj",sep="")) 

norm.beta.random <- meffil.normalize.samples(norm.objects.random, cpglist.remove=qc.summary$bad.cpgs$name)
save(norm.beta.random,file=paste("/data/epzjlm/ARIES/methylation/normalisation/aries_funnorm.randomeffect.pc",pc,".160915.Robj",sep=""))

#norm.summary.random <- meffil.normalization.summary(norm.beta.random, norm.objects.random)
#meffil.normalization.report(norm.summary.random, output.file=paste("/data/epzjlm/ARIES/methylation/normalisation/normalization-report160915.random.pc",pc,".html",sep=""))


