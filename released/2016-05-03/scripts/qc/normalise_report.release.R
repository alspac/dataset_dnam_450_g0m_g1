library(meffil)
pc=10
options(mc.cores=16)

load(paste("/projects/MRC-IEU/groups/ARIES/normalisation/aries_funnorm.randomeffect.pc",pc,".160915.Robj",sep=""))
gc()
load(paste("/projects/MRC-IEU/groups/ARIES/normalisation/ARIESall_norm.obj.randomeffect.pc",pc,".160915.Robj",sep=""))
gc()
norm.beta<-norm.beta.random
norm.objects<-norm.objects.random

load("./release/samplesheet.release3.4854indiv.inclreplicates_popstrat.Robj")
m<-match(colnames(norm.beta),samplesheet$Sample_Name)
samplesheet<-samplesheet[m,]

norm.parameters<-meffil.normalization.parameters(norm.objects,variables=c("sample_type2","time_code","Slide","Sex"),control.pcs=1:10,batch.pcs=1:10)

pcs <- meffil.methylation.pcs(norm.beta,probe.range=20000)
save(pcs,file="pcs.aries_funnorm.randomeffect.pc10.Robj")
load("pcs.aries_funnorm.randomeffect.pc10.Robj")
norm.summary <- meffil.normalization.summary(norm.objects,pcs=pcs,parameters=norm.parameters)
gc()

meffil.normalization.report(norm.summary, output.file=paste("/data/epzjlm/ARIES/methylation/normalisation/release/normalization-report.pc",pc,"160915.html"))

tp<-names(table(samplesheet$time_point))

for (i in 1:length(tp)){
cat(tp[i],"\n")
w<-which(samplesheet$time_point%in%tp[i])
length(w)
norm.objects.tp<-norm.objects[w]
length(norm.objects.tp)
norm.beta.tp<-norm.beta[,w]
dim(norm.beta.tp)
pcs <- meffil.methylation.pcs(norm.beta.tp,probe.range=20000)
save(pcs,file=paste("pcs.",tp[i],".Robj",sep=""))
norm.summary.tp <- meffil.normalization.summary(norm.objects.tp,pcs=pcs,parameters=norm.parameters)
meffil.normalization.report(norm.summary.tp, output.file=paste("/data/epzjlm/ARIES/methylation/normalisation/release/normalization-report.pc",pc,".",tp[i],".160915.html",sep=""))
}




