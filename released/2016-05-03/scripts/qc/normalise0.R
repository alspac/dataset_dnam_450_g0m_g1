library(meffil)
options(mc.cores=16)

samplesheet <- meffil.create.samplesheet("/panfs/panasas01/sscm/epzjlm/ARIES/methylation/normalisation150508/idat")
m<-read.table("/panfs/panasas01/sscm/epzjlm/ARIES/methylation/manifest/LIMS_Current_ARIES.only.tab",sep="\t",header=T)
w<-which(names(m)%in%c("Basename"))
m<-m[,-w]
w<-which(names(m)%in%c("gender"))
names(m)[w]<-c("Sex")
#change mums to female
w1<-which(m$time_code%in%c("FOM","antenatal","TF1-3")) #1994
m[w1,w]<-rep("f",length(w1))

table(m$time_code,m$Sex)
           
#               f    m
#  antenatal 1100    0
#  cord       574  553
#  F17        402  394
#  F7         557  529
#  FOM        894    0
#  TF1-3      189    0
#  TF3        147  130

samplesheet<-merge(samplesheet[,-2],m[,-2],by.x="Sample_Name",by.y="sentrix",all.x=T)
samplesheet$Sex<-gsub("m","M",samplesheet$Sex)
samplesheet$Sex<-gsub("f","F",samplesheet$Sex)

samplesheet$id2<-paste(samplesheet$ALN,samplesheet$QLET2,samplesheet$time_code)

#qc.objects <- meffil.qc(samplesheet, cell.type.reference="blood gse35069", verbose=TRUE)
#traceback()
#save(qc.objects,file="ARIESall_qc.Robj")

#snp.betas<-meffil.snp.betas(qc.objects)
#traceback()
#save(snp.betas,file="ARIESall_snp.betas.Robj")

load("/panfs/panasas01/sscm/epzjlm/ARIES/methylation/normalisation150508/ARIESall_qc.Robj")
load("/panfs/panasas01/sscm/epzjlm/ARIES/methylation/normalisation150508/ARIESall_snp.betas.Robj")

plinkfiles <- list.files("/panfs/panasas01/sscm/epzjlm/ARIES/methylation/genotypes/",pattern="*.raw")
plinkfiles<-paste("/panfs/panasas01/sscm/epzjlm/ARIES/methylation/genotypes/",plinkfiles,sep="")
genotypes <- meffil.extract.genotypes(plinkfiles)
genotypes <- genotypes[,match(samplesheet$id, colnames(genotypes))] #4904
colnames(genotypes) <- samplesheet$Sample_Name

#qc.parameters <-meffil.qc.parameters(beadnum.samples.threshold=0.1,detectionp.samples.threshold=0.1,snp.concordance.threshold=0.95,detectionp.cpgs.threshold = 0.1, beadnum.cpgs.threshold = 0.1,sex.outlier.sd=5,sample.genotype.concordance.threshold=0.8)
#qc.summary <- meffil.qc.summary(qc.objects,parameters=qc.parameters,genotypes=genotypes)
#traceback()
#meffil.qc.report(qc.summary, output.file="/panfs/panasas01/sscm/epzjlm/ARIES/methylation/normalisation160915/ARIESall_qc-report.html")
#traceback()
#save(qc.summary,file="/panfs/panasas01/sscm/epzjlm/ARIES/methylation/normalisation160915/qcsummary.Robj")

load("/panfs/panasas01/sscm/epzjlm/ARIES/methylation/normalisation160915/qcsummary.Robj")

outlier<-qc.summary$bad.samples
length(unique(outlier$sample.name)) #242
table(outlier$issue)

#            Control probe (dye.bias) Control probe (extension.G.31698466) 
#                                   9                                    8 
#Control probe (extension.G.74666473)             Control probe (oob.G.1%) 
#                                   9                                    2 
#    Control probe (spec2.G.17661470)     Control probe (spec2.G.29662396) 
#                                   2                                    2 
#    Control probe (spec2.G.34730329)          Control probe (spec2.ratio) 
#                                   5                                   91 
#                   Detection p-value                    Genotype mismatch
#                                 147  #166                                386 #411
#                    Low bead numbers           Methylated vs Unmethylated 
#                                   2                                   68 
#                        Sex mismatch                    X-Y ratio outlier 
#                                 161                                   30 


w<-which(outlier$issue%in%c("Control probe (dye.bias)","Methylated vs Unmethylated","X-Y ratio outlier","Low bead numbers","Detection p-value","Sex mismatch","Genotype mismatch"))
outlier<-outlier[w,] #150 #299 #256 #803 #847

#  Control probe (dye.bias)          Detection p-value 
#                         9                        147 #166
#         Genotype mismatch           Low bead numbers 
#                       386 #411                         2
#Methylated vs Unmethylated               Sex mismatch 
#                        68                        161 
#         X-Y ratio outlier 
#                        30 


##bad cpgs
bad.cpgs<-qc.summary$bad.cpgs
table(bad.cpgs$issue)

#Detection p-value Detection p-value, Low bead number
#4176                                  6
#Low bead number
#64

#23 detP probes on chrY

#sex mismatches
sex.sum<-qc.summary$sex.summary$tab[which(qc.summary$sex.summary$tab$sex.mismatch=="TRUE"),]
dim(sex.sum)
#[1] 161
sex.outl<-data.frame(sample.name=sex.sum$sample.name,issue="Sex mismatch")

w<-which(outlier$issue%in%c("Sex mismatch"))
length(which(outlier[w,1]%in%sex.outl[,1]))
#[1] 161

#outlier<-rbind(outlier,sex.outl)
length(unique(outlier$sample.name))
#[1] 263 #450 #465

m1<-merge(m,sex.sum,by.x="sentrix",by.y="sample.name",all.y=T)
table(m1$time_code)

#antenatal      cord       F17        F7       FOM     TF1-3       TF3 
#       37        57        12        17        29         0         9 


#Dye bias
w<-which(outlier$issue%in%"Control probe (dye.bias)")
dyebias<-samplesheet[which(samplesheet$Sample_Name%in%outlier$sample.name[w]),"Sample_Name"] #9

dyebias
#[1] "8784225091_R02C02" "8784225091_R05C01" "8784225091_R05C02"
#[4] "8784225091_R06C02" "8784225103_R01C01" "8784225103_R01C02"
#[7] "8784225103_R02C02" "8784225103_R03C01" "8963282048_R06C02"


dyebias<-samplesheet[samplesheet$Sample_Name%in%dyebias,]#9 #7 on BCD070- 12 samples in total on BCD plate

samplesheet[which(samplesheet$Slide%in%unique(dyebias$Slide)[1:2]),"Sample_Name"] #13 8784225091 & 8784225103

dyebias<-data.frame(samplesheet[samplesheet$Slide%in%unique(dyebias$Slide)[1:2],c("Sample_Name")],issue="Control probe (dye.bias)")
names(dyebias)<-c("sample.name","issue") #13
outlier<-rbind(outlier,dyebias) #816

#genotype QC
genoqc<-read.table("/panfs/panasas01/sscm/epzjlm/ARIES/methylation/normalisation150508/failed_GWAS_QC160915.txt",sep="\t",header=T,stringsAsFactors=F) #52
outlier<-rbind(outlier,genoqc) #928 #972
length(unique(outlier$sample.name)) #560 #575
write.table(unique(outlier),"/panfs/panasas01/sscm/epzjlm/ARIES/methylation/normalisation160915/outliers_ARIES160915.txt",sep="\t",col.names=T,row.names=F,quote=F) #860 #964


#genotype concordance
genotype<-data.frame(qc.summary$genotype.summary$tabs$samples)
snp<-data.frame(qc.summary$genotype.summary$tabs$snps)

pdf("/panfs/panasas01/sscm/epzjlm/ARIES/methylation/normalisation160915/scattersampleconcordance.pdf",width=6,height=6)
plot(genotype$concordance[order(genotype$concordance)],cex=0.8,pch=16,ylab="sample concordance (%)")
abline(h=0.9,col="red")
dev.off()

##
#w<-which(outlier$issue%in%"Genotype mismatch")

genotypefail<-genotype[which(genotype$concordance<0.90),]#470
dim(genotypefail)
genotypefail<-genotype[which(genotype$concordance<0.80),]#411
dim(genotypefail)
genotypefail<-data.frame(sample.name=genotypefail$sample.name,issue="GWA discordance 0.8") #411

#outlier<-rbind(outlier,genotypefail) #already added 
length(unique(outlier$sample.name)) #560 #575


#sample concordance-compare with old results

#m<-read.table("/panfs/panasas01/sscm/epzjlm/ARIES/methylation/meffil_results/sampleconcordance.txt",sep="\t",header=T) #4782
#SNPok<-read.table("/panfs/panasas01/sscm/epzjlm/ARIES/methylation/meffil_results/SNPswith99concordance.txt",sep="\t",header=T)
#w1<-match(colnames(snp.betas),m$sentrix)
#rm.samples<-colnames(snp.betas)[which(is.na(w1))]
#qc.objects.old <- meffil.remove.samples(qc.objects,rm.samples )
#qc.parameters <-meffil.qc.parameters(beadnum.samples.threshold=0.1,detectionp.samples.threshold=0.1,snp.concordance.threshold=0,detectionp.cpgs.threshold 0.1, beadnum.cpgs.threshold = 0.1)
#w1<-which(row.names(genotypes)%in%SNPok[,1])
#qc.summary.old <- meffil.qc.summary(qc.objects.old,parameters=qc.parameters,genotypes=genotypes[w1,])
#qc.summary.old2<- meffil.qc.summary(qc.objects.old,parameters=qc.parameters,genotypes=genotypes)

#genotype.old<-data.frame(qc.summary.old$genotype.summary$tabs$samples)
#w<-match(m$sentrix,genotype.old$sample.name)
#genotype.old<-genotype.old[w,]

#genotype.old2<-data.frame(qc.summary.old2$genotype.summary$tabs$samples)
#w<-match(m$sentrix,genotype.old2$sample.name)
#genotype.old2<-genotype.old2[w,]

#genotype.old<-read.table("/panfs/panasas01/sscm/epzjlm/ARIES/methylation/normalisation150508/SNP.concordance5469samples.txt",sep="\t",header=F)
#genotype.snp<-data.frame(qc.summary$genotype.summary$tabs$snps)
#o<-order(genotype.snp$snp.name)
#genotype.snp<-genotype.snp[o,]

#pdf("concordancecheck.pdf")
#plot(genotype.old[,6],genotype.snp$concordance,xlab="concordance Josine",ylab="concordance Matt",xlim=c(0.75,1),ylim=c(0.75,1))
#abline(a=0,b=1)
#abline(h=0.8)
#abline(v=0.8)
#dev.off()

#correlation across timepoints

a<-apply(snp.betas,1,function(x) y<-kmeans(x, centers=c(0.2, 0.5, 0.8))$cluster)
a<-t(a)
SNPok<-as.character(qc.summary$genotype.summary$tabs$snps[which(qc.summary$genotype.summary$tabs$snps$concordance>0.95),"snp.name"]) #27
a2<-a[which(row.names(a)%in%SNPok),]
#[1]   27 5469 #62 5469 #taking the best 50%

m<-match(colnames(a2),samplesheet$Sample_Name)
ind.cond2<-samplesheet[m,]

snp.ind<-as.character(ind.cond2[,"ALN"]) #indiv id indep of timepoint-only ALN
snp.ind2<-as.character(ind.cond2[,"id"]) #ALN QLET
snp.ind3<-as.character(ind.cond2[,"id2"]) #ALN QLET timepoint

x<-a2
y<-a2

dim(x)
#[1]   27 5469 #62 5469
dim(y)
#[1]   27 5469 #62 5469

flp <- function(x ,y){
results<-matrix(NA,ncol(x),ncol(y))
for (i in 1:ncol(x)){
	cat(i,"\n")
  for (j in 1:ncol(y)){
    r1<-x[,i]
    r2<-y[,j]
    indd<-r1==r2
    concordance<-round(100*sum(indd)/nrow(x),2)
    results[i,j]<-concordance  
  }
}
return(results)
}

#conc.matrix<-flp(x,y)
#colnames(conc.matrix)<-colnames(x)
#row.names(conc.matrix)<-snp.ind3
#save(conc.matrix,file="/panfs/panasas01/sscm/epzjlm/ARIES/methylation/normalisation150508/conc.matrix.Robj")

load("/panfs/panasas01/sscm/epzjlm/ARIES/methylation/normalisation150508/conc.matrix.Robj")
o<-order(snp.ind3)
conc.matrix<-conc.matrix[o,o]
x<-x[,o]
y<-y[,o]
snp.ind<-snp.ind[o] #ALN
snp.ind2<-snp.ind2[o] #ALN_QLET
snp.ind3<-snp.ind3[o] #ALN_QLET_TIMECODE
ind.cond2<-ind.cond2[o,]

t1<-data.frame(table(snp.ind))
t2<-data.frame(table(snp.ind2))
t3<-data.frame(table(snp.ind3))

#table(t3$Freq)

#1    2    3
#4723  361    8

#set non-ALN matches to NA
conc.matrix1<-conc.matrix
conc.matrix1b<-conc.matrix
for (aln in 1:dim(t1)[1]){
cat(aln,"\n")
w1<-which(snp.ind%in%t1$snp.ind[aln]==T)	#all ALNs
w2<-which(snp.ind%in%t1$snp.ind[aln]==F)	#all non-ALNs
conc.matrix1[w1,w2]<-NA

conc.matrix1b[w1,w1]<-NA
#w1.out<-append(w1.out,w1b)
#w2.out<-append(w2.out,w2b)

}


#conc.matrix1[upper.tri(conc.matrix1)] <- NA
diag(conc.matrix1)<-NA

conc.matrix2<-conc.matrix1
conc.matrix2b<-conc.matrix1
#conc.matrix2[upper.tri(conc.matrix2)] <- NA
row.names(conc.matrix)<-snp.ind3
#set non-ALN_QLET matches to NA

for (aln_qlet in 1:dim(t2)[1]){
cat(aln_qlet,"\n")
w1<-which(snp.ind2%in%t2$snp.ind2[aln_qlet]==T)	
w2<-which(snp.ind2%in%t2$snp.ind2[aln_qlet]==F)	
conc.matrix2[w1,w2]<-NA

}


pdf("/panfs/panasas01/sscm/epzjlm/ARIES/methylation/normalisation160915/concordancebetweentimepoints.pdf",height=6,width=6)
hist(na.omit(as.numeric(conc.matrix1)),xlab="concordance per ALN across timepoints", breaks=seq(20,100,2),main=NULL)
hist(na.omit(as.numeric(conc.matrix2)),xlab="concordance per ALN_QLET across timepoints", breaks=seq(20,100,2),main=NULL)
hist(na.omit(as.numeric(conc.matrix1b)),xlab="concordance without ALN",main=NULL)
dev.off()

pdf("/panfs/panasas01/sscm/epzjlm/ARIES/methylation/normalisation160915/concordancebetweentimepointsscatter.pdf",height=6,width=6)
plot(sort(na.omit(as.numeric(conc.matrix1))),ylab="concordance per ALN across timepoints", main=NULL,pch=16)
abline(h=90,col="red")
abline(h=80,col="red")
plot(sort(na.omit(as.numeric(conc.matrix2))),ylab="concordance per ALN_QLET across timepoints", main=NULL,pch=16)
abline(h=90,col="red")
abline(h=80,col="red")
dev.off()

#CONCORDANCE WITHOUT ALN pairs shows that 80% is the maximum concordance between random pairs.

a<-data.frame(ID=snp.ind3,ALN_QLET=snp.ind2,ALN=snp.ind,sentrix=colnames(conc.matrix2))
f<-which(a$sentrix%in%row.names(genotype))
a$genotyped<-NULL
a$genotyped[f]<-"genotyped"

#check if mums are different between mums and kids: 31145 M antenatal
a$failurerow_mum<-NULL
a$successrow_mum<-NULL

for (id in 1:dim(t1)[1]){
    cat(id,"\n")
    
    #id=53
    if(t1$Freq[id]>1){
    #is sample duplicate?
    w<-which(a[,c("ALN")]%in%t1[,1][id]) #[1] 280 281 282 283 284 285 ALN
    
    
    g1<-grep("M",snp.ind2[w])
    if(length(g1)>1){
    s<-w[g1] #283 284 285 mums
    s1<-w[-g1] #280 281 282 kids
    for (i in 1:length(g1)){
    
    s1<-s1[s1%in%s[-i]==F]	#kids only
    
    conc.tmp<-conc.matrix1[s[i],s1]
    mm<-which(row.names(conc.tmp)%in%genotypefail$sample.name)
				if(length(mm)>0){
				conc.tmp[mm,]<-NA
				conc.tmp[,mm]<-NA}

    fail.tmp<-length(which(conc.tmp<90&!is.na(conc.tmp)))
    success.tmp<-length(which(conc.tmp>90&!is.na(conc.tmp)))
                
                a[s[i],"successrow_mum"]<-fail.tmp
                a[s[i],"failurerow_mum"]<-success.tmp
    }
    }
    }
    }
    
fail.mum<-a[which(a$failurerow_mum>1&a$successrow_mum==0),]    #USE THIS 9- all checked manually
r<-which((fail.mum$sentrix%in%genotypefail$sample.name==F)&fail.mum$genotyped=="genotyped")

#check not genotyped samples
w<-which(snp.ind%in%c("43013"))
conc.matrix1[w,w]
a[w,]
test<-a[which(a$ALN%in%c("43013")),"sentrix"]
genotype[genotype$sample.name%in%test,]

w<-which(snp.ind%in%c("31145"))
conc.matrix1[w,w]
a[w,]
test<-a[which(a$ALN%in%c("31145")),"sentrix"]
genotype[genotype$sample.name%in%test,]




####check if kids are different between mums and kids

a$failurerow_kid<-NULL
a$successrow_kid<-NULL

for (id in 1:dim(t1)[1]){
    cat(id,"\n")
    
    #id=53
    if(t1$Freq[id]>1){
    #is sample duplicate?
    w<-which(a[,c("ALN")]%in%t1[,1][id]) #[1] 280 281 282 283 284 285 ALN
    
    
    g1a<-grep("A",snp.ind2[w])
    g1b<-grep("B",snp.ind2[w])
    g1<-c(g1a,g1b) #1 2 3
    
    if(length(g1)>1){
    s<-w[g1] #280 281 282 kids
    s1<-w[-g1] #283 284 285
    for (i in 1:length(g1)){
    
    s1<-s1[s1%in%s[-i]==F]	#kids only
    
    conc.tmp<-conc.matrix1[s[i],s1]
    mm<-which(row.names(conc.tmp)%in%genotypefail$sample.name)
				if(length(mm)>0){
				conc.tmp[mm,]<-NA
				conc.tmp[,mm]<-NA}

    fail.tmp<-length(which(conc.tmp<90&!is.na(conc.tmp)))
    success.tmp<-length(which(conc.tmp>90&!is.na(conc.tmp)))
                
                a[s[i],"successrow_kid"]<-fail.tmp
                a[s[i],"failurerow_kid"]<-success.tmp
    }
    }
    }
    }
    
fail.kid<-a[which(a$failurerow_kid>1&a$successrow_kid==0),]    #USE THIS 13
r<-which((fail.kid$sentrix%in%genotypefail$sample.name==F)&fail.kid$genotyped=="genotyped")
   
#check duplicates first before comparing across timepoints 

a$failurerow_dupl<-NULL
a$successrow_dupl<-NULL
for (id in 1:dim(t3)[1]){
    cat(id,"\n")
    
    #id=1258; 369 duplicates;id=1167 multiple dupl per ALN
    if(t3$Freq[id]>1){
    #is sample duplicate?
    w<-which(a[,c("ID")]%in%t3[,1][id]) #1251 1252
    w4<-NULL
    dupl<-NULL
    
    #id=1258; 369 duplicates;id=1167 multiple dupl per ALN; id658
   #are there duplicates in ALN
    aln_qlet<-unique(a[w,"ALN_QLET"]) #34930M
    aln<-unique(a[w,"ALN"]) #34930
    sampleid<-a[which(a$ALN%in%aln),"ID"] 
    # sampleid
#[1] 34930 A 15up      34930 A cord      34930 A F7        34930 M antenatal
#[5] 34930 M antenatal 34930 M FOM       34930 M FOM      

#identify all duplicates in ALN
dupl<-t3[t3[,1]%in%sampleid,]
    dupl<-dupl[dupl$Freq>1,1] #[1] 34930 M antenatal 34930 M FOM 
    dupl<-which(a$ID%in%dupl) # [1] 1249 1250 1251 1252
 
if(length(w)>1){
        
   
 w2<-which(a[,c("ALN")]%in%unique(a[w,"ALN"])) #[1] 1246 1247 1248 1249 1250 1251 1252
    
           
       for (i in 1:length(w)){
                 x2<-x
                 y2<-y
        
        #remove failed samples
         mm1<-match(fail.kid$sentrix,colnames(x))
         x2[,mm1]<-NA
         y2[,mm1]<-NA
		 mm2<-match(fail.mum$sentrix,colnames(x))
         x2[,mm2]<-NA
         y2[,mm2]<-NA

        
        
       #remove other duplicates than sample of interest
       if(length(dupl)>1){
       w4<-dupl[dupl%in%w==F] #[1] 1249 1250
       mm<-match(a[w4,"sentrix"],colnames(x))
       if(length(mm)>0){
       x2[,mm]<-NA
       y2[,mm]<-NA
       }
       }



#remove duplicate
                w3<-w[-i] #1252
                x2[,w3]<-NA
                y2[,w3]<-NA
               
               #set non ALN-QLET nmatches to NA
                conc.tmp<-flp(x2[,w2],y2[,w2])
                w_t<-which(a[w2,"ALN_QLET"]%in%aln_qlet==T)
                w_f<-which(a[w2,"ALN_QLET"]%in%aln_qlet==F)
                conc.tmp[w_t,w_f]<-NA
                conc.tmp[w_f,w_t]<-NA

                
                row.names(conc.tmp)<-a$sentrix[w2]
                colnames(conc.tmp)<-a$sentrix[w2]
                diag(conc.tmp)<-NA
                
                
               				
			    #mm<-which(row.names(conc.tmp)%in%genotypefail$sample.name)
				
				#conc.tmp[mm,]<-NA
				#conc.tmp[,mm]<-NA
                
                
                fail.tmp<-apply(conc.tmp,1,function(x) length(x[x<80&!is.na(x)]))
                success.tmp<-apply(conc.tmp,1,function(x) length(x[x>80&!is.na(x)]))
                m2<-w[i] #1251
                m<-match(a[m2,"sentrix"],names(fail.tmp))
                a[m2,"failurerow_dupl"]<-fail.tmp[m]
                a[m2,"successrow_dupl"]<-success.tmp[m]
            
            
             
            
            }
                        
            
            }
            
            }
        }


#first remove wrong duplicates and then calculate number of successes and failures 

fail.out2<-a[which(a$successrow_dupl==0&a$failurerow_dupl>1),] #213 ####USE THIS# 184 doesn't filter out mum replicates #189
length(intersect(fail.out2$sentrix,genotypefail$sample.name)) #318/342,318/386 #202 #175 #178

#all discordances are also found with GWA discordance if genotyped 
r<-which((fail.out2$sentrix%in%genotypefail$sample.name==F)&fail.out2$genotyped=="genotyped")
#0
r2<-which((fail.out2$sentrix%in%genotypefail$sample.name==F)&is.na(fail.out2$genotyped))
fail.out2[r2,]

fail.out2[r2,1:6] #all checked
#               ID ALN_QLET   ALN           sentrix genotyped successrow_mum
#103  30273 A 15up   30273A 30273 8221932070_R01C01      <NA>             NA
#293    31151 A F7   31151A 31151 8963291073_R02C02      <NA>             NA
#314  31311 A 15up   31311A 31311 8942342075_R03C01      <NA>             NA
#1957 37688 A cord   37688A 37688 9007217040_R04C01      <NA>             NA
#2276   38960 A F7   38960A 38960 8221932111_R05C02      <NA>             NA
#2573 40159 A 15up   40159A 40159 8221932106_R01C01      <NA>             NA
#2610 40363 A cord   40363A 40363 8221932111_R06C01      <NA>             NA #
#2612   40363 A F7   40363A 40363 8221932070_R04C02      <NA>             NA #
#3333 43117 A 15up   43117A 43117 8963311122_R05C01      <NA>             NA
#4028 45996 A cord   45996A 45996 9007225133_R06C01      <NA>             NA
#4676 50752 A cord   50752A 50752 8221932139_R06C02      <NA>             NA

#ID ALN_QLET   ALN           sentrix genotyped successrow_mum
#105   30273 A TF3   30273A 30273 8221932070_R01C01      <NA>             NA
#292    31151 A F7   31151A 31151 8963291073_R02C02      <NA>             NA
#315   31311 A F17   31311A 31311 8942342075_R03C01      <NA>             NA
#342  31399 A cord   31399A 31399 8784225040_R05C02      <NA>             NA #
#1956 37688 A cord   37688A 37688 9007217040_R04C01      <NA>             NA
#2276   38960 A F7   38960A 38960 8221932111_R05C02      <NA>             NA
#2531 40029 A cord   40029A 40029 8942289024_R01C02      <NA>             NA #
#2575  40159 A TF3   40159A 40159 8221932106_R01C01      <NA>             NA
#3334  43117 A F17   43117A 43117 8963311122_R05C01      <NA>             NA
#4027 45996 A cord   45996A 45996 9007225133_R06C01      <NA>             NA
#4675 50752 A cord   50752A 50752 8221932139_R06C02      <NA>             NA

#remove all incorrect samples with duplicates before recalculating counts
w<-which(a$sentrix%in%fail.out2$sentrix)
w2<-which(a$sentrix%in%genotypefail$sample.name)
w3<-which(a$sentrix%in%fail.mum$sentrix)
w4<-which(a$sentrix%in%fail.kid$sentrix)

w<-unique(c(w,w2,w3,w4)) #397

conc.matrix3<-conc.matrix2
conc.matrix3[w,]<-NA
conc.matrix3[,w]<-NA

faillow80<-apply(conc.matrix3,1,function(x) length(x[x<80&!is.na(x)]))
failhigh80<-apply(conc.matrix3,1,function(x) length(x[x>80&!is.na(x)]))
faillow90<-apply(conc.matrix3,1,function(x) length(x[x<90&!is.na(x)]))
failhigh90<-apply(conc.matrix3,1,function(x) length(x[x>90&!is.na(x)]))
a$failurerow_withoutdupl<-NULL
a$successrow_withoutdupl<-NULL
a[,"failurerow_withoutdupl"]<-faillow80
a[,"successrow_withoutdupl"]<-failhigh80


#check FOM and antenatal samples
w<-which(a$successrow_withoutdupl==0&a$failurerow_withoutdupl==1) # FOM and antenatal not matching with each other and no duplicates
fail1<-a[w,] #24 #USE THIS

#remove wrong mum duplicates
fail2a<-a[which(a$successrow_withoutdupl==0&a$failurerow_withoutdupl==0),] #441
length(intersect(fail2a$sentrix,a$sentrix[w]))#397
r<-which((fail2a$sentrix%in%genotypefail$sample.name==F)&fail2a$genotyped=="genotyped")
#0
r2<-which((fail2a$sentrix%in%genotypefail$sample.name==F)&is.na(fail2a$genotyped))

fail2<-a[which(a$successrow_dupl==0&a$failurerow_dupl==1&a$failurerow_withoutdupl>1),] #11 USE THIS
fail2
#                    ID ALN_QLET   ALN           sentrix genotyped
#1429 35670 M antenatal   35670M 35670 8221932139_R05C01      <NA>
#1574 36140 M antenatal   36140M 36140 7766148051_R03C02      <NA>
#1681 36574 M antenatal   36574M 36574 8221932127_R04C01      <NA>
#1986       37809 M FOM   37809M 37809 8221932106_R02C02      <NA>
#2125 38335 M antenatal   38335M 38335 8221932139_R04C01      <NA>
#2466 39770 M antenatal   39770M 39770 8221932106_R05C02      <NA>
#2610      40363 A cord   40363A 40363 8221932111_R06C01      <NA>
#2612        40363 A F7   40363A 40363 8221932070_R04C02      <NA>
#2761 40938 M antenatal   40938M 40938 8221932070_R04C01      <NA>
#4299 46885 M antenatal   46885M 46885 8942342143_R02C02      <NA>
#5441       54181 M FOM   54181M 54181 8221932127_R01C01      <NA>

w<-which(a$successrow_withoutdupl==0&a$failurerow_withoutdupl>1) # FOM and antenatal not matching with each other and no duplicates
fail6<-a[w,] #17 #overlap
fail6<-fail6[which(fail6$sentrix%in%fail2$sentrix==F),] #6

r<-which((fail6$sentrix%in%genotypefail$sample.name==F)&fail6$genotyped=="genotyped")
#0
r2<-which((fail6$sentrix%in%genotypefail$sample.name==F)&is.na(fail6$genotyped))
fail6[r2,]

#                    ID ALN_QLET   ALN           sentrix genotyped
#704  32632 M antenatal   32632M 32632 8221932092_R03C02      <NA>
#1032 34064 M antenatal   34064M 34064 8942300185_R02C01      <NA>
#1761        36882 A F7   36882A 36882 7766148062_R04C02      <NA>
#1765      36896 A cord   36896A 36896 7766148117_R02C01      <NA>
#2575      40159 A cord   40159A 40159 8655677014_R03C02      <NA>
#4298        46885 A F7   46885A 46885 8942300090_R06C02      <NA>

g<-grep("M",fail6$ALN_QLET)
fail6<-unique(rbind(a[which(a$ALN_QLET%in%fail6$ALN_QLET[g]),],fail6))

w<-which(snp.ind%in%c("30094"))
conc.matrix3[w,w]
a[w,]
test<-a[which(a$ALN%in%c("30094")),"sentrix"]
genotype[genotype$sample.name%in%test,]



#select duplicates with overall failure 
fail3<-a[which(a$successrow_dupl==0&a$failurerow_dupl==1&a$successrow_withoutdupl==0&a$failurerow_withoutdupl==0),] #122 all these samples are already removed #127
r<-which((fail3$sentrix%in%genotypefail$sample.name==F)&fail3$genotyped=="genotyped")
#0
r2<-which((fail3$sentrix%in%genotypefail$sample.name==F)&is.na(fail3$genotyped))
r2<-which((fail3$sentrix%in%a$sentrix[w]==F))
fail3[r2,]
#0
intersect(fail3$sentrix,a$sentrix[w]) #122

#mums with replicates but only concordance between replicates but not with other sample. 
fail4<-a[which(a$successrow_dupl==0&a$failurerow_dupl==1&a$successrow_withoutdupl==1&a$failurerow_withoutdupl>0),] #6 all these samples are already 
fail4<-a[which(a$ALN_QLET%in%fail4$ALN_QLET),]
r<-which((fail4$sentrix%in%genotypefail$sample.name==F)&fail4$genotyped=="genotyped")
#
w<-which(a$sentrix%in%fail.out2$sentrix)
w3<-which(a$sentrix%in%fail.mum$sentrix)
w4<-which(a$sentrix%in%fail.kid$sentrix)
w5<-which(a$sentrix%in%fail1$sentrix)
w6<-which(a$sentrix%in%fail2$sentrix)
w7<-which(a$sentrix%in%fail4$sentrix)
w8<-which(a$sentrix%in%fail6$sentrix)

w<-unique(c(w,w2,w3,w5,w6,w7,w8)) #442 #469

conc.matrix3<-conc.matrix2
conc.matrix3[w,]<-NA
conc.matrix3[,w]<-NA

faillow80<-apply(conc.matrix3,1,function(x) length(x[x<80&!is.na(x)]))
failhigh80<-apply(conc.matrix3,1,function(x) length(x[x>80&!is.na(x)]))
faillow90<-apply(conc.matrix3,1,function(x) length(x[x<90&!is.na(x)]))
failhigh90<-apply(conc.matrix3,1,function(x) length(x[x>90&!is.na(x)]))
a$failurerow_withoutdupl<-NULL
a$successrow_withoutdupl<-NULL
a[,"failurerow_withoutdupl2"]<-faillow80
a[,"successrow_withoutdupl2"]<-failhigh80

fail5<-a[which(a$successrow_withoutdupl2==0&a$failurerow_withoutdupl2>1),] #0
fail5<-a[which(a$successrow_withoutdupl2==0&a$failurerow_withoutdupl2==1),] #0 


fail5<-a[which(a$successrow_withoutdupl2==1&a$failurerow_withoutdupl2>0),]
g<-grep("M",fail5$ALN_QLET)
fail5<-unique(rbind(a[which(a$ALN_QLET%in%fail5$ALN_QLET[g]),])) #4

w<-which(snp.ind%in%c("30094"))
conc.matrix3[w,w]
a[w,]
test<-a[which(a$ALN%in%c("30094")),"sentrix"]
genotype[genotype$sample.name%in%test,]


##
outlier3<-data.frame(sample.name=fail.kid[,4],issue="genotype concordance (>90%) between mum-kid")
outlier4<-data.frame(sample.name=fail.mum[,4],issue="genotype concordance (>90%) between mum-kid")
outlier5<-data.frame(sample.name=fail.out2[,4],issue="genotype discordance (80%) across duplicates")#1
outlier6<-data.frame(sample.name=fail1[,4],issue="genotype discordance (80%) across timepoints with mums")#22
outlier7<-data.frame(sample.name=fail2[,4],issue="genotype discordance (80%) across duplicates")
outlier8<-data.frame(sample.name=fail4[,4],issue="genotype discordance (80%) across mums duplicates")#4
outlier9<-data.frame(sample.name=fail5[,4],issue="genotype discordance (80%) across mums duplicates")#4
outlier10<-data.frame(sample.name=fail6[,4],issue="genotype discordance (80%) across timepoints")

outlier2<-rbind(outlier,outlier3,outlier4,outlier5,outlier6,outlier7,outlier8,outlier9,outlier10) #534
outliertimepoint<-rbind(outlier3,outlier4,outlier5,outlier6,outlier7,outlier8,outlier9,outlier10)

check1<-genotypefail[which((genotypefail$sample.name%in%outliertimepoint$sample.name==F)),"sample.name"]
check1<-a[a$sentrix%in%check1,]

outliers<-merge(outlier2,a,all.x=T,by.x="sample.name",by.y="sentrix")
o<-order(outliers$ALN,outliers$sample.name)
outliers<-outliers[o,]
write.table(unique(outliers),"/panfs/panasas01/sscm/epzjlm/ARIES/methylation/normalisation160915/outliers_ARIES160915.txt",sep="\t",col.names=T,row.names=F,quote=F)

w<-which(a$sentrix%in%unique(outliers$sample.name))
a.clean<-a[-w,]
length(unique(a.clean$ID)) #4902 #4794
length(a.clean$ID) #4983 #4866

#[1] BCD070 BCD069 BCD070 BCD070 BCD070 BCD070 BCD070 BCD070 BCD119

