library('scran')
library('scater')
counts_aceto=read.csv('Aceto_countMatrix.csv')

counts_aceto=counts_aceto[!duplicated(counts_aceto$symbol),]
row.names(counts_aceto)=counts_aceto$symbol
counts_aceto=counts_aceto[,2:30]

cell_name=colnames(counts_aceto)[2:30]
cell_type=c('CL','SC','CL','CL','SC','CL','SC',
       'CL','SC','CL','CL','CL','SC','SC',
       'CL','SC','SC','CL','SC','CL',
       'SC','SC','CL','CL','SC','SC','SC',
       'CL','SC')                       
ids=data.frame(cell_name,cell_type)

sce=SingleCellExperiment(assays=list(counts=counts_aceto),
                         colData=DataFrame(label=cell_type))

df=perCellQCMetrics(sce)
qc.lib=isOutlier(df$sum,log=T,type='lower')
qc.nexprs=isOutlier(df$detected,log=T,type='lower')
qc.nexprs
qc.lib

discard=qc.lib|qc.nexprs

DataFrame(LibSize=sum(qc.lib),NExprs=sum(qc.nexprs),Total=sum(discard))
colData(sce)

reasons=quickPerCellQC(df)
filt_sce=sce[,!reasons$discard]
colData(filt)

clust=quickCluster(filt_sce,min.size=10)
clust

dconvfilt=computeSumFactors(filt_sce,cluster=clust,min.mean=0.1)
dconvfilt
summary(dconvfilt)

norm_sce=logNormCounts(dconvfilt)

logcounts(norm_sce)


dec.norm_sce=modelGeneVar(norm_sce)

dec.norm_sce=dec.norm_sce[order(dec.norm_sce$bio,decreasing = T),]
hvg.aceto=getTopHVGs(dec.norm_sce,var.threshold = 0)
str(hvg.aceto)

norm_sce_hvg=norm_sce[hvg.aceto,]


summed=aggregateAcrossCells(norm_sce,
                            id=colData(norm_sce)[,c('label')])
summed

