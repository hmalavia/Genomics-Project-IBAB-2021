library('scran')
library('scater')
library('batchelor')
counts_aceto=read.csv('Aceto_countMatrix.csv')

counts_aceto=counts_aceto[!duplicated(counts_aceto$symbol),]
row.names(counts_aceto)=counts_aceto$symbol
counts_aceto=counts_aceto[,2:30]

cell_name=colnames(counts_aceto)
cell_type=c('CL','SC','CL','CL','SC','CL','SC',
       'CL','SC','CL','CL','CL','SC','SC',
       'CL','SC','SC','CL','SC','CL',
       'SC','SC','CL','CL','SC','SC','SC',
       'CL','SC')
patient=c('Patient.1','Patient.1','Patient.10','Patient.10','Patient.10','Patient.2','Patient.2',
          'Patient.3','Patient.3','Patient.4','Patient.4','Patient.4','Patient.4','Patient.4',
          'Patient.5','Patient.5','Patient.5','Patient.6','Patient.6','Patient.7','Patient.7',
          'Patient.7','Patient.8','Patient.8','Patient.8','Patient.8','Patient.8','Patient.9',
          'Patient.9')
sample=c(1,1,10,10,10,2,2,3,3,4,4,4,4,4,5,5,5,6,6,7,7,7,8,8,8,8,8,9,9)
coldata=data.frame(CL_or_SC=cell_type,Patient=patient,sample=sample)
coldata$group=ifelse(coldata$CL_or_SC == 'CL',1,2)
coldata
sce=SingleCellExperiment(assays=list(counts=counts_aceto),
                         colData=coldata)

colData(sce)
df=perCellQCMetrics(sce)
reasons=quickPerCellQC(df)
sce[,reasons$discard]
filt_sce=sce[,!reasons$discard]
colData(filt_sce)

clust=quickCluster(filt_sce,min.size=10)
table(clust)

filt_sce <- computeSumFactors(filt_sce,cluster=clust)

filt_sce <- logNormCounts(filt_sce)

dec.sce=modelGeneVar(filt_sce,block=filt_sce$sample)

hvg.dec.sce <- getTopHVGs(dec.sce,var.threshold = 0)

sce.hvg=filt_sce[hvg.dec.sce,]


write.csv(logcounts(sce.hvg),file= 'aceto_HVGnormalizedcounts.csv',row.names = T,quote = F)

write.csv(counts(sce.hvg),file= 'aceto_HVGcounts.csv',row.names = T,quote = F)

set.seed(01001001)


summed=aggregateAcrossCells(sce.hvg,
                            id=colData(sce.hvg)[,c('CL_or_SC','sample')])
summed

current
library(edgeR)

y=DGEList(counts(sce.hvg),samples=colData(sce.hvg))

y$samples

keep <- filterByExpr(y,)

y <- y[keep,]
y$counts
y <- calcNormFactors(y)

y$samples

par(mfrow=c(3,2))

plotMD(y,column = 1)

plotMDS(cpm(y,log = T),
      col=ifelse(y$sample$CL_or_SC,"red","blue"))
CTCtype=factor(y$samples$CL_or_SC)
CTCtype <- relevel(CTCtype,ref='SC')
design <- model.matrix(~0+CTCtype,y$samples)
colnames(design) <- levels(CTCtype) 
design
y <- estimateDisp(y,design)
summary(y$trended.dispersion)
plotBCV(y)

fit <- glmQLFit(y,design,robust=T)
plotQLDisp(fit)


design
res <- glmQLFTest(fit,contrast = c(1,-1))
summary(decideTests(res))
topTags(res)

write.csv(res$table,file='aceto_CLvsSC_DGE_edgeR.csv',quote=F)

tr <- glmTreat(fit,coef = 2,lfc = 2)
summary(tr)
topTags(tr)
