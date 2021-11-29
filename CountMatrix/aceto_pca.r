library(scater)
library(scran)
library(edgeR)

df=read.csv('Aceto_countMatrix.csv')

df=df[!duplicated(df$symbol),]
row.names(df)=df$symbol

df=df[,2:30]
df
cell_name=colnames(df)
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
batch=c(1,1,10,10,10,2,2,3,3,4,4,4,4,4,5,5,5,6,6,7,7,7,8,8,8,8,8,9,9)
sce=SingleCellExperiment(assays=list(counts=df),
                         colData=DataFrame(CL_or_SC=cell_type,patient=patient,batch=batch))

sce <- addPerCellQC(sce)
colData(sce)

plotColData(sce,x='sum',y='detected', colour_by = 'CL_or_SC')
per.feat=perFeatureQCMetrics(sce)

summary(per.feat$mean)
summary(per.feat$detected)

ave=calculateAverage(sce)
summary(ave)

agg.sce=aggregateAcrossCells(sce,ids=sce$CL_or_SC)
head(assay(agg.sce))
agg.sce

stats <- perCellQCMetrics(sce)

qc <- quickPerCellQC(stats)

sce=sce[,!qc$discard]

set.seed(100)
clust=quickCluster(sce,min.size=10)
sce <- computeSumFactors(sce,cluster=clust)
sce <- logNormCounts(sce)

dec=modelGeneVar(sce,block=sce$patient)

chosenhvgs <- getTopHVGs(dec,prop=0.1)

library(scater)

sce <- runPCA(sce,subset_row=chosenhvgs)

reducedDimNames(sce)

dim(reducedDim(sce,'PCA'))
plotReducedDim(sce,dimred = 'PCA',colour_by = 'patient')
sce <- runTSNE(sce,dimred="PCA")

plotReducedDim(sce,dimred="TSNE",colour_by='CL_or_SC')
