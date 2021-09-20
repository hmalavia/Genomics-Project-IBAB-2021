if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scran")

BiocManager::install("scater")

library(SingleCellExperiment)
library(scran)

count_matrix=read.csv('Aceto_countMatrix.csv')

countsMatrix=count_matrix[!duplicated(count_matrix$symbol),]

row.names(countsMatrix)=countsMatrix$symbol

countsMatrix=countsMatrix[,2:30]

countsMatrix=as.matrix(countsMatrix)

sce=SingleCellExperiment(assays=list(counts=countsMatrix))

set.seed(100)
clust=quickCluster(sce,min.size=10)

table(clust)

sce <- computeSumFactors(sce, clusters=clust,min.mean=0.1)

sce<- logNormCounts(sce)

assayNames(sce)

write.csv(logcounts(sce),file = "aceto_normalizedcounts.csv",row.names = F,quote = F)


dec.aceto=modelGeneVar(sce)

dec.aceto

fit.aceto <- metadata(dec.aceto)

plot(fit.aceto$mean,fit.aceto$var,xlab='mean of log-expression',ylab='Variance of log-expression')
curve(fit.aceto$trend(x),col='dodgerblue',add=T,lwd=2)

dec.aceto=dec.aceto[order(dec.aceto$bio,decreasing = T),]

write.csv(dec.aceto,file="DGE/aceto_HVG.csv",row.names=T,quote=F)
