##Load libraries
library('scran')
library('scater')
library('edgeR')

##read in count data
counts <- read.csv('Gkountela_Patient_unprocessed_rawcounts.csv')

rownames(counts) <- counts$Geneid

counts <- counts[,2:ncol(counts)]

counts <- as.matrix(counts)

head(counts)

##load sample data

coldata <- read.csv('Gkountela_coldata.csv')

row.names(coldata) <- coldata$X.Sample_title

coldata <- coldata[,2:ncol(coldata)]

### Create singleCellExperimentObject ###

epigenes <- read.csv('../CountMatrix/Epigenes_unique.csv')

epigene_counts <- merge(counts,epigenes,by.x=0,by.y=1)

row.names(epigene_counts) <- epigene_counts$Row.names

epigene_counts <- epigene_counts[,2:ncol(epigene_counts)]

epigene_counts

sce.epi <- SingleCellExperiment(assays=list(counts=epigene_counts),colData=coldata)

sce.epi

## Adding spike in counts as altexp ##

spikes <- grep(rownames(counts),pattern = "^ERCC-",value = T)

spikecounts <- counts[spikes,]

spikein <- SummarizedExperiment(list(counts=spikecounts))

altExp(sce.epi,'spikes') <- spikein

counts(sce.epi)

colData(sce.epi)

altExp(sce.epi,'spikes')

## QC ##

qc <- perCellQCMetrics(sce.epi)
qc

reasons <- quickPerCellQC(qc,sub.fields='altexps_spikes_percent')

colSums(as.matrix(reasons))

unfiltered <-(sce.epi)

colData(unfiltered) <- cbind(colData(unfiltered),qc)
colData(unfiltered)

unfiltered$discard <- reasons$discard


gridExtra::grid.arrange(
  plotColData(unfiltered, x="Sample_Type", y="sum", 
              colour_by="discard") + scale_y_log10() + ggtitle("Total count"),
  plotColData(unfiltered, x="Sample_Type", y="detected", 
              colour_by="discard") + scale_y_log10() + ggtitle("Detected features"),
  plotColData(unfiltered, x="Sample_Type", y="altexps_spikes_percent", 
              colour_by="discard") + ggtitle("ERCC percent"),
  nrow=1,
  ncol=3
)

sce.epi <- sce.epi[,!reasons$discard]

stats <- perCellQCMetrics(sce.epi)

colData(sce.epi) <- cbind(colData(sce.epi),stats)

gridExtra::grid.arrange(
  plotColData(sce.epi, x="Sample_Type", y="sum", 
              colour_by="Donor") + scale_y_log10() + ggtitle("Total count"),
  plotColData(sce.epi, x="Sample_Type", y="detected", 
              colour_by="Donor") + scale_y_log10() + ggtitle("Detected features"),
  plotColData(sce.epi, x="Sample_Type", y="altexps_spikes_percent", 
              colour_by="Donor") + ggtitle("ERCC percent"),
  nrow=1,
  ncol=3
)

### Normalization ###
set.seed(100)
clust <- quickCluster(sce.epi,min.size=5)

summary(clust)

epi.deconv <- calculateSumFactors(sce.epi,cluster=clust)

summary(epi.deconv)

epi.hvgs <- computeSpikeFactors(sce.epi,spikes = 'spikes')

summary(sizeFactors(epi.hvgs))


to.plot <- data.frame(
  Deconv=epi.deconv,
  libsize=librarySizeFactors(sce.epi),
  spike=sizeFactors(epi.hvgs),
  CTC_type=epi.hvgs$Sample_Type
)

plot <- ggplot(to.plot, aes(x=Deconv, y=libsize)) +
  geom_point() + facet_wrap(~CTC_type) + scale_x_log10() + 
  scale_y_log10() + geom_abline(intercept=0, slope=1, color="red")

plot  

plot2 <- ggplot(to.plot, aes(x=spike, y=libsize)) +
  geom_point() + facet_wrap(~CTC_type) + scale_x_log10() + 
  scale_y_log10() + geom_abline(intercept=0, slope=1, color="red")

plot2

sce.epi <- logNormCounts(sce.epi,size.factors=to.plot$Deconv)

assays(sce.epi)
head(logcounts(sce.epi))

### Feature Selection ###

epi.hvgs <- modelGeneVarWithSpikes(sce.epi,'spikes')

epi.hvgs[order(epi.hvgs$bio,decreasing=T),]

plot(epi.hvgs$mean, epi.hvgs$total, xlab="Mean of log-expression",
     ylab="Variance of log-expression")

fit.spike <- metadata(epi.hvgs)

metadata(epi.hvgs)

points(fit.spike$mean, fit.spike$var, col="red", pch=16)
curve(fit.spike$trend(x), col="dodgerblue", add=TRUE, lwd=2)

dim(epi.hvgs)

chosen.hvgs <- getTopHVGs(epi.hvgs,var.field = 'bio',var.threshold = 1,row.names = T)

### PCA ###

epi.hvgs <- sce.epi[chosen.hvgs,] 
epi.hvgs

set.seed(100)

sce.epi <- fixedPCA(sce.epi,subset.row = NULL)

epi.hvgs <- fixedPCA(epi.hvgs,subset.row = NULL)

dim(reducedDim(sce.epi,'PCA'))

dim(reducedDim(epi.hvgs,'PCA'))

percent.var.sce <- attr(reducedDim(sce.epi),'percentVar')

plot(percent.var.sce,log='y',xlab='PC',ylab='Variance explained(%)')

percent.var.hvg <- attr(reducedDim(epi.hvgs),'percentVar')

plot(percent.var.hvg,log='y',xlab='PC',ylab='Variance explained(%)')

#BiocManager::install("PCAtools")

library('PCAtools')

chosen.elbow <- findElbowPoint(percent.var.hvg)

chosen.elbow

abline(v=chosen.elbow,col='red')

### Ploting PCA ###

plotReducedDim(epi.hvgs,dimred = 'PCA',colour_by = 'Sample_Type')

plotReducedDim(epi.hvgs,dimred = 'PCA',colour_by = 'Donor',shape_by = 'Sample_Type',point_size=2.5)

plotReducedDim(epi.hvgs,dimred = 'PCA',ncomponents = 4,colour_by = 'Donor',shape_by = 'Sample_Type',point_size=2.5)

### T-SNE ###
set.seed(100100)
epi.hvgs <- runTSNE(epi.hvgs,dimred='PCA')
plotReducedDim(epi.hvgs,dimred = 'TSNE',colour_by = 'Donor',shape_by = 'Sample_Type',point_size=3)

