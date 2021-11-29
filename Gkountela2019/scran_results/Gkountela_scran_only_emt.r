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

emt <- read.csv('Meschanchymal markers - Breast_cancer_only.csv',header = F)

colnames(emt)[1]='Geneid'

emt <- emt[,1]

emt <- as.data.frame(emt)

emt <- unique(emt)

head(emt)

emt_counts <- merge(counts,emt,by.x=0,by.y=1)

row.names(emt_counts) <- emt_counts$Row.names

emt_counts <- emt_counts[,2:ncol(emt_counts)]


sce.emt <- SingleCellExperiment(assays=list(counts=emt_counts),colData=coldata)

sce.emt

## Adding spike in counts as altexp ##

spikes <- grep(rownames(counts),pattern = "^ERCC-",value = T)

spikecounts <- counts[spikes,]

spikein <- SummarizedExperiment(list(counts=spikecounts))

altExp(sce.emt,'spikes') <- spikein

counts(sce.emt)

colData(sce.emt)

altExp(sce.emt,'spikes')

## QC ##

qc <- perCellQCMetrics(sce.emt)
qc

reasons <- quickPerCellQC(qc,sub.fields='altexps_spikes_percent')

colSums(as.matrix(reasons))

unfiltered <-(sce.emt)

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

sce.emt <- sce.emt[,!reasons$discard]

stats <- perCellQCMetrics(sce.emt)

colData(sce.emt) <- cbind(colData(sce.emt),stats)

gridExtra::grid.arrange(
  plotColData(sce.emt, x="Sample_Type", y="sum", 
              colour_by="Donor") + scale_y_log10() + ggtitle("Total count"),
  plotColData(sce.emt, x="Sample_Type", y="detected", 
              colour_by="Donor") + scale_y_log10() + ggtitle("Detected features"),
  plotColData(sce.emt, x="Sample_Type", y="altexps_spikes_percent", 
              colour_by="Donor") + ggtitle("ERCC percent"),
  nrow=1,
  ncol=3
)

### Normalization ###
set.seed(100)
clust <- quickCluster(sce.emt,min.size=5)

summary(clust)

emt.deconv <- calculateSumFactors(sce.emt,cluster=clust)

summary(emt.deconv)

emt.hvgs <- computeSpikeFactors(sce.emt,spikes = 'spikes')

summary(sizeFactors(emt.hvgs))


to.plot <- data.frame(
  Deconv=emt.deconv,
  libsize=librarySizeFactors(sce.emt),
  spike=sizeFactors(emt.hvgs),
  CTC_type=emt.hvgs$Sample_Type
)

plot <- ggplot(to.plot, aes(x=Deconv, y=libsize)) +
  geom_point() + facet_wrap(~CTC_type) + scale_x_log10() + 
  scale_y_log10() + geom_abline(intercept=0, slope=1, color="red")

plot  

plot2 <- ggplot(to.plot, aes(x=spike, y=libsize)) +
  geom_point() + facet_wrap(~CTC_type) + scale_x_log10() + 
  scale_y_log10() + geom_abline(intercept=0, slope=1, color="red")

plot2

sce.emt <- logNormCounts(sce.emt,size.factors=to.plot$Deconv)

assays(sce.emt)
head(logcounts(sce.emt))

### Feature Selection ###

emt.hvgs <- modelGeneVarWithSpikes(sce.emt,'spikes')

emt.hvgs[order(emt.hvgs$bio,decreasing=T),]

plot(emt.hvgs$mean, emt.hvgs$total, xlab="Mean of log-expression",
     ylab="Variance of log-expression")

fit.spike <- metadata(emt.hvgs)

metadata(emt.hvgs)

points(fit.spike$mean, fit.spike$var, col="red", pch=16)
curve(fit.spike$trend(x), col="dodgerblue", add=TRUE, lwd=2)

dim(emt.hvgs)

chosen.hvgs <- getTopHVGs(emt.hvgs,var.field = 'bio',var.threshold = 1,row.names = T)

### PCA ###

emt.hvgs <- sce.emt[chosen.hvgs,] 
emt.hvgs

set.seed(100)

sce.emt <- fixedPCA(sce.emt,subset.row = NULL)

emt.hvgs <- fixedPCA(emt.hvgs,subset.row = NULL)

dim(reducedDim(sce.emt,'PCA'))

dim(reducedDim(emt.hvgs,'PCA'))

percent.var.sce <- attr(reducedDim(sce.emt),'percentVar')

plot(percent.var.sce,log='y',xlab='PC',ylab='Variance explained(%)')

percent.var.hvg <- attr(reducedDim(emt.hvgs),'percentVar')

plot(percent.var.hvg,log='y',xlab='PC',ylab='Variance explained(%)')

#BiocManager::install("PCAtools")

library('PCAtools')

chosen.elbow <- findElbowPoint(percent.var.hvg)

chosen.elbow

abline(v=chosen.elbow,col='red')

plot(percent.var.hvg,log='y',xlab='PC',ylab='Variance explained(%)')

### Ploting PCA ###

plotReducedDim(emt.hvgs,dimred = 'PCA',colour_by = 'Sample_Type')

plotReducedDim(emt.hvgs,dimred = 'PCA',colour_by = 'Donor',shape_by = 'Sample_Type',point_size=2.5)

plotReducedDim(emt.hvgs,dimred = 'PCA',ncomponents = 4,colour_by = 'Donor',shape_by = 'Sample_Type',point_size=2.5)

### T-SNE ###
set.seed(100100)
emt.hvgs <- runTSNE(emt.hvgs,dimred='PCA')
plotReducedDim(emt.hvgs,dimred = 'TSNE',colour_by = 'Donor',shape_by = 'Sample_Type',point_size=3)

