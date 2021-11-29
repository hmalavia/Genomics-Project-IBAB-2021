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

## load sample data ##

coldata <- read.csv('Gkountela_coldata.csv')

row.names(coldata) <- coldata$X.Sample_title

coldata <- coldata[,2:ncol(coldata)]

## Create singleCellExperimentObject ##

sce <- SingleCellExperiment(assays = list(counts = counts),colData = coldata)

sce$Sample_Type <- factor(sce$Sample_Type)
sce$Donor <- factor(sce$Donor)
sce@colData

##Adding spike in counts as altexp

spikes <- grep(rownames(counts),pattern = "^ERCC-",value = T)

spikecounts <- counts[spikes,]

spikein <- SummarizedExperiment(list(counts=spikecounts))

altExp(sce,'spikes') <- spikein

## Adding epigene counts as altexps ##

epigenes <- read.csv('../CountMatrix/Epigenes_unique.csv')

head(epigenes)

epigene_counts <- merge(counts,epigenes,by.x=0,by.y=1)

row.names(epigene_counts) <- epigene_counts$Row.names

epigene_counts <- epigene_counts[,2:ncol(epigene_counts)]

epi_counts <- SummarizedExperiment(list(counts=epigene_counts))

altExp(sce,'epi') <- epi_counts

##Adding emtgene counts as altexps

emt <- read.csv('Meschanchymal markers - Breast_cancer_only.csv',header = F)

colnames(emt)[1]='Geneid'

emt <- emt[,1]

emt <- as.data.frame(emt)

emt <- unique(emt)

head(emt)

emt_counts <- merge(counts,emt,by.x=0,by.y=1)

row.names(emt_counts) <- emt_counts$Row.names

emt_counts <- emt_counts[,2:ncol(emt_counts)]

emt_counts <- SummarizedExperiment(list(counts=emt_counts))

altExp(sce,'emt') <- emt_counts

altExp(sce,'emt')

## QC ##

qc <- perCellQCMetrics(sce)
qc

reasons <- quickPerCellQC(qc,sub.fields='altexps_spikes_percent')

colSums(as.matrix(reasons))

unfiltered <- sce

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
  plotColData(unfiltered,x='Sample_Type',y='altexps_epi_percent')+scale_y_log10()+ggtitle('Epigenes_percentage'),
  plotColData(unfiltered,x='Sample_Type',y='altexps_emt_percent')+scale_y_log10()+ggtitle('EMTgenes_percentage'),
  nrow=2,
  ncol=3
)

sce <- sce[,!reasons$discard]

stats <- perCellQCMetrics(sce)

colData(sce) <- cbind(colData(sce),stats)


gridExtra::grid.arrange(
  plotColData(sce, x="Sample_Type", y="sum", 
              colour_by="Donor") + scale_y_log10() + ggtitle("Total count"),
  plotColData(sce, x="Sample_Type", y="detected", 
              colour_by="Donor") + scale_y_log10() + ggtitle("Detected features"),
  plotColData(sce, x="Sample_Type", y="altexps_spikes_percent", 
              colour_by="Donor") + ggtitle("ERCC percent"),
  plotColData(sce,x='Sample_Type',y='altexps_epi_percent',colour_by='Donor')+scale_y_log10()+ggtitle('Epigenes_percentage'),
  plotColData(sce,x='Sample_Type',y='altexps_emt_percent',colour_by='Donor')+scale_y_log10()+ggtitle('EMTgenes_percentage'),
  nrow=2,
  ncol=3
  )


### Normalization ###
set.seed(100)
clust <- quickCluster(sce,min.size=5)

summary(clust)

sce.deconv <- calculateSumFactors(sce,cluster=clust)

summary(sce.deconv)

sce.hvgs <- computeSpikeFactors(sce,spikes = 'spikes')

summary(sizeFactors(sce.hvgs))


to.plot <- data.frame(
  Deconv=sce.deconv,
  libsize=librarySizeFactors(sce),
  spike=sizeFactors(sce.hvgs),
  CTC_type=sce.hvgs$Sample_Type
)

plot <- ggplot(to.plot, aes(x=Deconv, y=libsize)) +
  geom_point() + facet_wrap(~CTC_type) + scale_x_log10() + 
  scale_y_log10() + geom_abline(intercept=0, slope=1, color="red")

plot  

plot2 <- ggplot(to.plot, aes(x=spike, y=libsize)) +
  geom_point() + facet_wrap(~CTC_type) + scale_x_log10() + 
  scale_y_log10() + geom_abline(intercept=0, slope=1, color="red")

plot2

sce <- logNormCounts(sce,size.factors=to.plot$Deconv)

assays(sce)
head(logcounts(sce))

### Feature Selection ###

sce.hvgs <- modelGeneVarWithSpikes(sce,'spikes')

sce.hvgs[order(sce.hvgs$bio,decreasing=T),]

plot(sce.hvgs$mean, sce.hvgs$total, xlab="Mean of log-expression",
     ylab="Variance of log-expression")

fit.spike <- metadata(sce.hvgs)

points(fit.spike$mean, fit.spike$var, col="red", pch=16)
curve(fit.spike$trend(x), col="dodgerblue", add=TRUE, lwd=2)

dim(sce.hvgs)

chosen.hvgs <- getTopHVGs(sce.hvgs,var.field = 'bio',var.threshold = 1,row.names = T)

### PCA ###

sce.hvgs <- sce[chosen.hvgs,] 
sce.hvgs

set.seed(100)

sce <- fixedPCA(sce,subset.row = NULL)

sce.hvgs <- fixedPCA(sce.hvgs,subset.row = NULL)

dim(reducedDim(sce,'PCA'))

dim(reducedDim(sce.hvgs,'PCA'))

percent.var.sce <- attr(reducedDim(sce),'percentVar')

plot(percent.var.sce,log='y',xlab='PC',ylab='Variance explained(%)')

percent.var.hvg <- attr(reducedDim(sce.hvgs),'percentVar')

plot(percent.var.hvg,log='y',xlab='PC',ylab='Variance explained(%)')

BiocManager::install("PCAtools")

library('PCAtools')

chosen.elbow <- findElbowPoint(percent.var.hvg)

abline(v=chosen.elbow,col='red')

### Ploting PCA ###

plotReducedDim(sce.hvgs,dimred = 'PCA',colour_by = 'Sample_Type')

plotReducedDim(sce.hvgs,dimred = 'PCA',colour_by = 'Donor',shape_by = 'Sample_Type')

plotReducedDim(sce.hvgs,dimred = 'PCA',ncomponents = 3,colour_by = 'Donor',shape_by = 'Sample_Type')

### T-SNE ###
set.seed(100100)
sce.hvgs <- runTSNE(sce.hvgs,dimred='PCA')
plotReducedDim(sce.hvgs,dimred = 'TSNE',colour_by = 'Donor',shape_by = 'Sample_Type')

### Varying Perplexity for T-SNE ###
set.seed(100)
sce.hvgs <- runTSNE(sce.hvgs, dimred="PCA", perplexity=5)
out5 <- plotReducedDim(sce.hvgs, dimred="TSNE",
                       colour_by="Donor",shape_by = 'Sample_Type') + ggtitle("perplexity = 5")

set.seed(100)
sce.hvgs <- runTSNE(sce.hvgs, dimred="PCA", perplexity=20)
out20 <- plotReducedDim(sce.hvgs, dimred="TSNE",
                        colour_by="Donor",shape_by = 'Sample_Type') + ggtitle("perplexity = 20")

set.seed(100)
sce.hvgs <- runTSNE(sce.hvgs, dimred="PCA", perplexity=80)
out80 <- plotReducedDim(sce.hvgs, dimred="TSNE", 
                        colour_by="Donor",shape_by = 'Sample_Type') + ggtitle("perplexity = 80")

gridExtra::grid.arrange(out5, out20, out80, ncol=3)

### Clustering ###

## Graph based Clustering ##

#BiocManager::install("bluster")

library('bluster')

graph.clusters <- clusterCells(sce.hvgs,use.dimred = 'PCA')

table(graph.clusters)

colLabels(sce.hvgs) <- graph.clusters

plotReducedDim(sce.hvgs,'PCA',colour_by = 'label',shape_by = 'Sample_Type',text_by = 'label')

graph.clusters2 <- clusterCells(sce.hvgs, use.dimred="PCA", 
                             BLUSPARAM=SNNGraphParam(k=5, type="rank", cluster.fun="walktrap"))
table(graph.clusters2)

colLabels(sce.hvgs) <- graph.clusters2

plotReducedDim(sce.hvgs,'TSNE',colour_by = 'label',shape_by = 'Sample_Type',text_by = 'Donor')
