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

##Create singleCellExperimentObject

sce <- SingleCellExperiment(assays = list(counts = counts),colData = coldata)

sce$Sample_Type <- factor(sce$Sample_Type)
sce$Donor <- factor(sce$Donor)
sce@colData

##Adding spike in counts as altexp

spikes <- grep(rownames(counts),pattern = "^ERCC-",value = T)

spikecounts <- counts[spikes,]

spikein <- SummarizedExperiment(list(counts=spikecounts))

altExp(sce,'spikes') <- spikein

##Adding epigene counts as altexps

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

##QC

qc <- perCellQCMetrics(sce)

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


##Normalization

clust <- quickCluster(sce,min.size=5)

summary(clust)

sce.deconv <- calculateSumFactors(sce,cluster=clust)

summary(sce.deconv)

sce.spike <- computeSpikeFactors(sce,spikes = 'spikes')

summary(sizeFactors(sce.spike))


to.plot <- data.frame(
  Deconv=sce.deconv,
  libsize=librarySizeFactors(sce),
  spike=sizeFactors(sce.spike),
  CTC_type=sce.spike$Sample_Type
)

plot <- ggplot(to.plot, aes(x=Deconv, y=libsize)) +
  geom_point() + facet_wrap(~CTC_type) + scale_x_log10() + 
  scale_y_log10() + geom_abline(intercept=0, slope=1, color="red")
plot  

logn
sce <- logNormCounts(sce,size.factors=to.plot$Deconv)

assays(sce)
logcounts(sce)
