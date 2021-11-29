library('SingleCellExperiment')
library('scran')
library('scater')
library('edgeR')
library('batchelor')
library('M3S')

counts=read.csv('GkountelaPatientCountMatrix.csv')

colnames(counts)[1] <- c('Geneid')

row.names(counts) <- counts$Geneid

counts <- counts[,2:ncol(counts)]

coldata <- sampleinfo

colnames(coldata)

coldata <- subset(coldata,select=c('!Sample_title','Sample_Type','Donor'))

rownames(coldata) <- coldata$`!Sample_title`

coldata <- coldata[,2:ncol(coldata)]

rownames(coldata)

colnames(counts)

keep <- rowSums(counts(sce)>0)>0

sce <- sce[keep,]

spikes <- grep("^ERCC-",rownames(counts),value = T)

spike_counts <- counts[spikes,]

notspikes <- setdiff(row.names(counts),spike_counts)

sce <- SingleCellExperiment(assays=list(counts=as.matrix(counts[notspikes,])),colData=coldata)

dim(counts(spike_counts))

SummarizedExperiment(spike_counts)
altExp(sce,'spike_ins') <- SummarizedExperiment(spike_counts,assays = list(counts=spike_counts))

counts(spike_counts)

df=perCellQCMetrics(sce)

df

reasons=quickPerCellQC(df,percent_subsets="altexps_spike_ins_percent",batch=sce$Donor)

discarded <- sce[,reasons$discard]@colData

discarded

write.table(discarded,'Gkountela_scran_discarded_samples.txt',sep='\t',quote = F)

filt_sce=sce[,!reasons$discard]

dim(colData(filt_sce))

colData(sce) <- cbind(colData(sce),df)

sce$discard <- reasons$discard

gridExtra::grid.arrange(
  plotColData(sce, x="Donor", y="sum", colour_by="discard") +
    scale_y_log10() + ggtitle("Total count") +
    theme(axis.text.x = element_text(angle = 90)),
  plotColData(sce, x="Donor", y="detected", colour_by="discard") +
    scale_y_log10() + ggtitle("Detected features") +
    theme(axis.text.x = element_text(angle = 90)),
  plotColData(sce, x="Donor", y="altexps_spike_ins_percent",
              colour_by="discard") + ggtitle("spike_ins percent") +
    theme(axis.text.x = element_text(angle = 90)),
  ncol=2
)

clust=quickCluster(filt_sce,min.size=10)
table(clust)

filt_sce <- computeSumFactors(filt_sce,cluster=clust)

filt_sce <- logNormCounts(filt_sce)

dec.sce=modelGeneVarWithSpikes(filt_sce,"spike_ins",block=filt_sce$Donor)

dec.sce

hvg <- dec.sce[dec.sce$bio>0,]

hvg <- hvg[hvg$p.value<0.05,]

hvg <- hvg[order(hvg$bio,decreasing = T),]

hvg
write.csv(hvg,'Gkountela_PatientSample_HVGs.csv',quote = F)

hvg.dec.sce <- getTopHVGs(dec.sce,var.threshold = 0)

sce.hvg=filt_sce[hvg.dec.sce,]

par(mfrow=c(3,3))

blocked.stats <- dec.sce$per.block

for (i in colnames(blocked.stats)) {
  current <- blocked.stats[[i]]
  plot(current$mean, current$total, main=i, pch=16, cex=0.5,
       xlab="Mean of log-expression", ylab="Variance of log-expression")
  curfit <- metadata(current)
  points(curfit$mean, curfit$var, col="red", pch=16)
  curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)
}
