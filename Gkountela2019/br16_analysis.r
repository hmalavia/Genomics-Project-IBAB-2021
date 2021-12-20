#rm(list=ls())

## Br 16 ###

## Load libraries ##
library('scran')
library('scater')
library('edgeR')

## Load Count Data ##

counts <- read.csv('Gkountela_Patient_unprocessed_rawcounts.csv')

rownames(counts) <- counts$Geneid

counts <- counts[,2:ncol(counts)]

counts <- as.matrix(counts)

head(counts)

samples <- grep('Br16',colnames(counts))

br16 <- counts[,samples]
br16

## load sample data ##

coldata <- read.csv('Gkountela_coldata.csv')

rownames(coldata) <- coldata$X

coldata <- coldata[,2:ncol(coldata)]

coldata

coldata16 <- coldata[samples,]

coldata16

coldata16 <- coldata16[order(coldata16$Sample_Type),]

generate_SampleID <- function(cd){
  k=1
  for (i in 1:nrow(cd)){
    if (i != nrow(cd)){
      if (cd$Sample_Type[i] == cd$Sample_Type[i+1]){
        type <- cd$Sample_Type[i]
        cd$Sample_ID[i] = paste(type,k,sep ="_" )
        k=k+1
      } else {
        cd$Sample_ID[i] = paste(type,k,sep ="_" )
        k=1
        next  
      }
    }
    else{
      cd$Sample_ID[i] = paste(type,k,sep ="_" )
    }
  }
  return(cd)
}

coldata16 <- generate_SampleID(coldata16)

#write.csv(coldata16,'coldata16.csv',quote = F)

#### Create singleCellExperimentObject ####

sce <- SingleCellExperiment(assays = list(counts = br16),colData = coldata16)

sce$Sample_Type <- factor(sce$Sample_Type)
sce$Donor <- factor(sce$Donor)
#sce@colData

#### Adding spike in counts as altexp ####

spikes <- grep(rownames(br16),pattern = "^ERCC-",value = T)

spikecounts <- br16[spikes,]

spikein <- SummarizedExperiment(list(counts=spikecounts))

altExp(sce,'spikes') <- spikein

## Adding epigene counts as altexps ##

epigenes <- read.csv('../CountMatrix/Epigenes_unique.csv')

head(epigenes)

epigene_counts <- merge(br16,epigenes,by.x=0,by.y=1)

row.names(epigene_counts) <- epigene_counts$Row.names

epigene_counts <- epigene_counts[,2:ncol(epigene_counts)]

epi_counts <- SummarizedExperiment(list(counts=epigene_counts))

altExp(sce,'epi') <- epi_counts

## Adding emtgene counts as altexps

emt <- read.csv('Meschanchymal markers - Breast_cancer_only.csv',header = F)

colnames(emt)[1]='Geneid'

emt <- emt[,1]

emt <- as.data.frame(emt)

emt <- unique(emt)

head(emt)

emt_counts <- merge(br16,emt,by.x=0,by.y=1)

row.names(emt_counts) <- emt_counts$Row.names

emt_counts <- emt_counts[,2:ncol(emt_counts)]

emt_counts <- SummarizedExperiment(list(counts=emt_counts))

altExp(sce,'emt') <- emt_counts

altExp(sce,'emt')

#### QC ####

qc <- perCellQCMetrics(sce)
#qc

reasons <- quickPerCellQC(qc,sub.fields='altexps_spikes_percent')

colSums(as.matrix(reasons))

unfiltered <- sce

colData(unfiltered) <- cbind(colData(unfiltered),qc)
colData(unfiltered)

unfiltered$discard <- reasons$discard

colData(unfiltered[,unfiltered$discard])[1,ncol(colData(unfiltered))] = F

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

#sce[,unfiltered$discard]

colData(sce[,unfiltered$discard])

sce <- sce[,!unfiltered$discard]

sce

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

##### Normalization #####
set.seed(100)

lib.sf <- librarySizeFactors(sce)

summary(lib.sf)

sce <- logNormCounts(sce,size.factors=lib.sf)

assays(sce)

##### Feature Selection #####

sce.hvgs <- modelGeneVarWithSpikes(sce,'spikes')

sce.hvgs[order(sce.hvgs$bio,decreasing=T),]

plot(sce.hvgs$mean, sce.hvgs$total, xlab="Mean of log-expression",
     ylab="Variance of log-expression")

fit.spike <- metadata(sce.hvgs)

points(fit.spike$mean, fit.spike$var, col="red", pch=16)
curve(fit.spike$trend(x), col="dodgerblue", add=TRUE, lwd=2)

dim(sce.hvgs)

chosen.hvgs <- getTopHVGs(sce.hvgs,var.field = 'bio',var.threshold = 1,row.names = T)

#length(chosen.hvgs)

##### PCA #####

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

#BiocManager::install("PCAtools")

library('PCAtools')

chosen.elbow <- findElbowPoint(percent.var.hvg)

abline(v=chosen.elbow,col='red')

##### Ploting PCA #####

plotReducedDim(sce.hvgs,dimred = 'PCA',colour_by = 'Sample_Type',point_size=2.5)

plotReducedDim(sce.hvgs,dimred = 'PCA',colour_by = 'Sample_Type',point_size=2.5,text_by = 'Sample_ID',text_size = 2.5)

plotReducedDim(sce.hvgs,dimred = 'PCA',ncomponents = 4,colour_by = 'Sample_Type',shape_by = 'Sample_Type', point_size=2.5)

plotReducedDim(sce.hvgs,dimred = 'PCA', ncomponents=c(1,3), colour_by='Sample_Type', shape_by = 'Sample_Type',point_size=2.5,text_by = 'Sample_ID',text_size = 2.5)

plotReducedDim(sce.hvgs,dimred = 'PCA', ncomponents=c(1,3), colour_by='Sample_Type', shape_by = 'Sample_Type',point_size=2.5)


##### T-SNE #####
set.seed(100100)

sce.hvgs <- runTSNE(sce.hvgs,dimred='PCA')

plotReducedDim(sce.hvgs,dimred = 'TSNE',colour_by = 'Sample_Type',shape_by = 'Sample_Type', point_size=2.5)

plotReducedDim(sce.hvgs,dimred = 'TSNE',colour_by = 'Sample_Type',shape_by = 'Sample_Type', point_size=2.5,text_by = 'Sample_ID',text_size = 2.5)


#### Varying Perplexity for T-SNE ####
set.seed(100)
sce.hvgs <- runTSNE(sce.hvgs, dimred="PCA", perplexity=5)
out5 <- plotReducedDim(sce.hvgs, dimred="TSNE",
                       colour_by="Sample_Type",shape_by = 'Sample_Type', point_size=2.5) + ggtitle("perplexity = 5")

set.seed(100)
sce.hvgs <- runTSNE(sce.hvgs, dimred="PCA", perplexity=20)
out20 <- plotReducedDim(sce.hvgs, dimred="TSNE",
                        colour_by="Sample_Type",shape_by = 'Sample_Type', point_size=2.5) + ggtitle("perplexity = 20")

set.seed(100)
sce.hvgs <- runTSNE(sce.hvgs, dimred="PCA", perplexity=80)
out80 <- plotReducedDim(sce.hvgs, dimred="TSNE", 
                        colour_by="Sample_Type",shape_by = 'Sample_Type', point_size=2.5) + ggtitle("perplexity = 80")

gridExtra::grid.arrange(out5, out20, out80, ncol=3)

#### Graph based Clustering ####

library(bluster)

clust.walktrap_rank <- clusterCells(sce.hvgs,use.dimred='PCA',
                                    BLUSPARAM = SNNGraphParam(k=4, type = 'rank', cluster.fun = 'walktrap')) 

table(clust.walktrap_rank)

clust.louvain_jaccard <- clusterCells(sce.hvgs,use.dimred='PCA',
                                      BLUSPARAM = SNNGraphParam(k=2, type = 'jaccard', cluster.fun = 'louvain')) 

table(clust.louvain_jaccard)

colLabels(sce.hvgs) <- clust.louvain_jaccard

clust.louvain_jaccard

plotReducedDim(sce.hvgs,'PCA',colour_by = 'label',shape_by='Sample_Type',point_size=2.5)

#### DE analysis using edgeR ####

#library('edgeR')

sce.hvgs$Sample_ID
groups <- c(1,1,1,2,2)

y <- DGEList(counts(sce.hvgs),samples=colData(sce.hvgs),group = groups )

y$counts

nrow(y$samples['Sample_Name'])

y <- calcNormFactors(y)

length(y$samples)

sce.hvgs$sizeFactor

ncol(y)

## MD plot ##
par(mfrow=c(3,4))
for (i in seq_len(ncol(y))) {
  plotMD(y, column=i)
}

y <- estimateDisp(y)

et <- exactTest(y)

toptags <- topTags(et)

toptags

summary(decideTestsDGE(et))

result <- et$table

result <- result[order(-result$logFC),]

write.csv(result,'Br16_edgeR_result.csv',quote = F)

# upreg <- result[result$logFC > 2,]
# 
# upregsig <- upreg[upreg$PValue <= 0.05,]
# 
# Downregs <- result[result$logFC < -2,]
# 
# Downregsig <- Downregs[Downregs$PValue <= 0.05,]

upreg <- result[result$logFC >= 1,]

upregsig <- upreg[upreg$PValue <= 0.05,]

Downregs <- result[result$logFC <= -1,]

Downregsig <- Downregs[Downregs$PValue <= 0.05,]

upreg <- upreg[order(-upreg$logFC),]

upregsig <- upregsig[order(-upregsig$logFC),]

write.csv(upregsig,'scran_results/Br16_upregulated_genes_ALL_edgeR.csv',quote = F)

Downregs <- Downregs[order(Downregs$logFC),]

Downregsig <- Downregsig[order(Downregsig$logFC),]

write.csv(Downregsig,'scran_results/Br16_downregulated_genes_ALL_edgeR.csv',quote = F)

logcpm <- cpm(counts(sce.hvgs),log=T)

logUp <- logcpm[rownames(upregsig),]

logDown <- logcpm[rownames(Downregsig),]

new_col_names <- sce.hvgs$Sample_ID

new_col_names

colnames(logUp) <- new_col_names 

colnames(logDown) <- new_col_names 

Updown <- rbind(logUp[0:10,],logDown)

## HeatMap ##

library('pheatmap')

pheatmap(as.matrix(Updown))

