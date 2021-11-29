#install.packages('Seurat')
library('Seurat')


counts <- read.csv('Gkountela_Patient_unprocessed_rawcounts.csv')

colnames(counts)

sampleinfo <- read.csv('Gkountela_coldata.csv')

head(sampleinfo)

rownames(counts) <- counts$Geneid

counts <- counts[2:ncol(counts)]

dups <- duplicated(counts)

counts <- counts[!dups,]

head(counts)


obj <- CreateSeuratObject(counts = counts,project = 'gkountela')

obj[['percent.ercc']] <- PercentageFeatureSet(obj,pattern = '^ERCC')

obj[['percent.ercc']]

obj[['orig.ident']] <- sampleinfo$Sample_Type

obj[['donor']] <- sampleinfo$Donor

obj@meta.data

VlnPlot(obj,features=c('nFeature_RNA','nCount_RNA','percent.ercc'),ncol=3,group.by = 'orig.ident')

#QC
#Removing cells with high spike ins low features

obj <- subset(obj,subset = nFeature_RNA > 1000 & nCount_RNA > 80000 & percent.ercc < 50 )

nrow(obj@meta.data)

VlnPlot(obj,features=c('nFeature_RNA','nCount_RNA','percent.ercc'),ncol=3,group.by = 'orig.ident')

#Normalization

obj <- NormalizeData(obj)

#Finding Variable Features(HVGs)

obj <- FindVariableFeatures(obj,selection.method = 'vst',nfeatures = 3000)

top10 <- head(VariableFeatures(obj,10)

top10

epigenes <- read.csv('../CountMatrix/Epigenes_unique.csv')

var_genes <- as.data.frame(VariableFeatures(obj))

var_epi <- merge(var_genes,epigenes,by.x = 1,by.y = 1)

write.csv(var_epi,'Gkountela_Epi_HVGs_seurat.csv',quote = F,row.names = F)

obj@meta.data

plot1 <- VariableFeaturePlot(obj)
plot2 <- LabelPoints(plot = plot1,points=c(var_epi[1:nrow(var_epi),]),repel = T,xnudge = 0.5,ynudge = 0.5)

plot2
#Scaling

obj <- ScaleData(obj,features = VariableFeatures(obj))

#PCA

obj <- RunPCA(obj,features = VariableFeatures(obj),approx=F)

VizDimLoadings(obj,dims=1:2,reduction = 'pca')

unique(obj@meta.data$donor)

pl1 <- DimPlot(obj,reduction = 'pca',cells=obj@meta.data$donor == 'Br37',group.by = 'orig.ident')
pl2 <- DimPlot(obj,reduction = 'pca',cells=obj@meta.data$donor == 'Br23',group.by = 'orig.ident')
pl3 <- DimPlot(obj,reduction = 'pca',cells=obj@meta.data$donor == 'Br39',group.by = 'orig.ident')
pl4 <- DimPlot(obj,reduction = 'pca',cells=obj@meta.data$donor == 'Br11',group.by = 'orig.ident')
pl5 <- DimPlot(obj,reduction = 'pca',cells=obj@meta.data$donor == 'Br53',group.by = 'orig.ident')
pl6 <- DimPlot(obj,reduction = 'pca',cells=obj@meta.data$donor == 'Br57',group.by = 'orig.ident')
pl7 <- DimPlot(obj,reduction = 'pca',cells=obj@meta.data$donor == 'Br61',group.by = 'orig.ident')
pl8 <- DimPlot(obj,reduction = 'pca',cells=obj@meta.data$donor == 'Br7',group.by = 'orig.ident')
pl9 <- DimPlot(obj,reduction = 'pca',cells=obj@meta.data$donor == 'Br16',group.by = 'orig.ident')

pl2+pl3+pl4+pl6+pl7+pl8+pl9

DimPlot(obj,dims=c(1,5),reduction = 'pca',group.by = c('donor','orig.ident'))
pl10
pl7

DimHeatmap(obj,dims = 1,balanced = T)
(obj <- JackStraw(obj,num.replicate = 100)(obj <- ScoreJackStraw(obj,dims = 1:20)

JackStrawPlot(obj

ElbowPlot(obj


