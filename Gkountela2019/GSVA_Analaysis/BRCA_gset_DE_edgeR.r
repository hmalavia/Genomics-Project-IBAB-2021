library('edgeR')

#### Load Count Data ####

counts <- read.csv('Gkountela_Patient_unprocessed_rawcounts.csv')

rownames(counts) <- counts$Geneid

counts <- counts[,2:ncol(counts)]

counts <- as.matrix(counts)

head(counts)

br61 <- grep(colnames(counts),'Br61')

br61 <- counts[,grep('Br61',colnames(counts))]
br61

#### Subsetting only EMT geneset Counts ####
counts <- br61[rownames(br61_common_genes_proliferation),colnames(br61_common_genes_proliferation)]

#### DE analysis using edgeR ####

sample_data <- read.csv('br61_scran_filtered_coldata.csv')

sample_data

y <- DGEList(counts,samples=sample_data,group = groups )

# y$counts
# 
# nrow(y$samples['Sample_Name'])

y <- calcNormFactors(y)

y <- estimateDisp(y, design)

group <- factor(c(1,1,1,1,1,1,2,2,2,2,2,2))
design <- model.matrix(~group)
fit <- glmQLFit(y, design)

qlf.2vs1 <- glmQLFTest(fit, coef=2)
toptags_qlf <- topTags(qlf.2vs1)

toptags_qlf

summary(decideTestsDGE(qlf.2vs1))

result <- qlf.2vs1$table

result <- result[order(-result$logFC),]

write.csv(result,'GSVA_Analaysis/Br61_Proliferation_gset_DEG.csv',quote = F)

# upreg <- result[result$logFC > 2 & result$PValue < 0.05,]
# 
# Downregs <- result[result$logFC < -2 &  result$PValue < 0.05,]

Sigup <- result[result$logFC > 1 & result$PValue < 0.05,]

SigDown <- result[result$logFC < -1 & result$PValue < 0.05,]

SigUpDown <- rbind(Sigup,SigDown)

SigUpDown <- SigUpDown[order(-SigUpDown$logFC),]

write.csv(SigUpDown,'GSVA_Analaysis/Br61_Proliferation_gset_SigUpDown_DEG.csv',quote=F)

logcpm <- read.csv('br61_scran_norm_counts.csv',row.names = "X")

logUpDown <- logcpm[rownames(SigUpDown),]

new_col_names <- y$samples$Sample_Name

colnames(logUpDown) <- new_col_names 

#### HeatMap ####

#library('pheatmap')

pheatmap(as.matrix(logUpDown),cluster_rows = T,cluster_cols = T,scale = 'row',main = "Differentialy expressed Proliferation Genes in Br61",xlab='Samples',ylab='Genes')


