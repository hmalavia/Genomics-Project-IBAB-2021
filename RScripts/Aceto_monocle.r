#BiocManager::install("monocle")
library('monocle')
library('scuttle')
library('scran')
library('scater')

####HVG_TPM_NORMALIZATION####

counts=read.csv('aceto_HVGcounts.csv')
hvg_ano=read.csv('Aceto_HVGs_anotation.csv')
hvg_ano$GeneLen=hvg_ano$Gene.end..bp. - hvg_ano$Gene.start..bp.
counts_ano=merge(counts,hvg_ano,by.x='X',by.y='Gene.name')
rowdata=subset(counts_ano,select = c('X','GeneLen'))

counts_sce=counts_ano[,1:25]
uniq.counts_sce=counts_sce[!duplicated(counts_sce$X),]
row.names(uniq.counts_sce)=uniq.counts_sce$X
uniq.counts_sce=uniq.counts_sce[,2:ncol(uniq.counts_sce)]

sce=SingleCellExperiment(assays=list(counts=uniq.counts_sce))
counts(sce)

genes <- row.names(uniq.counts_sce)
genes <- as.data.frame(genes)
colnames(genes)[1]='Gene.name'

lengths=merge(genes,hvg_ano,by.x='Gene.name',by.y='Gene.name')
lengths=lengths[!duplicated(lengths$Gene.name),]
lengths=subset(lengths,select=c('Gene.name','GeneLen'))
aceto_hvgfpkm=calculateFPKM(sce,lengths = lengths$GeneLen)
write.csv(aceto_hvgfpkm,file='Aceto_HVG_FPKMcounts.csv',quote=F)
####MONOCLE####
exprs <- aceto_hvgfpkm
cell=colnames(exprs)
cell
celltype=c(ifelse(grepl("CL",cell),"CL","SC"))
pd=data.frame(Cell=cell,Cell_type=celltype)
row.names(pd)=pd$Cell
pd=subset(pd,select=c('Cell_type'))
fd=as.data.frame(rownames(exprs))
row.names(fd)=fd$`rownames(exprs)`
fd$gene_short_name=rownames(exprs)
row.names(fd)=fd$`rownames(exprs)`
fd=subset(fd,select=c('gene_short_name'))
pd <- new('AnnotatedDataFrame',data=pd)
fd <- new('AnnotatedDataFrame',data=fd)

cds=newCellDataSet(exprs,phenoData=pd,featureData = fd,
                   expressionFamily = tobit())
cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds)

disp <- dispersionTable(cds)

diff_test_result=differentialGeneTest(cds,fullModelFormulaStr = '~Cell_type')
top_res=subset(diff_test_result,qval < 0.1)
top_res <- top_res[,c('pval','qval')]
top_res <- top_res[order(top_res$pval),]

write.csv(top_res,'Aceto_CL_VS_SC_ALL_DEG_monocle.csv',quote=F,row.names = F)

epigenes=read.csv('Epigenes_unique.csv')

epi_res=merge(epigenes,top_res,by.x='HGNC.approved.symbol',by.y=0)

epi_res <- epi_res[order(epi_res$pval),]

write.csv(epi_res,'Aceto_CL_VS_SC_EPI_DEG_monocle.csv',quote = F,row.names = F)
top_5=top_res[1:5,]
top5Gs <- row.names(top_5)
plot_genes_jitter(cds['MAPKAPK3',],grouping = 'Cell_type')
