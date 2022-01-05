# BiocManager::install("GSVA")
library(GSVA)
library(stringr)
library(ggplot2)


file_paths <- list.files('CancerSEA_markers',pattern = "BRCA_*",full.names = T)

file_names <- list.files('CancerSEA_markers',pattern = "BRCA_*",full.names = F)

file_names

set_names <- tools::file_path_sans_ext(file_names)

brca_sets <- lapply(file_paths,read.csv,header=F)

names(brca_sets) <- set_names

names(brca_sets)

brca_sets_up <- lapply(brca_sets,function(x) x[x$V4 != 'negative',"V2"])

brca_sets_up <- lapply(brca_sets_up, function(x) x[names(x) %in% c("V2")])

brca_sets_up1 <- lapply(brca_sets_up, function(x) as.list(x))

brca_sets_up1 <- lapply(brca_sets_up, function(x) x[x != "Symbol"])


expr11 <- read.csv('br11_scran_norm_counts.csv')

colnames(expr11)[1] <- "GeneName"

row.names(expr11) <- expr11$GeneName

expr11$GeneName <- NULL

expr11 <- as.matrix(expr11)

es11_ssgsea <- gsva(expr11,brca_sets_up1,kcdf="Gaussian",method='ssgsea')

es11_gsva <- gsva(expr11,brca_sets_up1,kcdf="Gaussian",mx.diff=T,abs.ranking=T)

sdata11 <- read.csv('br11_filtered_coldata.csv')

colnames(es11_gsva) <- sdata$Sample_ID

es11_gsva_trans <- t(es11_gsva)

write.csv(es11_gsva_trans,'GSVA_analysis_br11_BRCA_gsets',quote = F)

expr61 <- read.csv('br61_scran_norm_counts.csv')

colnames(expr61)[1] <- "GeneName"

row.names(expr61) <- expr61$GeneName

expr61$GeneName <- NULL

expr61 <- as.matrix(expr61)

es61_ssgsea <- gsva(expr61,brca_sets_up1,kcdf="Gaussian",method='ssgsea')

es61_gsva <- gsva(expr61,brca_sets_up1,kcdf="Gaussian",mx.diff=T,abs.ranking=T)

sdata61 <- read.csv('br61_scran_filtered_coldata.csv')

colnames(es61_gsva) <- sdata61$Sample_Name
colnames(es61_ssgsea) <- sdata61$Sample_Name

es61_gsva_trans <- t(es61_gsva)

write.csv(es61_gsva_trans,'GSVA_analysis_br61_BRCA_gsets',quote = F)


brca_sets_down <- lapply(brca_sets,function(x) x[x$V4 != 'positive',"V2"])

#brca_sets_down <- lapply(brca_sets_down, function(x) x[names(x) %in% c("V2")])

#brca_sets_down1 <- lapply(brca_sets_down, function(x) as.list(x))

brca_sets_down1 <- lapply(brca_sets_down, function(x) x[x != "Symbol"])

es11_ssgsea_down <- gsva(expr11,brca_sets_down1,kcdf="Gaussian",method='ssgsea')

es11_gsva_down <- gsva(expr11,brca_sets_down1,kcdf="Gaussian",mx.diff=F)

sdata11 <- read.csv('br11_filtered_coldata.csv')

colnames(es11_gsva_down) <- sdata11$Sample_ID
colnames(es11_ssgsea_down) <- sdata11$Sample_ID

expr61 <- read.csv('br61_scran_norm_counts.csv')

colnames(expr61)[1] <- "GeneName"

row.names(expr61) <- expr61$GeneName

expr61$GeneName <- NULL

expr61 <- as.matrix(expr61)

es61_ssgsea_down <- gsva(expr61,brca_sets_down1,kcdf="Gaussian",method='ssgsea')

es61_gsva_down <- gsva(expr61,brca_sets_down1,kcdf="Gaussian",mx.diff=F)

sdata61 <- read.csv('br61_scran_filtered_coldata.csv')

colnames(es61_gsva_down) <- sdata61$Sample_Name
colnames(es61_ssgsea_down) <- sdata61$Sample_Name

barplot(es11_gsva['BRCA_EMT',],xlab = 'Cells',ylab = 'Enrichment Score',main = 'Enrichment Score of EMT geneset')

barplot(es61_gsva['BRCA_EMT',],xlab = 'Cells',ylab = 'Enrichment Score',main = 'Enrichment Score of EMT geneset')


es11_gsva <- as.matrix(es11_gsva)

library(RColorBrewer)


br61_common_genes_emt <- expr61[rownames(expr61) %in% brca_sets_up$BRCA_EMT,]

br61_common_genes_emt <- br61_common_genes_emt[rowSums(br61_common_genes_emt) > 0,]

br61_common_genes_invasion <- expr61[rownames(expr61) %in% brca_sets_up$BRCA_invasion,]

br61_common_genes_invasion <- br61_common_genes_invasion[rowSums(br61_common_genes_invasion) > 0,]

br61_common_genes_metastasis <- expr61[rownames(expr61) %in% brca_sets_up$BRCA_metastasis,]

br61_common_genes_metastasis <- br61_common_genes_metastasis[rowSums(br61_common_genes_metastasis) > 0,]

br61_common_genes_proliferation <- expr61[rownames(expr61) %in% brca_sets_up$BRCA_proliferation,]

br61_common_genes_proliferation <- br61_common_genes_proliferation[rowSums(br61_common_genes_proliferation) > 0,]

br61_common_genes_stemness <- expr61[rownames(expr61) %in% brca_sets_up$BRCA_stemness,]

br61_common_genes_stemness <- br61_common_genes_stemness[rowSums(br61_common_genes_stemness) > 0,]


br11_common_genes_emt <- expr11[rownames(expr11) %in% brca_sets_up$BRCA_EMT,]

br11_common_genes_emt <- br11_common_genes_emt[rowSums(br11_common_genes_emt) > 0,]

br11_common_genes_invasion <- expr11[rownames(expr11) %in% brca_sets_up$BRCA_invasion,]

br11_common_genes_invasion <- br11_common_genes_invasion[rowSums(br11_common_genes_invasion) > 0,]

br11_common_genes_metastasis <- expr11[rownames(expr11) %in% brca_sets_up$BRCA_metastasis,]

br11_common_genes_metastasis <- br11_common_genes_metastasis[rowSums(br11_common_genes_metastasis) > 0,]

br11_common_genes_proliferation <- expr11[rownames(expr11) %in% brca_sets_up$BRCA_proliferation,]

br11_common_genes_proliferation <- br11_common_genes_proliferation[rowSums(br11_common_genes_proliferation) > 0,]

br11_common_genes_stemness <- expr11[rownames(expr11) %in% brca_sets_up$BRCA_stemness,]

br11_common_genes_stemness <- br11_common_genes_stemness[rowSums(br11_common_genes_stemness) > 0,]


