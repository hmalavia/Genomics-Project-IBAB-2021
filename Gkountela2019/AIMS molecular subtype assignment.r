#### AIMS molecular subtype assignment ####

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("AIMS")
# BiocManager::install("EnsDb.Hsapiens.v86")
# BiocManager::install("ensembldb")

library('AIMS')
library('ensembldb')
library('EnsDb.Hsapiens.v86')

#### Loading count Data #### 

counts <- read.csv('Gkountela_Patient_unprocessed_rawcounts.csv')

rownames(counts) <- counts$Geneid

counts <- counts[,2:ncol(counts)]

counts <- as.matrix(counts)

br61 <- counts[,grep('Br61',colnames(counts))]
br61

outliers <- c('Br61.CTC.24') ## Remove outlier

br61_new <- br61[,!colnames(br61) %in% outliers]

br61_new$sum <- rowSums(br61_new) 

### Remove spikes ###
spikes <- grep('ERCC-',rownames(br61_new),value = T)  

br61_new <- br61_new[!rownames(br61_new) %in% spikes,]

br61_new <- br61_new[order(-br61_new$sum),]

uniq <- br61_new[!duplicated(br61_new),]

edb <-  EnsDb.Hsapiens.v86

columns(edb)

my.symbols <- rownames(uniq)

entrez_ids <- select(edb,  
       keys = my.symbols, 
       columns = c("ENTREZID", "SYMBOL"), 
       keytype = "SYMBOL")

uniq$entrez_ids <- entrez_ids$ENTREZID

uniq <- merge(uniq,entrez_ids,by.x=0,by.y='SYMBOL')

uniq$Row.names <- NULL
ncol(uniq)
br61_subtypes <-applyAIMS(as.matrix(uniq[,1:14]),uniq$ENTREZID) 

write.csv(br61_subtypes$cl,'br61_subtypes.csv',quote = F)
