#### Molecular subtype assignment using genefu ####
# BiocManager::install("genefu")
# 
# BiocManager::install("ensembldb")
# 
# BiocManager::install("EnsDb.Hsapiens.v86")

library('genefu')
library('ensembldb')
library('EnsDb.Hsapiens.v86')

data(scmod2.robust)
data(pam50.robust)
data(scmgene.robust)
data(sig.ggi)
data(scmod1.robust)
data(sig.genius)
data("ssp2006.robust")
#### Loading count Data #### 
count_data <- read.csv('br61_scran_norm_counts.csv',row.names = 'X')

sum <- rowSums(count_data)

count_data$sum <- sum

my.symbols <- rownames(count_data)

edb <- EnsDb.Hsapiens.v86

entrez_ids <- select(edb,  
                     keys = my.symbols, 
                     columns = c("ENTREZID", "SYMBOL"), 
                     keytype = "SYMBOL")

count_data <- merge(count_data,entrez_ids,by.x=0,by.y='SYMBOL')

count_data <- count_data[order(-count_data$sum),]

count_data <- count_data[!duplicated(count_data$Row.names),]

annot <- subset(count_data,select=c('Row.names','ENTREZID'))

colnames(annot)[2] <- 'EntrezGene.ID'

colnames(annot)[1] <- 'Gene.Symbol'

rownames(annot) <- annot$Gene.Symbol

rownames(count_data) <- count_data$Row.names

count_data <- count_data[,2:(ncol(count_data)-2)]

count_data <- t(count_data)

subtypePredictions_scmod1 <- molecular.subtyping(sbt.model = 'scmod1',data = count_data,annot = annot,do.mapping = T,verbose = T)

subtypePredictions_scmod2 <- molecular.subtyping(sbt.model = 'scmod2',data = count_data,annot = annot,do.mapping = T,verbose = T)

subtypePredictions_pam50 <- molecular.subtyping(sbt.model = '',data = count_data,annot = annot,do.mapping = T,verbose = T)


subtypePredictions_scmod1$subtype

Basals<-names(which(subtypePredictions_scmod1$subtype == "ER-/HER2-"))
#Select samples pertaining to HER2 Subtype
HER2s<-names(which(subtypePredictions_scmod1$subtype == "HER2+"))
#Select samples pertaining to Luminal Subtypes
LuminalB<-names(which(subtypePredictions_scmod1$subtype == "ER+/HER2- High Prolif"))
LuminalA<-names(which(subtypePredictions_scmod1$subtype == "ER+/HER2- Low Prolif"))

Basals<-names(which(subtypePredictions_scmod2$subtype == "ER-/HER2-"))
#Select samples pertaining to HER2 Subtype
HER2s<-names(which(subtypePredictions_scmod2$subtype == "HER2+"))
#Select samples pertaining to Luminal Subtypes
LuminalB<-names(which(subtypePredictions_scmod2$subtype == "ER+/HER2- High Prolif"))
LuminalA<-names(which(subtypePredictions_scmod2$subtype == "ER+/HER2- Low Prolif"))


sampledata <- read.csv('br61_scran_filtered_coldata.csv')

rownames(sampledata) <- sampledata$X

sampledata <- sampledata[,2:ncol(sampledata)]

sampledata$scmod1 <- ""

sampledata[Basals,]$scmod1 <- 'basal'
sampledata[HER2s,]$scmod1 <- 'her2'
sampledata[LuminalA,]$scmod1 <- 'lumA'
sampledata[LuminalB,]$scmod1 <- 'lumB'

sampledata$scmod2 <- ""

sampledata[Basals,]$scmod2 <- 'basal'
sampledata[HER2s,]$scmod2 <- 'her2'
sampledata[LuminalA,]$scmod2 <- 'lumA'
sampledata[LuminalB,]$scmod2 <- 'lumB'


write.csv(sampledata,'br61_molecular_subtypes_genefu.csv',quote = F)
