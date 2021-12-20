## Finding Common Genes ##

#### CTC DEG list ####

br11UP <- read.csv('scran_results/Br11/Br11_upregulated_genes_ALL_edgeR.csv')
br11Down <- read.csv('scran_results/Br11_downregulated_genes_ALL_edgeR.csv')

br16Up <- read.csv('scran_results/Br16_upregulated_genes_ALL_edgeR.csv')
br16Down <- read.csv('scran_results/Br16_downregulated_genes_ALL_edgeR.csv')

br7Up <- read.csv('scran_results/Br7_upregulated_genes_ALL_edgeR.csv')
br7Down <- read.csv('scran_results/Br7_downregulated_genes_ALL_edgeR.csv')

br61Up <- read.csv('scran_results/Br61_upregulated_genes_ALL_edgeR.csv')
br61Down <- read.csv('scran_results/Br61_downregulated_genes_ALL_edgeR.csv')

#### Load Epigenes and EMT genes list ####

epigenes <- read.csv('../CountMatrix/Epigenes_unique.csv')

emt <- read.csv('EMT_genes_unique.csv')

#### Finding DE Epigenes ####

FindDEepi <- function(br,){
   epi<- merge(epigenes,br,by.x = 1,by.y = 1)
   return(epi)
   
}

#br11EpiUp <- merge(epigenes,br11UP,by.x = 1,by.y = 1)

br11EpiUp <- FindDEepi(br11UP)
br16EpiUp <- FindDEepi(br16Up)
br61EpiUp <- FindDEepi(br61Up)
br7EpiUp <- FindDEepi(br7Up)

br11EpiDown <- FindDEepi(br11Down)
br16EpiDown <- FindDEepi(br16Down)
br61EpiDown <- FindDEepi(br61Down)
br7EpiDown <- FindDEepi(br7Down)

#### Writing Output ####
writeCsv <- function(br){
  df <- deparse(substitute(br))
  write.csv(br,paste0("scran_results/",df,".csv",sep=""),row.names = F,quote = F)
}
 
writeCsv(br11EpiDown) 
writeCsv(br11EpiUp) 
writeCsv(br7EpiDown) 
writeCsv(br7EpiUp) 
writeCsv(br16EpiDown) 
writeCsv(br16EpiUp) 
writeCsv(br61EpiDown) 
writeCsv(br61EpiUp) 

#### Finding DE EMT genes ####
writeEMT <- function(br){
  df <- deparse(substitute(br))
  write.csv(br,paste0("scran_results/",df,".csv",sep=""),row.names = F,quote = F)
}

FindDEemt <- function(br){
  emt<- merge(emt,br,by.x = 1,by.y = 1)
  return(emt)
}

br11EMTUp <- FindDEemt(br11UP)
br16EMTUp <- FindDEemt(br16Up)
br61EMTUp <- FindDEemt(br61Up)
br7EMTUp <- FindDEemt(br7Up)

br11EMTDown <- FindDEemt(br11Down)
br16EMTDown <- FindDEemt(br16Down)
br61EMTDown <- FindDEemt(br61Down)
br7EMTDown <- FindDEemt(br7Down)

writeCsv(br11EMTDown) 
writeCsv(br11EMTUp) 
writeCsv(br7EMTDown) 
writeCsv(br7EMTUp) 
writeCsv(br16EMTDown) 
writeCsv(br16EMTUp) 
writeCsv(br61EMTDown) 
writeCsv(br61EMTUp) 

#### Venn Diagrams ####
#install.packages('ggvenn')              

library('ggvenn')       

#### Venn diagrams for Upregulated Epigenes ####

