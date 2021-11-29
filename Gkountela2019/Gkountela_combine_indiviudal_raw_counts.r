#install.packages('plyr')
library('plyr')
library('stringr')


path='Individual_Rawcounts/'

load_data <- function(path) { 
  files <- dir(path, pattern = '*.txt', full.names = TRUE)
  tables <- lapply(files, read.csv,sep="\t")
  do.call(cbind,tables)
}

rawmat <- load_data(path)

rawmat <- rawmat[,unique(colnames(rawmat))]

colnames(rawmat) <- str_replace_all(colnames(rawmat),'[_]+','.')

head(rawmat)

write.csv(rawmat,'Gkountela_Patient_unprocessed_rawcounts.csv',quote = F,row.names = F)
