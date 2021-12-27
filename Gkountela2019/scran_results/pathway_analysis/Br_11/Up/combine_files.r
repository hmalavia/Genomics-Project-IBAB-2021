file1 <- read.csv('Br11_up_ARCHS4_Kinases_Coexp_table.txt',sep='\t')
file2 <- read.csv('Br11_up_BioCarta_2016_table.txt',sep='\t')

load_data <- function(path) { 
  files <- dir(path, pattern = '*.txt', full.names = T)
  tables <- lapply(files, read.csv,sep='\t')
  do.call(rbind, tables)
  }

Br11_up <- load_data("~/Documents/Project/dataset/Gkountela2019/scran_results/pathway_analysis/Br_11/Up")

Br11_up <- Br11_up[order(Br11_up$Adjusted.P.value,decreasing = F),]

Br11_up_sig <- Br11_up[Br11_up$Adjusted.P.value < 0.05,]

write.csv(Br11_up,'Br11_up_allPathways',row.names = F,quote = F)

write.csv(Br11_up_sig,'Br11_up_SigPathways',row.names = F,quote = F)

Br11_down <- load_data("~/Documents/Project/dataset/Gkountela2019/scran_results/pathway_analysis/Br_11/down/")

Br11_down <- Br11_down[order(Br11_down$Adjusted.P.value,decreasing = F),]

Br11_down_sig <- Br11_down[Br11_down$Adjusted.P.value < 0.05,]

write.csv(Br11_down,'Br11_down_allPathways',row.names = F,quote = F)

write.csv(Br11_down_sig,'Br11_down_SigPathways',row.names = F,quote = F)

Br61_up <- load_data("~/Documents/Project/dataset/Gkountela2019/scran_results/pathway_analysis2/Br_61/up/")

Br61_up <- Br61_up[order(Br61_up$Adjusted.P.value,decreasing = F),]

Br61_up_sig <- Br61_up[Br61_up$Adjusted.P.value < 0.05,]

write.csv(Br61_up,'../../../pathway_analysis2/Br_61/up/Br61_up_allPathways',row.names = F,quote = F)

write.csv(Br61_up_sig,'../../../pathway_analysis2/Br_61/up/Br61_up_SigPathways',row.names = F,quote = F)

Br61_down <- load_data("~/Documents/Project/dataset/Gkountela2019/scran_results/pathway_analysis2/Br_61/down/")

Br61_down <- Br61_down[order(Br61_down$Adjusted.P.value,decreasing = F),]

Br61_down_sig <- Br61_down[Br61_down$Adjusted.P.value < 0.05,]

write.csv(Br61_down,'../../../pathway_analysis2/Br_61/down/Br61_down_allPathways',row.names = F,quote = F)

write.csv(Br61_down_sig,'../../../pathway_analysis2/Br_61/down/Br61_down_SigPathways',row.names = F,quote = F)
