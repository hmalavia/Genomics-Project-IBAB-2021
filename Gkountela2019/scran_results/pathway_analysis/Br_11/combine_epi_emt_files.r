file1 <- read.csv('Br11_up_ARCHS4_Kinases_Coexp_table.txt',sep='\t')
file2 <- read.csv('Br11_up_BioCarta_2016_table.txt',sep='\t')

load_data <- function(path) { 
  files <- dir(path, pattern = '*.txt', full.names = T)
  tables <- lapply(files, read.csv,sep='\t')
  do.call(rbind, tables)
}

#### Br11 ####

#### Epigenes ####

Br11_epi_up <- load_data("~/Documents/Project/dataset/Gkountela2019/scran_results/pathway_analysis/Br_11/Epigenes_up_Br11/")

Br11_epi_up <- Br11_epi_up[order(Br11_epi_up$Adjusted.P.value,decreasing = F),]

Br11_epi_up_sig <- Br11_epi_up[Br11_epi_up$Adjusted.P.value < 0.05,]

write.csv(Br11_epi_up,'Br11_epi_up_allPathways',row.names = F,quote = F)

write.csv(Br11_epi_up_sig,'Br11_epi_up_SigPathways',row.names = F,quote = F)

Br11_epi_down <- load_data("~/Documents/Project/dataset/Gkountela2019/scran_results/pathway_analysis/Br_11/Epigenes_Down_Br11/")

Br11_epi_down <- Br11_epi_down[order(Br11_epi_down$Adjusted.P.value,decreasing = F),]

Br11_epi_down_sig <- Br11_epi_down[Br11_epi_down$Adjusted.P.value < 0.05,]

write.csv(Br11_epi_down,'Br11_epi_down_allPathways',row.names = F,quote = F)

write.csv(Br11_epi_down_sig,'Br11_epi_down_SigPathways',row.names = F,quote = F)

#### EMT genes ####

Br11_emt_up <- load_data("~/Documents/Project/dataset/Gkountela2019/scran_results/pathway_analysis/Br_11/EMT_up_Br11/")

Br11_emt_up <- Br11_emt_up[order(Br11_emt_up$Adjusted.P.value,decreasing = F),]

Br11_emt_up_sig <- Br11_emt_up[Br11_emt_up$Adjusted.P.value < 0.05,]

write.csv(Br11_emt_up,'Br11_emt_up_allPathways',row.names = F,quote = F)

write.csv(Br11_emt_up_sig,'Br11_emt_up_SigPathways',row.names = F,quote = F)

Br11_emt_down <- load_data("~/Documents/Project/dataset/Gkountela2019/scran_results/pathway_analysis/Br_11/EMT_down_Br11/")

Br11_emt_down <- Br11_emt_down[order(Br11_emt_down$Adjusted.P.value,decreasing = F),]

Br11_emt_down_sig <- Br11_emt_down[Br11_emt_down$Adjusted.P.value < 0.05,]

write.csv(Br11_emt_down,'Br11_emt_down_allPathways',row.names = F,quote = F)

write.csv(Br11_emt_down_sig,'Br11_emt_down_SigPathways',row.names = F,quote = F)
###### Br61 #######

#### Epi genes ####

Br61_epi_up <- load_data("~/Documents/Project/dataset/Gkountela2019/scran_results/pathway_analysis2/Br_61/Epigenes_up_Br61/")

Br61_epi_up <- Br61_epi_up[order(Br61_epi_up$Adjusted.P.value,decreasing = F),]

Br61_epi_up_sig <- Br61_epi_up[Br61_epi_up$Adjusted.P.value < 0.05,]

write.csv(Br61_epi_up,'../../../pathway_analysis2/Br_61/Epigenes_up_Br61/Br61_epi_up_allPathways',row.names = F,quote = F)

write.csv(Br61_epi_up_sig,'../../../pathway_analysis2/Br_61/Epigenes_up_Br61/Br61_epi_up_SigPathways',row.names = F,quote = F)

#### EMT genes ####
Br61_emt_up <- load_data("~/Documents/Project/dataset/Gkountela2019/scran_results/pathway_analysis2/Br_61/EMT_up_Br61/")

Br61_emt_up <- Br61_emt_up[order(Br61_emt_up$Adjusted.P.value,decreasing = F),]

Br61_emt_up_sig <- Br61_emt_up[Br61_emt_up$Adjusted.P.value < 0.05,]

write.csv(Br61_emt_up,'../../../pathway_analysis2/Br_61/EMT_up_Br61/Br61_emt_up_AllPathways',row.names = F,quote = F)

write.csv(Br61_emt_up_sig,'../../../pathway_analysis2/Br_61/EMT_up_Br61/Br61_emt_up_SigPathways',row.names = F,quote = F)

Br61_emt_down <- load_data("~/Documents/Project/dataset/Gkountela2019/scran_results/pathway_analysis2/Br_61/EMT_down_Br61/")

Br61_emt_down <- Br61_emt_down[order(Br61_emt_down$Adjusted.P.value,decreasing = F),]

Br61_emt_down_sig <- Br61_emt_down[Br61_emt_down$Adjusted.P.value < 0.05,]

write.csv(Br61_emt_down,'../../../pathway_analysis2/Br_61/EMT_down_Br61/Br61_emt_down_allPathways',row.names = F,quote = F)

write.csv(Br61_emt_down_sig,'../../../pathway_analysis2/Br_61/EMT_down_Br61/Br61_emt_down_SigPathways',row.names = F,quote = F)


#### Find Common Pathways ####

Common_emt_up <- merge(Br11_emt_up_sig,Br61_emt_up_sig,by.x=1,by.y = 1)

Common_emt_down <- merge(Br11_emt_down_sig,Br61_emt_down_sig,by.x=1,by.y=1)

Common_epi_up <- merge(Br11_epi_up_sig,Br61_epi_up_sig,by.x = 1,by.y=1)

Common_emt_up <- Common_emt_up[order(Common_emt_up$Combined.Score.x,decreasing = T),]

Common_emt_down <- Common_emt_down[order(Common_emt_down$Combined.Score.x,decreasing = T),]

Common_epi_up <- Common_epi_up[order(Common_epi_up$Combined.Score.x,decreasing = T),]


write.csv(Common_emt_up,'../Common_emt_up_pathways.csv',row.names = F,quote = F)
write.csv(Common_emt_down,'../Common_emt_down_pathways.csv',row.names = F,quote = F)
write.csv(Common_epi_up,'../Common_epi_up_pathways.csv',row.names = F,quote = F)
