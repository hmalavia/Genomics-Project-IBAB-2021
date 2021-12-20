br7all <- read.csv('Br7_edgeR_result.csv')

row.names(br7all) <- br7all$X 

br7sigUP <- br7all[br7all$PValue < 0.1 & br7all$logFC >= 1,]

br7sigDown <- br7all[br7all$PValue < 0.1 & br7all$logFC <= -1,]

br7sig <- rbind(br7sigUP,br7sigDown)

br7sig <- br7sig[order(-br7sig$logFC),]

write.csv(br7sig,'scran_results/Br7/br7_pval0.1_edgeR',quote = F)

br11all <- read.csv('Br11_edgeR_result.csv')

row.names(br11all) <- br11all$X 

br11sigUP <- br11all[br11all$PValue < 0.1 & br11all$logFC >= 1,]

br11sigDown <- br11all[br11all$PValue < 0.1 & br11all$logFC <= -1,]

br11sig <- rbind(br11sigUP,br11sigDown)

br11sig <- br11sig[order(-br11sig$logFC),]

write.csv(br11sig,'scran_results/Br11/br11_pval0.1_edgeR',quote = F)

br16all <- read.csv('Br16_edgeR_result.csv')

row.names(br16all) <- br16all$X 

br16sigUP <- br16all[br16all$PValue < 0.1 & br16all$logFC >= 1,]

br16sigDown <- br16all[br16all$PValue < 0.1 & br16all$logFC <= -1,]

br16sig <- rbind(br16sigUP,br16sigDown)

br16sig <- br16sig[order(-br16sig$logFC),]

write.csv(br16sig,'scran_results/Br16/br16_pval0.1_edgeR',quote = F)

br61all <- read.csv('Br61_edgeR_result.csv')

row.names(br61all) <- br61all$X 

br61sigUP <- br61all[br61all$PValue < 0.1 & br61all$logFC >= 1,]

br61sigDown <- br61all[br61all$PValue < 0.1 & br61all$logFC <= -1,]

br61sig <- rbind(br61sigUP,br61sigDown)

br61sig <- br61sig[order(-br61sig$logFC),]

write.csv(br61sig,'scran_results/Br61/br61_pval0.1_edgeR',quote = F)

epigenes <- read.csv('../CountMatrix/Epigenes_unique.csv')

emt <- read.csv('EMT_genes_unique.csv')

br7epi <- merge(br7sig,epigenes,by.x=1,by.y=1)

write.csv(br7epi,'scran_results/OnlyEpi/br7_epi_pval0.1.csv',row.names = F,quote = F)

br11epi <- merge(br11sig,epigenes,by.x=1,by.y=1)

write.csv(br11epi,'scran_results/OnlyEpi/br11_epi_pval0.1.csv',row.names = F,quote = F)

br16epi <- merge(br16sig,epigenes,by.x=1,by.y=1)

write.csv(br16epi,'scran_results/OnlyEpi/br16_epi_pval0.1.csv',row.names = F,quote = F)

br61epi <- merge(br61sig,epigenes,by.x=1,by.y=1)

write.csv(br61epi,'scran_results/OnlyEpi/br61_epi_pval0.1.csv',row.names = F,quote = F)

br7emt <- merge(br7sig,emt,by.x=1,by.y=1)

br7emt <- br7emt[order(-br7emt$logFC),]

write.csv(br7epi,'scran_results/OnlyEMT/br7_emt_pval0.1.csv',row.names = F,quote = F)

br11emt <- merge(br11sig,emt,by.x=1,by.y=1)

br11emt <- br11emt[order(-br11emt$logFC),]

write.csv(br11epi,'scran_results/OnlyEMT/br11_emt_pval0.1.csv',row.names = F,quote = F)

br61emt <- merge(br61sig,emt,by.x=1,by.y=1)

br61emt <- br61emt[order(-br61emt$logFC),]

write.csv(br61epi,'scran_results/OnlyEMT/br61_emt_pval0.1.csv',row.names = F,quote = F)

br16emt <- merge(br16sig,emt,by.x=1,by.y=1)

br16emt <- br16emt[order(-br16emt$logFC),]

write.csv(br16epi,'scran_results/OnlyEMT/br16_emt_pval0.1.csv',row.names = F,quote = F)
