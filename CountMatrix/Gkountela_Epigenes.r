gCounts=read.delim('Gkountela_processed_normalized_matrix.csv',sep='\t')

ncol(gCounts)

gCounts['total']=rowSums(gCounts[,2:211],dims = 1)

gCounts_filtered=gCounts[gCounts$total>0,]

gCounts_filtered=gCounts_filtered[order(-gCounts_filtered$total),]

write.csv(gCounts_filtered,file = 'Gkountela_counts_filtered',row.names = F,quote = F)

epigenes_gkountela=merge(epigenes_uniq,gCounts_filtered,by.x='HGNC.approved.symbol',by.y = 'Geneid')

gkountela_epigenes_filtered=epigenes_gkountela[epigenes_gkountela$total>10,]

gkountela_epigenes_filtered=gkountela_epigenes_filtered[order(-gkountela_epigenes_filtered$total),]

write.csv(gkountela_epigenes_filtered,file = 'Gkountela_epigenesCount_sorted.csv',row.names = F,quote = F)
