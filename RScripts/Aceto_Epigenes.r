counts=read.csv('Aceto_anotated_readcounts.csv')
dim(counts)


epigenes_main=read.csv('Epifactorgenes.csv')

epigenes=epigenes_main['HGNC.approved.symbol']
head(epigenes)

epigenes_uniq=unique(epigenes)

write.csv(epigenes_uniq,file='Epigenes_unique.csv',row.names = F,quote = F)

aceto_epigenes=merge(epigenes_uniq,counts,by.x = 'HGNC.approved.symbol',by.y = 'symbol')

aceto_epigenes_sorted=aceto_epigenes[aceto_epigenes$total>10,]

aceto_epigenes_sorted=aceto_epigenes_sorted[order(-aceto_epigenes_sorted$total),]
write.csv(aceto_epigenes_sorted,file = 'EpigeneCounts_Aceto_FIlteredSorted.csv',row.names = F,quote = F)
