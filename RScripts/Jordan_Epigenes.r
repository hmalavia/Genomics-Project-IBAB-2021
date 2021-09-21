jCounts=read.csv('Jordan_readCountsFilteredSorted.csv')

jordan_epigenes=merge(epigenes_uniq,jCounts,by.x='HGNC.approved.symbol',by.y = 'symbol')

jordan_epigenes_filt=jordan_epigenes[jordan_epigenes$total>10,]

jordan_epigenes_filt=jordan_epigenes_filt[order(-jordan_epigenes_filt$total),]

write.csv(jordan_epigenes_filt,file='Jordan_epigenesCounts_sorted.csv',row.names = F,quote = F)
