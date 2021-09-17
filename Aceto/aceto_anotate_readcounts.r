ano=read.csv('Aceto_platform.csv')
head(ano)
counts=read.csv('Aceto_readCounts.csv')
head(counts)
colnames(counts)[1]='ID'
head(counts)


anofilt=subset(ano,select = c(ID,symbol))
dim(anofilt)

merged=merge(anofilt,counts)

anocounts=merged[,2:ncol(merged)]
head(anocounts)


s=sort(colnames(anocounts))


sortedanocounts=subset(anocounts,select = s) 

ncol(sortedanocounts)

sortedanocounts=sortedanocounts[,c(30,1:29)]
head(sortedanocounts)

total=rowSums(sortedanocounts[2:30],dims=1)

sortedanocounts['total']=total

head(sortedanocounts)

sortedanocounts=sortedanocounts[order(-total),]

dim(sortedanocounts[sortedanocounts$total == 0,])

sortedanocounts=sortedanocounts[!sortedanocounts$total == 0, ]

write.csv(sortedanocounts,file='Aceto_anotated_readcounts.csv',row.names =F ,quote = F)


