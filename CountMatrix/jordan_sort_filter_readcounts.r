jor=read.delim('Jordan_readCounts.csv',sep='\t',header=T)
dim(jor)
colnames(jor)[1]='symbol'
head(jor)

s=sort(colnames(jor))

jor=subset(jor,select = c(s))
ncol(jor)
jor=jor[,c(75,1:74)]
jor


jor['total']=rowSums(jor[,2:75],dims=1)
jor=jor[order(-jor$total),]
jor

jorFilt=jor[!jor$total == 0, ]
dim(jorFilt)

write.table(jorFilt,file='Jordan_readCountsFilteredSorted.csv',sep=',',row.names = F,quote=F)
