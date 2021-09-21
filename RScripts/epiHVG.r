hvg=read.csv('DGE/Gokentela_HVG.csv')
dim(hvg[hvg$bio >0,])
epigenes=read.csv('Epigenes_unique.csv')

epi_hvg=merge(epigenes,hvg,by.x='HGNC.approved.symbol',by.y='X')
epi_hvg=epi_hvg[order(epi_hvg$bio,decreasing = T),]
epi_hvg=epi_hvg[epi_hvg$bio >0,]
dim(epi_hvg)
write.csv(epi_hvg,file='gkountela_epiHVG.csv',row.names = F,quote=F)

top10=epi_hvg[1:10,]

write.csv(top10,file='jordan_top10_epiHVG.csv',row.names = F,quote = F)

####Anotation####

ano=read.csv('DGE/Epifactorgenes.csv')

Anotop10=merge(top10,ano,by.x='HGNC.approved.symbol',by.y='HGNC.approved.symbol')
Anotop10=Anotop10[order(-Anotop10$bio),]
write.csv(Anotop10,file='gkountela_Top10_HVEpigenes_anotated.csv',row.names = F,quote = F)
