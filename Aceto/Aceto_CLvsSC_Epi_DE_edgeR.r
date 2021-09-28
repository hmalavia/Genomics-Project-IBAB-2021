allde=read.csv('aceto_CLvsSC_DGE_edgeR.csv')
epifactors=read.csv('Epigenes_unique.csv')
colnames(allde)
epigenesDE=merge(epifactors,allde,by.x='HGNC.approved.symbol',by.y='X')
epigenesDE=epigenesDE[order(-epigenesDE$logFC),]
upregepi=epigenesDE[epigenesDE$logFC > 2,]
dim(upregepi)
signiup=upregepi[upregepi$PValue < 0.05,]

dim(signiup)
upregepi
downregepi=epigenesDE[epigenesDE$logFC < -2,]
dim(downregepi)
signidown=downregepi[downregepi$PValue<0.05,]
dim(signidown)
dim(downregepi)
downregepi

alldeUP=allde[allde$logFC>2,]
dim(alldeUP)

alldeDown=allde[allde$logFC < -2,]
dim(alldeDown)

allsigniUP=alldeUP[alldeUP$PValue < 0.05,]
dim(allsigniUP)

allsignidown=alldeDown[alldeDown$PValue < 0.05,]
dim(allsignidown)

write.csv(upregepi,file='Aceto_CL_vs_SC_EpiDE_AllUP_edgeR.csv',row.names=F,quote = F)
write.csv(downregepi,file='Aceto_CL_vs_SC_EpiDE_AllDown_edgeR.csv',row.names=F,quote = F)
write.csv(signiup[1:10,],file='Aceto_CL_vs_SC_EpiDE_Top10UP_edgeR.csv',row.names=F,quote = F)
write.csv(signidown[1:10,],file='Aceto_CL_vs_SC_EpiDE_Top10Down_edgeR.csv',row.names=F,quote = F)

ano=read.csv('DGE/Epifactorgenes.csv')
signiup=read.csv('Aceto_CL_vs_SC_EpiDE_Top10UP_edgeR.csv')

signiupano=merge(ano,signiup,by.x='HGNC.approved.symbol',by.y='HGNC.approved.symbol')
signiupano=signiupano[order(-signiupano$logFC),]
write.csv(signiupano,file = 'Aceto_CL_vs_SC_EpiDE_Top10UP_edgeR.csv',row.names=F,quote = F )
signidown=read.csv('Aceto_CL_vs_SC_EpiDE_Top10Down_edgeR.csv')
signidownano=merge(ano,signidown,by.x='HGNC.approved.symbol',by.y='HGNC.approved.symbol')
signidownano=signidownano[order(-signidownano$logFC),]
write.csv(signidownano[1:10,],file = 'Aceto_CL_vs_SC_EpiDE_Top10Down_edgeR.csv',row.names=F,quote = F )
