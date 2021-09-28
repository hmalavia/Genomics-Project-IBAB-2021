BiocManager::install("DEsingle")
library('DEsingle')
colData(sce.hvg)
group <- factor(sce.hvg$CL_or_SC)

results <- DEsingle(counts=sce.hvg, group = group)
results.classified <- DEtype(results = results,threshold = 0.05)
results.sig <-results.classified[results.classified$pvalue.adj.FDR<0.05,]

results.sig
results.classified
results.DEs <- results.sig[results.sig$Type == "DEs", ]
results.DEa <- results.sig[results.sig$Type == "DEa", ]
results.DEg <- results.sig[results.sig$Type == "DEg", ]

results.DEg.up=results.DEg[results.DEg$State == 'up',]
results.DEa.up=results.DEa[results.DEa$State == 'up',]

results.DEup=rbind.data.frame(results.DEg.up,results.DEa.up)

results.DEup$log2FC=log2(results.DEup$norm_foldChange)

results.DEup =subset(results.DEup,select=c('norm_foldChange','log2FC','Type'))

results.DEup=results.DEup[order(-results.DEup$log2FC),]
write.csv(results.DEup,file='Aceto_CL_VS_SC_DE_UP_DEsignal.csv',row.names = T,quote=F)

results.DEg.down=results.DEg[results.DEg$State == 'down',]
results.DEa.down=results.DEa[results.DEa$State == 'down',]

results.DEdown=rbind.data.frame(results.DEg.down,results.DEa.down)

results.DEdown$log2FC=log2(results.DEdown$norm_foldChange)

results.DEdown =subset(results.DEdown,select=c('norm_foldChange','log2FC','Type'))

results.DEdown=results.DEdown[order(-results.DEdown$log2FC),]
write.csv(results.DEdown,file='Aceto_CL_VS_SC_DE_Down_DEsignal.csv',row.names = T,quote=F)



results.DEup

write.csv(results.DEa,file='aceto_results.DEa.csv',row.names = T,quote = F)
write.csv(results.DEg,file='aceto_results.DEg.csv',row.names = T,quote = F)



epigenes=read.csv('Epigenes_unique.csv')

results.EpiDE.UP=merge(epigenes,results.DEup,by.x='HGNC.approved.symbol',by.y=0)
results.EpiDE.down=merge(epigenes,results.DEdown,by.x='HGNC.approved.symbol',by.y=0)

ano.EpiDE.UP=merge(ano,results.EpiDE.UP)
ano.EpiDE.UP=ano.EpiDE.UP[order(-ano.EpiDE.UP$log2FC),]
write.csv(ano.EpiDE.UP,file='Aceto_CL_Vs_SC_EpiDE_Up_DESignal.csv',row.names = F,quote = F)
ano.EpiDE.down=merge(ano,results.EpiDE.down)
ano.EpiDE.down=ano.EpiDE.down[order(-ano.EpiDE.down$log2FC),]
write.csv(ano.EpiDE.down,file='Aceto_CL_Vs_SC_EpiDE_Down_DESignal.csv',row.names = F,quote = F)
