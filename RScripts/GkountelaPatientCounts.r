library('stringr')
s_df=read.delim('GSE111065_series_matrix_sample.csv',sep = ' ',header = F )
s_df=as.data.frame(t(s_df))
colnames(s_df) <- s_df[1,]

countMatrix=read.csv('../CountMatrix/Gkountela_processed_normalized_matrix.csv',sep='\t')

s_df=s_df[2:nrow(s_df),]
s_df <- s_df[,s_df$`!Sample_characteristics_ch1 == origin: patient`]
patients <- s_df[s_df$`!Sample_characteristics_ch1` == 'origin: patient',]
patientSamples <- patients$`!Sample_title`

rownames(countMatrix) <- countMatrix$Geneid

colnames(countMatrix)

countMatrix <- countMatrix[,2:ncol(countMatrix)]

patientSamples
newSamples=str_replace(patientSamples,pattern = " #2",replacement = "")
newSamples

tmp <-intersect(colnames(countMatrix),newSamples)

tmp <-rbind(countMatrix$Geneid,tmp)

patientCounts <- countMatrix[,tmp]
rownames(patientCounts) <- row.names(countMatrix)
patientCounts


setdiff(newSamples,colnames(patientCounts))

colnames(countMatrix)
cols=str_replace_all(colnames(countMatrix),pattern = "[_]+",replacement = ".")
cols=str_replace_all(cols,pattern = "[././...]",replacement = ".")
cols
newSamples2=str_replace_all(newSamples,"[_]+",".")
newSamples2

tmp=intersect(cols,newSamples2)


cm2=read.csv('../CountMatrix/Gkountela_processed_normalized_matrix.csv',sep='\t')

row.names(cm2) <- cm2$Geneid
cm2=cm2[,2:ncol(cm2)]  
head(cm2)
colnames(cm2) <- str_replace_all(colnames(countMatrix),pattern = "[_]+",replacement = ".")

patientCounts <- cm2[,tmp]

write.csv(patientCounts,file='GkountelaPatientCountMatrix.csv',sep='\t',quote=F)
