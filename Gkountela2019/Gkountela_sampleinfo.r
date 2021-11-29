library('stringr')
patients=read.csv('Gkountela_patients.csv',header = F)

colnames(patients)=patients[1,]

colnames(patients)

patients <- patients[2:nrow(patients),]

grep("!Sample_characteristics_ch1",colnames(patients))

head(patients[,10])
colnames(patients)[10] <- "Origin"

head(patients[,11])
colnames(patients)[11] <- "Number_of_Cells"

head(patients[,12])
colnames(patients)[12] <- "Number_of_CTC"


head(patients[,14])
colnames(patients)[14] <- "Sample_Type"

head(patients[,15])
colnames(patients)[15] <- "Cell_Type"

head(patients[,16])
colnames(patients)[16] <- "CTC_Cluster_Size"

head(patients[,21])
colnames(patients)[21] <- "Donor"

sampleinfo=subset(patients,select=c('!Sample_title','Origin','Sample_Type','Number_of_Cells','CTC_Cluster_Size','Number_of_CTC','Donor'))
head(sampleinfo)

sampleinfo$`!Sample_title`=str_remove_all(sampleinfo$`!Sample_title`,pattern = ' #2')

head(sampleinfo)

sampleinfo$`!Sample_title`=str_replace_all(sampleinfo$`!Sample_title`,"[_]+",".")

head(sampleinfo)

sampleinfo$Origin=str_remove_all(sampleinfo$Origin,'origin: ')

head(sampleinfo)

sampleinfo$`Sample_Type`=str_remove_all(sampleinfo$`Sample_Type`,'sample type: ')

head(sampleinfo)

sampleinfo$`CTC_Cluster_Size`=str_remove_all(sampleinfo$`CTC_Cluster_Size`,'ctc-cluster size: ')

head(sampleinfo)

sampleinfo$Donor=str_remove_all(sampleinfo$Donor,'donor: ')

head(sampleinfo)

sampleinfo$`Number_of_CTC` <- str_remove_all(sampleinfo$`Number_of_CTC`,'number of ctc: ')

head(sampleinfo)

sampleinfo$`Number_of_Cells` <- str_remove_all(sampleinfo$`Number_of_Cells`,'number of cells: ')

head(sampleinfo)

write.csv(sampleinfo,'Gkountela_coldata.csv',row.names = F,quote = F)

patientcounts=read.csv('GkountelaPatientCountMatrix.csv')

colnames(patientcounts)[2:ncol(patientcounts)]

sampleinfo$`!Sample_title`

sampleinfo$`Number_of_Cells` <- as.integer(sampleinfo$`Number_of_Cells`)

sampleinfo$`CTC_Cluster_size` <- as.integer(sampleinfo$`CTC_Cluster_Size`)

head(sampleinfo)

sampleinfo$`Number_of_CTC` <- as.integer(sampleinfo$`Number_of_CTC`)

sampleinfo$`Sample_Type` <- factor(sampleinfo$`Sample_Type`)

levels(sampleinfo$Sample_Type)

str(sampleinfo$`Sample_Type`)
#install.packages('dplyr')
                 
library('dplyr')

glimpse(sampleinfo)

cluster_single_count <- sampleinfo %>%
  group_by(Sample_Type,Donor) %>%
  summarise(Samples=`!Sample_title`,count=n(),Nos_cells=Number_of_Cells)

write.csv(cluster_single_count,'Gkountela_CTC_Cluster_Single_counts.csv',quote = F,row.names = F)

sampleinfo

