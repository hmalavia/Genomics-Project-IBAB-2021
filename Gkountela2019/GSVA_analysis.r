# BiocManager::install("GSVA")
library(GSVA)
library(stringr)

file_names <- list.files('CancerSEA_markers/',full.names = T)

file_names

set_names <- tools::file_path_sans_ext(file_names)

csea_sets <- lapply(file_names,read.table,header=F) 

names(csea_sets) <- set_names

names(csea_sets)[[14]] <- "Stemness"

names(csea_sets)

csea_sets <- lapply(csea_sets, function(x) x[!names(x) %in% c("V1")])

csea_sets1 <- lapply(csea_sets, function(x) as.list(x))

csea_sets1 <- lapply(csea_sets, function(x) x[x != "GeneName"])


expr11 <- read.csv('br11_scran_norm_counts.csv')

colnames(expr11)[1] <- "GeneName"

row.names(expr11) <- expr11$GeneName

expr11$GeneName <- NULL

expr11 <- as.matrix(expr11)

es11_2 <- gsva(expr11,csea_sets1,kcdf="Gaussian",method='ssgsea')

sdata <- read.csv('br11_filtered_coldata.csv')

colnames(es11_2) <- sdata$Sample_ID











