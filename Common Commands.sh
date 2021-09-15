##Prefetch##
nohup /opt/ohpc/pub/apps/sratoolkit.2.11.1-centos_linux64/bin/prefetch --option-file <Accesion list.txt>&

##Fasterq-dump##
nohup /opt/ohpc/pub/apps/sratoolkit.2.11.1-centos_linux64/bin/fasterq-dump -v -e 12 SRR* &

##fastqc##
nohup fastqc *.fastq -o QCFiles/
