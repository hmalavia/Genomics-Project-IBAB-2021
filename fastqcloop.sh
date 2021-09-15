dir=patient_*

for subdir in $dir
do
	cd $subdir
	mkdir QCFiles
	fastqc *.fastq -o QCFiles/
	cd ../
done 	
