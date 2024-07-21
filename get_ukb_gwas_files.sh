#!/bin/sh

while read f
    do 
    
    # print some thing to the screen
    echo "Processing data for phenotype $f."
    
    # GET the file from AWS 
	wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/$f

	# rename the file
	mv $f tempFile.tsv.gz
	gunzip tempFile.tsv.gz
	
	echo "Now cleaning the data"

	# SELECTION ONLY A COLUMNS
	mv tempFile.tsv $f.txt
    
done < /required_files/ukb_epi_gwas_filenames.txt