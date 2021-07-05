# Distributed ChAMP pipeline - distChAMP

It`s simple tool to create myNorm file from *.idats using ChAMP package [1], but instead of process all files simultaneously it use distributed manner. This method is less computationally expensive and allows to process large amount of data in one run, however may consume much more time.

### How does it work

Sample sheet containing n-samples is shuffled and then split into k-batches containing ~ n/k elements. Then each batch is processed into temporary myNorm file using ChAMP pipeline. Next all temporary myNorms are merge into final object.


### How to run

To run **run_distributed_champ** function you need to install ChAMP package:

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    
    BiocManager::install("ChAMP")

then simply clone repository:

    git clone https://github.com/EpiGenMed/distChAMP
    
and use:

    source("path_to_file.R")

when function is correctly loaded type:
    
    run_distributed_champ(path_ss = path_to_sample_sheet_file, 
                          output = path_to_output_directory,
                          path_idats = path_to_idats_directory,
                          batch_size = number_of_samples_per_chunk)


**path_ss** -> path to sample sheet file [can not be in the same directory as idats].

**output** -> path to output directory.

**path_idats** -> path to directory containing *idats.

**batch_size** -> number of samples per batch.

##### ChAMP args
**array_type** -> EPIC / 450K.

**norm_type** -> BMIQ / SWAN / PBC / FunctionalNormliazation

**cores** -> number of cores to use [accelerate only normalization step]


Finally *output* directory contains myNorm file, ad QC/Norm files for each chunk. 


1. Tian Y, Morris TJ, Webster AP, Yang Z, Beck S, Andrew F, Teschendorff AE (2017). “ChAMP: updated methylation analysis pipeline for Illumina BeadChips.” Bioinformatics, btx513. doi: 10.1093/bioinformatics/btx513.


