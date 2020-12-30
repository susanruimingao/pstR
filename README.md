# pstR
Prophage Sequence Typing writen in R language; 

pstR version='0.3.0'

This pstR pipeline works in R

The pipleline is a continuous effort from previous PST bash pipeline: 
author='duceppemo' version='0.2.0'

In a terminal, conda create a "phage_typing" environment 
```
conda install -c bioconda fastp

conda install -c bioconda spades

conda install -c bioconda cd-hit
```

#If the python version is incompatible, the qiime environment can be created separately

conda install -c qiime2 qiime2

The scripts "checkPhasterServer.py" and "cdHitClstr2table.pl" which are writen by Marco, should be available in the path. This will be added as a parameter in the respective step.  

Download prophageTypingR_0.1.0.tar.gz from
https://github.com/susanruimingao/pstR.git


Go to work directory containing raw reads, and the reads name is preferred in "_R1.fastq.gz"/"_R2.fastq.gz"; otherwise specify in the command line

Conda activate phage_typing

Go to R environment


install.packages("../prophageTypingR_0.1.0.tar.gz", repos = NULL, type="source")

#To run the pipeline steps until submitting to PHASTER, using general function: including trim reads, spades assembly and submit the assemblies to PHASTER server

prophageTypingR::trimAssembleSubmit(inputDir = "rawdata",  suffixNameR1 = "_R1.fastq.gz",suffixNameR2 = "_R2.fastq.gz" );


#After submitting the assemlies to PHASTER, To check the status of PHASTER server running:

prophageTypingR::CheckPhasterServer(path = "YOUR_OWN_PATH_to../checkPhasterServer.py")


(Option 1) #After downloading zip files from PHASTER 
(including retrieve phage .fasta sequences, clustering and create phylogenetic tree with qiime)

prophageTypingR::ClusterRunQiime()


(Option 2) if the qiime environment is installed individually, the cd-hit clustering step is separated with the rest running qiime
under phage_typing environment

prophageTypingR::extract_fasta()

prophageTypingR::cluster_sequences(inputFile = "./extractFasta/all_phage.fasta", c = 0.99, s = 0.99, outputDir = "./extractFasta/clusterSeqs_99_99", path = "YOUR_OWN_PATH_to../cdHitClstr2table.pl")


Under qiime environment:

prophageTypingR::runQiime(inputFile = "./extractFasta/clusterSeqs_99_99/phage_clustered_c0.99_s0.99.fasta.clstr", sampleList = "./extractFasta/sampleList.txt", path = "YOUR_OWN_PATH_to../cdHitClstr2table.pl")

#The finall created .tree file is the desired output



###Notes: the individual functions are also available if you want to run each of them separately:

prophageTypingR::trimReads(inputDir = "rawdata",  suffixNameR1 = "_R1.fastq.gz",suffixNameR2 = "_R2.fastq.gz" );

prophageTypingR::assemblySpades()

prophageTypingR::submit_to_PHASTER()

prophageTypingR::CheckPhasterServer(path = "YOUR_OWN_PATH_to../checkPhasterServer.py")

prophageTypingR::extract_fasta()

prophageTypingR::cluster_sequences(inputFile = "./extractFasta/all_phage.fasta", c = 0.99, s = 0.99, outputDir = "./extractFasta/clusterSeqs_99_99")

prophageTypingR::create_biom_table(path = "YOUR_OWN_PATH_to../cdHitClstr2table.pl")

prophageTypingR::biom_convert()

prophageTypingR::beta_diversity()

prophageTypingR::neighbor_joining()
