# pstR
Prophage Sequence Typing writen in R language
pstR version='0.3.0'
This pstR pipeline works in R language

The pipleline is a continuous effort from previous PST bash pipeline: 
author='duceppemo' version='0.2.0'

In a terminal, conda create a "phage_typing" environment 

conda install -c bioconda fastp

conda install -c bioconda spades

conda install -c bioconda cd-hit

#If the python version is incompatible, the qiime environment can be created separately
conda install -c qiime2 qiime2


Download prophageTypingR_0.1.0.tar.gz from
https://github.com/susanruimingao/pstR.git

Go to workdirectory containing raw reads # the reads name is preferred in "_R1.fastq.gz"/"_R2.fastq.gz"; otherwise specify in the command line
Conda activate phage_typing
Go to R environment


install.packages("../prophageTypingR_0.1.0.tar.gz", repos = NULL, type="source")

#To run the pipeline steps after submitting to PHASTER, using general function: including trim reads, spades assembly and submit the assemblies to PHASTER server
prophageTypingR::trimAssembleSubmit(inputDir = "rawdata",  suffixNameR1 = "_R1.fastq.gz",suffixNameR2 = "_R2.fastq.gz" );


To check the status of PHASTER server running:
prophageTypingR::CheckPhasterServer()


(Option 1) #To start after downloading zip files from PHASTER (retrieve phage .fasta sequences, clustering and create phylogenetic tree with qiime)
prophageTypingR::ClusterRunQiime()


(Option 2) if the qiime environment is installed individually, the cd-hit clustering step is separated with the rest running qiime
under phage_typing environment
prophageTypingR::extract_fasta()

prophageTypingR::cluster_sequences(inputFile = "./extractFasta/all_phage.fasta", c = 0.99, s = 0.99, outputDir = "./extractFasta/clusterSeqs_99_99")


under qiime environment:
prophageTypingR::runQiime(inputFile = "./extractFasta/clusterSeqs_99_99/phage_clustered_c0.99_s0.99.fasta.clstr", sampleList = "./extractFasta/sampleList.txt")

The finall created .tree file is the final output


Notes: the individual functions are also available if you want to run each of them separately:

prophageTypingR::trimReads(inputDir = "rawdata",  suffixNameR1 = "_R1.fastq.gz",suffixNameR2 = "_R2.fastq.gz" );

prophageTypingR::assemblySpades()

prophageTypingR::submit_to_PHASTER()

prophageTypingR::CheckPhasterServer()

prophageTypingR::extract_fasta()

prophageTypingR::cluster_sequences(inputFile = "./extractFasta/all_phage.fasta", c = 0.99, s = 0.99, outputDir = "./extractFasta/clusterSeqs_99_99")

prophageTypingR::create_biom_table()

prophageTypingR::biom_convert()

prophageTypingR::beta_diversity()

prophageTypingR::neighbor_joining()
