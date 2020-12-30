# pstR
Prophage Sequence Typing pipeline is developed in R language and is a R wrapper for multiple functions; 
This pstR is based on a previous PST bash pipeline: 
author='duceppemo' version='0.2.0'


In a terminal, conda create a "phage_typing" environment and install three software fastp(v0.20.1), spades(v3.13.2) and cd-hit(v4.8.1) 
```
codna create –n phage_typing –c bioconda fastp spades cd-hit qiime2

```

ALTERNATIVELY: If qiime requires a different version python, which is incompatible with "phage_typing" environment, the qiime environment can be created separately
```
conda create -n qiime

conda install -c qiime2 qiime2
```

Two more scripts named "checkPhasterServer.py" and "cdHitClstr2table.pl" need to be downloaded from the folder named "scripts". The PATH containing this two scripts need to be called in R envrironment as puting as parameters in two steps of this pipeline. 


Download prophageTypingR_0.1.0.tar.gz from

https://github.com/susanruimingao/pstR.git


Go to work directory containing raw reads, and the reads name is preferred in "_R1.fastq.gz"/"_R2.fastq.gz"; otherwise specify in the command line.

Go to R environment and install the above downloaded prophageTyping package and four more other packages

```
R
install.packages("PATH/prophageTypingR_0.1.0.tar.gz", repos = NULL, type="source")
install.package("parallel"); 
install.package("crayon"); 
install.package("stringr"); 
install.package("seqinr");

```
To run the pipeline steps until submitting to PHASTER, using general function: including trim reads, spades assembly and submit the assemblies to PHASTER server

```
prophageTypingR::trimAssembleSubmit(inputDir = "rawdata",  suffixNameR1 = "_R1.fastq.gz",suffixNameR2 = "_R2.fastq.gz" );
```

After submitting the assemlies to PHASTER, to check the status of PHASTER server running:
```
prophageTypingR::CheckPhasterServer(path = "YOUR_OWN_PATH_to../checkPhasterServer.py")
```

After downloading zip files from PHASTER server 

##(Option 1): qiime is installed in the "phage_typing" environment
One step of running to get the final phylogenetic tree file (including retrieve phage .fasta sequences, clustering and create phylogenetic tree with qiime)

```
prophageTypingR::ClusterRunQiime()
```

##(Option 2): if the qiime environment is installed individually.

under phage_typing environment
```
prophageTypingR::extract_fasta()

prophageTypingR::cluster_sequences(inputFile = "./extractFasta/all_phage.fasta", c = 0.99, s = 0.99, outputDir = "./extractFasta/clusterSeqs_99_99", path = "YOUR_OWN_PATH_to../cdHitClstr2table.pl")
```

Under qiime environment:
```
prophageTypingR::runQiime(inputFile = "./extractFasta/clusterSeqs_99_99/phage_clustered_c0.99_s0.99.fasta.clstr", sampleList = "./extractFasta/sampleList.txt", path = "YOUR_OWN_PATH_to../cdHitClstr2table.pl")
```

***The finall created .tree file is the desired output



OPTIONAL: the individual functions are also available if you want to run each of them separately:
```
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
```
