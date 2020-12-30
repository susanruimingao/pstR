trimAssembleSubmit <- function(inputDir,
                               suffixNameR1 = "_1.fastq.gz", #must be these forms; _R1.fastq.gz, _1.fastq.gz, _R1.fq.gz, _1.fq.gz;,
                               suffixNameR2 = "_2.fastq.gz",
                               nThreads = 16,
                               kmer = 21, #kmer size must be 21, 33, 55, 77, 99, 127, odd numbers;
                               minLenth = 2000, #filter the minimum contig length;
                               nMemory = 30,
                               ...){
  trimReads(inputDir = inputDir,
            suffixNameR1 = suffixNameR1,
            suffixNameR2 = suffixNameR2,
            nThreads = nThreads);

  assemblySpades(trimmedReadsDir = "trimmedReads",
                            outputDir = "spadesOut",
                            kmer = kmer, #kmer size must be 21, 33, 55, 77, 99, 127, odd numbers;
                            minLenth = minLenth, #filter the minimum contig length;
                            nMemory = nMemory, #memory in GB;
                            nThreads = nThreads);

  submit_to_PHASTER(inputDir = "spadesOut",
                    suffix = "_assembly_filter2000.fasta",
                    outputDir = "PHASTEROut");
}


ClusterRunQiime <- function (inputDir = "PHASTEROut",
                             outputFile = "phage_clustered.tree",
                             c = 0.9,
                             s = 0.9,
                             n = 10,
                             nThreads = 8,
                             mem = 1000,
                             d = 0,
                     ...){
  
  extract_fasta(inputDir = inputDir,
                outputDir = "extractFasta");
  
  cluster_sequences(inputFile = "extractFasta/all_phage.fasta",
                    outputDir = "extractFasta/clusterSeqs",
                    c = c,
                    s = s,
                    n = n,
                    nThreads = nThreads,
                    mem = mem,
                    d = d);
  
  create_biom_table(inputFile = "extractFasta/clusterSeqs/phage_clustered.fasta.clstr",
                    outputFile = "extractFasta/clusterSeqs/phage_clustered.fasta.tsv",
                    sampleList = "extractFasta/sampleList.txt");
  
  biom_convert(inputFile = "extractFasta/clusterSeqs/phage_clustered.fasta.tsv",
               outputFile = "extractFasta/clusterSeqs/phage_clustered.fasta.biom");
  beta_diversity(inputFile = "extractFasta/clusterSeqs/phage_clustered.fasta.biom",
                 outputDir = "extractFasta/clusterSeqs/beta_div_nonPhylo");
  neighbor_joining(inputFile = "extractFasta/clusterSeqs/euclidean_phage_clustered.txt",
                   outputFile = outputFile);
}



runQiime <- function(inputFile = "phage_clustered.fasta.clstr",
                         outputFile = "phage_clustered.tree",
                      sampleList = "sampleList.txt",
                         ...){
  
  create_biom_table(inputFile = inputFile,
                    outputFile = file.path(dirname(inputFile), "phage_clustered.fasta.tsv"),
                    sampleList = sampleList);
  
  biom_convert(inputFile = file.path(dirname(inputFile), "phage_clustered.fasta.tsv"),
               outputFile = file.path(dirname(inputFile), "phage_clustered.fasta.biom"));
  
  beta_diversity(inputFile = file.path(dirname(inputFile), "phage_clustered.fasta.biom"),
                 outputDir = file.path(dirname(inputFile), "beta_div_nonPhylo"));

  
  neighbor_joining(inputFile = file.path(dirname(inputFile), "beta_div_nonPhylo/euclidean_phage_clustered.fasta.txt"),
                   outputFile = outputFile);
  
}
