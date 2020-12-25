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
