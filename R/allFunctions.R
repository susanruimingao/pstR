
trimReads = function(inputDir,
                     outputDir = "trimmedReads",
                     suffixNameR1 = "_1.fastq.gz", #must be these forms; _R1.fastq.gz, _1.fastq.gz, _R1.fq.gz, _1.fq.gz;,
                     suffixNameR2 = "_2.fastq.gz",
                     nThreads = 16, ...){
  #loading libraries;
  #if(!require("pacman")) install.packages("pacman");
  #pacman::p_load(parallel, crayon, magrittr, tidyverse);
  library(parallel); library(crayon); library(stringr);

  #check number of threads;
  nCores = detectCores();
  if(nCores < nThreads){
    nThreads = nCores;
    cat(format(Sys.time(), usetz = TRUE), yellow(" using ", nThreads), green("threads"), "\n");
  }

  #create output file;
  if(!dir.exists(outputDir)){
    cat(format(Sys.time(), usetz = TRUE), yellow(" creating the outputfile ", outputDir), "\n");
    dir.create(outputDir);
  } else {
    cat(format(Sys.time(), usetz = TRUE), red(" outputfile ", outputDir, " exits, is going to remove and create a new one"), "\n");
    unlink(outputDir, recursive = TRUE);
    dir.create(outputDir);
  }

  #get file infor
  fileNames = list.files(inputDir);
  fileNamesLeft = fileNames[str_detect(fileNames, suffixNameR1)] %>% str_remove(suffixNameR1);
  fileNamesRight = fileNames[str_detect(fileNames, suffixNameR2)] %>% str_remove(suffixNameR2);
  nSamples = 0;

  if(all(fileNamesLeft == fileNamesRight)){
    nSamples = length(fileNamesLeft);
    cat(format(Sys.time(), usetz = TRUE), yellow(" all read files (", nSamples, ") are paired\n"));
    fileLeftInput = file.path(inputDir, paste0(fileNamesLeft, suffixNameR1));
    fileRightInput = file.path(inputDir, paste0(fileNamesLeft, suffixNameR2));
    fileMegeredOutput = file.path(outputDir, paste0(fileNamesLeft, "_clean_merged.fastq.gz"));
    fileLeftOutput = file.path(outputDir, paste0(fileNamesLeft, "_clean_unmerged_R1.fastq.gz"));
    fileRightOutput = file.path(outputDir, paste0(fileNamesRight, "_clean_unmerged_R2.fastq.gz"));
    fileHtml = file.path(outputDir, paste0(fileNamesLeft, ".html"));
  } else {
    cat(format(Sys.time(), usetz = TRUE));
    stop(red(" pairend reads are not paired, please check and ensure pairedend reads"), "\n");
  }

  #using fastp to trim;
  for(i in 1:nSamples){
    cmd = paste("fastp",
                "-i", fileLeftInput[i],
                "-I", fileRightInput[i],
                "--correction",
                "--merge",
                "--merged_out", fileMegeredOutput[i],
                "-o", fileLeftOutput[i],
                "-O", fileRightOutput[i],
                "-w", nThreads,
                "--html", fileHtml[i],
                collapse = "\t");
    cat(format(Sys.time(), usetz = TRUE), yellow(" begin to process quality check using fastp for sample ", fileNamesLeft[i]), " ", green(i, " out of ", nSamples, "\n"));
    system(cmd);
    cat(format(Sys.time(), usetz = TRUE), yellow(" finish to process quality check using fastp for sample ", fileNamesLeft[i]), " ", green(i, " out of ", nSamples, "\n"));
  }
}



assemblySpades = function(trimmedReadsDir = "trimmedReads",
                          outputDir = "spadesOut",
                          kmer = 21, #kmer size must be 21, 33, 55, 77, 99, 127, odd numbers;
                          minLenth = 2000, #filter the minimum contig length;
                          nMemory = 30, #memory in GB;
                          nThreads = 16, ...){
  #loading libraries;
  #if(!require("pacman")) install.packages("pacman");
  #pacman::p_load(parallel, crayon, magrittr, tidyverse, seqinr);
  library(parallel); library(crayon); library(stringr); library(seqinr);

  #check number of threads;
  nCores = detectCores();
  if(nCores < nThreads){
    nThreads = nCores;
    cat(format(Sys.time(), usetz = TRUE), yellow(" using ", nThreads), green("threads"), "\n");
  }

  #create output file;
  if(!dir.exists(outputDir)){
    cat(format(Sys.time(), usetz = TRUE), yellow(" creating the outputfile ", outputDir), "\n");
    dir.create(outputDir);
  } else {
    cat(format(Sys.time(), usetz = TRUE), red(" outputfile ", outputDir, " exits, is going to remove and create a new one"), "\n");
    unlink(outputDir, recursive = TRUE);
    dir.create(outputDir);
  }

  #get file infor
  fileNames = list.files(trimmedReadsDir, pattern = "_clean_merged.fastq.gz") %>%
    str_remove("_clean_merged.fastq.gz");

  nSamples = length(fileNames);

  if(nSamples < 1){
    cat(format(Sys.time(), usetz = TRUE));
    stop(red(" detect ", nSamples, " , please check your ", trimmedReadsDir), "\n");
  } else {
    cat(format(Sys.time(), usetz = TRUE), yellow(" detect", nSamples, ") for spades assemble\n"));
    fileMergedInput = file.path(trimmedReadsDir, paste0(fileNames, "_clean_merged.fastq.gz"));
    fileLeftUnMerged = file.path(trimmedReadsDir, paste0(fileNames, "_clean_unmerged_R1.fastq.gz"));
    fileRightUnMerged = file.path(trimmedReadsDir, paste0(fileNames, "_clean_unmerged_R2.fastq.gz"));
  }

  for(i in 1 : nSamples){
    cmd = paste("spades.py",
                "--threads", nThreads,
                "--memory", nMemory,
                "-k", paste(kmer, collapse = ","),
                "--careful",
                "--s1", fileMergedInput[i],
                "--pe1-1", fileLeftUnMerged[i],
                "--pe1-2", fileRightUnMerged[i],
                "-o", file.path(outputDir, fileNames[i]));
    print(cmd);
    cat(format(Sys.time(), usetz = TRUE), yellow(" begin to spades assemble for sample ", fileNames[i]), " ", green(i, " out of ", nSamples, "\n"));
    system(cmd);

    for(j in list.files(file.path(outputDir, fileNames[i]))){
      if(j != "scaffolds.fasta"){
        j = file.path(outputDir, fileNames[i], j);
        if(file.exists(j) && !dir.exists(j)){
          file.remove(j);
        } else {
          unlink(j, recursive = TRUE);
        }
      }
    }

    fastqFile = read.fasta(file.path(outputDir, fileNames[i], "scaffolds.fasta"),
                           as.string = TRUE, seqtype = "DNA");

    fastqFile[getLength(fastqFile) >= minLenth] %>%
      write.fasta(sequences = .,
                  names = names(.),
                  as.string = TRUE,
                  file.out = file.path(outputDir, fileNames[i], paste0(fileNames[i], "_assembly_filter", minLenth, ".fasta")));

    file.rename(file.path(outputDir, fileNames[i], "scaffolds.fasta"),
                file.path(outputDir, fileNames[i], paste0(fileNames[i], "_assembly.fasta")));

    cat(format(Sys.time(), usetz = TRUE), yellow(" finish to spades assemble for sample ", fileNames[i]), " ", green(i, " out of ", nSamples, "\n"));
  }
}



submit_to_PHASTER = function(inputDir = "spadesOut",
                             suffix = "_assembly_filter2000.fasta",
                             outputDir = "PHASTEROut"
                             , ...){
  #loading libraries;
  #if(!require("pacman")) install.packages("pacman");
  #pacman::p_load(parallel, crayon, magrittr, tidyverse, seqinr);
  library(crayon);

  #create output file;
  if(!dir.exists(outputDir)){
    cat(format(Sys.time(), usetz = TRUE), yellow(" creating the outputfile ", outputDir), "\n");
    dir.create(outputDir);
  } else {
    cat(format(Sys.time(), usetz = TRUE), red(" outputfile ", outputDir, " exits, is going to remove and create a new one"), "\n");
    unlink(outputDir, recursive = TRUE);
    dir.create(outputDir);
  }

  fileNames = list.files(inputDir);

  for(i in 1: length(fileNames)){
    fileName = fileNames[i];
    fastaFileName = paste0(fileName, "_assembly_filter2000.fasta");
    inputFile = file.path(inputDir, fileName, fastaFileName);
    print(inputFile);
    outputFile = file.path(outputDir, fileName);
    print(outputFile);
    cmd = paste0('wget --post-file="',
                 inputFile,
                 '" "https://phaster.ca/phaster_api?contigs=1" ',
                 '-O ', paste0(outputFile, ".json"),
                 ' -o ', paste0(outputFile, "_wget.log"));

    print(cmd);
    cat(format(Sys.time(), usetz = TRUE), yellow(" begin to submit file to online server ", fileName), " ", green(i, " out of ", length(fileNames), "\n"));
    system(cmd);
  }

}



CheckPhasterServer = function (inputDir = "spadesOut",
                               outputDir = "PHASTEROut",
                               ...){
  library(crayon);
  if(!dir.exists(inputDir) || length(list.files(inputDir)) == 0){
    stop("please check the ", inputDir);
  }

  if(!dir.exists(outputDir) || length(list.files(outputDir)) == 0){
    stop("please check the ", outputDir);
  }

  fileNames = list.files(inputDir);

  for(i in 1: length(fileNames)){
    fileName = fileNames[i];
    fastaFileName = paste0(fileName, "_assembly_filter2000.fasta");
    inputFile = file.path(inputDir, fileName, fastaFileName);
    print(inputFile);
    outputFile = file.path(outputDir, fileName);
    print(outputFile);

    cmd = paste0("python3 /home/CFIA-ACIA/gaoru/bin/phage_typing/checkPhasterServer.py ",
                 " --check ",
                 " -i ", inputDir,
                 " -o ", outputDir);
    print(cmd);
    cat(format(Sys.time(), usetz = TRUE), yellow(" begin to check the checkPhasterServer.py"));
    system(cmd);
  }
}


extract_fasta = function (inputDir = "PHASTEROut",
                          outputDir = "extractFasta"
                          , ...){
  #loading libraries;
  #if(!require("pacman")) install.packages("pacman");
  #pacman::p_load(crayon, magrittr, tidyverse, seqinr);
  library(crayon);  library(stringr); library(seqinr); library(magrittr);

  #create output file;
  if(!dir.exists(outputDir)){
    cat(format(Sys.time(), usetz = TRUE), yellow(" creating the outputfile ", outputDir), "\n");
    dir.create(outputDir);
  } else {
    cat(format(Sys.time(), usetz = TRUE), red(" outputfile ", outputDir, " exits, is going to remove and create a new one"), "\n");
    unlink(outputDir, recursive = TRUE);
    dir.create(outputDir);
  }

  fileNames = list.files(inputDir, pattern = "_phaster.zip") %>%
    str_remove(., "_phaster.zip");

  write.table(fileNames,
              file.path(outputDir, "sampleList.txt"),
              sep = "\t",
              row.names = F,
              quote = F);


  fastaList <- list();
  for(i in 1:length(fileNames)){
    fileName = fileNames[i];
    cat(format(Sys.time(), usetz = TRUE), yellow(" begin to unzip file ", fileName), " ", green(i, " out of ", length(fileNames), "\n"));
    unzip(file.path(inputDir, paste0(fileName, "_phaster.zip")),
          exdir = file.path(inputDir, paste0(fileName, "_phaster")));

    cat(format(Sys.time(), usetz = TRUE), yellow(" begin to read fasta file ", fileName), " ", green(i, " out of ", length(fileNames), "\n"));
    tmp = read.fasta(file.path(inputDir, paste0(fileName, "_phaster"),  "phage_regions.fna"),
                     as.string = TRUE, seqtype = "DNA");
    names(tmp) <- paste(fileName, names(tmp), sep = "_");
    cat(format(Sys.time(), usetz = TRUE), yellow(" begin to write fasta file ", fileName), " ", green(i, " out of ", length(fileNames), "\n"));
    write.fasta(sequences = tmp,
                names=names(tmp),
                as.string = TRUE,
                file.out = file.path(outputDir, paste0(fileName, ".fasta")));
    fastaList <- c(fastaList, tmp);
  }

  cat(format(Sys.time(), usetz = TRUE), yellow(" begin to write all_phage.fasta"));
  write.fasta(sequences = fastaList,
              names=str_remove(names(fastaList), ".*\\."),
              as.string = TRUE,
              file.out = file.path(outputDir, "all_phage.fasta"));

}


cluster_sequences = function (inputFile = "all_phage.fasta",
                              outputDir = "clusterSeqs",
                              c = 0.9,
                              s = 0.9,
                              n = 10,
                              nThreads = 8,
                              mem = 1000,
                              d = 0,
                              ...){

  #loading libraries;
  #if(!require("pacman")) install.packages("pacman");
  #pacman::p_load(parallel, crayon, magrittr, tidyverse, seqinr);
  library(crayon);

  #create output file;
  if(!dir.exists(outputDir)){
    cat(format(Sys.time(), usetz = TRUE), yellow(" creating the outputfile ", outputDir), "\n");
    dir.create(outputDir);
  } else {
    cat(format(Sys.time(), usetz = TRUE), red(" outputfile ", outputDir, " exits, is going to remove and create a new one"), "\n");
    unlink(outputDir, recursive = TRUE);
    dir.create(outputDir);
  }


  cmd = paste0("cd-hit-est ",
               " -i ", inputFile,
               " -c ", c,
               " -s ", s,
               " -n ", n,
               " -o ", file.path(outputDir, paste0("phage_clustered_c", c, "_s", s, ".fasta")),
               " -T ", nThreads,
               " -M ", mem,
               " -d ", d);

  print(cmd);
  cat(format(Sys.time(), usetz = TRUE), yellow(" begin to cdhit"));
  system(cmd);
  cat(format(Sys.time(), usetz = TRUE), yellow(" finish the cdhit"));

}




# If qiime is installed under the phage_typing environment, ignore this comment.
# In my case, As the qiime requires Python2.7, thus create an individual new environment
# Under /home/CFIA-ACIA/gaoru/miniconda3/envs/qiime1 environment to call "PST_pipelineFunctions.R"


create_biom_table = function (inputFile = "phage_clustered.fasta.clstr",
                              outputFile = "phage_clustered.fasta.tsv",
                              sampleList = "sampleList.txt",
                              ...){


  #loading libraries;
  #if(!require("pacman")) install.packages("pacman");
  #pacman::p_load(parallel, crayon, magrittr, tidyverse, seqinr);
  library(crayon);
  library(seqinr);

  cmd = paste0("perl /home/CFIA-ACIA/gaoru/bin/phage_typing/cdHitClstr2table.pl ",
               " -s ", sampleList,
               " -i ", inputFile,
               " -o ", outputFile);
  print(cmd);
  cat(format(Sys.time(), usetz = TRUE), yellow(" begin to create biom table using qiime1"));
  system(cmd);
}



# To make run_qiime funtion (three parts)

biom_convert = function (inputFile = "phage_clustered.fasta.tsv",
                         outputFile = "phage_clustered.fasta.biom",
                         ...){

  #loading libraries;
  #if(!require("pacman")) install.packages("pacman");
  #pacman::p_load(parallel, crayon, magrittr, tidyverse, seqinr);
  library(crayon);

  cmd = paste0("biom convert ",
               " -i ", inputFile,
               " -o ", outputFile,
               " --table-type='OTU table'",
               " --to-json");
  print(cmd);
  cat(format(Sys.time(), usetz = TRUE), yellow(" begin to convert biom"));
  system(cmd);
  cat(format(Sys.time(), usetz = TRUE), yellow(" finish convert biom"));
}



beta_diversity = function (inputFile = "phage_clustered.fasta.biom",
                           outputDir = "beta_div_nonPhylo",
                           ...){

  #loading libraries;
  #if(!require("pacman")) install.packages("pacman");
  #pacman::p_load(parallel, crayon, magrittr, tidyverse, seqinr);
  library(crayon);

  #create output file;
  if(!dir.exists(outputDir)){
    cat(format(Sys.time(), usetz = TRUE), yellow(" creating the outputfile ", outputDir), "\n");
    dir.create(outputDir);
  } else {
    cat(format(Sys.time(), usetz = TRUE), red(" outputfile ", outputDir, " exits, is going to remove and create a new one"), "\n");
    unlink(outputDir, recursive = TRUE);
    dir.create(outputDir);
  }


  cmd = paste0("beta_diversity.py ",
               " -i ", inputFile,
               " -o ", outputDir,
               " -m euclidean");


  print(cmd);
  cat(format(Sys.time(), usetz = TRUE), yellow(" begin for performing beta-diversity"));
  system(cmd);
  cat(format(Sys.time(), usetz = TRUE), yellow(" finish beta-diversity"));

}



neighbor_joining = function (inputFile = "euclidean_phage_clustered.txt",
                             outputFile = "phage_clustered.tree",
                             ...){

  #loading libraries;
  #if(!require("pacman")) install.packages("pacman");
  #pacman::p_load(parallel, crayon, magrittr, tidyverse, seqinr);
  library(crayon);


  cmd = paste0("neighbor_joining.py ",
               " -i ", inputFile,
               " -o ", outputFile);


  print(cmd);
  cat(format(Sys.time(), usetz = TRUE), yellow(" begin to build tree"));
  system(cmd);
  cat(format(Sys.time(), usetz = TRUE), yellow(" finish neighbor_joining"));

}





