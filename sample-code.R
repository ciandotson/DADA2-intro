# General pipeline to run the Divisive Amplicon Denoising Algorithm 2 (dada2) with your microbiome data #  

## Installing dada2 package from benjjneb/dada2 ##
install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16") # change the ref argument to get other versions
library(dada2); packageVersion("dada2")

## Making a connection to the fastq files that have had their primers trimmed ##
### This can be performed using a series of trimming software, including trimmomatic, or can be done on the commandline using grep commands ### 

path <- "~/DADA2-intro/MiSeq_SOP"
list.files(path)

## Separate forward and reverse reads #
fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2_001.fastq", full.names = TRUE))

## Visualizing the general patterns of read quality through graphing average Phred score at each base pair ##
plotQualityProfile(fnFs[1:16])
plotQualityProfile(fnRs[1:16])

## Assigning each individual submission by its sample name found within the actual file name ##
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

## Copying file paths containing the reads into a subdirectory and adding "filt" to the end of each file name ##
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

## Actually filtering the reads within the files denoted by the file paths to trim forward reads to 250 bps and reverse reads to 160 (as inferred from the quality profiles of each sample), no Ns in the reads, removing forward and reverse reads with expected error greater than or equal to 2, truncating sequences to only include bases that have a Q score greater than 2, removing reads that match the PhiX genome, and performing everything in parallel ##
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)
head(out)

## Learning error rates that is specific to YOUR data ## 
### This will be important when the actual algorithm develops an error model when performing sequence comparisons
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)

# Plotting errors compared to the expected values generated a priori; outputs a 16 x N matrix of transistion probabilities from nucletiode base to another when considering quality scores #
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

## Constructing dada-class objects using the actual algorithm ##
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)
dadaFs[[1]]

## merging denoised forward and reversed reads ##
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

## Constructing a sequence table from the merged reads that will become the ASV table ##
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

## removing chimeras, or sequences that are exact copies of two abundant sequences spliced together ##
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

## determine ratio of non-chimeras to all reads to determine how usable your data is post all sequencing steps ##
sum(seqtab.nochim)/sum(seqtab)

## track reads through the pipeline to see if any steps filtered out a substantial number of reads ##
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
