# General pipeline to run the Divisive Amplicon Denoising Algorithm 2 (dada2) with your microbiome data #  

## Installing dada2 package from benjjneb/dada2 ##
install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16") # change the ref argument to get other versions
library(dada2); packageVersion("dada2")


path <- "~/Sequencing_Data/16SLegumeReads/trimmed_reads"
list.files(path)

# separate forward and reverse reads #
fnFs <- sort(list.files(path, pattern="_R1_001.trimmed.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.trimmed.fastq", full.names = TRUE))

# Assigning each individual submission by its sample name #
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# placing filtered files in subdirectory #
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filtering #
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

# Learning error rates # 
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Plotting Errors #
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Constructing dada-class objects #
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

# merging denoised forward and reversed reads #
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Constructing Sequence Table #
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# removing chimeras #
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# determine ratio of non-chimeras to all reads #
sum(seqtab.nochim)/sum(seqtab)

# track reads through the pipeline #
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
