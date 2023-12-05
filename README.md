# DADA2 Walkthrough #
In this repository, you will find a basic intro to the Divisive Amplicon Denoising Algorithm 2 (DADA2) pipeline based .

### Below you will see the basic functionalities of each command ###

## filterAndTrim() ##
The `filterAndTrim()` function allows you to clean up your reads before putting it through the dada algorthim. You have can clean it up using the parameters of:
1. `truncLen` trims your sequence to the desired length of your amplicon.
1. `maxN` how many ambiguous base pairs (N) you will allow in your reads; if 0, then reads with "N"s will be removed.
1. `maxEE` how low the expected error can be; if the read has a lower expected error, it is removed.
1. `truncQ` how low the Phred score for any single base call can be; if a base call has a lower score, the read is removed.
1. `rm.phix` if true, removes sequences that match the PhiX genome.

```
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)
```
## learnErrors() ##
The `learnErrors()` function takes your reads and learns the rate of error for each nucleotide comparison for all forward and reverse reads.
```
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)
```
## dada() ##
The `dada()` function constructs a dada-class object by denoising your sequences using the dada algorthim using the errors specific to your data.
```
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errF, multithread=FALSE)
```
## mergePairs() ##
After constructing dada-class objects from your forward and reverse reads, you can recombine the two to form one dada-class object 
```
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```
## makeSequenceTable() ##
Converts dada-class object to digestible form of an ASV table
```
seqtab <- makeSequenceTable(mergers)
```
## removeBimeraDenovo() ##
Finds and removes sequences that may have been produced through chimera sequencing (two or more templates used to produce one read)
```
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```
