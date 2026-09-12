library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)

path <- "rawData/"
list.files(path)

fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))

plotQualityProfile(fnFs[1], aggregate = T)
plotQualityProfile(fnRs[1], aggregate = T)

#remove Ns
fnFs.filtN <- file.path(path, "filtN", basename(fnFs))
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

#primers
MIFISH_U_F <- "GTCGGTAAAACTCGTGCCAGC"
MIFISH_E_F <- "GTTGGTAAATCTCGTGCCAGC"

MIFISH_U_R <- "CATAGTGGGGTATCTAATCCCAGTTTG"
MIFISH_E_R <- "CATAGTGGGGTATCTAATCCTAGTTTG"

##############################################################################################################
allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))
}

primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(
  U_FWD.ForwardReads = sapply(allOrients(MIFISH_U_F), primerHits, fn = fnFs.filtN[[1]]),
  U_FWD.ReverseReads = sapply(allOrients(MIFISH_U_F), primerHits, fn = fnRs.filtN[[1]]),
  E_FWD.ForwardReads = sapply(allOrients(MIFISH_E_F), primerHits, fn = fnFs.filtN[[1]]),
  E_FWD.ReverseReads = sapply(allOrients(MIFISH_E_F), primerHits, fn = fnRs.filtN[[1]]),
  U_REV.ForwardReads = sapply(allOrients(MIFISH_U_R), primerHits, fn = fnFs.filtN[[1]]),
  U_REV.ReverseReads = sapply(allOrients(MIFISH_U_R), primerHits, fn = fnRs.filtN[[1]]),
  E_REV.ForwardReads = sapply(allOrients(MIFISH_E_R), primerHits, fn = fnFs.filtN[[1]]),
  E_REV.ReverseReads = sapply(allOrients(MIFISH_E_R), primerHits, fn = fnRs.filtN[[1]])
)
##############################################################################################################

cutadapt <- "/home/admuser/miniconda3/envs/cutadapt/bin/cutadapt"
system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

R1.flags <- c("-g", MIFISH_U_F, "-g", MIFISH_E_F)
R2.flags <- c("-G", MIFISH_U_R, "-G", MIFISH_E_R)

for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(
    R1.flags, R2.flags,
    "-n", 2, 
    "--discard-untrimmed", 
    "--minimum-length", 50,
    "--revcomp", 
    "-o", fnFs.cut[i],
    "-p", fnRs.cut[i],
    fnFs.filtN[i], fnRs.filtN[i] 
  ))
}

rbind(
  U_FWD.ForwardReads = sapply(allOrients(MIFISH_U_F), primerHits, fn = fnFs.cut[[1]]),
  U_FWD.ReverseReads = sapply(allOrients(MIFISH_U_F), primerHits, fn = fnRs.cut[[1]]),
  E_FWD.ForwardReads = sapply(allOrients(MIFISH_E_F), primerHits, fn = fnFs.cut[[1]]),
  E_FWD.ReverseReads = sapply(allOrients(MIFISH_E_F), primerHits, fn = fnRs.cut[[1]]),
  U_REV.ForwardReads = sapply(allOrients(MIFISH_U_R), primerHits, fn = fnFs.cut[[1]]),
  U_REV.ReverseReads = sapply(allOrients(MIFISH_U_R), primerHits, fn = fnRs.cut[[1]]),
  E_REV.ForwardReads = sapply(allOrients(MIFISH_E_R), primerHits, fn = fnFs.cut[[1]]),
  E_REV.ReverseReads = sapply(allOrients(MIFISH_E_R), primerHits, fn = fnRs.cut[[1]])
)


# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq.gz", full.names = TRUE))


get.sample.name <- function(fname) {
  # remove the extension part first
  fname <- basename(fname)
  fname <- sub("\\_L1_R1\\.fastq\\.gz$", "", fname)
  return(fname)
}

sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

#quality profiles of the primer removed reads
fq <- plotQualityProfile(cutFs, aggregate = T)
fq + geom_vline(xintercept = 110) + geom_hline(yintercept = 30)

rq <- plotQualityProfile(cutRs, aggregate = T)
rq + geom_vline(xintercept = 105) + geom_hline(yintercept = 30)

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,
                     truncLen=c(110, 105),
                     maxN=0,
                     maxEE=c(2,2),
                     truncQ=2,
                     rm.phix=TRUE,
                     compress=TRUE,
                     multithread=TRUE)
head(out)

plotQualityProfile(filtFs[1], aggregate = T)
plotQualityProfile(filtRs[1], aggregate = T)

#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool="pseudo")

#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap=12, verbose=TRUE)
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
track <- cbind(track, retained_pct = round(track[,"nonchim"] / track[,"input"] * 100, 1))

#export data
rownames(seqtab.nochim) <- sub("_R1\\.fastq\\.gz$", "", rownames(seqtab.nochim))

seqs <- colnames(seqtab.nochim)
otab <- t(seqtab.nochim)
asv_headers <- paste0("ASV_", seq(nrow(otab)))
rownames(otab) <- asv_headers

# Export FASTA
#writeLines(c(rbind(paste0(">", asv_headers), seqs)), "mifish_asvs.fa")

#export format
otab_export <- cbind("#OTU ID" = rownames(otab), otab)

# Export table
#write.table(otab_export, "mifish_asv_table.txt", quote=FALSE, sep="\t", row.names=FALSE)

#write.csv(track, "mifish_track_reads.csv", quote=FALSE)



