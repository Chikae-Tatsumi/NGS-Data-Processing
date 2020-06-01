#DADA2 ITS Pipeline Workflow (1.8)
#See https://benjjneb.github.io/dada2/ITS_workflow.html

#Getting Ready
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")

DATABASE = "~/R/Database/sh_general_release_dynamic_02.02.2019.fasta"
setwd("~/R/Analysis/1_Test/ITS")  ## CHANGE ME to the directory containing the fastq files.
list.files()

# fnFs <- sort(list.files(getwd(), pattern = "_R1_001.fastq", full.names = TRUE))
# fnRs <- sort(list.files(getwd(), pattern = "_R2_001.fastq", full.names = TRUE)
fnFs <- sort(list.files(getwd(), pattern = ".fastq", full.names = TRUE))

#Identify primers
# FWD <- "CTTGGTCATTTAGAGGAAGTAA"  ## ITS1f
# REV <- "GCTGCGTTCTTCATCGATGC"  ## ITS2

FWD <-	"TCCGTAGGTGAACCTGCGG" ## ITS1
REV <-	"GCTGCGTTCTTCATCGATGC" ## ITS2

#to ensure we have the right primers, and the correct orientation of the primers on the reads
allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
        RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name ))
head(sample.names)
fnFs.filtN <- file.path(getwd(), "filtN", paste0(sample.names, ".fastq.gz"))

#to “pre-filter” the sequences just to remove those with  ambiguous bases (Ns)
# fnFs.filtN <- file.path(getwd(), "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
# fnRs.filtN <- file.path(getwd(), "filtN", basename(fnRs))
# filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)
filterAndTrim(fnFs, fnFs.filtN, maxN = 0, multithread = TRUE)

# Identifying and counting the primers on one set of paired end FASTQ files(Sufficient, because all the files were created using the same library preparation)
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
  #  FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
  #  REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
  #  REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]) )

#Remove Primers
cutadapt <- "/miniconda2/bin/cutadapt" 
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut <- file.path(getwd(), "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
# fnRs.cut <- file.path(path.cut, basename(fnRs))
fnRs.cut <- file.path(path.cut, basename(fnFs)) # dummy file

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
# for(i in seq_along(fnFs)) {
#  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
#                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
#                             fnFs.filtN[i], fnRs.filtN[i])) # input files
#}

for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, "-n", 1,
                            "-o", fnFs.cut[i],  # output file                             
                            fnFs.filtN[i])) # input files
}


# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
#    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
#    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
#    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
# cutFs <- sort(list.files(path.cut, pattern = "_1.fastq.gz", full.names = TRUE))
# cutRs <- sort(list.files(path.cut, pattern = "_2.fastq.gz", full.names = TRUE))
cutFs <- sort(list.files(path.cut, pattern = ".fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

#Inspect read quality profiles
plotQualityProfile(cutFs[1:2])
# plotQualityProfile(cutRs[1:2])

#Filter and trim
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
# filtRs <- file.path(path.cut, "filtered", basename(cutRs))

# out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
#    truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE

out <- filterAndTrim(cutFs, filtFs, maxN = 0, maxEE = 2,
    truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

#Learn the Error Rates
errF <- learnErrors(filtFs, multithread = TRUE)
# errR <- learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)

#Dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)
# derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
# names(derepRs) <- sample.names

#Sample Inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
# dadaRs <- dada(derepRs, err = errR, multithread = TRUE)


#Merge paired reads
# mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#Construct Sequence Table
# seqtab <- makeSequenceTable(mergers)

#see https://github.com/benjjneb/dada2/issues/384
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
table(nchar(getSequences(seqtab.nochim)))


#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
# track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
#    getN), rowSums(seqtab.nochim))

track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
# colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
#     "nonchim")

colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")

rownames(track) <- sample.names
head(track)
write.table(track,file="track.txt")

#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, DATABASE, multithread = TRUE, tryRC = TRUE)
taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
write.table(taxa, file="taxonomy.txt")
write.table(seqtab.nochim, file="seqtabnochim.txt")

samples.out<-rownames(seqtab.nochim)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                         tax_table(taxa))
#store the DNA sequences of our ASVs in the refseq slot of the phyloseq object

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#To output OTU table
otu_table.t<-t(ps@otu_table)
ps.t<-cbind(otu_table.t,ps@tax_table)
write.table(ps.t,  file="ASV_table.txt")
