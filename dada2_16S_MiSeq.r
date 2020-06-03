#DADA2 Pipeline Tutorial (1.12)
#See http://benjjneb.github.io/dada2/tutorial.html

#Getting ready
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")

DATABASE = "~/R/Database/tax/silva_nr_v132_train_set.fa.gz"
setwd("~/R/Analysis/1_Test/16S")  ## CHANGE ME to the directory containing the fastq files.
list.files()
# a set of Illumina-sequenced paired-end fastq files that have been split by sample and from which the barcodes/adapters have already been removed

# Define which is forward fastq and reverse fastq
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(getwd(), pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(getwd(), pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Inspect read quality profiles
# visualizing the quality profiles of the forward reads and reverse reads (only sample 1 and 2):
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
# You should judge how length is good quality, and truncate the bad quality part

# Filter and trim
filtFs <- file.path(getwd(), "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(getwd(), "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
# If too few reads are passing the filter, consider relaxing maxEE, perhaps especially on the reverse reads (eg. maxEE=c(2,5)), and reducing the truncLen to remove low quality tails
# Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE) #plot

# Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Merge paired reads
 mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) 
table(nchar(getSequences(seqtab))) # Inspect distribution of sequence lengths

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) #>90% is good?

# Track reads through the pipeline
# look at the number of reads that made it through each step in the pipeline:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track,file="track.txt")

#Assign taxonomy
#Install "silva_nr_v132_train_set.fa.gz" from: https://zenodo.org/record/1172783#.XUmvQ_ZFw2w
#Other taxonomic reference database: http://benjjneb.github.io/dada2/training.html
taxa <- assignTaxonomy(seqtab.nochim,DATABASE, multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
#If your reads do not seem to be appropriately assigned, for example lots of your bacterial 16S sequences are being assigned as Eukaryota NA NA NA NA NA, your reads may be in the opposite orientation as the reference database. Tell dada2 to try the reverse-complement orientation with assignTaxonomy(..., tryRC=TRUE) and see if this fixes the assignments. 

write.table(taxa,file="taxonomy.txt")
write.table(seqtab.nochim,file="seqtabnochim.txt")

taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
write.table(taxa, file="taxonomy.txt")
write.table(seqtab.nochim, file="seqtabnochim.txt")

samples.out<-rownames(seqtab.nochim)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                         tax_table(taxa))
#store the DNA sequences of our ASVs in the refseq slot of the phyloseq object

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#To output OTU table
otu_table.t<-t(ps@otu_table)
ps.t<-cbind(otu_table.t,ps@tax_table)
write.table(ps.t,  file="ASV_table.txt")