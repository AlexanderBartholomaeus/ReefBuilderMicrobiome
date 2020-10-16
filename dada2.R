# load lib
library(dada2)


###### dada2 standard pipeline -------------------------------------------------
# see https://benjjneb.github.io/dada2/tutorial.html

# path to fastq files
fnFs <- c()
fnRs <- c()
sample.names <- c()

path <- paste0('/home/abartho/shany/dada2/adapterCut/out/')  
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path=path, pattern="R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path=path, pattern="R2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- gsub("trimmedFinal_(.+)_R1.fastq","\\1",basename(fnFs))

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(paste0("filtered"), paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(paste0("filtered"), paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# filter and trim 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,150),
                     maxN=0, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE,minLen = 150) 
head(out)

# learn errors
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

# sample inference
derepFs <- derepFastq(filtFs)
ddsFs <- dada(derepFs, err=errF, multithread=TRUE, pool = TRUE)
derepRs <- derepFastq(filtRs)
ddsRs <- dada(derepRs, err=errR, multithread=TRUE, pool = TRUE)

# merge pairs
mergers <- mergePairs(ddsFs, derepFs, ddsRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# contruct sequence tab
seqtab <- makeSequenceTable(mergers)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
dim(seqtab)

# track reads
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(ddsFs, getN), sapply(ddsRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

## assign taxonomy
# SILVA 138
taxa <- assignTaxonomy(seqtab.nochim, "/opt/databases/DADA2/dada2_training_data/silva_nr_v138_train_set.fa.gz", multithread=TRUE, tryRC = TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
# SILVA 132
taxa2 <- assignTaxonomy(seqtab.nochim, "/opt/databases/DADA2/dada2_training_data/silva_nr_v132_train_set.fa.gz", multithread=TRUE, tryRC = TRUE)
taxa2.print <- taxa2 # Removing sequence rownames for display only
rownames(taxa2.print) <- NULL

# try to assign species (not in standard workflow)
spec <- assignSpecies(seqtab.nochim,"/opt/databases/DADA2/dada2_training_data/silva_species_assignment_v132.fa.gz",tryRC = TRUE, n = 10000)


###### data analysis -------------------------------------------------

### merge and write tables
write.table(seqtab.nochim,'table_ASV_noChimera_raw.csv',sep=';',quote=T,col.names = T,row.names = F)
write.table(t(seqtab),'table_ASV_raw.csv',sep=';',quote=T,col.names = F,row.names = F)
seqName <- cbind(colnames(seqtab.nochim),paste0('ASV_',1:ncol(seqtab.nochim)))
colnames(seqName) <- c('Sequence','Name')
seqTable <- t(seqtab.nochim)
rownames(seqTable) <- seqName[,2]
colnames(seqTable) <- sample.names
seqTable <- cbind(seqTable,taxa,taxa2)
write.table(seqName,paste0('table_sequences_asvNames.csv'),sep=';',quote=T,col.names = T,row.names = F)
write.table(seqTable,paste0('table_ASVs.csv'),sep=';',quote=T,col.names = T,row.names = F)
write.table(track,paste0('table_stats.csv'),sep=';',quote=T,col.names = T,row.names = F)


# save workspace [optional]
#save.image()