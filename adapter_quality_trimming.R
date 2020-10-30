##### library 
library(seqRFLP)

##### parameters

# cutadapt location
cutadaptPath <- '/home/abartho/.local/bin/cutadapt'

# output folder
outPath <- 'out'


# universal adapters
adapterUniF_R1 <- 'GTGCCAGCMGCCGCGGTAA' # 5'->3'
adapterUniF_R2 <- 'GGACTACHVGGGTWTCTAAT' # 5'->3'
adapterUniR_R1 <- revComp(adapterUniF_R2)
adapterUniR_R2 <- revComp(adapterUniF_R1)

# create output folder if not exists
outPath <- gsub("/$","",outPath) # remove tailing slash
outPath <- paste0(outPath,'/')
if(outPath != '' && !dir.exists(outPath)){
  dir.create(outPath)
}

# generate sample names
fnFs <- c()
fnRs <- c()
sample.names <- c()
path <- '/rawdata'
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs_buff <- sort(list.files(path=path, pattern="R1_001.fastq", full.names = TRUE))
fnRs_buff <- sort(list.files(path=path, pattern="R2_001.fastq", full.names = TRUE))
fnFs <- c(fnFs,fnFs_buff)
fnRs <- c(fnRs,fnRs_buff)
# Extract sample names
sample.names <- gsub("_L001_.+","",basename(fnFs_buff))

### status
print(paste0('  ',Sys.time()))
print('  Trimming adapters')

samples <- sample.names
### trim quality and universal adapter 
for(i in 1:length(fnFs)) {
  samp <- gsub("_L001_.+","",basename(fnFs[i]))
  system(paste0(
    cutadaptPath,' --cores=20 -e 0.2 -q 15,15 -m 150 --discard-untrimmed',
    ' -a "',adapterUniF_R1,'...',adapterUniR_R1,';e=0.2"',
    ' -A "',adapterUniF_R2,'...',adapterUniR_R2,';e=0.2"',
    ' -o ',outPath,'trimmedFinal_',samp,'_R1.fastq', # R1 output
    ' -p ',outPath,'trimmedFinal_',samp,'_R2.fastq ', # R2 output
    fnFs[i],' ', # R1 input
    fnRs[i],' > ', # R2 input
    outPath,'outputTrimmedAll_',samp,'.txt' # store STDOUT in file
  ))
}
