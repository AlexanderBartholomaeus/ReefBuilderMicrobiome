# load libs
library(readxl)
library(writexl)

# read sample
samp <- read_excel('data/table_sample_stats.xlsx', skip = 1)
asv <- read.table('data/table_ASVs.csv', sep = ";", header = T, stringsAsFactors = F)
rownames(asv) <- paste0('ASV_',1:nrow(asv))

### parse ASV table
asv <- asv[,1:44] # remove SILVA 132
# modify names
samp$SampleName[match(colnames(asv)[1:38],samp$SampleID)]
colnames(asv)[1:38] <- samp$SampleName[match(colnames(asv)[1:38],samp$SampleID)]
# !!!exchange S_E3 and S_TP1a (wrong labeling in samples for sequencing)
buff <- asv$S_TP1a
asv$S_TP1a <- asv$S_E3
asv$S_E3 <- buff

# remove chloroplast and mitochondria
sum(asv$Order == 'Chloroplast',na.rm=T)
asv <- asv[!is.element(asv$Order,'Chloroplast'), ]
sum(asv$Family == 'Mitochondria',na.rm=T)
asv <- asv[!is.element(asv$Family,'Mitochondria'), ]


# calculate relative abundance
apply(asv[,1:38],2,sum)
for(i in 1:38){
  asv[,i] <- round(asv[,i]/sum(asv[,i]), digits=6) 
}

# merge sample names
agg <- aggregate(t(asv[,1:38]),by=list(nam=colnames(asv)[1:38]),FUN =mean)
agg <- t(agg)
colnames(agg) <- agg[1,]
agg <- agg[-1,]
agg <- data.frame(agg, stringsAsFactors = F)
for(i in 1:25){
 agg[,i] <- as.numeric(agg[,i])  
}
agg[,(ncol(agg)+1):(ncol(agg)+6)] <- asv[,39:44]

# remove rare-biosphere
for(i in 1:25){
  agg[agg[,i] < 0.001,i] <- 0
}
apply(agg[,1:25],2,sum)

# check and remove 0 rows
sum(apply(agg[,1:25],1,sum)==0)
agg <- agg[apply(agg[,1:25],1,sum)!=0,]
# do relative abundance again to come to 1 again
for(i in 1:25){
  agg[,i] <- round(agg[,i]/sum(agg[,i]), digits=6) 
}
agg_write2 <- cbind(rownames(agg),agg)
colnames(agg_write2)[1] <- 'name'
write_xlsx(agg_write2,'data/rel_abundance.xlsx')



### rarefy
raw <- read.table('data/table_ASVs.csv', sep = ';', stringsAsFactors = F, header = T)
raw <- raw[,1:38] # remove taxonomy
# rarefy (just do once)
rare <- t(rrarefy(t(raw),3000))
write_xlsx(data.frame(rare, stringsAsFactors = F),'data/rarefied_3000.xlsx')


