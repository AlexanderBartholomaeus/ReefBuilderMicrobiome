# load libraries
library(readxl)
library(writexl)
library(ggplot2)
library(phyloseq)
library(vegan)


# load meta / sample data
meta <- data.frame(read_excel('data/table_sample_stats.xlsx', skip=1), stringsAsFactors = F)

# load rarefied 
rare <- data.frame(read_excel('data/rarefied_3000.xlsx'))

# check 0 rows
sum(apply(rare,1,sum)==0)

# do a first plot to get a feeling for the data
plot_richness(phyloseq(otu_table(rare,taxa_are_rows = T)))

# calculate alpha diversity indices (just once)
#aDiv <- estimate_richness(phyloseq(otu_table(rare,taxa_are_rows = T)))
#aDiv <- cbind(rownames(aDiv),aDiv)

# merge with sample data (just once)
#merged <- merge(aDiv,meta[,7:12], by.x=1, by.y=1)
#write_xlsx(data.frame(merged, stringsAsFactors = F),'data/alphaDiv_indices.xlsx')

# load 
merged <- data.frame(read_excel('data/alphaDiv_indices.xlsx'))

# without aquariaum 
merged_2 <- merged[merged$Location!='Aquarium',]
merged_2 <- merged_2[merged_2$Tissue_water=='Algal',]
merged_2$sampling_location <- paste(merged_2$Sampling,'-',merged_2$Location)
merged_3 <- merged[merged$Location!='Aquarium',]
merged_3$sampling_location <- paste(merged_3$Sampling,'-',merged_3$Location)
merged_3$sampling_location_tissue <- paste(merged_3$Tissue_water,'-',merged_3$Sampling,'-',merged_3$Location)


# plots -------------------------------------------------------------------

# path to store 
plotPath = 'rarefied/' # path to store plot files

# for each alpha diversity index 
for(i in c('Shannon','Observed','Chao1','ACE','Simpson','InvSimpson','Fisher')){
  # all + colored
  ggplot(merged,aes_string(x='Tissue_water',y=i)) + 
    geom_boxplot(lwd = 0.9, outlier.shape = NA) + 
    geom_jitter(aes(color = Location, shape = Sampling),position=position_jitter(0.1), size = 4) +
    xlab('') + 
    ylab(i) +
    theme_light()
  ggsave(paste0(plotPath,'/boxplot_algal_water_all_',i,'.svg'),width = 5, height = 4)
  ggsave(paste0(plotPath,'/boxplot_algal_water_all_',i,'.png'),width = 5, height = 4)

  # no aquarium
  ggplot(merged_3,aes_string(x='Tissue_water',y=i)) + 
    geom_boxplot(lwd = 0.9, outlier.shape = NA) + 
    geom_jitter(shape=16, position=position_jitter(0.1), size = 4, color = '#00000080') +
    xlab('') + 
    ylab(i) +
    theme_light()
  ggsave(paste0(plotPath,'/boxplot_algal_water_noAquarium_',i,'.svg'),width = 4, height = 4)
  ggsave(paste0(plotPath,'/boxplot_algal_water_noAquarium_',i,'.png'),width = 4, height = 4)
  # no aquarium + colored
  ggplot(merged_3,aes_string(x='Tissue_water',y=i)) + 
    geom_boxplot(lwd = 0.9, outlier.shape = NA) + 
    geom_jitter(aes(color = Location, shape = Sampling),position=position_jitter(0.1), size = 4) +
    xlab('') + 
    ylab(i) +
    theme_light()
  ggsave(paste0(plotPath,'/boxplot_algal_water_noAquarium_color_',i,'.svg'),width = 5, height = 4)
  ggsave(paste0(plotPath,'/boxplot_algal_water_noAquarium_color_',i,'.png'),width = 5, height = 4)
  
  # only algal and sampling
  ggplot(merged_2,aes_string(x='Sampling',y=i)) + 
    geom_boxplot(lwd = 0.9, outlier.shape = NA) + 
    geom_jitter(shape=16, position=position_jitter(0.1), size = 4, color = '#00000080') +
    xlab('') + 
    ylab(i) +
    theme_light()
  ggsave(paste0(plotPath,'/boxplot_algal_sampling_',i,'.svg'),width = 4, height = 4)
  ggsave(paste0(plotPath,'/boxplot_algal_sampling_',i,'.png'),width = 4, height = 4)
  
  # only algal and location
  ggplot(merged_2,aes_string(x='Location',y=i)) + 
    geom_boxplot(lwd = 0.9, outlier.shape = NA) + 
    geom_jitter(shape=16, position=position_jitter(0.1), size = 4, color = '#00000080') +
    xlab('') + 
    ylab(i) +
    theme_light()
  ggsave(paste0(plotPath,'/boxplot_algal_location_',i,'.svg'),width = 4, height = 4)
  ggsave(paste0(plotPath,'/boxplot_algal_location_',i,'.png'),width = 4, height = 4)
  
  # tide pool vs edge in each sampling
  ggplot(merged_2,aes_string(x='sampling_location',y=i)) + 
    geom_boxplot(lwd = 0.9, outlier.shape = NA) + 
    geom_jitter(shape=16, position=position_jitter(0.1), size = 4, color = '#00000080') +
    xlab('') + 
    ylab(i) +
    theme_light() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave(paste0(plotPath,'/boxplot_algal_sampling_location_',i,'.svg'),width = 6, height = 5)
  ggsave(paste0(plotPath,'/boxplot_algal_sampling_location_',i,'.png'),width = 6, height = 5)
  
  # algal and water with tide pool vs edge in each sampling
  ggplot(merged_3,aes_string(x='sampling_location_tissue',y=i)) + 
    geom_boxplot(lwd = 0.9, outlier.shape = NA) + 
    geom_jitter(shape=16, position=position_jitter(0.1), size = 4, color = '#00000080') +
    xlab('') + 
    ylab(i) +
    theme_light() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  ggsave(paste0(plotPath,'/boxplot_algal_water_sampling_location_',i,'.svg'),width = 7, height = 5)
  ggsave(paste0(plotPath,'/boxplot_algal_water_sampling_location_',i,'.png'),width = 7, height = 5)
}


# statistical tests --------------------------------------------------------

# perform statistical tests
summary(aov(merged_3$Shannon~merged_3$Tissue_water)) # ANOVA algal vs water
summary(aov(merged_2$Shannon~merged_2$Location)) # ANOVA tide pool vs edge (algal only)
summary(aov(merged_2$Shannon~merged_2$Sampling)) # ANOVA sampling1 vs sampling2 (algal only)
# more detailed test
TukeyHSD(aov(merged_2$Shannon~merged_2$sampling_location)) # sampling & location (algal only)
TukeyHSD(aov(merged_3$Shannon~merged_3$Tissue_water))
TukeyHSD(aov(merged_3$Shannon~merged_3$Tissue_water*merged_3$Sampling)) 
TukeyHSD(aov(merged_3$Shannon~merged_3$Tissue_water*merged_3$Location))
TukeyHSD(aov(merged_3$Shannon~merged_3$Tissue_water*merged_3$Location*merged_3$Sampling))

# remove winter tide
merged_3_1 <- merged_3[merged_3$sampling_location_tissue!='Algal - Winter - Edge',]
TukeyHSD(aov(merged_3_1$Shannon~merged_3_1$Tissue_water))
# add winter tide pool to water  
merged_3_2 <- merged_3
merged_3_2$Tissue_water[merged_3_2$sampling_location_tissue == 'Algal - Winter - Edge'] <- 'Water'
TukeyHSD(aov(merged_3_2$Shannon~merged_3_2$Tissue_water))
