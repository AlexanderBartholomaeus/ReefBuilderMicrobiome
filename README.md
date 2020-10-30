# The microbiome associated with the reef builder *Neogoniolithon sp.* in the East Mediterranean

This repository supports the submission of study **The microbiome associated with the reef builder** ***Neogoniolithon sp.*** **in the East Mediterranean**. We will provide code and data that was used to generate results and plots.

### Link to publication

Currently we are working on the revision.

### Data

The data folder holds different data files:

* **table_ASVs.csv**: ASV raw table including taxonomic assignment a result from the DADA2 pipeline. There is taxonomy assignmend based on SILVA 138 (columns 39 to 44) and also SILVA 132 (columns 45 to 50). Note this table still contains chloroplast and mitochondria ASVs which are removed during the data processing. 
* **rarefied_3000.xlsx**: Rarefied ASV t 
* **XXX**: Relative abundance ASV table we samples averaged by mean
* **alphaDiv_indices.xlsx**: Alpha diversity indices for each sample generated from the rarefied_3000.xlsx (see script [alpha_diversity.R lines 21-27](https://github.com/AlexanderBartholomaeus/ReefBuilderMicrobiome/blob/main/alpha_diversity.R#L21-L27). 
* ***table_sample_stats.xlsx***: This table contains the basic statistics of DADA2 pipeline and include a sample description.

The raw sequencing data files are stored at ENA  https://www.ebi.ac.uk/ena/browser/view/PRJEB38881.

### Code

We are providing the code for the data processing and generation of many of the main figures. Please note: some plots were generated with the help of a R/Shiny app using the packages described in our publication. The app will be open-source and hopefully published soon (manuscript is in preparation). The code here is logically divided into different steps:

* Adapter and quality trimming of the raw reads: adapter_quality_trimming.R
* DADA2: dada2.R
* Data processing: process_data.R
* NMDS and bubbleplots: nmds_bubbleplots.R
* Alpha diversity: alpha_diversity.R

In the following subsection the different files will be shortly presented

#### adapter_quality_trimming.R



#### dada2.R
The 

#### alpha_diversity.R

