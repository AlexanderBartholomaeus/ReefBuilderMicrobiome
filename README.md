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
* Alpha diversity: alpha_diversity.R
* NMDS and bubbleplots: nmds_bubbleplot.R

In the following subsection the different files will be shortly presented

#### adapter_quality_trimming.R

This script uses the raw sequencing files (download [here](https://www.ebi.ac.uk/ena/browser/view/PRJEB38881)) and perfoms adapter and quality trimming including a minimal length filtering. To execute the [cutadapt](https://cutadapt.readthedocs.io/en/stable/) tool and the R package [seqRFLP](https://github.com/helixcn/seqRFLP) are required required.

#### DADA2

The [dada2.R](https://github.com/AlexanderBartholomaeus/ReefBuilderMicrobiome/blob/main/dada2.R) script uses the [DADA2](https://benjjneb.github.io/dada2/index.html) pipeline published in [nature methods](https://www.nature.com/articles/nmeth.3869) to generate Amplicon Sequence Variants (ASVs). In addition taxonomic assigment of the ASVs is done with SILVA 138. This results in a [ASV table](https://github.com/AlexanderBartholomaeus/ReefBuilderMicrobiome/blob/main/data/table_ASV.csv). The script requires R packages [dada2](https://benjjneb.github.io/dada2/dada-installation.html).

#### process_data.R

The [process_data.R](https://github.com/AlexanderBartholomaeus/ReefBuilderMicrobiome/blob/main/process_data.R) script removes chloroplast and mitochondria ASVs, generates [rarefied data](https://github.com/AlexanderBartholomaeus/ReefBuilderMicrobiome/blob/main/data/rarefied_3000.xlsx) and generates the [relative abundance data](https://github.com/AlexanderBartholomaeus/ReefBuilderMicrobiome/blob/main/data/rel_abundance.xlsx). 


#### alpha_diversity.R

The [alpha_diversity.R](https://github.com/AlexanderBartholomaeus/ReefBuilderMicrobiome/blob/main/alpha_diversity.R) script uses rarefied data to generate alpha diversity plots and give statistical measures. [Lines 21-27](https://github.com/AlexanderBartholomaeus/ReefBuilderMicrobiome/blob/main/alpha_diversity.R#L21-L27) are outcommented as alpha diversity measures need to be calculated only once. Alpha diversity measures of each sample are available [here](https://github.com/AlexanderBartholomaeus/ReefBuilderMicrobiome/blob/main/data/alphaDiv_indices.xlsx). The script needs the R packages [readxl](https://readxl.tidyverse.org/), [writexl](https://github.com/ropensci/writexl), [ggplot2](https://ggplot2.tidyverse.org/), [phyloseq](https://joey711.github.io/phyloseq/), [vegan](https://github.com/vegandevs/vegan).

#### NMDS bubble

The [nmds_bubbleplot.R](https://github.com/AlexanderBartholomaeus/ReefBuilderMicrobiome/blob/main/nmds_bubbleplot.R) script generates NMDS and bubbleplot figures. Bubbleplot generation need code from the https://github.com/AlexanderBartholomaeus/BubblePlot is used. The script relies on various R packages (see [Lines 1-15](https://github.com/AlexanderBartholomaeus/ReefBuilderMicrobiome/blob/main/nmds_bubbleplot.R#L1-L15)).
