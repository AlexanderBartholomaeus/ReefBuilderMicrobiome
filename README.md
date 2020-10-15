# The microbiome associated with the reef builder Neogoniolithon sp. in the East Mediterranean

This repository support the submission of study performed by *The microbiome associated with the reef builder Neogoniolithon sp. in the East Mediterranean*. We will provide code and data that was used to generate results and plots.

## Link to publication 

... will be added as soon as possible 

## Structure of this repository

### Data

The data folder holds different data files:

* ASV raw table including taxonomic assignment a result from the DADA2 pipeline.
* Relative abundance ASV table we samples averaged by mean
* Rarefied ASV table

The raw sequencing files are stored at ENA ...link...

### Code

We are providing the code for the data processing and generation of many of the main figures. Please note: some plots were generated with the help of a R/Shiny app using the packages described in our publication. The app will be open-source and hopefully published soon (manuscript is in preparation). The code here is logically divided into different steps:

* DADA2 
* Data processing
* NMDS
* Alpha diversity
