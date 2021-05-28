#reset
rm(list=ls())

#loading libraries
#install_github("thomasp85/patchwork")
library(devtools)
library(tidyverse)
library(phyloseq)
library(vegan)
library(readxl)
library(ggplot2)


### Preparing data ###
#Data of crispr subtype abundances
load("~/data/CRISPRs.RData")

#data of bacteria abundances (mOTUs)
#load("~/data/mOTUs/motus2_phyloseq.RData")

#data of bacteria abundances (VAMB clusters)
load("~/data/VAMB/210415_VAMB_phyloseq_v02.RData")

#Checking the data:
#head(sample_data(final.physeq))
#sample_variables(final.physeq)
#sample_sums(final.physeq)
# 662 samples
#rank_names(final.physeq)
#get_variable(CRISPRs,"Supplement_oil")
#head(sample_data(CRISPRs))
#sample_variables(CRISPRs)


### Preparing sample data to be binary ###

# Extracting sample data from the phyloseq object:
samples_df <- get_variable(CRISPRs)

#Editing to be binary
samples_df <- samples_df %>% 
  #preserving the row names:
  rownames_to_column() %>% 
  #Having taken supplement is a 1, placebos are 0
  mutate(Supplement_oil = case_when(Supplement_oil == "Placebo" ~ 0,
                                    Supplement_oil == "Fiskeolie" ~ 1),
         Supplement_vitamin = case_when(Supplement_vitamin == "Placebo" ~ 0,
                                        Supplement_vitamin == "VitaminD" ~ 1),
         #Any type of sectio is 1, while normal is 0
         Delivery = case_when(Delivery == "Normal" ~ 0,
                              Delivery == "Planned sectio" ~ 1,
                              Delivery == "Acute sectio" ~ 1),
         #If the children spent 1 week or less in a year with animals, 
         #it will be categorized as "low" and therefore = 0
         Furred_animal_days = case_when(Furred_animal_days <= 7 ~ 0,
                                        TRUE ~ 1),
         #Pamcluster 2 becomes 1 and pamcluster 1 becomes 0
         pamcluster = as.numeric(pamcluster) - 1,
         #Rural becomes 0, urban becomes 1
         RuralUrbanStatus = case_when(RuralUrbanStatus == "Rural" ~ 0,
                                      RuralUrbanStatus == "Urban" ~ 1)
  )%>% 
  #2nd step for preserving row names
  column_to_rownames() %>% 
  #Deselect unnecessary variables not used for this analysis
  select(-"Birthdate", -"Sectiotype", -matches("days_with_"), -"maz")

#for checking if there are no variables left out of being 0/1:
#samples_df %>% group_by(Furred_animal_days) %>% summarise(N=n())

#Transform into phyloseq object
OTU = otu_table(CRISPRs)
TAX = tax_table(CRISPRs)
SAMPLES = sample_data(samples_df)

newCRISPRs <- phyloseq(OTU, TAX, SAMPLES)


### Relative abundance ###

#For doing relative abundance, relatively to each row's sum:
#install_github("JStokholm/Abundance")
library(abundance)

#taking relative abundance on bacteria table and CRISPR-Cas systems subtype table
#genus_mOTU_ab <- abundance(phylo_ob=mof1y, level="genus", id="Abcno",sample_id="SampleID", relative_abun=TRUE, remove_collapsed_taxa=FALSE, select_taxa=NULL,select_level=NULL)
species_VAMB_ab <- abundance(phylo_ob=final.physeq, level="species", id="abcno",sample_id="sampleid", relative_abun=TRUE, remove_collapsed_taxa=FALSE, select_taxa=NULL,select_level=NULL)
#genus_VAMB_ab <- abundance(phylo_ob=final.physeq, level="genus", id="abcno",sample_id="sampleid", relative_abun=TRUE, remove_collapsed_taxa=FALSE, select_taxa=NULL,select_level=NULL)
CRISPR_ab <- abundance(newCRISPRs,level="Subtype",id = "ABCno",sample_id="SampleName")



### Sorting and cleaning up data ###

#Ensuring we only have data in which ABC no. is in both tables
CRISPR_ab <- CRISPR_ab[CRISPR_ab$ABCno %in% species_VAMB_ab$abcno, ]
#removing duplicates of samples that belong to the same child
CRISPR_ab <- CRISPR_ab[!duplicated(CRISPR_ab$ABCno),]
#putting in order according to ABC number
CRISPR_ab <- CRISPR_ab[order(CRISPR_ab$ABCno),]

#doing the same for the bacterial data table
species_VAMB_ab <- species_VAMB_ab[species_VAMB_ab$abcno %in% CRISPR_ab$ABCno, ]
#removing duplicates of samples that belong to the same child
species_VAMB_ab <- species_VAMB_ab[!duplicated(species_VAMB_ab$abcno),]
#ordering
species_VAMB_ab <- species_VAMB_ab[order(species_VAMB_ab$abcno),]


### Transforming the abundance data ###

#Taking out the data with the bacterias' abundance (first 2 columns being ABC no. and sample ID)
species_abundOnly <- species_VAMB_ab[,3:ncol(species_VAMB_ab)]
#Ordering by median abundance to find the most abundant bacteria
species_abundOnly <- species_abundOnly[,order(-apply(species_abundOnly, 2, median))]
#Doing a log transformation with pseudocount 1 and multiplied by a million
species_abundOnly <- species_abundOnly %>%  apply(2, function(x) log(10**6 * x + 1))


#Doing the same for the CRISPR data
crispr_abundOnly <- CRISPR_ab[,3:ncol(CRISPR_ab)]
#Ordering by median
crispr_abundOnly <- crispr_abundOnly[,order(-apply(crispr_abundOnly, 2, median))] 
#Select all system subtypes and do log transformation of the abundance
crispr_abundOnly <- crispr_abundOnly %>%  apply(2, function(x) log(10**6 * x + 1))

#these should match:
nrow(species_abundOnly)
nrow(crispr_abundOnly)


# and now we're ready for analysis!