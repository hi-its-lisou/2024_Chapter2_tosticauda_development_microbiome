#load libraries
library(phyloseq)
library(tidyverse)
library(qiime2R)
library(decontam)
library(microbiome)
library(ggh4x)
library(ggtext)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(scales)

# Load in sequencing data
ps <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_4.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps

################################################################################
#####          calculating the % of chloroplasts and mitochondria           ####
################################################################################

# Function to computer percentage
compute_percentage <- function(taxon, ps) {
  counts_data <- as.data.frame(otu_table(ps))  #Extract the counts data from the phyloseq object
  taxon_row <- which(tax_table(ps)[, "Genus"] == taxon)  #Identify the row corresponding to the specific taxon in the taxonomic table
  taxon_counts <- sum(counts_data[taxon_row, ])  #Sum the counts of the specific taxon across all samples
  total_counts <- sum(counts_data)  #Calculate the percentage of reads attributed to the specific taxon
  percentage_taxon <- (taxon_counts / total_counts) * 100
  return(percentage_taxon)
}

# Sublist sample types
Prepupae <- subset_samples(ps, sample_type %in% c("Prepupae")& !(sample_data(ps)$AB_treatment %in% c("control", "antibiotic")))
Larva <- subset_samples(ps, sample_type %in% c("Larvae"))
Honey_bee <- subset_samples(ps, sample_type %in% c("Honey_bee"))
Food <- subset_samples(ps, sample_type %in% c("Food"))
AdultsFree <- subset_samples(ps, Env_exposure %in% c("Free-flying")) #seperate only natural adult samples

antibiotic <- subset_samples(ps, AB_treatment %in% c("antibiotic"))
control <- subset_samples(ps, AB_treatment %in% c("control"))
AdultsNest <- subset_samples(ps, Env_exposure %in% c("Nest_only"))
AdultsNone <- subset_samples(ps, Env_exposure %in% c("None"))

# Run function for each sublisted sample type
# For Chloroplast
percentage_Chloroplast_Adults <- compute_percentage("Chloroplast", AdultsFree) 
percentage_Chloroplast_Larva <- compute_percentage("Chloroplast", Larva)
percentage_Chloroplast_Prepupae <- compute_percentage("Chloroplast", Prepupae)
percentage_Chloroplast_Food <- compute_percentage("Chloroplast", Food)
percentage_Chloroplast_HB <- compute_percentage("Chloroplast", Honey_bee)

percentage_Chloroplast_Adults
percentage_Chloroplast_Food
percentage_Chloroplast_Larva
percentage_Chloroplast_Prepupae
percentage_Chloroplast_HB

percentage_Chloroplast_antibiotic <- compute_percentage("Chloroplast", antibiotic) 
percentage_Chloroplast_control <- compute_percentage("Chloroplast", control) 
percentage_Chloroplast_AdultsNest <- compute_percentage("Chloroplast", AdultsNest)
percentage_Chloroplast_AdultsNone <- compute_percentage("Chloroplast", AdultsNone)

percentage_Chloroplast_antibiotic
percentage_Chloroplast_control
percentage_Chloroplast_AdultsNest
percentage_Chloroplast_AdultsNone

# For Mitochondria
percentage_Mitochondria_Adults <- compute_percentage("Mitochondria", AdultsFree) 
percentage_Mitochondria_Larva <- compute_percentage("Mitochondria", Larva)
percentage_Mitochondria_Prepupae <- compute_percentage("Mitochondria", Prepupae)
percentage_Mitochondria_Food <- compute_percentage("Mitochondria", Food)
percentage_Mitochondria_HB <- compute_percentage("Mitochondria", Honey_bee)

percentage_Mitochondria_Adults
percentage_Mitochondria_Food
percentage_Mitochondria_Larva
percentage_Mitochondria_Prepupae
percentage_Mitochondria_HB

percentage_Mitochondria_antibiotic <- compute_percentage("Mitochondria", antibiotic) 
percentage_Mitochondria_control <- compute_percentage("Mitochondria", control) 
percentage_Mitochondria_AdultsNest <- compute_percentage("Mitochondria", AdultsNest)
percentage_Mitochondria_AdultsNone <- compute_percentage("Mitochondria", AdultsNone)

percentage_Mitochondria_antibiotic
percentage_Mitochondria_control
percentage_Mitochondria_AdultsNest
percentage_Mitochondria_AdultsNone

# Proportion of contaminate asvs (called using decontam at a threshold of 0.5)
#vector for contaminant asvs
contam_asvs <-c("e5701abae521c279d3c45df09d57b2db",
                              "57ef5006312a9662d87d57c1404c1a87",
                              "37e52baccd1ddbf9c3c7d5b2a4f97f98",
                              "7fae0a4914d95fd2426aa8a2ea9ef487",
                              "54580fb940f267b8e94283e28ad1f13f",
                              "2df836c908b2263e29287bd2a7594b02",
                              "a53c599ad2ea0934852c9c0a029bc2b6",
                              "9bfe08c51b7a63c8130c31cf1c4c5665",
                              "c53215f72a5f346870f7e81352cdef29",
                              "c709213e70b1275f0dd0d033f26f0c1a",
                              "7ea68727c4e0241bdd33c1ec09ef87bb",
                              "438dc2a2ea2405a782df5830ae62100e",
                              "ab6a8e2af52c17a984b49aa8fb6f60c8",
                              "6c891abaa8f2147383dd332e601800eb",
                              "8bc54d575f88d09b2ba16b7f571812ae",
                              "30554407069b2367b77ee68c1982db14",
                              "f6444ec21f62323b6f0c08c4a7e48a07",
                              "ff2bd29ff42e4dc25a31714e0b6c2dca",
                              "0345b451b2a01dbc32e694d5310d02aa",
                              "73bf8d1a5983e34a0cb84e3cae127815",
                              "a5189f77a2cfeab3bc1602ff5c8ac3e9",
                              "5337f7e5d0db185e5513265bf151a2ba",
                              "1a819f1e863b7b7ee819b90aa0f59a1e",
                              "a5d1217aa1cbc154c499796f9040d954",
                              "d23fbef2f31d48eda40876cdbc49933a",
                              "b39c338f5e964b6cb87e07f10badc6c4",
                              "68cdc2abe341d7a696241454566de4fd",
                              "6ba7d4f61fd8143e8aa856596c09b9dd",
                              "82c18677dafa3465b9436d4699bebaa5",
                              "dbac3cc914a4a00d24412eefd5c894d1",
                              "ee02b9d7321d96af270ccd201ba17b6b",
                              "c873733a1c472c123f8063297a02c773",
                              "081432aa78080cbd9570800e938e256c",
                              "588b5ccdde9b0d39d61568731d7e223f",
                              "50a3b56e0b7db75ee9daccf7a751ba41",
                              "1f878f615fcfc8d7bd381a7841ac1e41",
                              "b158e58dd489c672f709dcb0963abdb2",
                              "0dc8e0e89e89cf4f0c71808666f6c6e7",
                              "2f86aa90c73ac3f94160b89f1f225ea1",
                              "e27680d4009f98f30248d823bc17fb8e",
                              "1b158b8b2922d4fcad5d9cea607cbb7d")

all_asvs <- taxa_names(ps)
asvs_to_keep <- all_asvs[!(all_asvs %in% contam_asvs)]
ps_no_contam <- prune_taxa(asvs_to_keep, ps)
ps_no_contam #pyloseq object without any contaminant asvs

#Adult tostis contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, Env_exposure %in% c("Free-flying")))))) #Specifically for free-fying the adults
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, Env_exposure %in% c("Free-flying"))))))
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

#Pollen provision contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Food")))))) 
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, sample_type == "Food")))))
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

#Larvae contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Larvae"))))))
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, sample_type == "Larvae")))))
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

#Prepupae contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == "Prepupae" & !(sample_data(ps)$AB_treatment %in% c("control", "antibiotic")))))))
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, sample_type == "Prepupae" & !(sample_data(ps)$AB_treatment %in% c("control", "antibiotic"))))))) 
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

#Honey bee contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Honey_bee")))))) 
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, sample_type == "Honey_bee"))))) 
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

# Antibiotic prepupae contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, AB_treatment %in% c("antibiotic"))))))
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, AB_treatment %in% c("antibiotic"))))))
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

# Control prepupae contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, AB_treatment %in% c("control"))))))
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, AB_treatment %in% c("control"))))))
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

#Nest eclosed bees contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, Env_exposure %in% c("Nest_only"))))))
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, Env_exposure %in% c("Nest_only"))))))
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

#Controlled eclosion bees contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, Env_exposure %in% c("None"))))))
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, Env_exposure %in% c("None"))))))
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

#Working biomass after removing all non-target reads
#vector for chloroplast, mitochondrial, and contaminant asvs
# Remove contaminants, chloroplasts, and mitochondria ####
chloro_mito_decontam_asvs <- read_delim("input_files/chloro_mito_decontam_asvs.txt", 
                     delim = "\n", 
                     col_names = "asv")
chloro_mito_decontam_asvs <- chloro_mito_decontam_asvs$asv
all_asvs <- taxa_names(ps)
asvs_to_keep <- all_asvs[!(all_asvs %in% chloro_mito_decontam_asvs)]
ps_filtered <- prune_taxa(asvs_to_keep, ps)
ps_filtered #phyloseq object with all off-target reads removed

#Specifically for the adults
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, Env_exposure %in% c("Free-flying"))))))
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered,Env_exposure %in% c("Free-flying"))))))
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

#Larvae working biomass
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Larvae")))))) #replace with sample_type in question
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, sample_type == "Larvae"))))) #replace with sample_type in question
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

#Prepupae working biomass
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Prepupae") & !(sample_data(ps)$AB_treatment %in% c("control", "antibiotic"))))))) 
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, sample_type == "Prepupae" & !(sample_data(ps_filtered)$AB_treatment %in% c("control", "antibiotic")))))))
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

#Food working biomass
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Food")))))) #replace with sample_type in question
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, sample_type == "Food"))))) #replace with sample_type in question
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

#honey bee working biomass
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Honey_bee")))))) #replace with sample_type in question
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, sample_type == "Honey_bee"))))) #replace with sample_type in question
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining


# working biomass
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, AB_treatment %in% c("antibiotic"))))))
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, AB_treatment %in% c("antibiotic"))))))
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, AB_treatment %in% c("control"))))))
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, AB_treatment %in% c("control"))))))
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, Env_exposure %in% c("Nest_only"))))))
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, Env_exposure %in% c("Nest_only"))))))
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, Env_exposure %in% c("None"))))))
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, Env_exposure %in% c("None"))))))
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

