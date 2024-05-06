setwd("D:/Research/PhD/Manuscripts/Chapter 2/code/R project/Chapter2_analysis/data")

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
  features = "merged_table.qza",
  taxonomy = "merged_taxonomy.qza",
  tree = "merged_sepp_tree.qza",
  metadata = "combined_metadata_4.tsv")
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

# Run function for each sublisted sample type
# For Chloroplast
percentage_Chloroplast_Adults <- compute_percentage("Chloroplast", AdultsFree) 
percentage_Chloroplast_Larva <- compute_percentage("Chloroplast", Larva) #replace with sample_type in question
percentage_Chloroplast_Prepupae <- compute_percentage("Chloroplast", Prepupae)
percentage_Chloroplast_Food <- compute_percentage("Chloroplast", Food)
percentage_Chloroplast_HB <- compute_percentage("Chloroplast", Honey_bee)

percentage_Chloroplast_Adults
percentage_Chloroplast_Food
percentage_Chloroplast_Larva
percentage_Chloroplast_Prepupae
percentage_Chloroplast_HB

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
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, Env_exposure %in% c("Free-flying")))))) #Specifically for freelfying the adults
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, Env_exposure %in% c("Free-flying"))))))
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

#Pollen provision contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Food")))))) 
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, sample_type == "Food")))))
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

#Larvae contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Larvae")))))) #replace with sample_type in question
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, sample_type == "Larvae"))))) #replace with sample_type in question
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


#Working biomass after removing all non-target reads
#vector for chloroplast, mitochondrial, and contaminant asvs
chloro_mito_decontam_asvs <-c("fdf625d4c0cc6895a87a0023e3f52a18",
                              "ddbb3a2a7d6486a2201bb8770feb0a6d",
                              "d499ea4449d3cffa7a5acfbd226fc1bb",
                              "5502d5699b37790ac43b2b75532a621b",
                              "ce0fe9dd4a105852b029c2ab1d997024",
                              "f61bb013d7b0348b65975072ceb09619",
                              "1ab5f68f3e755878ca8afdb5563c218f",
                              "a3856381d1f19248dbe50d2c16724665",
                              "2ae8d0f7c8332e3df1366d60f0726123",
                              "713eae4c577dd5b3f8a4448c1173b23b",
                              "c3463ea4a604e0fd848ca2ceda06ca9e",
                              "ad4056eac2e7cbef811432e72eed6f52",
                              "3137676e72277173e6b9e9630ef18bee",
                              "54459e59e8344a06665a8e2d734bd918",
                              "b74fb77d1e3038ad7239c9c9a4dd36be",
                              "ed08849f037a13ec071c44366d8029a8",
                              "0a650376bddbd340d0494907ab3c892c",
                              "a877d39eed7fd090191019d61b4aa1ec",
                              "a17b7b0b059c0ff6f3a50c30272ec893",
                              "11be956e8b33ae5031de5d3be792c33e",
                              "d4e1bad156b1cac400e6d00c5a83512b",
                              "8f1cc7f72ac6c8542db82222104ce32e",
                              "794abc83beea78b6b6b44ee2043a3746",
                              "149ff2029c30eb2342d05344d9af29db",
                              "36bb717d13e68ed4719be14301b5e1d1",
                              "54e7b2139640ec28fbaf7f0491cad845",
                              "705d53c2b383d3d3a804afa8fbb2db8a",
                              "1bc1c2d08bb32d49f6095c52e6ac1876",
                              "2e74be53b6ee8ded1af833202c19931b",
                              "74eb42089635632919bff035677eacf9",
                              "da3a12b9d65930c33d447d14adf8dd1a",
                              "922d59e2ac30e517253bbcee836ceb16",
                              "e4a4f2899fd5a101eb3451507b24423f",
                              "0f8b275c51a0f7b7d9c53748b4daca03",
                              "966558565f78941a9f7329b84f3e02f1",
                              "e66741e5f1b594e738ef4ffa5a3d328d",
                              "7b5f0042ef3d5b32faee86f4b10578a6",
                              "bc87d429952f1b6c15580c6272af8d53",
                              "7312bb1dc386eb1e89c0a902eef94f47",
                              "d05d070d32c6213a6582f7915690c2fb",
                              "07b4d1ff39fc193ca91d28378707dd0e",
                              "d0bc2ccfe65e8cba6c1cf36cf84b12e8",
                              "af347ab82f76c2802339bcfcad64e876",
                              "ddf35962b88c833d12a42c0855b133b5",
                              "08578fc8fe4013ffae56edcb8aee5a4d",
                              "b8110a61ce2f6af058221eb370b215c2",
                              "630f779d78e36bc0ecbaaaccb1b4c844",
                              "d7ec4030ebd1eaf630f6cad5830e9662",
                              "a926e2dba6582506417b2b7ac1cad5cb",
                              "006447e914f38b29a86b92aab7b8145d",
                              "b84ee2364b00f6073b81ce44217504f5",
                              "0251c873b9e6e76973408b336adeb8b2",
                              "371a76e8c18bf7811c0f1299ce9fcbec",
                              "47b2b1223d7d4207e0340988d37ca900",
                              "1c9e5f2d1e7fd315c8346cd76b4bfdb6",
                              "7323f28da79c43975261bf628d9f4a8d",
                              "e5716aaf5d7bf569d21090c02de1ffe9",
                              "e9cbfc116f360488c1511eef4c7f19b9",
                              "b1eb9522d23662657cd5b9e02a96e347",
                              "2efbf01909c1b8c4b79ae8655c50b949",
                              "1a445cb7db38c2604dc792dbe91e059f",
                              "86b394556e3ff629b18106a842affe00",
                              "b220bf517e5078694409d128aa0fc5d9",
                              "fe8b52a854181e2e440d42ffa2cd83a8",
                              "88a927a190b1423b1406f3746eabd0b6",
                              "42d06841f5de3aca03238a5781452ebc",
                              "0236242c1bd991ed5ff2b072739c9a8f",
                              "636b1bad013fdad596fe39fd0216202f",
                              "784a38a560dab9f2dde36b735f3670b5",
                              "102c2423ad57d60605c7e3aa81128349",
                              "94ee97f7ac645cd3312e1401cfdbea00",
                              "0ebe3aa0ccdf468c4e6f2b0bda6ad67d",
                              "034e30578bf25dfb6a667d3db2758c1c",
                              "eb9d866bd6acfdd8096b50bacad65a8e",
                              "0767cb7cae16a55313f0707f602d4a00",
                              "977e31b3f58a7a3a705456ef4c63a0a1",
                              "1884fb0c3e9007ffc8af2625a5e7f87c",
                              "9d5aacb2cd576a8433985104b169fd46",
                              "6d226ab3a7ca4eaef04f5dbec99ac59d",
                              "c537e706d258bf535dbe923d4c13b750",
                              "55e9d859cd6814527e42d651718bd62b",
                              "52744a971602f4982249e89c38952713",
                              "226e953a1f4c8fcad650247be3807035",
                              "5323c6dd766d629cc750b0aabb4aaf72",
                              "59dda3111163eb7affe703f9eee3db9e",
                              "b0a6c67ebcf4f78f03ccf9c7bb037087",
                              "174be33d93aa22d41599ee9ab5663fad",
                              "e5701abae521c279d3c45df09d57b2db",
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