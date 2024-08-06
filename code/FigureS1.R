# Load the required libraries ####
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
library(extrafont)
library(forcats)
library(ggforce)
library(vegan)
library(xfun)

# Load phyloseq object ####
ps <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_4.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps

# Load custom colour palette ####
colours_df <- file.path(
  "D:/Research/PhD/Manuscripts/Chapter 2/Git/2024_Chapter2_tosticauda_development_microbiome/input_files/",
  "colour_list.csv"
) %>% read_csv

my_palette <- colours_df$colours
names(my_palette) <- colours_df$genera
my_palette

# Rarefy ####
ps_rar <- rarefy_even_depth(ps, sample.size = 1500, rngseed = 1337)

# Remove contaminants, chloroplasts, and mitochondria ####
contam <- read_delim("input_files/chloro_mito_decontam_asvs.txt", 
                     delim = "\n", 
                     col_names = "asv")
chloro_mito_decontam_asvs <- contam$asv
all_asvs <- taxa_names(ps_rar)
asvs_to_keep <- all_asvs[!(all_asvs %in% chloro_mito_decontam_asvs)]
ps_rar <- prune_taxa(asvs_to_keep, ps_rar)
ps_rar

# Filter out 0 abundance reads ####
zero_abundance_samples <- sample_sums(ps_rar) == 0
zero_abundance_samples <- sample_sums(ps) == 0

filtered_physeq <- subset_samples(ps_rar, !zero_abundance_samples)
filtered_physeq <- subset_samples(ps, !zero_abundance_samples)
filtered_physeq

# Subset the data for samples_types ####
food_through_time <- subset_samples(filtered_physeq, sample_type %in% c("Food") 
                                    & Nest_id %in% c("5", "13", "14", "27")
                                    & Cell_ID %in% c("J", "I", "H", "G", "F", "E", "D", "C", "B", "A"))
food_through_time@sam_data$Nest_id[which(food_through_time@sam_data$Nest_id == "5")] <- "1"
food_through_time@sam_data$Nest_id[which(food_through_time@sam_data$Nest_id == "13")] <- "2"
food_through_time@sam_data$Nest_id[which(food_through_time@sam_data$Nest_id == "14")] <- "3"
food_through_time@sam_data$Nest_id[which(food_through_time@sam_data$Nest_id == "27")] <- "4"

# Alpha diversity of pollen provisions as they are consumed ####
alpha_diversity <- alpha(food_through_time, index = "Shannon")
metadata <- meta(food_through_time)
metadata$name <- rownames(metadata)
alpha_diversity$name <- rownames(alpha_diversity)
alpha_diversity_metadata <- merge(alpha_diversity, metadata, by = "name")

# Box plot for alpha diversity ####
alpha_diversity_metadata %>%
  ggplot(aes(x = Cell_ID, y = diversity_shannon, colour = Cell_ID)) +
  geom_boxplot() +
  labs(x="Shannon diversity", y = "Cell positione (oldest to youngest)")+
  geom_point(size = 3)

ggsave("figures/S1.png", height=10, width=15)