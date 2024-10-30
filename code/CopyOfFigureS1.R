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
food_with_brood_through_time <- subset_samples(filtered_physeq, sample_type %in% c("Food") 
                                    & Nest_id %in% c("5", "13", "14", "27")
                                    & Cell_ID %in% c("J", "I", "H", "G", "F", "E", "D", "C", "B", "A"))
food_with_brood_through_time@sam_data$Nest_id[which(food_with_brood_through_time@sam_data$Nest_id == "5")] <- "1"
food_with_brood_through_time@sam_data$Nest_id[which(food_with_brood_through_time@sam_data$Nest_id == "13")] <- "2"
food_with_brood_through_time@sam_data$Nest_id[which(food_with_brood_through_time@sam_data$Nest_id == "14")] <- "3"
food_with_brood_through_time@sam_data$Nest_id[which(food_with_brood_through_time@sam_data$Nest_id == "27")] <- "4"



# Alpha diversity of pollen provisions as they are consumed ####
alpha_diversity <- alpha(food_with_brood_through_time, index = "Shannon")
metadata <- meta(food_with_brood_through_time)
metadata$name <- rownames(metadata)
alpha_diversity$name <- rownames(alpha_diversity)
alpha_diversity_metadata <- merge(alpha_diversity, metadata, by = "name")

# Set the factor levels for Cell_ID in the merged dataframe
custom_Cell_ID_order <- c("J", "I", "H", "G", "F", "E", "D", "C", "B", "A")
alpha_diversity_metadata$Cell_ID <- factor(alpha_diversity_metadata$Cell_ID, levels = custom_Cell_ID_order)

# Create the plot with regression line and R-squared annotation
boxplot <- alpha_diversity_metadata %>%
  ggplot(aes(x = Cell_ID, y = diversity_shannon, colour = Cell_ID)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(y = "Shannon diversity", x = "Cell position (youngest to oldest)") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14, face = 'bold'),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = 'bold'),
        legend.position = "none")

# Convert levels to numeric values for regression
numeric_cell_ids <- as.numeric(factor(alpha_diversity_metadata$Cell_ID, levels = custom_Cell_ID_order))

# Create the plot
ggplot(alpha_diversity_metadata, aes(x = numeric_cell_ids, y = diversity_shannon)) +
  geom_point(size = 5, color = "black") +  
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  labs(y = "Shannon diversity", x = "Cell position (youngest to oldest)") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14, face = 'bold'),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = 'bold')) +
  facet_wrap(~ Nest_id) +
  scale_x_continuous(breaks = 1:length(custom_Cell_ID_order), 
                     labels = custom_Cell_ID_order)


ggsave("figures/S1.png", height=8, width=12)


model_interaction <- lm(diversity_shannon ~ as.numeric(Cell_ID) * Nest_id, data = alpha_diversity_metadata)
summary(model_interaction)
