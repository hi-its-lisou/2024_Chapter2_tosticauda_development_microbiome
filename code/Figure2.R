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
filtered_physeq <- subset_samples(ps_rar, !zero_abundance_samples)
filtered_physeq

# Determine the number of reads per sample ####
sort(sample_sums(filtered_physeq))

# Subset the data for samples_types ####

frass_subset <- subset_samples(filtered_physeq, Description_of_sample %in% c("CocoonFrass", "Frass", "PinkCocoonFrass"))

food_through_time <- subset_samples(filtered_physeq, sample_type %in% c("Food") 
                                    & Nest_id %in% c("5", "13", "14", "27")
                                    & Cell_ID %in% c("J", "I", "H", "G", "F", "E", "D", "C", "B", "A"))
food_through_time@sam_data$Nest_id[which(food_through_time@sam_data$Nest_id == "5")] <- "1"
food_through_time@sam_data$Nest_id[which(food_through_time@sam_data$Nest_id == "13")] <- "2"
food_through_time@sam_data$Nest_id[which(food_through_time@sam_data$Nest_id == "14")] <- "3"
food_through_time@sam_data$Nest_id[which(food_through_time@sam_data$Nest_id == "27")] <- "4"

# Obtain top 20 genera ####
ps_Genus <- tax_glom(food_through_time, taxrank = "Genus", NArm = FALSE)
top20Genus = names(sort(taxa_sums(ps_Genus), TRUE)[1:20])
taxtab20 = cbind(tax_table(ps_Genus), Genus_20 = NA)
taxtab20[top20Genus, "Genus_20"] <- as(tax_table(ps_Genus)
                                       [top20Genus, "Genus"], "character")
tax_table(ps_Genus) <- tax_table(taxtab20)
ps_Genus_ra <- transform_sample_counts(ps_Genus, function(x) 100 * x/sum(x))
df_Genus <- psmelt(ps_Genus_ra)
df_Genus <- arrange(df_Genus, sample_type)
df_Genus$Genus_20[is.na(df_Genus$Genus_20)] <- c("Other")

# % of reads that make up the top 20 genera
mean(
  sample_sums(
    prune_taxa(top20Genus, ps_Genus_ra)
  )
)

# Plot the relative abundance ####
# Custom order for Cell_ID from youngest to oldest
custom_Cell_ID_order <- c("J", "I", "H", "G", "F", "E", "D", "C", "B", "A")
df_Genus$Cell_ID <- factor(df_Genus$Cell_ID, levels = custom_Cell_ID_order)

(fig3 <- df_Genus %>%
  mutate(Genus_20 = reorder(Genus_20, -Abundance)) %>%
  ggplot(aes(x = sample_type, y = Abundance, fill = Genus_20)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = my_palette) +
  facet_nested(~ Nest_id + Year + Cell_ID + sample_type, scales = "free", space = "free") +
  labs(x = "sample_type", y = "Relative abundance") +
  theme( 
    axis.text.y = element_text(size=16, face =, 'bold'),
    axis.title.y = element_text(size=16, face = 'bold'),
    axis.ticks.y = element_line(linewidth = 1),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_textbox_simple(
      padding = margin(5, 0, 5, 0),
      margin = margin(5, 5, 5, 5),
      size = 16,
      face = "bold",
      halign = 0.5,
      fill = "white",
      box.color = "grey",
      linewidth = 1.5,
      linetype = "solid",
      family = "Calibri"),
    panel.background = element_blank()
  ))

ggsave("figures/Figure2uneditted.png", height=10, width=20)

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
  geom_point(size = 3)

ggsave("figures/S1.png", height=10, width=15)

