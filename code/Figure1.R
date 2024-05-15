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
library(patchwork)

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
ps_rar

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

# Subset the data for samples_types 
prepupal_subset <- subset_samples(filtered_physeq, 
                                  sample_type %in% c("Prepupae") & 
                                    is.na(AB_treatment))
adults_subset <- subset_samples(filtered_physeq, 
                                sample_type == "Adults" & 
                                  Env_exposure == "Free-flying")
larva_subset <- subset_samples(filtered_physeq, sample_type %in% c("Larvae"))
egg_subset <- subset_samples(filtered_physeq, sample_type %in% c("Eggs"))
frass_subset <- subset_samples(filtered_physeq, Description_of_sample %in% c("CocoonFrass", "Frass", "PinkCocoonFrass"))
frass_subset@sam_data$sample_type[which(frass_subset@sam_data$sample_type == "Nest_contents")] <- "Frass contents"
HB_subset <- subset_samples(filtered_physeq, sample_type %in% c("Honey_bee"))

# Compare the alpha diversity of tosti samples with honey bees ####
combined_subset1 <- merge_phyloseq(egg_subset, prepupal_subset, adults_subset, larva_subset, frass_subset, HB_subset)

alpha_diversity <- alpha(combined_subset1, index = "Shannon")
metadata <- meta(combined_subset1)
metadata$name <- rownames(metadata)
alpha_diversity$name <- rownames(alpha_diversity)
alpha_diversity_metadata <- merge(alpha_diversity, metadata, by = "name")

#Plot the alpha diversity
alpha_diversity_metadata %>%
  ggplot(aes(x = sample_type, y = diversity_shannon, colour = sample_type)) +
  geom_boxplot() +
  geom_point(size = 3)

### Analysis development stages ###
# Compare the beta diversity of the different development stages ####
combined_subset2 <- merge_phyloseq(prepupal_subset, adults_subset, larva_subset)
unweighted_unifrac <- ordinate(combined_subset2, method = "PCoA", distance = "unifrac", weighted=F)

# Create plot
ordination_plot <- plot_ordination(physeq = combined_subset2,
                                   ordination = unweighted_unifrac,
                                   color = "sample_type",
                                   axes = c(1, 2)) +
  theme_minimal() +
  geom_point(size = 3, alpha = 0.6)

# Add ellipses
fig1b <- ordination_plot +
  stat_ellipse(geom = "polygon", type="norm", alpha=0) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 16))

# Adonis statistic to determine whether communities are significantly different
#Refactor metadata
md_combined_subset <- as(sample_data(combined_subset2), "data.frame")
#Run adonis with 9999 permutations
uw_adonis <- adonis2(distance(combined_subset2, method="unifrac") ~ sample_type, 
                     data = md_combined_subset, 
                     permutations = 9999)
uw_adonis

# Combine the subsetted data ####
combined_subset <- merge_phyloseq(egg_subset, prepupal_subset, adults_subset, larva_subset, frass_subset)

# Obtain top 20 genera ####
ps_Genus <- tax_glom(combined_subset, taxrank = "Genus", NArm = FALSE)
top20Genus = names(sort(taxa_sums(ps_Genus), TRUE)[1:20])
taxtab20 = cbind(tax_table(ps_Genus), Genus_20 = NA)
taxtab20[top20Genus, "Genus_20"] <- as(tax_table(ps_Genus)
                                       [top20Genus, "Genus"], "character")
tax_table(ps_Genus) <- tax_table(taxtab20)
ps_Genus_ra <- transform_sample_counts(ps_Genus, function(x) 100 * x/sum(x))
df_Genus <- psmelt(ps_Genus_ra)
df_Genus <- arrange(df_Genus, sample_type)
df_Genus$Genus_20[is.na(df_Genus$Genus_20)] <- c("Other")

# Plot the relative abundance ####
# Custom order for sample types
custom_order <- c("Adults", "Eggs", "Larvae", "Prepupae", "Frass contents")
df_Genus$sample_type <- factor(df_Genus$sample_type, levels = custom_order)

# Custom order for cell ID from youngest to oldest
custom_Cell_ID_order <- c("J", "I", "H", "G", "F", "E", "D", "C", "B", "A")
df_Genus$Cell_ID <- factor(df_Genus$Cell_ID, levels = custom_Cell_ID_order)

(fig1a <- df_Genus %>%
  mutate(Genus_20 = reorder(Genus_20, -Abundance)) %>%
  ggplot(aes(x = sample_type, y = Abundance, fill = Genus_20)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = my_palette) +
  facet_nested(~ sample_type + Cell_ID + sampleid, 
               scales = "free", 
               space = "free",
               ) +
  coord_cartesian(expand = F) +
  labs(x = "sampleid", y = "Relative abundance") +
  theme( 
    axis.text.y = element_text(size = 16, face = 'bold'),
    axis.title.y = element_text(size = 16, face = 'bold'),
    axis.ticks.y = element_line(linewidth = 1),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "bottom",
    strip.text = element_textbox_simple(
      padding = margin(2, 0, 2, 0),
#      margin = margin(5, 5, 5, 5),
      size = 14,
      face = "bold",
      halign = 0.5,
      fill = "grey",
      box.color = "black",
      linewidth = 0.5,
      linetype = "solid"),
    panel.background = element_blank()
  ))


#use patchwork to stich a/b
fig1a + fig1b +
  plot_layout(widths = c(3, 1))

#Save plot
ggsave("figures/Figure1.png", height=9, width=15)
