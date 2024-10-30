# Load the required libraries ####
library(phyloseq)
library(tidyverse)
library(qiime2R)
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
colours_df <- read_csv("input_files/colour_list.csv")
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

# Filter out samples with low abundance reads ####
low_abundance_samples <- sample_sums(ps_rar) <= 100
filtered_physeq <- subset_samples(ps_rar, !low_abundance_samples)
filtered_physeq

sort(sample_sums(filtered_physeq))

# Subset the data for samples_types 
prepupal_subset <- subset_samples(filtered_physeq, 
                                  sample_type %in% c("Prepupae") & 
                                    is.na(AB_treatment))
adults_subset <- subset_samples(filtered_physeq, 
                                sample_type == "Adults" & 
                                  Env_exposure == "Free-flying")
larva_subset <- subset_samples(filtered_physeq, sample_type %in% c("Larvae")
                               & Year %in% c("2022"))
frass_subset <- subset_samples(filtered_physeq, Description_of_sample %in% c("CocoonFrass", "Frass", "PinkCocoonFrass"))
frass_subset@sam_data$sample_type[which(frass_subset@sam_data$sample_type == "Nest_contents")] <- "Frass contents"
HB_subset <- subset_samples(filtered_physeq, sample_type %in% c("Honey_bee"))

### Analysis development stages ###
# Compare the beta diversity of the different development stages ####
development <- merge_phyloseq(prepupal_subset, adults_subset, larva_subset)
unweighted_unifrac <- ordinate(development, method = "PCoA", distance = "unifrac", weighted=F)

colours <- c("#f87970", "green", "#4B0092")

(Figure1B <- plot_ordination(physeq = development, 
                             ordination = unweighted_unifrac, 
                             color = "sample_type",
                             axes = c(1, 2)) +
    geom_point(size = 5, alpha = 0.6) +
    stat_ellipse(geom = "polygon", type="norm", alpha=0) +
    theme_minimal() +
    scale_color_manual(values = colours) +
    theme (axis.text.y = element_text(size=14, face = 'bold'),
           axis.title.y = element_text(size=14, face = 'bold'),
           axis.text.x = element_text(size=14, face = 'bold'),
           axis.title.x = element_text(size=14, face = 'bold'),
           legend.text = element_text(size = 14),
           legend.position = "top", 
           legend.title = element_blank()))

# Adonis statistic to determine whether communities are significantly different
#Refactor metadata
md_combined_subset <- as(sample_data(development), "data.frame")
#Run adonis with 9999 permutations
uw_adonis <- adonis2(distance(development, method="unifrac") ~ sample_type, 
                     data = md_combined_subset, 
                     permutations = 9999)
uw_adonis

# Combine the subsetted data ####
combined_subset <- merge_phyloseq(prepupal_subset, adults_subset, larva_subset, frass_subset)

# Determine the number of reads per sample ####
sort(sample_sums(combined_subset))

# Obtain top 20 genera ####
ps_Genus1 <- tax_glom(combined_subset, taxrank = "Genus", NArm = FALSE)

if ("Family" %in% colnames(tax_table(ps_Genus1))) {
  
  # Replace NA values in the Genus column with "Family_unknown"
  tax_table(ps_Genus1)[is.na(tax_table(ps_Genus1)[, "Genus"]), "Genus"] <- 
    paste(as.character(tax_table(ps_Genus1)[is.na(tax_table(ps_Genus1)[, "Genus"]), "Family"]), "unclassified", sep = "_")
  
} else {
  # If the Family column doesn't exist, just replace NA with "Unknown"
  tax_table(ps_Genus1)[is.na(tax_table(ps_Genus1)[, "Genus"]), "Genus"] <- "Unclassified"
}

top20Genus1 = names(sort(taxa_sums(ps_Genus1), TRUE)[1:19])
taxtab20 = cbind(tax_table(ps_Genus1), Genus_20 = NA)
taxtab20[top20Genus1, "Genus_20"] <- as(tax_table(ps_Genus1)
                                       [top20Genus1, "Genus"], "character")
tax_table(ps_Genus1) <- tax_table(taxtab20)
ps_Genus1_ra <- transform_sample_counts(ps_Genus1, function(x) 100 * x/sum(x))
df_Genus1 <- psmelt(ps_Genus1_ra)
df_Genus1 <- arrange(df_Genus1, sample_type)
df_Genus1$Genus_20[is.na(df_Genus1$Genus_20)] <- c("Other")

# View the result
print(unique(df_Genus1$Genus_20))

# % of reads that make up the top 20 genera
mean(sample_sums(prune_taxa(top20Genus1, ps_Genus1_ra)))


# Plot the relative abundance ####
# Custom order for sample types
custom_order <- c("Adults", "Larvae", "Prepupae", "Frass contents")
df_Genus1$sample_type <- factor(df_Genus1$sample_type, levels = custom_order)

# Custom order for cell ID from youngest to oldest
custom_Cell_ID_order <- c("J", "I", "H", "G", "F", "E", "D", "C", "B", "A")
df_Genus1$Cell_ID <- factor(df_Genus1$Cell_ID, levels = custom_Cell_ID_order)

(Figure1A <- df_Genus1 %>%
  mutate(Genus_20 = reorder(Genus_20, -Abundance)) %>%
  ggplot(aes(x = sample_type, y = Abundance, fill = Genus_20)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = my_palette) +
  facet_nested(~ sample_type + Year + Cell_ID + sampleid, 
               scales = "free", 
               space = "free") +
  labs(x = "sampleid", y = "Relative abundance") +
  theme( 
    axis.text.y = element_text(size = 14, face = 'bold'),
    axis.title.y = element_text(size = 14, face = 'bold'),
    axis.ticks.y = element_line(linewidth = 1),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = "bottom",
    strip.text = element_textbox_simple(
      padding = margin(2, 0, 2, 0),
      size = 12,
      face = "bold",
      halign = 0.5,
      fill = "grey",
      box.color = "black",
      linewidth = 0.5,
      linetype = "solid"),
    panel.background = element_blank()
  ))


#use patchwork to stich a/b
Figure1 <- cowplot::plot_grid(Figure1A,
                              Figure1B,
                              ncol = 2,
                              rel_heights = c(1,0.6),
                              rel_widths = c(1, 0.5))
Figure1

#Save plot
ggsave("figures/OctFigure1.png", height=10, width=18)

