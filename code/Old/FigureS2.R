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

# Load phyloseq object ####
ps <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_4.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps

# Subset the data for negative controls and adults in the acquisition experiment ####
supp_adults2023 <- subset_samples(ps, sample_type %in% c("Adults", "Negative_control"))
supp_adults2023

# Obtain top 20 genera ####
ps_Genus <- tax_glom(supp_adults2023, taxrank = "Genus", NArm = FALSE)

# Ensure the taxonomy table has a Family column
if ("Family" %in% colnames(tax_table(ps_Genus))) {
  
  # Replace NA values in the Genus column with "Family_unknown"
  tax_table(ps_Genus)[is.na(tax_table(ps_Genus)[, "Genus"]), "Genus"] <- 
    paste(as.character(tax_table(ps_Genus)[is.na(tax_table(ps_Genus)[, "Genus"]), "Family"]), "unclassified", sep = "_")
  
} else {
  # If the Family column doesn't exist, just replace NA with "Unknown"
  tax_table(ps_Genus)[is.na(tax_table(ps_Genus)[, "Genus"]), "Genus"] <- "Unclassified"
}

top20Genus = names(sort(taxa_sums(ps_Genus), TRUE)[1:19])
taxtab20 = cbind(tax_table(ps_Genus), Genus_20 = NA)
taxtab20[top20Genus, "Genus_20"] <- as(tax_table(ps_Genus)
                                       [top20Genus, "Genus"], "character")
tax_table(ps_Genus) <- tax_table(taxtab20)
ps_Genus_ra <- transform_sample_counts(ps_Genus, function(x) 100 * x/sum(x))
df_Genus <- psmelt(ps_Genus_ra)
df_Genus <- arrange(df_Genus, sample_type)
df_Genus$Genus_20[is.na(df_Genus$Genus_20)] <- c("Other")

# Load custom colour palette ####
colours_df <- read_csv("input_files/colour_list.csv")
my_palette <- colours_df$colours
names(my_palette) <- colours_df$genera
my_palette

# Plot the relative abundance ####
(FigureS1 <- df_Genus %>%
  mutate(Genus_20 = reorder(Genus_20, -Abundance)) %>%
  ggplot(aes(x = sample_type, y = Abundance, fill = Genus_20)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = my_palette) +
  facet_nested(~ sample_type + Env_exposure + sampleid, scales = "free", space = "free") +
  labs(x = "sample_type", y = "Relative abundance") +
  theme( 
    axis.text.y = element_text(size=16, face = 'bold'),
    axis.title.y = element_text(size=16, face = 'bold'),
    axis.ticks.y = element_line(linewidth = 1),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_textbox_simple(
      padding = margin(5, 0, 5, 0),
      margin = margin(5, 5, 5, 5),
      size = 10,
      face = "bold",
      halign = 0.5,
      fill = "white",
      box.color = "grey",
      linewidth = 1.5,
      linetype = "solid",),
    panel.background = element_blank()
  ))

#Save plot
ggsave("figures/S2.png", height=10, width=15)

