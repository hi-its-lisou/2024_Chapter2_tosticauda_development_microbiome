#####################################################################################################
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

ps <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_4.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps #1089 taxa and 208 samples

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
GenusList = unique(tax_table(ps)[,"Genus"])
GenusPalette = getPalette(length(GenusList))
names(GenusPalette) = GenusList
GenusPalette["Other"] <- "grey" 

new_color <- "#3D550C"
taxon_to_change <- "Erwinia"
GenusPalette[taxon_to_change] <- new_color
new_color2 <- "#81B622"
taxon_to_change2 <- "Arsenophonus"
GenusPalette[taxon_to_change2] <- new_color2
new_color3 <- "#ECF87F"
taxon_to_change3 <- "Lactobacillus"
GenusPalette[taxon_to_change3] <- new_color3
new_color4 <- "#59981A"
taxon_to_change4 <- "Dolosigranulum"
GenusPalette[taxon_to_change4] <- new_color4
new_color5 <- "#F51720"
taxon_to_change5 <- "Sodalis"
GenusPalette[taxon_to_change5] <- new_color5
new_color6 <- "#FA26A0"
taxon_to_change6 <- "Acinetobacter"
GenusPalette[taxon_to_change6] <- new_color6
new_color7 <- "#F8D210"
taxon_to_change7 <- "Tyzerella"
GenusPalette[taxon_to_change7] <- new_color7
new_color8 <- "#2FF3E0"
taxon_to_change8 <- "Sphingomonas"
GenusPalette[taxon_to_change8] <- new_color8
new_color9 <- "#B1D4E0"
taxon_to_change9 <- "Skermanella"
GenusPalette[taxon_to_change9] <- new_color9
new_color10 <- "#2E8BC0"
taxon_to_change10 <- "Gilliamella"
GenusPalette[taxon_to_change10] <- new_color10
new_color11 <- "#0C2D48"
taxon_to_change11 <- "Aquabacterium"
GenusPalette[taxon_to_change11] <- new_color11
new_color12 <- "#145DA0"
taxon_to_change12 <- "Thauera"
GenusPalette[taxon_to_change12] <- new_color12
new_color14 <- "#FFD4DB"
taxon_to_change14 <- "Massilia"
GenusPalette[taxon_to_change14] <- new_color14
new_color15 <- "#D3B5E5"
taxon_to_change15 <- "Modestobacter"
GenusPalette[taxon_to_change15] <- new_color15
new_color17 <- "#D8A7B1"
taxon_to_change17 <- "Asaia"
GenusPalette[taxon_to_change17] <- new_color17
new_color20 <- "#EF7C8E"
taxon_to_change20 <- "Carnimonas"
GenusPalette[taxon_to_change20] <- new_color20
new_color21 <- "#0023F5"
taxon_to_change21 <- "Bacillus"
GenusPalette[taxon_to_change21] <- new_color21
new_color22 <- "#A020F0"
taxon_to_change22 <- "Enterococcus"
GenusPalette[taxon_to_change22] <- new_color22
new_colorX <- "#000000"
taxon_to_changeX <- "Mitochondria"
GenusPalette[taxon_to_changeX] <- new_colorX
new_colorY <- "#008000"
taxon_to_changeY <- "Chloroplast"
GenusPalette[taxon_to_changeY] <- new_colorY

#rarefy at 1500
ps_rar <- rarefy_even_depth(ps, sample.size = 1500, rngseed = 1337)
ps_rar

#import contaminant asvs
contam <- read_delim("input_files/chloro_mito_decontam_asvs.txt", 
                     delim = "\n", 
                     col_names = "asv")

chloro_mito_decontam_asvs <- contam$asv

#vector for chloroplast, mitochondrial, and contaminant asvs
chloro_mito_decontam_asvs

#remove off-target reads
all_asvs <- taxa_names(ps_rar)
asvs_to_keep <- all_asvs[!(all_asvs %in% chloro_mito_decontam_asvs)]
ps_rar <- prune_taxa(asvs_to_keep, ps_rar)
ps_rar

#subset prepupae from antibiotic experiment only
AB_exp <- subset_samples(ps_rar, AB_treatment %in% c("control", "antibiotic"))
AB_exp #821 taxa and 39 samples

#remove 0 abundsance samples
zero_abundance_samples <- sample_sums(AB_exp) == 0
filtered_physeq <- subset_samples(AB_exp, !zero_abundance_samples)
filtered_physeq

#shannon diversity box plot
alpha_diversity <- alpha(filtered_physeq, index = "Shannon")

metadata <- meta(filtered_physeq)
metadata$name <- rownames(metadata)
alpha_diversity$name <- rownames(alpha_diversity)
alpha_diversity_metadata <- merge(alpha_diversity, metadata, by = "name")

fig5 <- alpha_diversity_metadata %>%
  ggplot(aes(x = AB_treatment, y = diversity_shannon, colour = AB_treatment)) +
  geom_boxplot() +
  geom_point(size = 3)

fig5

# relative abundance plot
#top 20 genera
ps_Genus <- tax_glom(filtered_physeq, taxrank = "Genus", NArm = FALSE)
top20Genus = names(sort(taxa_sums(ps_Genus), TRUE)[1:20])
taxtab20 = cbind(tax_table(ps_Genus), Genus_20 = NA)
taxtab20[top20Genus, "Genus_20"] <- as(tax_table(ps_Genus)
                                       [top20Genus, "Genus"], "character")
tax_table(ps_Genus) <- tax_table(taxtab20)
ps_Genus_ra <- transform_sample_counts(ps_Genus, function(x) 100 * x/sum(x))
df_Genus <- psmelt(ps_Genus_ra)
df_Genus <- arrange(df_Genus, sample_type)
df_Genus$Genus_20[is.na(df_Genus$Genus_20)] <- c("Other")

mean(
  sample_sums(
    prune_taxa(top20Genus, ps_Genus_ra)
  )
)


fig4 <- df_Genus %>%
  mutate(Genus_20 = reorder(Genus_20, -Abundance)) %>%
  ggplot(aes(x = sample_type, y = Abundance, fill = Genus_20)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = GenusPalette) +
  facet_nested(~ AB_treatment + sampleid + Year, scales = "free", space = "free") +
  labs(x = "sample_type", y = "Relative abundance") +
  theme( 
    axis.text.y = element_text(size=16, face = 'bold'),
    axis.title.y = element_text(size=16, face = 'bold'),
    axis.ticks.y = element_line(linewidth = 1),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
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
  )

fig4


ggsave("Figure3.png", height=15, width=20)



### Acinetobacter plot - Supplementary Figure 3 ###
# Phyloseq object with only Acinetobacter ASVs
ps_acinetobacter <- subset_taxa(ps, Genus=="Acinetobacter")
df_temp <- as.data.frame(ps_acinetobacter@tax_table) 
df_temp$asvs <- row.names(ps_acinetobacter@tax_table)
ps_acinetobacter@tax_table <- tax_table(as.matrix(df_temp))

# Subset samples
supp_AB_exp1 <- subset_samples(ps_acinetobacter, AB_treatment %in% c("control", "antibiotic"))
supp_AB_exp2 <- subset_samples(ps_acinetobacter, sample_type %in% c("Negative_control"))
supp_AB_exp <- merge_phyloseq(supp_AB_exp1, supp_AB_exp2)
supp_AB_exp

sample_sums(supp_AB_exp)

# Sort top 20 ASVs
top20ASV = names(sort(taxa_sums(supp_AB_exp), TRUE)[1:20])
taxtab20 = cbind(tax_table(supp_AB_exp), ASV_20 = NA)
taxtab20[top20ASV, "ASV_20"] <- as(tax_table(supp_AB_exp)
                                   [top20ASV, "asvs"], "character")

tax_table(supp_AB_exp) <- tax_table(taxtab20)
ps_ASV_ra <- transform_sample_counts(supp_AB_exp, function(x) 100 * x/sum(x))
df_ASV <- psmelt(ps_ASV_ra)
df_ASV <- arrange(df_ASV, sample_type)
df_ASV$ASV_20[is.na(df_ASV$ASV_20)] <- c("Other")

# Colour palette to distinguish true Acinetobacter from contam
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
ASVList = unique(tax_table(ps_acinetobacter)[,"asvs"])
ASVPalette = getPalette(length(ASVList))
names(ASVPalette) = ASVList

new_color_Acinetobacter <- "#FA26A0"
Acinetobacter_to_change <- "c624f2e4228eea7296b2a77e2d4b7e50"
ASVPalette[Acinetobacter_to_change] <- new_color_Acinetobacter

#Plot relative abundance of Acinetobacter ASVs
Acinetobacter_plot_supp3 <- df_ASV %>%
  filter(Abundance > 0) %>%
  ggplot(aes(x = sample_type, y = Abundance, fill = ASV_20)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = ASVPalette) +
  facet_nested(~ sample_type + AB_treatment + Sample, scales = "free", space = "free") +
  labs(x = "sample_type", y = "Relative abundance") +
  theme( 
    axis.text.y = element_text(size=16, face = 'bold'),
    axis.title.y = element_text(size=16, face = 'bold'),
    axis.ticks.y = element_line(linewidth = 1),
    axis.ticks.x = element_blank(),
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
      linetype = "solid",),
    panel.background = element_blank())


Acinetobacter_plot_supp3
ggsave("FigureS3.png", height=15, width=25)
