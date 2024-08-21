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
library(cowplot)

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

# Identify zero abundance samples
zero_abundance_samples <- sample_sums(ps) == 0
ps <- subset_samples(ps, !zero_abundance_samples)
ps

# Remove primary isolates
ps <- subset_samples(ps, !sample_type %in% c("Primary_isolate")
                     & !AB_treatment %in% c("control", "antibiotic"))

#Collapse to Genus level and pull our relative abundance of top 20 common genuses.
ps_Genus5 <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)

if ("Family" %in% colnames(tax_table(ps_Genus5))) {
  
  # Replace NA values in the Genus column with "Family_unknown"
  tax_table(ps_Genus5)[is.na(tax_table(ps_Genus5)[, "Genus"]), "Genus"] <- 
    paste(as.character(tax_table(ps_Genus5)[is.na(tax_table(ps_Genus5)[, "Genus"]), "Family"]), "unclassified", sep = "_")
  
} else {
  # If the Family column doesn't exist, just replace NA
  tax_table(ps_Genus5)[is.na(tax_table(ps_Genus5)[, "Genus"]), "Genus"] <- "Unclassified"
}

# Replace "NA" (as a string) in the Genus column with "NA_unclassified"
tax_table(ps_Genus5)[tax_table(ps_Genus5)[, "Genus"] == "NA", "Genus"] <- "NA_unclassified"
tax_table(ps_Genus5)[tax_table(ps_Genus5)[, "Genus"] == "<NA>", "Genus"] <- "NA_unclassified"

top20Genus5 = names(sort(taxa_sums(ps_Genus5), TRUE)[1:19])
taxtab20 = cbind(tax_table(ps_Genus5), Genus_20 = NA)
taxtab20[top20Genus5, "Genus_20"] <- as(tax_table(ps_Genus5)
                                       [top20Genus5, "Genus"], "character")

tax_table(ps_Genus5) <- tax_table(taxtab20)
ps_Genus5_ra <- transform_sample_counts(ps_Genus5, function(x) 100 * x/sum(x))
df_Genus5 <- psmelt(ps_Genus5_ra)
df_Genus5 <- arrange(df_Genus5, sample_type)
df_Genus5$Genus_20[is.na(df_Genus5$Genus_20)] <- c("Other")

# View the result
print(unique(df_Genus5$Genus_20))

# % of reads that make up the top 20 genera ####
mean(sample_sums(prune_taxa(top20Genus5, ps_Genus5_ra))
)


custom_order <- c("Mitochondria", 
                  "Chloroplast", 
                  "Acinetobacter",
                  "Pseudomonas", 
                  "Tyzzerella",
                  "Lactobacillus",
                  "Sodalis", 
                  "Sphingobium", 
                  "Bacillus",
                  "Enterobacteriaceae_unclassified", 
                  "Erwiniaceae_unclassified", 
                  "Escherichia-Shigella",  
                  "Commensalibacter",
                  "Gilliamella", 
                  "Snodgrassella", 
                  "Bifidobacterium",
                  "Bartonella", 
                  "Brevibacterium",
                  "NA_unclassified", 
                  "Other")
                  
df_Genus5$Genus_20 <- factor(df_Genus5$Genus_20, levels = custom_order)

df_Genus5 <- df_Genus5 %>%
  dplyr::mutate(sample_type2 = dplyr::if_else(sample_type == "Adults",
                                              stringr::str_c("Tosti_", Env_exposure),
                                              sample_type), .after = sample_type)

order <- c("Negative_control", "Tosti_Free-flying", "Tosti_Nest_only", "Tosti_None", "Food", "Honey_bee")
df_Genus5$sample_type2 <- factor(df_Genus5$sample_type2, levels=order)


(absolutePlot <- df_Genus5 %>%
    dplyr::filter(complete.cases(logDNA)) %>%
    dplyr::distinct(sampleid, Genus_20, .keep_all = TRUE) %>%
    dplyr::mutate(logDNA_weighted = logDNA*(Abundance/100)) %>%
    ggplot(aes(x = sampleid)) +
    geom_col(width = 1, aes(fill = Genus_20, y = logDNA_weighted)) +
    scale_fill_manual(values = my_palette) +
    facet_nested(~ sample_type2, scales = "free", space = "free") +
    labs(x = "sampleid", y = "Absolute abundance (ng)") +
    theme( 
      axis.text.y = element_text(size=14, face = "bold"),
      axis.title.y = element_text(size=14, face = "bold"),
      axis.ticks.y = element_line(linewidth = 1),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = "none",
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
      panel.background = element_blank()
    ))

(relativeBar <- df_Genus5 %>%
    dplyr::filter(complete.cases(logDNA)) %>%
    ggplot(aes(x = sampleid, y = Abundance, fill = Genus_20)) +
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values = my_palette) +
    facet_nested(~ sample_type2 , scales = "free", space = "free") +
    labs(x = "sampleid", y = "Relative abundance") +
    theme( 
      axis.text.y = element_text(size=14, face = "bold"),
      axis.title.y = element_text(size=14, face = "bold"),
      axis.ticks.y = element_line(linewidth = 1),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_blank(),
      panel.background = element_blank()
    ))



(legendPlot <- df_Genus5 %>%
    dplyr::filter(complete.cases(logDNA)) %>%
    dplyr::distinct(sampleid, Genus_20, .keep_all = TRUE) %>%
    mutate(Genus_20 = reorder(Genus_20, -Abundance)) %>%
    dplyr::mutate(logDNA_weighted = logDNA*(Abundance/100)) %>%
    ggplot(aes(x = sampleid)) +
    geom_col(width = 1, aes(fill = Genus_20, y = logDNA_weighted)) +
    scale_fill_manual(values = my_palette) +
    facet_nested(~ sample_type2 , scales = "free", space = "free") +
    #facet_nested(~ sample_type , scales = "free", space = "free") +
    labs(x = "sampleid", y = "Relative abundance") +
    theme( 
      axis.text.y = element_text(size=12, face =, 'bold'),
      axis.title.y = element_text(size=12, face = 'bold'),
      axis.ticks.y = element_line(linewidth = 1),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
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
      panel.background = element_blank()
    ))


(abundancesplot_nolegend <- plot_grid(absolutePlot, relativeBar, ncol = 1, align = "hv"))

(Abundances_plot <- cowplot::plot_grid(abundancesplot_nolegend,
                                       cowplot::get_legend(legendPlot),
                                       cols = 1,
                                       rel_heights = c(1.3, 0.2)
                                       ))

ggsave("figures/Figure5.png", height=10, width=15)

