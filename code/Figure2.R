#load required libraries
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

#test

#load phyloseq object
ps <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_4.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps

ps_acinetobacter <- subset_taxa(ps, Genus=="Acinetobacter")
df_temp <- as.data.frame(ps_acinetobacter@tax_table) 
df_temp$asvs <- row.names(ps_acinetobacter@tax_table)
ps_acinetobacter@tax_table <- tax_table(as.matrix(df_temp))
supp_adults2023 <- subset_samples(ps_acinetobacter, sample_type %in% c("Adults", "Negative_control"))
supp_adults2023

top20ASV = names(sort(taxa_sums(supp_adults2023), TRUE)[1:20])
taxtab20 = cbind(tax_table(supp_adults2023), ASV_20 = NA)
taxtab20[top20ASV, "ASV_20"] <- as(tax_table(supp_adults2023)
                                   [top20ASV, "asvs"], "character")

tax_table(supp_adults2023) <- tax_table(taxtab20)
ps_ASV_ra <- transform_sample_counts(supp_adults2023, function(x) 100 * x/sum(x))
df_ASV <- psmelt(ps_ASV_ra)
df_ASV <- arrange(df_ASV, sample_type)
df_ASV$ASV_20[is.na(df_ASV$ASV_20)] <- c("Other")


Acinetobacter_plot_supp <- df_ASV %>%
  filter(Abundance > 0) %>%
  ggplot(aes(x = sample_type, y = Abundance, fill = ASV_20)) +
geom_bar(stat = "identity") +
facet_nested(~ sample_type + Env_exposure + Sample, scales = "free", space = "free") +
labs(x = "sample_type", y = "Relative abundance") +
theme( 
axis.text.y = element_text(size=16, face = 'bold'),
axis.title.y = element_text(size=16, face = 'bold'),
axis.ticks.y = element_line(linewidth = 1),
axis.ticks.x = element_blank(),
axis.text.x = element_blank(),
axis.title.x = element_blank(),
legend.title = element_blank(),
legend.text = element_text(size = 18),
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
panel.background = element_blank())

#Load costum colour palette 
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

#rarefy at 1500
ps_rar <- rarefy_even_depth(ps, sample.size = 1500, rngseed = 1337)
ps_rar

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
#remove off-target reads
all_asvs <- taxa_names(ps_rar)
asvs_to_keep <- all_asvs[!(all_asvs %in% chloro_mito_decontam_asvs)]
ps_rar <- prune_taxa(asvs_to_keep, ps_rar)
ps_rar

#filter out 0 abundance reads
zero_abundance_samples <- sample_sums(ps_rar) == 0
filtered_physeq <- subset_samples(ps_rar, !zero_abundance_samples)
filtered_physeq

sort(sample_sums(filtered_physeq))

# Subset the data for samples_types 
food_through_time <- subset_samples(filtered_physeq, sample_type %in% c("Food") 
                                    & Nest_id %in% c("5", "13", "14", "27")
                                    & Cell_ID %in% c("J", "I", "H", "G", "F", "E", "D", "C", "B", "A"))
#obtain top 20 genera
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
mean(
  sample_sums(
    prune_taxa(top20Genus, ps_Genus_ra)
  )
)

# custom order for Cell_ID from youngest to oldest
custom_Cell_ID_order <- c("J", "I", "H", "G", "F", "E", "D", "C", "B", "A")
df_Genus$Cell_ID <- factor(df_Genus$Cell_ID, levels = custom_Cell_ID_order)
                                       
fig3 <- df_Genus %>%
  mutate(Genus_20 = reorder(Genus_20, -Abundance)) %>%
  ggplot(aes(x = sample_type, y = Abundance, fill = Genus_20)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = GenusPalette) +
  facet_nested(~ Nest_id + Cell_ID + sample_type, scales = "free", space = "free") +
  labs(x = "sample_type", y = "Relative abundance") +
  theme( 
    axis.text.y = element_text(size=16, face =, 'bold'),
    axis.title.y = element_text(size=16, face = 'bold'),
    axis.ticks.y = element_line(linewidth = 1),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
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
      linetype = "solid",
      family = "Calibri"),
    panel.background = element_blank()
  )


fig3

#alpha diversity if pollen provisions as they are consumed
alpha_diversity <- alpha(food_through_time, index = "Shannon")
metadata <- meta(food_through_time)
metadata$name <- rownames(metadata)
alpha_diversity$name <- rownames(alpha_diversity)
alpha_diversity_metadata <- merge(alpha_diversity, metadata, by = "name")

alpha_diversity_metadata %>%
  ggplot(aes(x = Cell_ID, y = diversity_shannon, colour = Cell_ID)) +
  geom_boxplot() +
  geom_point(size = 3)
