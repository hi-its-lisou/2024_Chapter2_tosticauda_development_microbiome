# Load the required libraries ####
library(phyloseq)
library(tidyverse)
library(qiime2R)
library(decontam)
library(microbiome)


# Load phyloseq object ####
ps_raw <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_5.tsv")
ps_raw@sam_data$sampleid = rownames(ps_raw@sam_data)
ps_raw

# Total reads in the raw (unfiltered) dataset
total_reads_raw <- sum(sample_sums(ps_raw))
sort(sample_sums(ps_raw))

# Filter out samples with low abundance reads ####
high_read_samples <- sample_sums(ps_raw) >= 300 
ps <- subset_samples(ps_raw, high_read_samples)
ps

# Remove contaminants, chloroplasts, and mitochondria ####
contam <- read_delim("input_files/chloro_mito_decontam_asvs.txt", 
                     delim = "\n", 
                     col_names = "asv")
chloro_mito_decontam_asvs <- contam$asv
all_asvs <- taxa_names(ps)
asvs_to_keep <- all_asvs[!(all_asvs %in% chloro_mito_decontam_asvs)]
ps_filtered <- prune_taxa(asvs_to_keep, ps)
ps_filtered

sort(sample_sums(ps_filtered))

# Total reads in the filtered dataset
total_reads_filtered <- sum(sample_sums(ps_filtered))

# Reads per sample in the raw dataset
reads_per_sample_raw <- sample_sums(ps_raw)
mean_reads_raw <- mean(reads_per_sample_raw)

# Reads per sample in the filtered dataset
reads_per_sample_filtered <- sample_sums(ps_filtered)
mean_reads_filtered <- mean(reads_per_sample_filtered)

# Summary
cat("Sequencing Summary:\n")
cat("Total reads generated (raw):", total_reads_raw, "\n")
cat("Total reads remaining after filtering:", total_reads_filtered, "\n")
cat("Mean reads per sample (raw):", round(mean_reads_raw), "\n")
cat("Mean reads per sample (filtered):", round(mean_reads_filtered), "\n")


# Extract the OTU table and taxonomy table from the phyloseq object
otu_table_df <- as.data.frame(as(otu_table(ps_raw), "matrix"))
tax_table_df <- as.data.frame(as(tax_table(ps_raw), "matrix"))

# Ensure the taxa names are consistent in both tables
otu_table_df$TaxaID <- rownames(otu_table_df)
tax_table_df$TaxaID <- rownames(tax_table_df)

# Merge OTU and taxonomy tables based on TaxaID
merged_data <- merge(otu_table_df, tax_table_df, by = "TaxaID")

# Identify chloroplast and mitochondrial ASVs/OTUs in the taxonomy
chloroplast_mitochondria_taxa <- merged_data %>%
  filter(grepl("chloroplast", Genus, ignore.case = TRUE) | 
           grepl("mitochondria", Family, ignore.case = TRUE))

# Select only the OTU count columns (excluding non-numeric columns like taxonomic information)
# Assuming taxonomic columns are at the end of the data frame after merging
otu_columns <- select(chloroplast_mitochondria_taxa, where(is.numeric))

# Calculate total reads for chloroplast and mitochondria ASVs/OTUs
chloroplast_mitochondria_reads <- rowSums(otu_columns)

# Calculate total reads in the entire dataset
total_reads <- sum(rowSums(otu_table(ps)))

# Calculate the proportion of reads that are chloroplasts or mitochondria
proportion_chloroplast_mitochondria <- sum(chloroplast_mitochondria_reads) / total_reads

# Print the results
cat("Total reads:", total_reads, "\n")
cat("Chloroplast/Mitochondria reads:", sum(chloroplast_mitochondria_reads), "\n")
cat("Proportion of reads that are chloroplast or mitochondria:", round(proportion_chloroplast_mitochondria * 100, 2), "%\n")





# Extract the OTU table from the phyloseq object as a matrix
otu_table_matrix <- as(otu_table(filtered_physeq), "matrix") 

if (taxa_are_rows(otu_table(filtered_physeq))) {
  otu_table_matrix <- t(otu_table_matrix)
}
otu_table_df <- as.data.frame(otu_table_matrix)
otu_table_df <- tibble::rownames_to_column(otu_table_df, var = "SampleID")
write.csv(otu_table_df, file = "otu_table_filtered.csv", row.names = FALSE)



# Extract the taxonomy table from the phyloseq object
taxonomy_table_matrix <- as(tax_table(ps), "matrix") 

# Convert to a data frame for better formatting
taxonomy_table_df <- as.data.frame(taxonomy_table_matrix)

# Add ASV/OTU names as a column for identification
taxonomy_table_df <- tibble::rownames_to_column(taxonomy_table_df, var = "ASV")

# Save the taxonomy table to a CSV file
write.csv(taxonomy_table_df, file = "taxonomy_table.csv", row.names = FALSE)

