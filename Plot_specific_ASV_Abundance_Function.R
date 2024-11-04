### Function to plot the abundance of a given ASV across all samples.
### This is useful for identifying ASVs that show evidence of cross-sample contamination, i.e. high abundance
### in some samples, but very low abundance in others. I have based this on phyloseq's "plot_bar" function.

plot.specific.ASV <- function(physeq, asv, facet_grid=NULL) {
  
  #Sort the phyloseq object by most abundant features, then select only the feature of interest based on its
  #abundance. E.g. 1 = most abundant, 5 = 5th most abundant
  psASV = prune_taxa(taxa = asv, x = physeq)
  
  #Calculate the total abundance of selected ASV in dataset
  total_abundance <- sum(sample_sums(physeq))
  asv_abundance <- sum(sample_sums(psASV))
  
  #Melt the phyloseq data, as ggplot likes long format
  mdf = psmelt(psASV)
  
  #Label is used as input data for geom_text
  plot = ggplot(mdf, aes(x = Sample, y = Abundance, label = Abundance)) +
    
    #Columns are cleaner than bars, and we're only looking at 1 ASV
    geom_col(fill="BLACK") +
    
    #Add text to the top of each column -- easier to see for low feature count samples
    geom_text(position = position_dodge(width = 1), vjust = -0.5, size = 3) +
    
    #Make pretty
    ggtitle(paste("ASV:", asv, "       ",
                  "Genus:", as.data.frame(psASV@tax_table) %>% slice(1) %>% select(Genus), "\n",
                  "Total relative abundance of ASV:", scales::percent(asv_abundance/total_abundance, accuracy = 0.1))) +
    
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  #If the user provides a metadata category to facet data, use it here. Note that the paste(". ~", facet_grid)
  #was essential, otherwise facet_grid() does not read the variable!
# Modify the faceting part
if( !is.null(facet_grid) ){
    # Modify the formula for facet_grid to include both variables
    plot <- plot + facet_grid(as.formula(paste(facet_grid, "~ .")), scales = "free_x", space = "free_x")
}

  
  return(plot)  
  
}


##############################################################################
