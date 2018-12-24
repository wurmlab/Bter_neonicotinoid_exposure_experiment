################################################################################################
##
## Author Joe Colgan http://wurmlab.com
## Copyright Queen Mary University London
## All rights reserved
##
################################################################################################
## The following scripts are used by plot_DEGs_heatmap.Rmd to generate a combined heatmap
## for both bumblebee castes.  

## Define function for plotting of heatmap:
caste_plotter <- function(caste){
        ## Step One: Load data:
        ## Load raw count estimates:
        raw_count_directory <- paste("input/count_estimates_",
                                     caste, "/", sep = "")
        load(paste(raw_count_directory, "txi_count_estimates.Rdata",
                   sep = ""))
        ## Load significant genes for plotting:
        load(paste("input/", pesticide, "_de_genes_output_",
                   caste, "/significant_genes.Rdata",
                   sep = ""))
        ## Step Two: Organise data for generating a heatmap:
        ## Normalise counts:
        normalised_counts_df <- t(scale(t(txi_counts$counts),
                                        center = TRUE,
                                        scale = TRUE))
        ## Select significant gene counts
        normalised_counts_df_sig_df <- as.data.frame(normalised_counts_df[de_genes_CLO, ])
        ## Step Three: Extract annotations for significant genes
        ## Using custom script, load gene annotations from NCBI for each
        ## significant gene:
        gene_descriptions <- description_from_LOCid(de_genes_CLO)
        ## Subtract the first three characters of the string and capitalize:
        pesticide_abbrev <- toupper(substr(pesticide,
                                           start = 1,
                                           stop = 3))
        control_abbrev <- toupper(substr(control,
                                         start = 1,
                                         stop = 3))
        ## Generate new sample names for plotting:
        sample_names <- c(paste(pesticide_abbrev, "_", 1:4, sep = ""),
                          paste(control_abbrev, "_", 1:4, sep = ""))
        ## Extract columns of interest:
        column_nums <- c(grep(pesticide_abbrev,
                              colnames(normalised_counts_df_sig_df)),
                         grep(control_abbrev,
                              colnames(normalised_counts_df_sig_df)))
        ## Subset columns of interest:
        normalised_counts_subset <- normalised_counts_df_sig_df[, column_nums]
        ## Rename columns:
        colnames(normalised_counts_subset) <- sample_names
        ## Add rownames with gene descriptions:
        rownames(normalised_counts_subset) <- gene_descriptions$id
        normalised_counts_subset$description <- gsub(":", "", gene_descriptions$both)
        ## Step Four: Cluster samples by euclidean distance for plotting:
        ## Calculate eucleidean distance between rows:
        dd <- hclust(dist(as.matrix(normalised_counts_subset[, 1:8])))
        ## Reorder based on cluster order:
        normalised_counts_ordered <- as.data.frame(normalised_counts_subset[dd$order, ])
        ## Reshape ("melt") for plotting with ggplot2:
        normalised_counts_ordered_melt <- melt(normalised_counts_ordered)
        ## Rename columns for plotting:
        colnames(normalised_counts_ordered_melt) <- c("gene_name",
                                                      "sample",
                                                      "normalised_counts")
        ## Step Five: Cluster samples by euclidean distance for plotting:
        ## Reorder "gene_name" for plotting:
        normalised_counts_ordered_melt$gene_name <-
        factor(normalised_counts_ordered_melt$gene_name,
               levels = rev(unique(c(as.character(unlist(normalised_counts_ordered_melt$gene_name))))))
        ## Create a treatment column to plot:
        normalised_counts_ordered_melt$treatment <- rep(c(pesticide, control),
                                                              each = (nrow(normalised_counts_ordered_melt) / 2))
        ## Depending on which caste, call a specific function to plot:
        if (caste == "workers"){
                plot <- heatmap_plot_worker(normalised_counts_ordered_melt)
                return(plot)
        } else {
                plot <- heatmap_plot_queen(normalised_counts_ordered_melt)
                return(plot)
        }
}

## Function for generating worker heatmap:
heatmap_plot_worker <- function(input){
                        ggplot(input, aes(sample, gene_name)) +
                                geom_tile(aes(fill = normalised_counts)) +
                                geom_vline(xintercept = c(0.5, 4.5, 8.5)) +
                                scale_fill_gradient2(low = "black",
                                                     mid = "white",
                                                     high = "red") +
                                ylab("") +
                                xlab("") +
                                theme(legend.title = element_text(size = 12,
                                                        face = "bold"),
                                      legend.text = element_text(size = 12,
                                                        face = "bold"),
                                      axis.text.x = element_blank(),
                                      axis.text.y = element_text(size = 12,
                                                        face = "plain",
                                                        colour = "black"),
                                      axis.title = element_text(size = 12,
                                                        face = "bold"),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.background = element_blank(),
                                      legend.position = "top") +
                                scale_y_discrete(position = "right") +
                                labs(fill = "Log fold change") +
                                theme(axis.text.y = element_text(face = c(rep("plain", 4),
                                                                            "bold.italic",
                                                                            rep("plain", 10),
                                                                            "bold.italic",
                                                                            rep("plain", 6),
                                                                            "bold",
                                                                            rep("plain", 31),
                                                                            "italic")))
}

## Function for plotting heatmap for queens:
heatmap_plot_queen <- function(input) {
                        ggplot(input, aes(sample, gene_name)) +
                                geom_tile(aes(fill = normalised_counts)) +
                                geom_vline(xintercept = c(0.5, 4.5, 8.5)) +
                                scale_fill_gradient2(low = "black",
                                                     mid = "white",
                                                     high = "red") +
                                ylab("") +
                                xlab("") +
                                theme(legend.title = element_text(size = 12,
                                                                face = "bold"),
                                legend.text = element_text(size = 12,
                                                        face = "bold"),
                                axis.text.x = element_text(angle = 45,
                                                            size = 10,
                                                            hjust = 1,
                                                            face = "bold"),
                                axis.text.y = element_text(size = 12,
                                                            face = "plain",
                                                            colour = "black"),
                                axis.title = element_text(size = 12,
                                                           face = "bold"),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.background = element_blank(),
                                legend.position = "top") +
                                scale_y_discrete(position = "right") +
                                labs(fill = "Log fold change") +
                                theme(axis.text.y = element_text(face = c(rep("plain", 2),
                                                                            "bold",
                                                                            rep("plain", 20))))
}
