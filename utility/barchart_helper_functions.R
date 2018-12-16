################################################################################################
##
## Author Joe Colgam
## Copyright Queen Mary University London
## All rights reserved
##
################################################################################################
## The following scripts are used by plot_go_barchart.Rmd to generate a combined barhart.


## Define function for plotting of heatmap:
caste_barchart_plotter <- function(caste){
        
                                ### Read in biological processes GO terms:
                                pesticide_caste_BP <- read.table(file=paste("input/", caste, "/input_", pesticide_abbrev, "/BP_sig.txt", sep=""), 
                                                                           header = TRUE)
                                ### Read in molecular function GO terms:
                                pesticide_caste_MF <- read.table(file=paste("input/", caste, "/input_", pesticide_abbrev, "/MF_sig.txt", sep=""), 
                                                                           header = TRUE)
                                
                                ### Read in cellular component GO terms:
                                pesticide_caste_CC <- read.table(file=paste("input/", caste, "/input_", pesticide_abbrev, "/CC_sig.txt", sep=""), 
                                                                           header = TRUE)
                                
                                ## Combine GO terms:
                                pesticide_caste_combined <- rbind(pesticide_caste_BP, 
                                                                  pesticide_caste_MF,
                                                                  pesticide_caste_CC)
                                
                                pesticide_caste_combined$category <- factor(pesticide_caste_combined$category)
                                
                                # Set order of GO categories to plot
                                levels(pesticide_caste_combined$category) <- c("BP", "MF", "CC")
                                
                                ## Round p-values for plot:
                                pesticide_caste_combined$weight_ks <- round(pesticide_caste_combined$weight_ks, digits = 3)
                                
                                # Log transform p values for plotting:
                                pesticide_caste_combined$weight_ks <- -log(pesticide_caste_combined$weight_ks)
                                
                                # Reorder terms for plotting:
                                pesticide_caste_combined$Term <- factor(pesticide_caste_combined$Term, levels = pesticide_caste_combined$Term[order(pesticide_caste_combined$category, pesticide_caste_combined$weight_ks)])
                                
                                ## For plotting, update GO term names to include total number of annotated terms:
                                pesticide_caste_combined$updated_terms <- paste(pesticide_caste_combined$Term, " ", "(", pesticide_caste_combined$Annotated, ")", sep="")
                                
                                ## Remove underscore:
                                pesticide_caste_combined$updated_terms <- gsub("_", " ", pesticide_caste_combined$updated_terms)
                                
                                # Depending on which caste, call a specific function to plot:
                                if (caste=="worker"){
                                        plot <- worker_barchart_plot(pesticide_caste_combined)
                                        return(plot)
                                } else {
                                        plot <- queen_barchart_plot(pesticide_caste_combined)
                                        return(plot)
                                }
}

# Plot for workers:
worker_barchart_plot<- function(input){
                        pesticide_worker_plot <- ggbarplot(input, x = "updated_terms", y = "weight_ks",
                                                      position = position_dodge(0.1),
                                                      fill = "category",# change fill color by mpg_level
                                                      color = NULL, #"white"            # Set bar border colors to white
                                                      palette = "jco",# jco journal color palett. see ?ggpar
                                                      sort.val = "asc",          # Sort the value in descending order
                                                      sort.by.groups = TRUE,     # Don't sort inside each group
                                                      ylab = "-log10(p)",
                                                      xlab = "Gene ontology term",
                                                      legend.title = "Gene ontology",
                                                      lab.col = "black",
                                                      lab.size = 4,
                                                      lab.vjust = 0.5,
                                                      lab.hjust = 1,
                                                      legend = "top",
                                                      rotate = TRUE,
                                                      ggtheme = theme_minimal())
        
                        pesticide_worker_plot <- pesticide_worker_plot +
                                                        scale_y_continuous(expand = c(0, 0)) +
                                                        theme(axis.text=element_text(size=15),
                                                              axis.title.x = element_text(size=15,face="bold"),
                                                              axis.title.y = element_text(size=15,face="bold"),
                                                              axis.text.y = element_text(size=12, face="bold"),
                                                              axis.text.x = element_text(size=12),
                                                              legend.position="none") +
                                                        expand_limits(y = 10) +
                                                        geom_hline(yintercept = 1.301, linetype="dashed", colour="black")
                        
                        ## Update colours for plotting
                        pesticide_worker_plot<- pesticide_worker_plot + 
                                                                scale_fill_manual(values = c("orange", 
                                                                                             "light blue", 
                                                                                             "light grey"))
}


## Plot queen:
queen_barchart_plot<- function(input){
                        pesticide_queen_plot <- ggbarplot(input, x = "updated_terms", y = "weight_ks",
                                                      position = position_dodge(0.1),
                                                      fill = "category",# change fill color by mpg_level
                                                      color = NULL, #"white"            # Set bar border colors to white
                                                      palette = "jco",# jco journal color palett. see ?ggpar
                                                      sort.val = "asc",          # Sort the value in descending order
                                                      sort.by.groups = TRUE,     # Don't sort inside each group
                                                      ylab = "-log10(p)",
                                                      xlab = "Gene ontology term",
                                                      legend.title = "Gene ontology",
                                                      lab.col = "black",
                                                      lab.size = 4,
                                                      lab.vjust = 0.5,
                                                      lab.hjust = 1,
                                                      legend = "top",
                                                      rotate = TRUE,
                                                      ggtheme = theme_minimal())

                        pesticide_queen_plot <- pesticide_queen_plot +
                                                        scale_y_continuous(expand = c(0, 0)) +
                                                        theme(axis.text=element_text(size=15),
                                                              axis.title.x = element_text(size=15,face="bold"),
                                                              axis.title.y = element_text(size=15,face="bold"),
                                                              axis.text.y = element_text(size=12, face="bold"),
                                                              axis.text.x = element_text(size=12),
                                                              legend.position="none") +
                                                        expand_limits(y = 10) +
                                                        geom_hline(yintercept = 1.301, linetype="dashed", colour="black")

                        ## Update colours for plotting
                        pesticide_queen_plot <- pesticide_queen_plot + 
                                                         scale_fill_manual(values = c("orange", 
                                                                                     "light blue", 
                                                                                     "light grey"))
}


