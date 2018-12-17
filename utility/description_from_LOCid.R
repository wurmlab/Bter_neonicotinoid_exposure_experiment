#
# Usage: 
# genes <-  c("LOC100647169", "LOC100650804", "LOC100644325", "LOC100642741")
# description_from_LOCid(genes)
# returns data frame of id, description, both (helpful for figures) 
#
# Author Yannick Wurm
# Copyright Queen Mary University London
# All rights reserved
description_from_LOCid <- function(gene_ids) {

  library(rentrez)

  gene_ncbi_info <- entrez_fetch(db = "gene",
                                 id = gene_ids,
                                 rettype = "xml",
                                 parsed = TRUE)

  tmp_ids <- xpathSApply(doc = gene_ncbi_info,
                         path = "//Entrezgene_gene/Gene-ref/Gene-ref_locus",
                         xmlValue)

  tmp_description <- xpathSApply(doc = gene_ncbi_info,
                                 path = "//Entrezgene_gene/Gene-ref/Gene-ref_desc",
                                 xmlValue)

  tmp_description_clean <- gsub(pattern = " LOC[0-9]+",
                                replacement = "",
                                x = tmp_description)

  desc_for_fig <- paste(tmp_ids, tmp_description_clean, sep = ": ")

  return(data.frame(id = tmp_ids,
                    description = tmp_description,
                    both = desc_for_fig))
}
