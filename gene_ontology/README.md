#### Gene Ontology Enrichment Analysis

Gene ontology enrichment analysis performed using the R package, <a href="https://bioconductor.org/packages/release/bioc/html/topGO.html">topGO.</a>  

### Running the code:  
For running the code, two input files are required.   
1. A <i>genelist</i> file.  
The genelist file is a tab-delimited file containing two columns:
Column 1: Locus (contains gene or transcript name of interest).  
Column 2: Rank value of interest (e.g. p-values or log fold changes).  

2. A <i>GO database</i> file.  
The GO database file is a tab-delimited file containing two columns:
Column 1: Locus (contains gene or transcript name of interest).
Column 2: Comma separated GO terms (e.g. GO:0000001, GO:0000002, etc.).  


