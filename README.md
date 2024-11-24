# Target Enrichment mini-pipeline

<!-- [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://2ly.link/216Di) -->

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](
https://colab.research.google.com/github/ststevanovic/ni_victoria/blob/main/notebook.ipynb
)

# About
Provided collection of protein/gene identifiers build  
associations and apply scoring to answer following questions:  
1. What are important pathways having group of targets or some of the targets from my group?
2. What biological processes depend on my target group?

# Results
The `processed` (default) contains results from queries:
* gene names, kegg to uniprot id mappings, patwhays data.
The main `results` (default) contains:
* Top 10 bio-processes and pathways and enrichment based on target group.

# Run notebook sample code
Open notebook via colab link and run all code cells.  
To create custom dataset as an input -   
Use UniProt's primary accessions as a list of identifiers.
`e.g. ['P35462', 'Q9H244', 'P21452' ]`

# (Advanced) Run cross-profiling
Supports cross-profiling providing additional gene group 
usually from a specific cell type. Best use to support narrow
expression analysis. To apply the option:   
Define the cross-profile as `cross_gene_list`.  

## Notes
Best used for smaller dataset (400 targets);  
KEGG server overloads dictates retrieval rates.  

## Reference
DOI

