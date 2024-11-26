# Target-based Enrichment Data Pipeline

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://2ly.link/216Di)

# About
This data pipeline processes a collection of protein/gene identifiers (using UniProt accession codes) by querying relevant databases. It generates rankings and visualizations to address the following key questions:
* Which important pathways involve my group of targets?
* What biological processes depend on my target group?

# Results directory
By default, the pipeline generates two subdirectories:
The subdir `processed` contains curated results from database queries, including:
* gene names, KEGG gene IDs, KEGG to UniProt ID mappings and associated patwhays data.
The main `results` subdirectory contains:
* Includes enrichment data and ranked (Top 10 by default) biological processes and pathways.  

# How to use
Run all cells in colab notebook to test the suite OR create custom dataset as an input -   
Use UniProt's primary accessions as a list of identifiers.
`e.g. ['P35462', 'Q9H244', 'P21452' ]`

## Notes
Best used for smaller dataset (400 targets);  
KEGG server overloads dictates retrieval rates.  
