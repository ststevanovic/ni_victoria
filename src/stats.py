import Bio.KEGG.KGML.KGML_pathway as Pathway
from Bio.KEGG.KGML.KGML_parser import read
import matplotlib.patches as mpatches
from collections import defaultdict
import matplotlib.pyplot as plt
import Bio.UniProt as uniprot
from Bio.KEGG import REST
from io import StringIO
import networkx as nx
import pandas as pd 
import requests 
import re
import gseapy as gp
from collections import Counter


def run_enricher(assoc_to_gene):
    Enricher = gp.enrichr
    # Extract all gene names
    gene_list = list(filter(None, list(set([l[0] for l in list(filter(None, [li for subli in list(filter(None, [[list(k.values()) for k in v.values()] for v in assoc_to_gene.values()])) for li in subli]))]))))
    gene_sets = "GO_Biological_Process_2021" 
    my_enrichment = Enricher(gene_list=gene_list, gene_sets=gene_sets, organism="human")
    return my_enrichment

def create_mapping_hsa_to_access(kegg_identifiers):
    hsa_to_uniprot = {}

    for k, v in kegg_identifiers.items():
        try:
            for hsa in v.values():
                # many to 1, keep all
                for h in hsa:
                    hsa_to_uniprot[h] = k 
        except:
            continue
    return hsa_to_uniprot

def produce_stats_inner_comparison(hsa_pathways_cleaned, hsa_to_uniprot, path_kgmls):
    all_pathways_counts = Counter([l for sl in [h for h in hsa_pathways_cleaned.values()] for l in sl])
    all_pathways_counts_stats = {}
    all_pathways_counts_targets = {}
    len_targets = len(hsa_pathways_cleaned)

    # Collect all targets per pathways
    for pathway, _ in all_pathways_counts.items():
        pathway_ts = []
        for k,v in hsa_pathways_cleaned.items():
            if pathway in v:
                pathway_ts.append(k)
        
        all_pathways_counts_stats[pathway] = round(len(pathway_ts)/len_targets *100, 2)
        all_pathways_counts_targets[pathway] = pathway_ts

    # Sort pathways per amount of targets present, descending
    data_stats = dict(sorted(all_pathways_counts_stats.items(), key=lambda item: item[1], reverse=True))
    data_stats

    # Extract top 10
    x = list(data_stats.keys())[:10]
    y = list(data_stats.values())[:10]

    # Label targets in pathway
    n_to_show = 4
    x_annot = [sorted(all_pathways_counts_targets[i])[:n_to_show] + [f'{len(sorted(all_pathways_counts_targets[i]))-n_to_show} more'] for i in x]
    x_annot_uniprot = [[hsa_to_uniprot[m] for m in i[:-1]] + [i[-1]] for i in x_annot]
    x = [path_kgmls[i].title for i in x]

    # Setup plot
    plt.figure(figsize=(20, 16))
    bars = plt.bar(x, y)

    # Annotating each bar with some target names...
    for en, (bar, value) in enumerate(zip(bars, y)):
        # Adjust the placement of text based on the bar's height
        plt.text(
            bar.get_x() + bar.get_width() /2,      # X-coordinate: center of the bar
            bar.get_height() - 5,                  # Y-coordinate
            f"{"\n".join(x_annot_uniprot[en])}",   # Text to display         
            ha='center',                           # Align horizontally
            va='bottom',                           # Align vertically
            color="white",                         # Font color
            fontsize=14                            # Font size
        )

    plt.xlabel("Pathways", fontsize=18)
    plt.ylabel(f"Percent of targets (total {len_targets})", fontsize=18)
    plt.title("Top 10: Targets per Pathway (annotated)", fontsize=24)
    plt.xticks(rotation=45, ha='right', fontsize=18)
    plt.tight_layout()
    plt.show()

def process_stats_from_BO_enrichment(my_enrichment):
    ddt = my_enrichment.results

    def get_occupancy_rank(row, gene_list):
        coll = row.Genes.split(";")
        count = 0 
        icount = 0
        for gene in coll:
            if gene in gene_list :
                count+=1
            else:
                icount+=1
        return count / len(coll) * 100

    ddt.loc[:, 'Target_Gene_Occupancy'] = ddt.apply(lambda row: get_occupancy_rank(row, gene_list), axis=1)

    # p values < 0.05 drop insignificant BP...
    ddt = ddt[ddt['Adjusted P-value'] < 0.05].copy()

    # Overlap: The ratio (number of your genes in this term / total genes in the term). 
    # 9/17 means 9 of your genes are found in a process that includes 17 genes in total, indicating substantial overlap.
    ddt.loc[:, 'Overlap_perc'] = ddt.Overlap.apply(lambda x: round(int(x.split("/")[0]) / int(x.split("/")[1]) * 100, 2))

    # Combining all to compensate log2Fold change
    sort_cols = ['Combined Score', "Overlap_perc", "Target_Gene_Occupancy"]

    sorted_ddt = ddt.sort_values(by=sort_cols, ascending=False).reset_index(drop=True)

    # Higher the combined score the better ( higher enrichment )
    sorted_ddt['Combined Score'].iloc[:10].plot(title="Top 10 Combined Scores")
    plt.show()
    sorted_ddt['Combined Score'].plot(title="All BP Combined Scores")

    print(f"Top 10 BO process Terms:\n==========\n {"\n".join(sorted_ddt.Term.to_list()[:10])}\n ------")

    
# Run target cross path
# e.g.
# path_id    #targets
# path:      20/50

# Run swa620 cross target path
# swa620 path id  #targets

# Find common swa620 targets BO from enricher
# swa620 BP  target BP 
