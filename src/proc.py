from Bio.KEGG.KGML.KGML_parser import read
from collections import defaultdict
import Bio.UniProt as uniprot
import pandas as pd
import requests 


def read_target_fishing_xls(path   = "Targetfishing_res.xlsx"):
    dfs     = pd.read_excel(path, sheet_name=None)
    sheets  = [sheet for sheet in dfs.keys()]
    target_ids = set()
    for sheet in sheets:
        if sheet == "SwissTarget":
            df = pd.read_excel(path, sheet_name=sheet, skiprows=1)
        else:
            df = pd.read_excel(path, sheet_name=sheet)
        identifier_col = [col for col in df.columns if col in ['Uniprot ID', 'Query', 'Uniprot_ID']][0]
        assert len([identifier_col]) == 1
        for t in df[identifier_col].to_list():
            target_ids.add(t)
    return target_ids

def flatten(all_t):
    all_t = [l.split(" ") for l in all_t]
    return set([li for subli in all_t for li in subli])

# targets = read_target_fishing_xls("Targetfishing_res.xlsx")
# targets = flatten(targets)
# print(f"There are {len(targets)} potential napthalimide targets.")

def get_kegg_id(uniprot_id="P05129") -> dict.values:
    """Retrieve kegg ids from uniprot"""
    res = uniprot.search(uniprot_id)
    kegg_ids_global = defaultdict(dict)

    if len(res.results_cache) > 1:
        print(f"Multiple results found for {uniprot_id}")

    for en, res_results_cache in enumerate(res.results_cache):
        try:
            kegg_data = [
                x for x in res_results_cache['uniProtKBCrossReferences']
                if x['database'] == "KEGG"
            ]
        
            if len(kegg_data) > 1:
                print(f"Multiple KEGG cross-references found for {uniprot_id}")
            
            kegg_ids_global[uniprot_id][en] = [k["id"] for k in kegg_data]
        except:
            continue 
    
    return kegg_ids_global[uniprot_id]

def filter_none_hsa(data):
    # Iterate through the dictionary
    cleaned_data = {}
    
    for key, value in data.items():
        # Filter only "hsa" entries in the lists of values
        cleaned_values = {k: [v for v in vs if 'hsa' in v] for k, vs in value.items()}
        
        # Remove empty lists
        cleaned_values = {k: vs for k, vs in cleaned_values.items() if vs}
        
        # If the value is not empty, add to the cleaned data
        if cleaned_values:
            cleaned_data[key] = cleaned_values
    
    return cleaned_data

kegg_identifiers_hsa = filter_none_hsa(kegg_identifiers) # Keep as mapping
kegg_identifiers_hsa_list_all = set([li for subli in [[hsa[0] for hsa in val.values()] for val in kegg_identifiers_hsa.values()] for li in subli])
print(f"There are {len(kegg_identifiers_hsa)} human protein KO identifier groups and a total of {len(kegg_identifiers_hsa_list_all)} hsa identifiers.")


def get_gene_names(uniprot_id):
    
    gene_names_glob = {}
    
    res = uniprot.search(uniprot_id)
    if len(res.results_cache) > 1:
        print(f"Has multiple results cache. {uniprot_id}")

    total_len = len(res.results_cache)
    for en, res_results_cache in enumerate(res.results_cache):
        try:
            gene_data = res_results_cache['genes']
        except:
            assert 'genes' not in res_results_cache.keys(), "There are genes info in a record."
            gene_data = None
        
        gene_names = {}

        if bool(gene_data):
            for g_en, g in enumerate(gene_data) :
                try:
                    gene_names[g_en] = g['geneName']['value']
                except:
                    assert "geneName" not in g.keys(), "There is geneName attribute."
                    gene_names[g_en] = None
        else:
            print(f"No gene name in {en}-th cache - {uniprot_id}. Total check {total_len}")

        gene_names_glob[en] = gene_names

    return gene_names_glob
        
# assoc_to_gene = {key: get_gene_names(key) for key in kegg_identifiers_hsa.keys()}


def kegg_id_to_pathways(kegg_id):
    url = f"https://rest.kegg.jp/link/pathway/{kegg_id}"
    response = requests.get(url)
    if response.ok:
        # Extract pathway IDs from the response
        try:
            pathways = [line.split('\t')[1] for line in response.text.strip().split('\n')]
            return pathways
        except:
            # Some of the paths will return empty response text,e.g '\n'
            # Probably due to lack of annotations...
            print(f"No pathway found: {kegg_id}")
            return []
    else:
        print(f"Error retrieving pathways for KEGG ID {kegg_id}")
        return []
    
# hsa_pathways = {hsa_id : kegg_id_to_pathways(hsa_id) for hsa_id in kegg_identifiers_hsa_list_all}


def clean_hsa_pathways(hsa_pathways):
    # Filter None 
    hsa_pathways_cleaned = { k : v for k, v in hsa_pathways.items() if bool(v)}
    # Spread
    hsa_pathways_cleaned_list = [li for subli in hsa_pathways_cleaned.values() for li in subli]
    # Make set
    hsa_pathways_cleaned_set = set(hsa_pathways_cleaned_list)
    print(
        "Number of pathway groups per hsa:", len(hsa_pathways_cleaned), "\n",
        "Total number of pathways:", len(hsa_pathways_cleaned_list), "\n"
        "Total number of unique pathways:", len(hsa_pathways_cleaned_set))
    
    return hsa_pathways_cleaned_set


def retrieve_kegg_pathway_kgml(pathway_id):
    # Fetch KGML file from KEGG REST API
    url = f"https://rest.kegg.jp/get/{pathway_id}/kgml"
    response = requests.get(url)
    if response.status_code != 200:
        # raise ValueError(f"Failed to retrieve KGML for {pathway_id}")
        print(f"Failed to retrieve KGML for {pathway_id}")
        pathway = None
    else:
        # Parse KGML content using BioPython
        pathway = read(response.text)
    
    return pathway

# path_kgmls = {query: retrieve_kegg_pathway_kgml(query) for query in hsa_pathways_cleaned_set}

