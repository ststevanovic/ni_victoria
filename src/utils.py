import os 
import json 
import pickle



def save_data(targets, kegg_identifiers, hsa_pathways, assoc_to_gene, path_kgmls, out="data_processed"):
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "targets.txt"), "w") as ofile:
        ofile.write("\n".join(targets))
    with open(os.path.join(out, "kegg_identifiers.json"), "w") as ofile:
        json.dump(kegg_identifiers, ofile)

    with open(os.path.join(out, "hsa_pathways.json"), "w") as ofile:
        json.dump(hsa_pathways, ofile)

    with open(os.path.join(out, "assoc_to_gene.json"), "w") as ofile:
        json.dump(assoc_to_gene, ofile)

    with open(os.path.join(out, "hsa_pathways_kgmls.pkl"), "wb") as ofile:
        pickle.dump(path_kgmls, ofile)


# Load input data
def load_data(out="data_processed"):
    targets = json.load(open(os.path.join(out, "targets.txt"), "r"))
    kegg_identifiers = json.load(open(os.path.join(out, "kegg_identifiers.json"), "r"))
    hsa_pathways = json.load(open(os.path.join(out, "hsa_pathways.json"), "r"))
    assoc_to_gene = json.load(open(os.path.join(out, "assoc_to_gene.json"), "r"))
    path_kgml = pickle.load(open(os.path.join(out, "hsa_pathways_kgmls.pkl"), "rb"))

    return targets, kegg_identifiers, hsa_pathways, assoc_to_gene, path_kgml

