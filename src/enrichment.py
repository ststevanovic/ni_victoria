import matplotlib.pyplot as plt
import gseapy as gp
from collections import defaultdict
from itertools import chain
import os

class Utils:
    @staticmethod
    def flatten_list(nested_list):
        return [li for sli in nested_list for li in sli]

    @staticmethod
    def unpack_dl(nested_dict):
        return list(chain.from_iterable(
            v for lst in nested_dict.values()
            for d in lst if isinstance(d, dict)
            for v in d.values()
        ))
    @staticmethod
    def p2a_mapping(nested_dict):
        di_recon=defaultdict(list)
        for p_code, pathways in nested_dict.items():
            for _, paths in pathways.items():
                for path in paths:
                    di_recon[path].append(p_code)
        return di_recon
    
class Enrichment:
    def __init__(self, mapping_object, pathways_name_mapping, organism="human", view_top_n=10, output=None):
        self.pth   = self._ex_map(mapping_object)['PATHWAYS']
        self.gen   = self._ex_map(mapping_object)['GENENAME']
        self.ptu   = self._ex_map(mapping_object)['PATHWAYS_UNQ']
        self.kid   = self._ex_map(mapping_object)['KEGGID']
        self.p2n   = pathways_name_mapping
        self.organism = organism
        self.view_top_n = view_top_n 
        self.output = output
        # Default
        self.k2g_map = {}
        self.enr = None
        self._make_out()

    def _make_out(self):
        if bool(self.output):
            os.makedirs(self.output, exist_ok=True)

    def _save_stats(self, fn, data):
        if bool(self.output):
            output = os.path.join(self.output, fn)
            with open(output, "w") as ifile:
                ifile.write(data)
    
    def _save_fig(self, fn, plt_obj):
        if bool(self.output):
            output = os.path.join(self.output, fn)
            plt_obj.savefig(output)

    def _ex_map(self, mobj):
        return getattr(mobj, "mapping", None) or mobj

    def enricher(self, gene_list):
        Enricher  = gp.enrichr
        gene_sets = "GO_Biological_Process_2021" 
        my_enrichment = Enricher(gene_list=gene_list, gene_sets=gene_sets, organism=self.organism)
        self.enr = my_enrichment.results
        self.enr.to_csv(os.path.join(self.output, "enrichment.csv"), index=False)

    def create_kegg_to_gene_name_mapper(self):
        for k, v in self.kid.items():
            for org in v.values():
                for h in org: # many to 1, keep all
                    self.k2g_map[h] = k 

    def produce_stats_inner_comparison(self):

        targets_len = len(self.pth)

        # Restructure input to map pathways to P-codes
        p2a_map = Utils.p2a_mapping(self.pth)
        
        occupancy_stats = {pth : round(len(targets_list)/targets_len *100, 2) for pth, targets_list in p2a_map.items()}
        
        # Sort pathways per amount of targets present, descending
        data_stats = dict(sorted(occupancy_stats.items(), key=lambda item: item[1], reverse=True))

        # Extract top 10
        x = list(data_stats.keys())[:self.view_top_n]   # str-int
        y = list(data_stats.values())[:self.view_top_n] # score

        # Label targets on each x point
        n_to_show = 4
        x_annot = [
            sorted(p2a_map[i])[:n_to_show] + 
            [f'{max(0, len(sorted(p2a_map[i]))-n_to_show)} more'] for i in x]
        
        # Pathways names on x
        x = [self.p2n[i] for i in x][:self.view_top_n]

        # Setup plot
        fig, ax = plt.subplots(figsize=(20, 16))
        bars = ax.bar(x, y)
        annot_text = '\n'.join(x_annot[en])
        # Annotating each bar with some target names...
        for en, (bar, _) in enumerate(zip(bars, y)):
            # Adjust the placement of text based on the bar's height
            ax.text(
                bar.get_x() + bar.get_width() /2,      # X-coordinate: center of the bar
                bar.get_height() - 5,                  # Y-coordinate
                f"{annot_text}",                       # Text to display         
                ha='center',                           # Align horizontally
                va='bottom',                           # Align vertically
                color="white",                         # Font color
                fontsize=14                            # Font size
            )

        ax.set_xlabel("Pathways", fontsize=18)
        ax.set_ylabel(f"Percent of targets (total {targets_len})", fontsize=18)
        ax.set_title("Top 10: Targets per Pathway (annotated)", fontsize=24)
        ax.set_xticks(range(len(x)))
        ax.set_xticklabels(x, rotation=45, ha='right', fontsize=18)
        fig.tight_layout()
        self._save_fig("Top10_Pathways.png", fig)
        plt.show()

    def get_occupancy_rank(self, row, gene_list):
        coll = row.Genes.split(";")
        count = 0 
        icount = 0
        for gene in coll:
            if gene in gene_list :
                count+=1
            else:
                icount+=1
        return count / len(coll) * 100

    def process_stats_from_BP_enrichment(self, gene_list):
        df = self.enr
        
        df.loc[:, 'Target_Gene_Occupancy'] = df.apply(
            lambda row: self.get_occupancy_rank(row, gene_list), axis=1
            )

        # p values < 0.05 drop insignificant BP...
        df = df[df['Adjusted P-value'] < 0.05].copy()

        # Overlap: The ratio (number of your genes in this term / total genes in the term). 
        # 9/17 means 9 of your genes are found in a process that includes 17 genes in total, indicating substantial overlap.
        df.loc[:, 'Overlap_perc'] = df.Overlap.apply(
            lambda x: round(int(x.split("/")[0]) / int(x.split("/")[1]) * 100, 2)
            )

        # Combining all to compensate log2Fold change
        sort_cols = ['Combined Score', "Overlap_perc", "Target_Gene_Occupancy"]

        df_sorted = df.sort_values(by=sort_cols, ascending=False).reset_index(drop=True)

        # Higher the combined score the better ( higher enrichment )
        fig = plt.figure()
        df_sorted['Combined Score'].iloc[:10].plot(title="Top 10 BioProc (Combined Scores)")
        self._save_fig("Top10_BioProcess.png", fig)
        plt.show()

        fig = plt.figure()
        df_sorted['Combined Score'].plot(title="All BioProc (Combined Scores)")
        self._save_fig("All_BioProcess.png", fig)
        plt.show()

        stxt1 = "Top 10 Biological Process Terms:\n==========\n  ------"
        stxt2 = "\n ".join(df_sorted.Term.to_list()[:self.view_top_n])
        stxt3 = "\n ------"

        stats = f"{stxt1} {stxt2} {stxt3}"
        print(stats)
        self._save_stats("rated_bioprocess.txt", stats)

    def process(self):
        gene_list = Utils.flatten_list(self.gen.values())

        self.enricher(gene_list)
        self.create_kegg_to_gene_name_mapper()
        self.produce_stats_inner_comparison()
        self.process_stats_from_BP_enrichment(gene_list)
