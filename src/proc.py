from Bio.KEGG.KGML.KGML_parser import read
from requests.exceptions import ReadTimeout, RequestException
from http.client import IncompleteRead
from collections import defaultdict
from itertools import chain
import Bio.UniProt as uniprot
from tqdm import tqdm
import pandas as pd
import requests 
import pickle
import json
import os
import time
from threading import Lock


# checkpoints TODO

class Utils:
    # Lock to synchronize requests and enforce rate limit
    _last_request_time = 0
    _rate_limit_lock = Lock()

    @staticmethod
    def fetch_with_retries(request_obj, retries=3, delay=2,  rate_limit=None):
        """
        Utility function to handle retries for HTTP requests.
        
        Parameters:
        - request_func: Callable function to execute the request.
        - retries: Number of retry attempts (default: 3).
        - delay: Delay between retries in seconds (default: 2).
        - args, kwargs: Arguments to pass to the request function.
        
        Returns:
        - The result of the request function if successful.
        
        Raises:
        - The last encountered exception if all retries fail.
        """
        request_func, params = request_obj
        for attempt in range(retries):
            try:
                # Enforce rate limit
                if rate_limit:
                    with Utils._rate_limit_lock:
                        current_time = time.time()
                        elapsed = current_time - Utils._last_request_time
                        if elapsed < rate_limit:
                            time.sleep(rate_limit - elapsed)
                        Utils._last_request_time = time.time()


                # Execute the request
                return request_func(params)

            except Exception as e:
                print(f"Error on attempt {attempt + 1}/{retries}: {e}")
                if attempt < retries - 1:
                    time.sleep(delay)
                else:
                    raise # Re-raise the exception if out of retries

    @staticmethod
    def check_list_len(input_list, input_id):
        if len(input_list) > 1:
            print(f"Multiple results found for {input_id}")

    @staticmethod 
    def check_dict_val(input_dict):
        if any(not bool(v) for v in input_dict.values()):
            raise Exception("dict has False val")

    @staticmethod
    def unpack_dl(nested_dict):
        return list(chain.from_iterable(
            v for lst in nested_dict.values()
            for d in lst if isinstance(d, dict)
            for v in d.values()
        ))
    
    @staticmethod
    def unpack_dl(nested_dicts_list):
        return list(chain.from_iterable(
            v for d in nested_dicts_list
            for v in d.values()
        ))

    @staticmethod
    def flatten(input_list):
        return [li for subli in list(input_list) for li in subli]
    
    @staticmethod 
    def drop_duplicates(input_list):
        return set(input_list)
   
class KeggOperations:
    def __init__(self, output:str, uniprot_id: str =None, uniprot_ids: list=None, organism: str="hsa", download: bool=True, verbose: bool = False):
        self.output      = output 
        self.uniprot_id  = uniprot_id
        self.uniprot_ids = uniprot_ids
        self.organism    = organism 
        self.download    = download
        self.mapping     = defaultdict(dict)
        self.records     = None
        self.verbose     = verbose
        self._prep_out()

    def _prep_out(self):
        os.makedirs(self.output, exist_ok=True)

    def get_taxid(self, taxid_mapping={"hsa": 9606}):
        """Translate organism name to Tax ID."""
        return taxid_mapping[self.organism]

    def fetch_uniprot_data(self, uniprot_id=None):
        """Use BioPython UniProt module and retrieve data for a given Uniprot ID"""
        uniprot_id = uniprot_id or self.uniprot_id
        if uniprot_id:
            unirepr = (uniprot.search, f"(organism_id:{self.get_taxid()}) and ({uniprot_id})")
            res  = Utils.fetch_with_retries(unirepr)
            if self.verbose:
                Utils.check_list_len(res.results_cache, uniprot_id)
            if bool(res): # none-matching organism
                self.records = res
        else:
            raise ValueError("Uniprot ID not provided.")

    def get_kegg_id(self, uniprot_id="P05129"):
        """Retrieves KEGG ID providing UniProt Accession ID"""
        
        records = self.records.results_cache 
        
        collector = defaultdict(dict)

        for en, rc in enumerate(records):
            try:
                data = [
                    x for x in rc['uniProtKBCrossReferences']
                    if x['database'] == "KEGG"
                ]

                if len(data) > 1 and self.verbose:
                    print(f"Multiple KEGG cross-references found for {uniprot_id}")
                    # Avoid ambiguity TODO?

                result = [k["id"] for k in data if self.organism in k['id']]
                if bool(result):
                    collector[uniprot_id][en] = result
            
            except Exception as e:
                if self.verbose:
                    print(f"Error processing record for {uniprot_id}: {e}")
                continue
        
        self.mapping['KEGGID'][uniprot_id] = collector[uniprot_id]

    def get_gene_names(self, uniprot_id, cutoff=3): # less greedy? set cutoff to 5 or more.    
        """Retrives gene names provinding uniprot accession. Note - first one is primary"""
        records = self.records.results_cache
        
        collector = list(
            filter(
                lambda x: x, 
                [
                    g.get('geneName', {}).get('value') for r in records[:cutoff]
                    for g in r.get('genes', []) 
                    if g.get('geneName', {}).get('value') is not None
                ]
            )
        )
        if bool(collector):
            self.mapping['GENENAME'][uniprot_id] = collector

    def kegg_id_to_pathways(self, uniprot_id, kegg_id):
        """
        Retrieves pathway IDs. 
        Some paths will return empty response text, e.g '\n'
        Most likely lacking annotation.
        """
        request_obj = (
            lambda params : requests.get(url=params['url'], timeout=params['timeout']), 
            {
                "url": f"https://rest.kegg.jp/link/pathway/{kegg_id}", 
                "timeout": 10
            }
        )

        response = Utils.fetch_with_retries(request_obj, rate_limit=0.333)

        if response.ok:
            try:
                pathways = [
                    line.split('\t')[1] for line in response.text.strip().split('\n')
                    ] 
                if bool(pathways):
                    self.mapping['PATHWAYS'][uniprot_id] = {kegg_id : pathways}
            except:
                if self.verbose:
                    print(f"No pathway found: {kegg_id}")
        else:
            print(f"Error retrieving pathways for KEGG ID {kegg_id}")
        
    def get_pathways(self, uniprot_id, greedy=True):
        kegg_ids = self.mapping['KEGGID'][uniprot_id]            
        kegg_ids = Utils.flatten(kegg_ids.values())

        for i, kegg_id in enumerate(kegg_ids, 1):
            if self.verbose:
                print("Querying KEGG_ID -", {kegg_id})
            self.kegg_id_to_pathways(uniprot_id, kegg_id)
            if greedy and self.mapping['PATHWAYS'].get(uniprot_id):  # Stop early if greedy
                if self.verbose:
                    print(f"Found pathways - {i} retries")
                break

    def clean_pathways(self):
        """Deduplicates and updates unique pathways in the mapping."""
        pathways = list(self.mapping['PATHWAYS'].values())
        pathways_list = Utils.unpack_dl(pathways)
        pathways_set  = Utils.drop_duplicates(pathways_list)
        self.mapping['PATHWAYS_UNQ'] = list(pathways_set)

        print(
            "-- STATS --","\n",
            "Total pathways avilable:", len(pathways_list), "\n"
            "Total unique pathways:", len(pathways_set)
            )
        
    def retrieve_kegg_pathway_kgml(self, pathway_id):
        """Fetch KGML file from KEGG REST API"""

        request_obj = (
            lambda params : requests.get(url=params['url'], timeout=params['timeout']), 
            {
                "url": f"https://rest.kegg.jp/get/{pathway_id}/kgml", 
                "timeout": 10
            }
        )

        response = Utils.fetch_with_retries(request_obj, rate_limit=0.333)

        if response.status_code != 200:
            if self.verbose:
                print(f"Failed to retrieve KGML for {pathway_id}")
            pathway = None
        else:
            pathway = read(response.text)
        
        return pathway
    
    def retrieve_kegg_pathway_name(self, pathway_id):
        """Fetch pathway name from KEGG REST API"""
        
        def make_request(params):
            response = requests.get(params["url"], timeout=params["timeout"])
            if response.status_code == 200:
                second_line = response.text.splitlines()[1]
                if second_line.startswith("NAME"):
                    return second_line.split("  ", 1)[1].strip().split(" - ")[0]
                    # return second_line.split("  ", 1)[1].strip()
                else:
                    raise ValueError("Unexpected response format")
            else:
                raise ValueError(f"HTTP Error {response.status_code}")

        request_obj = (
            make_request, 
                {
                "url": f"https://rest.kegg.jp/get/{pathway_id}", 
                "timeout": 10
                }
            )
        
        try:
            return Utils.fetch_with_retries(request_obj, rate_limit=0.333)
        except Exception as e:
            if self.verbose:
                print(f"Failed to retrieve pathway name for {pathway_id}: {e}")
            return None

    def kegg_batch(self, func, fn, ext="pkl"):
        """Batch collect and save KGML data"""
        tools = {
            "pkl": pickle,
            "json": json,
        }
        typ = {
            "pkl" : "wb",
            "json": "w"
        }
        results = {pathway: func(pathway) for pathway in tqdm(self.mapping['PATHWAYS_UNQ'])}
        with open(os.path.join(self.output, f"{fn}.{ext}"), typ[ext]) as ofile:
            tools[ext].dump(results, ofile)

    def check(self, locvar):
        """Input ID processing from locvar"""
        locvar.pop('self')
        ky, ve = next(iter(locvar.items()))
        key = ve or getattr(self, ky)
        if not key:
            raise ValueError(f"Not provided - {ky}")
        return key

    def no_records(self):
        """Checks if there are records"""
        return not bool(getattr(self.records, "results_cache", None))

    def query(self, uniprot_id=None):
        """Runs the query on Uniprot ID"""
        upid = self.check(locals())
        self.fetch_uniprot_data(upid)
        if self.no_records():
            if self.verbose:
                print(f"Filtered - {uniprot_id}")
            return

        self.get_kegg_id(upid)
        self.get_gene_names(upid)
        self.get_pathways(upid)
        
        self.records = None  # Reset records
        if self.verbose:
            print(f"Query {uniprot_id} - completed.")

    def query_batch(self, uniprot_ids=None):
        """Batch runs the query on a list with Uniprot IDs"""
        for uid in tqdm(self.check(locals())):
            self.query(uid)
        self.clean_pathways()
        self.kegg_batch(self.retrieve_kegg_pathway_name, "pathways_names", "json")