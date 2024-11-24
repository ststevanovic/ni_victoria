import os
import json
import pickle
import pandas as pd
from collections import defaultdict


class DatasetUtils:
    def __init__ (self, schema={
            "KEGGID"    : "kegg_ids.json",
            "GENENAME"  : "gene_names.json",
            "PATHWAYS"  : "pathways.json",
            "PATHWAYS_UNQ" : "pathways_unq.json"
        }):
        self.schema = schema
        self.mapping = defaultdict(dict)

    def load_json(self, path, fn=None):
        if fn:
            path = os.path.join(path, fn)
        return json.load(open(path, "r"))
    
    def save_data(self, data, path, fn=None):
        """
        Save data to a file based on the file extension.
        
        :param data: Data to be saved.
        :param path: Path to the output file (supports .txt, .json, .pkl).
        """
        os.makedirs(path, exist_ok=True)

        if fn:
            path = os.path.join(path, fn)

        if path.endswith(".txt"):
            with open(path, "w") as ofile:
                ofile.write("\n".join(data))
        elif path.endswith(".json"):
            with open(path, "w") as ofile:
                json.dump(data, ofile, indent=4)
        elif path.endswith(".pkl"):
            with open(path, "wb") as ofile:
                pickle.dump(data, ofile)
        else:
            raise ValueError("Unsupported file format. Use .txt, .json, or .pkl")

    def save_mapping(self, raw_object, path):
        """
        Save data to a file based on the file extension.
        
        :param data: mapping to be saved.
        :param path: Path to the output file (base dir name).
        """
        mapping_object = getattr(raw_object, "mapping", None) or raw_object
        for k, v in self.schema.items():
            try:
                self.save_data(mapping_object[k], path, v)
            except:
                raise TypeError(f"Unsupported structure - Can't save mapped data : {k} - {v }")

    def load_mapping(self, base_dir):
        for k, fn in self.schema.items():
            self.mapping[k] = self.load_json(base_dir, fn)
        return self.mapping
            # try:
            # except:
            #     print(k, fn)



class Dataset:
    def __init__(self, path, usecols=[], skiprows=0, skip_sheet=None):
        """
        DatasetUtils to load and save datasets.
        
        :param path: Path to the file (supports .csv, .xlsx, .json, .pkl).
        :param usecols: Column names expected in dataframes (if applicable).
        :param skiprows: Rows to skip (Excel only).
        :param skip_sheet: Specific sheet name for special handling (Excel only).
        """
        self.path = path
        self.usecols = usecols
        self.skiprows = skiprows
        self.skip_sheet = skip_sheet
        self._data = set()

    @property
    def data(self):
        """Wrapper for accessing the data."""
        return self._data
    
    @data.setter
    def data(self, new_data):
        """Intercept updates to self.data."""
        if isinstance(new_data, (set, list)) and isinstance(new_data[0], (str, int)):
            flattened_data = set([li for subli in [l.split(" ") for l in new_data] for li in subli])
            self._data.update(flattened_data)
        else:
            # print("Warning! Unsupported data type for data transformation.")
            self._data.update(new_data)

    def load(self):
        """Load data based on the file type."""
        if self.path.endswith(".csv"):
            self._load_csv()
        elif self.path.endswith(".xlsx"):
            self._load_excel()
        elif self.path.endswith(".json"):
            self._load_json()
        elif self.path.endswith(".pkl"):
            self._load_pickle()
        else:
            raise ValueError("Unsupported file format. Use .csv, .xlsx, .json, or .pkl")
        return self.data

    def _load_csv(self):
        """Load data from a CSV file."""
        if bool(self.usecols):
            data = pd.read_csv(
                self.path, 
                usecols=lambda col: any(col == uc for uc in self.usecols)
                ).iloc[:, 0].drop_duplicates().dropna().tolist()
            self.data = data
        else:
            self.data = pd.read_csv(self.path)

    def _load_excel(self):
        """Load data from an Excel file."""
        sheets = pd.read_excel(self.path, sheet_name=None)
        if bool(self.usecols):
            for sheet in sheets.keys():
                skip = self.skiprows if sheet == self.skip_sheet else 0
                data = pd.read_excel(
                    self.path, 
                    sheet_name=sheet, 
                    skiprows=skip, 
                    usecols=lambda col: any(col == uc for uc in self.usecols)
                    ).iloc[:, 0].drop_duplicates().dropna().tolist()
                self.data = data
        else:
            self.data = sheets

    def _load_json(self):
        """Load data from a JSON file."""
        data = json.load(open(self.path, "r"))
        self.data = data
    
    def _load_pickle(self):
        """Load data from a Pickle file."""
        data = pickle.load(open(self.path, "rb"))
        self.data = data


   
        
   
    
   