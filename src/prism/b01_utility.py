import re
from prism.paths import *
import json
import pandas as pd
import os

def custom_parse(fasta_filepath):
    """
    Parses FASTA file custom-made for this program
    \n Custom Format: chr'X' | Type (Ref, Alt, Flank_1/2) | Name / Identifier
    \n Where X is a valid chromosome number with X = 23 and Y = 24
    :param fasta_filepath:
    :return:
    """
    dataframe = []
    with open(fasta_filepath, 'r') as outfile:
        FASTAFile = [l.strip() for l in outfile.readlines()]
    FASTAstrings = {}
    FASTALabel = ''

    # Separate the fasta formatted DNA strings into labels and contents
    # This is so we can further sort them by the contents of their labels
    for line in FASTAFile:
        if line.startswith('>'):
            FASTALabel = line.lower()
            FASTAstrings[FASTALabel] = ''
        else:
            FASTAstrings[FASTALabel] += line

    current_chr = None
    current_name = None
    df_dict = {}

    for label, string in FASTAstrings.items():
        match = re.match(r'^>chr(\w+)\|(\w+)\|(\w+)$', label)
        if not match:
            raise ValueError(
                f"Invalid Format Detected: {label}. DNA headers must follow the following format: chr'X'| Type [Ref, Alt, Flank_1/2] | [Name / ID]")
        chrom, tag, name = match.groups()

        if (current_chr is not None and current_name is not None
                and (chrom != current_chr or name != current_name)):
            # add chromosome and clinical significance - reset dictionary
            df_dict['Name'] = str(current_name)
            df_dict['Chromosome'] = int(current_chr)
            df_dict['ClinicalSignificance'] = 'N/A'
            dataframe.append(df_dict)
            df_dict = {}

        if 'ref' in tag:
            df_dict['ReferenceAlleleVCF'] = str(string)
        elif 'alt' in tag:
            df_dict['AlternateAlleleVCF'] = str(string)
        elif 'flank_1' in tag or 'flank1' in tag:
            df_dict['Flank_1'] = str(string)
        elif 'flank_2' in tag or 'flank2' in tag:
            df_dict['Flank_2'] = str(string)
        else:
            raise ValueError(f"Improper label: {tag}")
        current_chr = chrom
        current_name = name

    if df_dict:
        df_dict['Name'] = str(current_name)
        df_dict['Chromosome'] = int(current_chr)
        df_dict['ClinicalSignificance'] = 'N/A'
        dataframe.append(df_dict)

    return dataframe


class SiftControl:

    def __init__(self):
        self.folder_path = DATASIFT_CONFIGS
        self.filepath = None
        self.config = None

    def LoadConfig(self, model_name):
        self.filepath = self.folder_path / f"{model_name}.json"
        try:
            file_size = self.filepath.stat().st_size
            print(f"Loading {model_name} Sift: {file_size} bytes")
        except FileNotFoundError:
            print(f"Error: File not found at {self.filepath}")

        with open(self.filepath, 'r') as confile:
            self.config = json.load(confile)

        return self.config

    def check(self):
        if self.config is None:
            raise FileNotFoundError(f"Config not detected: call LoadConfig before utilizing Config-dependent functions")
        else:
            return True

    def LoadSift(self):

        if self.check:
            return self.config['features']
        else:
            raise ValueError("Config diagnostic failed")

    def SiftData(self, X_dataframe):
        features = self.LoadSift

        missing_features = set(features) - set(X_dataframe.columns)
        if missing_features:
            raise ValueError(f"Missing {len(missing_features)}: {list(missing_features)}")

        # extra features are fine, the Sift will ignore it anyway
        return X_dataframe[features]


# file navigation - going to use this across prism so its useful to keep here
def file_access(branch, filename):
    if branch == 'gns':
        datapath = GENE_DATABANK / f"{filename}.fasta"
        txt_parser(datapath)
    elif branch == 'scr':
        folder = SCREEN_RESULTS / filename
        for file in os.listdir(folder):
            if '.csv' in str(file) and '_df' not in str(file):
                datapath = folder / file
                csv_parser(datapath)
    elif branch == 'rpr':
        folder = REGEN_CANDIDATES / filename
        for file in os.listdir(folder):
            if '.txt' in str(file):
                datapath = folder / file
                txt_parser(datapath)


def txt_parser(datapath, verbose=True):
    try:
        with open(datapath, 'r') as infile:
            content = infile.read()

            if verbose:
                print(content)

            return content

    except FileNotFoundError:
        print(f"Error: '{datapath}' was not found.")
    except PermissionError:
        print(f"Error: Permission denied when trying to read the file '{datapath}'.")

def csv_parser(datapath, verbose=True):
    try:
        with open(datapath, 'r') as infile:
            content = pd.read_csv(infile)
            content = content.astype(str)
            transposed = content.T


            if verbose:
                print(content)

            return transposed.to_string(index=True)

    except FileNotFoundError:
        print(f"Error: '{datapath}' was not found.")
    except PermissionError:
        print(f"Error: Permission denied when trying to read the file '{datapath}'.")


