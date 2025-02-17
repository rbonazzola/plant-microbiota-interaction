import os
import pysam
import pickle as pkl
import numpy as np
import pandas as pd
from tqdm import tqdm
import re

import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import norm
import statsmodels.api as sm



def read_annotation_data(gff_file, type="CDS", parse_attributes=False):

    col_names = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    
    gff_data = pd.read_csv(
        gff_file,
        sep="\t",
        comment="#",  # Skip comment lines
        names=col_names,
        dtype={"seqid": str, "source": str, "type": str, "start": int, "end": int, "score": str, "strand": str, "phase": str, "attributes": str},
    )
    
    genes = gff_data[gff_data['type'] == type]
    
    if parse_attributes:
        gff_data.attributes = gff_data.attributes.apply(lambda x: x.split(";"))
        gff_data.attributes = gff_data.attributes.apply(lambda xx: {x.split("=")[0]:x.split("=")[1] for x in xx})
    
    return gff_data

def get_genome_metadata(metadata_file="./data/genomes/metadata_whole_genome.xlsx", as_dataframe=False):

    '''
    return: 
    '''

    data = pd.read_excel(metadata_file, engine='openpyxl')

    treatment = data.media.apply(lambda x: re.sub(pattern="[^A-Z]", repl="", string=x))
    replica = data.rep.astype(int)
    generation = data.Generation.apply(lambda x: x if isinstance(x, int) else 0).astype(int)

    batch_mapping = dict(zip(data.BGI_ID, tuple(zip(treatment, replica, generation))))
    """
    { 
      '1_L1': ('K', 1, 4),
      '2_L1': ('K', 3, 74), 
      ...
    }
    """
    
    DUPLICATED_SAMPLE = "19_L6"

    if as_dataframe:
        batch_mapping_df = pd.DataFrame(batch_mapping).T.set_axis(axis=1, labels=["treatment", "replica", "generation"])
        batch_mapping_df = batch_mapping_df.sort_values(["treatment", "replica", "generation"])        
        batch_mapping_df = batch_mapping_df.drop_duplicates()
        return batch_mapping_df
    else:
        batch_mapping.pop(DUPLICATED_SAMPLE)
        return batch_mapping


def compute_af(record):
    
    dp4 = record.info['DP4']
    if sum(dp4) == 0:
        return None
    allele_freq = (dp4[2]+dp4[3]) / sum(dp4)
    return allele_freq


def inverse_rank_normalization(data):
    ranks = data.rank(method="average", na_option="keep")  # Get ranks of the data
    normalized = (ranks - 0.5) / len(data)  # Scale ranks to [0, 1]
    irn = norm.ppf(normalized)  # Map to standard normal
    return irn


def process_plant_phenotype(file="data/phenotype/Phenotyping_data.csv", agg_function="median"):
    
    phenotypes_df_raw = pd.read_csv(file)
    phenotypes_df = phenotypes_df_raw.copy()

    phenotypes_df['n_replica'] = phenotypes_df.Treatment.apply(lambda x: x.split(".")[1] if "." in x else 1).astype(int)
    phenotypes_df['Treatment'] = phenotypes_df.Treatment.apply(lambda x: x.split(".")[0] if "." in x else x)

    phenotypes_df_reds = {
        (treatment, replica): phenotypes_df.query("n_replica == @replica and Treatment == @treatment").sort_values("Batch")
        for treatment in ('K', 'MS')
        for replica in range(1, 4)
    }

    phenotypes_df_nb = phenotypes_df.query("n_replica == 1 and Treatment == 'NB'").sort_values("Batch") # if add_nb else None
    
    phenotypes_df_red_nb = phenotypes_df_nb.groupby(['Treatment', 'n_replica', 'Batch']).agg(agg_function).reset_index()
    phenotypes_df["id"] = phenotypes_df.apply(lambda row: (row.Treatment, row.n_replica, row.Batch), axis=1)

    return phenotypes_df, phenotypes_df_reds, phenotypes_df_nb, phenotypes_df_red_nb