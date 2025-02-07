import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import glob
import ipywidgets as widgets
from ipywidgets import interact
from datetime import timedelta, time

# Function to convert time strings to datetime.time format
def convert_to_float(time_str):
    parts = time_str.split()
    hours = int(parts[0])
    minutes = int(parts[2]) if len(parts) > 2 else 0
    return hours + minutes / 60

get_treatment = lambda x: x.split(".")[0] if "." in x else x
get_replica   = lambda x: int(x.split(".")[1]) if "." in x else 1

def process_df(df):
    # print(fitness_df.shape)
    df = df.set_index(["Well", "generation number", "treatment"])    
    n_times = (df.shape[1] // 2)
    columns = [("raw")] * n_times + ["blank"] * n_times   

    times = df.iloc[0, :].apply(convert_to_float)
    df.columns = pd.MultiIndex.from_tuples([ (columns[i], time) for i, time in enumerate(times) ])
    df = df.iloc[1:,:]

    df = df.reset_index()    
    
    df['n_replica'] = df.treatment.apply(get_replica)
    df['Treatment'] = df.treatment.apply(get_treatment)

    df = df.drop("treatment", axis=1, level=0).rename({"Treatment": "treatment"}, axis=1)
    df = df.set_index(["Well", "generation number", "treatment", "n_replica"])
    df = df.sort_index()
    return df