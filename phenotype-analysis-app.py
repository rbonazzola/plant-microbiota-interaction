import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
from scipy.ndimage import gaussian_filter1d


st.set_page_config(layout="wide")

fitness_dict = {}

fitness_files = glob.glob("data/fitness/*xlsx")
fitness_files = sorted(fitness_files)

suffixes = sorted([ file[-8:-5] for file in fitness_files ])

for i, (file, suffix) in enumerate(zip(fitness_files, suffixes)):
    fitness_dict[suffix] = pd.read_excel(file)

from datetime import timedelta, time

# from plant_plots import (
#     convert_to_float,
#     get_treatment,
#     get_replica,
#     process_df
# )

# Function to convert time strings to datetime.time format
def convert_to_float(time_str):
    parts = time_str.split()
    hours = int(parts[0])
    minutes = int(parts[2]) if len(parts) > 2 else 0
    return hours + minutes / 60

def process_df(df):
    # print(fitness_df.shape)
    df = df.set_index(["Well", "generation number", "treatment"])    
    n_times = (df.shape[1] // 2)
    columns = [("raw")] * n_times + ["blank"] * n_times   

    times = df.iloc[0, :].apply(convert_to_float)
    df.columns = pd.MultiIndex.from_tuples([ (columns[i], time) for i, time in enumerate(times) ])
    df = df.iloc[1:,:]

    df = df.reset_index()
    
    get_treatment = lambda x: x.split(".")[0] if "." in x else x
    get_replica   = lambda x: x.split(".")[1] if "." in x else 1
    
    df['n_replica'] = df.treatment.apply(get_replica)
    df['Treatment'] = df.treatment.apply(get_treatment)
    df = df.drop("treatment", axis=1).rename({"Treatment": "treatment"}, axis=1)
    df = df.set_index(["Well", "generation number", "treatment", "n_replica"])
    df = df.sort_index()
    return df

def plot_opt_density(generation, max_time=None):
    
    # replica = str(replica)
    
    _fitness_1_df = fitness_1_df.copy()    
    # _fitness_1_df = _fitness_1_df[(_fitness_1_df.treatment == treatment)] #  & (_fitness_1_df.n_replica == replica)]
    
    fig, ax = plt.subplots(1, 1, figsize=(5, 2))
    # display(_fitness_1_df)
    
    print(_fitness_1_df[(_fitness_1_df.treatment == 'MS') & (_fitness_1_df.n_replica == '1')].shape) # .iloc[generation-1, 4:]

    ax.plot(_fitness_1_df[(_fitness_1_df.treatment == 'MS') & (_fitness_1_df.n_replica == '1')].iloc[generation-1, 4:], label='MS-1', color='lavender')
    ax.plot(_fitness_1_df[(_fitness_1_df.treatment == 'MS') & (_fitness_1_df.n_replica == '2')].iloc[generation-1, 4:], label='MS-2', color='blue')
    ax.plot(_fitness_1_df[(_fitness_1_df.treatment == 'MS') & (_fitness_1_df.n_replica == '3')].iloc[generation-1, 4:], label='MS-3', color='darkblue')

    ax.plot(_fitness_1_df[(_fitness_1_df.treatment == 'K') & (_fitness_1_df.n_replica == '1')].iloc[generation-1, 4:], label='K-1', color='lime')
    ax.plot(_fitness_1_df[(_fitness_1_df.treatment == 'K') & (_fitness_1_df.n_replica == '2')].iloc[generation-1, 4:], label='K-2', color='limegreen')
    ax.plot(_fitness_1_df[(_fitness_1_df.treatment == 'K') & (_fitness_1_df.n_replica == '3')].iloc[generation-1, 4:], label='K-3', color='forestgreen')
    
    ax.set_title(f"Generation {generation}")
    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Optical density")

    ax.legend(loc='best', fontsize=6)
    ax.set_ylim(0, 1.7)

    if max_time is not None:
        ax.set_xlim(0, max_time)
    
    return fig

fitness_dict = { k: process_df(v) for k, v in fitness_dict.items() }

plant_pheno_tab, fitness_tab = st.tabs(["Plant phenotypes", "Bacterial fitness"])

def iqr(x):
    return x.quantile(0.75) - x.quantile(0.25)

from wfork_streamlit_profiler import Profiler

with Profiler():

    with plant_pheno_tab:
    
        st.write('## Plant phenotypes')
        
        phenotypes_df_raw = pd.read_csv("data/phenotype/Phenotyping_data.csv")
        phenotypes_df = phenotypes_df_raw.copy()
        phenotypes_df['n_replica'] = phenotypes_df.Treatment.apply(lambda x: x.split(".")[1] if "." in x else 1).astype(int)
        phenotypes_df['Treatment'] = phenotypes_df.Treatment.apply(lambda x: x.split(".")[0] if "." in x else x)
        
        with st.sidebar:
        
            st.title("Plant phenotypes")
            treatment = st.select_slider("Tratamiento", ["K", "MS"])
            n_replicas = phenotypes_df.query("Treatment == @treatment").n_replica.nunique()
            replica = st.slider("Réplica", min_value=1, max_value=n_replicas) if n_replicas > 1 else 1
            smoothing = st.slider("Smoothing sigma", min_value=0, max_value=20)
    
            st.title("Bacterial fitness")
            generation = st.slider(label="Generation", min_value=1, max_value=78)
    
        plot_type = st.radio("Plot ratio?", options=["Absolute value", "Ratio to NB"])        
        
        # scatter = st.radio("Scatter", ["scatter", "line (avg)"])
        phenotypes_df_red_ = phenotypes_df.query("n_replica == @replica and Treatment == @treatment").sort_values("Batch")
        phenotypes_df_red_nb_ = phenotypes_df.query("n_replica == 1 and Treatment == 'NB'").sort_values("Batch")
        phenotypes_df_red = phenotypes_df_red_.groupby(['Treatment', 'n_replica', 'Batch']).agg("median").reset_index()
        phenotypes_df_red_nb = phenotypes_df_red_nb_.groupby(['Treatment', 'n_replica', 'Batch']).agg("median").reset_index()
        
        # phenotypes_df_red_std = phenotypes_df_red_.groupby(['Treatment', 'n_replica', 'Batch']).agg("std").reset_index()
        # phenotypes_df_red_nb_std = phenotypes_df_red_nb_.groupby(['Treatment', 'n_replica', 'Batch']).agg("std").reset_index()
        phenotypes_df_red_std = phenotypes_df_red_.groupby(['Treatment', 'n_replica', 'Batch']).agg(iqr).reset_index()
        phenotypes_df_red_nb_std = phenotypes_df_red_nb_.groupby(['Treatment', 'n_replica', 'Batch']).agg(iqr).reset_index()
    
        phenotypes_df_red_min = phenotypes_df_red_.groupby(['Treatment', 'n_replica', 'Batch']).agg("min").reset_index()
        phenotypes_df_red_max = phenotypes_df_red_.groupby(['Treatment', 'n_replica', 'Batch']).agg("max").reset_index()
        phenotypes_df_red_nb_min = phenotypes_df_red_nb_.groupby(['Treatment', 'n_replica', 'Batch']).agg("min").reset_index()
        phenotypes_df_red_nb_max = phenotypes_df_red_nb_.groupby(['Treatment', 'n_replica', 'Batch']).agg("max").reset_index()
        
    
        phenotypes = [('PR_Length', 8), ('LR_number', 30), ('LR_Density', 9)]
        # phenotypes = [('LR_number', 30)]
    
        fig, ax = plt.subplots(len(phenotypes), 1, figsize=(25, 15)) 
    
        if plot_type == "Absolute value":        
            
            add_nb = st.checkbox("Add NB?")        
            T = phenotypes_df_red.Batch
    
            for i, (phenotype, ylim) in enumerate(phenotypes):
                
                try:
                    ax_i = ax[i]
                except:
                    ax_i = ax
    
                Y           = phenotypes_df_red[phenotype]
                Y_ctrl      = phenotypes_df_red_nb[phenotype]
                stdY        = phenotypes_df_red_std[phenotype]
                stdY_ctrl   = phenotypes_df_red_nb_std[phenotype]
                Ymin        = phenotypes_df_red_min[phenotype]
                Ymax        = phenotypes_df_red_max[phenotype]
                Yctrl_min   = phenotypes_df_red_nb_min[phenotype]
                Yctrl_max   = phenotypes_df_red_nb_max[phenotype]
            
                if smoothing != 0:
                    Y = gaussian_filter1d(Y, smoothing)
                    Y_ctrl = gaussian_filter1d(Y_ctrl, smoothing)
    
                print(stdY)
                ax_i.plot(range(len(Y)), Y, color="blue")               
                ax_i.fill_between(range(len(Y)), Y - stdY, Y + stdY, alpha=0.2, label="Std Dev")
                # ax_i.fill_between(range(len(Y)), Ymin, Ymax, alpha=0.2, label="Std Dev")
                # ax_i.scatter(range(len(Y)), Y, color="blue")               
    
                if add_nb: 
                    ax_i.plot(range(len(Y_ctrl)), Y_ctrl, color="red")
                    # ax_i.fill_between(range(len(Y_ctrl)), Yctrl_min, Yctrl_max, alpha=0.2, label="Std Dev", color="red")
                    ax_i.fill_between(range(len(Y_ctrl)), Y_ctrl-stdY_ctrl, Y_ctrl+stdY_ctrl, alpha=0.2, label="Std Dev", color="red")
    
        
                ax_i.set_title(phenotype)
                ax_i.set_xlabel('Batch (generation)')
                # ax_i.set_ylim(0, ylim)
                
                ax_i.set_ylim(min(Ymin.min(), Yctrl_min.min()), max(Yctrl_max.max(), Ymax.max()))
                
            
        elif plot_type == "Ratio to NB":
    
            T = phenotypes_df_red.Batch
            for i, (phenotype, ylim) in enumerate([('PR_Length', 1.5), ('LR_number', 3), ('LR_Density', 4.5)]):
                
                if smoothing == 0:
                    Y = phenotypes_df_red[phenotype] / phenotypes_df_red_nb[phenotype]
                else:
                    Y = gaussian_filter1d(phenotypes_df_red[phenotype] / phenotypes_df_red_nb[phenotype], smoothing)
        
                ax[i].plot(range(len(Y)), Y, color="blue")
                ax[i].set_title(phenotype)
                ax[i].set_xlabel('Batch (generation)')
                ax[i].set_ylim(0, ylim)
            
        fig.suptitle(f"Treatment: {treatment}, replica: {replica}\n\n\n", fontsize=16)
        st.pyplot(fig)
    
    with fitness_tab:
    
        st.write('## Bacterial fitness')
        M = st.slider("Select N in Fitness_Plate_N.M.xlsx", min_value=1, max_value=5)
        st.write(f"Files: Fitness_Plate_N.{M}.xlsx")
    
        fitness_1_df = pd.concat([fitness_dict[f'{N}.{M}'] for N in range(1, 6)])
        fitness_1_df = fitness_1_df['blank'].reset_index()
        fig = plot_opt_density(generation, max_time=10)
        st.pyplot(fig)
    