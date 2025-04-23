import pickle as pkl
import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

import plotly.express as px
import plotly.graph_objects as go
import statsmodels.api as sm

from genome_helpers import (
    inverse_rank_normalization
)

def adj_phenotypes_for_gwas(phenotype_df, control_df, phenotypes, sample_ids=None):

    phenotype_df = pd.merge(phenotype_df, control_df, on=["Batch"], suffixes=("", "_NB"))
    
    for phenotype in phenotypes:
        phenotype_df[phenotype + "_adj"] = phenotype_df[phenotype] - phenotype_df[phenotype + "_NB"]
        phenotype_df[phenotype + "_adj_irn"] = inverse_rank_normalization(phenotype_df[phenotype + "_adj"])
        
    phenotype_df = ( phenotype_df
        .drop(["Treatment", "n_replica", "Batch"], axis=1)
        .drop(phenotypes, axis=1)
        .drop([ f"{p}_adj" for p in phenotypes], axis=1)
        .drop([ f"{p}_NB" for p in phenotypes], axis=1)
        .drop(["Treatment_NB", "n_replica_NB"], axis=1)
        .pipe(lambda df: df[df.id.isin(sample_ids)] if sample_ids else df)
        .set_index("id")
        .rename({f"{p}_adj_irn": p for p in phenotypes}, axis=1)
    )
    
    return phenotype_df


def run_gwas(data, snps, phenotypes):
    
    results = []

    for snp in snps:
        
        for phenotype in phenotypes:
            
            X = data[snp]
            valid = ~X.isna()
            X = X[valid]
            X = sm.add_constant(X)
            
            # Y si permutamos los datos? Da lo esperado
            #X = X.sample(len(X))
            
            y = data[phenotype][valid]
            y = inverse_rank_normalization(y)
            
            model = sm.OLS(y, X)
            result = model.fit()
            
            results.append({
                'SNP': snp,
                'phenotype': phenotype,
                'p_value': result.pvalues[snp],
                'beta': result.params[snp],
                'r_squared': result.rsquared
            })
    
    results_df = (
        pd.DataFrame(results)
        .pipe(lambda df: df.assign(contig=df.SNP.apply(lambda x: x[0])))
        .pipe(lambda df: df.assign(position=df.SNP.apply(lambda x: x[1])))
        .sort_values("p_value")
    )

    return results_df


def annotation_to_dict(annotation): 
    
    if annotation is None:
        return dict(ID="", locus="")
    else:
        return dict(item.split('=') for item in d.strip(';').split(';'))


def annotate_results(results_df, annotation_df):

    annotations = list()

    count = 0
    assert len(results_df.SNP.iloc[0]) == 2, "There should be a SNP columns with (contig, position)"
    assert "contig" in annotation_df.columns, "There must be a contig column in the annotation file"
    
    for i, row in results_df.iterrows():
    
        variant_contig = row.SNP[0] # row.contig
        variant_position = row.SNP[1] # row.position
    
        matching_rows = annotation_df[
            # (gff_data["contig"] == variant_contig) & 
            (annotation_df["seqid"] == variant_contig) & 
            (annotation_df["start"] <= variant_position) & 
            (annotation_df["end"] >= variant_position)
        ]
        
        annotations.append(matching_rows.attributes if len(matching_rows) > 0 else None)
    
        # Result
        if matching_rows.empty:        
            # print("Variant belongs to the following rows:")
            # print(matching_rows)
            count += 1
            print("Variant does not belong to any row.")
            
    print(count)
    print(f"{count} out of {results_df.shape[0]} did not match any annotation.")
    annotations_df = pd.Series(annotations).apply(lambda x: None if x is None else x.to_list()[0])
    return annotations_df


def compute_distance(signed_dist1, signed_dist2):

    if signed_dist1 * signed_dist2 < 0:
        return 0
    else:
        return min(abs(signed_dist1), abs(signed_dist2))


def get_closest_gene(df):

    dist_to_start = df.start - df.position
    dist_to_end   = df.end   - df.position
    pp = pd.DataFrame([dist_to_start, dist_to_end]).T
    
    min_index = pp.apply(lambda row: compute_distance(row[0], row[1]), axis=1).argmin()
    return df.iloc[min_index]


def manhattan_static(results_df):

    # Sort chromosomes naturally (1, 2, ..., X)
    results_df['Chromosome'] = pd.Categorical(results_df['contig'], 
                                        categories=sorted(results_df['contig'].unique(), key=lambda x: (not x.isdigit(), x)),
                                        ordered=True)
    results_df = results_df.sort_values(['Chromosome', 'position'])
    
    # Add a column for -log10(P-value)
    results_df['-log10(P-value)'] = -np.log10(results_df['p_value'])
    
    # Generate chromosome-specific positions
    chrom_offsets = {}
    cumulative_position = 0
    
    for chrom in results_df['Chromosome'].cat.categories:
        chrom_data = results_df[results_df['Chromosome'] == chrom]
        chrom_offsets[chrom] = cumulative_position
        results_df.loc[results_df['Chromosome'] == chrom, 'Cumulative Position'] = chrom_data['position'] + cumulative_position
        cumulative_position += chrom_data['position'].max()
    
    # Create the plot
    plt.figure(figsize=(12, 6))
    
    # Alternate colors for chromosomes
    colors = ['#1f77b4', '#ff7f0e']
    for i, chrom in enumerate(results_df['Chromosome'].cat.categories):
        chrom_data = results_df[results_df['Chromosome'] == chrom]
        plt.scatter(chrom_data['Cumulative Position'], chrom_data['-log10(P-value)'], 
                    c=colors[i % len(colors)], s=10, label=f"Chr {chrom}")
    
    # Add a genome-wide significance line (optional)
    significance_threshold = -np.log10(5e-8)
    plt.axhline(y=significance_threshold, color='red', linestyle='--', label='Genome-wide significance')
    
    labels = [x[10:-2] for x in chrom_offsets.keys()]
    
    # Customize the plot
    plt.xlabel('Chromosome')
    plt.ylabel('-log10(P-value)')
    plt.title('Manhattan Plot')
    plt.xticks(ticks=[chrom_offsets[chrom] + (results_df[results_df['Chromosome'] == chrom]['position'].max() / 2) for chrom in chrom_offsets],
               labels=labels, rotation=0)
    # plt.legend(loc='upper right')
    plt.tight_layout()
    
    # Show the plot
    plt.show()


def manhattan_interactive(results_df, annotation="annotation"):

    assert annotation in results_df.columns, f"A column named {annotation} must be present"
    assert isinstance(results_df[annotation].iloc[0], str), f"The {annotation} column must be a str"
        
    # Assuming results_df is already defined and contains columns: contig, position, p_value, annotation
    # Sort chromosomes naturally (1, 2, ..., X)
    results_df['Chromosome'] = pd.Categorical(
        results_df['contig'], 
        categories=sorted(results_df['contig'].unique(), key=lambda x: (not x.isdigit(), x)),
        ordered=True)
    
    results_df = results_df.sort_values(['Chromosome', 'position'])
    
    # Add a column for -log10(P-value)
    results_df['-log10(P-value)'] = -np.log10(results_df['p_value'])
    
    # Generate chromosome-specific positions
    chrom_offsets = {}
    cumulative_position = 0
    ticks = []
    
    for chrom in results_df['Chromosome'].cat.categories:
        chrom_data = results_df.query('Chromosome == @chrom')
        chrom_offsets[chrom] = cumulative_position
        results_df.loc[results_df['Chromosome'] == chrom, 'Cumulative Position'] = (
            chrom_data['position'] + cumulative_position
        )
        ticks.append(cumulative_position + (chrom_data['position'].max() / 2))
        cumulative_position += chrom_data['position'].max()
    
    # Create interactive scatter plot
    fig = px.scatter(
        results_df,
        x='Cumulative Position',
        y='-log10(P-value)',
        color='Chromosome',
        hover_data={
            'Chromosome': True,
            'position': True,
            'p_value': True,
            'annotation': True,  # Include annotation in hover data
            'Cumulative Position': False  # Hide cumulative position from hover
        },
        title='Interactive Manhattan Plot',
        labels={'Cumulative Position': 'Chromosome', '-log10(P-value)': '-log10(P-value)'},
        color_discrete_sequence=px.colors.qualitative.Set3
    )
    
    # Add genome-wide significance threshold
    significance_threshold = -np.log10(5e-2/600)
    fig.add_hline(
        y=significance_threshold,
        line_dash="dash",
        line_color="red",
        annotation_text="Genome-wide significance",
        annotation_position="top left"
    )
    
    
    labels = [x[10:-2] for x in  results_df['Chromosome'].cat.categories]
    
    # Customize axis
    fig.update_layout(
        xaxis=dict(
            title='Contig',
            tickvals=ticks,
            ticktext=labels
        ),
        yaxis=dict(title='-log10(P-value)'),
        legend_title='Contig',
        template='plotly_white'
    )
    
    # Show the plot
    fig.show()


def qqplot(p_values):
    
    """
    Generates a QQ-plot comparing the observed -log10(p-values) to the expected -log10(p-values) 
    under the null hypothesis.

    Parameters:
    p_values (array-like): A sequence of p-values to be plotted.
    
    The function creates a scatter plot of sorted observed vs. expected -log10(p-values) and 
    overlays a reference line (y=x) representing the null hypothesis. This is useful for assessing 
    whether the distribution of p-values deviates from the expected uniform distribution.
    """

    observed = -np.log10(p_values)
    expected = -np.log10(np.linspace(1/len(p_values), 1, len(p_values)))
    
    plt.figure(figsize=(6, 6))
    plt.scatter(np.sort(expected), np.sort(observed), c="blue", s=10, label="SNPs")
    plt.plot([0, max(expected)], [0, max(expected)], color="red", linestyle="--", label="Expected under the null (y=x)")
    plt.xlabel("Expected -log10(p)")
    plt.ylabel("Observed -log10(p)")
    plt.title("QQ-plot")
    plt.legend()
    plt.grid(alpha=0.3)