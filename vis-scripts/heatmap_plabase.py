import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

def generate_heatmap(df, output_path, title, normalized=True):
    plt.rcParams.update({'font.size': 20})
    fig, ax = plt.subplots(figsize=(max(20, df.shape[1]*2), max(10, df.shape[0]*0.5)))
    sns.heatmap(df, cmap="viridis", ax=ax)
    ax.set_title(title)
    fig.savefig(output_path, bbox_inches='tight')
    plt.clf()

def generate_facet_heatmap(data, level, output_path, normalized=True):
    lv2_levels = data['Lv2'].unique()
    for lv2 in lv2_levels:
        subset = data[data['Lv2'] == lv2]
        df_pivot = subset.pivot_table(index=level, columns="Sample", values="Count", aggfunc='sum', fill_value=0)

        if normalized:
            total_sum = df_pivot.sum().sum()
            df_pivot = df_pivot / total_sum * 100

        if df_pivot.size > 0:
            generate_heatmap(df_pivot, os.path.join(output_path, f'{level}_heatmap_facet_{lv2}.svg'),
                             f"{'Normalized' if normalized else 'Non-normalized'} gene counts by sample and {level} with Lv2={lv2}", normalized)
        else:
            print(f"No data for {level} with Lv2={lv2}. Skipping heatmap.")

if __name__ == "__main__":
    import sys
    gene_counts_path = sys.argv[1]
    out_dir = sys.argv[2]
    pathways_file = sys.argv[3]

    gene_counts = pd.read_csv(gene_counts_path, sep="\t")
    pathways = pd.read_csv(pathways_file, sep="\t")
    gene_counts['ID'] = gene_counts['ID'].str.split("_", expand=True)[0]

    df = gene_counts.merge(pathways, on="ID", how="outer")

    tables_dir = os.path.join(out_dir, 'tables')
    figures_dir = os.path.join(out_dir, 'figures')
    os.makedirs(tables_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)

    tables_normalized_dir = os.path.join(tables_dir, 'normalized')
    tables_nonnormalized_dir = os.path.join(tables_dir, 'non-normalized')
    figures_normalized_dir = os.path.join(figures_dir, 'normalized')
    figures_nonnormalized_dir = os.path.join(figures_dir, 'non-normalized')

    os.makedirs(tables_normalized_dir, exist_ok=True)
    os.makedirs(tables_nonnormalized_dir, exist_ok=True)
    os.makedirs(figures_normalized_dir, exist_ok=True)
    os.makedirs(figures_nonnormalized_dir, exist_ok=True)

    # Generate the table with the sum of IDs and PGPT_IDs
    df_sum = df.pivot_table(index=['ID', 'PGPT_ID'], columns='Sample', values='Count', aggfunc='sum', fill_value=0)
    df_sum.to_csv(os.path.join(tables_dir, 'gene_counts_sum.txt'), sep='\t')

    for level in range(5, 2, -1):
        lv = f"Lv{level}"

        df_pivot = df.pivot_table(index=lv, columns="Sample", values="Count", aggfunc='sum', fill_value=0)
        df_pivot.to_csv(os.path.join(tables_nonnormalized_dir, f'gene_counts_{lv}.txt'), sep='\t')

        total_sum = df_pivot.sum().sum()
        df_pivot_normalized = df_pivot / total_sum * 100
        df_pivot_normalized.to_csv(os.path.join(tables_normalized_dir, f'normalized_gene_counts_{lv}.txt'), sep='\t')

        if lv != "Lv5":
            generate_heatmap(df_pivot, os.path.join(figures_nonnormalized_dir, f'{lv}_heatmap.svg'),
                             f"Non-normalized gene counts by sample and {lv}", normalized=False)

            generate_heatmap(df_pivot_normalized, os.path.join(figures_normalized_dir, f'{lv}_heatmap.svg'),
                             f"Normalized gene counts by sample and {lv}", normalized=True)

        generate_facet_heatmap(df, lv, figures_normalized_dir, normalized=True)
        generate_facet_heatmap(df, lv, figures_nonnormalized_dir, normalized=False)

