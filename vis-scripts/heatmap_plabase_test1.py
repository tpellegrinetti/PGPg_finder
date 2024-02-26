import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import subprocess

def generate_heatmap(df, output_path, title, normalized=True):
    """Generate a heatmap from a DataFrame."""
    plt.rcParams.update({'font.size': 20})
    fig, ax = plt.subplots(figsize=(max(20, df.shape[1]*2), max(10, df.shape[0]*0.5)))
    sns.heatmap(df, cmap="Spectral_r", ax=ax, annot=True)
    ax.set_title(title)
    fig.savefig(output_path, bbox_inches='tight')
    plt.clf()

def generate_facet_heatmap(data, level, output_path, normalized=True):
    """Generate facet heatmaps based on different levels."""
    lv2_levels = data['Lv2'].unique()
    for lv2 in lv2_levels:
        subset = data[data['Lv2'] == lv2]
        df_pivot = subset.pivot_table(index=level, columns="Sample", values="Count", aggfunc='sum', fill_value=0)

        if normalized:
            total_sum = df_pivot.sum().sum()
            df_pivot = df_pivot / total_sum * 100

        if df_pivot.size > 0:
            heatmap_title = f"{'Normalized' if normalized else 'Non-normalized'} gene counts by sample and {level} with Lv2={lv2}"
            generate_heatmap(df_pivot, os.path.join(output_path, f'{level}_heatmap_facet_{lv2}.svg'), heatmap_title, normalized)
        else:
            print(f"No data for {level} with Lv2={lv2}. Skipping heatmap.")

if __name__ == "__main__":
    import sys

    # Parse command line arguments
    gene_counts_path = sys.argv[1]
    out_dir = sys.argv[2]
    pathways_file = sys.argv[3]
    summary = sys.argv[4]

    # Load data
    gene_counts = pd.read_csv(gene_counts_path, sep="\t")
    pathways = pd.read_csv(pathways_file, sep="\t")
    gene_counts['ID'] = gene_counts['ID'].str.split("_", expand=True)[0]

    # Merge gene counts with pathways
    df = gene_counts.merge(pathways, on="ID", how="outer")

    # Create output directories
    tables_dir = os.path.join(out_dir, 'tables')
    figures_dir = os.path.join(out_dir, 'figures')
    os.makedirs(tables_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)

    # Subdirectories for normalized and non-normalized data
    tables_normalized_dir = os.path.join(tables_dir, 'normalized')
    tables_nonnormalized_dir = os.path.join(tables_dir, 'non-normalized')
    figures_normalized_dir = os.path.join(figures_dir, 'normalized')
    figures_nonnormalized_dir = os.path.join(figures_dir, 'non-normalized')

    os.makedirs(tables_normalized_dir, exist_ok=True)
    os.makedirs(tables_nonnormalized_dir, exist_ok=True)
    os.makedirs(figures_normalized_dir, exist_ok=True)
    os.makedirs(figures_nonnormalized_dir, exist_ok=True)

    # Generate aggregated data
    df_sum = df.pivot_table(index='ID', columns='Sample', values='Count', aggfunc='sum', fill_value=0)
    df_sum.reset_index(inplace=True)

    # Merge with pathways dataframe
    df_sum_with_pathways = pd.merge(df_sum, pathways[['ID', 'PGPT_ID', 'Lv1', 'Lv2', 'Lv3', 'Lv4', 'Lv5']], on='ID', how='left')

    # Reorder columns
    cols_to_move = ['ID', 'Lv1', 'Lv2', 'Lv3', 'Lv4', 'Lv5', 'PGPT_ID']
    df_sum_with_pathways = df_sum_with_pathways[cols_to_move + [col for col in df_sum_with_pathways if col not in cols_to_move]]

    # Save the dataframe with pathways
    df_sum_with_pathways.to_csv(os.path.join(tables_dir, 'gene_counts_sum_with_pathways.txt'), sep='\t', index=False)

    # Combine various levels into a single column
    df_sum_with_pathways['Pathway'] = df_sum_with_pathways[['Lv1', 'Lv2', 'Lv3', 'Lv4', 'Lv5', 'PGPT_ID']].apply(lambda row: ';'.join(row.values.astype(str)), axis=1)
    df_sum_with_pathways.drop(columns=['Lv1', 'Lv2', 'Lv3', 'Lv4', 'Lv5', 'PGPT_ID'], inplace=True)

    # Reorder columns again
    df_sum_with_pathways = df_sum_with_pathways[['ID', 'Pathway'] + [col for col in df_sum_with_pathways if col not in ['ID', 'Pathway']]]

    # Rename columns
    df_sum_with_pathways = df_sum_with_pathways.rename(columns={'ID': '#OTU_ID', df_sum_with_pathways.columns[-1]: 'taxonomy'})

    # Save the combined dataframe
    output_file = os.path.join(tables_dir, 'gene_counts_sum_with_combined_pathways.txt')
    df_sum_with_pathways.to_csv(output_file, sep='\t', index=False)

    # Modify header for biom file
    with open(output_file, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write("# Constructed from biom file\n")
        f.write(content.replace("Pathway", "#OTU ID", 1))

    # Convert to biom format
    output_biom_file = os.path.join(out_dir, 'tables/table.json.biom')
    command = f"biom convert -i {output_file} -o {output_biom_file} --to-json --table-type='OTU table' --process-obs-metadata taxonomy"
    subprocess.run(command, shell=True)

    # Generate heatmaps for different levels
    for level in range(5, 2, -1):
        lv = f"Lv{level}"
        df_pivot = df.pivot_table(index=lv, columns="Sample", values="Count", aggfunc='sum', fill_value=0)

        # Save the pivot table
        df_pivot.to_csv(os.path.join(tables_nonnormalized_dir, f'gene_counts_{lv}.txt'), sep='\t')

        # Normalize the data
        total_sum = df_pivot.sum().sum()
        df_pivot_normalized = df_pivot / total_sum * 100
        df_pivot_normalized.to_csv(os.path.join(tables_normalized_dir, f'normalized_gene_counts_{lv}.txt'), sep='\t')

        # Generate heatmaps
        if lv != "Lv5":
            generate_heatmap(df_pivot, os.path.join(figures_nonnormalized_dir, f'{lv}_heatmap.svg'), f"Non-normalized gene counts by sample and {lv}", normalized=False)
            generate_heatmap(df_pivot_normalized, os.path.join(figures_normalized_dir, f'{lv}_heatmap.svg'), f"Normalized gene counts by sample and {lv}", normalized=True)

        # Generate facet heatmaps
        generate_facet_heatmap(df, lv, figures_normalized_dir, normalized=True)
        generate_facet_heatmap(df, lv, figures_nonnormalized_dir, normalized=False)

    # Generate summary info
    summary_df = pd.read_csv(summary, sep="\t")
    df_summary = gene_counts.merge(summary_df, on="ID", how="outer")

    # Create directories for summary tables and figures
    tables_summary_dir = os.path.join(tables_dir, 'summary')
    figures_summary_dir = os.path.join(figures_dir, 'summary')
    os.makedirs(tables_summary_dir, exist_ok=True)
    os.makedirs(figures_summary_dir, exist_ok=True)

    # Normalization and heatmap generation for summary data
    df_summary = df_summary.dropna(subset=['LV_SUM'])
    sum_lv_sum = df_summary.groupby('LV_SUM')['Count'].sum()
    df_summary['Normalized_Count'] = df_summary.apply(lambda row: (row['Count'] / sum_lv_sum[row['LV_SUM']]) * 100, axis=1)

    # Create normalized pivot table
    df_normalized = df_summary.pivot_table(index='LV_SUM', columns='Sample', values='Normalized_Count', aggfunc='sum', fill_value=0)

    # Save the normalized table
    df_normalized.to_csv(os.path.join(tables_summary_dir, 'normalized_summary_table.txt'), sep='\t')

    # Generate heatmap for the normalized summary table
    heatmap_output_path_normalized = os.path.join(figures_summary_dir, 'normalized_summary_heatmap.svg')
    generate_heatmap(df_normalized, heatmap_output_path_normalized, "Normalized Summary Heatmap")

