import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import subprocess

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

    # Generate the table with the sum of IDs
    df_sum = df.pivot_table(index='ID', columns='Sample', values='Count', aggfunc='sum', fill_value=0)
    df_sum.reset_index(inplace=True)

    # Merge with pathways dataframe
    df_sum_with_pathways = pd.merge(df_sum, pathways[['ID', 'PGPT_ID', 'Lv1', 'Lv2', 'Lv3', 'Lv4', 'Lv5']], on='ID', how='left')

    # Reorder the columns
    df_sum_with_pathways = df_sum_with_pathways[['ID', 'Lv1', 'Lv2', 'Lv3', 'Lv4', 'Lv5', 'PGPT_ID'] + [c for c in df_sum_with_pathways.columns if c not in ['ID', 'Lv1', 'Lv2', 'Lv3', 'Lv4', 'Lv5', 'PGPT_ID']]]

    # Save the dataframe with pathways as a separate file
    df_sum_with_pathways.to_csv(os.path.join(tables_dir, 'gene_counts_sum_with_pathways.txt'), sep='\t', index=False)

    # Combine 'Lv1', 'Lv2', 'Lv3', 'Lv4', 'Lv5' and 'PGPT_ID' into a single column
    df_sum_with_pathways['Pathway'] = df_sum_with_pathways[['Lv1', 'Lv2', 'Lv3', 'Lv4', 'Lv5', 'PGPT_ID']].apply(lambda row: ';'.join(row.values.astype(str)), axis=1)
    df_sum_with_pathways.drop(columns=['Lv1', 'Lv2', 'Lv3', 'Lv4', 'Lv5', 'PGPT_ID'], inplace=True)

    # Rearrange the columns
    columns_order = ['ID', 'Pathway'] + [col for col in df_sum_with_pathways.columns if col not in ['ID', 'Pathway']]

    # Mover a segunda coluna para o final
    second_column = columns_order.pop(1)
    columns_order.append(second_column)

    df_sum_with_pathways = df_sum_with_pathways[columns_order]
    
    # Determine the last column's name
    last_column = df_sum_with_pathways.columns[-1]
    
    # Rename the first column to '#OTU_ID'
    df_sum_with_pathways = df_sum_with_pathways.rename(columns={'ID': '#OTU_ID'})

    # Rename the last column
    df_sum_with_pathways = df_sum_with_pathways.rename(columns={last_column: 'taxonomy'})
  
    # Save the dataframe with combined columns as a new file
    output_file = os.path.join(tables_dir, 'gene_counts_sum_with_combined_pathways.txt')
    df_sum_with_pathways.to_csv(output_file, sep='\t', index=False)

    # Modify header and add lines
    with open(output_file, 'r+') as f:
       content = f.read()
       f.seek(0, 0)
       f.write("# Constructed from biom file\n")
       f.write(content.replace("Pathway", "#OTU ID", 1))
       
    output_biom_file = os.path.join(out_dir, 'tables/table.json.biom')
    command = f"biom convert -i {os.path.join(tables_dir, 'gene_counts_sum_with_combined_pathways.txt')} -o {output_biom_file} --to-json --table-type='OTU table' --process-obs-metadata taxonomy"
    subprocess.run(command, shell=True)

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
