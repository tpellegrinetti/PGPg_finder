#!/bin/bash

# Function to display help
display_help() {
    echo "Usage: $0 -i <input_directory> -o <output_directory> -t <threads>"
    echo
    echo "This script performs gene search with DIAMOND on a collection of genomes,"
    echo "generates a gene count table, and creates gene presence/absence heatmaps."
    echo "Results are saved to the output directory."
    echo
    echo "Arguments:"
    echo "  -i    Path to the directory containing genome files (.fasta, .fna, .fa)."
    echo "  -o    Path to the directory where results will be saved."
    echo "  -t    Number of CPU threads to be used in the pipeline."
    echo
}

# Check if no argument was provided
if [ $# -eq 0 ]; then
    display_help
    exit 1
fi

# Read command line options
while getopts "i:o:t:h" option; do
    case $option in
        i) genomes_dir=$OPTARG ;;
        o) out_dir=$OPTARG ;;
        t) threads=$OPTARG ;;
        h) display_help
           exit 0 ;;
        *) display_help
           exit 1 ;;
    esac
done

# Path to the DIAMOND database
diamond_db="database/genome.dmnd"
        
# Check if all required options were provided
if [ -z "$genomes_dir" ] || [ -z "$out_dir" ] || [ -z "$threads" ]; then
    echo "Error: All options are required."
    display_help
    exit 1
fi

# Find all genome files with .fasta, .fna, or .fa extensions
genome_files=($(find $genomes_dir -type f \( -name "*.fasta" -o -name "*.fna" -o -name "*.fa" \)))

# Create a file to store gene counts
echo -e "Sample\tID\tCount" > ${out_dir}/gene_counts.txt

# Iterate over each genome file
for genome in "${genome_files[@]}"; do
    sample=$(basename ${genome} .fasta)
    sample=$(basename ${sample} .fna)
    sample=$(basename ${sample} .fa)
    echo "Processing sample ${sample}..."
 
    # Run DIAMOND
    diamond blastx -d ${diamond_db} -q ${genome} -o ${out_dir}/${sample}_diamond.txt -k 1 -e 0.0001 -p ${threads}
 
    # Count the number of genes that match and print the occurrence count
    awk '{print $2}' ${out_dir}/${sample}_diamond.txt | sort | uniq -c | while read count gene; do
        echo -e "${sample}\t${gene}\t${count}" >> ${out_dir}/gene_counts.txt
    done
done


# Python script to generate the heatmaps
python vis-scripts/heatmap_plabase.py ${out_dir}/gene_counts.txt ${out_dir} database/pathways_plabase.txt

