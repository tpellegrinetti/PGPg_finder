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

# Function to log messages
log() {
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] $1" | tee -a "$log_file"
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
    log "Error: All options are required."
    display_help
    exit 1
fi

# Create output directory if it doesn't exist
if [ ! -d "$out_dir" ]; then
    mkdir -p "$out_dir"
    log "Created output directory: $out_dir"
fi

# Create log file
log_file="${out_dir}/log.txt"
touch "$log_file"

# Find all genome files with .fasta, .fna, or .fa extensions
genome_files=($(find "$genomes_dir" -type f \( -name "*.fasta" -o -name "*.fna" -o -name "*.fa" \)))

# Create a file to store gene counts
gene_counts_file="${out_dir}/gene_counts.txt"
echo -e "Sample\tID\tCount" > "$gene_counts_file"
log "Created gene counts file: $gene_counts_file"

# Iterate over each genome file
for genome in "${genome_files[@]}"; do
    sample=$(basename "${genome}" .fasta)
    sample=$(basename "${sample}" .fna)
    sample=$(basename "${sample}" .fa)
    log "Processing sample ${sample}..."

    # Run prodigal to generate protein sequences
    prodigal -i "${genome}" -a "${out_dir}/${sample}_proteins.fa" -p single
    log "Generated protein sequences for sample ${sample}"

    # Run DIAMOND
    diamond blastp -d "${diamond_db}" -q "${out_dir}/${sample}_proteins.fa" -o "${out_dir}/${sample}_diamond.txt" -k 1 -e 0.0001 -p "${threads}"
    log "Completed DIAMOND search for sample ${sample}"

    # Count the number of genes that match and print the occurrence count
    awk '{print $2}' "${out_dir}/${sample}_diamond.txt" | sort | uniq -c | while read count gene; do
        echo -e "${sample}\t${gene}\t${count}" >> "${gene_counts_file}"
    done
    log "Generated gene counts for sample ${sample}"
done

# Python script to generate the heatmaps
heatmap_script="vis-scripts/heatmap_plabase.py"
python "$heatmap_script" "${gene_counts_file}" "${out_dir}" database/pathways_plabase.txt
log "Generated heatmaps"

log "Gene search and heatmap generation completed. Check ${out_dir}/gene_counts.txt for the results."

