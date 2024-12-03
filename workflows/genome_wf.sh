#!/bin/bash

# Function to display help
display_help() {
    echo "PGPg_finder v1.1.0 | by Thierry Pellegrinetty <thierry.pellegrinetti@hotmail.com>"
    echo "Check https://github.com/tpellegrinetti/PGPg_finder for updates"
    echo "Usage: $0 -i <input_directory> -o <output_directory> -t <threads> [--piden <min_identity>] [--qcov <min_query_cover>] [--extra <extra_args>] [--bitscore <min_score>] [--evalue <evalue>] [--mode <mode>] [-h]"
    echo
    echo "This script performs gene search with DIAMOND on a collection of genomes,"
    echo "generates a gene count table, and creates gene presence/absence heatmaps."
    echo "Results are saved to the output directory."
    echo
    echo "Arguments:"
    echo "  -i        Path to the directory containing genome files (.fasta, .fna, .fa)."
    echo "  -o        Path to the directory where results will be saved."
    echo "  -t        Number of CPU threads to be used in the pipeline."
    echo
    echo "Optional Arguments:"
    echo "  --piden     Minimum identity for DIAMOND (default: 30)."
    echo "  --qcov      Minimum query coverage for DIAMOND (default: 50)."
    echo "  --extra     Additional DIAMOND arguments (optional)."
    echo "  --bitscore  Minimum bit score to report alignments."
    echo "  --evalue    Maximum e-value to report alignments (default: 1e-5)."
    echo "  --dmode      DIAMOND mode for sequence search (e.g., fast, sensitive, very-sensitive)."
    echo "  -h          Display this help message."
}

# Log function
log() {
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] $1" | tee -a "$log_file"
}

# Parse long and short options using `getopt`
ARGS=$(getopt -o i:o:t:h --long piden:,qcov:,extra:,bitscore:,evalue:,dmode:,help -n "$0" -- "$@")
if [ $? -ne 0 ]; then
    display_help
    exit 1
fi
eval set -- "$ARGS"

# Default values
min_identity=30
min_query_cover=50
diamond_extra=""
min_score=""
evalue="1e-5"
diamond_mode=""

# Parse options
while true; do
    case "$1" in
        -i) genomes_dir=$2; shift 2 ;;
        -o) out_dir=$2; shift 2 ;;
        -t) threads=$2; shift 2 ;;
        --piden) min_identity=$2; shift 2 ;;
        --qcov) min_query_cover=$2; shift 2 ;;
        --extra) diamond_extra=$2; shift 2 ;;
        --bitscore) min_score=$2; shift 2 ;;
        --evalue) evalue=$2; shift 2 ;;
        --dmode) diamond_mode=$2; shift 2 ;;
        -h|--help) display_help; exit 0 ;;
        --) shift; break ;;
        *) display_help; exit 1 ;;
    esac
done

# Validation of required arguments
if [ -z "$genomes_dir" ] || [ -z "$out_dir" ] || [ -z "$threads" ]; then
    log "Error: Missing required arguments."
    display_help
    exit 1
fi

# Check if DIAMOND database exists
script_dir=$(dirname "$(dirname "$(readlink -f "$0")")")
diamond_db="$script_dir/database/genome.dmnd"

if [ ! -f "$diamond_db" ]; then
    log "Error: DIAMOND database not found at $diamond_db."
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$out_dir"
log_file="${out_dir}/log.txt"
touch "$log_file"

# Find genome files
genome_files=($(find "$genomes_dir" -type f \( -name "*.fasta" -o -name "*.fna" -o -name "*.fa" \)))

if [ ${#genome_files[@]} -eq 0 ]; then
    log "Error: No genome files found in $genomes_dir."
    exit 1
fi

# Process each genome file
gene_counts_file="${out_dir}/gene_counts.txt"
echo -e "Sample\tID\tCount" > "$gene_counts_file"

for genome in "${genome_files[@]}"; do
    sample=$(basename "${genome%.*}")
    log "Processing sample ${sample}..."

    # Run prodigal
    prodigal -i "${genome}" -a "${out_dir}/${sample}_proteins.fa" -p single
    if [ $? -ne 0 ]; then
        log "Error: Prodigal failed for ${sample}."
        exit 1
    else
        log "Generated protein sequences for sample ${sample}"
    fi
    
    # Run DIAMOND
    diamond blastp -d "$diamond_db" \
        -q "${out_dir}/${sample}_proteins.fa" \
        -o "${out_dir}/${sample}_diamond.txt" \
        -p "$threads" \
        -k 1 \
        -e "$evalue" \
        --id "$min_identity" \
        --query-cover "$min_query_cover" \
        $( [ -n "$min_score" ] && echo "--min-score $min_score" ) \
        $( [ -n "$diamond_mode" ] && echo "--mode $diamond_mode" ) \
        $diamond_extra

    if [ $? -ne 0 ]; then
        log "Error: DIAMOND failed for ${sample}."
        exit 1    
    else
        log "Completed DIAMOND search for sample ${sample}"
    fi

    # Process results
    awk '{print $2}' "${out_dir}/${sample}_diamond.txt" | sort | uniq -c | while read count gene; do
        echo -e "${sample}\t${gene}\t${count}" >> "$gene_counts_file"
    done
    log "Completed processing for ${sample}."
done

# Python script to generate the heatmaps
heatmap_script="$script_dir/vis-scripts/heatmap_plabase.py"
python "$heatmap_script" "${gene_counts_file}" "${out_dir}" "$script_dir/database/pathways_plabase.txt" "$script_dir/database/summary.txt"
log "Generated heatmaps"
log "Pipeline completed. Results are saved in ${out_dir}."
