#!/bin/bash

# Function to display help
display_help() {
    echo "PGPg_finder v1.1.0 | by Thierry Pellegrinetty <thierry.pellegrinetti@hotmail.com>"
    echo "Check https://github.com/tpellegrinetti/PGPg_finder for updates"
    echo "Usage: $0 -i <input_directory> -o <output_directory> -t <threads> [--piden <min_identity>] [--qcov <min_query_cover>] [--extra <extra_args>] [--bitscore <min_score>] [--evalue <evalue>] [--mode <mode>] [-h]"
    echo
    echo "Read-Based alignment against PLaBAse on a collection of metagenomes,"
    echo "Results are saved to the output directory."
    echo
    echo "Arguments:"
    echo "  -i        Path to the directory containing genome files (.fasta, .fna, .fa)."
    echo "  -o        Path to the directory where results will be saved."
    echo "  -t        Number of CPU threads to be used in the pipeline."
    echo
    echo "Optional Arguments:"
    echo "  --piden     Minimum identity for DIAMOND (default: 30)."
    echo "  --qcov      Minimum query coverage for DIAMOND (default: 30)."
    echo "  --extra     Additional DIAMOND arguments (optional)."
    echo "  --bitscore  Minimum bit score to report alignments."
    echo "  --evalue    Maximum e-value to report alignments (default: 1e-5)."
    echo "  --dmode     DIAMOND mode for sequence search (e.g., fast, sensitive, very-sensitive)."
    echo "  -h          Display this help message."
}

# Function for logging
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
min_query_cover=30
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

# Create log file
log_file="${out_dir}/log.txt"
touch "$log_file"

# Validation of required arguments
if [ -z "$genomes_dir" ] || [ -z "$out_dir" ] || [ -z "$threads" ]; then
    log "Error: Missing required arguments."
    display_help
    exit 1
fi

# Check if DIAMOND database exists
script_dir=$(dirname "$(dirname "$(readlink -f "$0")")")
diamond_db="$script_dir/database/metagenome.dmnd"

if [ ! -f "$diamond_db" ]; then
    log "Error: DIAMOND database not found at $diamond_db."
    exit 1
else
    log "DIAMOND database found at $diamond_db"
fi

# Create output directory if it doesn't exist
if [ ! -d "$out_dir" ]; then
    mkdir -p "$out_dir"
    log "Created output directory: $out_dir"
fi

# Create a file to store gene counts
gene_counts_file="${out_dir}/gene_counts.txt"
echo -e "Sample\tID\tCount" > "$gene_counts_file"
log "Created gene counts file: $gene_counts_file"

# Loop to iterate over each pair of read files
for reads_1 in "$genomes_dir"/*_*1.*; do
    sample=$(basename "$reads_1")
    sample=${sample%_1.*}
    log "Processing sample $sample..."
    
    read_file_1="$reads_1"
    read_file_2="${reads_1/_1./_2.}"
    trimmed_file_1="${out_dir}/${sample}_trimmed_1.fq"
    trimmed_file_2="${out_dir}/${sample}_trimmed_2.fq"
    trimmed_se="${out_dir}/${sample}_trimmed_se.fq"

    # Quality trimming with Trimmomatic
    log "Running Trimmomatic for quality trimming"
    if [ -f "$read_file_2" ]; then
        # Paired-end reads
        trimmomatic PE -threads "$threads" "$read_file_1" "$read_file_2" "$trimmed_file_1" "$trimmed_se" "$trimmed_file_2" "$trimmed_se" SLIDINGWINDOW:4:20 MINLEN:36
    else
        # Single-end reads
        trimmomatic SE -threads "$threads" "$read_file_1" "$trimmed_file_1" SLIDINGWINDOW:4:20 MINLEN:36
    fi

    # Run DIAMOND
    log "Running DIAMOND for PLaBAse alignment for ${sample}"
    diamond blastx -d "$diamond_db" \
        -q "$trimmed_file_1" \
        -o "${out_dir}/${sample}_diamond.txt" \
        -k 1 \
        -p "$threads" \
        -e "$evalue" \
        --id "$min_identity" \
        --query-cover "$min_query_cover" \
        $( [ -n "$min_score" ] && echo "--min-score $min_score" ) \
        $diamond_mode \
        $diamond_extra

    if [ $? -ne 0 ]; then
        log "Error: DIAMOND failed for ${sample}."
        exit 1
    else
        log "Completed DIAMOND search for sample $sample"
    fi

    # Parse DIAMOND output and get gene counts
    awk '{print $2}' "${out_dir}/${sample}_diamond.txt" | sort | uniq -c | while read count id; do
        echo -e "${sample}\t${id}\t${count}" >> "$gene_counts_file"
    done
    log "Generated gene counts for sample $sample"
done

log "Gene search is completed. Check ${gene_counts_file} for the results."

# Python script to generate the heatmaps
heatmap_script="$script_dir/vis-scripts/heatmap_plabase.py"
python "$heatmap_script" "$gene_counts_file" "$out_dir" "$script_dir/database/pathways_plabase.txt" "$script_dir/database/summary.txt"
log "Generated heatmaps"
