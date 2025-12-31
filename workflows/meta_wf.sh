#!/bin/bash

###############################################################################
# PGPg_finder meta_wf
# Assembly-based metagenomic workflow
###############################################################################

# Function to display help
display_help() {
    echo "PGPg_finder v1.1.0 | by Thierry Pellegrinetty <thierry.pellegrinetti@hotmail.com>"
    echo "Check https://github.com/tpellegrinetti/PGPg_finder for updates"
    echo
    echo "Usage:"
    echo "  $0 -i <reads_directory> -o <output_directory> -t <threads> [-a <assembly_directory>]"
    echo "     [--piden <min_identity>] [--qcov <min_query_cover>] [--extra <extra_args>]"
    echo "     [--bitscore <min_score>] [--evalue <evalue>] [--dmode <mode>] [-h]"
    echo
    echo "Assembly-based metagenomic workflow for detection and quantification of"
    echo "plant growth–promoting genes using PLaBAse."
    echo
    echo "This workflow performs read trimming, optional metagenome assembly,"
    echo "gene prediction, functional annotation, and read mapping to estimate"
    echo "gene abundance across samples."
    echo
    echo "Required arguments:"
    echo "  -i        Directory containing metagenomic read files."
    echo "            Paired-end reads must follow the *_1 / *_2 naming convention."
    echo "  -o        Output directory."
    echo "  -t        Number of CPU threads."
    echo
    echo "Optional arguments:"
    echo "  -a        Directory containing pre-assembled contigs."
    echo "            If provided, the assembly step is skipped."
    echo
    echo "DIAMOND parameters:"
    echo "  --piden     Minimum identity percentage (default: 30)."
    echo "  --qcov      Minimum query coverage percentage (default: 30)."
    echo "  --bitscore  Minimum bit score to report alignments."
    echo "  --evalue    Maximum e-value (default: 1e-5)."
    echo "  --dmode     DIAMOND search mode (fast, sensitive, very-sensitive)."
    echo "  --extra     Additional DIAMOND arguments."
    echo
    echo "Other options:"
    echo "  -h, --help  Display this help message."
}

# Logging function
log() {
    local timestamp
    timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] $1" | tee -a "$log_file"
}

###############################################################################
# Argument parsing
###############################################################################

ARGS=$(getopt -o i:o:t:a:h --long piden:,qcov:,extra:,bitscore:,evalue:,dmode:,help -n "$0" -- "$@")
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

# Parse arguments
while true; do
    case "$1" in
        -i) reads_dir=$2; shift 2 ;;
        -o) out_dir=$2; shift 2 ;;
        -t) threads=$2; shift 2 ;;
        -a) assembly_dir=$2; shift 2 ;;
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

###############################################################################
# Checks and setup
###############################################################################

if [ -z "$reads_dir" ] || [ -z "$out_dir" ] || [ -z "$threads" ]; then
    echo "ERROR: Missing required arguments."
    display_help
    exit 1
fi

mkdir -p "$out_dir"
log_file="${out_dir}/log.txt"
touch "$log_file"

script_dir=$(dirname "$(dirname "$(readlink -f "$0")")")
diamond_db="$script_dir/database/genome.dmnd"

if [ ! -f "$diamond_db" ]; then
    log "ERROR: DIAMOND database not found at $diamond_db"
    exit 1
fi

###############################################################################
# Main loop
###############################################################################

for reads_1 in "$reads_dir"/*_*1.*; do
    sample=$(basename "$reads_1")
    sample=${sample%_1.*}
    log "Processing sample ${sample}"

    read_file_1="$reads_1"
    read_file_2="${reads_1/_1./_2.}"

    trimmed_1="${out_dir}/${sample}_trimmed_1.fq"
    trimmed_2="${out_dir}/${sample}_trimmed_2.fq"
    trimmed_se="${out_dir}/${sample}_trimmed_se.fq"

    log "Running Trimmomatic"
    if [ -f "$read_file_2" ]; then
        trimmomatic PE -threads "$threads" \
            "$read_file_1" "$read_file_2" \
            "$trimmed_1" "$trimmed_se" \
            "$trimmed_2" "$trimmed_se" \
            SLIDINGWINDOW:4:20 MINLEN:36
    else
        trimmomatic SE -threads "$threads" \
            "$read_file_1" "$trimmed_1" \
            SLIDINGWINDOW:4:20 MINLEN:36
    fi

    if [ -z "$assembly_dir" ]; then
        log "Assembling metagenome with MEGAHIT"
        megahit -1 "$trimmed_1" -2 "$trimmed_2" -t "$threads" -o "${out_dir}/${sample}_assembly"
        assembly="${out_dir}/${sample}_assembly/final.contigs.fa"
    else
        log "Using provided assembly"
        assembly=$(ls "${assembly_dir}/${sample}".*)
    fi

    log "Running Prodigal"
    prodigal -i "$assembly" -q -a "${out_dir}/${sample}_proteins.faa" \
             -o "${out_dir}/${sample}_genes.gbk" \
             -d "${out_dir}/${sample}_nucleotide.ffn" -p meta

    log "Running DIAMOND"
    diamond blastp -d "$diamond_db" \
        -q "${out_dir}/${sample}_proteins.faa" \
        -o "${out_dir}/${sample}_diamond.txt" \
        -p "$threads" -k 1 -e "$evalue" \
        --id "$min_identity" \
        --query-cover "$min_query_cover" \
        $( [ -n "$min_score" ] && echo "--min-score $min_score" ) \
        $( [ -n "$diamond_mode" ] && echo "--mode $diamond_mode" ) \
        $diamond_extra

    log "Building Bowtie2 index"
    bowtie2-build "${out_dir}/${sample}_nucleotide.ffn" "${out_dir}/${sample}_bt2"

    log "Mapping reads back to genes"
    bowtie2 -x "${out_dir}/${sample}_bt2" \
        -1 "$read_file_1" -2 "$read_file_2" \
        -S "${out_dir}/${sample}.sam" -p "$threads"

    log "Calculating coverage"
    pileup.sh usejni=t in="${out_dir}/${sample}.sam" out="${out_dir}/${sample}.pileup"

    python "$script_dir/vis-scripts/gene_relative_abundance.py" \
        -p "${out_dir}/${sample}.pileup" -b "$sample" -o "$out_dir"

    python "$script_dir/vis-scripts/merge_blastp.py" \
        -b "${out_dir}/${sample}_diamond.txt" \
        -o "${out_dir}/${sample}_diamond_table.txt"

    python "$script_dir/vis-scripts/merge_abund_blastp.py" \
        -a "${out_dir}/${sample}.abundance" \
        -b "${out_dir}/${sample}_diamond_table.txt" \
        -o "${out_dir}/diamond_merged.txt"

    log "Cleaning temporary files"
    rm -f "${out_dir}/${sample}.sam" "${out_dir}/${sample}.pileup" \
          "${out_dir}/${sample}_diamond_table.txt"
done

log "Generating heatmaps"
python "$script_dir/vis-scripts/heatmap_plabase.py" \
    "${out_dir}/diamond_merged.txt" "$out_dir" \
    "$script_dir/database/pathways_plabase.txt" \
    "$script_dir/database/summary.txt"

log "meta_wf completed successfully"
