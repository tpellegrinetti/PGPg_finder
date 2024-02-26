#!/bin/bash

# Function to display help
display_help() {
    echo "Usage: $0 -i <input_directory> -o <output_directory> -t <number_of_threads> [-m <mode>]"
    echo
    echo "This script performs gene search with DIAMOND on a collection of metagenomes,"
    echo "after quality trimming with Trimmomatic, merging paired reads with PEAR, and gene prediction with Prodigal."
    echo "Results are saved to the output directory."
    echo
    echo "Arguments:"
    echo "  -i   Path to the directory containing read files (_1.fq and _2.fq)."
    echo "  -o   Path to the directory where results will be saved."
    echo "  -t   Number of threads to be used by the tools."
    echo "  -m   DIAMOND mode (optional)."
    echo
}

# Function for logging
log() {
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] $1" | tee -a "$log_file"
}

# Check if no argument was provided
if [ $# -eq 0 ]; then
    display_help
    exit 1
fi

# Parse command-line arguments
diamond_mode=""
while getopts i:o:t:m:h opt; do
  case $opt in
    i) reads_dir=$OPTARG ;;
    o) out_dir=$OPTARG ;;
    t) threads=$OPTARG ;;
    m) diamond_mode=$OPTARG ;;
    h) display_help
       exit 0 ;;
    *) display_help
       exit 1 ;;
  esac
done

# Obtenha o caminho do diretório onde o script está sendo executado
script_dir=$(dirname "$(dirname "$(readlink -f "$0")")")

# Path to the DIAMOND database
diamond_db="$script_dir/database/metagenome.dmnd"

# Verify if required arguments were provided
if [ -z "$reads_dir" ] || [ -z "$out_dir" ] || [ -z "$threads" ]; then
    log "Missing arguments. See usage below:"
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

# Create a file to store gene counts
gene_counts_file="${out_dir}/gene_counts.txt"
echo -e "Sample\tID\tCount" > "$gene_counts_file"
log "Created gene counts file: $gene_counts_file"

# Loop to iterate over each pair of read files
for reads_1 in $reads_dir/*_1.*; do
    sample=$(basename "$reads_1")
    sample=${sample%_1.*}
    log "Processing sample $sample..."
    
    read_file_1=$reads_1
    read_file_2=${reads_1/_1./_2.}
    trimmed_file_1="${out_dir}/${sample}_trimmed_1.fq"
    trimmed_file_2="${out_dir}/${sample}_trimmed_2.fq"
    trimmed_se="${out_dir}/${sample}_trimmed_se.fq"

    # Quality trimming with Trimmomatic
    log "Running Trimmomatic for quality trimming"
    if [ -f "$read_file_2" ]; then
        # Paired-end reads
        trimmomatic PE -threads $threads $read_file_1 $read_file_2 $trimmed_file_1 $trimmed_se $trimmed_file_2 $trimmed_se SLIDINGWINDOW:4:20 MINLEN:36
    else
        # Single-end reads
        trimmomatic SE -threads $threads $read_file_1 $trimmed_file_1 SLIDINGWINDOW:4:20 MINLEN:36
    fi

    # Check if diamond_mode is set and add it to the DIAMOND command
    if [ ! -z "$diamond_mode" ]; then
        diamond_mode_flag="--mode $diamond_mode"
    else
        diamond_mode_flag=""
    fi

    if [ -f "$read_file_2" ]; then
        # Paired-end reads
        # Merge paired reads with PEAR
        pear -f $trimmed_file_1 -r $trimmed_file_2 -o ${out_dir}/${sample}_merged -j $threads -q 20
        log "Merged paired reads for sample $sample"
    
        # Run DIAMOND
        diamond blastx -d $diamond_db -q ${out_dir}/${sample}_merged.assembled.fastq -o ${out_dir}/${sample}_diamond.txt -k 1 -e 0.0001 -p $threads $diamond_mode_flag
        log "Completed DIAMOND search for sample $sample"

        # Parse DIAMOND output and get gene counts
        awk '{print $2}' ${out_dir}/${sample}_diamond.txt | sort | uniq -c | while read count id; do
            echo -e "${sample}\t${id}\t${count}" >> $gene_counts_file
        done
        log "Generated gene counts for sample $sample"

        # Clean up
        rm -f ${out_dir}/${sample}_merged*
    else
        # Single-end reads
        # Run DIAMOND directly on the single-end reads
        diamond blastx -d $diamond_db -q $trimmed_file_1 -o ${out_dir}/${sample}_diamond.txt -k 1 -e 0.0001 -p $threads $diamond_mode_flag
        log "Completed DIAMOND search for sample $sample"

        # Parse DIAMOND output and get gene counts
        awk '{print $2}' ${out_dir}/${sample}_diamond.txt | sort | uniq -c | while read count id; do
            echo -e "${sample}\t${id}\t${count}" >> $gene_counts_file
        done
        log "Generated gene counts for sample $sample"
    fi
done

log "Gene search is completed. Check ${gene_counts_file} for the results."

# Python script to generate the heatmaps
heatmap_script="$script_dir/vis-scripts/heatmap_plabase.py"
python "$heatmap_script" "${gene_counts_file}" "${out_dir}" "$script_dir/database/pathways_plabase.txt" "$script_dir/database/summary.txt"
log "Generated heatmaps"

