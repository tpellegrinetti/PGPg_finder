#!/bin/bash

# Function to display help
display_help() {
    echo "Usage: $0 -i <input_directory> -o <output_directory> -t <number_of_threads>"
    echo
    echo "This script performs gene search with DIAMOND on a collection of metagenomes,"
    echo "after merging paired reads with PEAR and gene prediction with Prodigal."
    echo "Results are saved to the output directory."
    echo
    echo "Arguments:"
    echo "  -i   Path to the directory containing read files (_1.fq and _2.fq)."
    echo "  -o   Path to the directory where results will be saved."
    echo "  -t   Number of threads to be used by the tools."
    echo
}

# Function for logging
log() {
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] $1"
}

# Parse command-line arguments
while getopts i:o:t:a:h opt; do
  case $opt in
    i) reads_dir=$OPTARG ;;
    o) out_dir=$OPTARG ;;
    t) threads=$OPTARG ;;
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
if [ -z "$reads_dir" ] || [ -z "$out_dir" ] || [ -z "$threads" ] || [ -z "$diamond_db" ]; then
    log "Missing arguments. See usage below:"
    display_help
    exit 1
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
for reads_1 in $reads_dir/*_1.*; do
    sample=$(basename "$reads_1")
    sample=${sample%_1.*}
    log "Processing sample $sample..."
    
    read_file_1=$reads_1
    read_file_2=${reads_1/_1./_2.}

    if [ -f "$read_file_2" ]; then
        # Paired-end reads
        # Merge paired reads with PEAR
        pear -f ${read_file_1} -r ${read_file_2} -o ${out_dir}/${sample}_merged -j ${threads} -q 20
        log "Merged paired reads for sample $sample"

        # Run DIAMOND
        diamond blastx -d ${diamond_db} -q ${out_dir}/${sample}_merged.assembled.fastq -o ${out_dir}/${sample}_diamond.txt -k 1 -e 0.0001 -p ${threads}
        log "Completed DIAMOND search for sample $sample"

        # Parse DIAMOND output and get gene counts
        awk '{print $2}' ${out_dir}/${sample}_diamond.txt | sort | uniq -c | while read count id; do
            echo -e "${sample}\t${id}\t${count}" >> ${gene_counts_file}
        done
        log "Generated gene counts for sample $sample"

        # Clean up
        rm -f ${out_dir}/${sample}_merged*

    else
        # Single-end reads
        # Run DIAMOND directly on the single-end reads
        diamond blastx -d ${diamond_db} -q ${read_file_1} -o ${out_dir}/${sample}_diamond.txt -k 1 -e 0.0001 -p ${threads}
        log "Completed DIAMOND search for sample $sample"

        # Parse DIAMOND output and get gene counts
        awk '{print $2}' ${out_dir}/${sample}_diamond.txt | sort | uniq -c | while read count id; do
            echo -e "${sample}\t${id}\t${count}" >> ${gene_counts_file}
        done
        log "Generated gene counts for sample $sample"
    fi
done

log "Gene search is completed. Check ${gene_counts_file} for the results."

# Python script to generate the heatmaps
heatmap_script="$script_dir/vis-scripts/heatmap_plabase.py"
python "$heatmap_script" "${gene_counts_file}" "${out_dir}" "$script_dir/database/pathways_plabase.txt" "$script_dir/database/summary.txt"
log "Generated heatmaps"

