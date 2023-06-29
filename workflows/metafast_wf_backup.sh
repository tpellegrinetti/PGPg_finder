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

# Path to the DIAMOND database
diamond_db="database/metagenome.dmnd"

# Verify if required arguments were provided
if [ -z "$reads_dir" ] || [ -z "$out_dir" ] || [ -z "$threads" ] || [ -z "$diamond_db" ]; then
    echo "Missing arguments. See usage below:"
    display_help
    exit 1
fi

# Create output directory if it doesn't exist
if [ ! -d "$out_dir" ]; then
    mkdir -p "$out_dir"
fi

# Create a file to store gene counts
echo -e "Sample\tGene\tCount" > ${out_dir}/gene_counts.txt

# Loop to iterate over each pair of read files
for reads in $reads_dir/*_1.*; do
    sample=$(basename "$reads")
    sample=${sample%_1.*}
    echo "Processing sample $sample..."
    
    read_file_1=$reads
    read_file_2=${reads_dir}/${sample}_2.*

    # Merge paired reads with PEAR
    pear -f ${read_file_1} -r ${read_file_2} -o ${out_dir}/${sample}_merged -j ${threads}

    # Run prodigal to generate protein sequences
    prodigal -i ${out_dir}/${sample}_merged.assembled.fastq -a ${out_dir}/${sample}_proteins.fa -p meta

    # Run DIAMOND
    diamond blastx -d ${diamond_db} -q ${out_dir}/${sample}_merged.assembled.fastq -o ${out_dir}/${sample}_diamond.txt -k 1 -e 0.0001 -p ${threads}

    # Parse DIAMOND output and get gene counts
    awk '{print $2}' ${out_dir}/${sample}_diamond.txt | sort | uniq -c | while read count gene; do
        echo -e "${sample}\t${gene}\t${count}" >> ${out_dir}/gene_counts.txt
    done

    # Clean up
    rm -f ${out_dir}/${sample}_merged*
    rm -f ${out_dir}/${sample}_proteins.fa
done

echo "Gene search is completed. Check ${out_dir}/gene_counts.txt for the results."

# Python script to generate the heatmaps
python vis-scripts/heatmap_plabase.py ${out_dir}/gene_counts.txt ${out_dir} database/pathways_plabase.txt

