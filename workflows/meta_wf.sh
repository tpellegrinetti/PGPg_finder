#!/bin/bash

# Function to display help
display_help() {
    echo "Usage: $0 -i <input_directory> -o <output_directory> -t <number_of_threads> [-a <assembly_directory>] [-m <mode>]"
    echo
    echo "This script performs gene search with HMM profiles on a collection of metagenomes,"
    echo "after quality trimming with Trimmomatic, generates a gene count table, and creates gene presence/absence heatmaps."
    echo "Results are saved to the output directory."
    echo
    echo "Arguments:"
    echo "  -i   Path to the directory containing read files (_1.fq.gz and _2.fq.gz)."
    echo "  -o   Path to the directory where results will be saved."
    echo "  -t   Number of threads to be used by the tools."
    echo "  -a   (Optional) Path to the directory containing assembly files (.fasta). If provided, assembly step is skipped."
    echo "  -m   (Optional) DIAMOND mode for sequence search."
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
while getopts i:o:t:a:m:h opt; do
  case $opt in
    i) reads_dir=$OPTARG ;;
    o) out_dir=$OPTARG ;;
    t) threads=$OPTARG ;;
    a) assembly_dir=$OPTARG ;;
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

# Check if diamond_mode is set and add it to the DIAMOND command
if [ ! -z "$diamond_mode" ]; then
    diamond_mode_flag="--mode $diamond_mode"
else
    diamond_mode_flag=""
fi

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

    # Quality trimming with Trimmomatic
    log "Running Trimmomatic for quality trimming"
    if [ -f "$read_file_2" ]; then
        # Paired-end reads
        trimmomatic PE -threads $threads $read_file_1 $read_file_2 $trimmed_file_1 $trimmed_se $trimmed_file_2 $trimmed_se SLIDINGWINDOW:4:20 MINLEN:36
    else
        # Single-end reads
        trimmomatic SE -threads $threads $read_file_1 $trimmed_file_1 SLIDINGWINDOW:4:20 MINLEN:36
    fi

    if [ -z "$assembly_dir" ]; then
        # If no assembly file is provided, assemble the genome with MEGAHIT
        log "Assembling metagenome for sample $sample"
        megahit -1 $trimmed_file_1 -2 $trimmed_file_2 -t $threads -o ${out_dir}/${sample}_assembly
        assembly=${out_dir}/${sample}_assembly/final.contigs.fa
    else
        # Use the provided assembly
        log "Using provided assembly for sample $sample"
        assembly=$(ls ${assembly_dir}/*.fasta)
    fi
    
    # Run prokka to generate protein and annotate sequences
    log "Generating protein sequences for sample $sample"
    prodigal -i ${assembly} -q -a ${out_dir}/${sample}_proteins.faa -o ${out_dir}/${sample}_genes.gbk -d ${out_dir}/${sample}_nucleotide.ffn -p meta
   
    # Run DIAMOND
    log "Running DIAMOND search for sample $sample"
    diamond blastp -d ${diamond_db} -q ${out_dir}/${sample}_proteins.faa -o ${out_dir}/${sample}_diamond.txt -k 1 -e 0.0001 -p ${threads} $diamond_mode_flag
   
    # Build bowtie2 index
    log "Building bowtie2 index for sample $sample"
    bowtie2-build ${out_dir}/${sample}_nucleotide.ffn ${out_dir}/${sample}_bowtie2_index

    # Map reads back to the assembly with bowtie2
    log "Mapping reads to the assembly for sample $sample"
    bowtie2 -x ${out_dir}/${sample}_bowtie2_index -1 ${read_file_1} -2 ${read_file_2} -S ${out_dir}/${sample}_mapped_reads.sam -p ${threads}

    #Calculate gene abundance with BBMap pileup
    log "Calculate coverage for sample $sample"
    pileup.sh usejni=t in=${out_dir}/${sample}_mapped_reads.sam out=${out_dir}/${sample}.pileup

    #Calculate gene abundance coverage
    log "Generate contig coverage $sample"
    python "$script_dir/vis-scripts/gene_relative_abundance.py" -p ${out_dir}/${sample}.pileup -b ${sample} -o ${out_dir}
    
    #Get diamond out results with contigs
    log "Generate contig and gene for $sample"
    python "$script_dir/vis-scripts/merge_blastp.py" -b ${out_dir}/${sample}_diamond.txt -o ${out_dir}/${sample}_diamond_table.txt

    #Merge table to obtain Sample - ID - Count
    log "Merge contig and coverage with blastp result $sample"
    python "$script_dir/vis-scripts/merge_abund_blastp.py" -a ${out_dir}/${sample}.abundance -b ${out_dir}/${sample}_diamond_table.txt -o ${out_dir}/diamond_merged.txt

    # Moving files
    log "Moving Files"
    if [ ! -d "${out_dir}/prodigal" ]; then
        mkdir -p "${out_dir}/prodigal"
        log "Created output directory: ${out_dir}/prodigal"
    fi

    mv ${out_dir}/${sample}_proteins.fa ${out_dir}/prodigal/${sample}_proteins.fa
    log "${sample}_proteins.fa moved to ${out_dir}/prodigal"
    
    if [ ! -d "${out_dir}/coverage" ]; then
        mkdir -p "${out_dir}/coverage"
        log "Created output directory: ${out_dir}/coverage"
    fi
    
    mv ${out_dir}/${sample}_pileup.txt ${out_dir}/coverage/${sample}_pileup.txt

    rm ${out_dir}/${sample}_bowtie2_index* 
    rm -f ${out_dir}/${sample}_bowtie2_index*
    rm -f ${out_dir}/${sample}_mapped_reads.sam
    rm -f ${out_dir}/${sample}_mapped_reads.bam
    rm -f ${out_dir}/${sample}_mapped_reads_sorted.bam
    rm -f ${out_dir}/${sample}.abundance
    rm -f ${out_dir}/${sample}_diamond_table.txt
    
    log "All mapping files were removed"
    
    if [ ! -d "${out_dir}/assembly" ]; then
        mkdir -p "${out_dir}/assembly"
        log "Created output directory: ${out_dir}/assembly"
    fi
    
    mv ${out_dir}/${sample}_assembly ${out_dir}/assembly/${sample}_assembly
    log "${sample}_assembly was moved to ${out_dir}/assembly"
     
done

log "Gene search is completed. Check ${gene_counts_file} for the results."

# Python script to generate the heatmaps
log "Generating heatmaps"
python "$script_dir/vis-scripts/heatmap_plabase.py" ${out_dir}/diamond_merged.txt ${out_dir} "$script_dir/database/pathways_plabase.txt" "$script_dir/database/summary.txt"
log "Generated heatmaps"

