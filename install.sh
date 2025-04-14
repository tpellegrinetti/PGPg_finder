#!/bin/bash

# Welcome message
echo "Hello, welcome to the PGPB_finder!"
echo "This pipeline finds plant growth-promoting genes in your metagenomes or genomes."
echo "This tool was developed by Dr. Thierry Pellegrinetti."
echo "If this tool proves useful to you, please cite: PGPg_finder."
echo ""

# Create a new conda environment
echo "First, we need to create a separate environment for PGPg_finder..."
conda create -n PGPg_finder -y
source activate PGPg_finder
echo "Environment 'PGPg_finder' is now active!"
echo ""

# Install necessary programs
echo "Now, let's install the dependencies..."
conda install -c bioconda prodigal diamond megahit bowtie2 samtools gawk pear trimmomatic -y
conda install pandas seaborn matplotlib
echo "Dependencies installed successfully!"
echo ""

# Verify that all programs were installed correctly
echo "Verifying if all dependencies were installed correctly..."
prokka --version
megahit --version
bwa
samtools --version
awk --version
echo "All dependencies are installed and working correctly!"
echo ""

echo "To use PGPg_finder, activate the conda environment using:"
echo "    conda activate PGPg_finder"
echo "To deactivate the environment when you're done, use:"
echo "    conda deactivate"
echo ""

echo "Thank you for installing PGPg_finder. Happy analyzing!"

echo "Lets starting the download of databases"
echo "The first is the genome PGPTdb"

# Download databases
cd database/
echo "Downloading genome database..."
wget https://plabase.cs.uni-tuebingen.de/pb/tools/PGPTblhm/data/factors/PGPT_BASE_nr_Aug2021n_ul_1.fasta.gz

echo "The second is the mgPGPT-db database"
echo "Downloading metagenome database..."
wget https://plabase.cs.uni-tuebingen.de/pb/tools/PGPTblhm/data/factors/mgPGPT/mgPGPT-db_Feb2022_ul_dwnld.fasta.gz

echo "Databases downloaded successfully."

echo "Decompressing database"

gzip -d *.gz

###Making diamond databases

echo "Making diamond db for genome and metagenome"

diamond makedb --in PGPT_BASE_nr_Aug2021n_ul_1.fasta --db genome
diamond makedb --in mgPGPT-db_Feb2022_ul_dwnld.fasta --db metagenome

echo "remove fasta files"

rm PGPT_BASE_nr_Aug2021n_ul_1.fasta
rm mgPGPT-db_Feb2022_ul_dwnld.fasta
