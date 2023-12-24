# PGPg_finder

![My Image](static/logo.png)
PGPg Finder is a versatile command-line tool designed for the annotation of Plant-Growth Promotion Traits (PGPT) in both genomes and metagenomes. Rooted in the comprehensive PLant-associated BActeria web resource database (PLaBAse), this tool boasts flexibility and user-friendliness. PGPg Finder offers a range of pre-configured workflows, enabling users to select the tools and processes best suited to their specific research needs.


# Features

The pipeline allows you to choose from several workflows based on the tools and steps you need. The workflows currently available are:

## Genome Analysis

genome_wf: A genome annotation pipeline that uses Prodigal for gene prediction and DiAMOND against PLaBAse.

## Metagenome Analysis

meta_wf: A metagenome workflow with megahit assembly, Prodigal gene prediction, DIAMOND annotation, Bowtie2 mapping, samtools quantification.  

metafast_wf: A metagenome workflow with PEAR assembly and direct DIAMOND annotation.


Each workflow is designed to be run with a specific set of input files and produces a set of output files. The tool also supports running on multiple threads to improve performance.


# Installation

To use the tool, you just need to clone the repository and choose between conda or nextflow installation (see below).

## Conda (recomended)

We assume users have conda installed. If not, please follow the instructions on the [conda website](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

All environment dependencies are listed in the [enviroment_pgpg_finder.yml file](installation/enviroment_pgpg_finder.yml). To install the dependencies and create the conda environment ("PGPb_finder"), run the following commands:

```bash
git clone https://github.com/tpellegrinetti/PGPg_finder/
cd PGPg_finder/installation
conda env create -f enviroment_pgpg_finder.yml
conda activate PGPg_finder
```

TODO: Testing 

### Downloading the database

To download the database, you just need to run command below.

TODO: add database versions
TODO: rethink how to download the database. Adding hard links is usually not a good idea because if the link does not work, the user will not be able to download the database. Also, users may want to download a different version of the database.

```bash
git clone https://github.com/tpellegrinetti/PGPg_finder/
cd PGPg_finder/installation
bash download_database.sh
```

## Nextflow (optional)

TODO: add nextflow installation instructions (including database download and tests)

# Usage:

TODO: adapt usage to the new version. Also, usage will change in the Nextflow pipeline.

TODO: add instructions on how to add set environment variables (e.g. PGPg_finder scripts and database paths)

You can run the PGPg_finder tool using the following command:

`python PGPb_finder.py -w workflow -i input_directory -o output_directory -t threads`

the -a argument is optional and works if you want to provide assembly files


1) Running PGPb_finder with genomes:

`python PGPb_finder.py -w genome_wf -i /path/to/fasta/folder/ -o /path/to/your/desired/out/ -t 12`


2) Running PGPb_finder with metagenomes:
 
a) Fast way (less accurated)

`python PGPb_finder.py -w metafast_wf -i /path/to/fasta/folder -o /path/to/your/desired/out/ -t 12`

b) Slow way (more accurated) 

`python PGPb_finder.py -w meta_wf -i /path/to/fasta/folder -o /path/to/your/desired/out/ -t 12`

Here you can provide your metagenome assemblies with -a option

# Citation

If PGPb_finder was useful for you please cite us:

Pellegrinetti, TA; Monteiro, G; Lemos, LN; Tsai, SM; Mendes, L. PGP_finder: PGPg_finder: A Comprehensive and User-friendly Pipeline for Identifying Plant Growth-Promoting Genes in Genomic and Metagenomic Data. 

This pipeline was developed based on a curated database called "PLant-associated BActeria web resource (PLaBAse)". Access PLaBAse website: https://plabase.cs.uni-tuebingen.de/

We encourage you to cite the PLaBAse:

Patz S, Gautam A, Becker M, Ruppel S, Rodríguez-Palenzuela P, Huson DH. PLaBAse: A comprehensive web resource for analyzing the plant growth-promoting potential of plant-associated bacteria. (submitted 2021, meanwhile you can read the preprint)