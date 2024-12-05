#!/usr/bin/env python3

import argparse
import os
import subprocess

# Get the path of the directory where the script is being executed
dir_path = os.path.dirname(os.path.realpath(__file__))

workflows = {
    "genome_wf": "workflows/genome_wf.sh",
    "metafast_wf": "workflows/metafast_wf.sh",
    "meta_wf": "workflows/meta_wf.sh",
}

def call_workflow(workflow, input_dir, output_dir, threads, assembly, mode):
    command = ["bash", os.path.join(dir_path, workflows[workflow]), "-i", input_dir, "-o", output_dir, "-t", str(threads)]
    if assembly:
        command.extend(["-a", assembly])
    if mode:
        command.extend(["--dmode", mode])
    subprocess.call(command)

def print_workflows():
    # Códigos ANSI para cores
    GREEN = "\033[32m"
    BLUE = "\033[36m"
    RED = "\033[31m"
    RESET = "\033[0m"

    print(f'''
                {BLUE}🌱🦠 PGPg_finder v1.0.1 🦠🌱{RESET}

PGPg_finder is a tool that offers different workflows for genome and metagenome annotation of Plant-Growth Promotion Genes.

{GREEN}Annotation of Genomes or Metagenome-Assembled Genomes:{RESET}

  genome_wf        ->  Performs genome annotation with PGPT-db database of PLaBAse.

{GREEN}Annotation of Short-Reads Metagenomes:{RESET}

  metafast_wf      ->  Rapid metagenome read alignment with mgPGPT-db database of PLaBAse
  meta_wf          ->  Accurate metagenome analysis with assembling and alignment with mgPGPT-db database of PLaBAse

{BLUE}Usage:{RESET}
  PGPg_finder <command> -h for command-specific help


PGPg_finder v1.1.0 | by Thierry Pellegrinetty <thierry.pellegrinetti@hotmail.com>
Check https://github.com/tpellegrinetti/PGPg_finder for updates



{BLUE}We encourage you to cite the PLaBAse and PGPg_finder:{RESET}

{GREEN}1){RESET} Patz S, Rauh M, Gautam A, Huson DH. mgPGPT: Metagenomic analysis of plant growth-promoting traits. (submitted, 2024, preprint)
{GREEN}2){RESET} Patz S, Gautam A, Becker M, Ruppel S, Rodríguez-Palenzuela P, Huson DH. PLaBAse: A comprehensive web resource for analyzing the plant growth-promoting potential of plant-associated bacteria. (submitted 2021, preprint)
{GREEN}3){RESET} Pellegrinetti, TA; Monteiro, G; Lemos, LN; RAC, Santos; Barros, A; Mendes, L. (2024) PGPg_finder: A Comprehensive and User-friendly Pipeline for Identifying Plant Growth-Promoting Genes in Genomic and Metagenomic Data. Rhizosphere.
''')


def print_workflow_help(workflow):
# Códigos ANSI para cores
    GREEN = "\033[32m"
    BLUE = "\033[36m"
    RED = "\033[31m"
    RESET = "\033[0m"
    if workflow == "genome_wf":
        print(f'''
{GREEN} 🧬 Genome workflow 🧬: {RESET}

Performs genome annotation with DIAMOND against the PGPT-db database of PLaBAse.

{GREEN} Required arguments: {RESET}
  -i <input_dir>:               Directory containing assemblies as nucleotide sequences (.fasta, .fa, .fna)
  -o <output_dir>:              Output directory to save results

{BLUE} Optional arguments: {RESET}
  -t <threads>:                 Number of threads used for analysis (default: 1)
  --dmode <diamond-mode>:       DIAMOND mode for alignment (e.g., faster, sensitive, very-sensitive)
  --piden <min_identity>:       Minimum identity for alignment (%) (default: 30)
  --qcov <query_cover>:         Minimum query coverage (%) (default: 30)
  --bitscore <bit_score>:       Minimum bit score to report alignments
  --extra <extra>:              Extra arguments for DIAMOND
  --evalue <evalue>:            Maximum e-value to report alignments"
  
  {GREEN}Usage:{RESET}
  PGPg_finder -w genome_wf -i input_folder/ -o output_folder/ -t 12
''')
    elif workflow == "metafast_wf":
        print(f'''
{GREEN} 🧬 Metagenome Read-Based workflow 🧬: {RESET}

Performs metagenome annotation using reads with DIAMOND against the PGPT-db database of PLaBAse.

{GREEN} Required arguments: {RESET}
  -i <input_dir>:               Directory containing metagenomes in FASTQ format (.fastq, .fq, .fastq.gz, .fq.gz)
  -o <output_dir>:              Output directory to save results

{BLUE} Optional arguments: {RESET}
  -t <threads>:                 Number of threads used for analysis (default: 1)
  --dmode <diamond-mode>:       DIAMOND mode for alignment (e.g., fast, sensitive, very-sensitive)
  --piden <min_identity>:       Minimum identity for alignment (%) (default: 30)
  --qcov <query_cover>:         Minimum query coverage (%) (default: 30)
  --bitscore <bit_score>:       Minimum bit score to report alignments
  --extra <extra>:              Extra arguments for DIAMOND
  --evalue <evalue>:            Maximum e-value to report alignments"
  
  {GREEN}Usage:{RESET}
  PGPg_finder -w metafast_wf -i input_folder/ -o output_folder/ -t 12
''')
    elif workflow == "meta_wf":
        print(f'''
{GREEN} 🧬 Metagenome Assembly-Based workflow 🧬: {RESET}
 
Performs metagenome annotation using contigs with DIAMOND against the PGPT-db database of PLaBAse.

{GREEN} Required arguments: {RESET}
  -i <input_dir>:              Directory containing metagenomes in FASTQ format (.fastq, .fq, .fastq.gz, .fq.gz)
  -o <output_dir>:             Output directory to save results

{BLUE} Optional arguments: {RESET}
  -t <threads>:                Number of threads used for analysis (default: 1)
  -a <assembly>:               Pre-assembled contigs provided by the user (in FASTA format)
  --dmode <diamond-mode>:        DIAMOND mode for alignment (e.g., fast, sensitive, very-sensitive)
  --piden <min_identity>:       Minimum identity for alignment (%) (default: 30)
  --qcov <query_cover>:         Minimum query coverage (%) (default: 30)
  --bitscore <bit_score>:       Minimum bit score to report alignments
  --extra <extra>:              Extra arguments for DIAMOND
  --evalue <evalue>:            Maximum e-value to report alignments"
  
  {GREEN}Usage:{RESET}
  PGPg_finder -w meta_wf -i input_folder/ -o output_folder/ -t 12

 
''')
    else:
        print("Invalid workflow specified.")

def main():
    parser = argparse.ArgumentParser(
        description=(
            "PGPg_finder v1.1.0 | by Thierry Pellegrinetty <thierry.pellegrinetti@hotmail.com>\n"
            "Check https://github.com/tpellegrinetti/PGPg_finder for updates.\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=False
    )

    parser.add_argument('-w', '--workflow', choices=workflows.keys(), help="Choose the workflow to execute.")
    parser.add_argument('--list-workflows', action='store_true', help="List all available workflows.")
    parser.add_argument('-h', '--help', action='store_true', help="Show this help message and exit.")

    args, remaining_args = parser.parse_known_args()

    if args.help:
        if args.workflow:
            print_workflow_help(args.workflow)
        else:
            print_workflows()
        return

    if args.list_workflows:
        print_workflows()
        return

    if args.workflow:
        print_workflow_help(args.workflow)
        subparser = argparse.ArgumentParser()
        subparser.add_argument('-i', '--input', required=True, help='Path to the input directory.')
        subparser.add_argument('-o', '--output', required=True, help='Path to the output directory.')
        subparser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use. Default is 1.')
        if args.workflow == "genome_wf":
            subparser.add_argument('--dmode', '--mode', help='DIAMOND mode for alignment.')
            subparser.add_argument('--piden', '--min_identity', type=float, help='Minimum identity for DIAMOND alignment.')
            subparser.add_argument('--qcov', '--query_cover', type=float, help='Minimum query coverage for DIAMOND alignment.')
            subparser.add_argument('--bitscore', '--min_score', type=float, help='Minimum bit score to report alignments.')
            subparser.add_argument('--extra', '--diamond_extra', help='Additional arguments for DIAMOND.')
            subparser.add_argument('--evalue', '--evalue', help='Maximum e-value to report alignments.')
        elif args.workflow == "metafast_wf":
            subparser.add_argument('--dmode', '--mode', help='DIAMOND mode for alignment.')
            subparser.add_argument('--piden', '--min_identity', type=float, help='Minimum identity for DIAMOND alignment.')
            subparser.add_argument('--qcov', '--query_cover', type=float, help='Minimum query coverage for DIAMOND alignment.')
            subparser.add_argument('--bitscore', '--min_score', type=float, help='Minimum bit score to report alignments.')
            subparser.add_argument('--extra', '--diamond_extra', help='Additional arguments for DIAMOND.')
            subparser.add_argument('--evalue', '--evalue', help='Maximum e-value to report alignments.')
        elif args.workflow == "meta_wf":
            subparser.add_argument('-a', '--assembly', help='Path to the assembly file.')
            subparser.add_argument('--dmode', '--mode', help='DIAMOND mode for alignment.')
            subparser.add_argument('--piden', '--min_identity', type=float, help='Minimum identity for DIAMOND alignment.')
            subparser.add_argument('--qcov', '--query_cover', type=float, help='Minimum query coverage for DIAMOND alignment.')
            subparser.add_argument('--bitscore', '--min_score', type=float, help='Minimum bit score to report alignments.')
            subparser.add_argument('--extra', '--diamond_extra', help='Additional arguments for DIAMOND.')
            subparser.add_argument('--evalue', '--evalue', help='Maximum e-value to report alignments.')
        parsed_args = subparser.parse_args(remaining_args)
        call_workflow(args.workflow, parsed_args.input, parsed_args.output, parsed_args.threads, getattr(parsed_args, 'assembly', None), getattr(parsed_args, 'mode', None))

if __name__ == "__main__":
    main()
