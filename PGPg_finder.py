#!/usr/bin/env python3

import argparse
import os
import subprocess

# Caminho do diretório onde este script está localizado
dir_path = os.path.dirname(os.path.realpath(__file__))

# Dicionário com os caminhos dos workflows
workflows = {
    "genome_wf": "workflows/genome_wf.sh",
    "metafast_wf": "workflows/metafast_wf.sh",
    "meta_wf": "workflows/meta_wf.sh",
}

def call_workflow(workflow, args):
    command = ["bash", os.path.join(dir_path, workflows[workflow]),
               "-i", args.input, "-o", args.output, "-t", str(args.threads)]

    if hasattr(args, "assembly") and args.assembly:
        command.extend(["-a", args.assembly])
    if args.dmode:
        command.extend(["--dmode", args.dmode])
    if args.piden:
        command.extend(["--piden", str(args.piden)])
    if args.qcov:
        command.extend(["--qcov", str(args.qcov)])
    if args.bitscore:
        command.extend(["--bitscore", str(args.bitscore)])
    if args.extra:
        command.extend(["--extra", args.extra])
    if args.evalue:
        command.extend(["--evalue", str(args.evalue)])

    subprocess.call(command)

def print_workflows():
    GREEN = "\033[32m"
    BLUE = "\033[36m"
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
  PGPg_finder -w (genome_wf or metafast_wf or meta_wf) -h for command-specific help

PGPg_finder v1.1.0 | by Thierry Pellegrinetty <thierry.pellegrinetti@hotmail.com>
Check https://github.com/tpellegrinetti/PGPg_finder for updates

{BLUE}We encourage you to cite the PLaBAse and PGPg_finder:{RESET}
{GREEN}1){RESET} Patz S et al. mgPGPT: Metagenomic analysis of plant growth-promoting traits. (submitted, 2024)
{GREEN}2){RESET} Patz S et al. PLaBAse: A comprehensive web resource... (submitted 2021)
{GREEN}3){RESET} Pellegrinetti TA et al. (2024) PGPg_finder: A Comprehensive Pipeline... Rhizosphere.
''')

def print_workflow_help(workflow):
    GREEN = "\033[32m"
    BLUE = "\033[36m"
    RESET = "\033[0m"

    if workflow == "genome_wf":
        print(f'''
{GREEN} 🧬 Genome workflow 🧬: {RESET}

Performs genome annotation with DIAMOND against the PGPT-db database of PLaBAse.

{GREEN} Required arguments: {RESET}
  -i <input_dir>         Directory with assemblies (.fasta, .fa, .fna)
  -o <output_dir>        Output directory

{BLUE} Optional arguments: {RESET}
  -t <threads>           Number of threads (default: 1)
  --dmode                DIAMOND mode (fast, sensitive, very-sensitive)
  --piden                Minimum identity (%) (default: 30)
  --qcov                 Minimum query coverage (%) (default: 30)
  --bitscore             Minimum bit score
  --evalue               Max e-value (default: 1e-5)
  --extra                Extra DIAMOND options

{GREEN}Usage:{RESET}
  PGPg_finder -w genome_wf -i input_dir -o output_dir -t 12
''')
    elif workflow == "metafast_wf":
        print(f'''
{GREEN} 🧬 Metagenome Read-Based workflow 🧬: {RESET}

Annotates short-read metagenomes directly using DIAMOND and PGPT-db.

{GREEN} Required arguments: {RESET}
  -i <input_dir>         Directory with reads (.fastq or .fastq.gz)
  -o <output_dir>        Output directory

{BLUE} Optional arguments: {RESET}
  -t <threads>           Number of threads (default: 1)
  --dmode                DIAMOND mode (fast, sensitive, very-sensitive)
  --piden                Minimum identity (%)
  --qcov                 Minimum query coverage (%)
  --bitscore             Minimum bit score
  --evalue               Max e-value
  --extra                Extra DIAMOND options

{GREEN}Usage:{RESET}
  PGPg_finder -w metafast_wf -i input_dir -o output_dir -t 12
''')
    elif workflow == "meta_wf":
        print(f'''
{GREEN} 🧬 Metagenome Assembly-Based workflow 🧬: {RESET}

Annotates assembled contigs using DIAMOND and PGPT-db.

{GREEN} Required arguments: {RESET}
  -i <input_dir>         Directory with reads or assemblies
  -o <output_dir>        Output directory

{BLUE} Optional arguments: {RESET}
  -t <threads>           Number of threads (default: 1)
  -a <assembly>          Use pre-assembled contigs (FASTA)
  --dmode                DIAMOND mode (fast, sensitive, very-sensitive)
  --piden                Minimum identity (%)
  --qcov                 Minimum query coverage (%)
  --bitscore             Minimum bit score
  --evalue               Max e-value
  --extra                Extra DIAMOND options

{GREEN}Usage:{RESET}
  PGPg_finder -w meta_wf -i input_dir -o output_dir -t 12
''')
    else:
        print("Invalid workflow specified.")

def main():
    parser = argparse.ArgumentParser(
        description=(
            "PGPg_finder v1.1.0 | by Thierry Pellegrinetti\n"
            "Check https://github.com/tpellegrinetti/PGPg_finder for updates."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=False
    )

    parser.add_argument('-w', '--workflow', choices=workflows.keys())
    parser.add_argument('--list-workflows', action='store_true')
    parser.add_argument('-h', '--help', action='store_true')

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
        subparser.add_argument('-i', '--input', required=True)
        subparser.add_argument('-o', '--output', required=True)
        subparser.add_argument('-t', '--threads', type=int, default=1)
        subparser.add_argument('--dmode')
        subparser.add_argument('--piden', type=float)
        subparser.add_argument('--qcov', type=float)
        subparser.add_argument('--bitscore', type=float)
        subparser.add_argument('--evalue')
        subparser.add_argument('--extra')

        if args.workflow == "meta_wf":
            subparser.add_argument('-a', '--assembly')

        parsed_args = subparser.parse_args(remaining_args)
        call_workflow(args.workflow, parsed_args)

if __name__ == "__main__":
    main()
