#!/usr/bin/env python3

import argparse
import subprocess
import os  # Import the os library

# Get the path of the directory where the script is being executed
dir_path = os.path.dirname(os.path.realpath(__file__))

def call_workflow(args):
    command = ["bash", os.path.join(dir_path, f"workflows/{args.workflow}.sh"), "-i", args.input, "-o", args.output, "-t", str(args.threads)]
    if args.assembly:
        command.extend(["-a", args.assembly])
    if args.mode:
        command.extend(["-m", args.mode])
    subprocess.call(command)

def main(args):
    call_workflow(args)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='PGPg_finder is a tool that offers different workflows for genome and metagenome annotation. Select the workflow based on the tools and steps needed.', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-w', '--workflow', required=True, help='''Choose the workflow:\n
                        genome_wf - Genome annotation with diamond against PGPT-db database of PLaBAse.\n
                        metafast_wf - Fast metagenome direct read annotation with diamond against mgPGPT-db database of PLaBAse.\n
                        meta_wf - Metagenome assembly annotation (diamond) and coverage (Bowtie2, samtools bedtools) with mgPGPT-db database of PLaBAse.''')
                        
    parser.add_argument('-i', '--input', required=True, help='Path to the input directory.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output directory.')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use. Default is 1.')
    parser.add_argument('-a', '--assembly', help='Path to the assembly file. Use if you prefer to provide your own assembly.')
    parser.add_argument('-m', '--mode', help = '(faster, fast, mid-sensitive, sensitive, more-sensitive, very-sensitive, ultra-sensitive) - [defaul=fast]', default=None)

    args = parser.parse_args()
    
    main(args)

