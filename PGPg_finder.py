#!/usr/bin/env python3

import argparse
import subprocess

def call_workflow(args):
    if args.workflow == 'genome_wf':
        command = ["bash", "workflows/genome_wf.sh", "-i", args.input, "-o", args.output, "-t", str(args.threads)]
        if args.assembly:
            command.extend(["-a", args.assembly])
        subprocess.call(command)
    elif args.workflow == 'metafast_wf':
        command = ["bash", "workflows/metafast_wf.sh", "-i", args.input, "-o", args.output, "-t", str(args.threads)]
        if args.assembly:
            command.extend(["-a", args.assembly])
        subprocess.call(command)
    elif args.workflow == 'meta_wf':
        command = ["bash", "workflows/meta_wf.sh", "-i", args.input, "-o", args.output, "-t", str(args.threads)]
        if args.assembly:
            command.extend(["-a", args.assembly])
        subprocess.call(command)
    else:
        print("Invalid workflow. Choose one of the following: 'genome_wf','metafast_wf', 'meta_wf'.")

def main(args):
    call_workflow(args)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='PGPg_finder is a tool that offers different workflows for genome and metagenome annotation. Select the workflow based on the tools and steps needed.', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-w', '--workflow', required=True, 
                        help='''Choose the workflow:\n
                        genome_wf - Genome annotation with diamond against PGPT-db database of PLaBAse.\n
                        metafast_wf - Fast metagenome direct read annotation with diamond against mgPGPT-db database of PLaBAse.\n
                        meta_wf - Metagenome assembly annotation (diamond) and coverage (Bowtie2, samtools bedtools) with mgPGPT-db database of PLaBAse.''')
                        
    parser.add_argument('-i', '--input', required=True, help='Path to the input directory.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output directory.')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use. Default is 1.')
    parser.add_argument('-a', '--assembly', help='Path to the assembly file. Use if you prefer to provide your own assembly.')
    
    args = parser.parse_args()
    
    main(args)

