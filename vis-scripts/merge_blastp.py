import os
import argparse

def merge_blastp_annotations(blastout_file, output):
    print('Merge blastp annotations')
    with open(output, 'w') as fo:
        fo.write('#sample\tgene_id\taccession\n')
    if blastout_file.endswith('.txt'):
        basename = os.path.basename(blastout_file).split('_diamond')[0]
        with open(blastout_file, 'r') as fi:
            for line in fi:
                elements = line.split()
                gene_id = elements[0]
                accession = elements[1]
                with open(output, 'a') as fo:
                    fo.write(basename + '\t' + gene_id + '\t' + accession + '\n')

def main():
    parser = argparse.ArgumentParser(description="Merge blastp annotations into a single file.")
    parser.add_argument("-b", "--blastout_file", required=True, help="Path to the blastp output file.")
    parser.add_argument("-o", "--output", required=True, help="Output file to write merged annotations.")
    args = parser.parse_args()

    merge_blastp_annotations(args.blastout_file, args.output)

if __name__ == "__main__":
    main()

