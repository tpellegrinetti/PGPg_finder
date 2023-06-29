import os
import argparse
import logging

def merge_abun_tab(abundance_file, blastp, output):
    logging.info('Merge abundance table with blastp table')

    # Check if the file exists. If not, writes the header.
    if not os.path.exists(output):
        with open(output, 'w') as fo:
            fo.write('Sample\tID\tCount\n')  

    abun_tab_dict = {}
    basename = os.path.basename(abundance_file).rsplit('.', 1)[0]
    with open(abundance_file) as fi:
        next(fi)
        for line in fi:
            line = line.strip('\n')
            gene_id, abundance = line.split('\t')
            if gene_id != "#ID":
                key = str(basename) + '+' + str(gene_id)
                abun_tab_dict[key] = abundance

    a = []
    with open(blastp, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                a.append(line.strip())

    for item in sorted(a):
        basename = item.split('\t')[0]
        gene_id = item.split('\t')[1]
        accession = item.split('\t')[2]
        key = basename + '+' + gene_id
        if key in abun_tab_dict:
            abundance = abun_tab_dict[key]
            with open(output, 'a') as fo:
                fo.write(basename + '\t' + accession + '\t' + abundance + '\n')

def main():
    parser = argparse.ArgumentParser(description='Merge abundance table with blastp table.')
    parser.add_argument('-a', '--abundance_file', help='Path to the abundance file', required=True)
    parser.add_argument('-b', '--blastp', help='Filepath of the merged blastp table', required=True)
    parser.add_argument('-o', '--output', help='Output file', required=True)
    args = parser.parse_args()
    merge_abun_tab(args.abundance_file, args.blastp, args.output)

if __name__ == '__main__':
    main()

