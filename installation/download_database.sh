#!/bin/bash
cd ../database/
wget https://plabase.cs.uni-tuebingen.de/pb/tools/PGPTblhm/data/factors/PGPT_BASE_nr_Aug2021n_ul_1.fasta.gz
wget https://plabase.cs.uni-tuebingen.de/pb/tools/PGPTblhm/data/factors/mgPGPT/mgPGPT-db_Feb2022_ul_dwnld.fasta.gz
gzip -d *.gz
diamond makedb --in PGPT_BASE_nr_Aug2021n_ul_1.fasta --db genome
diamond makedb --in mgPGPT-db_Feb2022_ul_dwnld.fasta --db metagenome
rm PGPT_BASE_nr_Aug2021n_ul_1.fasta
rm mgPGPT-db_Feb2022_ul_dwnld.fasta
cd ../installattion