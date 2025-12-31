# Tutorial

Thank you very much for your interest in our pipeline, **PGPg_finder**, which was developed to address the need for identifying specific genes related to plant growth promotion in genomes and metagenomes. We are grateful to the PLaBAse team for creating a dedicated database for this purpose, which served as the foundation for the step-by-step pipeline we developed.

The first requirement is to have **Conda or Mamba** installed on your computer. If you do not have Conda installed, please follow the short instructions to install Miniforge or Conda for your operating system:
```bash
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m | sed 's/x86_64/x86_64/;s/aarch64/aarch64/').sh | bash

```

Once Conda or Mamba is available, you can clone the PGPg_finder repository using the command below:

```bash
git clone https://github.com/tpellegrinetti/PGPg_finder/
```

After cloning the repository, install PGPg_finder using the provided installation script:

```bash
bash install.sh
```

If the installation completes successfully, the pipeline is ready to use.

---

## Manual database installation (when PLaBAse is under maintenance)

Occasionally, PLaBAse may be temporarily unavailable due to maintenance. In such cases, the databases must be downloaded and installed manually by following the steps below.

First, download one or both databases, depending on your analysis.

PGPT_BASE_nr
[https://unitc-my.sharepoint.com/:u:/r/personal/iijag01_cloud_uni-tuebingen_de/Documents/PLaBAse/PGPT_BASE_nr_Aug2021n_ul_1.fasta.gz?csf=1&web=1&e=jhBIlL](https://unitc-my.sharepoint.com/:u:/r/personal/iijag01_cloud_uni-tuebingen_de/Documents/PLaBAse/PGPT_BASE_nr_Aug2021n_ul_1.fasta.gz?csf=1&web=1&e=jhBIlL)

mgPGPT-db
[https://unitc-my.sharepoint.com/:u:/r/personal/iijag01_cloud_uni-tuebingen_de/Documents/PLaBAse/mgPGPT-db_Feb2022_ul_dwnld.fasta.gz?csf=1&web=1&e=yZ2XCT](https://unitc-my.sharepoint.com/:u:/r/personal/iijag01_cloud_uni-tuebingen_de/Documents/PLaBAse/mgPGPT-db_Feb2022_ul_dwnld.fasta.gz?csf=1&web=1&e=yZ2XCT)

Place the downloaded files in the `database` directory of PGPg_finder. Then, with the PGPg_finder Conda environment activated, create the DIAMOND databases using the commands below:

```bash
diamond makedb --in PGPT_BASE_nr_Aug2021n_ul_1.fasta.gz --db genome
diamond makedb --in mgPGPT-db_Feb2022_ul_dwnld.fasta.gz --db metagenome
```

These commands will generate two `.dmnd` files. If desired, the original `.fasta.gz` files can be removed afterward.

After this step, PGPg_finder will function normally.

---

## Example usage of PGPg_finder

### Running PGPg_finder with genomes or MAGs

First, create a directory to store your genomes or MAGs and move into it:

```bash
mkdir PGPG
cd PGPG
```

Activate the PGPg_finder environment:

```bash
conda activate PGPg_finder
```

### Executing the genome workflow

The general command for running PGPg_finder with genomes is:

```bash
python PGPg_finder.py -w genome_wf -i input_directory -o output_directory -t threads
```

The `input_directory` should contain your genome files in FASTA format (`.fasta`, `.fa`, or `.fna`).

An example using the provided test genomes is shown below:

```bash
python PGPg_finder.py -w genome_wf -i genome_example/ -o genomeresult -t 22
```

---

## Additional parameters

PGPg_finder provides several optional parameters to customize DIAMOND searches:

```
--piden     Minimum identity percentage for DIAMOND (default: 30)
--qcov      Minimum query coverage percentage for DIAMOND (default: 30)
--extra     Additional DIAMOND arguments (optional)
--bitscore  Minimum bit score to report alignments
--evalue    Maximum e-value to report alignments (default: 1e-5)
--dmode     DIAMOND search mode (e.g., fast, sensitive, very-sensitive)
-h          Display the help message
```

You can adjust the identity threshold using `--piden`, the coverage threshold using `--qcov`, or modify the DIAMOND behavior by providing additional arguments with `--extra`. The alignment stringency can also be controlled using `--bitscore`, `--evalue`, and the DIAMOND search mode via `--dmode`.


## Analysis using raw reads

PGPg_finder also supports the identification of plant growth–promoting genes directly from metagenomic reads, without requiring pre-assembled genomes. This functionality is useful when working with large metagenomic datasets or when assemblies are not available or not desired.

Before starting, make sure the PGPg_finder Conda environment is activated:

```bash
conda activate PGPg_finder
```

### Read-based analysis without assembly (metafast_wf)

The `metafast_wf` workflow performs a read-based alignment against the PLaBAse databases using DIAMOND. In this mode, reads are directly queried against the reference database, allowing fast screening of plant growth–promoting genes across multiple metagenomes.

The general command is:

```bash
python PGPg_finder.py -w metafast_wf -i input_directory -o output_directory -t threads
```

The `input_directory` should contain your metagenomic read files. Paired-end reads must follow the `_1` and `_2` naming convention, while single-end reads are also supported. Results are written to the specified output directory.

Optional parameters such as minimum identity, query coverage, e-value, DIAMOND mode, and additional DIAMOND arguments can be adjusted using `--piden`, `--qcov`, `--evalue`, `--dmode`, and `--extra`, respectively. These options allow users to control the stringency and performance of the read-based search.

---

## Analysis using reads with assembly (meta_wf)

PGPg_finder also provides a workflow that combines metagenomic assembly with downstream gene detection and quantification. This approach is recommended when users want gene-level resolution, coverage estimates, and contig-based annotation.

In the `meta_wf` workflow, reads are first quality-trimmed, optionally assembled using MEGAHIT, and then processed through gene prediction, functional annotation, and read mapping to estimate gene abundance.

The general command is:

```bash
python PGPg_finder.py -w meta_wf -i reads_directory -o output_directory -t threads
```

By default, the workflow assembles the metagenomes internally. If precomputed assemblies are available, they can be provided using the `-a` option, allowing PGPg_finder to skip the assembly step and use existing contigs instead.

During execution, the pipeline performs quality trimming, assembly when required, gene prediction with Prodigal, functional annotation against the PLaBAse database using DIAMOND, and read mapping back to predicted genes to estimate coverage and relative abundance. The final output includes merged tables linking samples, genes, and abundance values, as well as summary visualizations such as heatmaps.

As with the read-based workflow, alignment parameters can be customized using `--piden`, `--qcov`, `--bitscore`, `--evalue`, and `--dmode` to adjust the sensitivity and specificity of the analysis.


