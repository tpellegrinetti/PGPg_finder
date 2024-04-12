# Tutorial
## 1. First steps
* Create a directory for PGPg_finder
```bash
mkdir pgpg
```
* Enter in directory maked
```bash
cd pgpg
```
### 1.1. Instalation
* Copy the archives for your computer
```bash
git clone https://github.com/tpellegrinetti/PGPg_finder.git
```
* Enter in directory which copy from github
```bash
cd PGPg_finder
```
* Install the PGPg_finder in your computer
```bash
bash install.sh
```
## 2. Executing the PGPg for genomes
* General comand
```bash
python PGPb_finder.py -w genome_wf -i input_directory -o output_directory -t threads
```

* Genomes example
```bash
python PGPg_finder.py -w genome_wf -i genome_example/ -o genomeresult -t 22
```
* **The "-a" argument is optional and works if you want to provide assembly files**
* **Is not necessary create output folder before run the workflow**
* 
## 3. Executing the PGPg for metagenomes
### 3.1 greater accuracy

* General comand
  ```bash
  python PGPb_finder.py -w meta_wf -i input_directory -o output_directory -t threads
  ```
