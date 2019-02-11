# snplord
Snplord is a Snakemake based pipeline that automates SNP-based phylogenetic analyses using snippy (https://github.com/tseemann/snippy), Gubbins (https://github.com/sanger-pathogens/gubbins), snp-sites (https://github.com/sanger-pathogens/snp-sites) and FastTree2.

This analysis is most suited to comparing closely related genomes for which there is a good quality reference genome available, however contig assemblies may also be used as a reference.

It requires Illumina short read pairs with filenames in format foo.R1.fastq.gz/foo.R2.fastq.gz and reference genome with filename foobar.fa (NOT .fasta).

## Installation
You will need to install Miniconda3 and create an environment 
```python
conda create -n snakemake -c bioconda snakemake
```
Create a working directory for your analysis called "snake-snippy" or similar and clone this repository into it
```bash
mkdir snake-snippy
cd snake-snippy
git clone 
```
