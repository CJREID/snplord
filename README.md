# snplord
Snplord is a Snakemake based pipeline that automates SNP-based phylogenetic analyses using snippy (https://github.com/tseemann/snippy), Gubbins (https://github.com/sanger-pathogens/gubbins), snp-sites (https://github.com/sanger-pathogens/snp-sites) and FastTree2.

This analysis is most suited to comparing closely related genomes for which there is a good quality reference genome available, however contig assemblies may also be used as a reference.

It requires Illumina short read pairs with filenames in format foo.R1.fastq.gz/foo.R2.fastq.gz and reference genome with filename foobar.fa (NOT .fasta).

## Installation
You will need to install Miniconda3 and create an environment 
```python
conda create -n snakemake -c bioconda snakemake
```
Clone the snplord repository into your home or data directory and enter this directory, this will be your working directory for analysis
```bash
git clone git@github.com:CJREID/snplord.git
cd snplord
```
Next you will need to create a directory called 'data' with two subdirectories called 'reads' and 'ref'
```
mkdir data
cd data
mkdir reads
mkdir ref
```
Copy your paired Illumina reads into 'reads' and your reference fasta file into 'ref'.

## Running snplord
Return to the snplord directory and activate your snakemake environment
```
source activate snakemake
```
In order to check the intended workflow that snplord will run, execute:
```
snakemake -np 
```
If you have no errors then run
```
snakemake -p --use-conda --resource mem_mb=64000 -s Snakefile_snplord
```

## Output
Output files will be present in four subdirectories of 'data'
1. snippyout

   This contains folders for each of your samples aligned to the reference genome,  
2. core

   This contains the outputs from snippy-core including full and core alignments. You may wish to use the 'core.aln' file to generate a    tree, however this workflow takes the full alignment 'full.core.aln' and runs it through recombination filtering with Gubbins before    identifying SNPs and inferring the tree.  
3. gubbins

   This folder contains all the output of the Gubbins analysis.  
4. fasttree

   This folder contains the final recombination filtered SNP tree,  
   
## Acknowledgments
I'd like to thank Max Cummins (github.com/maxlcummins) for his help with creating this pipeline. I'd also like to acknowledge Torsten Seemann (github.com/tseemann) and the sanger-pathogens group (github.com/sanger-pathogens) for developing programs used in this pipeline.

## Citations
Price MN, Dehal PS, Arkin AP. FastTree 2--approximately maximum-likelihood trees for large alignments. PLoS One 2010;5(3):e9490.   

