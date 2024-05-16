Wants to annotate ncRNA in a genome, but is having trouble navigating the dozens of different tools and databases out there? Trying to find functions for lncRNAs, but finding almost nothing? 

This software may help you solve those problems.

# RNA Gatherer

RNA Gatherer is a software with ready to use pipelines for:

- explorer.py: Annotation and prediction of ncRNA in genomes, taking into account transcriptome data, covariance models, reference sequences, reference annotations and data from public APIs;
- prophet.py: Computational prediction of lncRNA functions using gene coexpression;

## Installation

RNA Gatherer requires some databases and software in order to run. It was developed for Linux x64 environments and uses a command line interface.

First of all, you should clone (or [download](https://github.com/pentalpha/rna_gatherer/archive/master.zip)) this repository:

```sh
git clone https://github.com/pentalpha/rna_gatherer.git
cd rna_gatherer
```
### Databases

|File|Is it mandatory?|Download Link|
|:-|:-:|-:|
|Gene Ontology Graph | Yes | [http://purl.obolibrary.org/obo/go.obo](http://purl.obolibrary.org/obo/go.obo) |
|RFAM Covariance Models | Yes | [ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz](ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz) |
|Non-Redundant Proteins | Only if you want to remove known protein's mRNA from lncRNA data | [ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz](ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz) |
|ncRNA Database FASTAs | Only if you want to look for known ncRNA through alignment | It can be ANY .fasta file. We suggest using RNA Central's database: [ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_active.fasta.gz](ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_active.fasta.gz) |

After downloading them, edit the config.json file to include the full paths. If the file does not already exist, create it:

```sh
cp config.dummy.json config.json
```

Now, open config.json with your favorite text editor. Fill in the empty fields with the path to the downloaded files:

```json
[...]
    
    "rna_dbs": {'DB Name': db_path, 
        'DB Name 2': db_path_2, ...},
    "non_redundant": "path/to/nr.fasta",
    "go_obo": "path/to/go.obo",
    "rfam_cm": "path/to/Rfam.cm"
}
```

Non-mandatory fields can be left empty.

### Required software

The required software are listed in the [environment.yml](environment.yml) file. Using [conda](https://docs.conda.io/en/latest/miniconda.html), you can create the environment in one command:

```sh
conda env create -f environment.yml
```

Now activate the fresh new environment in order to use the software:

```sh
conda activate rna
```

## How To Use

### explorer.py

This is an extensive pipeline for detecting ncRNA in a given genome. Given a genome (and maybe some optional inputs), it will give you a non-redundant .GFF annotation file and a .TSV file with functional annotations, based on RFAM and other databases.

A basic command would be:

```sh
python explorer.py -g [genome.fasta] \
    -tx [taxonomic ID for species] \
    -o [output directory]
```

These are the only required input arguments. But other inputs can be passed in order to make the annotation a lot better!

This enables the annotation of lncRNA transcripts:
```sh
python explorer.py -g [genome.fasta]\
    -tx [taxonomic ID for species] \
    -tr [transcriptome.fasta] \
    -o [output directory]
```

This includes a ncRNA reference annotation file (.gff format):
```sh
python explorer.py -g [genome.fasta] \
    -tx [taxonomic ID for species] \
    -gff [reference.gff] \
    -o [output directory]
```
You can find reference files like these for many species [here](ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/genome_coordinates/gff3/). Please note that inclusing mRNA in the reference annotation can mess things up a little bit...

Many species have reference ncRNA sequences out there with no position in the genome. RNA Gatherer can map them for you:
```sh
python explorer.py -g [genome.fasta]\
    -tx [taxonomic ID for species] \
    -ref [reference.fasta] \
    -o [output directory]
```

For more detailed description of the command line arguments, use --help:
```sh
python explorer.py --help
```

### prophet.py

Given a count reads table, a list of lncRNA names and a annotation of coding genes, this tool enables you to predict the functions (Gene Ontology terms) of lncRNA.

An example command:
```sh
python prophet.py -cr test_data/counts/mus_musculus_tpm.tsv \
    -reg test_data/lnc_list/mus_musculus_lncRNA.txt \
    -ann test_data/annotation/mgi_genes_annotation.tsv \
    -o output_directory
```

-cr: The count reads table is a simple .TSV table where the first row is the sample names and the following rows start with a gene name, followed by the read counts at each sample. The counts must be normalized, preferably with TPM. It must include counts for both lncRNA and mRNA ([example](test_data/counts/mus_musculus_tpm.tsv)).

-reg: The lncRNA list specifies which ones of the genes in the count reads table are lncRNA. It's a simple .TXT file where every line is a lncRNA name ([example](test_data/lnc_list/mus_musculus_lncRNA.txt)).

-ann: The functional annotation for the coding genes, another .TSV table. Each line contains a gene name, a GO term and the respective ontology - molecular_function, biological_process or cellular_component ([example](test_data/annotation/mgi_genes_annotation.tsv)).

For more detailed description of the command line arguments, use --help:
```sh
python prophet.py --help
```

## Project structure



## To-do list

- Create an utility to download the databases for the user;
