# LoReMINE: Long Read-based Microbial genome mining pipeline

LoreMiNE (Long Read-based Microbial genome mining pipeline) is an end-to-end pipeline for microbial natural product discovery directly from long-read sequencing data. LoreMiNE integrates multiple modules into a unified pipeline: (i) de-novo genome assembly using multiple long-read assemblers with automated selection of the highest-quality assembly, (ii) taxonomic classification, (iii) Biosynthetic gene clusters (BGC) detection, and (iv) clustering of identified BGCs to facilitate comparative analysis across different databases and taxa. By combining these steps into a unified pipeline, LoreMiNE facilitates comprehensive exploration of natural product diversity which has the potential to yield novel drug candidates.


## Installation

LoReMINE can be used by installing it via [conda](https://anaconda.org/kalininalab/loremine) 

To install it into the new environment:

`````shell
conda create -n loremine -c conda-forge -c kalininalab -c bioconda loremine
conda activate loremine
`````

Or to install it into the already existing environment (for e.g: loremine):

`````shell
conda create -n loremine
conda activate
conda install -c conda-forge -c kalininalab -c bioconda loremine
`````
We still need to install [BiG-SCAPE](https://github.com/medema-group/BiG-SCAPE) before using the LoReMINE pipeline. You can install it anywhere, but we prefer to install in the same conda environment directory (````` ~/anaconda3/envs/loremine/ `````) so that all tools are in one place. To install [BiG-SCAPE](https://github.com/medema-group/BiG-SCAPE), follow the instructions given below

`````shell
cd ~/anaconda3/envs/loremine/
git clone https://github.com/medema-group/BiG-SCAPE
cd BiG-SCAPE
pip install .
`````

## Setting up local databases

We still need to setup the local databases for before using the LoReMINE pipeline. For DFAST_QC, we will download the compact version (````` DQC_REFERENCE_COMPACT `````, <1.5GB) of the database as the full version (````` DQC_REFERENCE_FULL `````) is quite huge (>100 GB). If you wnat to still download the full version, refer [here](https://github.com/nigyta/dfast_qc#quick-set-up-recommended). To download the compact version of the database, run the below command
`````shell
dqc_ref_manager.py download
`````
Once that is done, the [BiG-SLiCE](https://github.com/medema-group/bigslice) databases can be retrieved using the following command:
`````shell
download_bigslice_hmmdb
`````
Afterward, we need to set up the antiSMASH database, which will be used to identify BGCs. antiSMASH allows you to specify a custom location for the database (recommended option). Make sure to remember this location, since it will be needed later when running the pipeline. To set up the database in a custom directory, use:
`````shell
download-antismash-databases --database-dir custom_directory_to_download_antismash_database
`````
Finally, we need to set up the Pfam database, which will be used for clustering with [BiG-SCAPE](https://github.com/medema-group/BiG-SCAPE). You may install it in any custom directory, but be sure to remember this location as it will be required later to run the pipeline. To download and extract the Pfam database, run:
`````shell
cd custom_directory_to_download_Pfam_database
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
`````
All required databases are now set up, and the **LoReMINE** pipeline is ready to use.

## Usage

**LoReMINE** pipleine is now installed as a command-line tool. To get information about all the sub-modules, run:
`````shell
loremine --help
`````
this will produce the following output
`````shell

LoReMINE: Long read-based microbial genome mining pipeline

positional arguments:
  {assemble,taxonomy,identify_bgcs,bgc_clustering,all_submodules}
    assemble            run automated genome assembly pipeline
    taxonomy            identify the taxonomy of the genome
    identify_bgcs       identify the BGCs in the genome
    bgc_clustering      Cluster BGCs to identify Gene cluster families (GCFs)
    all_submodules      Run all the submodules (assemble, taxonomy, identify_bgcs, bgc_clustering) together in one run
`````

- ````` assemble: run automated genome assembly pipeline ````` this sub-module is used to perform the genome assembly. It supports long reads from both Pacbio as well as Nanopore. Depending on the type of reads (Pacbio CLR vs Pacbio HiFi vs Nanopore long reads), it uses different assemblers to perform the assembly and automatically selects the best assembly out of all the assemblies

- ````` taxonomy: identify the taxonomy of the genome ````` this sub-module is used to identify the taxonomy of the input genome. It uses both NCBI & GTDB databases to identify the taxonomy and gives a combined taxonomy output in a single text file

- ````` identify_bgcs: identify the BGCs in the genome ````` this sub-module is used to identify the BGCs in the input genome using antiSMASH

- ````` bgc_clustering: Cluster BGCs to identify Gene cluster families (GCFs) ````` this sub-module is used for clustering the input BGCs using both BiG-SLiCE & BiG-SCAPE to identify the GCFs. It also allows the user to include the [MiBiG](https://mibig.secondarymetabolites.org/) BGCs for comparision

- ````` all_submodules: Run all the submodules (assemble, taxonomy, identify_bgcs, bgc_clustering) together in one run ````` this sub-module is used to run all the above 4 submodules together 
