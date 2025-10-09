# LoReMINE: Long Read-based Microbial genome mining pipeline

LoreMiNE (Long Read-based Microbial genome mining pipeline) is an end-to-end pipeline for microbial natural product discovery directly from long-read sequencing data. LoreMiNE integrates multiple modules into a unified pipeline: (i) de-novo genome assembly using multiple long-read assemblers with automated selection of the highest-quality assembly, (ii) taxonomic classification, (iii) Biosynthetic gene clusters (BGC) detection, and (iv) clustering of identified BGCs to facilitate comparative analysis across different databases and taxa. By combining these steps into a unified pipeline, LoreMiNE facilitates comprehensive exploration of natural product diversity which has the potential to yield novel drug candidates.


# Installation

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

# Setting up local databases

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

# Usage

**LoReMINE** pipeline is now installed as a command-line tool. To get information about all the sub-modules, run:
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

# Usage of individual submodules
## Assemble
This submodule is used to perform the genome assembly using the long-reads from both Pacbio (raw as well as HiFi) and Nanopore. To see all the available options, run the following command
`````shell
loremine assemble --help
`````
this will produce the following output
`````shell

usage: loremine assemble [-h] [--reads READS] [--reads_type READS_TYPE] [--pacbio-raw PACBIO_RAW] [--pacbio-hifi PACBIO_HIFI] [--batch_run BATCH_RUN] [-g GENOME_SIZE] -o OUTPUT [-t THREADS] --prefix PREFIX [--alt_param ALT_PARAM]

options:
  -h, --help            show this help message and exit
  --reads READS         path to the input reads (.fastq format). If ".bam" file is available instead of ".fastq" file, then use the "--pacbio-raw" or "--pacbio-hifi"
  --reads_type READS_TYPE
                        type of reads in the ".fastq" file. Possible inputs are "raw_pacbio", "raw_nanopore" and "hifi_pacbio"
  --pacbio-raw PACBIO_RAW
                        path to the input Pacbio raw reads (.bam file)
  --pacbio-hifi PACBIO_HIFI
                        path to the input Pacbio HiFi reads (.bam file)
  --batch_run BATCH_RUN
                        path to the .tsv (tab seperated) file which contains 4 columns in the following order (Location of raw reads (.fastq format), type of reads (raw_pacbio, raw_nanopore or hifi_pacbio), genome size (deafult = 5000000 bp), prefix). No header is assumed, so
                        start from first line itself
  -g GENOME_SIZE, --genome-size GENOME_SIZE
                        estimated genome size (default = 5000000 bp (5Mbp))
  -o OUTPUT, --output OUTPUT
                        path to the save the output of the assembly
  -t THREADS, --threads THREADS
                        number of threads to use, default = 1
  --prefix PREFIX       Prefix for the output. If you use "batch_run" parameter, then write "NA" in this field
  --alt_param ALT_PARAM
                        Run the assembly using pacbio/nanopore raw reads with alternate parameters (True or False). Use this parameter only when the assembly using default parameters in not satisfactory. Can only be used with Pacbio/Nanopore "raw" reads and not with Pacbio "hifi"
                        reads
`````

- ````` --reads `````: Use this option if you have `````.fastq````` or `````.fastq.gz````` file for the input reads. If you still have the `````.bam````` file for the reads from Pacbio, then use the `````--pacbio-raw````` or `````--pacbio-hifi````` option to input the reads file depending on whether the reads are raw or HiFi
  
- ````` --reads_type `````: Use this option to input the type of reads you have in `````.fastq````` or `````.fastq.gz````` file. Possible option for this parameter are `````raw_pacbio`````, `````raw_nanopore````` or `````hifi_pacbio`````. This parameter is important because depending on the type of input reads, the pipeline uses different assemblers to perform the assembly. If you have `````hifi_pacbio````` reads, then the assembly is performed using 3 different assemblers and the best assembly is automatically selected by the pipeline. The best assembly selected by the pipeline can be found in `````chosen_best_assembly.txt````` file which can be found inside the `````given_output_folder_path/assembly/given_prefix/````` folder

- ````` pacbio-raw `````: Use this option to input the `````.bam````` file, if your input reads are Pacbio raw reads

- ````` pacbio-hifi `````: Use this option to input the `````.bam````` file, if your input reads are Pacbio HiFi reads

-  ````` batch_run `````: Use this option to assemble multiple strains simultaneously instead of performing individual assemblies for each strain. As an input, it requires a **".tsv"** file which should contain 4 columns
    - 1st column is the path to raw reads (`````.fastq````` or `````.fastq.gz`````) file. It does not accept reads in `````.bam````` format. To convert the reads from `````.bam````` to `````.fastq````` format, you can use the following command `````samtools fastq input_filename.bam > output_filename.fastq`````
    - 2nd column should contain the type of reads (`````raw_pacbio`````, `````raw_nanopore````` or `````hifi_pacbio`````)
    - 3rd column should contain the genome size (if you already know, else the deafult is 5000000 bp)
    - 4th column should specify the **"prefix"** to use while saving the strain's assembly
    
- ````` -g `````: Use this option to input the estimated genome size of the strain whose assembly you are performing. It is used to estimate coverage and guide assembly algorithms while performing the assembly. For e.g: if your organism's genome size is `````10 Mbp (10 mega base pairs)`````, it should be given as `````10000000`````. Default value is `````5 Mbp (5000000 bp)`````
