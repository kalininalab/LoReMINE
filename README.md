# LoReMINE: Long Read-based Microbial genome mining pipeline

LoreMiNE (Long Read-based Microbial genome mining pipeline) is an end-to-end pipeline for microbial natural product discovery directly from long-read sequencing data. It integrates multiple sub-modules into a unified pipeline: (i) de-novo genome assembly using multiple long-read assemblers with automated selection of the highest-quality assembly, (ii) taxonomic classification, (iii) Biosynthetic gene clusters (BGC) detection, and (iv) clustering of BGCs to idenitfy gene cluster families (GCFs) which facilitates comparative analysis across different databases and taxa. By combining these steps into a unified pipeline, LoreMiNE facilitates comprehensive exploration of natural product diversity which has the potential to yield novel drug candidates.


# Installation

LoReMINE can be used by installing it via [conda](https://anaconda.org/kalininalab/loremine) 

To install it into the new environment:

`````shell
conda create -n loremine -c conda-forge -c kalininalab -c bioconda loremine
conda activate loremine
`````
Or you can also use mamba instead of conda (usually faster):

`````shell
mamba create -n loremine -c conda-forge -c kalininalab -c bioconda loremine
mamba activate loremine
`````

Or to install it into the already existing environment (for e.g: loremine):

`````shell
conda create -n loremine python=3.11.9
conda activate loremine
conda install -c conda-forge -c kalininalab -c bioconda loremine
`````
Or you can also use mamba instead of conda (usually faster):
`````shell
mamba create -n loremine python=3.11.9
mamba activate loremine
mamba install -c conda-forge -c kalininalab -c bioconda loremine
`````
We still need to install [BiG-SCAPE](https://github.com/medema-group/BiG-SCAPE) before using the LoReMINE pipeline. You can install it anywhere, but we prefer to install in the same conda environment directory (````` ~/anaconda3/envs/loremine/ `````) so that all tools are in one place. To install [BiG-SCAPE](https://github.com/medema-group/BiG-SCAPE), follow the instructions given below

`````shell
cd ~/anaconda3/envs/loremine/
git clone https://github.com/medema-group/BiG-SCAPE
cd BiG-SCAPE
pip install .
`````

# Setting up local databases

We still need to setup the local databases before using the LoReMINE pipeline. For DFAST_QC, we will download the compact version (````` DQC_REFERENCE_COMPACT `````, <1.5GB) of the database as the full version (````` DQC_REFERENCE_FULL `````) is quite huge (>100 GB). If you want to still download the full version, refer [here](https://github.com/nigyta/dfast_qc#quick-set-up-recommended). To download the compact version of the database, run the below command:
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

LoReMINE: Long Read-based Microbial genome mining pipeline

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

- ````` all_submodules: Run all the submodules (assemble, taxonomy, identify_bgcs, bgc_clustering) together in one run ````` this sub-module is used to run all the above 4 submodules together in a single run

# Usage of individual sub-modules
## assemble
This sub-module is used to perform the genome assembly using the long-reads from both Pacbio (raw as well as HiFi) and Nanopore. To see all the available options, run the following command
`````shell
loremine assemble --help
`````
this will produce the following output
`````shell

usage: loremine assemble [-h] [--reads READS] [--reads_type READS_TYPE] [--pacbio-raw PACBIO_RAW] [--pacbio-hifi PACBIO_HIFI] [--batch_run BATCH_RUN] [-g GENOME_SIZE] -o OUTPUT [-t THREADS] --prefix PREFIX [--alt_param ALT_PARAM]

options:
  -h, --help            show this help message and exit
  --reads READS         path to the input reads (.fastq or .fastq.gz format). If ".bam" file is available instead of ".fastq" file, then use the "--pacbio-raw" or "--pacbio-hifi"
  --reads_type READS_TYPE
                        type of reads in the ".fastq" or ".fastq.gz" file. Possible inputs are "raw_pacbio", "raw_nanopore" or "hifi_pacbio"
  --pacbio-raw PACBIO_RAW
                        path to the input Pacbio raw reads (.bam file)
  --pacbio-hifi PACBIO_HIFI
                        path to the input Pacbio HiFi reads (.bam file)
  --batch_run BATCH_RUN
                        path to the .tsv (tab seperated) file which contains 4 columns in the following order (Location of raw reads (.fastq format), type of reads (raw_pacbio, raw_nanopore or hifi_pacbio), genome size (default = 5000000 bp), prefix). No header is assumed, so
                        start from first line itself
  -g GENOME_SIZE, --genome-size GENOME_SIZE
                        estimated genome size (default = 5000000 bp (5Mbp))
  --weights WEIGHTS     weights of parameters (a,b & c) for calculating the assembly score while selecting the best assembly among different candidate assemblies. The formula to calculate assembly score is: 10*a + 2*b - 2*c + d/1e-6, where a = has_circular_chromosome, b =
                        no_of_circular_contigs, c = no_of_contigs & d = n50 (default = 10,2,2 as seen in formula)
  -o OUTPUT, --output OUTPUT
                        path to the save the output of the assembly
  -t THREADS, --threads THREADS
                        number of threads to use, default = 1
  --prefix PREFIX       Prefix for the output. If you use "batch_run" parameter, then provide "NA" as an input for this parameter
  --asm-coverage ASM_COVERAGE
                        reduced coverage for initial disjointig assembly (Used only for raw Pacbio and ONT reads) incase there is a high coverage of reads. Default value is not set, so that it uses all reads to perform the assembly. Incase, the initial disjointigs doesn't get
                        assembled due to very high coverage, then suggested value is "50", so that it uses longest 50x reads for initial disjointigs assembly
  --alt_param ALT_PARAM
                        Run the assembly using pacbio/nanopore raw reads with alternate parameters. Possible inputs are "True" or "False" (default = False). Use this parameter only when the assembly using default parameters in not satisfactory. Can only be used with
                        Pacbio/Nanopore "raw" reads and not with Pacbio "hifi" reads
  --force               Override the output in the existing output directory
  --verbose VERBOSE     Verbosity level of the output. Possible inputs are "0" or "1", where 0 = only prints status of the pipeline, 1 = prints status of the pipeline + output of each tools (default = 0)
`````

- ````` --reads `````: Use this option if you have `````.fastq````` or `````.fastq.gz````` file for the input reads. If you still have the `````.bam````` file for the reads from Pacbio, then use the `````--pacbio-raw````` or `````--pacbio-hifi````` option to input the reads file depending on whether the reads are raw or HiFi
  
- ````` --reads_type `````: Use this option to input the type of reads you have in `````.fastq````` or `````.fastq.gz````` file. Possible options for this parameter are `````raw_pacbio`````, `````raw_nanopore````` or `````hifi_pacbio`````. This parameter is important because depending on the type of input reads, the pipeline uses different assemblers to perform the assembly

- ````` --pacbio-raw `````: Use this option to input the `````.bam````` file, if your input reads are Pacbio raw reads

- ````` --pacbio-hifi `````: Use this option to input the `````.bam````` file, if your input reads are Pacbio HiFi reads

-  ````` --batch_run `````: Use this option to assemble multiple strains simultaneously instead of performing individual assemblies for each strain. An example file named `````input_batch_run.tsv````` is provided in the repository and can be used as a template for this parameter. As an input, it requires a **".tsv"** file which should contain 4 columns
    - 1st column is the path to raw reads (`````.fastq````` or `````.fastq.gz`````) file. It does not accept reads in `````.bam````` format. To convert the reads from `````.bam````` to `````.fastq````` format, you can use the following command `````samtools fastq input_filename.bam > output_filename.fastq`````
    - 2nd column should contain the type of reads (`````raw_pacbio`````, `````raw_nanopore````` or `````hifi_pacbio`````)
    - 3rd column should contain the genome size (if you already know, else the default is 5000000 bp)
    - 4th column should specify the **"prefix"** to use while saving the strain's assembly
    
- ````` -g `````: Use this option to input the estimated genome size of the strain whose assembly you are performing. It is used to estimate coverage and guide assembly algorithms while performing the assembly. For e.g: if your organism's genome size is `````10 Mbp (10 mega base pairs)`````, it should be given as `````10000000`````. Default value is `````5 Mbp (5000000 bp)`````

- ````` --weights `````: Use this option to input the weights of parameters (a,b & c) for calculating the assembly score while selecting the best assembly among different candidate assemblies. To calculate the assembly score, we use multiple parameters like whether the assembly has a circular chromosome **(highest priority)**, number of circular contigs **(bonus for other circular contigs)**, number of contigs **(strong penalty for fragmented assemblies)** & N50 **(small bonus for continuity)**. The formula to calculate assembly score is `````Assembly score = 10*a + 2*b - 2*c + d/1e-6`````, where a = has_circular_chromosome, b = no_of_circular_contigs, c = no_of_contigs & d = n50. Default value is `````10,2,2````` as we can observe in the formula that these are the weights of a, b & c
  
- ````` -o `````: Path to the output directory where you want to save the assembly output. After assembly completion, the assembly statistics and contig circularity information (whether a contig is circular or linear) are provided in the `````report.tsv````` and `````circularity.tsv````` files within the `````given_output_folder_path/assembly/given_prefix/quast_outputs/````` folder, while assembly completeness and contamination are reported in `````checkm_summary.tsv````` file inside the `````given_output_folder_path/assembly/given_prefix/checkm_output/````` folder. If you have Pacbio HiFi reads or you used alternate parameters (`````--alt_param`````) for performing the assembly using Pacbio & Nanopore raw reads, the pipeline generates multiple assemblies using different assemblers or parameter settings and automatically selects the best assembly. The best selected assembly, along with its completeness, contamination, and summary statistics, is provided in the `````chosen_best_assembly.txt````` file, while the best selected assembly is saved as `````given_prefix.fasta`````. Both files are located inside the `````given_output_folder_path/assembly/given_prefix/````` folder 

- ````` -t `````: Number of threads to use while running the assembly. Default value is `````1`````

- ````` --prefix `````: Specifies the prefix to be used for naming the assembly output files. **This parameter is mandatory**. If you are performing batch assemblies instead of a single-strain assembly, use **"NA"** as the value for this parameter

- ````` --asm-coverage `````: Use this parameter to limit read coverage during the initial disjointig assembly. In cases of extremely high read coverage, the initial disjointig assembly may fail.This parameter should be used only with **raw PacBio and Oxford Nanopore (ONT) reads** and must not be used with **HiFi reads**. By default, the parameter is not set, meaning all reads are used for assembly. If the initial disjointigs fail to assemble due to excessive coverage, a recommended value is `````50`````, which restricts the initial disjointigs assembly to the longest 50× reads

- ````` --alt_param `````: Use this parameter only if the assembly results obtained with the default settings are unsatisfactory. It is applicable only for PacBio or Nanopore **raw reads** and should not be used with PacBio **HiFi reads**. Possible values are `````True````` or `````False`````. Default value is `````False`````

- ````` --force `````: Use this parameter to allow writing output to an existing output directory. By default, the pipeline will terminate with an error if the specified output directory already exists. When `````--force````` is enabled, the existing output directory is overwritten and the pipeline proceeds normally

- ````` --verbose `````: Use this parameter to control the level of output printed on the console. Possible values are `````0````` or `````1`````, where `````0````` prints only high-level pipeline status messages (e.g., current step and tool being executed), while `````1````` prints pipeline status messages and the full stdout/stderr output from all tools. Default value is `````0`````
  
## taxonomy
This sub-module is responsible for determining the taxonomy of the input genome. To see all the available options, run the following command
`````shell
loremine taxonomy --help
`````
this will produce the following output
`````shell

usage: loremine taxonomy [-h] [-i INPUT_FASTA] [--input_dir INPUT_DIR] -o OUTPUT [-t THREADS]

options:
  -h, --help            show this help message and exit
  -i INPUT_FASTA, --input_fasta INPUT_FASTA
                        path to the input fasta file (Use this when you want to identify the taxonomy for single genome)
  --input_dir INPUT_DIR
                        path to the input directory containing multiple fasta files (Use this option to identify the taxonomy for multiple genomes)
  -o OUTPUT, --output OUTPUT
                        path to the save the output of the taxonomy
  -t THREADS, --threads THREADS
                        number of threads to use, default = 1
  --force               Override the output in the existing output directory
  --verbose VERBOSE     Verbosity level of the output. Possible inputs are "0" or "1", where 0 = only prints status of the pipeline, 1 = prints status of the pipeline + output of each tools (default = 0)
`````

- ````` -i `````: Path to the input FASTA file of the genome for which you want to identify the taxonomy. Use this option to determine the taxonomy of a single genome
  
- ````` --input_dir `````: Path to the input directory containing FASTA files of the genomes whose taxonomy you want to identify. Use this option to determine the taxonomy of multiple genomes in a single run

- ````` -o `````: Path to the output directory where you want to save the taxonomy output. Suppose, if your input FASTA file is named `````strainame.fasta`````, then the taxonomy results will be available in the `````given_output_folder_path/taxonomy/strainame/````` directory. You can find a summary of the taxonomy identified using both the NCBI and GTDB databases in the `````identified_taxonomy.txt````` file within the output directory (described in previous sentence), rather than checking individual output files

- ````` -t `````: Number of threads to use while identifying the taxonomy. Default value is `````1`````

- ````` --force `````: Use this parameter to allow writing output to an existing output directory. By default, the pipeline will terminate with an error if the specified output directory already exists. When `````--force````` is enabled, the existing output directory is overwritten and the pipeline proceeds normally

- ````` --verbose `````: Use this parameter to control the level of output printed on the console. Possible values are `````0````` or `````1`````, where `````0````` prints only high-level pipeline status messages (e.g., current step and tool being executed), while `````1````` prints pipeline status messages and the full stdout/stderr output from all tools. Default value is `````0`````

## identify_bgcs
This sub-module is responsible for identifying the BGCs in the input genome. To see all the available options, run the following command
`````shell
loremine identify_bgcs --help
`````
this will produce the following output
`````shell

usage: loremine identify_bgcs [-h] [-i INPUT_FASTA] [--input_dir INPUT_DIR] [--db_path DB_PATH] -o OUTPUT [-t THREADS]

options:
  -h, --help            show this help message and exit
  -i INPUT_FASTA, --input_fasta INPUT_FASTA
                        path to the input fasta file (Use this when you want to identify the BGCs for single genome)
  --input_dir INPUT_DIR
                        path to the input directory containing multiple fasta files (Use this option to identify the BGCs for multiple genomes)
  --db_path DB_PATH     path to the directory where you downloaded antismash databases (should point to directory which includes clusterblast, knownclusterblast, pfam etc as sub-directories). Use this option only when you downloaded databases at a custom location
  -o OUTPUT, --output OUTPUT
                        path to the output directory where you want to save the identified BGCs
  -t THREADS, --threads THREADS
                        number of threads to use, default = 1
  --force               Override the output in the existing output directory
  --verbose VERBOSE     Verbosity level of the output. Possible inputs are "0" or "1", where 0 = only prints status of the pipeline, 1 = prints status of the pipeline + output of each tools (default = 0)
`````

- ````` -i `````: Path to the input FASTA file of the genome for which you want to identify the BGCs. Use this option to idenitfy the BGCs of a single genome
  
- ````` --input_dir `````: Path to the input directory containing FASTA files of the genomes for which you want to identify BGCs. Use this option to detect BGCs across multiple genomes in a single run

- ````` --db_path `````: Path to the directory where the antiSMASH databases are downloaded (provide the full path to the directory containing subdirectories such as clusterblast, knownclusterblast, pfam, etc...). **Use this option only if you downloaded the antiSMASH databases at a custom location**. If you did not specify a custom location during downloading the antiSMASH database, **you can omit this parameter**, as the databases will be automatically fetched from the default directory within your conda environment

- ````` -o `````: Path to the output directory where you want to save the output of identified BGCs. Suppose, if your input FASTA file is named `````strainame.fasta`````, then the identified BGCs will be available in the `````given_output_folder_path/identified_bgcs/strainame/````` directory

- ````` -t `````: Number of threads to use while identifying the BGCs. Default value is `````1`````

- ````` --force `````: Use this parameter to allow writing output to an existing output directory. By default, the pipeline will terminate with an error if the specified output directory already exists. When `````--force````` is enabled, the existing output directory is overwritten and the pipeline proceeds normally

- ````` --verbose `````: Use this parameter to control the level of output printed on the console. Possible values are `````0````` or `````1`````, where `````0````` prints only high-level pipeline status messages (e.g., current step and tool being executed), while `````1````` prints pipeline status messages and the full stdout/stderr output from all tools. Default value is `````0`````

## bgc_clustering
This sub-module is responsible for clustering the BGCs to identify gene cluster families (GCFs). To see all the available options, run the following command
`````shell
loremine bgc_clustering --help
`````
this will produce the following output
`````shell

usage: loremine bgc_clustering [-h] [--input_dir INPUT_DIR] [--mibig]
                               [-o OUTPUT] --clustering_type CLUSTERING_TYPE
                               [-t THREADS] --pfam_dir PFAM_DIR
                               [--bigslice_cutoff BIGSLICE_CUTOFF]
                               [--bigscape_cutoff BIGSCAPE_CUTOFF]

options:
  -h, --help            show this help message and exit
  --input_dir INPUT_DIR
                        path to the input directory which contains all the bgcs for clustering
  --mibig               Use this option when you want to include MiBiG BGCs for clustering
  -o OUTPUT, --output OUTPUT
                        path to the output directory which will contain the clustering output
  --clustering_type CLUSTERING_TYPE
                        tool to use for clustering BGCs into GCFs. Possible inputs are "bigslice", "bigscape" or "both" (default = both)
  -t THREADS, --threads THREADS
                        number of threads to use, default = 1
  --pfam_dir PFAM_DIR   path to the directory where you have extracted the Pfam database. Please provide the complete path to the "Pfam-A.hmm" file
  --bigslice_cutoff BIGSLICE_CUTOFF
                        BiG-SLiCE cutoff value (default = 0.4)
  --bigscape_cutoff BIGSCAPE_CUTOFF
                        BiG-SCAPE cutoff value (default = 0.5)
  --force               Override the output in the existing output directory
  --verbose VERBOSE     Verbosity level of the output. Possible inputs are "0" or "1", where 0 = only prints status of the pipeline, 1 = prints status of the pipeline + output of each tools (default = 0)
`````

- ````` --input_dir `````: Path to the input directory containing all BGC (.gbk) files to be clustered into Gene Cluster Families (GCFs)

- ````` --mibig `````: This parameter indicates whether to include MiBiG (v4.0) BGCs for clustering or not. Using this parameter will **include MiBiG BGCs** along with input BGCs for clustering

- ````` --clustering_type `````: This parameter specifies the tool(s) to use for clustering BGCs into Gene Cluster Families (GCFs). Possible inputs are `````bigscape`````, `````bigslice`````, or `````both`````. Default value is `````both`````

- ````` -o `````: Path to the output directory where you want to save the clustering output. Suppose, If you have used both bigscape and bigslice for clustering, then the clustering output can be found in `````given_output_folder_path/clustering/bigscape/output/````` and `````given_output_folder_path/clustering/bigslice/output/````` directories respectively. Please refer the file `````output_clusters.tsv````` in both the output directories for the final clustering output. If 2 or more bgcs have **same id** in **"GCF_id"** column for bigslice output and **"Family"** column in bigscape output, it means that they belong to same GCF

- ````` -t `````: Number of threads to use while clustering the BGCs. Default value is `````1`````

- ````` --pfam_dir `````: Complete path to the directory where you have extracted the Pfam database as it will be used for BiG-SCAPE clustering. For e.g: I have extracted my Pfam database in a loremine environment directory, so I provide this `````~/anaconda3/envs/loremine/BiG-SCAPE/Pfam_dir/Pfam-A.hmm````` as as input

- ````` --bigslice_cutoff `````: Cut-off to use for clustering using BiG-SliCE. Default value is `````0.4`````. **Increasing** the cut-off will result in smaller & more clusters, while **decreasing** the cut-off will result in larger & fewer clusters

- ````` --bigscape_cutoff `````: Cut-off to use for clustering using BiG-SCAPE. Default value is `````0.5`````. **Decreasing** the cut-off will result in smaller & more clusters, while **increasing** the cut-off will result in larger & fewer clusters

- ````` --force `````: Use this parameter to allow writing output to an existing output directory. By default, the pipeline will terminate with an error if the specified output directory already exists. When `````--force````` is enabled, the existing output directory is overwritten and the pipeline proceeds normally

- ````` --verbose `````: Use this parameter to control the level of output printed on the console. Possible values are `````0````` or `````1`````, where `````0````` prints only high-level pipeline status messages (e.g., current step and tool being executed), while `````1````` prints pipeline status messages and the full stdout/stderr output from all tools. Default value is `````0`````

## all_submodules
This sub-module executes the complete pipeline in a single run. It takes raw sequencing reads (from one or multiple strains) as input, performs genome assembly, taxonomic identification, BGC detection, and BGC clustering to identify Gene Cluster Families (GCFs).To see all the available options, run the following command
`````shell
loremine all_submodules --help
`````
this will produce the following output
`````shell

usage: loremine all_submodules [-h] [--reads READS] [--reads_type READS_TYPE] [--pacbio-raw PACBIO_RAW] [--pacbio-hifi PACBIO_HIFI] [--batch_run BATCH_RUN] [-g GENOME_SIZE] -o OUTPUT [-t THREADS] --prefix PREFIX [--alt_param ALT_PARAM] [--db_path DB_PATH] [--mibig]
                               --clustering_type CLUSTERING_TYPE --pfam_dir PFAM_DIR [--bigslice_cutoff BIGSLICE_CUTOFF] [--bigscape_cutoff BIGSCAPE_CUTOFF]

options:
  -h, --help            show this help message and exit
  --reads READS         path to the input reads (.fastq or .fastq.gz format). If ".bam" file is available instead of ".fastq" file, then use the "--pacbio-raw" or "--pacbio-hifi"
  --reads_type READS_TYPE
                        type of reads in the ".fastq" or ".fastq.gz" file. Possible inputs are "raw_pacbio", "raw_nanopore" or "hifi_pacbio"
  --pacbio-raw PACBIO_RAW
                        path to the input Pacbio raw reads (.bam file)
  --pacbio-hifi PACBIO_HIFI
                        path to the input Pacbio HiFi reads (.bam file)
  --batch_run BATCH_RUN
                        path to the .tsv (tab seperated) file which contains 4 columns in the following order (Location of raw reads (.fastq format), type of reads (raw_pacbio, raw_nanopore or hifi_pacbio), genome size (default = 5000000 bp), prefix). No header is assumed, so
                        start from first line itself
  -g GENOME_SIZE, --genome-size GENOME_SIZE
                        estimated genome size (default = 5000000 bp (5Mbp))
  --weights WEIGHTS     weights of parameters (a,b & c) for calculating the assembly score while selecting the best assembly among different candidate assemblies. The formula to calculate assembly score is: 10*a + 2*b - 2*c + d/1e-6, where a = has_circular_chromosome, b =
                        no_of_circular_contigs, c = no_of_contigs & d = n50 (default = 10,2,2 as seen in formula)
  -o OUTPUT, --output OUTPUT
                        path to the save the output of the pipeline
  -t THREADS, --threads THREADS
                        number of threads to use, default = 1
  --prefix PREFIX       Prefix for the output. If you use "batch_run" parameter, then provide "NA" as an input for this parameter
  --asm-coverage ASM_COVERAGE
                        reduced coverage for initial disjointig assembly (Used only for raw Pacbio and ONT reads) incase there is a high coverage of reads. Default value is not set, so that it uses all reads to perform the assembly. Incase, the initial disjointigs doesn't get
                        assembled due to very high coverage, then suggested value is "50", so that it uses longest 50x reads for initial disjointigs assembly
  --alt_param ALT_PARAM
                        Run the assembly using pacbio/nanopore raw reads with alternate parameters. Possible inputs are "True" or "False" (default = False). Use this parameter only when the assembly using default parameters in not satisfactory. Can only be used with
                        Pacbio/Nanopore "raw" reads and not with Pacbio "hifi" reads
  --db_path DB_PATH     path to the directory where you downloaded antismash databases (should point to directory which includes clusterblast, knownclusterblast, pfam etc as sub-directories). Use this option only when you downloaded databases at a custom location
  --mibig               Use this option when you want to include MiBiG BGCs for clustering
  --clustering_type CLUSTERING_TYPE
                        tool to use for clustering BGCs into GCFs. Possible inputs are "bigslice", "bigscape" or "both" (default = both)
  --pfam_dir PFAM_DIR   Path to the directory where you have extracted the Pfam database. Please provide the complete path to the "Pfam-A.hmm" file
  --bigslice_cutoff BIGSLICE_CUTOFF
                        BiG-SLiCE cutoff value (default = 0.4)
  --bigscape_cutoff BIGSCAPE_CUTOFF
                        BiG-SCAPE cutoff value (default = 0.5)
  --force               Override the output in the existing output directory
  --verbose VERBOSE     Verbosity level of the output. Possible inputs are "0" or "1", where 0 = only prints status of the pipeline, 1 = prints status of the pipeline + output of each tools (default = 0)
`````

- ````` -o `````: Path to the output directory where you want to save the output of the pipeline. The output of all sub-modules will be saved in individual sub-module directories (`````given_output_folder_path/assembly/`````, `````given_output_folder_path/taxonomy/`````, `````given_output_folder_path/identified_bgcs/````` and `````given_output_folder_path/clustering/`````)

All other parameters used in this submodule have been thoroughly described in the respective sub-module sections above. Please refer to the corresponding documentation for detailed explanations.
