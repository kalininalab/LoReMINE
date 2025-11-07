"""This is the entry point of the LoReMINE pipeline"""

import argparse
import sys
from src.convert_reads_to_fastq import *
from src.run_pacbio_raw_assembly import *
from src.run_nanopore_raw_assembly import *
from src.run_pacbio_hifi_assembly import *
from src.evaluate_best_assembly import *
from src.identify_taxonomy import *
from src.identify_bgcs import *
from src.run_bgc_clustering import *
from src.run_all_submodules import *
import shutil
import os
import glob
from pathlib import Path

def main():
    # Create the parser
    parser = argparse.ArgumentParser(prog="loremine", description="LoReMINE: Long Read-based Microbial genome mining pipeline")

    #################################################  Parser for running the genome assembly ############################################

    subparsers = parser.add_subparsers(dest='command', required=True)
    parser_assemble = subparsers.add_parser('assemble', help='run automated genome assembly pipeline')
    parser_assemble.add_argument('--reads', type=str, help='path to the input reads (.fastq or .fastq.gz format). If \".bam\" file is available instead of \".fastq\" file, then use the \"--pacbio-raw\" or \"--pacbio-hifi\"')
    parser_assemble.add_argument('--reads_type',  type=str, help='type of reads in the \".fastq\" or \".fastq.gz\" file. Possible inputs are \"raw_pacbio\", \"raw_nanopore\" or \"hifi_pacbio\"')
    parser_assemble.add_argument('--pacbio-raw',type=str, help='path to the input Pacbio raw reads (.bam file)')
    parser_assemble.add_argument('--pacbio-hifi',type=str, help='path to the input Pacbio HiFi reads (.bam file)')
    parser_assemble.add_argument('--batch_run', type=str, help='path to the .tsv (tab seperated) file which contains 4 columns in the following order (Location of raw reads (.fastq format), type of reads (raw_pacbio, raw_nanopore or hifi_pacbio), genome size (default = 5000000 bp), prefix). No header is assumed, so start from first line itself')
    parser_assemble.add_argument('-g', '--genome-size', type=str, help='estimated genome size (default = 5000000 bp (5Mbp))',default="5000000")
    parser_assemble.add_argument('-o', '--output',type=str, help='path to the save the output of the assembly', required=True)
    parser_assemble.add_argument('-t', '--threads', type=str, help='number of threads to use, default = 1', default = "1")
    parser_assemble.add_argument('--prefix',type=str, help='Prefix for the output. If you use "batch_run" parameter, then provide "NA" as an input for this parameter', required=True)
    parser_assemble.add_argument('--alt_param', type=str, help='Run the assembly using pacbio/nanopore raw reads with alternate parameters. Possible inputs are \"True\" or \"False\" (default = False). Use this parameter only when the assembly using default parameters in not satisfactory. Can only be used with Pacbio/Nanopore "raw" reads and not with Pacbio "hifi" reads', default=False)

    #################################################  Parser for identifying the taxonomy of the genome ############################################

    parser_taxonomy = subparsers.add_parser('taxonomy', help='identify the taxonomy of the genome')
    parser_taxonomy.add_argument('-i', '--input_fasta', type=str, help='path to the input fasta file (Use this when you want to identify the taxonomy for single genome)')
    parser_taxonomy.add_argument('--input_dir', type=str, help='path to the input directory containing multiple fasta files (Use this option to identify the taxonomy for multiple genomes)')
    parser_taxonomy.add_argument('-o', '--output', type=str, help='path to the save the output of the taxonomy', required=True)
    parser_taxonomy.add_argument('-t', '--threads', type=str, help='number of threads to use, default = 1', default = "1")

    #################################################  Parser for identifying the BGCs in the genome ############################################

    parser_bgc_identification = subparsers.add_parser('identify_bgcs', help='identify the BGCs in the genome')
    parser_bgc_identification.add_argument('-i', '--input_fasta', type=str, help='path to the input fasta file (Use this when you want to identify the BGCs for single genome)')
    parser_bgc_identification.add_argument('--input_dir', type=str, help='path to the input directory containing multiple fasta files (Use this option to identify the BGCs for multiple genomes)')
    parser_bgc_identification.add_argument('--db_path', type=str, help='path to the directory where you downloaded antismash databases (should point to directory which includes clusterblast, knownclusterblast, pfam etc as sub-directories). Use this option only when you downloaded databases at a custom location')
    parser_bgc_identification.add_argument('-o', '--output', type=str, help='path to the output directory where you want to save the identified BGCs', required=True)
    parser_bgc_identification.add_argument('-t', '--threads', type=str, help='number of threads to use, default = 1', default = "1")

    #################################################  Parser for running the BGC clustering ############################################

    parser_bgc_clustering = subparsers.add_parser('bgc_clustering', help='Cluster BGCs to identify Gene cluster families (GCFs)')
    parser_bgc_clustering.add_argument('--input_dir', type=str, help='path to the input directory which contains all the bgcs for clustering')
    parser_bgc_clustering.add_argument('--mibig', action='store_true', help='Use this option when you want to include MiBiG BGCs for clustering')
    parser_bgc_clustering.add_argument('-o', '--output', type=str, help='path to the output directory which will contain the clustering output')
    parser_bgc_clustering.add_argument('--clustering_type', type=str, help='tool to use for clustering BGCs into GCFs. Possible inputs are \"bigslice\", \"bigscape\" or \"both\" (default = both)', required=True, default = "both")
    parser_bgc_clustering.add_argument('-t', '--threads', type=str, help='number of threads to use, default = 1', default = "1")
    parser_bgc_clustering.add_argument('--pfam_dir', type=str, help='path to the directory where you have extracted the Pfam database. Please provide the complete path to the "Pfam-A.hmm" file', required=True)
    parser_bgc_clustering.add_argument('--bigslice_cutoff',type=float, help='BiG-SLiCE cutoff value (default = 0.4)', default = 0.4)
    parser_bgc_clustering.add_argument('--bigscape_cutoff', type=float, help='BiG-SCAPE cutoff value (default = 0.5)', default = 0.5)

    #################################################  Parser for running the whole pipeline together ############################################

    parser_all_submodules = subparsers.add_parser('all_submodules', help='Run all the submodules (assemble, taxonomy, identify_bgcs, bgc_clustering) together in one run')
    parser_all_submodules.add_argument('--reads', type=str, help='path to the input reads (.fastq or .fastq.gz format). If \".bam\" file is available instead of \".fastq\" file, then use the \"--pacbio-raw\" or \"--pacbio-hifi\"')
    parser_all_submodules.add_argument('--reads_type', type=str, help='type of reads in the \".fastq\" or \".fastq.gz\" file. Possible inputs are \"raw_pacbio\", \"raw_nanopore\" or \"hifi_pacbio\"')
    parser_all_submodules.add_argument('--pacbio-raw', type=str, help='path to the input Pacbio raw reads (.bam file)')
    parser_all_submodules.add_argument('--pacbio-hifi', type=str, help='path to the input Pacbio HiFi reads (.bam file)')
    parser_all_submodules.add_argument('--batch_run', type=str, help='path to the .tsv (tab seperated) file which contains 4 columns in the following order (Location of raw reads (.fastq format), type of reads (raw_pacbio, raw_nanopore or hifi_pacbio), genome size (default = 5000000 bp), prefix). No header is assumed, so start from first line itself')
    parser_all_submodules.add_argument('-g', '--genome-size', type=str, help='estimated genome size (default = 5000000 bp (5Mbp))', default="5000000")
    parser_all_submodules.add_argument('-o', '--output', type=str, help='path to the save the output of the pipeline', required=True)
    parser_all_submodules.add_argument('-t', '--threads', type=str, help='number of threads to use, default = 1', default = "1")
    parser_all_submodules.add_argument('--prefix', type=str, help='Prefix for the output. If you use "batch_run" parameter, then provide "NA" as an input for this parameter', required=True)
    parser_all_submodules.add_argument('--alt_param', type=str, help='Run the assembly using pacbio/nanopore raw reads with alternate parameters. Possible inputs are \"True\" or \"False\" (default = False). Use this parameter only when the assembly using default parameters in not satisfactory. Can only be used with Pacbio/Nanopore "raw" reads and not with Pacbio "hifi" reads', default=False)
    parser_all_submodules.add_argument('--db_path', type=str, help='path to the directory where you downloaded antismash databases (should point to directory which includes clusterblast, knownclusterblast, pfam etc as sub-directories). Use this option only when you downloaded databases at a custom location')
    parser_all_submodules.add_argument('--mibig', action='store_true', help='Use this option when you want to include MiBiG BGCs for clustering')
    parser_all_submodules.add_argument('--clustering_type', type=str, help='tool to use for clustering BGCs into GCFs. Possible inputs are \"bigslice\", \"bigscape\" or \"both\" (default = both)', required=True, default = "both")
    parser_all_submodules.add_argument('--pfam_dir', type=str, help='Path to the directory where you have extracted the Pfam database. Please provide the complete path to the "Pfam-A.hmm" file', required=True)
    parser_all_submodules.add_argument('--bigslice_cutoff', type=float, help='BiG-SLiCE cutoff value (default = 0.4)', default = 0.4)
    parser_all_submodules.add_argument('--bigscape_cutoff', type=float, help='BiG-SCAPE cutoff value (default = 0.5)', default = 0.5)

    args = parser.parse_args()

    if args.command == 'assemble':

        if args.batch_run == None:

            """Performing the assembly for single strain"""

            if args.pacbio_raw != None:  # If the input reads are in .bam file, then converting it to fastq
                print('Pacbio raw reads provided at this location:' + os.path.abspath(args.pacbio_raw))
                convert_reads(args)
                print("\n\nStaring the assembly using pacbio raw reads now.....\n\n")
                perform_assembly_raw_reads_pacbio(args)

            elif args.pacbio_hifi != None:  # If the input reads are in .bam file, then converting it to fastq
                print('Pacbio HiFi reads provided at this location:' + os.path.abspath(args.pacbio_hifi))
                convert_reads(args)
                print("\n\nStaring the assembly using pacbio HiFi reads now.....\n\n")
                perform_assembly_hifi_reads(args)

            elif args.pacbio_raw != None and args.pacbio_hifi != None:
                print("Sorry, we need only single type of reads (Raw or HiFi)")
                sys.exit(1)

            elif args.reads != None and args.reads_type == "raw_pacbio":
                print('Pacbio raw reads provided at this location:' + os.path.abspath(args.reads))
                os.makedirs(args.output + "/raw_reads")
                suffixes = Path(args.reads).suffixes
                raw_reads_filename = args.prefix + "".join(suffixes)
                shutil.copyfile(args.reads, args.output + "/raw_reads/" + raw_reads_filename)
                print("\n\nStaring the assembly using pacbio raw reads now.....\n\n")
                perform_assembly_raw_reads_pacbio(args)

            elif args.reads != None and args.reads_type == "hifi_pacbio":
                print('Pacbio HiFi reads provided at this location:' + os.path.abspath(args.reads))
                os.makedirs(args.output + "/raw_reads")
                suffixes = Path(args.reads).suffixes
                raw_reads_filename = args.prefix + "".join(suffixes)
                shutil.copyfile(args.reads, args.output + "/raw_reads/" + raw_reads_filename)
                print("\n\nStaring the assembly using pacbio HiFi reads now.....\n\n")
                perform_assembly_hifi_reads(args)

            elif args.reads != None and args.reads_type == "raw_nanopore":
                print('Nanopore raw reads provided at this location:' + os.path.abspath(args.reads))
                os.makedirs(args.output + "/raw_reads")
                suffixes = Path(args.reads).suffixes
                raw_reads_filename = args.prefix + "".join(suffixes)
                shutil.copyfile(args.reads, args.output + "/raw_reads/" + raw_reads_filename)
                print("\n\nStaring the assembly using nanopore raw reads now.....\n\n")
                perform_assembly_raw_reads_nanopore(args)

        else:

            """Performing the assembly for multiple strains"""

            output_path = os.path.abspath(args.output)
            with open(args.batch_run) as batch_file:
                for line in batch_file.readlines():
                    split_array = line.split('\t')
                    raw_reads_location = split_array[0].strip()
                    reads_type = split_array[1].strip()
                    if split_array[2].strip() == "":
                        genome_size = str(5000000)
                    else:
                        genome_size = split_array[2].strip()
                    prefix = split_array[3].strip()

                    if raw_reads_location.endswith('.bam'):
                        print('\n\nSorry, we only accept reads in ".fastq" format for batch run. You can convert your reads into ".fastq" format using "samtools fastq your_input_bam_file > output_prefix.fastq". Then you can provide the paths to these ".fastq" files for running the pipeline\n\n')
                        sys.exit(1)

                    if reads_type == "hifi_pacbio":
                        print('Pacbio HiFi reads provided at this location:' + raw_reads_location)
                        if not os.path.exists(output_path + "/raw_reads"):
                            os.makedirs(output_path + "/raw_reads")
                        suffixes = Path(raw_reads_location).suffixes
                        raw_reads_filename = prefix + "".join(suffixes)
                        shutil.copyfile(raw_reads_location, output_path + "/raw_reads/" + raw_reads_filename)
                        print("\n\nStaring the assembly using pacbio HiFi reads now.....\n\n")
                        perform_assembly_hifi_reads_batch_run(args, output_path + "/raw_reads/" + prefix + ".fastq", genome_size, prefix, output_path)

                    elif reads_type == "raw_pacbio":
                        print('Pacbio raw reads provided at this location:' + raw_reads_location)
                        if not os.path.exists(output_path + "/raw_reads"):
                            os.makedirs(output_path + "/raw_reads")
                        suffixes = Path(raw_reads_location).suffixes
                        raw_reads_filename = prefix + "".join(suffixes)
                        shutil.copyfile(raw_reads_location, output_path + "/raw_reads/" + raw_reads_filename)
                        print("\n\nStaring the assembly using pacbio raw reads now.....\n\n")
                        perform_assembly_raw_reads_pacbio_batch_run(args, output_path + "/raw_reads/" + raw_reads_filename, genome_size, prefix, output_path)

                    elif reads_type == "raw_nanopore":
                        print('Nanopore raw reads provided at this location:' + raw_reads_location)
                        if not os.path.exists(output_path + "/raw_reads"):
                            os.makedirs(output_path + "/raw_reads")
                        suffixes = Path(raw_reads_location).suffixes
                        raw_reads_filename = prefix + "".join(suffixes)
                        shutil.copyfile(raw_reads_location, output_path + "/raw_reads/" + raw_reads_filename)
                        print("\n\nStaring the assembly using nanopore raw reads now.....\n\n")
                        perform_assembly_raw_reads_nanopore_batch_run(args, output_path + "/raw_reads/" + raw_reads_filename, genome_size, prefix, output_path)


    elif args.command == "taxonomy":
        if args.input_fasta != None:
            taxonomy_single_genome(args)
        elif args.input_dir != None:
            taxonomy_multiple_genomes(args)

    elif args.command == "identify_bgcs":
        if args.input_fasta != None:
            identify_bgcs_single_genome(args)
        elif args.input_dir != None:
            identify_bgcs_multiple_genomes(args)

    elif args.command == "bgc_clustering":
        if args.clustering_type == "bigslice":
            run_bigslice_clustering(args)
        elif args.clustering_type == "bigscape":
            run_bigscape_clustering(args)
        elif args.clustering_type == "both":
            run_bigslice_clustering(args)
            run_bigscape_clustering(args)

    elif args.command == "all_submodules":
        run_all_submodules(args)


    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()