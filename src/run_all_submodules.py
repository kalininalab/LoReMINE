"""This file is used to run all the submodules of LoReMINE pipeline together"""

from src.convert_reads_to_fastq import *
from src.run_pacbio_raw_assembly import *
from src.run_nanopore_raw_assembly import *
from src.run_pacbio_hifi_assembly import *
from src.evaluate_best_assembly import *
from src.identify_taxonomy import *
from src.identify_bgcs import *
from src.run_bgc_clustering import *
from pathlib import Path
import shutil
import os
import glob
from Bio import SeqIO
import shutil
import sys

def run_all_submodules(args):

    final_assembly_path = ""
    basepath = os.path.abspath(args.output)

    if args.batch_run == None:

        """Running the whole pipeline for a single strain"""

        if args.pacbio_raw != None: # If the input reads are in .bam file, then converting it to fastq
            print('Pacbio raw reads provided at this location:' + os.path.abspath(args.pacbio_raw))
            convert_reads(args)
            print("\n\nStaring the assembly using pacbio raw reads now.....\n\n")
            final_assembly_path, final_circularity_file_path = perform_assembly_raw_reads_pacbio(args)

        elif args.pacbio_hifi != None: # If the input reads are in .bam file, then converting it to fastq
            print('Pacbio HiFi reads provided at this location:' + os.path.abspath(args.pacbio_hifi))
            convert_reads(args)
            print("\n\nStaring the assembly using pacbio HiFi reads now.....\n\n")
            final_assembly_path, final_circularity_file_path = perform_assembly_hifi_reads(args)

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
            final_assembly_path, final_circularity_file_path = perform_assembly_raw_reads_pacbio(args)

        elif args.reads != None and args.reads_type == "hifi_pacbio":
            print('Pacbio HiFi reads provided at this location:' + os.path.abspath(args.reads))
            os.makedirs(args.output + "/raw_reads")
            suffixes = Path(args.reads).suffixes
            raw_reads_filename = args.prefix + "".join(suffixes)
            shutil.copyfile(args.reads, args.output + "/raw_reads/" + raw_reads_filename)
            print("\n\nStaring the assembly using pacbio HiFi reads now.....\n\n")
            final_assembly_path, final_circularity_file_path = perform_assembly_hifi_reads(args)

        elif args.reads != None and args.reads_type == "raw_nanopore":
            print('Nanopore raw reads provided at this location:' + os.path.abspath(args.reads))
            os.makedirs(args.output + "/raw_reads")
            suffixes = Path(args.reads).suffixes
            raw_reads_filename = args.prefix + "".join(suffixes)
            shutil.copyfile(args.reads, args.output + "/raw_reads/" + raw_reads_filename)
            print("\n\nStaring the assembly using nanopore raw reads now.....\n\n")
            final_assembly_path, final_circularity_file_path = perform_assembly_raw_reads_nanopore(args)

        if not os.path.exists(basepath + '/assembly/' + args.prefix + '/best_assembly'):
            os.makedirs(basepath + '/assembly/' + args.prefix + '/best_assembly')
        readjusted_assembly_header_path = basepath + '/assembly/' + args.prefix + '/best_assembly/' + args.prefix + '.fasta'
        reheader_fasta(final_assembly_path, readjusted_assembly_header_path, final_circularity_file_path, basepath + '/assembly/' + args.prefix + '/best_assembly/', args.prefix)

        taxonomy_single_genome_all_submodules(args, readjusted_assembly_header_path, basepath)
        bgc_output_path = identify_bgcs_single_genome_all_submodules(args, readjusted_assembly_header_path, basepath)


        if args.clustering_type == "bigslice":
            run_bigslice_clustering_all_submodules(args, bgc_output_path, basepath)
        elif args.clustering_type == "bigscape":
            run_bigscape_clustering_all_submodules(args, bgc_output_path, basepath)
        elif args.clustering_type == "both":
            run_bigslice_clustering_all_submodules(args, bgc_output_path, basepath)
            run_bigscape_clustering_all_submodules(args, bgc_output_path, basepath)

    else:

        """Running the whole pipeline for a batch of strains"""

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
                    if not os.path.exists(basepath + "/raw_reads"):
                        os.makedirs(basepath + "/raw_reads")
                    suffixes = Path(raw_reads_location).suffixes
                    raw_reads_filename = prefix + "".join(suffixes)
                    shutil.copyfile(raw_reads_location, basepath + "/raw_reads/" + raw_reads_filename)
                    print("\n\nStaring the assembly using pacbio HiFi reads now.....\n\n")
                    final_assembly_path, final_circularity_file_path = perform_assembly_hifi_reads_batch_run(args, basepath + "/raw_reads/" + prefix + ".fastq", genome_size, prefix, basepath)

                elif reads_type == "raw_pacbio":
                    print('Pacbio raw reads provided at this location:' + raw_reads_location)
                    if not os.path.exists(basepath + "/raw_reads"):
                        os.makedirs(basepath + "/raw_reads")
                    suffixes = Path(raw_reads_location).suffixes
                    raw_reads_filename = prefix + "".join(suffixes)
                    shutil.copyfile(raw_reads_location, basepath + "/raw_reads/" + raw_reads_filename)
                    print("\n\nStaring the assembly using pacbio raw reads now.....\n\n")
                    final_assembly_path, final_circularity_file_path = perform_assembly_raw_reads_pacbio_batch_run(args, basepath + "/raw_reads/" + raw_reads_filename, genome_size, prefix, basepath)

                elif reads_type == "raw_nanopore":
                    print('Nanopore raw reads provided at this location:' + raw_reads_location)
                    if not os.path.exists(basepath + "/raw_reads"):
                        os.makedirs(basepath + "/raw_reads")
                    suffixes = Path(raw_reads_location).suffixes
                    raw_reads_filename = prefix + "".join(suffixes)
                    shutil.copyfile(raw_reads_location, basepath + "/raw_reads/" + raw_reads_filename)
                    print("\n\nStaring the assembly using nanopore raw reads now.....\n\n")
                    final_assembly_path, final_circularity_file_path = perform_assembly_raw_reads_nanopore_batch_run(args, basepath + "/raw_reads/" + raw_reads_filename, genome_size, prefix, basepath)

                if not os.path.exists(basepath + '/assembly/best_assemblies'):
                    os.makedirs(basepath + '/assembly/best_assemblies')
                readjusted_assembly_header_path = basepath + '/assembly/best_assemblies/' + prefix + '.fasta'
                reheader_fasta(final_assembly_path, readjusted_assembly_header_path, final_circularity_file_path, basepath + '/assembly/best_assemblies/', prefix)

        taxonomy_multiple_genomes_all_submodules(args, basepath + '/assembly/best_assemblies/', basepath)
        all_bgcs_path = identify_bgcs_multiple_genomes_all_submodules(args, basepath + '/assembly/best_assemblies/', basepath)



        if args.clustering_type == "bigslice":
            run_bigslice_clustering_all_submodules(args, basepath + all_bgcs_path, basepath)
        elif args.clustering_type == "bigscape":
            run_bigscape_clustering_all_submodules(args, basepath + all_bgcs_path, basepath)
        elif args.clustering_type == "both":
            run_bigslice_clustering_all_submodules(args, basepath + all_bgcs_path, basepath)
            run_bigscape_clustering_all_submodules(args, basepath + all_bgcs_path, basepath)


def reheader_fasta(input_fasta, output_fasta, final_circularity_file_path, best_assemblies_path, prefix):

    """Rename the contig headers of a best selected genome assembly file so that all contigs are numbered sequentially with "prefix" also added to each contig header"""

    mapping = {}
    with open(output_fasta, "w") as output_file:
        for i, record in enumerate(SeqIO.parse(input_fasta, "fasta"), start=1):
            new_id = f"contig_{i}_{prefix}"
            mapping[record.id] = new_id
            record.id = new_id
            record.description = ""
            SeqIO.write(record, output_file, "fasta")

    if not os.path.exists(best_assemblies_path + "/quast_output"):
        os.makedirs(best_assemblies_path + "/quast_output")

    os.system("quast -o " + best_assemblies_path + "/quast_output/" + prefix + "/ " + output_fasta)

    updated_circularity_file = open(best_assemblies_path + "/quast_output/" + prefix + "/circularity.tsv", 'x')

    counter = 0
    with open(final_circularity_file_path) as circularity_file:
        for line in circularity_file.readlines():
            counter += 1
            if counter == 1:
                updated_circularity_file.write(line.strip() + '\n')
                continue
            else:
                split_array = line.split('\t')
                updated_circularity_file.write(mapping[split_array[0].strip()] + '\t' + split_array[1].strip() + '\t' + split_array[2].strip() + '\n')