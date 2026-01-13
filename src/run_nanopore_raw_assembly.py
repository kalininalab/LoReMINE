import argparse
import os
from pathlib import Path
import shutil
from Bio import SeqIO
import subprocess
from src.evaluate_best_assembly import *

def perform_assembly_raw_reads_nanopore(args):

    """Perform the genome assembly using nanopore raw reads"""

    basepath = os.path.abspath(args.output)
    os.makedirs(basepath + "/assembly")
    os.makedirs(basepath + "/assembly/" + args.prefix)
    assembly_basepath = os.path.join(basepath, "assembly", args.prefix)
    os.makedirs(assembly_basepath + "/flye")
    os.makedirs(assembly_basepath + "/quast_outputs")

    if args.reads == None:
        args.reads = basepath + "/raw_reads/" + args.prefix + ".fastq"

    suffixes = Path(args.reads).suffixes
    raw_reads_filename = args.prefix + "".join(suffixes)

    if args.alt_param == False:   # Performing the assembly using default parameters
        if args.asm_coverage == None:
            flye_command = "flye --nano-raw " + basepath + "/raw_reads/" + raw_reads_filename + " -o " + assembly_basepath + "/flye/ -t " + args.threads + " -i 3 --scaffold -g " + args.genome_size
        else:
            flye_command = "flye --nano-raw " + basepath + "/raw_reads/" + raw_reads_filename + " -o " + assembly_basepath + "/flye/ -t " + args.threads + " -i 3 --scaffold -g " + args.genome_size + " --asm-coverage " + str(args.asm_coverage)

        print("Performing the genome assembly using Flye for " + args.prefix + "...")
        run_command(flye_command, verbosity=args.verbose)
        os.system("mv " + assembly_basepath + "/flye/assembly.fasta " + assembly_basepath + "/flye/assembly_temp.fasta")
        filter_low_quality_contigs(assembly_basepath + "/flye/assembly_temp.fasta", assembly_basepath + "/flye/assembly.fasta", args)
        os.remove(assembly_basepath + "/flye/assembly_temp.fasta")
        if os.path.isfile(assembly_basepath + "/flye/assembly.fasta"):
            print("\n\nAssembly completed. Collecting the assembly stats...\n\n")
            quast_command = "quast -o " + assembly_basepath + "/quast_outputs " + assembly_basepath + "/flye/assembly.fasta"   # Collect the assembly stats using quast
            run_command(quast_command, verbosity=args.verbose)
            generate_contig_circularity_info(assembly_basepath + "/quast_outputs", assembly_basepath + "/flye")
            print("Assembly stats can be found here: " + assembly_basepath + "/quast_outputs/report.pdf" + "\n")
            print("Contigs circularity stats can be found here: " + assembly_basepath + "/quast_outputs/circularity.tsv" + "\n\n")
            print("If you are not satisfied with the default assembly, please re-run the assembly with \"alt_param\" flag as \"True\"")
        else:
            print("There was an error while performing the assembly. Please look if you can solve it by yourself or send a log file to amayajaykumar.agrawal@helmholtz-hips.de")

    else:  # Performing the assembly using alternate parameters as the assembly using default parameters was not satisfactory
        min_overlap = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
        files_present = {}
        for ele in min_overlap:
            os.makedirs(assembly_basepath + "/flye/" + str(ele))
            if args.asm_coverage == None:
                flye_command = "flye --nano-raw " + basepath + "/raw_reads/" + raw_reads_filename + " -o " + assembly_basepath + "/flye/" + str(ele) + "/ -t " + args.threads + " -i 3 --scaffold -m " + str(ele) + " -g " + args.genome_size
            else:
                flye_command = "flye --nano-raw " + basepath + "/raw_reads/" + raw_reads_filename + " -o " + assembly_basepath + "/flye/" + str(ele) + "/ -t " + args.threads + " -i 3 --scaffold -m " + str(ele) + " -g " + args.genome_size + " --asm-coverage " + str(args.asm_coverage)

            print("Performing the genome assembly using Flye for " + args.prefix + " at min-overlap: " + str(ele) + "...")
            run_command(flye_command, verbosity=args.verbose)
            os.system("mv " + assembly_basepath + "/flye/" + str(ele) + "/assembly.fasta " + assembly_basepath + "/flye/" + str(ele) + "/assembly_temp.fasta")
            filter_low_quality_contigs(assembly_basepath + "/flye/" + str(ele) + "/assembly_temp.fasta", assembly_basepath + "/flye/" + str(ele) + "/assembly.fasta", args)
            os.remove(assembly_basepath + "/flye/" + str(ele) + "/assembly_temp.fasta")
            if os.path.isfile(assembly_basepath + "/flye/" + str(ele)+  "/assembly.fasta"):
                files_present[ele] = "yes"
                quast_command = "quast -o " + assembly_basepath + "/quast_outputs/" + str(ele) + "/ " + assembly_basepath + "/flye/" + str(ele) + "/assembly.fasta"   # Collect the assembly stats using quast
                run_command(quast_command, verbosity=args.verbose)
                generate_contig_circularity_info(assembly_basepath + "/quast_outputs/" + str(ele), assembly_basepath + "/flye/" + str(ele))
            else:
                continue

        print("\n\nGenome assembly completed using all alternate parameters. The assembly stats can be found here:\n\n")

        for key in files_present.keys():
            if files_present[key] == "yes":
                print(os.path.abspath(assembly_basepath + "/quast_outputs/" + str(key) + "/report.pdf"))

        print("Selecting the best assembly out of all the assemblies....")
        best_assembly = evaluate_assemblies(glob.glob(assembly_basepath + "/quast_outputs/*"), args.genome_size, args.weights, assembly_basepath + "/")   # Choosing the best assembly out of all the assemblies created using alternate parameters
        print("The best assembly (among all assemblies) according to us can be found here: " + assembly_basepath + "/flye/" + str(best_assembly) + "/assembly.fasta")

    if args.alt_param == False:
        final_assembly_path = assembly_basepath + "/flye/assembly.fasta"
        final_circularity_file_path = assembly_basepath + "/quast_outputs/circularity.tsv"
    else:
        final_assembly_path = assembly_basepath+ "/flye/" + str(best_assembly) + "/assembly.fasta"
        final_circularity_file_path = assembly_basepath + "/quast_outputs/" + str(best_assembly) + "/circularity.tsv"

    shutil.copyfile(final_assembly_path, assembly_basepath + "/" + args.prefix + ".fasta")
    print("Final assembly can be found here: " + assembly_basepath + "/" + args.prefix + ".fasta" + "\n\n")
    os.makedirs(assembly_basepath + "/checkm_output")
    os.makedirs(assembly_basepath + "/checkm_output/best_selected_assembly")
    shutil.copyfile(final_assembly_path, assembly_basepath + "/checkm_output/best_selected_assembly/" + args.prefix + ".fasta")
    print("\n\nRunning the contamination and completeness check on the best selected assembly\n\n")
    checkm_command = "checkm lineage_wf -x fasta -t " + args.threads + " " + assembly_basepath + "/checkm_output/best_selected_assembly/" + " " + assembly_basepath + "/checkm_output"
    run_command(checkm_command, verbosity=args.verbose)
    os.system("checkm qa " + assembly_basepath + "/checkm_output/lineage.ms " + assembly_basepath + "/checkm_output -o 1 -t " + args.threads + " --tab_table --file " + assembly_basepath + "/checkm_output/checkm_summary.tsv")
    print("\n\nYou can find the contamination and completeness stats of the best selected assembly here: " + assembly_basepath + "/checkm_output/checkm_summary.tsv\n\n")

    if args.alt_param == True:
        with open(assembly_basepath + "/checkm_output/checkm_summary.tsv") as checkm_file:
            next(checkm_file)
            for line in checkm_file:
                split_array = line.split('\t')
                completeness = str(split_array[-3].strip())
                contamination = str(split_array[-2].strip())

        with open(assembly_basepath + "/chosen_best_assembly.txt", 'a') as best_assembly_file:
            best_assembly_file.write("-" * len("The contamination and completeness stats of the best chosen assembly is given below:") + "\n")
            best_assembly_file.write("The contamination and completeness stats of the best chosen assembly is given below:\n")
            best_assembly_file.write("-" * len("The contamination and completeness stats of the best chosen assembly is given below:") + "\n")
            best_assembly_file.write("Completeness: " + completeness + '\n')
            best_assembly_file.write("Contamination: " + contamination + '\n')

    return final_assembly_path, final_circularity_file_path


def perform_assembly_raw_reads_nanopore_batch_run(args, raw_reads_path, genome_size, prefix, output_path):

    """Perform the genome assembly using nanopore raw reads (while running the whole pipeline together)"""

    basepath = output_path
    if not os.path.exists(basepath + "/assembly"):
        os.makedirs(basepath + "/assembly")
    os.makedirs(basepath + "/assembly/" + prefix)
    assembly_basepath = os.path.join(basepath, "assembly", prefix)
    os.makedirs(assembly_basepath + "/flye")
    os.makedirs(assembly_basepath + "/quast_outputs")

    if args.alt_param == False:  # Performing the assembly using default parameters
        if args.asm_coverage == None:
            flye_command = "flye --nano-raw " + raw_reads_path + " -o " + assembly_basepath + "/flye/ -t " + args.threads + " -i 3 --scaffold -g " + genome_size
        else:
            flye_command = "flye --nano-raw " + raw_reads_path + " -o " + assembly_basepath + "/flye/ -t " + args.threads + " -i 3 --scaffold -g " + genome_size + " --asm-coverage " + str(args.asm_coverage)

        print("Performing the genome assembly using Flye for " + prefix + "...")
        run_command(flye_command, verbosity=args.verbose)
        os.system("mv " + assembly_basepath + "/flye/assembly.fasta " + assembly_basepath + "/flye/assembly_temp.fasta")
        filter_low_quality_contigs(assembly_basepath + "/flye/assembly_temp.fasta", assembly_basepath + "/flye/assembly.fasta", args)
        os.remove(assembly_basepath + "/flye/assembly_temp.fasta")
        if os.path.isfile(assembly_basepath + "/flye/assembly.fasta"):
            print("\n\nAssembly completed. Collecting the assembly stats...\n\n")
            quast_command = "quast -o " + assembly_basepath + "/quast_outputs " + assembly_basepath + "/flye/assembly.fasta"   # Collect the assembly stats using quast
            run_command(quast_command, verbosity=args.verbose)
            generate_contig_circularity_info(assembly_basepath + "/quast_outputs", assembly_basepath + "/flye")
            print("Assembly stats can be found here: " + assembly_basepath + "/quast_outputs/report.pdf" + "\n")
            print("Contigs circularity stats can be found here: " + assembly_basepath + "/quast_outputs/circularity.tsv" + "\n\n")
            print("If you are not satisfied with the default assembly, please re-run the assembly with \"alt_param\" flag as \"True\"")
        else:
            print("There was an error while performing the assembly. Please look if you can solve it by yourself or send a log file to amayajaykumar.agrawal@helmholtz-hips.de")

    else:  # Performing the assembly using alternate parameters as the assembly using default parameters was not satisfactory
        min_overlap = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
        files_present = {}
        for ele in min_overlap:
            os.makedirs(assembly_basepath + "/flye/" + str(ele))
            if args.asm_coverage == None:
                flye_command = "flye --nano-raw " + raw_reads_path + " -o " + assembly_basepath + "/flye/" + str(ele) + "/ -t " + args.threads + " -i 3 --scaffold -m " + str(ele) + " -g " + genome_size
            else:
                flye_command = "flye --nano-raw " + raw_reads_path + " -o " + assembly_basepath + "/flye/" + str(ele) + "/ -t " + args.threads + " -i 3 --scaffold -m " + str(ele) + " -g " + genome_size + " --asm-coverage " + str(args.asm_coverage)

            print("Performing the genome assembly using Flye for " + prefix + " at min-overlap: " + str(ele) + "...")
            run_command(flye_command, verbosity=args.verbose)
            os.system("mv " + assembly_basepath + "/flye/" + str(ele) + "/assembly.fasta " + assembly_basepath + "/flye/" + str(ele) + "/assembly_temp.fasta")
            filter_low_quality_contigs(assembly_basepath + "/flye/" + str(ele) + "/assembly_temp.fasta", assembly_basepath + "/flye/" + str(ele) + "/assembly.fasta", args)
            os.remove(assembly_basepath + "/flye/" + str(ele) + "/assembly_temp.fasta")
            if os.path.isfile(assembly_basepath + "/flye/" + str(ele)+  "/assembly.fasta"):
                files_present[ele] = "yes"
                quast_command = "quast -o " + assembly_basepath + "/quast_outputs/" + str(ele) + "/ " + assembly_basepath + "/flye/" + str(ele) + "/assembly.fasta"   # Collect the assembly stats using quast
                run_command(quast_command, verbosity=args.verbose)
                generate_contig_circularity_info(assembly_basepath + "/quast_outputs/" + str(ele), assembly_basepath + "/flye/" + str(ele))
            else:
                continue

        print("\n\nGenome assembly completed using all alternate parameters. The assembly stats can be found here:\n\n")

        for key in files_present.keys():
            if files_present[key] == "yes":
                print(os.path.abspath(assembly_basepath + "/quast_outputs/" + str(key) + "/report.pdf"))

        print("Selecting the best assembly out of all the assemblies....")
        best_assembly = evaluate_assemblies(glob.glob(assembly_basepath + "/quast_outputs/*"), genome_size, args.weights, assembly_basepath + "/")  # Choosing the best assembly out of all the assemblies created using alternate parameters
        print("The best assembly (among all assemblies) according to us can be found here: " + assembly_basepath + "/flye/" + str(best_assembly) + "/assembly.fasta")

    if args.alt_param == False:
        final_assembly_path = assembly_basepath + "/flye/assembly.fasta"
        final_circularity_file_path = assembly_basepath + "/quast_outputs/circularity.tsv"
    else:
        final_assembly_path = assembly_basepath+ "/flye/" + str(best_assembly) + "/assembly.fasta"
        final_circularity_file_path = assembly_basepath + "/quast_outputs/" + str(best_assembly) + "/circularity.tsv"

    shutil.copyfile(final_assembly_path, assembly_basepath + "/" + prefix + ".fasta")
    print("Final assembly can be found here: " + assembly_basepath + "/" + prefix + ".fasta" + "\n\n")
    os.makedirs(assembly_basepath + "/checkm_output")
    os.makedirs(assembly_basepath + "/checkm_output/best_selected_assembly")
    shutil.copyfile(final_assembly_path, assembly_basepath + "/checkm_output/best_selected_assembly/" + prefix + ".fasta")
    print("\n\nRunning the contamination and completeness check on the best selected assembly\n\n")
    checkm_command = "checkm lineage_wf -x fasta -t " + args.threads + " " + assembly_basepath + "/checkm_output/best_selected_assembly/" + " " + assembly_basepath + "/checkm_output"
    run_command(checkm_command, verbosity=args.verbose)
    os.system("checkm qa " + assembly_basepath + "/checkm_output/lineage.ms " + assembly_basepath + "/checkm_output -o 1 -t " + args.threads + " --tab_table --file " + assembly_basepath + "/checkm_output/checkm_summary.tsv")
    print("\n\nYou can find the contamination and completeness stats of the best selected assembly here: " + assembly_basepath + "/checkm_output/checkm_summary.tsv\n\n")

    if args.alt_param == True:
        with open(assembly_basepath + "/checkm_output/checkm_summary.tsv") as checkm_file:
            next(checkm_file)
            for line in checkm_file:
                split_array = line.split('\t')
                completeness = str(split_array[-3].strip())
                contamination = str(split_array[-2].strip())

        with open(assembly_basepath + "/chosen_best_assembly.txt", 'a') as best_assembly_file:
            best_assembly_file.write("-" * len("The contamination and completeness stats of the best chosen assembly is given below:") + "\n")
            best_assembly_file.write("The contamination and completeness stats of the best chosen assembly is given below:\n")
            best_assembly_file.write("-" * len("The contamination and completeness stats of the best chosen assembly is given below:") + "\n")
            best_assembly_file.write("Completeness: " + completeness + '\n')
            best_assembly_file.write("Contamination: " + contamination + '\n')

    return final_assembly_path, final_circularity_file_path


def generate_contig_circularity_info(quast_path, output_path):

    """Checking the assemblies to identify the length of the contigs along with their circularity info (whether the contig is circular or linear)"""

    write_file = open(quast_path + '/circularity.tsv', 'x')
    write_file.write("contig_name\tcircular\tLength\n")

    headers = [record.id for record in SeqIO.parse(output_path + "/assembly.fasta", "fasta")]

    with open(output_path + "/assembly_info.txt") as assembly_info_file:
        counter = 0
        for line in assembly_info_file.readlines():
            counter += 1
            if counter == 1:
                continue
            split_array = line.split()

            if split_array[0].strip() not in headers:
                continue

            if split_array[3].strip() == "Y":
                write_file.write(split_array[0].strip() + '\tyes' + '\t' + str(split_array[1].strip()) + '\n')
            elif split_array[3].strip() == "N":
                write_file.write(split_array[0].strip() + '\tno' + '\t' + str(split_array[1].strip()) + '\n')

def filter_low_quality_contigs(input_fasta, output_fasta, args):

    """Checking the assemblies to filter out low quality contigs like homopolymer runs and SSRs"""

    masked_path = str(Path(input_fasta).resolve().parent)
    masked_file = masked_path + '/' + str(Path(input_fasta).resolve().stem) + "_masked.fasta"
    run_command("dustmasker -in " + str(input_fasta) + " -outfmt fasta -out " + str(masked_file), verbosity=args.verbose)
    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(masked_file, "fasta"):
            low_record_quality = is_low_quality(record)
            if low_record_quality == "No":
                for input_fasta_record in SeqIO.parse(input_fasta, "fasta"):
                    if input_fasta_record.id == record.id:
                        SeqIO.write(input_fasta_record, out, "fasta")
                    else:
                        continue
            else:
                continue
    os.remove(masked_file)

def is_low_quality(record):

    """Checking if >=80% of the contig is marked as low quality by dustmasker. If yes, then we will discard it"""

    low_quality = ""
    seq = str(record.seq)
    total = len(seq)
    lowercase = sum(1 for b in seq if b.islower())
    lower_fraction = lowercase / total

    if lower_fraction >= 0.8:
        low_quality = "Yes"
    else:
        low_quality = "No"

    return low_quality

def run_command(cmd, verbosity=0):
    """
    Run a shell command with controlled verbosity
    """
    if verbosity == 0:
        subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    else:
        subprocess.run(cmd, shell=True, check=True)