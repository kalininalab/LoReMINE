import argparse
import os
from Bio import SeqIO
import shutil
import glob
from pathlib import Path
import subprocess
from src.evaluate_best_assembly import *

def perform_assembly_hifi_reads(args):

    """Perform the genome assembly using pacbio hifi reads"""

    basepath = os.path.abspath(args.output)
    os.makedirs(basepath + "/assembly")
    os.makedirs(basepath + "/assembly/" + args.prefix)
    assembly_basepath = os.path.join(basepath, "assembly", args.prefix)
    os.makedirs(assembly_basepath + "/quast_outputs")
    os.makedirs(assembly_basepath + "/flye")
    print("Performing the genome assembly using Flye for " + args.prefix + "...")

    flye_command = "flye --pacbio-hifi " + basepath + "/raw_reads/" + args.prefix + ".fastq -o " + assembly_basepath + "/flye/ -t " + args.threads + " -i 3 --scaffold --asm-coverage 50 -g " + args.genome_size  # Performing the assembly using flye
    run_command(flye_command, verbosity=args.verbose)
    os.system("mv " + assembly_basepath + "/flye/assembly.fasta " + assembly_basepath + "/flye/assembly_temp.fasta")
    filter_low_quality_contigs(assembly_basepath + "/flye/assembly_temp.fasta", assembly_basepath + "/flye/assembly.fasta", args)
    os.remove(assembly_basepath + "/flye/assembly_temp.fasta")
    if os.path.isfile(assembly_basepath + "/flye/assembly.fasta"):
        quast_command = "quast -o " + assembly_basepath + "/quast_outputs/flye " + assembly_basepath + "/flye/assembly.fasta"   # Collect the assembly stats using quast
        run_command(quast_command, verbosity=args.verbose)
        generate_contig_circularity_info("flye", assembly_basepath + "/quast_outputs/flye/", assembly_basepath + "/flye/")
        print("\n\nAssembly using flye completed. Next, we will perform assembly using HiFiasm for " + args.prefix + "...")
    else:
        print("There was an error while performing the assembly. Please look if you can solve it by yourself or send a log file to amayajaykumar.agrawal@helmholtz-hips.de")

    os.makedirs(assembly_basepath + "/Hifiasm")
    os.chdir(assembly_basepath + "/Hifiasm/")
    genome_size_intermediate = int(args.genome_size)/1000000
    updated_genome_size = str(genome_size_intermediate) + 'm'
    hifiasm_command = "hifiasm -o " + args.prefix + " -t " + args.threads + " -a 8 -f 0 --hg-size " + str(updated_genome_size) + " " + basepath + "/raw_reads/" + args.prefix + ".fastq"  # Performing the assembly using Hifiasm
    run_command(hifiasm_command, verbosity=args.verbose)
    os.system('''awk '/^S/{print \">\"$2;print $3}\' ''' + args.prefix + ".bp.p_ctg.gfa" + " > final_assembly_temp.fasta")
    filter_low_quality_contigs(assembly_basepath + "/Hifiasm/final_assembly_temp.fasta", assembly_basepath + "/Hifiasm/final_assembly.fasta", args)
    os.remove(assembly_basepath + "/Hifiasm/final_assembly_temp.fasta")
    if os.path.isfile(assembly_basepath + "/Hifiasm/final_assembly.fasta"):
        quast_command = "quast -o " + assembly_basepath + "/quast_outputs/Hifiasm " + assembly_basepath + "/Hifiasm/final_assembly.fasta"  # Collect the assembly stats using quast
        run_command(quast_command, verbosity=args.verbose)
        generate_contig_circularity_info("Hifiasm", assembly_basepath + "/quast_outputs/Hifiasm/", assembly_basepath + "/Hifiasm/")
        print("\n\nAssembly using HiFiasm completed. Next, we will perform assembly using IPA for " + args.prefix + "...")
    else:
        print("There was an error while performing the assembly. Please look if you can solve it by yourself or send a log file to xxxx@gmail.com")

    os.makedirs(assembly_basepath + "/IPA")
    ipa_command = "ipa local -i " + basepath + "/raw_reads/" + args.prefix + ".fastq" + " --run-dir " + assembly_basepath + "/IPA/" + " --tmp-dir "  + assembly_basepath + "/IPA/" + " --njobs 1 --nthreads " + args.threads  # Performing the assembly using IPA
    run_command(ipa_command, verbosity=args.verbose)
    shutil.copyfile(assembly_basepath + "/IPA/19-final/final.p_ctg.fasta", assembly_basepath + "/IPA/final_assembly_temp.fasta")
    filter_low_quality_contigs(assembly_basepath + "/IPA/final_assembly_temp.fasta", assembly_basepath + "/IPA/final_assembly.fasta", args)
    os.remove(assembly_basepath + "/IPA/final_assembly_temp.fasta")
    if os.path.isfile(assembly_basepath + "/IPA/final_assembly.fasta"):
        quast_command = "quast -o " + assembly_basepath + "/quast_outputs/IPA " + assembly_basepath + "/IPA/final_assembly.fasta"   # Collect the assembly stats using quast
        run_command(quast_command, verbosity=args.verbose)
        generate_contig_circularity_info("ipa", assembly_basepath + "/quast_outputs/IPA/", assembly_basepath + "/IPA/")
        print("\n\nAssembly using IPA completed\n\n")
    else:
        print("There was an error while performing the assembly. Please look if you can solve it by yourself or send a log file to xxxx@gmail.com")

    print("Assembly stats for the flye assembly can be found here: " + assembly_basepath + "/quast_outputs/flye/report.pdf")
    print("Assembly stats for the HiFiasm assembly can be found here: " + assembly_basepath + "/quast_outputs/Hifiasm/report.pdf")
    print("Assembly stats for the IPA assembly can be found here: " + assembly_basepath + "/quast_outputs/IPA/report.pdf")
    print("\n\n")
    print("Selecting the best assembly out of all the assemblies....")
    best_assembly = evaluate_assemblies(glob.glob(assembly_basepath + "/quast_outputs/*"), args.genome_size, args.weights, assembly_basepath + '/')   # Choosing the best assembly out of all the 3 assemblies created using flye, Hifiasm and IPA

    final_assembly_path = ""
    final_circularity_file_path = ""
    if best_assembly == "flye":
        print("The best assembly (among all assemblies) according to us can be found here: " + assembly_basepath + '/' + str(best_assembly) + "/assembly.fasta")
        final_assembly_path = assembly_basepath + '/' + str(best_assembly) + "/assembly.fasta"
        final_circularity_file_path = assembly_basepath + "/quast_outputs/" + str(best_assembly) + "/circularity.tsv"
    else:
        print("The best assembly (among all assemblies) according to us can be found here: " + assembly_basepath + '/' + str(best_assembly) + "/final_assembly.fasta")
        final_assembly_path = assembly_basepath + '/' + str(best_assembly) + "/final_assembly.fasta"
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

    with open(assembly_basepath + "/checkm_output/checkm_summary.tsv") as checkm_file:
        next(checkm_file)
        for line in checkm_file:
            split_array = line.split('\t')
            completeness = str(split_array[-3].strip())
            contamination = str(split_array[-2].strip())

    with open(assembly_basepath + "/chosen_best_assembly.txt", 'a') as best_assembly_file:
        best_assembly_file.write("\n\nThe contamination and completeness stats of the best chosen assembly is given below:\n\n")
        best_assembly_file.write("Completeness: " + completeness + '\n')
        best_assembly_file.write("Contamination: " + contamination + '\n')

    return final_assembly_path, final_circularity_file_path


def perform_assembly_hifi_reads_batch_run(args, raw_reads_path, genome_size, prefix, output_path):

    """Perform the genome assembly using pacbio hifi reads (while running the whole pipeline together)"""

    basepath = output_path
    if not os.path.exists(basepath + "/assembly"):
        os.makedirs(basepath + "/assembly")
    os.makedirs(basepath + "/assembly/" + prefix)
    assembly_basepath = os.path.join(basepath, "assembly", prefix)
    os.makedirs(assembly_basepath + "/quast_outputs")
    os.makedirs(assembly_basepath + "/flye")
    print("Performing the genome assembly using Flye for " + prefix + "...")

    flye_command = "flye --pacbio-hifi " + raw_reads_path + " -o " + assembly_basepath + "/flye/ -t " + args.threads + " -i 3 --scaffold --asm-coverage 50 -g " + genome_size  # Performing the assembly using flye
    run_command(flye_command, verbosity=args.verbose)
    os.system("mv " + assembly_basepath + "/flye/assembly.fasta " + assembly_basepath + "/flye/assembly_temp.fasta")
    filter_low_quality_contigs(assembly_basepath + "/flye/assembly_temp.fasta", assembly_basepath + "/flye/assembly.fasta", args)
    os.remove(assembly_basepath + "/flye/assembly_temp.fasta")
    if os.path.isfile(assembly_basepath + "/flye/assembly.fasta"):
        quast_command = "quast -o " + assembly_basepath + "/quast_outputs/flye " + assembly_basepath + "/flye/assembly.fasta"  # Collect the assembly stats using quast
        run_command(quast_command, verbosity=args.verbose)
        generate_contig_circularity_info("flye", assembly_basepath + "/quast_outputs/flye/", assembly_basepath + "/flye/")
        print("\n\nAssembly using flye completed. Next, we will perform assembly using HiFiasm for " + prefix + "...")
    else:
        print("There was an error while performing the assembly. Please look if you can solve it by yourself or send a log file to xxxx@gmail.com")

    os.makedirs(assembly_basepath + "/Hifiasm")
    os.chdir(assembly_basepath + "/Hifiasm/")
    genome_size_intermediate = int(genome_size)/1000000
    updated_genome_size = str(genome_size_intermediate) + 'm'
    hifiasm_command = "hifiasm -o " + prefix + " -t " + args.threads + " -a 8 -f 0 --hg-size " + str(updated_genome_size) + " " + raw_reads_path  # Performing the assembly using Hifiasm
    run_command(hifiasm_command, verbosity=args.verbose)
    os.system('''awk '/^S/{print \">\"$2;print $3}\' ''' + prefix + ".bp.p_ctg.gfa" + " > final_assembly_temp.fasta")
    filter_low_quality_contigs(assembly_basepath + "/Hifiasm/final_assembly_temp.fasta", assembly_basepath + "/Hifiasm/final_assembly.fasta", args)
    os.remove(assembly_basepath + "/Hifiasm/final_assembly_temp.fasta")
    if os.path.isfile(assembly_basepath + "/Hifiasm/final_assembly.fasta"):
        quast_command = "quast -o " + assembly_basepath + "/quast_outputs/Hifiasm " + assembly_basepath + "/Hifiasm/final_assembly.fasta"  # Collect the assembly stats using quast
        run_command(quast_command, verbosity=args.verbose)
        generate_contig_circularity_info("Hifiasm", assembly_basepath + "/quast_outputs/Hifiasm/", assembly_basepath + "/Hifiasm/")
        print("\n\nAssembly using HiFiasm completed. Next, we will perform assembly using IPA for " + prefix + "...")
    else:
        print("There was an error while performing the assembly. Please look if you can solve it by yourself or send a log file to amayajaykumar.agrawal@helmholtz-hips.de")

    os.makedirs(assembly_basepath + "/IPA")
    ipa_command = "ipa local -i " + raw_reads_path + " --run-dir " + assembly_basepath + "/IPA/" + " --tmp-dir "  + assembly_basepath + "/IPA/" + " --njobs 1 --nthreads " + args.threads  # Performing the assembly using IPA
    run_command(ipa_command, verbosity=args.verbose)
    shutil.copyfile(assembly_basepath + "/IPA/19-final/final.p_ctg.fasta", assembly_basepath + "/IPA/final_assembly_temp.fasta")
    filter_low_quality_contigs(assembly_basepath + "/IPA/final_assembly_temp.fasta", assembly_basepath + "/IPA/final_assembly.fasta", args)
    os.remove(assembly_basepath + "/IPA/final_assembly_temp.fasta")
    if os.path.isfile(assembly_basepath + "/IPA/final_assembly.fasta"):
        quast_command = "quast -o " + assembly_basepath + "/quast_outputs/IPA " + assembly_basepath + "/IPA/final_assembly.fasta"  # Collect the assembly stats using quast
        run_command(quast_command, verbosity=args.verbose)
        generate_contig_circularity_info("ipa", assembly_basepath + "/quast_outputs/IPA/", assembly_basepath + "/IPA/")
        print("\n\nAssembly using IPA completed\n\n")
    else:
        print("There was an error while performing the assembly. Please look if you can solve it by yourself or send a log file to xxxx@gmail.com")

    print("Assembly stats for the flye assembly can be found here: " + assembly_basepath + "/quast_outputs/flye/report.pdf")
    print("Assembly stats for the HiFiasm assembly can be found here: " + assembly_basepath + "/quast_outputs/Hifiasm/report.pdf")
    print("Assembly stats for the IPA assembly can be found here: " + assembly_basepath + "/quast_outputs/IPA/report.pdf")
    print("\n\n")
    print("Selecting the best assembly out of all the assemblies....")
    best_assembly = evaluate_assemblies(glob.glob(assembly_basepath + "/quast_outputs/*"), genome_size, args.weights, assembly_basepath + '/')   # Choosing the best assembly out of all the 3 assemblies created using flye, Hifiasm and IPA

    final_assembly_path = ""
    final_circularity_file_path = ""
    if best_assembly == "flye":
        print("The best assembly (among all assemblies) according to us can be found here: " + assembly_basepath + '/' + str(best_assembly) + "/assembly.fasta")
        final_assembly_path = assembly_basepath + '/' + str(best_assembly) + "/assembly.fasta"
        final_circularity_file_path = assembly_basepath + "/quast_outputs/" + str(best_assembly) + "/circularity.tsv"
    else:
        print("The best assembly (among all assemblies) according to us can be found here: " + assembly_basepath + '/' + str(best_assembly) + "/final_assembly.fasta")
        final_assembly_path = assembly_basepath + '/' + str(best_assembly) + "/final_assembly.fasta"
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

    with open(assembly_basepath + "/checkm_output/checkm_summary.tsv") as checkm_file:
        next(checkm_file)
        for line in checkm_file:
            split_array = line.split('\t')
            completeness = str(split_array[-3].strip())
            contamination = str(split_array[-2].strip())

    with open(assembly_basepath + "/chosen_best_assembly.txt", 'a') as best_assembly_file:
        best_assembly_file.write("\n\nThe contamination and completeness stats of the best chosen assembly is given below:\n\n")
        best_assembly_file.write("Completeness: " + completeness + '\n')
        best_assembly_file.write("Contamination: " + contamination + '\n')

    return final_assembly_path, final_circularity_file_path

def generate_contig_circularity_info(assembler, quast_path, output_path):

    """Checking the assemblies to identify the length of the contigs along with their circularity info (whether the contig is circular or linear)"""

    write_file = open(quast_path + 'circularity.tsv', 'x')
    write_file.write("contig_name\tcircular\tLength\n")

    if assembler == "flye":
        headers = [record.id for record in SeqIO.parse(output_path + "/assembly.fasta", "fasta")]

        with open(output_path + "assembly_info.txt") as assembly_info_file:
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

    elif assembler == "Hifiasm":
        fasta_file = output_path + "final_assembly.fasta"
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id[-1] == 'c':
                write_file.write(record.id + '\tyes' + '\t' + str(len(record.seq)) + '\n')
            elif record.id[-1] == 'l':
                write_file.write(record.id + '\tno' + '\t' + str(len(record.seq)) + '\n')

    elif assembler == "ipa":
        fasta_file = output_path + "final_assembly.fasta"
        for record in SeqIO.parse(fasta_file, "fasta"):
            split_array = record.id.split('/')
            if split_array[2].strip() == 'c':
                write_file.write(record.id + '\tyes' + '\t' + str(len(record.seq)) + '\n')
            elif split_array[2].strip() == 'l':
                write_file.write(record.id + '\tno' + '\t' + str(len(record.seq)) + '\n')


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