import argparse
import os
from Bio import SeqIO
import shutil
import glob
from src.evaluate_best_assembly import *

def perform_assembly_hifi_reads(args):
    basepath = os.path.abspath(args.output)
    os.makedirs(basepath + "/assembly")
    os.makedirs(basepath + "/assembly/" + args.prefix)
    assembly_basepath = os.path.join(basepath, "assembly", args.prefix)
    os.makedirs(assembly_basepath + "/quast_outputs")
    os.makedirs(assembly_basepath + "/flye")
    print("Performing the genome assembly using Flye...")

    flye_command = "flye --pacbio-hifi " + basepath + "/raw_reads/" + args.prefix + ".fastq -o " + assembly_basepath + "/flye/ -t " + args.threads + " -i 3 --scaffold --asm-coverage 50 -g " + args.genome_size
    os.system(flye_command)

    if os.path.isfile(assembly_basepath + "/flye/assembly.fasta"):
        quast_command = "quast -o " + assembly_basepath + "/quast_outputs/flye " + assembly_basepath + "/flye/assembly.fasta"
        os.system(quast_command)
        generate_contig_circularity_info("flye", assembly_basepath + "/quast_outputs/flye/", assembly_basepath + "/flye/")
        print("\n\nAssembly using flye completed. Next, we will perform assembly using HiFiasm...")
    else:
        print("There was an error while performing the assembly. Please look if you can solve it by yourself or send a log file to xxxx@gmail.com")

    os.makedirs(assembly_basepath + "/Hifiasm")
    os.chdir(assembly_basepath + "/Hifiasm/")
    genome_size_intermediate = int(args.genome_size)/1000000
    updated_genome_size = str(genome_size_intermediate) + 'm'
    hifiasm_command = "hifiasm -o " + args.prefix + " -t " + args.threads + " -a 8 -f 0 --hg-size " + str(updated_genome_size) + " " + basepath + "/raw_reads/" + args.prefix + ".fastq"
    os.system(hifiasm_command)
    os.system('''awk '/^S/{print \">\"$2;print $3}\' ''' + args.prefix + ".bp.p_ctg.gfa" + " > final_assembly.fasta")
    if os.path.isfile(assembly_basepath + "/Hifiasm/final_assembly.fasta"):
        quast_command = "quast -o " + assembly_basepath + "/quast_outputs/Hifiasm " + assembly_basepath + "/Hifiasm/final_assembly.fasta"
        os.system(quast_command)
        generate_contig_circularity_info("Hifiasm", assembly_basepath + "/quast_outputs/Hifiasm/", assembly_basepath + "/Hifiasm/")
        print("\n\nAssembly using HiFiasm completed. Next, we will perform assembly using IPA...")
    else:
        print("There was an error while performing the assembly. Please look if you can solve it by yourself or send a log file to xxxx@gmail.com")

    os.makedirs(assembly_basepath + "/IPA")
    ipa_command = "ipa local -i " + basepath + "/raw_reads/" + args.prefix + ".fastq" + " --run-dir " + assembly_basepath + "/IPA/" + " --tmp-dir "  + assembly_basepath + "/IPA/" + " --njobs 1 --nthreads " + args.threads
    os.system(ipa_command)
    shutil.copyfile(assembly_basepath + "/IPA/19-final/final.p_ctg.fasta", assembly_basepath + "/IPA/final_assembly.fasta")
    if os.path.isfile(assembly_basepath + "/IPA/final_assembly.fasta"):
        quast_command = "quast -o " + assembly_basepath + "/quast_outputs/IPA " + assembly_basepath + "/IPA/final_assembly.fasta"
        os.system(quast_command)
        generate_contig_circularity_info("ipa", assembly_basepath + "/quast_outputs/IPA/", assembly_basepath + "/IPA/")
        print("\n\nAssembly using IPA completed\n\n")
    else:
        print("There was an error while performing the assembly. Please look if you can solve it by yourself or send a log file to xxxx@gmail.com")

    print("Assembly stats for the flye assembly can be found here: " + assembly_basepath + "/quast_outputs/flye/report.pdf")
    print("Assembly stats for the HiFiasm assembly can be found here: " + assembly_basepath + "/quast_outputs/Hifiasm/report.pdf")
    print("Assembly stats for the IPA assembly can be found here: " + assembly_basepath + "/quast_outputs/IPA/report.pdf")
    print("\n\n")
    print("Selecting the best assembly out of all the assemblies....")
    best_assembly = evaluate_assemblies(glob.glob(assembly_basepath + "/quast_outputs/*"), args.genome_size, assembly_basepath + '/')

    final_assembly_path = ""
    if best_assembly == "flye":
        print("The best assembly (among all assemblies) according to us can be found here: " + assembly_basepath + '/' + str(best_assembly) + "/assembly.fasta")
        final_assembly_path = assembly_basepath + '/' + str(best_assembly) + "/assembly.fasta"
    else:
        print("The best assembly (among all assemblies) according to us can be found here: " + assembly_basepath + '/' + str(best_assembly) + "/final_assembly.fasta")
        final_assembly_path = assembly_basepath + '/' + str(best_assembly) + "/final_assembly.fasta"

    return final_assembly_path


def perform_assembly_hifi_reads_batch_run(args, raw_reads_path, genome_size, prefix, output_path):
    basepath = output_path
    if not os.path.exists(basepath + "/assembly"):
        os.makedirs(basepath + "/assembly")
    os.makedirs(basepath + "/assembly/" + prefix)
    assembly_basepath = os.path.join(basepath, "assembly", prefix)
    os.makedirs(assembly_basepath + "/quast_outputs")
    os.makedirs(assembly_basepath + "/flye")
    print("Performing the genome assembly using Flye...")

    flye_command = "flye --pacbio-hifi " + raw_reads_path + " -o " + assembly_basepath + "/flye/ -t " + args.threads + " -i 3 --scaffold --asm-coverage 50 -g " + genome_size
    os.system(flye_command)

    if os.path.isfile(assembly_basepath + "/flye/assembly.fasta"):
        quast_command = "quast -o " + assembly_basepath + "/quast_outputs/flye " + assembly_basepath + "/flye/assembly.fasta"
        os.system(quast_command)
        generate_contig_circularity_info("flye", assembly_basepath + "/quast_outputs/flye/", assembly_basepath + "/flye/")
        print("\n\nAssembly using flye completed. Next, we will perform assembly using HiFiasm...")
    else:
        print("There was an error while performing the assembly. Please look if you can solve it by yourself or send a log file to xxxx@gmail.com")

    os.makedirs(assembly_basepath + "/Hifiasm")
    os.chdir(assembly_basepath + "/Hifiasm/")
    genome_size_intermediate = int(genome_size)/1000000
    updated_genome_size = str(genome_size_intermediate) + 'm'
    hifiasm_command = "hifiasm -o " + prefix + " -t " + args.threads + " -a 8 -f 0 --hg-size " + str(updated_genome_size) + " " + raw_reads_path
    os.system(hifiasm_command)
    os.system('''awk '/^S/{print \">\"$2;print $3}\' ''' + prefix + ".bp.p_ctg.gfa" + " > final_assembly.fasta")
    if os.path.isfile(assembly_basepath + "/Hifiasm/final_assembly.fasta"):
        quast_command = "quast -o " + assembly_basepath + "/quast_outputs/Hifiasm " + assembly_basepath + "/Hifiasm/final_assembly.fasta"
        os.system(quast_command)
        generate_contig_circularity_info("Hifiasm", assembly_basepath + "/quast_outputs/Hifiasm/", assembly_basepath + "/Hifiasm/")
        print("\n\nAssembly using HiFiasm completed. Next, we will perform assembly using IPA...")
    else:
        print("There was an error while performing the assembly. Please look if you can solve it by yourself or send a log file to xxxx@gmail.com")

    os.makedirs(assembly_basepath + "/IPA")
    ipa_command = "ipa local -i " + raw_reads_path + " --run-dir " + assembly_basepath + "/IPA/" + " --tmp-dir "  + assembly_basepath + "/IPA/" + " --njobs 1 --nthreads " + args.threads
    os.system(ipa_command)
    shutil.copyfile(assembly_basepath + "/IPA/19-final/final.p_ctg.fasta", assembly_basepath + "/IPA/final_assembly.fasta")
    if os.path.isfile(assembly_basepath + "/IPA/final_assembly.fasta"):
        quast_command = "quast -o " + assembly_basepath + "/quast_outputs/IPA " + assembly_basepath + "/IPA/final_assembly.fasta"
        os.system(quast_command)
        generate_contig_circularity_info("ipa", assembly_basepath + "/quast_outputs/IPA/", assembly_basepath + "/IPA/")
        print("\n\nAssembly using IPA completed\n\n")
    else:
        print("There was an error while performing the assembly. Please look if you can solve it by yourself or send a log file to xxxx@gmail.com")

    print("Assembly stats for the flye assembly can be found here: " + assembly_basepath + "/quast_outputs/flye/report.pdf")
    print("Assembly stats for the HiFiasm assembly can be found here: " + assembly_basepath + "/quast_outputs/Hifiasm/report.pdf")
    print("Assembly stats for the IPA assembly can be found here: " + assembly_basepath + "/quast_outputs/IPA/report.pdf")
    print("\n\n")
    print("Selecting the best assembly out of all the assemblies....")
    best_assembly = evaluate_assemblies(glob.glob(assembly_basepath + "/quast_outputs/*"), args.genome_size, assembly_basepath + '/')

    final_assembly_path = ""
    if best_assembly == "flye":
        print("The best assembly (among all assemblies) according to us can be found here: " + assembly_basepath + '/' + str(best_assembly) + "/assembly.fasta")
        final_assembly_path = assembly_basepath + '/' + str(best_assembly) + "/assembly.fasta"
    else:
        print("The best assembly (among all assemblies) according to us can be found here: " + assembly_basepath + '/' + str(best_assembly) + "/final_assembly.fasta")
        final_assembly_path = assembly_basepath + '/' + str(best_assembly) + "/final_assembly.fasta"

    return final_assembly_path


def generate_contig_circularity_info(assembler, quast_path, output_path):
    write_file = open(quast_path + 'circularity.tsv', 'x')
    write_file.write("contig_name\tcircular\tLength\n")

    if assembler == "flye":
        with open(output_path + "assembly_info.txt") as assembly_info_file:
            counter = 0
            for line in assembly_info_file.readlines():
                counter += 1
                if counter == 1:
                    continue
                split_array = line.split()
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
