import argparse
import os
from src.evaluate_best_assembly import *

def perform_assembly_raw_reads_nanopore(args):
    basepath = os.path.abspath(args.output)
    os.makedirs(basepath + "/assembly")
    os.makedirs(basepath + "/assembly/" + args.prefix)
    assembly_basepath = os.path.join(basepath, "assembly", args.prefix)
    os.makedirs(assembly_basepath + "/flye")
    os.makedirs(assembly_basepath + "/quast_outputs")

    if args.alt_param == False:
        flye_command = "flye --nano-raw " + basepath + "/raw_reads/" + args.prefix + ".fastq -o " + assembly_basepath + "/flye/ -t " + args.threads + " -i 3 --scaffold"
        os.system(flye_command)
        if os.path.isfile(assembly_basepath + "/flye/assembly.fasta"):
            print("\n\nAssembly completed. Collecting the assembly stats...\n\n")
            quast_command = "quast -o " + assembly_basepath + "/flye/quast " + assembly_basepath + "/flye/assembly.fasta"
            os.system(quast_command)
            print("Assembly stats can be found here: " + assembly_basepath + "/flye/quast/report.pdf" + "\n\n")
            print("If you are not satisfied with the default assembly, please re-run the assembly with \"alt_param\" flag as \"True\"")
        else:
            print("There was an error while performing the assembly. Please look if you can solve it by yourself or send a log file to xxxx@gmail.com")
    else:
        min_overlap = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
        files_present = {}
        for ele in min_overlap:
            os.makedirs(assembly_basepath + "/flye/" + str(ele))
            flye_command = "flye --nano-raw " + basepath + "/raw_reads/" + args.prefix + ".fastq -o " + assembly_basepath + "/flye/" + str(ele) + "/ -t " + args.threads + " -i 3 --scaffold -m " + str(ele)
            os.system(flye_command)
            if os.path.isfile(assembly_basepath + "/flye/" + str(ele)+  "/assembly.fasta"):
                files_present[ele] = "yes"
                quast_command = "quast -o " + assembly_basepath + "/quast_outputs/" + str(ele) + "/ " + assembly_basepath + "/flye/" + str(ele) + "/assembly.fasta"
                os.system(quast_command)
                generate_contig_circularity_info(assembly_basepath + "/quast_outputs/" + str(ele), assembly_basepath + "/flye/" + str(ele))
            else:
                continue

        print("\n\nGenome assembly completed using all alternate parameters. The assembly stats can be found here:\n\n")

        for key in files_present.keys():
            if files_present[key] == "yes":
                print(os.path.abspath(assembly_basepath + "/quast_outputs/" + str(key) + "/report.pdf"))

        print("Selecting the best assembly out of all the assemblies....")
        best_assembly = evaluate_assemblies(glob.glob(assembly_basepath + "/quast_outputs/*"), args.genome_size, assembly_basepath + "/")
        print("The best assembly (among all assemblies) according to us can be found here: " + assembly_basepath + "/flye/" + str(best_assembly) + "/assembly.fasta")

    if args.alt_param == False:
        final_assembly_path = assembly_basepath + "/flye/assembly.fasta"
    else:
        final_assembly_path = assembly_basepath+ "/flye/" + str(best_assembly) + "/assembly.fasta"

    return final_assembly_path


def perform_assembly_raw_reads_nanopore_batch_run(args, raw_reads_path, genome_size, prefix, output_path):
    basepath = output_path
    if not os.path.exists(basepath + "/assembly"):
        os.makedirs(basepath + "/assembly")
    os.makedirs(basepath + "/assembly/" + prefix)
    assembly_basepath = os.path.join(basepath, "assembly", prefix)
    os.makedirs(assembly_basepath + "/flye")
    os.makedirs(assembly_basepath + "/quast_outputs")

    if args.alt_param == False:
        flye_command = "flye --nano-raw " + raw_reads_path + " -o " + assembly_basepath + "/flye/ -t " + args.threads + " -i 3 --scaffold"
        os.system(flye_command)
        if os.path.isfile(assembly_basepath + "/flye/assembly.fasta"):
            print("\n\nAssembly completed. Collecting the assembly stats...\n\n")
            quast_command = "quast -o " + assembly_basepath + "/flye/quast " + assembly_basepath + "/flye/assembly.fasta"
            os.system(quast_command)
            print("Assembly stats can be found here: " + assembly_basepath + "/flye/quast/report.pdf" + "\n\n")
            print("If you are not satisfied with the default assembly, please re-run the assembly with \"alt_param\" flag as \"True\"")
        else:
            print("There was an error while performing the assembly. Please look if you can solve it by yourself or send a log file to xxxx@gmail.com")
    else:
        min_overlap = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]
        files_present = {}
        for ele in min_overlap:
            os.makedirs(assembly_basepath + "/flye/" + str(ele))
            flye_command = "flye --nano-raw " + raw_reads_path + " -o " + assembly_basepath + "/flye/" + str(ele) + "/ -t " + args.threads + " -i 3 --scaffold -m " + str(ele)
            os.system(flye_command)
            if os.path.isfile(assembly_basepath + "/flye/" + str(ele)+  "/assembly.fasta"):
                files_present[ele] = "yes"
                quast_command = "quast -o " + assembly_basepath + "/quast_outputs/" + str(ele) + "/ " + assembly_basepath + "/flye/" + str(ele) + "/assembly.fasta"
                os.system(quast_command)
                generate_contig_circularity_info(assembly_basepath + "/quast_outputs/" + str(ele), assembly_basepath + "/flye/" + str(ele))
            else:
                continue

        print("\n\nGenome assembly completed using all alternate parameters. The assembly stats can be found here:\n\n")

        for key in files_present.keys():
            if files_present[key] == "yes":
                print(os.path.abspath(assembly_basepath + "/quast_outputs/" + str(key) + "/report.pdf"))

        print("Selecting the best assembly out of all the assemblies....")
        best_assembly = evaluate_assemblies(glob.glob(assembly_basepath + "/quast_outputs/*"), args.genome_size, assembly_basepath + "/")
        print("The best assembly (among all assemblies) according to us can be found here: " + assembly_basepath + "/flye/" + str(best_assembly) + "/assembly.fasta")

    if args.alt_param == False:
        final_assembly_path = assembly_basepath + "/flye/assembly.fasta"
    else:
        final_assembly_path = assembly_basepath+ "/flye/" + str(best_assembly) + "/assembly.fasta"

    return final_assembly_path



def generate_contig_circularity_info(quast_path, output_path):
    write_file = open(quast_path + '/circularity.tsv', 'x')
    write_file.write("contig_name\tcircular\tLength\n")


    with open(output_path + "/assembly_info.txt") as assembly_info_file:
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