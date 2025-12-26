import os
import subprocess

def taxonomy_single_genome(args):

    """Identify the taxonomy of a single genome"""

    basepath = os.path.abspath(args.output)
    if not os.path.exists(basepath + "/taxonomy"):
        os.makedirs(basepath + "/taxonomy")
    filename = os.path.splitext(os.path.basename(args.input_fasta))[0]
    command = "dfast_qc -i " + args.input_fasta + " -o " + basepath + "/taxonomy/" + filename + "/ -n " + args.threads + " --enable_gtdb --disable_cc --force > /dev/null 2>&1"
    print("Running taxonomy for " + args.input_fasta + '\n')
    run_command(command, verbosity=args.verbose)
    print("Taxonomy identification completed for " + args.input_fasta + '\n')
    extract_taxonomy_info(basepath + "/taxonomy/" + filename + '/')   # Collect the taxonomy outputs
    print("For the final identified taxonomy output, refer to this file:" + basepath + "/taxonomy/" + filename + "/identified_taxonomy.txt\n\n")

def taxonomy_single_genome_all_submodules(args, input_fasta, basepath):

    """Identify the taxonomy of a single genome (while running the whole pipeline together)"""

    if not os.path.exists(basepath + "/taxonomy"):
        os.makedirs(basepath + "/taxonomy")
    filename = args.prefix
    command = "dfast_qc -i " + input_fasta + " -o " + basepath + "/taxonomy/" + filename + "/ -n " + args.threads + " --enable_gtdb --disable_cc --force > /dev/null 2>&1"
    print("Running taxonomy for " + input_fasta + '\n')
    run_command(command, verbosity=args.verbose)
    print("Taxonomy identification completed for " + input_fasta + '\n')
    extract_taxonomy_info(basepath + "/taxonomy/" + filename + '/')   # Collect the taxonomy outputs
    print("For the final identified taxonomy output, refer to this file:" + basepath + "/taxonomy/" + filename + "/identified_taxonomy.txt\n\n")


def taxonomy_multiple_genomes(args):

    """Identify the taxonomy of multiple genomes"""

    basepath = os.path.abspath(args.output)
    if not os.path.exists(basepath + "/taxonomy"):
        os.makedirs(basepath + "/taxonomy")
    for file in os.listdir(args.input_dir):
        filename = os.path.splitext(os.path.basename(file))[0]
        print("Running taxonomy for " + file + '\n')
        command = "dfast_qc -i " + args.input_dir + '/' + file + " -o " + basepath + "/taxonomy/" + filename + "/ -n " + args.threads + " --enable_gtdb --disable_cc --force > /dev/null 2>&1"
        run_command(command, verbosity=args.verbose)
        print("Taxonomy identification completed for " + file  + '\n')
        extract_taxonomy_info(basepath + "/taxonomy/" + filename + '/')   # Collect the taxonomy outputs
        print("For the final identified taxonomy output, refer to this file:" + basepath + "/taxonomy/" + filename + "/identified_taxonomy.txt\n\n")


def taxonomy_multiple_genomes_all_submodules(args, input_dir, basepath):

    """Identify the taxonomy of multiple genomes (while running the whole pipeline together)"""

    if not os.path.exists(basepath + "/taxonomy"):
        os.makedirs(basepath + "/taxonomy")
    for file in os.listdir(input_dir):
        if os.path.isdir(input_dir + '/' + file):
            continue
        filename = os.path.splitext(os.path.basename(file))[0]
        print("Running taxonomy for " + file + '\n')
        command = "dfast_qc -i " + input_dir + '/' + file + " -o " + basepath + "/taxonomy/" + filename + "/ -n " + args.threads + " --enable_gtdb --disable_cc --force > /dev/null 2>&1"
        run_command(command, verbosity=args.verbose)
        print("Taxonomy identification completed for " + file  + '\n')
        extract_taxonomy_info(basepath + "/taxonomy/" + filename + '/')   # Collect the taxonomy outputs
        print("For the final identified taxonomy output, refer to this file:" + basepath + "/taxonomy/" + filename + "/identified_taxonomy.txt\n\n")


def extract_taxonomy_info(output_path):

    """Collect both NCBI & GTDB taxonomy outputs and summarize in one file """

    write_file = open(output_path + 'identified_taxonomy.txt', 'x')

    #################################### Calculate the taxonomy output from NCBI ##############################

    write_file.write("\n\n##################### NCBI database taxonomy result #####################\n\n")

    status_ncbi = {}
    accession_ncbi = {}
    taxid_ncbi = {}

    status_array_ncbi = []

    with open(output_path + 'tc_result.tsv') as ncbi_taxonomy_file:
        counter = 0
        for line in ncbi_taxonomy_file.readlines():
            counter += 1
            if counter == 1:
                continue
            split_array = line.split('\t')
            status_ncbi[split_array[0].strip()] = split_array[-1].strip()
            status_array_ncbi.append(split_array[-1].strip())
            accession_ncbi[split_array[0].strip()] = split_array[2].strip()
            taxid_ncbi[split_array[0].strip()] = split_array[3].strip()

        if len(status_array_ncbi) == 0:
            write_file.write("Sorry, we were not able to identify the taxonomy using NCBI database\n\n")

        elif "conclusive" in status_array_ncbi:
            write_file.write("The identified taxonomy using NCBI is:\n\n")
            keys_with_conclusive_ncbi = [k for k, v in status_ncbi.items() if v == "conclusive"]
            key_counter_ncbi = 1
            for key_ncbi in keys_with_conclusive_ncbi:
                write_file.write("\t" + str(key_counter_ncbi) + ". Organism: " + key_ncbi + ", NCBI Accession: " + accession_ncbi[key_ncbi] + ", NCBI Taxonomy ID: " + taxid_ncbi[key_ncbi] + "\n")
                key_counter_ncbi += 1

        elif "inconclusive" in status_array_ncbi:
            write_file.write("Sorry, we weren't able to conclusively identify the taxonomy using NCBI. However, the top candidates (in decreasing order of priority) are listed below:\n\n")
            keys_with_inconclusive_ncbi = [k for k, v in status_ncbi.items() if v == "inconclusive"]
            key_counter_ncbi = 1
            for key_ncbi in keys_with_inconclusive_ncbi:
                write_file.write("\t" + str(key_counter_ncbi) + ". Organism: " + key_ncbi + ", NCBI Accession: " + accession_ncbi[key_ncbi] + ", NCBI Taxonomy ID: " + taxid_ncbi[key_ncbi] + "\n")
                key_counter_ncbi += 1

        else:
            write_file.write("Sorry, we weren't able to identify the exact taxonomy using NCBI as all the hits are below the defined species threshold. However, the identified candidates (in decreasing order of priority) below the threshold are listed below:\n\n")
            keys_with_below_threshold_ncbi = [k for k, v in status_ncbi.items() if v == "below_threshold"]
            key_counter_ncbi = 1
            for key_ncbi in keys_with_below_threshold_ncbi:
                write_file.write("\t" + str(key_counter_ncbi) + ". Organism: " + key_ncbi + ", NCBI Accession: " + accession_ncbi[key_ncbi] + ", NCBI Taxonomy ID: " + taxid_ncbi[key_ncbi] + "\n")
                key_counter_ncbi += 1

    #################################### Calculate the taxonomy output from GTDB ##############################

    write_file.write("\n\n\n\n##################### GTDB database taxonomy result #####################\n\n")

    status_gtdb = {}
    accession_gtdb = {}
    high_taxonomy_gtdb = {}

    status_array_gtdb = []

    with open(output_path + 'result_gtdb.tsv') as gtdb_taxonomy_file:
        counter = 0
        for line in gtdb_taxonomy_file.readlines():
            counter += 1
            if counter == 1:
                continue
            split_array = line.split('\t')
            status_gtdb[split_array[1].strip()] = split_array[-1].strip()
            status_array_gtdb.append(split_array[-1].strip())
            accession_gtdb[split_array[1].strip()] = split_array[0].strip()
            high_taxonomy_gtdb[split_array[1].strip()] = split_array[5].strip()

        if len(status_array_gtdb) == 0:
            write_file.write("Sorry, we were not able to identify the taxonomy using GTDB database\n\n")

        elif "conclusive" in status_array_gtdb:
            write_file.write("The identified taxonomy using GTDB database is:\n\n")
            keys_with_conclusive_gtdb = [k for k, v in status_gtdb.items() if v == "conclusive"]
            key_counter_gtdb = 1
            for key_gtdb in keys_with_conclusive_gtdb:
                write_file.write("\t" + str(key_counter_gtdb) + ". " + high_taxonomy_gtdb[key_gtdb] + ";" + key_gtdb + " (NCBI accession: " + accession_gtdb[key_gtdb] + ")\n")
                key_counter_gtdb += 1

        elif "inconclusive" in status_array_gtdb:
            write_file.write("Sorry, we weren't able to conclusively identify the taxonomy using GTDB database. However, the top candidates (in decreasing order of priority) are listed below:\n\n")
            keys_with_inconclusive_gtdb = [k for k, v in status_gtdb.items() if v == "inconclusive"]
            key_counter_gtdb = 1
            for key_gtdb in keys_with_inconclusive_gtdb:
                write_file.write("\t" + str(key_counter_gtdb) + ". " + high_taxonomy_gtdb[key_gtdb] + ";" + key_gtdb + " (NCBI accession: " + accession_gtdb[key_gtdb] + ")\n")
                key_counter_gtdb += 1

        else:
            write_file.write("Sorry, we weren't able to identify the exact taxonomy using GTDB database as all the hits are below the defined species threshold. However, the identified candidates (in decreasing order of priority) below the threshold are listed below:\n\n")
            keys_with_below_threshold_gtdb = [k for k, v in status_gtdb.items() if v == "below_threshold"]
            key_counter_gtdb = 1
            for key_gtdb in keys_with_below_threshold_gtdb:
                write_file.write("\t" + str(key_counter_gtdb) + ". " + high_taxonomy_gtdb[key_gtdb] + ";" + key_gtdb + " (NCBI accession: " + accession_gtdb[key_gtdb] + ")\n")
                key_counter_gtdb += 1

def run_command(cmd, verbosity=0):
    """
    Run a shell command with controlled verbosity
    """
    if verbosity == 0:
        subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    else:
        subprocess.run(cmd, shell=True, check=True)