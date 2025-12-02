import os
import shutil

def identify_bgcs_single_genome(args):

    """Identify the bgcs in a single genome using antiSMASH"""

    basepath = os.path.abspath(args.output)
    os.makedirs(basepath + "/identified_bgcs")
    filename = os.path.splitext(os.path.basename(args.input_fasta))[0]
    if args.db_path == None:
        command = "antismash --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees --genefinding-tool prodigal --output-dir " + basepath + "/identified_bgcs/" + filename + "/ --output-basename " + filename + " -c " + args.threads + " " + args.input_fasta
    else:
        command = "antismash --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees --genefinding-tool prodigal --output-dir " + basepath + "/identified_bgcs/" + filename + "/ --output-basename " + filename + " -c " + args.threads + " --databases " + args.db_path + " " + args.input_fasta
    print("Identifying BGCs for " + args.input_fasta + ' with ' + args.threads + ' threads\n')
    os.system(command)
    print("BGC identification completed for " + args.input_fasta + '\n')
    print("The identified BGCs can be found here:" + basepath + "/identified_bgcs/" + filename + "/\n\n")

def identify_bgcs_single_genome_all_submodules(args, input_fasta, basepath):

    """Identify the bgcs in a single genome using antiSMASH (while running the whole pipeline together)"""

    os.makedirs(basepath + "/identified_bgcs")
    filename = args.prefix
    if args.db_path == None:
        command = "antismash --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees --genefinding-tool prodigal --output-dir " + basepath + "/identified_bgcs/" + filename + "/ --output-basename " + filename + " -c " + args.threads + " " + input_fasta
    else:
        command = "antismash --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees --genefinding-tool prodigal --output-dir " + basepath + "/identified_bgcs/" + filename + "/ --output-basename " + filename + " -c " + args.threads + " --databases " + args.db_path + " " + input_fasta
    print("Identifying BGCs for " + input_fasta + ' with ' + args.threads + ' threads\n')
    os.system(command)
    print("BGC identification completed for " + input_fasta + '\n')
    print("The identified BGCs can be found here:" + basepath + "/identified_bgcs/" + filename + "/\n\n")
    bgc_output_path = basepath + "/identified_bgcs/" + filename
    return bgc_output_path


def identify_bgcs_multiple_genomes(args):

    """Identify the bgcs in multiple genomes using antiSMASH"""

    basepath = os.path.abspath(args.output)
    if not os.path.exists(basepath + "/identified_bgcs"):
        os.makedirs(basepath + "/identified_bgcs")
    for file in os.listdir(args.input_dir):
        filename = os.path.splitext(os.path.basename(file))[0]
        print("Identifying BGCs for " + file + ' with ' + args.threads + ' threads\n')
        if args.db_path == None:
            command = "antismash --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees --genefinding-tool prodigal --output-dir " + basepath + "/identified_bgcs/" + filename + "/ --output-basename " + filename + " -c " + args.threads + " " + args.input_dir + '/' + file
        else:
            command = "antismash --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees --genefinding-tool prodigal --output-dir " + basepath + "/identified_bgcs/" + filename + "/ --output-basename " + filename + " -c " + args.threads + " --databases " + args.db_path + " " + args.input_dir + '/' + file
        os.system(command)
        print("BGC identification completed for " + file  + '\n')
        print("The identified BGCs can be found here:" + basepath + "/identified_bgcs/" + filename + "/\n\n")

def identify_bgcs_multiple_genomes_all_submodules(args, input_dir, basepath):

    """Identify the bgcs in multiple genomes using antiSMASH (while running the whole pipeline together)"""

    if not os.path.exists(basepath + "/identified_bgcs"):
        os.makedirs(basepath + "/identified_bgcs")
    for file in os.listdir(input_dir):
        if os.path.isdir(input_dir + '/' + file):
            continue
        filename = os.path.splitext(os.path.basename(file))[0]
        print("Identifying BGCs for " + file + ' with ' + args.threads + ' threads\n')
        if args.db_path == None:
            command = "antismash --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees --genefinding-tool prodigal --output-dir " + basepath + "/identified_bgcs/" + filename + "/ --output-basename " + filename + " -c " + args.threads + " " + input_dir + '/' + file
        else:
            command = "antismash --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees --genefinding-tool prodigal --output-dir " + basepath + "/identified_bgcs/" + filename + "/ --output-basename " + filename + " -c " + args.threads + " --databases " + args.db_path + " " + input_dir + '/' + file
        os.system(command)
        print("BGC identification completed for " + file  + '\n')
        print("The identified BGCs can be found here:" + basepath + "/identified_bgcs/" + filename + "/\n\n")

        bgc_output_path = basepath + "/identified_bgcs/" + filename

        if not os.path.exists(basepath + '/identified_bgcs/identified_bgcs_all_strains'):
            os.makedirs(basepath + '/identified_bgcs/identified_bgcs_all_strains')

        for bgc_file in os.listdir(bgc_output_path):
            if "region" in bgc_file and bgc_file.endswith(".gbk"):
                shutil.copy(bgc_output_path + '/' + bgc_file, basepath + '/identified_bgcs/identified_bgcs_all_strains/')

    return basepath + '/identified_bgcs/identified_bgcs_all_strains/'