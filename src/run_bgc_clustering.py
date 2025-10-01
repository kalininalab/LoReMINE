import os
import shutil
import sqlite3
from sqlite3 import Error

def run_bigslice_clustering(args):

    """Run the bgc clustering using BiG-SLiCE"""

    basepath = os.path.abspath(args.output)
    if not os.path.exists(basepath + "/clustering"):
        os.makedirs(basepath + "/clustering")
    basepath = basepath + "/clustering"
    os.makedirs(basepath + "/bigslice")
    os.makedirs(basepath + "/bigslice/taxonomy")
    os.makedirs(basepath + "/bigslice/dataset")
    os.makedirs(basepath + "/bigslice/dataset/input_bgc")

    # Creating necessary input files and directory structure for bigslice run
    with open(basepath + "/bigslice/taxonomy/dataset_taxonomy.tsv", 'x') as taxonomy_file:
        taxonomy_file.write("Genome folder\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tOrganism\n")
        taxonomy_file.write(basepath + "/bigslice/dataset/input_bgc/" + "\tBacteria\tMyxococcota\tDeltaproteobacteria\tMyxococcales\tMyxococcaceae\tMyxococcus\tMyxococcus xanthus\tMyxococcus xanthus")

    # Creating necessary input files and directory structure for bigslice run
    with open(basepath + "/bigslice/datasets.tsv", 'x') as dataset_file:
        dataset_file.write("# Dataset name\tPath to folder\tPath to taxonomy\tDescription\n")
        dataset_file.write("dataset\t" + basepath + "/bigslice/dataset/\t" + basepath + "/bigslice/taxonomy/dataset_taxonomy.tsv\tNot available")

    for file in os.listdir(args.input_dir):
        if file.endswith(".gbk"):
            shutil.copy(args.input_dir + '/' + file, basepath + "/bigslice/dataset/input_bgc/")


    # Check if MiBiG is activated or not. If yes, then download mibig bgcs and convert them into proper format such that it can be accepted by bigslice
    if args.mibig == True:
        download_mibig(basepath)
        for file in os.listdir(basepath + "/mibig_gbk_4.0"):
            shutil.copy(basepath + "/mibig_gbk_4.0" + '/' + file, basepath + "/bigslice/dataset/input_bgc/")
        os.system(f"sed -i 's/Version      :: False/Version      :: 5.0.0/' {basepath}/bigslice/dataset/input_bgc/BGC*.gbk")
        shutil.rmtree(basepath + "/mibig_gbk_4.0")

    bigslice_command = "bigslice -i " + basepath + "/bigslice/ --threshold " + str(args.bigslice_cutoff) + " -t " + str(args.threads) + " " + basepath + "/bigslice/output"

    os.system(bigslice_command)

    extract_bigslice_output(basepath + "/bigslice/output/")

    print("Clustering using BiGSLiCE completed at threshold " + str(args.bigslice_cutoff) + ". You can find the final clustering output here:" + basepath + "/bigslice/output/output_clusters.tsv")


def run_bigslice_clustering_all_submodules(args, input_dir, basepath):

    """Run the bgc clustering using BiG-SLiCE (while running the whole pipeline together)"""

    if not os.path.exists(basepath + "/clustering"):
        os.makedirs(basepath + "/clustering")
    basepath = basepath + "/clustering"
    os.makedirs(basepath + "/bigslice")
    os.makedirs(basepath + "/bigslice/taxonomy")
    os.makedirs(basepath + "/bigslice/dataset")
    os.makedirs(basepath + "/bigslice/dataset/input_bgc")

    # Creating necessary input files and directory structure for bigslice run
    with open(basepath + "/bigslice/taxonomy/dataset_taxonomy.tsv", 'x') as taxonomy_file:
        taxonomy_file.write("Genome folder\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tOrganism\n")
        taxonomy_file.write(basepath + "/bigslice/dataset/input_bgc/" + "\tBacteria\tMyxococcota\tDeltaproteobacteria\tMyxococcales\tMyxococcaceae\tMyxococcus\tMyxococcus xanthus\tMyxococcus xanthus")

    # Creating necessary input files and directory structure for bigslice run
    with open(basepath + "/bigslice/datasets.tsv", 'x') as dataset_file:
        dataset_file.write("# Dataset name\tPath to folder\tPath to taxonomy\tDescription\n")
        dataset_file.write("dataset\t" + basepath + "/bigslice/dataset/\t" + basepath + "/bigslice/taxonomy/dataset_taxonomy.tsv\tNot available")

    for file in os.listdir(input_dir + '/'):
        if 'region' in file and file.endswith('.gbk'):
            shutil.copy(input_dir + '/' + file, basepath + "/bigslice/dataset/input_bgc/")

    # Check if MiBiG is activated or not. If yes, then download mibig bgcs and convert them into proper format such that it can be accepted by bigslice
    if args.mibig == True:
        download_mibig(basepath)
        for file in os.listdir(basepath + "/mibig_gbk_4.0"):
            shutil.copy(basepath + "/mibig_gbk_4.0" + '/' + file, basepath + "/bigslice/dataset/input_bgc/")
        os.system(f"sed -i 's/Version      :: False/Version      :: 5.0.0/' {basepath}/bigslice/dataset/input_bgc/BGC*.gbk")
        shutil.rmtree(basepath + "/mibig_gbk_4.0")

    bigslice_command = "bigslice -i " + basepath + "/bigslice/ --threshold " + str(args.bigslice_cutoff) + " -t " + str(args.threads) + " " + basepath + "/bigslice/output"

    os.system(bigslice_command)

    extract_bigslice_output(basepath + "/bigslice/output/")

    print("Clustering using BiGSLiCE completed at threshold " + str(args.bigslice_cutoff) + ". You can find the final clustering output here:" + basepath + "/bigslice/output/output_clusters.tsv")


def run_bigscape_clustering(args):

    """Run the bgc clustering using BiG-SCAPE"""

    basepath = os.path.abspath(args.output)
    if not os.path.exists(basepath + "/clustering"):
        os.makedirs(basepath + "/clustering")
    basepath = basepath + "/clustering"
    os.makedirs(basepath + "/bigscape")
    os.makedirs(basepath + "/bigscape/input_bgc")

    for file in os.listdir(args.input_dir):
        if file.endswith(".gbk"):
            shutil.copy(args.input_dir + '/' + file, basepath + "/bigscape/input_bgc/")


    bigscape_command = "bigscape cluster -i " + basepath + "/bigscape/input_bgc/ -o " + basepath + "/bigscape/output/ --include-gbk '*' --pfam-path " + args.pfam_dir + " --gcf-cutoffs " + str(args.bigscape_cutoff) + " --mix --include-singletons -c " + args.threads

    # Check if the MiBiG is activated or not. If yes, then use it for bigscape clustering
    if args.mibig == True:
        bigscape_command += " --mibig-version 4.0"

    os.system(bigscape_command)

    filename = "mix_clustering_c" + str(args.bigscape_cutoff) + ".tsv"

    for root, dirs, files in os.walk(basepath + "/bigscape/output/"):
        if filename in files:
            src_path = os.path.join(root, filename)

    clean_bigscape_output(src_path, basepath + "/bigscape/output/output_clusters.tsv")  # Clean the bigscape output to remove GCFs with only MiBiG BGCs

    print("Clustering using BiGSCAPE completed at threshold " + str(args.bigscape_cutoff) + ". You can find the final clustering output here:" + basepath + "/bigscape/output/output_clusters.tsv")


def run_bigscape_clustering_all_submodules(args, input_dir, basepath):

    """Run the bgc clustering using BiG-SCAPE (while running the whole pipeline together)"""

    if not os.path.exists(basepath + "/clustering"):
        os.makedirs(basepath + "/clustering")
    basepath = basepath + "/clustering"
    os.makedirs(basepath + "/bigscape")
    os.makedirs(basepath + "/bigscape/input_bgc")

    for file in os.listdir(input_dir + '/'):
        if 'region' in file and file.endswith('.gbk'):
            shutil.copy(input_dir + '/' + file, basepath + "/bigscape/input_bgc/")


    bigscape_command = "bigscape cluster -i " + basepath + "/bigscape/input_bgc/ -o " + basepath + "/bigscape/output/ --include-gbk '*' --pfam-path " + args.pfam_dir + " --gcf-cutoffs " + str(args.bigscape_cutoff) + " --mix --include-singletons -c " + args.threads

    # Check if the MiBiG is activated or not. If yes, then use it for bigscape clustering
    if args.mibig == True:
        bigscape_command += " --mibig-version 4.0"

    os.system(bigscape_command)

    filename = "mix_clustering_c" + str(args.bigscape_cutoff) + ".tsv"

    for root, dirs, files in os.walk(basepath + "/bigscape/output/"):
        if filename in files:
            src_path = os.path.join(root, filename)

    clean_bigscape_output(src_path, basepath + "/bigscape/output/output_clusters.tsv")  # Clean the bigscape output to remove GCFs with only MiBiG BGCs

    print("Clustering using BiGSCAPE completed at threshold " + str(args.bigscape_cutoff) + ". You can find the final clustering output here:" + basepath + "/bigscape/output/output_clusters.tsv")


def download_mibig(path):

    """Download and extract the MiBiG bgcs for BiG-SLiCE run"""

    import requests
    import tarfile
    import os

    # URL of the MIBiG 4.0 GenBank archive
    url = "https://dl.secondarymetabolites.org/mibig/mibig_gbk_4.0.tar.gz"


    base_dir = path
    extract_dir = base_dir
    tar_path = os.path.join(base_dir, "mibig_gbk_4.0.tar.gz")

    os.makedirs(base_dir, exist_ok=True)

    # Download the MiBiG BGCs
    print("Downloading MIBiG BGC dataset...")
    response = requests.get(url, stream=True)
    with open(tar_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
    print(f"Downloaded to {tar_path}")

    # Extract the bgcs from zip file
    print(f"Extracting into {extract_dir}/mibig_gbk_4.0/")
    os.makedirs(extract_dir, exist_ok=True)
    with tarfile.open(tar_path, "r:gz") as tar:
        tar.extractall(path=extract_dir)
    print(f"Extraction complete. Files are in {extract_dir}/mibig_gbk_4.0/")

    os.remove(base_dir + "/mibig_gbk_4.0.tar.gz")



def extract_bigslice_output(output_path):

    """Extract the BiG-SLiCE output from a database into a tsv file"""

    database = r"" + output_path + "result/data.db"

    # create a database connection
    conn = create_connection(database)
    cur = conn.cursor()

    cur.execute("SELECT bgc.orig_filename, gcf_membership.gcf_id, gcf_membership.membership_value, bgc.on_contig_edge, bgc.length_nt FROM gcf_membership,bgc WHERE gcf_membership.bgc_id = bgc.id")

    rows = cur.fetchall()
    write_file = open(output_path + 'output_clusters_with_length.tsv', 'x')
    write_file.write('BGC\tGCF_id\tDistance\ton_contig_edge\tlength\n')
    counter = 0

    for row in rows:
        counter += 1
        split_array = str(row).split(',')
        write_file.write(split_array[0].strip()[2:-5] + '\t' + split_array[1].strip() + '\t' + split_array[2].strip() + '\t' + split_array[3].strip() + '\t' + split_array[4].strip()[:-1] + '\n')

    write_file.close()

    mibig_bgc_with_gcf_id = {}
    counter = 0
    other_gcfs = set()

    with open(output_path + 'output_clusters_with_length.tsv') as clustering_file:
        for line in clustering_file.readlines():
            counter += 1
            if counter == 1:
                continue
            split_array = line.split('\t')
            bgc_name = split_array[0].strip()
            gcf_id = split_array[1].strip()

            if "BGC00" in bgc_name:
                mibig_bgc_with_gcf_id[bgc_name] = gcf_id
            else:
                other_gcfs.add(gcf_id)

    write_file = open(output_path + 'output_clusters.tsv', 'x')
    write_file.write('BGC\tGCF_id\tDistance\ton_contig_edge\tlength\n')


    # Clean the BiG-SLiCE output to remove GCFs which has only MiBiG bgcs
    counter = 0
    with open(output_path + 'output_clusters_with_length.tsv') as clustering_file:
        for line in clustering_file.readlines():
            counter += 1
            if counter == 1:
                continue
            split_array = line.split('\t')

            if split_array[0].strip() in mibig_bgc_with_gcf_id.keys():
                mibig_gcf_id = mibig_bgc_with_gcf_id[split_array[0].strip()]
                if mibig_gcf_id in other_gcfs:
                    write_file.write(line)
            else:
                write_file.write(line)

    os.remove(output_path + 'output_clusters_with_length.tsv')


def create_connection(db_file):
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except Error as e:
        print(e)

    return conn

def clean_bigscape_output(input_file, output_path):

    """Clean the BiG-SCAPE output to remove GCFs which has only MiBiG bgcs"""

    mibig_bgc_with_gcf_id = {}
    counter = 0
    other_gcfs = set()

    with open(input_file) as clustering_file:
        for line in clustering_file.readlines():
            counter += 1
            if counter == 1:
                continue
            split_array = line.split('\t')
            bgc_name = split_array[1].strip()
            gcf_id = split_array[5].strip()

            if "BGC00" in bgc_name:
                mibig_bgc_with_gcf_id[bgc_name] = gcf_id
            else:
                other_gcfs.add(gcf_id)

    write_file = open(output_path, 'x')
    write_file.write('Record\tGBK\tRecord_Type\tRecord_Number\tCC\tFamily\n')

    counter = 0
    with open(input_file) as clustering_file:
        for line in clustering_file.readlines():
            counter += 1
            if counter == 1:
                continue
            split_array = line.split('\t')

            if split_array[1].strip() in mibig_bgc_with_gcf_id.keys():
                mibig_gcf_id = mibig_bgc_with_gcf_id[split_array[1].strip()]
                if mibig_gcf_id in other_gcfs:
                    write_file.write(line)
            else:
                write_file.write(line)