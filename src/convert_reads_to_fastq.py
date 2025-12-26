import argparse
import os
import subprocess

def convert_reads(args):

    """Convert the .bam file of the reads into the .fastq format for performing the genome assembly"""

    basepath = os.path.abspath(args.output)
    os.makedirs(basepath + "/raw_reads")
    print("Started converting the raw reads for " + args.prefix + " to a .fastq file")
    raw_reads_command = "samtools fastq "


    if args.pacbio_raw != None:
        raw_reads_command += args.pacbio_raw + " > " + basepath + "/raw_reads/" + args.prefix + '.fastq'
    elif args.pacbio_hifi != None:
        raw_reads_command += args.pacbio_hifi + " > " + basepath + "/raw_reads/" + args.prefix + '.fastq'

    run_command(raw_reads_command, verbosity=args.verbose)

    print("Raw reads can be found at: " + basepath + "/raw_reads")

def run_command(cmd, verbosity=0):
    """
    Run a shell command with controlled verbosity
    """
    if verbosity == 0:
        subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    else:
        subprocess.run(cmd, shell=True, check=True)