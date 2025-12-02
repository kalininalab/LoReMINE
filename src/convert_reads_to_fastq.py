import argparse
import os

def convert_reads(args):

    """Convert the .bam file of the reads into the .fastq format for performing the genome assembly"""

    basepath = os.path.abspath(args.output)
    os.makedirs(basepath + "/raw_reads")
    print("Started converting the raw reads to a .fastq file")
    raw_reads_command =  "samtools fastq "


    if args.pacbio_raw != None:
        raw_reads_command += args.pacbio_raw + " > " + basepath + "/raw_reads/" + args.prefix + '.fastq'
    elif args.pacbio_hifi != None:
        raw_reads_command += args.pacbio_hifi + " > " + basepath + "/raw_reads/" + args.prefix + '.fastq'

    os.system(raw_reads_command)

    print("Raw reads can be found at: " + basepath + "/raw_reads")