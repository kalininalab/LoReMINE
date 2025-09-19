import argparse
import os

def convert_reads(args):
    os.makedirs(args.output + "/raw_reads")
    print("Started converting the raw reads to a .fastq file")
    raw_reads_command =  "samtools fastq "


    if args.pacbio_raw != None:
        raw_reads_command += args.pacbio_raw + " > " + args.output + "/raw_reads/" + args.prefix + '.fastq'
    elif args.pacbio_hifi != None:
        raw_reads_command += args.pacbio_hifi + " > " + args.output + "/raw_reads/" + args.prefix + '.fastq'

    os.system(raw_reads_command)

    print("Raw reads can be found at:" + os.path.abspath(args.output  + "/raw_reads/"))