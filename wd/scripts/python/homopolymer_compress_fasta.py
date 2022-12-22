import argparse
import pdb
import logging
import sys
import itertools

import pysam


##### Logging

Log_Format = "%(levelname)s %(asctime)s - %(message)s"

# logging.basicConfig(#filename = "logfile.log",
#                     stream = sys.stdout,
#                     filemode = "w",
#                     format = Log_Format,
#                     level = logging.ERROR)

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


##### Parsing

parser = argparse.ArgumentParser(description='Homopolymer Compress Fasta.')

parser.add_argument('--input', nargs=1, required=True)
parser.add_argument('--output', nargs=1, required=True)

args = parser.parse_args()
logger.info(args)

#### Functions

def rle(x):
    return [(key, sum(1 for i in group)) for key,group in itertools.groupby(x)]

def homopolymer_compress(x):
    runs = rle(x)
    compressed_read=''
    for base, length in runs:
        assert(base in ['A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n'])
        if base in ['N', 'n']:
            compressed_read += str(base) * length
        else:
            compressed_read += str(base)
    return compressed_read

# Only Used in Testing

def maxrun(x):
    run_lengths = rle(x)
    maxrun = max(run_lengths, key=lambda x:x[1])
    return maxrun


###### Main

with pysam.FastxFile(args.input[0]) as fin,  open(args.output[0], mode='w') as fout:
    for read in fin:
        # print(read.name)
        logger.info(read.name)
        hc = homopolymer_compress(read.sequence)
        read.sequence = hc
        fout.write(str(read) + '\n')
