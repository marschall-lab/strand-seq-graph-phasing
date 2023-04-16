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

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
logger = logging.getLogger(__name__)


##### Parsing

parser = argparse.ArgumentParser(description='Split Fasta One Seq per file.')

parser.add_argument('--input', nargs=1, required=True)
parser.add_argument('--output_prefix', nargs=1, required=True)

args = parser.parse_args()
logger.info(args)

###### Main
# 
# input = 'input/references/T2Tv11_hg002Yv2_chm13_hc.fasta'
# out_prefix = 'input/references/T2Tv11_hg002Yv2_chm13_hc.fasta'

with pysam.FastxFile(args.input[0]) as fin:
    n = 0
    for read in fin:
        n += 1
        with open(args.output_prefix[0] + str(n), mode='w') as fout:
            logging.info(read.name)
            fout.write(str(read) + '\n')
            
