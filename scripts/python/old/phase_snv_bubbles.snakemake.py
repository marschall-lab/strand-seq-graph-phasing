import argparse
import pdb
import logging
import sys


from module.phase_snv_bubbles import *
from module.parsing import *

parser = argparse.ArgumentParser(description='Phase SNV Bubbles.')

parser.add_argument('--phased-strand-states', nargs=1, required=True)
parser.add_argument('--ss-clust', nargs=1, required=True)
parser.add_argument('--bubbles', nargs=1, required=True)
parser.add_argument('--map', nargs=1, required=True)
parser.add_argument('--output', nargs=1, required=True)

args = parser.parse_args()

print(args)

print("getting haplo strand states ...")
lib_clust_to_haplo = read_strandphaser_strand_states(args.phased_strand_states[0])

print("getting ss clusters ...")
ss_to_clust = get_ss_clust(args.ss_clust[0])

print('reading bubbles\' fasta file')
bubbles = get_bubbles_from_json(args.bubbles[0]) #, with_km=False, with_unitig_name=True)

print("phasing the bubbles and writing the phase information in the output file ...")
phase_bubbles(args.map[0], bubbles, ss_to_clust, lib_clust_to_haplo, args.output[0])
