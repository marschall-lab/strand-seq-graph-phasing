#from output_valid_maps import *
import os
import pdb
import json

import argparse
from module.parsing import *




# unitig_to_bubble_allele = {}

	# with pysam.FastxFile(bubble_fasta_file) as fastafile:
	# 	for read in fastafile:
	# 		sp = read.name.split("_")
	# 		bubble_id, bubble_allele, unitig = sp[1], str(int(sp[3])-1), sp[-1]
	# 		unitig_to_bubble_allele[unitig] = (bubble_id, bubble_allele)
	#
	# return unitig_to_bubble_allele

def is_empty_file(file):
    with open(file) as my_file:
        my_file.seek(0, os.SEEK_END) # go to end of file
        if my_file.tell(): # if current position is truish (i.e != 0)
            my_file.seek(0) # rewind the file for later use
            return False
        else:
            return True

def map_unitig_to_bubble_allele_json(bubble_path):
    with open(bubble_path) as json_data:
        chains = json.load(json_data)
    unitig_to_bubble_allele = {}
    for chain in chains.values():
        # chain_id = chain['chain_id']
        for bubble in chain['bubbles']:
            if bubble['type'] == 'simple':
                bubble_id = str(bubble['id'])
                bubble_tigs = bubble['inside']
                unitig_to_bubble_allele[bubble_tigs[0]] = (bubble_id, '0')
                unitig_to_bubble_allele[bubble_tigs[1]] = (bubble_id, '1')
            else:
                continue
    return unitig_to_bubble_allele


#libs = snakemake.params['libs']
#
#libs_true_order = True
#for i in range(len(libs)):
#	print(libs[i], snakemake.input['ss_reads'][i])
#	print(libs[i], snakemake.input['map'][i])
#
#	if snakemake.input['ss_reads'][i].find(libs[i])==-1 or snakemake.input['map'][i].find(libs[i])==-1:
#		libs_true_order = False
#
#print('libs_true_order:', libs_true_order)

#if snakemake.params['input_type']=='unitig':
#	unitig_to_bubble_allele = map_unitig_to_bubble_allele(snakemake.input['bubbles'])

#clust_pair = snakemake.wildcards["clust_pair"].split('_')

#output_valid_maps(snakemake.input['ss_reads'], snakemake.input['ss_clust_file'], clust_pair, snakemake.input['unitigs'], snakemake.input['map'], snakemake.output[0], snakemake.log[0], snakemake.params['libs'], snakemake.params['input_type'], unitig_to_bubble_allele)


parser = argparse.ArgumentParser(description='Output Valid Maps.')

parser.add_argument('--bubbles', nargs=1, required=True)
parser.add_argument('--map', nargs=1, required=True)
parser.add_argument('--output', nargs=1, required=True)

args = parser.parse_args()

print(args)
# print(args.bubbles[0])
# pdb.set_trace()

# unitig_to_bubble_allele = map_unitig_to_bubble_allele(args.bubbles[0])
if not(is_empty_file(args.bubbles[0])):
    unitig_to_bubble_allele = map_unitig_to_bubble_allele_json(args.bubbles[0])
    output_bwa_fastmap_matches(unitig_to_bubble_allele, args.map[0], args.output[0])

# unitig_to_bubble_allele = map_unitig_to_bubble_allele(snakemake.input['bubbles'])
# output_bwa_fastmap_matches(unitig_to_bubble_allele, snakemake.input['map'], snakemake.output[0])
