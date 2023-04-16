import pdb

from .parsing import *


def phase_bubbles(ss_bubble_map_files, bubbles, ss_to_clust, lib_clust_to_haplo, bubble_phase_file, min_h0_frac=0.25, max_h0_frac=0.75):
	'''
	Computes the bubble allele corresponding to the first haplotype (h0) and writes it in the output file
	
	Args:
		ss_bubble_map_files: a list of name of files containing the exact mapping information of SS reads to SNV bubbles
		lib_clust_to_haplo:  A dictionaty that maps a pair (ss_lib, clust) to a haplotype (only for wc ss libs)		
		bubble_phase_file: The output file containing two columns for the bubble name and the bubble allele corresponding to the first haplotype, respectively
	'''

	unitig_haplo_cov = {}
	# for testing
	bubble_unitig_haplo_cov = {}

	if type(ss_bubble_map_files)==str:
		ss_bubble_map_files=[ss_bubble_map_files]

	for input_file in ss_bubble_map_files:
		with open(input_file) as f:
			# skip header
			next(f)
			
			for line in f:
				if line=="":
					break
				
				# skip the header line
				if line.startswith("SSname"):
					continue
					
				sp = line.split()
				
				ss_name, ss_lib, unitig_name, bubble_id, bubble_allele_id, is_reverse = sp[0], sp[1], sp[2], sp[3], sp[4], sp[5]
				
				ss_clust = ss_to_clust[ss_name]

				if (ss_lib, ss_clust) not in lib_clust_to_haplo:
					continue

				haplo = lib_clust_to_haplo[(ss_lib, ss_clust)]
				
				if bubble_id=="None":
					if unitig_name not in unitig_haplo_cov:
						unitig_haplo_cov[unitig_name]=[0,0]
						
					unitig_haplo_cov[unitig_name][haplo] +=1
					continue
					
				if unitig_name not in bubble_unitig_haplo_cov:
					bubble_unitig_haplo_cov[unitig_name]=[0,0]
						
					bubble_unitig_haplo_cov[unitig_name][haplo] +=1

				assert(bubble_allele_id=='0' or bubble_allele_id=='1'), 'error in bubble '+bubble_id+': bubble allele is not binary'
				
				bubble = bubbles[int(bubble_id)]
				bubble_allele = bubble.allele0 if bubble_allele_id=='0' else bubble.allele1
				
				bubble_allele.haplo_coverage[haplo] += 1


	with open(bubble_phase_file, 'w') as out:
		print("#unitigName\thaplotype\tbubbleName\tbubbleAllele", file=out)
		# printing bubble phase info
		for bubble_id, bubble in bubbles.items():
			
			h0_al0 = bubble.allele0.haplo_coverage[0]+bubble.allele1.haplo_coverage[1]
			h0_al1 = bubble.allele0.haplo_coverage[1]+bubble.allele1.haplo_coverage[0]

			if h0_al0+h0_al1==0 or min_h0_frac < h0_al0/(h0_al0+h0_al1) < max_h0_frac:
				continue

			al0_haplo = 'H1' if h0_al0 > h0_al1 else 'H2'
			al1_haplo = 'H1' if h0_al0 < h0_al1 else 'H2'

			print(bubble.allele0.unitig_name + "\t" + al0_haplo + "\t" + str(bubble.id) + "\t0", file=out)
			print(bubble.allele1.unitig_name + "\t" + al1_haplo + "\t" + str(bubble.id) + "\t1", file=out)
			
		# printing unitig phase info
		
		for unitig_name, haplo_cov in unitig_haplo_cov.items():
			
			if sum(haplo_cov)==0 or min_h0_frac < haplo_cov[0]/sum(haplo_cov) < max_h0_frac:
				#print(unitig_name + "\tNone\tNone\tNone", file=out)
				continue
			
			haplo = 'H1' if haplo_cov[0] > haplo_cov[1] else 'H2'
			print(unitig_name + "\t" + haplo + "\tNone\tNone", file=out)

