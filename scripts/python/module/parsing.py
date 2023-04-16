import time
import pdb
import gzip
import json

import pysam

from .bubble_long_read_alignment import *

global valid_chroms
valid_chroms = ['chr' + str(i) for i in range(1,23)] + ['chrX'] + ['chrY']


def get_seq_len(fasta_files):
	if type(fasta_files)==str:
		fasta_files = [fasta_files]

	seq_to_len = {}

	for fasta_file in fasta_files:
		with pysam.FastxFile(fasta_file) as fastafile:
			for read in fastafile:
				seq_to_len[read.name] = len(read.sequence)

	return (seq_to_len)


def map_unitig_to_bubble_allele(bubble_fasta_file):
	unitig_to_bubble_allele = {}

	with pysam.FastxFile(bubble_fasta_file) as fastafile:
		for read in fastafile:
			sp = read.name.split("_")
			bubble_id, bubble_allele, unitig = sp[1], str(int(sp[3])-1), sp[-1]
			unitig_to_bubble_allele[unitig] = (bubble_id, bubble_allele)

	return unitig_to_bubble_allele


def output_bwa_fastmap_matches(unitig_to_bubble_allele, fastmap_file, output_file):
	with open(output_file, 'w') as out:
		with open(fastmap_file) as f:
			prev_record = ""
			unique_match = True
			print("SSname\tSSlib\tunitig_name\tbubbleName\tbubbleAllele\tisReverseMapped", file=out)

			for line in f:
				sp = line.split()
				#TODO remove //SQ
				if len(sp)==0 or sp[0] not in ['SQ', 'EM',"//SQ"]:
					continue

				if sp[0]=="SQ" or sp[0]=="//SQ": # SS read
					if unique_match and prev_record=="EM":
						# print the previous match
						bubble_id, bubble_allele = "None", "None"
						# pdb.set_trace()
						if unitig_name in unitig_to_bubble_allele:
							bubble_id, bubble_allele = unitig_to_bubble_allele[unitig_name]
							#pdb.set_trace()

						print(ss_name+"\t"+ss_lib+"\t"+unitig_name+"\t"+bubble_id+"\t"+bubble_allele+"\t"+is_reverse, file=out)

					# Split only on first underscore
					ss_name, *ss_lib = sp[1].split('_')
					ss_lib = '_'.join(ss_lib)

					unique_match = True
					num_EM_lines = 0
					prev_record="SQ"

				else: # unitigs line
					num_EM_lines += 1

					if num_EM_lines > 1 or len(sp) > 5 or sp[-1]=="*":
						unique_match = False
						continue

					unitig_name, map_info = sp[-1].split(":")
					strand = map_info[0]
					is_reverse = "True" if strand=="-" else "False"
					prev_record="EM"


def read_strand_states(strand_state_files):
	'''
	Reads phased strand states from the input file

	Args:
		strand_state_file: the file containing phased strand states
	'''

	lib_clust_to_haplo = {}

	for f in strand_state_files:
		file_name = f.split("/")[-1]
		lib_name=file_name.split("_haplo_strand_states")[0]
		with open(f) as states:
			for line in states:
				# process the line if it is not empty nor the header line
				if line!="" and line[0]!="V":
					sp = line.split()
					lib_clust_to_haplo[(lib_name, sp[0])]=int(sp[1])

	return lib_clust_to_haplo


def read_strandphaser_strand_states(strand_states_file):
	'''
	Reads phased strand states from the input file

	Args:
		strand_states_file: the file containing phased strand states
	'''

	lib_clust_to_haplo = {}

	with open(strand_states_file) as states:
		#skip the header line
		next(states)
		for line in states:
			if line!="":
				sp = line.split()
				lib_clust_to_haplo[(sp[0], sp[1])]=int(sp[2])
	return lib_clust_to_haplo

def get_ss_clust(ss_clust_file):
	print('file', ss_clust_file)
	ss_to_clust = {}
	with open(ss_clust_file) as f:
		next(f) # skip the header line
		for line in f:
			if line == "":
				continue
			sp = line.split()
			ss_name, ss_clust = sp[0], sp[1]
			ss_name = ss_name.split('_')[0]
			ss_to_clust[ss_name] = ss_clust
	return ss_to_clust

def get_bubbles(bubble_fasta_file, with_km=True, with_unitig_name=False):

	'''
	Reads bubbles fasta file, creates bubble and bubbleAllele objects, and returns the list of bubbles

	Parameters (str):
		bubble_fasta_file: The path of bubbles fasta file

	Returns (dict(int): {bubble_id -> bubble}):
		A dictionary that maps each bubble_id to its bubble object

	'''

	start_time = time.time()
	print('getting bubbles from', bubble_fasta_file)

	bubbles = {}

	with pysam.FastxFile(bubble_fasta_file) as fastafile:
		for read in fastafile:
			bubble_name = read.name
			seq = read.sequence

			bubble_name_sp = bubble_name.split('_')
			bubble_id, bubble_allele_id = int(bubble_name_sp[1]), int(bubble_name_sp[3])-1

			if bubble_id in bubbles:
				bubble = bubbles[bubble_id]

			else:
				bubble = Bubble(bubble_id)
				bubbles[bubble_id] = bubble


			unitig_name = bubble_name_sp[-1] if with_unitig_name else None

			bubble_allele = BubbleAllele(bubble_allele_id, bubble_name, bubble, seq, unitig_name)
			bubble.add_allele(bubble_allele)

			if with_km:
				bubble_allele.km = float(bubble_name_sp[7])

	print('elapsed time =', time.time()-start_time)

	return bubbles



def get_bubbles_from_json(bubble_json_file):

  print('getting bubbles from', bubble_json_file)

  with open(bubble_json_file) as json_data:
    chains = json.load(json_data)

  bubbles = {}
  for chain in chains.values():
    for bubble in chain['bubbles']:
      if bubble['type'] != 'simple':
        continue
      else:
        bubble_id = int(bubble['id'])
        bubble_tigs = bubble['inside']

        bubble = Bubble(bubble_id)
        for bubble_allele_id in [0, 1]:
          unitig_name = bubble_tigs[bubble_allele_id]
          bubble_name = 'bubble_' + str(bubble_id) + '_allele_' + str(bubble_allele_id) + '_unitig_' + unitig_name
          # No seq info in BubbleGun Json
          bubble_allele = BubbleAllele(bubble_allele_id, bubble_name, bubble, None, unitig_name)
          bubble.add_allele(bubble_allele)

        bubbles[bubble_id] = bubble

  return bubbles



def get_bubbles_from_bam(bubble_haplotagged_bam_file):

	'''
	Reads bubbles happlotagged bam file, creates bubble and bubbleAllele objects, and returns the list of bubbles

	Parameters (str):
		bubble_happlotagged bam_file: The path of bubbles happlotagged bam file

	Returns (dict(int): {bubble_id -> bubble}):
		A dictionary that maps each bubble_id to its bubble object

	'''

	start_time = time.time()
	print('getting bubbles from', bubble_haplotagged_bam_file)

	bubbles = {}

	samfile = pysam.AlignmentFile(bubble_haplotagged_bam_file, 'rb')

	for read in samfile.fetch(until_eof=True): # until_eof=True allows to read also unmapped reads in sam file
		bubble_name = read.query_name
		bubble_name_sp = bubble_name.split('_')
		bubble_id, bubble_allele_id = int(bubble_name_sp[1]), int(bubble_name_sp[3])-1

		if bubble_id in bubbles:
			bubble = bubbles[bubble_id]

		else:
			bubble = Bubble(bubble_id)
			bubbles[bubble_id] = bubble

		bubble_allele = BubbleAllele(bubble_allele_id, bubble_name, bubble)
		bubble.add_allele(bubble_allele)

		if read.is_unmapped:
			bubble.actual_type = 'unmapped'
			continue

		chrom = read.reference_name

		if not chrom in valid_chroms:
			# the chromosome name is not valid
			bubble.actual_type = 'invalid_chrom'
			continue

		bubble.actual_chrom = chrom

		if not read.has_tag("HP"):
			# the read is not haplotagged
			bubble.actual_type = 'untagged'
			continue

		bubble.actual_type = 'tagged' + str(read.get_tag("HP")-1)
		bubble_allele.actual_haplo = read.get_tag("HP")-1

	print('num bubbles =', len(bubbles))

	print('elapsed time =', time.time()-start_time)

	return bubbles


def get_clust_to_chrom(clust_to_chrom_file):

	'''
	'''

	start_time = time.time()
	print('getting cluster to chrom_dir mapping from', clust_to_chrom_file)
	clust_to_chrom = {}
	with open(clust_to_chrom_file) as f:
		for line in f:
			if not line.startswith('chr'):
				# it is a hedear line
				continue

			chrom, clust_forward, clust_backward = line.split()

			if clust_forward not in clust_to_chrom:
				# remove the direction from chrom name
				clust_to_chrom[clust_forward] = chrom.split('_')[0]

	print('elapsed time =', time.time()-start_time)

	return clust_to_chrom


def add_bubble_clust(bubble_clust_file, bubbles):

	'''
	Adds SaaRclust (chrom+dir) clusters to bubble objects

	Parameters (str):
		bubble_clust_file (str): Path of bubbles clusters file

		bubbles (dict(int): {bubble_id -> bubble}):
		A dictionary that maps each bubble_id to its bubble object
	'''

	start_time = time.time()
	print('adding SaaRclust chromosome clusters to bubbles from the file', bubble_clust_file)

	with open(bubble_clust_file) as f:
		# skip the header line
		next(f)

		for line in f:
			sp = line.split()

			bubble_id, bubble_clust = int(sp[1].split('_')[1]), sp[0]

			assert (bubble_id in bubbles), 'bubble ' + str(bubble_id) + ' is not present in the bubbles'

			bubble = bubbles[bubble_id]
			bubble.clust = bubble_clust

	print('elapsed time =', time.time()-start_time)


def add_bubble_allele_pred_haplo(bubble_phase_file, bubbles):

	'''
	Adds Haploclust predicted haploypes to bubble allele objects

	Parameters (str):
		bubble_phase_file (str): Path of bubbles phase file

		bubbles (dict(int): {bubble_id -> bubble}):
		A dictionary that maps each bubble_id to its bubble object

	'''

	start_time = time.time()
	print('adding haplotype clusers to bubble alleles from the file', bubble_phase_file)

	#for phase_file in bubble_phase_file:
	with open(bubble_phase_file) as f:
		# skip the header line
		next(f)

		for line in f:
			sp = line.split()

			if sp[1]=="none":
				# the bubble is not phased
				continue

			bubble_id, al0_haplo = int(sp[0]), sp[1]
			al0_haplo = 0 if al0_haplo=="H1" else 1

			if bubble_id not in bubbles:
				# the bubble is not in this cluster being]
				continue

			bubble = bubbles[bubble_id]
			bubble.allele0.pred_haplo = al0_haplo
			bubble.allele1.pred_haplo = 1-al0_haplo

	print('elapsed time =', time.time()-start_time)


def get_long_reads(long_reads_fasta_file):

	'''
	Reads long reads from fasta files and creates LongRead objects

	Parameters (str):
		long_reads_fasta_files: List of the paths of long reads fasta files

	Returns (dict(int): {long_read_name -> long_read}):
		A dictionary that maps each long_read_name to its LongRead object
	'''

	start_time = time.time()
	print('getting long reads')

	long_reads = {}

	#for f in long_reads_fasta_files:
	print('getting long reads from', long_reads_fasta_file, '...')

	with pysam.FastxFile(long_reads_fasta_file) as fastafile:
		for read in fastafile:
			read_name = read.name
			seq = read.sequence

#			read_name_sp = read_name.split('/ccs')
#			read_name = read_name_sp[0]+'/ccs'

			long_read = LongRead(read_name, seq)
			long_reads[read_name] = long_read

	print('elapsed time =', time.time()-start_time)

	return long_reads


def get_long_reads_from_bam(long_reads_haplotagged_bam_files):

	'''
	Reads long reads from haplotagged bam files and creates LongRead objects

	Parameters (str):
		long_reads_fasta_files: List of the paths of long reads haplotagged bam files

	Returns (dict(int): {long_read_name -> long_read}):
		A dictionary that maps each long_read_name to its LongRead object
	'''

	start_time = time.time()
	print('getting long reads from', long_reads_haplotagged_bam_files)

	if type(long_reads_haplotagged_bam_files)!=list:
		long_reads_haplotagged_bam_files = list(long_reads_haplotagged_bam_files)

	long_reads = {}

	for alignmentfile in long_reads_haplotagged_bam_files:
		print('processing', alignmentfile, '...')
		samfile = pysam.AlignmentFile(alignmentfile, 'rb')
		for read in samfile.fetch(until_eof=True):
			read_name = read.query_name

			long_read = LongRead(read_name)
			long_reads[read_name] = long_read

			if read.is_unmapped:
				long_read.actual_type = 'unmapped'
				continue

			chrom = read.reference_name
			long_read.actual_chrom = chrom

			if not chrom in valid_chroms:
				# the chromosome name is not valid
				long_read.actual_type = 'invalid_chrom'
				continue

			if not read.has_tag("HP"):
				# the read is not haplotagged
				long_read.actual_type = 'untagged'
				continue

			haplotype = read.get_tag("HP")-1
			long_read.actual_haplo = haplotype
			long_read.actual_type = 'tagged' + str(haplotype)

	print('elapsed time =', time.time()-start_time)

	return long_reads


def add_long_reads_clust(long_reads_clust_files, long_reads):

	'''
	Adds SaaRclust (chrom+dir) clusters to LongRead objects

	Parameters (str):
		long_reads_clust_files (str): Path of long_reads clusters file

		long_reads (dict(str): {long_read_name -> long_read}):
		A dictionary that maps each read_name to its long_read object

	'''

	start_time = time.time()
	for clust_file in long_reads_clust_files:
		print('reading clusters from', clust_file)
		with open(clust_file) as f:
			for line in f:
				sp = line.split()
				read_name, clust = sp[0], sp[1]
				read_name_sp = read_name.split('/ccs')
				read_name = read_name_sp[0]+'/ccs'

				assert (read_name in long_reads), read_name + ' is not present in the long_reads'
				# to be removed
				#if read_name not in long_reads:
				#	continue

				long_reads[read_name].clust=clust


		print('elapsed time =', time.time()-start_time)


def add_long_reads_pred_haplotype(long_reads_phase_file_list, long_reads):

	'''
	Adds Haploclust predicted haploypes to bubble allele objects

	Parameters (str):
		long_reads_phase_file_list (str): A list of paths of bubbles haplotagged bam files

		long_reads (dict(str): {long_read name -> long_read}):
		A dictionary that maps each read name to its LongRead object

	'''

	start_time = time.time()

	for long_reads_phase_file in long_reads_phase_file_list:
		print('adding pred haplo to long reads from', long_reads_phase_file)
		with open(long_reads_phase_file) as f:
			# skip the header line
			next(f)
			for line in f:
				sp = line.split()
				read_name, haplo = sp[0], sp[-1]
#				read_name_sp = read_name.split('/ccs')
#				read_name = read_name_sp[0]+'/ccs'

				assert (read_name in long_reads), 'read ' + read_name + ' should be present in long_reads'
				# to be removed
				#if read_name not in long_reads:
				#	continue

				long_read = long_reads[read_name]

				if haplo!="none":
					long_read.pred_haplo = 0 if haplo=="H1" else 1

	print('elapsed time =', time.time()-start_time)




def set_alignments_from_minimap_file(minimap_file, bubbles, long_reads):

	'''
	Given a list of minimap alignment files, and bubbles and long_reads, creates alignment objects

	Parameters:
		minimap_files_list: A list of paths to minimap alignment files
		bubbles			 : A dictionary {bubble_id -> bubbles}
		long_reads		 : A dictionary {long_read_name -> long_read}
	'''

	start_time = time.time()

	#for minimap_file in minimap_files_list:
	print('reading alignments from file', minimap_file)
	with gzip.open(minimap_file) as minimap:

		for line in minimap:
			line = line.decode("utf-8")
			sp = line.split()

			bubble_name, bubble_len, bubble_start, bubble_end, strand, \
			read_name, long_read_len, long_read_start, long_read_end = \
			sp[0], int(sp[1]), int(sp[2]), int(sp[3])-1, sp[4], \
			sp[5], int(sp[6]), int(sp[7]), int(sp[8])-1

			bubble_name_sp = bubble_name.split('_')

			bubble_id, bubble_allele_id = int(bubble_name_sp[1]), int(bubble_name_sp[3])-1

			assert(bubble_allele_id == 0 or bubble_allele_id == 1), 'bubble ' + str(bubble_id) + ': allele should be 0 or 1'
			# remove the beginning part of the cigar string
			cigar = sp[-1] # assuming that always the last tag is cigar
			cigar = cigar.split('cg:Z:')[1]

			read_name_sp = read_name.split('/ccs')
			read_name = read_name_sp[0]+'/ccs'

			if bubble_id not in bubbles:
				# the bubble is not from the same clust_pair (chromosome) as the long read
				continue

			assert (read_name in long_reads), 'long read ' + read_name + ' is not present in long reads'

			bubble = bubbles[bubble_id]

			if bubble.allele0==None or bubble.allele1==None:
				# Not both alleles exist or not both clustered in the same chromosome
				del bubbles[bubble_id]
				continue

			long_read = long_reads[read_name]

			bubble_allele = bubble.allele0 if bubble_allele_id == 0 else bubble.allele1

			bubble_allele_seq = bubble_allele.seq

			assert (bubble_allele != None), 'bubble ' + str(bubble_id) + ' allele ' + str(bubble_al) + ' is None'

			aln = Alignment(long_read=long_read, bubble_allele=bubble_allele, long_read_start=long_read_start, long_read_end=long_read_end, \
									bubble_start=bubble_start, bubble_end=bubble_end, strand=strand, cigar=cigar)

	print('elapsed time =', time.time()-start_time)


def set_alignments_from_kmers_file(kmers_files_list, bubbles, phased_long_reads):

	'''
	Given a list of bubble/long_read alignment kmer files, and bubbles and phased_long_reads, creates alignment objects

	Parameters:
		kmers_files_list: A list of paths to bubble/long_read alignment kmer files
		bubbles			 : A dictionary {bubble_id -> bubbles}
		phased_long_reads		 : A dictionary {long_read_name -> long_read}
	'''

	start_time = time.time()
	print('getting alignments')

	for kmers_files in kmers_files_list:
		print('reading kmers from file', kmers_files)
		with open(kmers_files) as f:
			# skip the header line
			next(f)
			for line in f:
				sp = line.split()
				if sp[0]=='none':
					continue # the long read is not mapped to any bubble

				bubble_id, bubble_al, read_name, bubble_kmer, long_read_kmer, edit_dist = \
				int(sp[0]), int(sp[1]), sp[2], sp[3], sp[4], int(sp[5])

				read_name_sp = read_name.split('/ccs')
				read_name = read_name_sp[0]+'/ccs'

				assert (bubble_id in bubbles), 'bubble ' + str(bubble_id) + ' is not present in the bubbles'
				assert (read_name in phased_long_reads), 'long read ' + read_name + ' is not present in long reads'

				bubble = bubbles[bubble_id]
				long_read = phased_long_reads[read_name]

				bubble_allele = bubble.allele0 if bubble_al==0 else bubble.allele1

				assert (bubble_allele != None), 'bubble ' + str(bubble_id) + ' allele ' + str(bubble_al) + ' is None'

				aln = Alignment(long_read=long_read, bubble_allele=bubble_allele, bubble_kmer=bubble_kmer, long_read_kmer=long_read_kmer, edit_dist=edit_dist)

	print('elapsed time =', time.time()-start_time)
