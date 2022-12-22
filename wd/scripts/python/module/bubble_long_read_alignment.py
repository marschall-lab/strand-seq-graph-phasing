import pdb

from whatshap.align import edit_distance

from .string_functions import *


class Bubble:

	def __init__(self, id=None):
		self.id = id
		self.actual_chrom = None
		self.clust = None
		self.allele0, self.allele1 = None, None
		self.actual_type = None # in ["unmapped", "invalid_chrom", "untagged"]
		self.pred_type = None # in ["true_haplo", "false_haplo", "haploclust_false_pos", "not_chrom_clust", "garbage_clust", "not_haplo_clust"]
		self.het_positions = []
		self.num_al0_h0_reads = 0
		self.num_al0_h1_reads = 0

	def print(self):
		print('id =', self.id)
		print('actual_chrom =', self.actual_chrom)
		print('clust =', self.clust)
		print('actual_type =', self.actual_type)
		print('pred_type =', self.pred_type)
		print('****************************\nallele0\n****************************')
		if self.allele0 != None:
			self.allele0.print()
		else:
			print('allele0: None')

		print('****************************\nallele1\n****************************')
		if self.allele1 != None:
			self.allele1.print()
		else:
			print('allele1: None')

		print('****************************')


	def add_het_positions(self):
		assert (len(self.allele0.seq)==len(self.allele1.seq)), 'the lengths of the two bubble chains should be equal'
		self.het_positions = [i for i in range(len(self.allele0.seq)) if self.allele0.seq[i]!=self.allele1.seq[i]]
		assert (len(self.het_positions)>0), 'there should be at least one heterozygous position in the bubble'

	def add_allele(self, bubble_allele):
		assert (bubble_allele.id==0 or bubble_allele.id==1), 'bubble allele =, ' + str(bubble_allele.id) + 'should be 0 or 1'
		if bubble_allele.id == 0:
			self.allele0 = bubble_allele
		else:
			self.allele1 = bubble_allele

	def phase(self, ratio_h0_reads=(0.25, 0.75)):
		self.allele0.haplo_coverage = [0,0]
		self.allele1.haplo_coverage = [0,0]

		self.allele0.pred_haplo = None
		self.allele1.pred_haplo = None

		for aln in self.allele0.alignments:
			long_read = aln.long_read
			if long_read.pred_haplo == None:
				continue

			if aln.alt_bubble_allele_aln == None:
				continue

			al0_edit_dist = aln.edit_dist
			al1_edit_dist = aln.alt_bubble_allele_aln.edit_dist

			if al0_edit_dist < al1_edit_dist:
				# long read is better aligned to allele0
				self.allele0.haplo_coverage[long_read.pred_haplo] += 1

			if al0_edit_dist > al1_edit_dist:
				# long read is better aligned to allele1
				self.allele1.haplo_coverage[long_read.pred_haplo] += 1

		self.num_al0_h0_reads = self.allele0.haplo_coverage[0]+self.allele1.haplo_coverage[1]
		self.num_al0_h1_reads = self.allele0.haplo_coverage[1]+self.allele1.haplo_coverage[0]

		# filter out the set of bubbles that don't pass the phasing criteria
		if self.num_al0_h0_reads + self.num_al0_h1_reads == 0:
			return

		if ratio_h0_reads[0] < self.num_al0_h0_reads / (self.num_al0_h0_reads + self.num_al0_h1_reads) < ratio_h0_reads[1]:
			return

		if self.num_al0_h0_reads > self.num_al0_h1_reads:
			self.allele0.pred_haplo = 0
			self.allele1.pred_haplo = 1

		elif self.num_al0_h0_reads < self.num_al0_h1_reads:
			self.allele0.pred_haplo = 1
			self.allele1.pred_haplo = 0


class BubbleAllele:

	def __init__(self, id=None, name=None, bubble=None, seq=None, unitig_name=None):
		self.id = id
		self.name = name
		self.bubble = bubble
		self.seq = seq
		self.unitig_name = unitig_name
		self.rc_seq = None # reverse complement sequence
		self.actual_haplo, self.pred_haplo = None, None
		self.haplo_coverage = [0,0] # number of supporting (strand-seq or long) reads from haplo0 and haplo1
		self.km = 0
		self.alignments = [] #{}

	def print(self):
		print('allele id =', self.id)
		print('bubble_id =', self.bubble.id)
		print('km =', self.km)
		print('actual_haplo =', self.actual_haplo)
		print('pred_haplo =', self.pred_haplo)
		print('alignments =', self.alignments)

	def set_rc_seq(self):
		if self.rc_seq == None:
			self.rc_seq = reversecomp(self.seq)


	def get_haplotypes_edit_dist(self):

		haplo_edit_dist = [0,0]

		for aln in self.alignments:
			long_read_haplo = aln.long_read.pred_haplo

			if long_read_haplo == None:
				continue

			assert (long_read_haplo==0 or long_read_haplo==1), 'the predicted haplotype for long read ' + \
				aln.long_read.name + ' is ' + str(long_read_haplo) + ', should be 0 or 1'

			haplo_edit_dist[long_read_haplo] += aln.edit_dist

		return haplo_edit_dist


class LongRead:

	def __init__(self, name, seq=None):
		self.name = name
		self.seq = seq
		self.actual_haplo, self.pred_haplo = None, None
		self.actual_chrom, self.clust = None, None
		self.haplo0_edit_dist, self.haplo1_edit_dist = 0, 0
		self.num_haplo_bubbles = [0, 0] # num_haplo0_bubbles, num_haplo1_bubbles
		self.actual_type, self.pred_type = None, None
		self.alignments = []

	def print(self):
		print('name =', self.name)
		print('actual_chrom=', self.actual_chrom)
		print('clust =', self.clust)
		print('actual_haplo =', self.actual_haplo)
		print('pred_haplo =', self.pred_haplo)
		print('actual_type=', self.actual_type)
		print('alignments =', self.alignments)

	def link_alt_alignments(self):
		self.alignments.sort(key=lambda x: x.bubble_allele.bubble.id)
		for i in range(len(self.alignments)-1):
			if self.alignments[i].bubble_allele.bubble.id == self.alignments[i+1].bubble_allele.bubble.id and self.alignments[i].bubble_allele.id != self.alignments[i+1].bubble_allele.id:
				self.alignments[i].alt_bubble_allele_aln = self.alignments[i+1]
				self.alignments[i+1].alt_bubble_allele_aln = self.alignments[i]

	def set_haplotypes_edit_dist(self):

		haplo_edit_dist = [0, 0] # h0_dist, h1_dist, respectively

		for aln in self.alignments:
			bubble_allele_haplo = aln.bubble_allele.pred_haplo

			if bubble_allele_haplo == None:
				continue

			assert (bubble_allele_haplo==0 or bubble_allele_haplo==1), 'error in alignment ' + self.name + \
				', bubble ' + str(aln.bubble_allele.bubble.id) + ', allele ' + \
				str(aln.bubble_allele.id) + ': bubble allele predicted haplotype is ' + \
				str(bubble_allele_haplo) + ', should be 0 or 1'

			haplo_edit_dist[bubble_allele_haplo] += aln.edit_dist

		self.haplo0_edit_dist, self.haplo1_edit_dist = haplo_edit_dist[0], haplo_edit_dist[1]

	def set_alignments_edit_dist(self, q):
		for aln in self.alignments:
			aln.set_edit_dist(q)

	# this function is not used apparently!
	def set_num_haplo_bubbles(self):
		self.num_haplo_bubbles = [0, 0]

		for aln in self.alignments:
			bubble_allele = aln.bubble_allele
			if bubble_allele.id != 0:
				# we process the bothe alleles in one call of this function, so it suffices if we process only allele0
				continue

			bubble_allele0_haplo = bubble_allele.pred_haplo
			bubble_allele1 = bubble_allele.bubble.allele1
			bubble_allele1_haplo = bubble_allele1.pred_haplo

			if aln.alt_bubble_allele_aln == None:
				# this long read is only aligned to allele 0!
				continue

			if bubble_allele0_haplo == None:
				continue

			assert (bubble_allele1_haplo == 1-bubble_allele0_haplo), 'haplotypes of bubble' + str(bubble_allele.bubble.id) + 'should be 0 and 1'

			al0_edit_dist = aln.edit_dist
			al1_edit_dist = aln.alt_bubble_allele_aln.edit_dist #self.alignments[bubble_allele1].edit_dist

			if al0_edit_dist == al1_edit_dist:
				continue

			if al0_edit_dist < al1_edit_dist:
				self.num_haplo_bubbles[bubble_allele0_haplo] += 1

			else:
				self.num_haplo_bubbles[bubble_allele1_haplo] += 1

	def phase(self, ratio_h0_bubbles=(0.25, 0.75), min_haplotagged_bubbles=1):
		self.set_num_haplo_bubbles()

		self.pred_haplo = None

		# filter out if the long read doesn't pass the phasing criteria
		if self.num_haplo_bubbles[0] + self.num_haplo_bubbles[1] < min_haplotagged_bubbles:
			return

		if ratio_h0_bubbles[0] < self.num_haplo_bubbles[0] / (self.num_haplo_bubbles[0] + self.num_haplo_bubbles[1]) < ratio_h0_bubbles[1]:
			return

		if self.num_haplo_bubbles[0] > self.num_haplo_bubbles[1]:
			self.pred_haplo = 0
		elif self.num_haplo_bubbles[0] < self.num_haplo_bubbles[1]:
			self.pred_haplo = 1


class Alignment:
	def __init__(self, long_read, bubble_allele, edit_dist=0, long_read_start=None, long_read_end=None, bubble_start=None, bubble_end=None, strand = None, cigar=None):
		self.long_read, self.bubble_allele = long_read, bubble_allele
		self.long_read.alignments.append(self)
		self.bubble_allele.alignments.append(self)
		self.alt_bubble_allele_aln = None
		self.edit_dist = edit_dist
		self.long_read_start = long_read_start
		self.long_read_end = long_read_end
		self.bubble_start = bubble_start
		self.bubble_end = bubble_end
		self.strand=strand
		self.cigar = cigar
		self.bubble_allele_kmers=[]
		self.long_read_kmers=[]

	def set_edit_dist(self, q):

		for h in range(len(self.bubble_allele.bubble.het_positions)):
			het_pos = self.bubble_allele.bubble.het_positions[h]

			bubble_allele_seq = self.bubble_allele.seq
			if self.strand == "-":
				het_pos = len(self.bubble_allele.seq)-1-het_pos
				self.bubble_allele.set_rc_seq()
				bubble_allele_seq = self.bubble_allele.rc_seq

			assert (q <= het_pos <= len(bubble_allele_seq)-q-1), 'het position should be at least ' + str(q) + ' base pairs far from the start and end points of the bubble, ' \
				+ str(self.bubble_allele.bubble.id) + ', het_pos: ' + str(het_pos) + ', allele0: ' + str(self.bubble_allele.bubble.allele0.seq) + ', allele1: ' + str(self.bubble_allele.bubble.allele1.seq)

			### Apparently the query_start and end pos in paf file are also in the original strand!
			## look at bubble=1648373 ccs=m54329U_190827_173812/7079975/ccs in chunk 000
			if self.strand == '-':
				self.bubble_start, self.bubble_end = len(self.bubble_allele.seq)-1-self.bubble_end, len(self.bubble_allele.seq)-1-self.bubble_start
			###

			if not self.bubble_start+q <= het_pos <= self.bubble_end-q:
				# het pos is not fully covered in the alignment
				continue

			bubble_allele_kmer = bubble_allele_seq[het_pos-q:het_pos+q+1]
			self.bubble_allele_kmers.append(bubble_allele_kmer)

			long_read_kmer = get_reference_aln_substr(self.long_read.seq, bubble_allele_seq, self.long_read_start, self.bubble_start, self.cigar, het_pos-q, het_pos+q)
			self.long_read_kmers.append(long_read_kmer)

			# FIXME: some het kmers might have overlap with each other, and the edit distances will be computed wrongly (counting the overlapping mismatches twice)
			self.edit_dist += edit_distance(bubble_allele_kmer, long_read_kmer)

	def output_kmers(self):
		output_str = ""
		for i in range(len(self.bubble_allele_kmers)):
			output_str += str(self.bubble_allele.bubble.id) + '\t' \
							+ str(self.bubble_allele.id) + '\t' \
							+ self.long_read.name + '\t' \
							+ self.bubble_allele_kmers[i] + '\t' \
							+ self.long_read_kmers[i] + '\t' \
							+ str(self.edit_dist) + '\t' \
							+ str(self.bubble_allele.km)
			if i < len(self.bubble_allele_kmers)-1:
				 output_str += '\n'

		return output_str
