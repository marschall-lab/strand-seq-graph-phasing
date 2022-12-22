from whatshap.align import edit_distance

def reversecomp(seq):
	'''
	returns the reverse complement of string seq
	'''
	
	revcomp = {'a':'t', 'c':'g', 'g':'c', 't':'a', 'A':'T', 'C':'G', 'G':'C', 'T':'A'}
	rc = ''
	for i in range(len(seq)):
		rc = revcomp[seq[i]] + rc
	return rc

def print_dict_head(dic, num):
	for i in range(num):
		key = list(dic.keys())[i]
		print(key, ':', dic[key])

def print_dp_table(dp_table, str1=None, str2=None, max_dist_digits=2):
	def adjust_space(char, max_digits):
		n_digits = len(str(char))
		add_space = " "*(max_digits - n_digits)
		return add_space+str(char)
	
	if str2 != None:
		print("-", [adjust_space(char, max_dist_digits) for char in '-'+str2])
		
	for i in range(len(dp_table)):
		char = ""
		if str1 != None:
			char = "-" if i==0 else str1[i-1]			
			
		print(char, [adjust_space(dist, max_dist_digits) for dist in dp_table[i]])
		
def needleman_wunsch(str1, str2):
	m, n = len(str1), len(str2)
	edit_dist = [[j for j in range(n+1)]]
	
	for i in range(1, m+1):
		edit_dist.append([i])
		
		for j in range(1, n+1):				
			d = 0 if str1[i-1]==str2[j-1] else 1
			
			edit_dist[i].append(min(d+edit_dist[i-1][j-1], \
											1+edit_dist[i-1][j], # str1[i] aligned with gap \
											1+edit_dist[i][j-1]
											))
	assert (len(edit_dist)==m+1), 'number of rows in edit distance should be ' + str(m+1)
	for i in range(m+1):
		assert (len(edit_dist[i])==n+1), 'number of columns in edit distance should be ' + str(n+1)
		
	# backtracking
	str1_aln, aln, str2_aln = "", "", ""
	
	i, j = m, n
	
	while i > 0 or j > 0:
		#print('i, j =', i, j)
		#print(str1_aln+"\n"+aln+"\n"+str2_aln)
	
		back_dist = []
		if j != 0:
			back_dist.append(1+edit_dist[i][j-1])
			
		if i != 0:
			back_dist.append(1+edit_dist[i-1][j])
		
		if i != 0 and j != 0:
			d = 0 if str1[i-1]==str2[j-1] else 1
			back_dist.append(d+edit_dist[i-1][j-1])
			
		min_back_dist = min(back_dist)
		
		if i != 0 and j != 0 and min_back_dist == d+edit_dist[i-1][j-1]:
			# pos i in str1 is aligned with pos j in str2
			str1_aln = str1[i-1] + str1_aln
			str2_aln = str2[j-1] + str2_aln
			aln_char = "|" if str1[i-1]==str2[j-1] else "."
			aln = aln_char + aln
			i, j = i-1, j-1
			continue
		
		elif j != 0 and min_back_dist == 1+edit_dist[i][j-1]:
			# pos j in str2 is aligned with gap
			str1_aln = "-" + str1_aln
			str2_aln = str2[j-1] + str2_aln
			aln = " " + aln
			j = j-1
			continue
			
		else: # i!=0 and min_back_dist == 1+edit_dist[i-1][j]:
			# pos i in str1 is aligned with gap
			str1_aln = str1[i-1] + str1_aln
			str2_aln = "-" + str2_aln
			aln = " " + aln
			i = i-1
			
	print(str1_aln)
	print(aln)
	print(str2_aln)

	return edit_dist[m][n]
	
	
def get_alignment_from_cigar(ref, query, aln_ref_start_pos, aln_query_start_pos, cigar):
	'''
	This function returns the alignment from a cigar string
	
	Args:
		ref: the reference sequence
		query: the query sequence
		aln_ref_start_pos: 0-based alignment start position in the reference sequence
		aln_query_start_pos: 0-based alignment start position in the query sequence
		cigar: A string describing how the query sequence aligns to the reference sequence. It has integers followed by characters \in "MIDNSHP=X".
			The characters have the following meaning:
				M: alignment match (consumes_query=yes, consumes_reference=yes)
				I: insertion to the reference (consumes_query=yes, consumes_reference=no)
				D: deletion from the reference (consumes_query=no, consumes_reference=yes)
				N: skipped region from the reference (consumes_query=no, consumes_reference=yes)
				S: soft clipping (consumes_query=yes, consumes_reference=no)
				H: hard clipping (consumes_query=no, consumes_reference=no)
				P: padding: silent deletion from padded reference (consumes_query=no, consumes_reference=no)
				=: sequence match (consumes_query=yes, consumes_reference=yes)
				X: sequence mismatch (consumes_query=yes, consumes_reference=yes)
					
			consumes_query and consumes_reference indicate whether the CIGAR operation causes the alignment to step along the query sequence and the reference sequence respectively
				
	Returns:
		a tuple containing aligned reference, alignment sequence, and aligned query
	'''

	cigar_operations = 'MIDNSHP=X'	
	
	# defining consumes_query and consumes_ref for all cigar operations
	
	consumes_query = {'M': True, 'I': True , 'D': False, 'N': False, 'S': True , 'H': False, 'P': False, '=': True, 'X': True}
	consumes_ref =   {'M': True, 'I': False, 'D': True , 'N': True , 'S': False, 'H': False, 'P': False, '=': True, 'X': True}
		
	ref_pos = aln_ref_start_pos
	query_pos = aln_query_start_pos
	
	aln = ''
	ref_aln = ''
	query_aln = ''
		
	length = ''
	for i in range(len(cigar)):
		if cigar[i].isdigit():
			length = length + cigar[i]
			
		else:
			assert(len(length) > 0), 'there should not be two cigar operation characters next to each other'
			assert(cigar[i] in cigar_operations), "cigar " + cigar + " does not have valid characters"
			
			length = int(length)
			
			aln_char = '|'
			
			if consumes_ref[cigar[i]]:
				ref_aln += ref[ref_pos : ref_pos+length]
				ref_pos += length
			else:
				ref_aln += "-"*length
				aln_char = " "
				
			if consumes_query[cigar[i]]:
				query_aln += query[query_pos : query_pos+length]
				query_pos += length
			else:
				query_aln += "-"*length
				aln_char = " "
				
			aln = aln + aln_char * length
			length = ''

	assert(len(ref_aln)==len(query_aln)), 'the lengths of the aligned sequences should be the same'
	
	return ref_aln, aln, query_aln
	
	
def get_reference_aln_substr(ref, query, aln_ref_start_pos, aln_query_start_pos, cigar, query_start, query_end):
	'''
	This function returns the reference subsequence aligned to the given interval of the query sequence (query start and end are both included in the interval)
	
	Args:
		ref: 						@inherited from get_alignment_from_cigar
		query: 					@inherited from get_alignment_from_cigar
		aln_ref_start_pos: 	@inherited from get_alignment_from_cigar
		aln_query_start_pos: @inherited from get_alignment_from_cigar
		cigar: 					@inherited from get_alignment_from_cigar
		query_start: 0-based start of the interval of interest in the query sequence
		query_end  : 0-based end   of the interval of interest in the query sequence
	
	Returns:
		the reference subsequence aligned to the given interval of the query sequence
	'''

	# getting the alignments from the cigar string
	
	ref_aln, aln, query_aln = get_alignment_from_cigar(ref, query, aln_ref_start_pos, aln_query_start_pos, cigar)
	
	# computing the reference substr of interest
	
	ref_substr = ""
	query_pos = aln_query_start_pos
	for i in range(len(ref_aln)):
		if query_start <= query_pos <= query_end and ref_aln[i] != "-":
			ref_substr += ref_aln[i]
		
		if query_aln[i] != "-":
			# we read one char from query
			query_pos += 1
			
	return ref_substr
	

####################################################################################