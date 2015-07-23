#!/Users/panix/anaconda/bin/python

from genes import *
# Useful dicts

codon_table = dict(UUU='F', UUC='F', UUA='L', UUG='L', UCU='S', UCC='S', UCA='S', UCG='S', UAU='Y', UAC='Y', UAA='*', UAG='*', UGU='C', UGC='C', UGA='*', UGG='W', CUU='L', CUC='L', CUA='L', CUG='L', CCU='P', CCC='P', CCA='P', CCG='P', CAU='H', CAC='H', CAA='Q', CAG='Q', CGU='R', CGC='R', CGA='R', CGG='R', AUU='I', AUC='I', AUA='I', AUG='M', ACU='T', ACC='T', ACA='T', ACG='T', AAU='N', AAC='N', AAA='K', AAG='K', AGU='S', AGC='S', AGA='R', AGG='R', GUU='V', GUC='V', GUA='V', GUG='V', GCU='A', GCC='A', GCA='A', GCG='A', GAU='D', GAC='D', GAA='E', GAG='E', GGU='G', GGC='G', GGA='G', GGG='G')



class rna_seq_data:
    def __init__(self, file_path, organism, condition):
        self.file_path = file_path
        self.organism = organism
        self.condition = condition

    description = "This is a RNA sequence record. It contains information of the location and experimental significance"
    author = "To be set at a later stage. Probably the person who owns the data set will be author"



class genome():
	def __init__(self, organism, isolate):
		self.isolate = isolate
		self.organism = organism
		
	def get_snps(self, ref_isolate):
		con = lite.connect('Jon' + '_' + self.isolate + '.db')
		formatted_result = []
		with con:
			cur = con.cursor()
			cur.execute("SELECT Pos FROM " + ref_isolate + " WHERE Pos >= " + "0")
			result = cur.fetchall()	
			for i in result:
				not_tupp = str(i[0])
				formatted_result.append(not_tupp)
		return formatted_result
	
	def gene_list(self):
		con = lite.connect('Jon' + '_' + self.isolate + '.db')
		pretty_list = []
		with con:
			cur = con.cursor()
			cur.execute("SELECT " + "LOCUS" + " FROM " + "Gene")
			info = cur.fetchall()
		for i in info:
			this = str(i[0])
			pretty_list.append(this)
		return pretty_list

	
class bim_file:
    entry_count = 0

    def __init__(self, input_file):
        self.input_file = input_file

    def get_unique_elements(col):
        # return a unique set based on column
        subset_data = set()
        this_line = []
        for line in self.input_file:
            this_line = line.split()
            subset_data.add(this_line[col])
        print subset_data

    #This section of code has a bug, but im tired now and will just make it worse. so...
    #bim_var = bim_file(input_parser("/home/user/Exercises/plinktut/small_edited.bim"))
    #print bim_var.get_unique_elements()



class gene():
	def __init__(self, gene_ID, origin_organism):
		self.gene_ID = gene_ID
		self.origin_organism = origin_organism.organism
		if hasattr(origin_organism, 'species_ID'):
			self.species_ID = origin_organism.species_ID
		self.origin_isolate = origin_organism.isolate
		self.start = get_gene_item(session_info, self.origin_isolate, self.gene_ID, 'START')
		self.stop = get_gene_item(session_info, self.origin_isolate, self.gene_ID, 'STOP')
		self.strand = get_gene_item(session_info, self.origin_isolate, self.gene_ID, 'STRAND')
		self.layer = "gene"
		
	def ID(self):
		print self.gene_ID
		print self.origin_organism
		print self.origin_isolate
		print self.start, self.stop
		print self.strand
		print self.layer
		
	def layer(self):
		return self.layer

	def set_origin(self, origin_organism):
		self.origin_organism = origin_organism.organism
		self.origin_isolate = origin_organism.isolate

	def __str__ (self):
		return self.gene_ID
	
	def set_mRNA_seq(self, mRNA_seq):
		self.mRNA_seq = mRNA_seq

	def set_SNP_list(self, SNP_list):
		self.SNP_list = SNP_list

	def set_SNP_graph(self, SNP_graph):
		self.SNP_graph = SNP_graph

	def mRNA(self, *ref_isolate):
		"looking good, could add a check to make sure correct snp is selected"
		
		temp_gene_seq = get_seq_from_coords(session_info, self.origin_isolate, self.start, self.stop)
		
		for isolate in ref_isolate:
			con = lite.connect('Jon' + '_' + self.origin_isolate + '.db')
			formatted_result = []
			with con:
				cur = con.cursor()
				cur.execute("SELECT Pos FROM " + isolate + " WHERE Pos >= " + self.start + " AND Pos <=" + self.stop)
				result = cur.fetchall()	
				for i in result:
					not_tupp = str(i[0])
					formatted_result.append(not_tupp)
			#print formatted_result
			locus_alt_dict = {}

			for locus in formatted_result:
				cur = con.cursor()
				cur.execute("SELECT Alt FROM " + isolate + " WHERE Pos = " + locus)
				result = cur.fetchall()	
				for i in result:
					not_tupp = str(i[0])
					locus_alt_dict[locus] = not_tupp
			#print locus_alt_dict
			for locus in locus_alt_dict:
				#print locus, locus_alt_dict[locus]
				local_locus = int(locus) - int(self.start)
				#print local_locus
				temp_gene_seq = temp_gene_seq[0:(local_locus)] + locus_alt_dict[locus] + temp_gene_seq[(local_locus + 1):]

				
		if self.strand == '+':
			self.mRNA_seq = temp_gene_seq.replace('T','U')
		else:
			temp_gene_seq = reverse_compliment(temp_gene_seq)
			self.mRNA_seq = temp_gene_seq.replace('T','U')

		return self.mRNA_seq
	
	def set_start(self, start):
		self.start = start

	def set_stop(self, stop):
		self.stop = stop

	def set_REACTION(self, REACTION):
		self.REACTION = REACTION

	def seq(self):
		self.gene_seq = get_seq_from_coords(session_info, self.origin_isolate, self.start, self.stop)
		return self.gene_seq

	def snps(self, ref_isolate):
		con = lite.connect('Jon' + '_' + self.origin_isolate + '.db')
		formatted_result = []
		with con:
			cur = con.cursor()
			cur.execute("SELECT Pos FROM " + ref_isolate + " WHERE Pos >= " + self.start + " AND Pos <=" + self.stop)
			result = cur.fetchall()	
			for i in result:
				not_tupp = str(i[0])
				formatted_result.append(not_tupp)
		return formatted_result

	def aa(self, *snp_list):
		# this function will return an amino acid sequence for the given rna sequence assuming the
		# start position is at the first base
		self.aa = rna_to_aa(self.mRNA_seq)
		return self.aa
		
	def reaction(self):
		self.reaction = get_gene_item(session_info, self.origin_isolate, self.gene_ID, 'REACTION').split(',')
		return self.reaction
		
	def string_interaction(self):
		return get_STRING_interactions(self.gene_ID, self.species_ID)
		
	def get_base_count(self):
		# returns a dictionary counting the nucleotide content of a sequence
		base_dict = {}
		seq = self.gene_seq
		seq = seq.upper()
		seq = seq.strip()
		for base in seq:
			if base not in base_dict:
				base_dict[base] = 1
			else:
				base_dict[base] += 1

		return base_dict

	def get_GC_content(self):
		# returns a dictionary counting the nucleotide content of a sequence
		base_dict = {}
		seq = self.gene_seq
		seq = seq.upper()
		seq = seq.strip()
		for base in seq:
			if base not in base_dict:
				base_dict[base] = 1
			else:
				base_dict[base] += 1

		if 'G' not in base_dict:
			base_dict['G'] = 0
		if 'C' not in base_dict:
			base_dict['C'] = 0

		GC = (float(base_dict['G'] + base_dict['C']))/(len(seq))
		return GC
		
	def get_expression(self, experiment, comp_isolate, session_info):
		# Returns the expression data from an Array experiment
		result = get_gene_expression(session_info, comp_isolate, self.gene_ID.upper(), experiment)
		return result
	
	def get_average_expression(self, experiment, comp_isolate, session_info):
		result_dict = get_gene_expression(session_info, comp_isolate, self.gene_ID.upper(), experiment)
		total = 0
		count = 0
		if len(result_dict) > 0:
			for sample in result_dict:
				if sample['log_val'] != 'NA':
					value = float(sample['log_val'])
					total = total + value
					count = count + 1
			result = total / count
			return result
		else:
			return "NA"
	
	
class gff_feature():
	def __init__(self, seqname, start, stop):
		self.seqname = seqname
		self.start = start
		self.stop = stop

#    def get_uniprot_info(self):
        # Fills in info using ID. Takes from Uniprot.
        # Requires uniprot parser


	