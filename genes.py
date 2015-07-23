# The modules used in the nucleus are defined here
# The abreviation GAN refers to the Genbank Accession number, and acts
# to unify refords with a single identifier.

# Dev notes: Remember tp go back and close files where appropriate

# Import required modules

session_info = {}
session_info = 'Jon'

# For SQLite3
import sqlite3 as lite

import sys
import datetime
import csv
import numpy
import operator
import networkx as nx
import urllib2
import math

from subprocess import call, Popen, PIPE
import subprocess

from networkx import compose
from networkx import write_gml
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
from mpl_toolkits.mplot3d import Axes3D
import graphviz
from difflib import SequenceMatcher

# For importing classes
sys.path.append('/Users/panix/Dropbox/Programs/tools/Cell/Cell_core/')
#from rna import *
sys.setrecursionlimit(10000)
# gff and gene classes were removed from here

#-------------------------------------------------- Classes

class rna_seq_data:
    def __init__(self, file_path, organism, condition):
        self.file_path = file_path
        self.organism = organism
        self.condition = condition

    description = "This is a RNA sequence record. It contains information of the location and experimental significance"
    author = "To be set at a later stage. Probably the person who owns the data set will be author"



class genome():
	def __init__(self, organism, isolate, species_ID):
		self.isolate = isolate
		self.organism = organism
		self.species_ID = species_ID
		
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

#-------------------------------------------------- Database functions

def create_user_table(session_info):
	con = lite.connect(session_info['user_name'] + '.db')

	with con:
		cur = con.cursor()
		cur.execute("CREATE TABLE Organisms(Organism TEXT, Taxa_ID TEXT, Taxa_Heirarchy TEXT)")	
	

def create_table(session_info, organism):
	con = lite.connect(session_info['user_name'] + '_' + organism + '.db')
	
	with con:
		cur = con.cursor()
		cur.execute("CREATE TABLE Gene(Locus TEXT, Symbol TEXT, Synonym TEXT, Length INT, Start INT, Stop INT, Strand TEXT, Name TEXT, Chromosome TEXT, Genome_ontology TEXT, Enzyme_code TEXT, Kegg TEXT, Pathway TEXT, Reaction TEXT, Cog TEXT, Pfam TEXT, Operon TEXT)")
	
	with con:
		cur = con.cursor()
		cur.execute("CREATE TABLE Vcf(Chrom TEXT, Pos INT, Id TEXT, Ref TEXT, Alt TEXT, Qual TEXT, Filter TEXT, Info TEXT, Format TEXT, Fiveoseven TEXT, Origin TEXT)")
	
	with con:
		cur = con.cursor()
		cur.execute("CREATE TABLE Genome(Identity TEXT, Path TEXT)")


def add_to_user_database(session_info, table_name, data):
	# Where data is a dictionary. Generally.
	con = lite.connect(session_info['user_name'] + '.db') 
	# create a Session
	with con:
		cur = con.cursor()	
	# Create the entry
	# This needs to become dynamic for any dictionary, or specific to gene imports
		if table_name == 'Organism':
			cur.execute("SELECT Organism FROM Organisms")
			extant_organism = cur.fetchall()
			extant_organism_list = []
			for i in extant_organism:
				print i[0]
				j = i[0]
				extant_organism_list.append(j)
			print extant_organism_list
			if data['Organism'] not in extant_organism_list:
				cur.execute("INSERT INTO Organisms VALUES(" + "'" + data['Organism'] + "' , '" + data['Taxa_ID'] + "' , '" + data['Taxa_Heirarchy'] + "')")
				print "Importing: " + data['Organism']
			else:
				print 'Already in the database'

def drop_table(session_info, organism, table):
	con = lite.connect(session_info['user_name'] + '_' + organism + '.db')
	with con:
		cur = con.cursor()
		cur.execute("DROP TABLE " + table)
	

def add_to_biological_database(session_info, organism, table_name, dictionary):
	# Simple add entry to a database function
	# To be included in the import function some how
	con = lite.connect(session_info['user_name'] + '_' + organism + '.db')
	
	# create a Session
	with con:
		cur = con.cursor()

	# Create the entry
	# This needs to become dynamic for any dictionary, or specific to gene imports
		if table_name == 'Gene':
			for row in dictionary:
				cur.execute("INSERT INTO Gene VALUES(" + "'" + row['LOCUS'] + "' , '" + row['SYMBOL'] + "' , '" + row['SYNOYM'] + "' , '" + row['LENGTH'] + "' , '" + row['START'] + "' , '" + row['STOP'] + "' , '" + row['STRAND'] + "' , '" + row['NAME'] + "' , '" + row['CHROMOSOME'] + "' , '" + row['GENOME ONTOLOGY'] + "' , '" + row['ENZYME CODE'] + "' , '" + row['KEGG'] + "' , '" + row['PATHWAY'] + "' , '" + row['REACTION'] + "' , '" + row['COG'] + "' , '" + row['PFAM'] + "' , '" + row['OPERON'] + "')")
				print "Importing: ", row['LOCUS'], "\r",
			
		if table_name == 'VCF':
			for row in dictionary:
				#new_entry  = VCF(CHROM = row['CHROM'], POS = row['POS'], ID = row['ID'], REF = row['REF'], ALT = row['ALT'], QUAL = row['QUAL'], FILTER = row['FILTER'], INFO = row['INFO'], FORMAT = row['FORMAT'], FIVEOSEVEN = row['FIVEOSEVEN'], ORIGIN = row['ORIGIN'])
				cur.execute("INSERT INTO Vcf VALUES(" + "'" + row['CHROM'] + "' , '" + row['POS'] + "' , '" + row['ID'] + "' , '" + row['REF'] + "' , '" + row['ALT'] + "' , '" + row['QUAL'] + "' , '" + row['FILTER'] + "' , '" + row['INFO'] + "' , '" + row['FORMAT'] + "' , '" + row['FIVEOSEVEN'] + "' , '" + row['ORIGIN']  + "')")
				print "Importing: ", row['POS'], "\r",
				
		if table_name == 'Genome':
			cur.execute("INSERT INTO Genome VALUES(" + "'" + organism + "' , '" + dictionary + "')")
		
		if table_name == 'TF_targets':
			# from_gene_locus	from_symbol	from_product	to_gene_locus	to_symbol	to_product	count_links	cvg_over_mean	max_pos	type	subtype	dist2stt	fc_control	zscore	by_operon
			cur.execute("DROP TABLE " + table_name)
			cur.execute("CREATE TABLE " + table_name + "(locus TEXT, from_symbol TEXT, from_product TEXT, to_gene_locus TEXT, to_symbol TEXT, to_product TEXT, count_links TEXT, cvg_over_mean TEXT, max_pos TEXT, type TEXT, subtype TEXT, dist2stt TEXT, fc_control TEXT, zscore TEXT, by_operon TEXT)")			
			for row in dictionary:
				cur.execute("INSERT INTO TF_targets VALUES(" + "'" + row['from_gene_locus'] + "' , '" + row['from_symbol'] + "' , '" + row['from_product'] + "' , '" + row['to_gene_locus'] + "' , '" + row['to_symbol'] + "' , '" + row['to_product'] + "' , '" + row['count_links'] + "' , '" + row['cvg_over_mean'] + "' , '" + row['max_pos'] + "' , '" + row['type'] + "' , '" + row['subtype'] + "' , '" + row['dist2stt'] + "' , '" + row['fc_control'] + "' , '" + row['zscore'] + "' , '" + row['by_operon']  + "')")
				print "Importing: ", row['from_gene_locus'], "\r",
			
			
		else:
			SQL_table_string = "CREATE TABLE " + table_name + " ("
			
			for key in dictionary[0]:
				SQL_table_string = SQL_table_string + key + " TEXT, "
	
			SQL_table_string = SQL_table_string[:-2] + ")"
			
			cur.execute(SQL_table_string)
			
			for row in dictionary:
				SQL_add_string = 'INSERT INTO ' + table_name + ' VALUES('
				for key in row:
					SQL_add_string = SQL_add_string + "'" + row[key] + "' , "	
				SQL_add_string = SQL_add_string[:-3] + ')'
				
				print SQL_add_string
				cur.execute(SQL_add_string)
				#print "Importing: ", row, "\r"
			

def get_gene_details(session_info, organism, gene):
	# Pretty much get all we have on a particular gene. Return it as a large dictionary
	if type(session_info) is not str:
		session_info = session_info['user_name']
	else:
		session_info = session_info
	con = lite.connect(session_info + '_' + organism + '.db')
	
	gene_data = get_db_subset(session_info, organism, 'Gene', 'Locus', gene)
	gene_start = gene_data[0]['START']
	gene_stop = gene_data[0]['STOP']
	with con:
		cur = con.cursor()
		cur.execute("SELECT Pos FROM Vcf WHERE Pos >= " + gene_start + " AND Pos <=" + gene_stop)
		result = cur.fetchall()
		
	return result

def import_array_data(gds_code, gpl_code, user_name):
	"Import mocroarray data to a SQL table"
	
	sample_organism = "Mycobacterium_tuberculosis"
	table_name_exp = "gene_expression"
	print user_name
	
	script = '/Users/panix/Dropbox/Programs/tools/Cell/Cell_core/R_scripts/array_import.R'
	
	result = Popen(["Rscript", script, gds_code, gpl_code], stdout=PIPE)
	result = result.stdout.read()
	sample_block = True
	
	print "................Getting sample values per gene................"
	
	sample_value_dict = {}
	
	gene_loaded = False
	
	for line in result.split('\n'): 
		line = line.split()
		
		if '">>"' in line:
			sample_block = False
		
		# Get the gene name
		if "[1]" in line and sample_block == True:
			gene = line[-1]
			gene_loaded = True
			block_count = 1
			
		# Get the sample values
		if 	"[1]" not in line and gene_loaded == True:
			print "gene"
			print gene
			print line
			print block_count
			if block_count %2 == 0:
				values = line
			else:
				sample = line
				
			if block_count % 2 == 0:
				temp_dict = dict(zip(sample, values))
				print temp_dict
				gene = str(gene)
				gene = gene.replace('"', '')
				if gene in sample_value_dict.keys():
					sample_value_dict[gene] = dict(sample_value_dict[gene].items() + temp_dict.items())
				else:
					sample_value_dict[gene] = temp_dict
			
			block_count = block_count + 1
			
	
	#print sample_value_dict['RV1511']['GSM71989']
	
	print "................Getting isolate sample dict................"
	
	
	sample_block = False
	isolate_sample_dict = {}
	
	for line in result.split('\n'):
		#print line
		
		if sample_block == True and len(line.split()) > 1:
			line = line.split()
			#print line[-1]
			isolate_sample = line[-1].split("_")
			
			isolate_sample_dict[isolate_sample[1].replace('"', '')] = isolate_sample[0].replace('"', '')
			
			
		else:
			1 == 1
		
		if ">>" in line:
			sample_block = True
	
	print isolate_sample_dict
	
	# Combining the two dicts into one for the SQL database
	
	print "-------------------------------"
	
	sql_import_list_of_dict = []
	
	for gene in sample_value_dict:
		print gene
		print gds_code
		for sample in sample_value_dict[gene]:
			sql_import_temp_dict = {}
			print sample + sample_value_dict[gene][sample]
			sql_import_temp_dict['experiment'] = gds_code
			sql_import_temp_dict['gene'] = gene
			isolate = isolate_sample_dict[sample]
			sql_import_temp_dict['isolate'] = isolate
			sql_import_temp_dict['sample'] = sample
			sql_import_temp_dict['log_exp'] = sample_value_dict[gene][sample]
			# To be replaced with actual condition later
			sql_import_temp_dict['condition'] = 'NA'
			sql_import_list_of_dict.append(sql_import_temp_dict)
		
	print sql_import_list_of_dict[1116]
	
	# Adding the data to a SQL database
	
	#add_to_biological_database(session_info, sample_organism, table_name_exp, sql_import_list_of_dict)
	
	con = lite.connect(session_info['user_name'] + '_' + sample_organism + '.db')
	print "connected to db"
	
	# create a Session
	with con:
		cur = con.cursor()
		print "Creating table"
		cur.execute("CREATE TABLE expression(gene TEXT, isolate TEXT, sample TEXT, log_val TEXT, condition TEXT, experiment TEXT)")
		print "Adding data"
		
		for row in sql_import_list_of_dict:
			cur.execute("INSERT INTO expression VALUES('" + row['gene'] + "' , '" + row['isolate'] + "' , '" + row['sample'] + "' , '" + row['log_exp'] + "' , '" + row['condition'] + "' , '" + row['experiment'] + "')")
			print "Importing: ", row['gene'], "\r",
	
	return "Import complete \n"


def get_gene_item(session_info, organism, gene, info):
	# Pretty much get all we have on a particular gene. Return it as a large dictionary
	if type(session_info) is not str:
		session_info = session_info['user_name']
	else:
		session_info = session_info
	con = lite.connect(session_info + '_' + organism + '.db')
	
	gene_data = get_db_subset(session_info, organism, 'Gene', 'Locus', gene)
	gene_info = gene_data[0][info]
	return gene_info

def db_status(session_info, organism):
	con = lite.connect(session_info['user_name'] + '_' + organism + '.db')
	with con:
		cur = con.cursor()
		cur.execute("SELECT Name FROM sqlite_master WHERE type='table'")
		tables = cur.fetchall()
		cur.execute("SELECT Count(*) FROM Vcf")
		rows = cur.fetchall()
		print "Tables:"
		for l in tables:
			cur.execute("SELECT Count(*) FROM " + l[0])
			rows = cur.fetchall()
			print l[0], rows[0][0], ' entries'

def get_all_from(session_info, organism, table, category):
	con = lite.connect(session_info['user_name'] + '_' + organism + '.db')
	pretty_list = []
	with con:
		cur = con.cursor()
		cur.execute("SELECT " + category + " FROM " + table)
		info = cur.fetchall()
	for i in info:
		this = str(i[0])
		pretty_list.append(this)
	return pretty_list
	
	
def db_status_user(session_info):
	print 'User: ' + session_info['user_name']
	# Get summary of user db contents
	con = lite.connect(session_info['user_name'] + '.db')
	
	with con:
		cur = con.cursor()
		cur.execute("SELECT Name FROM sqlite_master WHERE type='table'")
		tables = cur.fetchall()
		cur.execute("SELECT * FROM Organisms")
		result = cur.fetchall()

	print 'Tables:'
	for item in tables:
		print item[0]

	print 'Organisms: '
	for item in result:
		print item[0]

def get_available_organisms(session_info):
	# Get available organisms and return as a list
	orgnaisms = []
	con = lite.connect(session_info['user_name'] + '.db')
	with con:
		cur = con.cursor()
		cur.execute("SELECT * FROM Organisms")
		result = cur.fetchall()
	return result
		
def get_db_subset(session_info, organism, table, column, key, *spec_attr):
	# Returns a dictionary containing entires where the column meets the constraint
	if type(session_info) is dict:
		session_info = session_info['user_name']
	
	con = lite.connect(session_info + '_' + organism + '.db')
	formatted_result = []
	with con:
		cur = con.cursor()
		if len(spec_attr) == 0:
			cur.execute("SELECT * FROM " + table + " WHERE " + column +"='" + key + "'")
		else:
			cur.execute("SELECT " + spec_attr[0] + " FROM " + table + " WHERE " + column +"='" + key + "'")
		result = cur.fetchall()

		if table == "Gene":
			for i in result:
				temp_dict = {'LOCUS': str(i[0]), 'SYMBOL': str(i[1]), 'SYNONYM': str(i[2]), 'LENGTH': str(i[3]), 'START': str(i[4]), 'STOP': str(i[5]), 'STRAND': str(i[6]), 'NAME': str(i[7]), 'CHROMOSOME': str(i[8]), 'GENOME_ONTOLOGY': str(i[9]), 'ENZYME_CODE': str(i[10]), 'KEGG': str(i[11]), 'PATHWAY': str(i[12]), 'REACTION': str(i[13]), 'COG': str(i[14]), 'PFAM': str(i[15]), 'OPERON': str(i[16])}
				formatted_result.append(temp_dict)
				
		if table == "Vcf":
			for i in result:
				temp_dict = {'CHROM': str(i[0]), 'POS': str(i[1]), 'ID': str(i[2]), 'REF': str(i[3]), 'ALT': str(i[4]), 'QUAL': str(i[5]), 'FILTER': str(i[6]), 'INFO': str(i[7]), 'FORMAT': str(i[8]), 'FIVEOSEVEN': str(i[9]), 'ORIGIN': str(i[10])}
				formatted_result.append(temp_dict)
		
		if table == "TFs":
			for i in result:
				if len(spec_attr) == 0:
					temp_dict = {'from_gene_locus': str(i[10]), 'from_symbol': str(i[13]), 'from_product': str(i[7]), 'to_gene_locus': str(i[12]), 'to_symbol': str(i[6]), 'to_product': str(i[8]), 'count_links': str(i[3]), 'cvg_over_mean': str(i[0]), 'max_pos': str(i[11]), 'type': str(i[14]), 'subtype': str(i[2]), 'dist2stt': str(i[5]), 'fc_control': str(i[4]), 'zscore': str(i[1]), 'by_operon': str(i[9])}
					formatted_result.append(temp_dict)
				else:
					temp_dict = {spec_attr[0]:str(i[0])}
					formatted_result.append(temp_dict)
		else:
			for i in result:
				entry = []
				for value in i:
					entry.append(str(value))
				formatted_result.append(entry)
		return formatted_result

def get_edge_value(session_info, organism, table, column1, key1, column2, key2, edge_val):
	# Returns a dictionary containing entires where the column meets the constraint
	if type(session_info) is dict:
		session_info = session_info['user_name']
	
	con = lite.connect(session_info + '_' + organism + '.db')
	formatted_result = []
	with con:
		cur = con.cursor()
		exec_string = "SELECT " + edge_val + " FROM " + table + " WHERE " + column1 + "='" + key1 + "' AND " + column2 + "='" + key2 + "'"
		#print exec_string
		cur.execute(exec_string)
		result = cur.fetchall()
	
	for i in result:
		formatted_result.append(str(i[0]))
	return formatted_result

def get_genome_path(session_info, organism):
	#change to attribute
	con = lite.connect(session_info + '_' + organism + '.db')
	with con:
		cur = con.cursor()
		cur.execute("SELECT Path FROM Genome")
		result = cur.fetchall()
	return str(result[0][0])
		
def get_annotations_from_locus_range(session_info, organism, start, stop):
	# Returns a list of genes and features found between the given start and stop locations
	con = lite.connect(session_info + '_' + organism + '.db')
	formatted_result = []
	with con:
		cur = con.cursor()
		cur.execute("SELECT Locus FROM Gene WHERE Start <=" + stop + " AND Stop >=" + start)
		result = cur.fetchall()
		for i in result:
			not_tupp = str(i[0])
			formatted_result.append(not_tupp)
	return formatted_result

def get_SNPs(session_info, organism, gene):
	con = lite.connect(session_info['user_name'] + '_' + organism + '.db')
	
	gene_data = get_db_subset(session_info, organism, 'Gene', 'Locus', gene.gene_ID)
	gene_start = gene_data[0]['START']
	gene_stop = gene_data[0]['STOP']
	with con:
		cur = con.cursor()
		cur.execute("SELECT Pos FROM Vcf WHERE Pos >= " + gene_start + " AND Pos <=" + gene_stop)
		result = cur.fetchall()
		
	return result

def get_SNPs_vcf(vcf_filepath, start, stop):
	"Returns the SNPs from a vcf file from between a given range"
	input_list = input_parser(vcf_filepath)
	
	result_list = []
	for entry in input_list:
		if (int(entry['POS']) > int(start) and int(entry['POS']) < int(stop)):
			result_list.append(entry)
	
	return result_list

def get_gene_expression(session_info, isolate, gene, experiment):
	# Establish connection
	con = lite.connect(session_info['user_name'] + '_' + session_info['organism'] + '.db')
	
	with con:
		cur = con.cursor()
		cur.execute("SELECT * FROM expression WHERE gene='" + gene + "' AND isolate='" + isolate + "' AND experiment='" + experiment + "'")
		result = cur.fetchall()
		formatted_result = []
		for i in result:
				temp_dict = {'gene': str(i[0]), 'isolate': str(i[1]), 'sample': str(i[2]), 'log_val': str(i[3]), 'condition': str(i[4]), 'experiment': str(i[5])}
				formatted_result.append(temp_dict)
		return formatted_result


# TESTING
#add_to_biological_database(user_name,database_title,'NC_007779.1', 'gtgaaacgat','982-1190')
#this = retrieveFromDatabase('hi')
#print this

#-------------------------------------------------- Online functions

def ensemble_query(GBA):
	print 'To be done'

def get_STRING_interactions(gene, species_ID):
	# requires urllib2
	
	try:
		urllib2.urlopen('http://string-db.org/api/tsv/interactors?identifier=' + gene + '&species=' + species_ID)
	except:
		print gene + ' not found in string' 
		result_list = []
		return result_list
	else:	
		response = urllib2.urlopen('http://string-db.org/api/tsv/interactors?identifier=' + gene + '&species=' + species_ID)
		html = response.read()
		html = html.split()[1:]
		result_list = []
		for entry in html:
			result_list.append(entry.split('.')[1])
		return result_list
	


#-------------------------------------------------- Special case functions
def convert_strings_to_class(list, isolate):
	classy_list = []
	for gene_string in list:
		classy_list.append(gene(gene_string, isolate))
	return classy_list


def convert_obj_to_dict_list(list_of_objects):
	# Converts a list of objects to a list of dictionaries
	new_list = []
	for i in list_of_objects:
		new = i.__dict__
		new_list.append(new)
	return new_list

# Running R Scripts

def get_expression_data(gds_soft_file_path, gpl_soft_file_path, gene):
	"This function will return the expression of a gene in a set of experiments and uses R"
	
	# Test R script running
	script = '/Users/panix/Dropbox/Programs/tools/Cell/Cell_core/R_scripts/get_gene_expression.R'
	gene = gene.upper()
	result = Popen(["Rscript", script, gds_soft_file_path, gpl_soft_file_path, gene], stdout=PIPE)
	result = result.stdout.read()
	
	# Format into something useful
	i = 0
	headder_list = []
	data_list = []
	
	for line in result.split('\n'):
		if i % 2 == 0:
			line = line.split()
			for headder in line:
				headder_list.append(headder)
		else:
			line = line.split()
			for value in line:
				data_list.append(value)	
		i = i + 1	
		
	expression_dict = {}
	
	# Combine the two lists into a dict
	i = 0
	for header in headder_list:
		expression_dict[header] = data_list[i]
		i = i + 1
			
	return expression_dict
	
def get_GSM_isolate(gds_soft_file_path, gpl_soft_file_path, GSM_number):
	"Gets the isolate from which a gsm entry was isolated from"
	script = '/Users/panix/Dropbox/Programs/tools/Cell/Cell_core/R_scripts/get_GSM_strain.R'
	result = Popen(["Rscript", script, gds_soft_file_path, gpl_soft_file_path, GSM_number], stdout=PIPE)
	result = result.stdout.read()
	return result[4:]
	
def import_array_data(gds_code, gpl_code, user_name):
	"Import mocroarray data to a SQL table"
	
	print user_name
	
	script = '/Users/panix/Dropbox/Programs/tools/Cell/Cell_core/R_scripts/array_import.R'
	
	result = Popen(["Rscript", script, gds_code, gpl_code], stdout=PIPE)
	result = result.stdout.read()
	for line in result.split('\n'):
		print line 
	return "done"

def ncRNA_import():
	ncRNA_add_data = csv.reader(open('/Volumes/HDD/Users/Admin/Work/Genomes/M_tuberculosis/Paida_ncRNA_Data/Final_Mtb_ncRNA/ncRNA-Table 1.csv','r'))

	# headder_file = "/Volumes/HDD/Users/Admin/Work/Genomes/M_tuberculosis/H37Rv/ncRNA/RNAPredator_predictions/headders.csv"
	working_directory = "/Volumes/HDD/Users/Admin/Work/Genomes/M_tuberculosis/H37Rv/ncRNA/RNAPredator_predictions/"
	file_names = "ncRNA_Pred_tar_"

	# Filtering parameters
	energy_cutoff = 0


	# Final list that will contain all the data
	all_ncRNA_list = []

	# First, using the list of the predicted ncRNA

	for item in ncRNA_add_data:
		# Now working with element n in the list
		#print item[0]
	
		ncRNA_dict = {}
	
		# Now working with current ncRNA target file
		if item[0] != 'NUMBER' and len(item[0]) > 0:
			if float(item[3].replace(',','')) > float(item[4].replace(',','')):
				# add info to ncRNA dict
				ncRNA_dict['number'] = item[0]
				ncRNA_dict['strain'] = item[1]
				ncRNA_dict['start'] = item[5]
				ncRNA_dict['stop'] = item[6]
				ncRNA_dict['orientation'] = item[8]
			
				# open the ncRNA target file
				ncRNA_target_file = open('/Volumes/HDD/Users/Admin/Work/Genomes/M_tuberculosis/H37Rv/ncRNA/ncRNA_tar_' + item[0] + '.txt','r')
			
				# Convert file object to list
				this_ncRNA_list = []
			
				for line in ncRNA_target_file:
					this_ncRNA_list.append(line)
			
				# Extract the predicted target info
				i = 0
				block = False
				ncRNA_target_list = []
			
			
				while i < len(this_ncRNA_list):
					if this_ncRNA_list[i][0:4] == "Rank":
						block = True
					if block == True and this_ncRNA_list[i][0:1] != "\n":
						# Dictionary for each target for this ncRNA
						ncRNA_target_dict = {}
						ncRNA_target_line = this_ncRNA_list[i].split()
						ncRNA_target_dict['Rank'] = ncRNA_target_line[0]
						ncRNA_target_dict['Gene'] = ncRNA_target_line[1]
						ncRNA_target_dict['Synonym'] = ncRNA_target_line[2]
						ncRNA_target_dict['Energy'] = ncRNA_target_line[3]
						ncRNA_target_dict['Pvalue'] = ncRNA_target_line[4]
						ncRNA_target_dict['sRNA_start'] = ncRNA_target_line[5]
						ncRNA_target_dict['sRNA_stop'] = ncRNA_target_line[6]
						ncRNA_target_dict['mRNA_start'] = ncRNA_target_line[7]
						ncRNA_target_dict['mRNA_stop'] = ncRNA_target_line[8]
						ncRNA_target_list.append(ncRNA_target_dict)
					else:
						block = False
				
					i = i + 1
			ncRNA_dict['TargetRNA2_predictions'] = ncRNA_target_list[1:]
		
			# importing the RNAPredator data
			if item[1] == "Mycobacterium tuberculosis H37Rv" and len(item[6]) > 0:
				Pred_target_file = csv.reader(open(working_directory + file_names + item[0] + ".csv", 'r'))
			
				RNA_pred_list = []
			
				for line in Pred_target_file:
					if float(line[9]) < energy_cutoff:
						Pred_target_dict = {}
						Pred_target_dict['Rank'] = line[0]
						Pred_target_dict['NC'] = line[1]
						Pred_target_dict['Accession'] = line[2]
						# Split into start and stop
						Pred_target_dict['Co-ordinates'] = line[3]
						Pred_target_dict['Interaction'] = line[4]
						Pred_target_dict['mRNA_start'] = line[5]
						Pred_target_dict['mRNA_stop'] = line[6]
						Pred_target_dict['sRNA_start'] = line[7]
						Pred_target_dict['sRNA_stop'] = line[8]
						Pred_target_dict['Energy'] = line[9]
						Pred_target_dict['zscore'] = line[13]
						Pred_target_dict['Annotation'] = line[14]
						# Appended gene name bugfix here
						gene_formatted = "Rv" + line[15].split("Rv")[1]
						Pred_target_dict['Gene'] = gene_formatted
					
						RNA_pred_list.append(Pred_target_dict)
					
		
			ncRNA_dict['RNAPredator_predictions'] = RNA_pred_list
	
		all_ncRNA_list.append(ncRNA_dict)
	
	#print all_ncRNA_list[12]['RNAPredator_predictions']
	#print all_ncRNA_list[12]['number']
	return all_ncRNA_list

#-------------------------------------------------- Maths functions


def list_mean(list):
	"returns the mean of a list of values"
	mean = sum(list) / float(len(list))
	return mean


def list_SD(list):
	"Returns the variance of a list"
	val_sum = sum(list)
	sum_diff = 0
	for value in list:
		print value
		sum_diff = sum_diff + ((value - list_mean(list)) ** 2)
		STD = math.sqrt(sum_diff/(len(list)-1))
		
	return STD
	

#-------------------------------------------------- General functions

def compliment_Base(base_in):
	# Returns the compliment of the given base. Assuning DNA
	# Non-DNA bases are returned as a empty string

	base_in = base_in.upper()
	if base_in == "A":
		return "T"
	if base_in == "T":
		return "A"
	if base_in == "C":
		return "G"
	if base_in == "G":
		return "C"
	if base_in == "N":
		return "N"
	else:
		return ""

def compliment_DNA(seq_in):
	# Returns the string with the complimentry base pairs
	seq_out = ""
	for char in seq_in:
		seq_out += compliment_Base(char)
	return seq_out

def reverse(text):
	# returns a reversed string using recursion
	if len(text) <= 1:
		return text
	return reverse(text[1:]) + text[0]

def reverse_compliment(seq_in):
	"Returns the reverse compliment of a sequence"
	comp_seq = compliment_DNA(seq_in)
	rev_comp_seq = reverse(comp_seq)
	return rev_comp_seq

def rna_to_aa(seq):
	# this function will return an amino acid sequence for the given rna sequence assuming the
	# start position is at the first base
	codon_table = dict(UUU='F', UUC='F', UUA='L', UUG='L', UCU='S', UCC='S', UCA='S', UCG='S', UAU='Y', UAC='Y', UAA='*', UAG='*', UGU='C', UGC='C', UGA='*', UGG='W', CUU='L', CUC='L', CUA='L', CUG='L', CCU='P', CCC='P', CCA='P', CCG='P', CAU='H', CAC='H', CAA='Q', CAG='Q', CGU='R', CGC='R', CGA='R', CGG='R', AUU='I', AUC='I', AUA='I', AUG='M', ACU='T', ACC='T', ACA='T', ACG='T', AAU='N', AAC='N', AAA='K', AAG='K', AGU='S', AGC='S', AGA='R', AGG='R', GUU='V', GUC='V', GUA='V', GUG='V', GCU='A', GCC='A', GCA='A', GCG='A', GAU='D', GAC='D', GAA='E', GAG='E', GGU='G', GGC='G', GGA='G', GGG='G')
	i = 0
	aa_out = ''
	seq = seq.upper().strip()
	reading_frame = True

	while i < len(seq) and reading_frame == True:
		codon = seq[0+i:3+i]
		aa_out = aa_out + codon_table[codon]
		i = i + 3

		# deal with rna seq that is not a multiple of 3
		if len(seq[0+i:3+i]) < 3:
			reading_frame = False

	return aa_out

def file_string_replace(file_path, from_str, to_str):
	infile = open(file_path)
	out_file_path = file_path[:-4] + '_replaced.txt'
	print out_file_path
	outfile = open(out_file_path, 'w')

	replacements = {from_str:to_str}

	for line in infile:
		for src, target in replacements.iteritems():
			line = line.replace(src, target)
		outfile.write(line)
	infile.close()
	outfile.close()

def input_parser(file_path):
	if file_path[-3:] == ".fa" or file_path[-6:] == ".fasta":
		input_file = open(file_path, "r")

		# set variables
		sequence_details = ""
		sequence = ""

		for line in input_file:
			if line[0] == ">":
				sequence_details = line
			else:
				sequence += line
		sequence_details = sequence_details.strip()
		sequence = sequence.strip()
		gene_ID_dict = {"gene_details" : sequence_details[1:], "DNA_seq" : sequence}

		return gene_ID_dict

	if file_path[-4:] == ".bim":
		input_file = open(file_path, "r")
		return input_file
		
	if file_path[-4:] == ".csv":
		data_table = csv.reader(open(file_path, 'r'), delimiter=',')
		data_matrix = list(data_table)
		result=numpy.array(data_matrix).astype("float")
		return result
		
	if file_path[-4:] == ".txt":
		list_of_dicts  = []
		reader = csv.DictReader(open(file_path, 'r'), delimiter='\t')
		for row in reader:
			list_of_dicts.append(row)
		return list_of_dicts
		
	if file_path[-4:] == ".vcf":
		list_of_dicts  = []
		# Deal with random info at start
		in_file = open(file_path, 'r')
		entry_label = file_path
		for line in in_file:
			if not line.startswith('#'):
				entries = line.split('\t')
				entry_dict = {'CHROM':entries[0], 'POS':entries[1], 'ID':entries[2], 'REF':entries[3], 'ALT':entries[4], 'QUAL':entries[5], 'FILTER':entries[6], 'INFO':entries[7], 'FORMAT':entries[8], 'FIVEOSEVEN':entries[9], 'ORIGIN':entry_label} 
				list_of_dicts.append(entry_dict)
		return list_of_dicts
	
	if file_path[-5:] == ".diff":
		list_of_dicts  = []
		reader = csv.DictReader(open(file_path, 'r'), delimiter='\t')
		for row in reader:
			list_of_dicts.append(row)
		return list_of_dicts
	
def csv_to_dict_importer(file_path):
	"Importing a variety of csv files and adding them to a dictionary while settign the headders"
	print file_path
	print "to be done"
	
def parse_soft_file(filepath):
	"Parses a soft file to be used "
	print filepath
	
# SOFT parse testing
#parse_soft_file('/Users/Admin/Work/Genomes/M_tuberculosis/Expression_data/GDS1552/GDS1552_full.soft')


def get_unique_col(in_file_path):
	input_data = input_parser(in_file_path)
	output_data = []
	subset_data = set()
	this_line = []
	for line in input_data:
		output_data.append(line.split())
		this_line = line.split()
		subset_data.add(this_line[1])

	print output_data
	print subset_data
	print "Original entries: " + str(len(output_data))
	print "Unique SNP IDs: " + str(len(subset_data))


def sliding_window(sequence, window_size, step_size, function):
	# Analises a sequence in windows, with defined steps and window sizes.
	# Returns a list of dictionaries
	returned_feature_list = []
	
	if function == 'GC':
		# generic bit
		i = 0
		frames = (len(sequence) - window_size)/step_size
		while i < (len(sequence) - window_size + 1):	
			feature = {}
			feature['start'] = i
			feature['stop'] = i + step_size
			seq_slice = sequence[i:i+window_size]
			i = i + step_size
			#print seq_slice
			feature['attribute'] = 'GC=' + str(get_GC_content(seq_slice)) + ';'
			returned_feature_list.append(feature)
	
	return returned_feature_list

def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

#seq_temp = 'GACTGACGATCAGCGAGCGACGACTATCTATCGACGTACGGACTAGCGAGCATCGAGGCGCGGCATTATATTTAAACCCGG'
#window_size = 3
#Steppppa = 3
#funk = 'GC'

#heythere =  sliding_window(seq_temp, window_size, Steppppa, funk)
#print heythere
		
#-------------------------------------------------- Sequence function
# These functions deal with sequence composition

def get_codon_usage():
	#Get the codon usage of a sequence, assuming sequence is in frame
	#Return the result as a dictionary, per codon, absolute count
	print 'TBD'

def get_base_count(seq):
	# returns a dictionary counting the nucleotide content of a sequence
	base_dict = {}
	seq = seq.upper()
	seq = seq.strip()
	for base in seq:
		if base not in base_dict:
			base_dict[base] = 1
		else:
			base_dict[base] += 1

	return base_dict

def get_GC_content(seq):
	# Get the GC content of a sequence
	base_dict = get_base_count(seq)
	if 'G' not in base_dict:
		base_dict['G'] = 0
	if 'C' not in base_dict:
		base_dict['C'] = 0

	GC = (float(base_dict['G'] + base_dict['C']))/(len(seq))
	return GC
	

def k_mer_comp(seq, k_mer):
	"Calculates the k-mer composition of sequence. Returns a dict - count"
	k_mer = int(k_mer)
	seq = seq.upper()
	result_dict = {}
	i = 0
	while i < (len(seq) - k_mer + 1):
		if result_dict.get(seq[i:i+k_mer]) == None:
			result_dict[(seq[i:i+k_mer])] = 1
		else:
			result_dict[(seq[i:i+k_mer])] = result_dict[(seq[i:i+k_mer])] + 1
		
		i = i + 1
	
	return result_dict


def generate_profile(list_of_sequences, k_mer):
	"Generates a profile for the given sequences and returns as a dict in percentage composition"
	
	result_profile_dict = {}
	# This depends on the input format...
	for sequence in list_of_sequences:
		seq_dict = k_mer_comp(sequence, k_mer)
		total = 0
		for key in seq_dict:
			seq_len = (len(sequence) - k_mer + 1)
			#print key, seq_dict[key]
			#print key, 'comp = ', str((1.0 / seq_len*seq_dict[key])*100)
			perc_comp = (1.0 / seq_len*seq_dict[key])*100
			total = total + perc_comp
			
			if result_profile_dict.get(key) != None:
				 result_profile_dict[key].append(perc_comp)
			else:
				result_profile_dict[key] = [perc_comp]
			
		if str(total) != "100.0":
			print "ERROR IN THE PERCENTAGE COMPOSITION OF THE generate_profile FUNCTION"
			print total

	
	return result_profile_dict
	

def compare_sequence(profile1, profile2, comparison_type, parameter):
	"Compares the composition of two different sequnces"
	# The parameter is going to be made flexible / optional
	# Handle the comparison of kmers that are not found in both profiles
	
	if comparison_type =="Kolmogorov_Smirnov":
		# Perform a Kolmogorov-Smirnov two sample test that two data samples come from the same distribution. Note that we are not specifying what that common distribution is		
		# Works with R. Duh.
		
		# First do the calculation to get the values
		query_pro = generate_profile(profile1, parameter)
		reference_pro = generate_profile(profile2, parameter)
		
		result_list = []

		
		# Then preform the test for each kmer
		for k_mer in query_pro:
			query_pro_string = ''
			reference_pro_string = ''
			
			print query_pro[k_mer]
			for value in query_pro[k_mer]:
				query_pro_string = query_pro_string + str(value) + ','
				
				
			if k_mer in reference_pro:
				for value in reference_pro[k_mer]:
					query_pro_string = query_pro_string + str(value) + ','
				print reference_pro[k_mer]
			else:
				reference_pro[k_mer] = '0'
				reference_pro_string = '0'
				print reference_pro[k_mer]
			
			# R script 
			script = '/Users/Admin/Dropbox/Programs/tools/Cell/Cell_core/R_scripts/Kolmogorov_Smirnov_test.R'
			
			query_pro_string = query_pro_string[:-1]
			reference_pro_string = reference_pro_string[:-1]
			result = Popen(["Rscript", script, query_pro_string, reference_pro_string], stdout=PIPE)
			print ''
			print 'Analysis output for ', k_mer
			result = result.stdout.read()
			print result
			result_list.append(result)
			#result = "pie"
		return result_list
		
		
	if comparison_type == 'SD':
		# This is pretty much TEAR. Profile 1 is considered the query set / sequence
		query_pro = generate_profile(profile1, parameter)
		reference_pro = generate_profile(profile2, parameter)
		result_list = []
		
		for k_mer in query_pro:
			if len(query_pro[k_mer]) > 1:
				print list_mean(query_pro[k_mer])
			
		for k_mer in reference_pro:
			if len(reference_pro[k_mer]) > 1:
				print list_SD(reference_pro[k_mer])
		
	return result_list

#test_list_1 = ['CACACCCACACACACACCCACACCACACACACCCACA', 'CCCACACACCACACCCACAAACACACACCAC', 'CAACAAACA', 'CAACAACCCACA', 'ACACCCAAAAAAAAAA','aAAAAAAC', 'CCAACCAAAAAAAAA']
#test_list_2 = ['GACTACGGCATACCACCCGATCGACT','AGCTAGCHACTTAGCACAACCAGCATCATC','ACGCACAAACCCACACT', 'GCTAGCGCGCGGCGGCGACACCCACACAACACACAGGTTAACCGCGATCGAGCATCTCGGCGCGGATCTAGCGACT', 'GCATCGGCGCCGGATCTACGACTCTGCTCTGCTCGCTCGAGCAGCTACGATCAGCTAGTACAAAAAAAAAAA', 'GCATCGGCGCCGGATCTCACACAAACAACGACTCTGCTCTGCTCGCTCGAGCAGCTACGATCAGCTAGTACAAAAAAAAAAA', 'GCATCGGCAAAAAAAGCCGGATCTACGACTCTGCTCTGCTCGCCCCCCCCTCGAGCAGCTACGATCAGCTAGTACAAAAAAAAAAA', 'GCATCGGCGCCGGCAATCTACGACTCTGCTCTGCTCGCTCGAGCAGCTACGATCAGCTAGTACAAAAAAAAAAA', 'GCATCGGCGCCGGATCTACGACTCTGCTCTGCTCGCTCGAGCAGCTACGATCAGCTAGTACAAAAAAAAAAA']
#
#answer_to_life = compare_sequence(test_list_1, test_list_2, "Kolmogorov_Smirnov", 2)
#for answer in answer_to_life:
#	print answer

#-------------------------------------------------- Genome functions 

def get_seq_from_coords(session_info, organism, start, stop, *contig):
	"Returns a sequence from the given position in an organisms genome. Contig required if not supercontig"
	genome_path = get_genome_path(session_info, organism)
	genome_dict = input_parser(genome_path)
	seq = genome_dict["DNA_seq"]
	seq = seq.replace('\n', '')
	seq_slice = seq[(int(start)-1):(int(stop))]
	seq_slice = seq_slice
	return seq_slice

#print 'get seq test'
#seq_start = '3483886'
#seq_stop = '3484356'
#true_seq = 'GTGTATTCTGGTTGTTGGATAAACAACCAGAATGGGGAGACGCGGGTGGGCGAGGACTCGCTGGAGGATCTGGAGCAGCGGCGAGCGCGACTGTATGACCAGTTGGCCGCGACCGGCGATTTCCGGCGCGGCTCGATCAGTGAGAACTATCGCCGCTGCGGCAAGCCCAATTGTGTGTGCGCGCAAGAGGGTCACCCCGGGCATGGGCCGCGATATTTGTGGACGCGCACGGTGGCCGGGCGGGGTACCAAGGGGCGGCAGCTCTCGGTCGAGGAGGTGGACAAGGTGCGCGCCGAGTTGGCCAACTATCACCGTTTCGCGCAGGTCAGTGAGCAGATCGTGGCGGTCAACGAGGCGATCTGCGAGGCCCGCCCACCGAACCCGGCGGCCACGGCGCCCCCGGCCGGCACAACGGGGCACAAAAAAGGGGGCTCTGCGACCAGATCGCGGCGGAGTTCACCGCCGAGGTAG'
#test_seq = get_seq_from_coords("Jon", "CDC1551", seq_start, seq_stop, 'contig')

#print 'Observed length: ' + str(len(test_seq))
#print 'Expected lenght: ' + str(3484356-3483886)
#print 'True seq length: ' + str(len(true_seq))
#print test_seq



def create_blast_db(session_info, organism):
	# Create a blast database for a genome
	con = lite.connect(session_info + '_' + organism + '.db')
	# create a Session
	with con:
		cur = con.cursor()	
		cur.execute("SELECT Path FROM Genome")
	result = cur.fetchall()
	genome_path = str(result[0][0])
	database_name = organism + "_BlastDB"
	return call(["makeblastdb", "-in", genome_path, "-out", database_name, "-dbtype", "nucl"])
	
# print create_blast_db('Jon', 'F11')

def conduct_blast_search(session_info, organism, sequence, e_val):
	# Search the selected database for the sequence
	# The input sequence must be in fasta format... which may prove a problem
	# Return a ordered result as a list of dictionaries... to be decided
	blast_db_directory = "./"
	blast_database = organism + "_BlastDB"
	
	# Make a fasta file with the sequence in it
	fasta_dir = "/Users/Admin/Dropbox/Programs/tools/Cell/Cell_core/"
	temp_fasta = open(fasta_dir + 'blasta.fa','w')
	temp_fasta.write(">temp_fasta_file\n")
	temp_fasta.write(sequence + "\n")
	temp_fasta.close()
	
	# Doing the blast 
	#  blastn -db F11_BlastDB -query /Users/Admin/Dropbox/Programs/tools/Cell/Cell_core/blasta.fa
	#print "blasting"
	full_directory = fasta_dir + "blasta.fa"
	blast_result_raw = Popen(["blastn", "-db", blast_database, "-query", full_directory, "-evalue", e_val, "-outfmt", "7"], stdout=PIPE)
	blast_result = blast_result_raw.stdout.read()
	#print len(blast_result)
	
	# Formatting into a list
	new_list = blast_result.split('\n')
	
	# Remove empty list items
	new_list = filter(None, new_list)
	
	# make into dict & add to list
	List_of_result_lists = []
	for entry in new_list:
		if entry[0] != "#":
			split_entry = entry.split('\t')
			List_of_result_lists.append(split_entry)
			
	# Lists are in the form:
	# query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
	return List_of_result_lists


#test_seq = "GTGTATTCTGGTTGTTGGATAAACAACCAGAATGGGGAGACGCGGGTGGGCGAGGACTCGCTGGAGGATCTGGAGCAGCGGCGAGCGCGACTGTATGACCAGTTGGCCGCGACCGGCGATTTCCGGCGCGGCTCGATCAGTGAGAACTATCGCCGCTGCGGCAAGCCCAATTGTGTGTGCGCGCAAGAGGGTCACCCCGGGCATGGGCCGCGATATTTGTGGACGCGCACGGTGGCCGGGCGGGGTACCAAGGGGCGGCAGCTCTCGGTCGAGGAGGTGGACAAGGTGCGCGCCGAGTTGGCCAACTATCACCGTTTCGCGCAGGTCAGTGAGCAGATCGTGGCGGTCAACGAGGCGATCTGCGAGGCCCGCCCACCGAACCCGGCGGCCACGGCGCCCCCGGCCGGCACAACGGGGCACAAAAAAGGGGGCTCTGCGACCAGATCGCGGCGGAGTTCACCGCCGAGGTAG"
#print conduct_blast_search("Jon", "CDC1551", test_seq, '0.000000001')

# temp testing
#print get_unique_col("/home/user/Exercises/plinktut/small_edited.bim")
#print input_parser('/Users/Admin/Dropbox/Programs/tools/Cell/Modules/Input_files/data.csv')

#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1)
#plt.show()

def output_parser(GAN, out_format, file_path, optional_string):
	# Outputs a object in the desired format to a location
	# GAN is a dictionary, later to be a object
	if len(optional_string) == 0:
		optional_string = ""

	if out_format == "fasta":
		output_file = open(file_path+GAN["GAN"]+ optional_string + ".fa", "w")
		first_fasta_ID_line = ">" + GAN['GAN'] + optional_string + "\n"
		output_file.write(first_fasta_ID_line)
		fasta_sequence_line = GAN["DNA_seq"] + "\n"
		output_file.write(fasta_sequence_line)

	else:
		print "Unknown format. "

def search_for_seq(seq, search_type):
	# Use different searches to find a string
	print "will do soon"

def aligh_sequences(auery, reference, e_val):
	"Aligns two sequences returning the details of the alignment"
	temp_dir = ""
	# Save sequences to temp fasta files
	
	# Align the two sequences
	
	# Parse the results
	

#-------------------------------------------------- Testing area

#sequence_1 = DNA_seq("TTTTTTTTTTTTTTTTTT","G45")
#sequence_2 = DNA_seq('GCGATCGAGCGACGACTAGC', 'H553')

#sequence_1.displayInfo()

#seq = "GCATCAGCAGGCTAGC"
#print get_GC_content(sequence_1.sequence)
#print get_GC_content(sequence_2.sequence)

#DNA = "CTGAAa"
#print compliment_DNA(DNA)
#print reverse("Hello")

#-------------------------------------------------- Annotation functions

def get_feature_count(dict, field, term):
	count = 0 
	for row in dict:
		if term in row[field]:
			count += 1
	return count

def cluster_on_funct(dict, field):
	cluster_count = {}
	temp_list = []
	for row in dict:
		row_string = row[field]
		temp_list = row_string.split(',')
		for item in temp_list:
			if item in cluster_count:
				cluster_count[item] += 1
			else:
				cluster_count[item] = 1
	return cluster_count

def get_homologues(session_info, gene_list, source, target):
	# Return the initial list of dictionaries for each gene WITH the homologues as a new key.
	# Homologues selected based on BLAST alignment - E-value must be adjusted

	# Get the sequence based on gene name and location
	results_list = []
	count = 0
	for gene_x in gene_list:
		gene_x_full = get_db_subset(session_info, source, "Gene", "Locus", gene_x)
		#print gene_x_full
		
		gene_x_seq = get_seq_from_coords(session_info, source, gene_x_full[0]['START'], gene_x_full[0]['STOP'])
		
		# Blast the sequences
		
		BLAST_hit_list = conduct_blast_search(session_info, target, gene_x_seq, '0.0000001')
		
		
		# Get genes at hit locations
		for hit in BLAST_hit_list:
			homologues_list = []
			homologues_list.append(get_annotations_from_locus_range(session_info, target, hit[8], hit[9]))
		
		# Add the sequences homologues
		gene_x_full[0]["HOMOLOGUES"] = homologues_list
		
		results_list.append(gene_x_full)
		count = count + 1
		print str(count) + " / " + str(len(gene_list))
	return results_list

#test_gene_list = ['Rv0397', 'Rv0398c', 'Rv0399c']

#get_sum_list = get_homologues("Jon", test_gene_list, "H37Rv", "F11")

def add_annotation(organism, ID, start, stop, keywords_list, *meta_dict):
	"Adds a custom annotation to the add_annotation table of the organism"
	# Meta dict is optional. Maybe keywords_list too
	print "tic"
	
	

#-------------------------------------------------- Grouping

# What is the input? is it a list?? come now past Jon, document.
# Ok so its a list

def intersect(a, b):
	# Check if string
	if isinstance(a, list) == True:
		c = []
		d = []
		for i in a:
			if isinstance(i, str) == False:
				i = i.gene_ID
				c.append(i)
			else:
				c.append(i)
		for j in b:
			if isinstance(j, str) == False:
				j = j.gene_ID
				d.append(j)
			else:
				d.append(j)
		return list(set(c) & set(d))
	else:
		# Treat as a list of dicts
		return 'th'



def create_homologue_graph(list_of_homo_gene_dicts):
	# Takes in a list of gene dictionaries and returns a graph object with the edges based on homology
	Graph=nx.DiGraph()
	
	# Add nodes
	for gene_list in list_of_homo_gene_dicts:
		for gene in gene_list:
			#print gene[0]['LOCUS']
			Graph.add_node(gene[0]['LOCUS'])

		# Add edges
		for gene in gene_list:
			if len(gene[0]['HOMOLOGUES'][0]) != 0:
				for homologue in gene[0]['HOMOLOGUES'][0]:
					Graph.add_edges_from([(gene[0]['LOCUS'], homologue)])
	return Graph


# Get lists of unique, connected nodes ect
# Relative to a particular genome
# Needs to work for variable genome combinations

def get_edge_lists(Graph, gene_dict):
	# Returns a dictionary where keys are the number of edges and the values are lists of genes with the amount of edges
	edge_nodes_result = {}

	for node in gene_dict:
		node_degree = Graph.out_degree(node[0]['LOCUS'])
		#print node[0]['LOCUS'] + ' ' + str(node_degree)
		
		if str(node_degree) in edge_nodes_result:
			#print str(node_degree)  + "found, appending"
			edge_nodes_result[str(node_degree)].append(node[0]['LOCUS'])
		else:
			#print str(node_degree) + 'not found, adding'
			edge_nodes_result[str(node_degree)] = [node[0]['LOCUS']]
			
	return edge_nodes_result



def union(a, b):
	c = []
	d = []
	for i in a:
		if isinstance(i, str) == False:
			i = i.gene_ID
			c.append(i)
		else:
			c.append(i)
	for j in b:
		if isinstance(j, str) == False:
			j = j.gene_ID
			d.append(j)
		else:
			d.append(j)
	return list(set(c) | set(d))
	
def difference(a, b):
	c = []
	d = []
	for i in a:
		if isinstance(i, str) == False:
			i = i.gene_ID
			c.append(i)
		else:
			c.append(i)
	for j in b:
		if isinstance(j, str) == False:
			j = j.gene_ID
			d.append(j)
		else:
			d.append(j)
	return list(set(c) - set(d))
	
def symDifference(a, b):
	c = []
	d = []
	for i in a:
		if isinstance(i, str) == False:
			i = i.gene_ID
			c.append(i)
		else:
			c.append(i)
	for j in b:
		if isinstance(j, str) == False:
			j = j.gene_ID
			d.append(j)
		else:
			d.append(j)
	return list(set(c) ^ set(d))
	
def unique_list(list):
	new_list = []
	for i in list:
		if i not in new_list:
			new_list.append(i)
	return new_list

#-------------------------------------------------- Graph related functions 

def graph_on_feature_old(list_of_obj, feature):
	# Given a list of objects for each gene? Or is it a list of gene objects? looks like objects...
	new_dict = {}
	awesome_graph_object = nx.Graph(label=feature)
	# BRUTE FORCE!!!!
	for item in list_of_obj:
		if feature == 'REACTION':	
			classifiers = item.REACTION.split(',')
			new_dict[item.gene_ID] = classifiers
			
	# all against all comparison
	
	
	for item in list_of_obj:
		awesome_graph_object.add_node(item.gene_ID)
	
	for gene in new_dict:
		for other_gene in new_dict:
			for reaction in new_dict[gene]:
				if reaction in new_dict[other_gene]:
					awesome_graph_object.add_edge(gene, other_gene)
					
			

	return awesome_graph_object

def organism_graph(session_info, organism):
	# Creates a base graph with all the needed gene info
	
	org_network = nx.MultiDiGraph()
	
	classy_list = []
	for gene_string in organism.gene_list():
		classy_list.append(gene(gene_string, organism))
	
	for gene_obj in classy_list:
		#print gene_obj.gene_ID
		org_network.add_node(gene_obj.gene_ID, type='gene', origin_organism = gene_obj.origin_organism, origin_isolate = gene_obj.origin_isolate, start = gene_obj.start, stop = gene_obj.stop, strand = gene_obj.strand)
	
	return org_network
		

def organism_TF_graph(session_info, organism, TF_table, list_of_obj, *threshold):
	if len(threshold) > 0:
		thres = float(threshold[0])
	else:
		thres = 0
	
	TF_list = get_all_from(session_info, organism.isolate, TF_table, "from_gene_locus")
	TF_list = unique_list(TF_list)
	
	# print 'point 2'
	TF_graph = nx.DiGraph(label='TFs')
	for gene in list_of_obj:
		print gene.gene_ID
		if gene.gene_ID in TF_list:
			#print 'point 3'
			TF_graph.add_node(gene.gene_ID)
			TF_graph.node[gene.gene_ID]["type"] = "TF"
			isolate = organism.isolate
			gene_details = get_db_subset(session_info['user_name'], isolate, TF_table, "from_gene_locus", gene.gene_ID)
			for target in  gene_details:
				#print target['to_gene_locus'], target['zscore']
				if target['zscore'] != "NULL":
					if abs(float(target['zscore'])) > thres:
						TF_graph.add_edge(target['from_gene_locus'],target['to_gene_locus'],weight=float(target['zscore']),type="TF")
					
	return TF_graph


def add_TF_to_graph(in_graph, session_info, organism, TF_table, list_of_obj, *threshold):
	if len(threshold) > 0:
		thres = float(threshold[0])
	else:
		thres = 0
	
	TF_list = get_all_from(session_info, organism.isolate, TF_table, "from_gene_locus")
	TF_list = unique_list(TF_list)
	
	isolate = organism.isolate
	
	# print 'point 2'
	# TF_graph = in_graph
	for gene in list_of_obj:
		print gene.gene_ID
		if gene.gene_ID in TF_list:
			# Assess whether node exists or not
			if gene.gene_ID in in_graph:
				appended_type = in_graph.node[gene.gene_ID]["type"] + ",TF"
				print appended_type
				in_graph.node[gene.gene_ID]["type"] = appended_type
				gene_details = get_db_subset(session_info['user_name'], isolate, TF_table, "from_gene_locus", gene.gene_ID)			
				# Adding the new edges
				for target in  gene_details:
					#print target['to_gene_locus'], target['zscore']
					if target['zscore'] != "NULL":
						if abs(float(target['zscore'])) > thres:
							in_graph.add_edge(target['from_gene_locus'],target['to_gene_locus'],weight=float(target['zscore']),type="TF")
			
			else:
				in_graph.add_node(gene.gene_ID)
				in_graph.node[gene.gene_ID]["type"] = "TF"
				gene_details = get_db_subset(session_info['user_name'], isolate, TF_table, "from_gene_locus", gene.gene_ID)
				for target in  gene_details:
					#print target['to_gene_locus'], target['zscore']
					if target['zscore'] != "NULL":
						if abs(float(target['zscore'])) > thres:
							in_graph.add_edge(target['from_gene_locus'],target['to_gene_locus'],weight=float(target['zscore']),type="TF")
					
	return in_graph

def organism_ncRNA_graph(organism, ncRNA_dict_list, *threshold):
	'''This function exists to import data generated from the RnaPredator ncRNA target prediction site and incorporate it into a network'''
	# Returns a graph object
	
	print organism.isolate
	print ncRNA_dict_list
	
	# Create graph object
	ncRNA_graph = nx.DiGraph(label='ncRNA')
	
	# Replacing the given ncRNA blank with data ---> CHANGE LATER <---
	ncRNA_dict_list = ncRNA_import()
	
	
	for ncRNA in ncRNA_dict_list:
		#print ncRNA.keys()
		if 'number' in ncRNA.keys():
			print ncRNA['number']
			ncRNA_node = 'ncRNA_' + str(ncRNA['number'])
			for target in ncRNA['RNAPredator_predictions']:
				ncRNA_graph.add_edge(ncRNA_node, target['Gene'], energy=target['Energy'], zscore=target['zscore'], type='ncRNA')
			ncRNA_graph.node[ncRNA_node]['type'] = 'ncRNA'
			
	return ncRNA_graph

def add_ncRNA_edges(graph, ncRNA_dict_list, *threshold):
	
	# Replacing the given ncRNA blank with data ---> CHANGE LATER <---
	ncRNA_dict_list = ncRNA_import()
	
	for ncRNA in ncRNA_dict_list:
		#print ncRNA.keys()
		if 'number' in ncRNA.keys():
			#print ncRNA['number']
			ncRNA_node = 'ncRNA_' + str(ncRNA['number'])
			for target in ncRNA['RNAPredator_predictions']:
				graph.add_edge(ncRNA_node, target['Gene'], energy=target['Energy'], zscore=target['zscore'], type='ncRNA')
			graph.node[ncRNA_node]['type'] = 'ncRNA'
			
	return graph
			
def add_snps_to_graph(graph, gene_class_list, query_organism):
	"""Adds the number of SNPs in the genes to a graph"""
	for gene in gene_class_list:
		graph.node[gene.gene_ID]['SNPs'] = len(gene.snps(query_organism))

def add_string_interaction_edges_web(in_graph, isolate, gene_obj_list):
	'''Adds edges to existing gene nodes based on string database interactions'''
	# Graph based on STRING interactions with genes not included in the list
	
	for gene in gene_obj_list:
		print gene.gene_ID
		edge_list = get_STRING_interactions(gene.gene_ID, isolate.species_ID)
		print edge_list
		for target in edge_list:
			if target != gene.gene_ID:
				in_graph.add_edges_from([(gene.gene_ID,target,{"type":"PPI"})])
	return in_graph
	
def add_string_interaction_edges(session_info, in_graph, isolate, edge_type, gene_obj_list):
	'''Adds edges to existing gene nodes based on string database interactions from a local sql database'''
	# Graph based on STRING interactions with genes not included in the list
	
	for gene in gene_obj_list:
		print gene.gene_ID
		edge_list = get_db_subset(session_info, "H37Rv", "STRING_83332", 'protein1', gene.gene_ID, 'protein2')
		#print edge_list
		for target in edge_list:
			if target != gene.gene_ID:
				weight_val = get_edge_value(session_info, gene.origin_isolate, "STRING_83332", 'protein1', gene.gene_ID, 'protein2', target[0], edge_type)
				weight_val = weight_val[0]
				target = target[0]
				in_graph.add_edge(gene.gene_ID, target, type="edge_type", weight=weight_val)
	return in_graph

def add_node_attribute(session_info, graph, attribute):
	print 'test'
	for this_node in graph:
		# print this_node, graph.node[this_node]['type']
		if graph.node[this_node]['type'] == 'gene':
			data = get_gene_item(session_info, graph.node[this_node]['origin_isolate'], str(this_node), attribute).split(',')
			graph.node[this_node][attribute] = data
	
def graph_on_dir_string(list_of_obj):
	# Graph based on STRING interactions
	string_graph = nx.DiGraph(label='STRING')
	for gene in list_of_obj:
		string_graph.add_node(gene.gene_ID)
		print species_identifier[gene.origin_organism]
		print get_STRING_interactions(gene.gene_ID, species_identifier[gene.origin_organism])

	for gene in list_of_obj:
		edge_list = get_STRING_interactions(gene.gene_ID, species_identifier[gene.origin_organism])
		for target in edge_list:
			if target != gene.gene_ID:
				if target in string_graph:
					string_graph.add_edges_from([(gene.gene_ID,target)])
	return string_graph

def graph_on_all_string(list_of_obj):
	# Graph based on STRING interactions with genes not included in the list
	ext_string_graph = nx.DiGraph(label='STRING')
	for gene in list_of_obj:
		ext_string_graph.add_node(gene.gene_ID, type='gene')
		#print get_STRING_interactions(gene.gene_ID, gene.species_ID)
	
	for gene in list_of_obj:
		print gene.gene_ID
		edge_list = get_STRING_interactions(gene.gene_ID, gene.species_ID)
		print edge_list
		for target in edge_list:
			if target != gene.gene_ID:
				ext_string_graph.add_edges_from([(gene.gene_ID,target,{"type":"PPI"})])
	return ext_string_graph

def get_edge(self, v1, v2):
    try:
        e = self[v1][v2] # order shouldn't matter
        print("edge exists")
        return e
    except KeyError:
        print("edge does not exist")
        return None

def test_edge_exists(self, v1, v2):
    try:
        e = self[v1][v2] # order shouldn't matter
        return True
    except KeyError:
        return False

def graph_on_reaction(list_of_obj):
	"""Create graph based on common functions"""
	# Use a multigraph so multiple edges can exist between nodes
	reaction_graph = nx.MultiGraph(label='REACTION')
	for gene in list_of_obj:
		print gene.gene_ID
		reaction_graph.add_node(gene.gene_ID)
	
	# Create edge dictionary
	edge_dict = {}
	for gene in list_of_obj:
		if len(gene.reaction()) > 0:
			for pred_reaction in gene.reaction:
				if len(pred_reaction) > 0: 
					print "pred_reaction: " + pred_reaction
					if pred_reaction not in edge_dict:
						temp_gene_list = []
						temp_gene_list.append(gene.gene_ID)
						edge_dict[pred_reaction] = temp_gene_list
					else:
						edge_dict[pred_reaction].append(gene.gene_ID)
	
	# Convert edge dictionary to edges with labels
	for k in edge_dict:
		print k, edge_dict[k]
		if len(edge_dict[k]) > 1:
			for reacting_gene in edge_dict[k]:
				i = 0
				while i < len(edge_dict[k]):
					if reacting_gene != edge_dict[k][i]:
						if test_edge_exists(reaction_graph, reacting_gene, edge_dict[k][i]) == False:
							reaction_graph.add_edges_from([(reacting_gene,edge_dict[k][i])], reaction=k)
					i = i + 1
	print reaction_graph.edges()
	
	#print test_edge_exists(reaction_graph, 'Rv2228c', 'Rv0054')
	
	return reaction_graph

def graph_on_feature(list_of_obj):
	# v2.0
	# Given a list of gene/annotation/node objects, we need to create a network linking these nodes
	graph_dict = {}
	
	# Graph based on STRING interactions	
	graph_dict['string'] = graph_on_dir_string(list_of_obj)

	# Graph based on STRING interactions with genes not included in the list
	graph_dict['ext_String'] = graph_on_all_string(list_of_obj)
	
	
	# Graph on reaction data 
	
	graph_dict['reaction'] = graph_on_reaction(list_of_obj)
			
	# all against all comparison

	return graph_dict

def merge_graph(graphA, graphB, *feature):
	# Combines two graphs, based on common nodes.
	# This functions exists as compose(G1,G2)
	print "use compose(G1,G2)"

def create_gene_object_list(session_info, list_of_diff_genes_string, organism):
	# This function takes a list of genes based on whats in the gene file and converts them 
	# to gene objects. These gene objects are populated with info like their start and stop 
	# locations and interactions. 
	newnew = convert_strings_to_class(list_of_diff_genes_string)

	newnew[0].set_start(get_gene_item(session_info, organism, newnew[0].gene_ID, 'START'))

	# Add needed info 
	for i in newnew:
		i.set_start(get_gene_item(session_info, organism, i.gene_ID, 'START'))
		i.set_stop(get_gene_item(session_info, organism, i.gene_ID, 'STOP'))

	for item in newnew:
		item.set_REACTION(get_gene_item(session_info, organism, item.gene_ID, 'REACTION'))

	return newnew

def rough_GSA_on_nodes(graph, node_label, edge_label, set):
	#'''A Work in progress'''
	#print '-----------------------'
	
	# First, lets get an idea of what the superset looks like (the set of every gene in the organism)
	
	profile_result_dict = {}
	
	superset_count_dict = {}
	
	# Make the dict with keys and 0 values
	for test_node in graph.nodes():
		if graph.node[test_node]['type'] == 'gene':
			for unique_set in graph.node[test_node][set]:
				superset_count_dict[unique_set] = 0
	#print superset_count_dict
	
	set_total = 0
	
	# Do the superset count
	for test_node in graph.nodes():
		if graph.node[test_node]['type'] == 'gene':
			for unique_set in graph.node[test_node][set]:
				superset_count_dict[unique_set] = superset_count_dict[unique_set] + 1
				# For values in superset_profile_dict, this is the total
				set_total = set_total + 1
	#print superset_count_dict
	
	superset_profile_dict = {}
	for superset in superset_count_dict.keys():
		superset_profile_dict[superset] = float(superset_count_dict[superset]) / float(set_total)
				
	# Do the per node analysis
	for test_node in graph.nodes():
		if 'type' in graph.node[test_node].keys():
			#print graph.node[test_node]['type']
			if graph.node[test_node]['type'] == node_label:
				
				# From here we work on THIS ncRNA / node
				# print test_node
				test_node_name = str(test_node)
				# print len(graph.neighbors(test_node)), 'attached nodes'
				
				# Calculating at varying thresholds
				window_size = 0.2
				start_threshold = 0
				stop_threshold = 2.0
				
				threshold = start_threshold
				
				gene_set_dict = {}
				
				while threshold < stop_threshold:
				
					target_gene_subset = []
					
					for connected_node in graph.neighbors(test_node):
						'''Here is working with each connected node'''
						#print connected_node, graph.edge[test_node][connected_node]
						
						for target_edge in graph.edge[test_node][connected_node]:
							'''Here is working with each connected node'''
							#print target_edge
							#print graph.edge[test_node][connected_node][target_edge]
							if graph.edge[test_node][connected_node][target_edge]['type'] == node_label:
								#print 'this far'
								if abs(float(graph.edge[test_node][connected_node][target_edge][edge_label])) > threshold:
									#print 'even further'
									#print graph.edge[test_node][connected_node][target_edge][edge_label]
									target_gene_subset.append(connected_node)
					
					gene_set_dict[threshold] = target_gene_subset
					# Update window
					threshold = threshold + window_size
				
				# We now have a dict containing thresholds as keys, and the genes at those thresholds for this node
				
				# The following dict is the profile composition for this node
				node_profile_dict = {}
				
				for test_node in graph.nodes():
					if graph.node[test_node]['type'] == 'gene':
						for unique_set in graph.node[test_node][set]:
							node_profile_dict[unique_set] = 0
			
				#print node_profile_dict
				
				threshold_profile_result_dict = {}
				
				for current_threshold in gene_set_dict.keys():
					'''sdf'''
					#print current_threshold, len(gene_set_dict[current_threshold])
					# Generate empty profile dict
					
					node_threshold_profile_dict = {}
					for test_node in graph.nodes():
						if graph.node[test_node]['type'] == 'gene':
							for unique_set in graph.node[test_node][set]:
								node_threshold_profile_dict[unique_set] = 0
					
					# Populate with values
					for tar_node in gene_set_dict[current_threshold]:
						for unique_set in graph.node[tar_node][set]:
							node_threshold_profile_dict[unique_set] = node_threshold_profile_dict[unique_set] + 1
					
					#print node_threshold_profile_dict
					
					node_threshold_profile_difference_dict = {}
					
					diff_ratio_dict = {}
					
					# Calculate ratios for the sets
					tot_targets = 0
					for pred_target in node_threshold_profile_dict.keys():
						tot_targets = tot_targets + node_threshold_profile_dict[pred_target]
					
					target_node_ratio_dict = {}
					for pred_target in node_threshold_profile_dict.keys():
						target_node_ratio_dict[pred_target] = float(node_threshold_profile_dict[pred_target]) / float(tot_targets)
					
					# Calculate the difference between the sets
					for category in node_threshold_profile_dict.keys():
						node_threshold_profile_difference_dict[category] = superset_count_dict[category] - node_threshold_profile_dict[category]
					
					# Calculate the differences between the profiles
					for pred_target in target_node_ratio_dict.keys():
						diff_ratio_dict[pred_target] = float(target_node_ratio_dict[pred_target]) - float(superset_profile_dict[pred_target])
										
					'''
					print "full dict"
					print superset_count_dict
					print "diff dict"
					print node_threshold_profile_difference_dict
					
					print "full profile dict"
					print superset_profile_dict
					print 'profile ratio dict'
					print target_node_ratio_dict
					print 'profile difference dict'
					print diff_ratio_dict
					'''
					threshold_profile_result_dict[current_threshold] = diff_ratio_dict
				
				profile_result_dict[test_node_name] = threshold_profile_result_dict
				
	# Done. I think
	return profile_result_dict
	
#-------------------------------------------------- Plotting / Visualisation functions

def RNA_s_structure_plot(input_file, outfile_name, bind_start, bind_stop):
	"draw a plot of the secondary structure of a RNA molecule given the Vienna notation. Requires the path to the VARNA program to be given"
	# use -highlightRegion "10-16:fill=#FF0000" to highlight binding site
	high_region = str(bind_start) + '-' + str(bind_stop) + ':fill=#FF0000'
	path_to_VARNA = '/Users/Admin/Dropbox/Programs/tools/VARNAv3-91.jar'
	result = call(['java', '-cp', path_to_VARNA, 'fr.orsay.lri.varna.applications.VARNAcmd', '-i', input_file, '-o', outfile_name, '-highlightRegion', high_region, '-algorithm', 'naview'])


#RNA_s_structure_plot('/Users/Admin/Work/Genomes/M_tuberculosis/S507/ncRNA/predicted_structures/507_ncRNA_16.vienna', 'ncRNA_507_16.png', '1', '29')

def plot_on_genome(genome, genes_list_o_dict):
	# plot the given data from the gff object to the genome
	
	# Labels and looks
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	
	plt.xlabel('Genome location')
	plt.ylabel('Value')
	plt.title('Genome plot')
	
	# parsing input data
	for gene in genes_list_o_dict:
		plt.plot([gene['start'], gene['stop']], [0,0], [1,1], linewidth=float(gene['attribute'][3:6])*10)
	
	# Actual plotting
	plt.plot([0,len(genome)],[0,0], [0,0], linewidth=4)
	plt.ylim([-10,20])
	plt.xlim([-10,len(genome)+10])
	plt.plot([1,1,3,1],[2,4,6,2],[6,8,9,6])
	# Boom
	plt.show()
	
#A_gene_dict = {'start':10, 'stop':50, 'name':'protease'}
#B_gene_dict = {'start':70, 'stop':90, 'name':'protease'}
#C_gene_dict = {'start':110, 'stop':150, 'name':'protease', 'expression':20}
#genes = [A_gene_dict, B_gene_dict, C_gene_dict]

#plot_on_genome(seq_temp, heythere)

def plot_dict_bar(aDictionary):
	"""Plot a bar chart based on a dictionary where the keys are the axis label and the values are, well, values."""
	# Convert strings to float
	for key in aDictionary:
		aDictionary[key] = float(aDictionary[key])
		
	# Plot the result
	plt.bar(range(len(aDictionary)), aDictionary.values(), align='center')
	plt.xticks(range(len(aDictionary)), aDictionary.keys(), rotation=90)
	
	plt.show()

#-------------------------------------------------- Useful dictionaries

COG_funct_dict = {
"A":"RNA processing and modification",
"B":"Chromatin Structure and dynamics",
"C":"Energy production and conversion",
"D":"Cell cycle control and mitosis",
"E":"Amino Acid metabolis and transport",
"F":"Nucleotide metabolism and transport",
"G":"Carbohydrate metabolism and transport",
"H":"Coenzyme metabolysm",
"I":"Lipid metabolysm",
"J":"Tranlsation",
"K":"Transcription",
"L":"Replication and repair",
"M":"Cell wall/membrane/envelope biogenesis",
"N":"Cell motility",
"O":"Post-translational modification, protein turnover, chaperone functions",
"P":"Inorganic ion transport and metabolism",
"Q":"Secondary Structure",
"T":"Signal Transduction",
"U":"Intracellular trafficing and secretion",
"Y":"Nuclear structure",
"Z":"Cytoskeleton",
"R":"General Functional Prediction only",
"S":"Function Unknown"
}

species_identifier = {
"Mycobacterium tuberculosis":"1773"
}		

print 'Genes loaded'