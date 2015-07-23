#!/Users/panix/anaconda/bin/python

# Importing modules
import sys
sys.path.append('/Users/panix/Dropbox/Programs/tools/Cell/Cell_core/')

from rna import *
from genes import *
import sqlite3
# Linking to database

# Variables required at startup/login
# Variables for user and database name

print "Welcome to Cell v0.1"

# Make into a dictionary to be passed around


# user_name = raw_input("User Name: ") # Jon for now
user_name = 'Jon'

conn = sqlite3.connect(user_name + ".db")


# Load user info

# This wont do, change/ move/ intergrate

session_info = {'user_name':user_name}



# Creating the database - See where and if this needs to be done and how many times

# create_user_table(session_info)


#print create_table(session_info)


# Functions


# Main loop of program


# TESTING
main_loop = True

while (main_loop == True):
	
	print 'This is the main loop'
	
	print "Menu"
	print "1: Run test"
	print "2: Import data"
	print "3: NA"
	print "4: Run Pipeline"
	print '5: Database Overview'
	print 'o: Open Script'
	print "Q: Quit"
	option = raw_input("Please select an option: ")
	print '\n'
	
	if option == '1':
		# Start with a list of genes
		test_gene_list = ['Rv0397', 'Rv0398c', 'Rv0399c', 'Rv0400c', 'Rv0401', 'Rv0402c', 'Rv0403c', 'Rv0404', 'Rv0405', 'Rv0406c', 'Rv0407', 'Rv0408', 'Rv0409', 'Rv0410c', 'Rv0411c', 'Rv0412c', 'Rv0980c', 'Rv0981', 'Rv0982', 'Rv0983', 'Rv0984', 'Rv0985c', 'Rv0986', 'Rv0987', 'Rv0988', 'Rv0989c', 'Rv0990c', 'Rv0991c', 'Rv0992c', 'Rv0993', 'Rv0994', 'Rv0995', 'Rv0996', 'Rv0997', 'Rv0998', 'Rv1308', 'Rv1309', 'Rv1310', 'Rv1311', 'Rv1312', 'Rv1313c', 'Rv1314c', 'Rv1315', 'Rv1316c', 'Rv1317c', 'Rv1318c', 'Rv1319c', 'Rv1320c', 'Rv1321', 'Rv1322', 'Rv1322A', 'Rv1323', 'Rv1324', 'Rv1325c', 'Rv1493', 'Rv1494', 'Rv1495', 'Rv1496', 'Rv1497', 'Rv1498A', 'Rv1498c', 'Rv1499', 'Rv1500', 'Rv1501', 'Rv1502', 'Rv1503c', 'Rv1504c', 'Rv1505c', 'Rv1506c', 'Rv1507A', 'Rv1507c', 'Rv1508A', 'Rv1508c', 'Rv1509', 'Rv1510', 'Rv1511', 'Rv1512', 'Rv1513', 'Rv1514c', 'Rv1515c', 'Rv1516c', 'Rv1658', 'Rv1659', 'Rv1660', 'Rv1661', 'Rv1662', 'Rv1663', 'Rv1664', 'Rv1665', 'Rv1666c', 'Rv1667c', 'Rv1668c', 'Rv2484c', 'Rv2485c', 'Rv2486', 'Rv2487c', 'Rv2488c', 'Rv2489c', 'Rv2490c', 'Rv2491', 'Rv2492', 'Rv2493', 'Rv2494', 'Rv2495c', 'Rv2496c', 'Rv2497c', 'Rv2498c', 'Rv2499c', 'Rv2500c', 'Rv2812', 'Rv2813', 'Rv2814c', 'Rv2815c', 'Rv2816c', 'Rv2817c', 'Rv2818c', 'Rv2819c', 'Rv2820c', 'Rv2821c', 'Rv2822c', 'Rv2823c', 'Rv2824c', 'Rv2825c', 'Rv2826c', 'Rv2827c', 'Rv2828c', 'Rv2942', 'Rv2943', 'Rv2943A', 'Rv2944', 'Rv2945c', 'Rv2946c', 'Rv2947c', 'Rv2948c', 'Rv2949c', 'Rv2950c', 'Rv2951c', 'Rv2952', 'Rv2953', 'Rv2954c', 'Rv2955c', 'Rv3106', 'Rv3107c', 'Rv3108', 'Rv3109', 'Rv3110', 'Rv3111', 'Rv3112', 'Rv3113', 'Rv3114', 'Rv3115', 'Rv3116', 'Rv3117', 'Rv3118', 'Rv3119', 'Rv3120', 'Rv3121', 'Rv3122', 'Rv3123', 'Rv3124', 'Rv3125c', 'Rv3126c', 'Rv3127']
		
		results_list = get_homologues("Jon", test_gene_list, "H37Rv", "F11")
		print results_list
		for entry in results_list:
			print entry[0]['LOCUS'], entry[0]['HOMOLOGUES'][0]
		print len(test_gene_list)
		print len(results_list)
		
		unique_genes = []
		for gene in results_list:
			if len(gene[0]['HOMOLOGUES'][0]) == 0:
				unique_genes.append(gene[0]['LOCUS'])
		print unique_genes
		
		# Convert gene string to dict
		unique_gene_dict = []
		for gene in unique_genes:
			gene_dict = get_db_subset("Jon", "H37Rv", "Gene", "Locus", gene)
			unique_gene_dict.append(gene_dict[0])
					
		print unique_gene_dict
		
		print cluster_on_funct(unique_gene_dict, "GENOME_ONTOLOGY")
		
	if option == 'o':
		# Load external script
		print 'Available scripts:'
		available_operons = subprocess.Popen(['ls', './operons/'], stdout=subprocess.PIPE).communicate()[0]
		available_operons = available_operons.split()
		i = 0
		for script in available_operons:
			i = i + 1
			print i, ':', script
			
		selected_script = raw_input("Select a script \n")
		
		script_result = subprocess.Popen(['./operons/' + available_operons[int(selected_script)-1]], stdout=subprocess.PIPE).communicate()[0]
		print available_operons[int(selected_script)-1] + " loaded"
		print ""
		print script_result
		
		
	
	if option == 't':
		# Quick test option
		create_blast_db(session_info['user_name'], 'C')
		create_blast_db(session_info['user_name'], 'H37Ra')
		create_blast_db(session_info['user_name'], 'RGTB327')
	
	if option == '6':
		# Pipeline for the analysis of ncRNA
		
		# Thresholds and parameters
		energy_threshold = -20.0
		zscore_threshold = 1
		blast_threshold = '0.1'
		session_info = 'Jon'
		
		# Isolate list to work with:
		# These isolates must have blast databases available, genome sequences and gene summaries loaded (last one may be optional)
		
		target_isolates = ['CDC1551', 'F11', 'H37Rv']
		
		# The reference isolate is the isolate to which all alignments take place, 
		# except later where binding sites are aligned directly
		
		reference_isolate = 'H37Rv'
		
		
		# Input files containing the gene expression data 
		
		GDS_Filepath = '~/Work/Genomes/M_tuberculosis/Expression_data/GDS1552/GDS1552_full.soft'
		GPL_Filepath = '/Users/Admin/Work/Genomes/M_tuberculosis/Expression_data/GDS1552/GPL2787.annot'
		
		# Output file paths
		
		out_dir = "/Users/Admin/Work/Genomes/M_tuberculosis/Expression_data/analysis_results/"
		
		from ncRNA_data_import import *
		
		# Temp hard coded dict of where samples originate from:
		sample_origins = {}
		
		sample_origins['CDC1551'] = ['GSM71949','GSM71953','GSM71957','GSM71984']
		
		sample_origins['H37Rv'] = ['GSM71958','GSM71988','GSM71989','GSM71990']
		
		# THIS IS A LIE because there is no F11 data
		sample_origins['F11'] = ['GSM71963','GSM71964','GSM71968','GSM71976']

		# This is where the data from the target prediction files is imported and stored as a large dict. 
		ncRNA_target_data = ncRNA_import()
		
		for ncRNA in ncRNA_target_data:
			# This loop moves over one ncRNA at a time, and in the end produces a csv file of the results.
			
			# The following dictionary is a way of ro-organising the data. 
			ncRNA_results_dict = {}
			
			# Here we make sure that we are working with non-empty entries
			if 'start' in ncRNA and len(ncRNA['start'] ) > 0:
				
				# Base information
				ncRNA_results_dict['ncRNA'] = ncRNA['number']
				
				# Get the co-ordinates of the ncRNA
				print ncRNA['number']
				print ncRNA['start'] + '-' + ncRNA['stop']
				
				# Get the genomic sequence from the reference for the ncRNA (NB coordinates are for a specific isolate!!)
				ncRNA_seq = get_seq_from_coords(session_info, reference_isolate, ncRNA['start'], ncRNA['stop'])
				print ncRNA_seq
				

				# Blast the sequence against the other genomes and return the score
				sim_dict = {}
				for isolate in target_isolates:
					blast_result = conduct_blast_search(session_info, isolate, ncRNA_seq, '0.1')
					print isolate
					print blast_result[0]

					sim_dict[isolate] = blast_result[0][2]
				ncRNA_results_dict['strain'] = sim_dict
					
				# Get the targets of the ncRNA that have a binding energy above the threshold and add those to a list (Just names)
				# First for RNAPredator
				# Additionally we are going to try use this info to get the ncRNA binding site for this ncRNA-gene
				# ---------------------------------
				
				ncRNA_pred_targets = []
				print ""
				
				for Pred_target in ncRNA['RNAPredator_predictions']:
					if float(Pred_target['Energy']) < energy_threshold:
						print Pred_target['Gene'], Pred_target['Annotation']
						ncRNA_pred_targets.append(Pred_target['Gene'])
						
						#print "ncRNA seq: ", ncRNA_seq
						print "ncRNA binding site"
						print ncRNA_seq[int(Pred_target['sRNA_start']):int(Pred_target['sRNA_stop'])]
						print Pred_target['Interaction']
						
						print "Current ncRNA target gene: ", Pred_target['Gene']
						
						# Now, we need to get the sequence of the gene from the reference genome
						reference_ncRNA_target_seq = get_seq_from_coords(session_info, reference_isolate, get_gene_item(session_info, reference_isolate, Pred_target['Gene'], "START"), get_gene_item(session_info, reference_isolate, Pred_target['Gene'], "STOP"))
						# print reference_ncRNA_target_seq
						print "Reference sequence binding_site"
						if Pred_target['Co-ordinates'][0:1] == 'c':
							# This compensates for binding in the opposite direction. And an indexing error that requires a +1
							ref_target_binding_site =  reference_ncRNA_target_seq[(-1)*int(Pred_target['mRNA_stop'])+1:(-1)*int(Pred_target['mRNA_start'])+1]
						else:	
							ref_target_binding_site = reference_ncRNA_target_seq[int(Pred_target['mRNA_start']):int(Pred_target['mRNA_stop'])]
						
						print ref_target_binding_site	
						
						print "---------------------------------------------"
						# Then align that sequence to the query genomes 
						for isolate in target_isolates:
							print ""
							homologous_target_blast = conduct_blast_search(session_info, isolate, reference_ncRNA_target_seq, blast_threshold)
							
							# Then get the sequence from the gene it aligns too
							homologous_target_seq = get_seq_from_coords(session_info, isolate, homologous_target_blast[0][8], homologous_target_blast[0][9])
							# print homologous_target_seq
							print "Target binding site"
							if Pred_target['Co-ordinates'][0:1] == 'c':
								homologous_tar_binding_site =  homologous_target_seq[(-1)*int(Pred_target['mRNA_stop'])+1:(-1)*int(Pred_target['mRNA_start'])+1]
							else:	
								homologous_tar_binding_site =  homologous_target_seq[int(Pred_target['mRNA_start']):int(Pred_target['mRNA_stop'])]
							print homologous_tar_binding_site
							print "Reference - Homologue simmilarity"
							print homologous_target_blast[0][2]
							print "Simmilarity of the reference target bindong site to the homologue target binding site"
							# Comparing length doesn't work because the homology based search can return incorrect co-ordinates...
							
							print len(reference_ncRNA_target_seq), " vs ", len(homologous_target_seq)
							print similar(ref_target_binding_site, homologous_tar_binding_site)
							
							# Then align the ncRNA binding site sequence to the target gene from the query genome
							
							
						# return the similarity of the two
						# possible return the binding energy of the alignment
				
				# Get the expression data for the targets

				expression_data_dict = {}
				for gene in ncRNA_pred_targets:
					gene_expression_data = get_expression_data(GDS_Filepath, GPL_Filepath, gene)
					print "Gene = " + gene
					print gene_expression_data
					
					expression_data_dict[gene] = gene_expression_data
				
				ncRNA_results_dict['target_expression'] = expression_data_dict
				
				# Get the sequence similarity between ncRNA targets
				
				# Label as a dict
				label_dict = {}
				
				for gene in ncRNA_pred_targets:
					# Get the sequence of the gene
					
					pred_gene_seq = get_seq_from_coords(session_info, reference_isolate, get_gene_item("Jon", reference_isolate, gene, 'START'), get_gene_item("Jon", reference_isolate, gene, 'STOP'))

					
					# Conduct blast search from reference isolate to other isolates
					target_homologue_dict = {}
					for isolate in target_isolates:
						isolate_ncRNA_target_homologue = conduct_blast_search(session_info, isolate, pred_gene_seq, blast_threshold)
						if len(isolate_ncRNA_target_homologue) == 0:
							# For genes with no homologues
							isolate_ncRNA_target_homologue = [['NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA']]
							
						print isolate_ncRNA_target_homologue[0]
						
						target_homologue_dict[isolate] = isolate_ncRNA_target_homologue[0][2]
				
					label_dict[gene] = target_homologue_dict
				
				ncRNA_results_dict['target_homology'] = label_dict
				
				print "##################------ ncRNA Complete ------##################"
				print "ncRNA: " + ncRNA_results_dict['ncRNA']
				print ncRNA_results_dict['strain']
				print "----------------------------------------------------------------"
				print ncRNA_results_dict['target_homology']
				print "----------------------------------------------------------------"
				print ncRNA_results_dict['target_expression']
				
				# Now to format and analyse the data 
				
				# Concept test
				
				print ncRNA_results_dict['ncRNA']
				print ncRNA_results_dict['strain']['CDC1551']
				#print ncRNA_results_dict['target_homology']['Rv0746']['CDC1551']
				print ncRNA_results_dict['target_expression']
				
				# For tomorrow, gets the isolate where the GSM sample came from
				#print get_GSM_isolate(GDS_Filepath, GPL_Filepath, 'GSM71970')
				
				print "Fun part"
				# With exporting to a csv file
				
				out_csv_file = open(out_dir + "ncRNA_"+ncRNA_results_dict['ncRNA'] + ".csv" , 'w')
				
				for input_isolate in target_isolates:
					# Get per isolate stats for this ncRNA
					print "# Results for" + input_isolate
					
					for gene in ncRNA_results_dict['target_homology']:
						print gene
						expression_values = []
					
						for isolate in ncRNA_results_dict['target_expression'][gene]:
							if isolate in sample_origins[input_isolate]:
								print isolate, ncRNA_results_dict['target_expression'][gene][isolate]
							
								# For calculating the mean
								if ncRNA_results_dict['target_expression'][gene][isolate] != 'NA':
									expression_values.append(float(ncRNA_results_dict['target_expression'][gene][isolate]))
						
						if len(expression_values) > 0:
							the_mean = str(sum(expression_values) / float(len(expression_values)))
							print "Mean = " + str(sum(expression_values) / float(len(expression_values)))
							the_variance = str(numpy.var(expression_values))
							print "Variance = " + str(numpy.var(expression_values))
						else:
							print "Mean = NA"
							the_mean = 'NA'
							print "Variance = NA"
							the_variance = 'NA'
							
						# Here we add the binding site data
						print "=================The adding the binding site part ========================"
						print int(ncRNA_results_dict['ncRNA'])

						print ncRNA_target_data[int(ncRNA_results_dict['ncRNA'])]['RNAPredator_predictions'][0]['sRNA_start']
						print ncRNA_target_data[int(ncRNA_results_dict['ncRNA'])]['RNAPredator_predictions'][0]['Gene']
						print ncRNA_target_data[int(ncRNA_results_dict['ncRNA'])]['RNAPredator_predictions'][1]
						
						# This part correctly returns the binding details for gene X of the current ncRNA
						print gene
						print '>'
						for this_ncRNA_target in ncRNA_target_data[int(ncRNA_results_dict['ncRNA'])]['RNAPredator_predictions']:
							if this_ncRNA_target['Gene'] == gene:
								print this_ncRNA_target
							
				
						out_string = input_isolate + ',' + gene + ',' + ncRNA_results_dict['strain'][input_isolate] + ',' + ncRNA_results_dict['target_homology'][gene][input_isolate] + ',' + str(the_mean) + ',' + str(the_variance) + '\n'
						out_csv_file.write(out_string)
				
				out_csv_file.close()	
				
				
	
	if option == '7':
		#add_to_biological_database(session_info, 'F11', 'Genome', '/Users/Admin/Work/Genomes/M_tuberculosis/F11/mycobacterium_tuberculosis_f11__finished__4_supercontigs.fasta')
		h37Rv_genes = []
		h37Rv_genes = h37Rv_genes + get_annotations_from_locus_range(session_info["user_name"], 'TB', '476000', '499599')
		h37Rv_genes = h37Rv_genes + get_annotations_from_locus_range(session_info["user_name"], 'TB', '1095500', '1115099')
		h37Rv_genes = h37Rv_genes + get_annotations_from_locus_range(session_info["user_name"], 'TB', '1464500', '1488599')
		h37Rv_genes = h37Rv_genes + get_annotations_from_locus_range(session_info["user_name"], 'TB', '1685500', '1708099')
		h37Rv_genes = h37Rv_genes + get_annotations_from_locus_range(session_info["user_name"], 'TB', '1871500', '1895099')
		h37Rv_genes = h37Rv_genes + get_annotations_from_locus_range(session_info["user_name"], 'TB', '2792000', '2814099')
		h37Rv_genes = h37Rv_genes + get_annotations_from_locus_range(session_info["user_name"], 'TB', '3117000', '3136099')
		h37Rv_genes = h37Rv_genes + get_annotations_from_locus_range(session_info["user_name"], 'TB', '3286000', '3308099')
		h37Rv_genes = h37Rv_genes + get_annotations_from_locus_range(session_info["user_name"], 'TB', '3474000', '3492599')
		print 'h37Rv_genes' , unique_list(h37Rv_genes)
		print 'length = ' , len(unique_list(h37Rv_genes))
		
		# Get the homologues of the genes
		H37Rv_genes_F11_homologues = get_homologues("Jon", unique_list(h37Rv_genes), "H37Rv", "F11")
		H37Rv_genes_CDC1551_homologues = get_homologues("Jon", unique_list(h37Rv_genes), "H37Rv", "CDC1551")
		
		F11_genes = []
		F11_genes = F11_genes + get_annotations_from_locus_range(session_info["user_name"], 'F11', '479500', '499099')
		F11_genes = F11_genes + get_annotations_from_locus_range(session_info["user_name"], 'F11', '830000', '852099')
		F11_genes = F11_genes + get_annotations_from_locus_range(session_info["user_name"], 'F11', '1099000', '1119599')
		F11_genes = F11_genes + get_annotations_from_locus_range(session_info["user_name"], 'F11', '1468500', '1493099')
		F11_genes = F11_genes + get_annotations_from_locus_range(session_info["user_name"], 'F11', '1630000', '1652099')
		F11_genes = F11_genes + get_annotations_from_locus_range(session_info["user_name"], 'F11', '1690000', '1712099')
		F11_genes = F11_genes + get_annotations_from_locus_range(session_info["user_name"], 'F11', '1867000', '1890599')
		F11_genes = F11_genes + get_annotations_from_locus_range(session_info["user_name"], 'F11', '2308000', '2330099')
		F11_genes = F11_genes + get_annotations_from_locus_range(session_info["user_name"], 'F11', '2806000', '2828099')
		F11_genes = F11_genes + get_annotations_from_locus_range(session_info["user_name"], 'F11', '3128500', '3148099')
		F11_genes = F11_genes + get_annotations_from_locus_range(session_info["user_name"], 'F11', '3298000', '3320099')
		F11_genes = F11_genes + get_annotations_from_locus_range(session_info["user_name"], 'F11', '3486000', '3504099')
		F11_genes = F11_genes + get_annotations_from_locus_range(session_info["user_name"], 'F11', '3796000', '3818099')
		print 'F11_genes' , unique_list(F11_genes)
		print 'length = ' , len(unique_list(F11_genes))
		
		# Get the homologues of the genes
		F11_genes_H37Rv_homologues = get_homologues("Jon", unique_list(F11_genes), "F11", "H37Rv")
		F11_genes_CDC1551_homologues = get_homologues("Jon", unique_list(F11_genes), "F11", "CDC1551")
		
		CDC1551_genes = []
		CDC1551_genes = CDC1551_genes + get_annotations_from_locus_range(session_info["user_name"], 'CDC1551', '474000', '496099')
		CDC1551_genes = CDC1551_genes + get_annotations_from_locus_range(session_info["user_name"], 'CDC1551', '831500', '849099')
		CDC1551_genes = CDC1551_genes + get_annotations_from_locus_range(session_info["user_name"], 'CDC1551', '1095500', '1115099')
		CDC1551_genes = CDC1551_genes + get_annotations_from_locus_range(session_info["user_name"], 'CDC1551', '1464000', '1488599')
		CDC1551_genes = CDC1551_genes + get_annotations_from_locus_range(session_info["user_name"], 'CDC1551', '1626000', '1648099')
		CDC1551_genes = CDC1551_genes + get_annotations_from_locus_range(session_info["user_name"], 'CDC1551', '1685500', '1708599')
		CDC1551_genes = CDC1551_genes + get_annotations_from_locus_range(session_info["user_name"], 'CDC1551', '1862000', '1886099')
		CDC1551_genes = CDC1551_genes + get_annotations_from_locus_range(session_info["user_name"], 'CDC1551', '2788000', '2810099')
		CDC1551_genes = CDC1551_genes + get_annotations_from_locus_range(session_info["user_name"], 'CDC1551', '3111000', '3129599')
		CDC1551_genes = CDC1551_genes + get_annotations_from_locus_range(session_info["user_name"], 'CDC1551', '3285000', '3302599')
		CDC1551_genes = CDC1551_genes + get_annotations_from_locus_range(session_info["user_name"], 'CDC1551', '3466000', '3488099')
		CDC1551_genes = CDC1551_genes + get_annotations_from_locus_range(session_info["user_name"], 'CDC1551', '3776000', '3798099')
		print 'CDC1551_genes' , unique_list(CDC1551_genes)
		print 'length = ' , len(unique_list(CDC1551_genes))
		
		# Get the homologues of the genes
		CDC1551_genes_H37Rv_homologues = get_homologues("Jon", unique_list(CDC1551_genes), "CDC1551", "H37Rv")
		CDC1551_genes_F11_homologues = get_homologues("Jon", unique_list(CDC1551_genes), "CDC1551", "F11")
		
		print CDC1551_genes_H37Rv_homologues
		print H37Rv_genes_CDC1551_homologues
		print H37Rv_genes_F11_homologues
		
		genome_list = [H37Rv_genes_F11_homologues, H37Rv_genes_CDC1551_homologues, F11_genes_H37Rv_homologues, F11_genes_CDC1551_homologues, CDC1551_genes_H37Rv_homologues, CDC1551_genes_F11_homologues]
		
		
		G = create_homologue_graph(genome_list) 
		
		edge_list_H37Rv_genes_F11 = get_edge_lists(G, H37Rv_genes_F11_homologues)
		edge_list_H37Rv_genes_CDC1551 = get_edge_lists(G, H37Rv_genes_CDC1551_homologues)
		
		edge_list_F11_genes_H37Rv = get_edge_lists(G, F11_genes_H37Rv_homologues)
		edge_list_F11_genes_CDC1551 = get_edge_lists(G, F11_genes_CDC1551_homologues)
		
		edge_list_CDC1551_genes_H37Rv = get_edge_lists(G, CDC1551_genes_H37Rv_homologues)
		edge_list_CDC1551_genes_F11 = get_edge_lists(G, CDC1551_genes_F11_homologues)
		
		print 'H37Rv'
		print edge_list_H37Rv_genes_F11['0']
		print edge_list_H37Rv_genes_CDC1551['0']
		
		print 'F11'
		print edge_list_F11_genes_H37Rv['0']
		print edge_list_F11_genes_CDC1551['0']
		
		print 'CDC1551'
		print edge_list_CDC1551_genes_H37Rv['0']
		print edge_list_CDC1551_genes_F11['0']
		
		nx.draw(G)
		plt.show()

	if option == '2':
		print "Current DB Status:"
		print db_status(session_info, 'TB')

		# Input data options
		print 'What would you like to import?'
		print '1: Gene data'
		print '2: VCF data'
		print '3: Back'
		print '4: New Organism'
		print '5: TF Data'
		
		import_option = raw_input('Selection: ')
		
		print 'Adding the extra stuff \n \n '
		# Do this only once 
		
		if import_option == '1':
			database_table = 'Gene'
			organism = raw_input("For which organism?: ")
			gene_filepath = raw_input("Gene summary filepath: ")
			parsed_genes = input_parser(gene_filepath)
			add_to_biological_database(session_info, organism, database_table, parsed_genes)
			print '\n Import complete'
			
		if import_option == '2':
			gene_filepath = raw_input("VCF summary filepath: ")
			organism = raw_input("Organism: ")
			database_table = raw_input("Reference organism: ")
			parsed_vcf = input_parser(gene_filepath)
			add_to_biological_database(session_info, organism, database_table, parsed_vcf)
			print '\n Import complete'
		
		if import_option == '3':
			print "Main:"

		if import_option == '4':
			#removed requirement for genome filepath
			new_organism = {}
			organism_descriptor = raw_input('Organism descriptor: ')
			new_organism['Organism'] = organism_descriptor
			tax_ID = raw_input('Taxonomic identifier: ')
			new_organism['Taxa_ID'] = tax_ID
			heir_data = raw_input('Taxomonic data: ')
			new_organism['Taxa_Heirarchy'] = heir_data
			
			# New db and tables for organism
			create_table(session_info, new_organism['Organism'])
			
			add_genome_opt = raw_input("Add genome sequence? (y/n)")
			if add_genome_opt == "y" or add_genome_opt == "Y":
				genome_path = raw_input('Genome file path: ')
				new_organism['genome_path'] = genome_path
				# Populate with info
				add_to_user_database(session_info, 'Organism', new_organism)
				# Add genome path to Genome table
				add_to_biological_database(session_info, new_organism['Organism'], 'Genome', new_organism['genome_path'])
				# make blast db
				create_blast_db(session_info['user_name'], organism_descriptor)
			else:
				# Populate with info
				add_to_user_database(session_info, 'Organism', new_organism)
				# Add genome path to Genome table
				new_organism['genome_path'] = "None"
				add_to_biological_database(session_info, new_organism['Organism'], 'Genome', new_organism['genome_path'])
		
		if import_option == '5':
			# Add TF data
			database_table = 'TF_targets'
			organism = raw_input("For which organism?: ")
			gene_filepath = raw_input("TF interactions filepath: ")
			gene_filepath = "/Users/Admin/Work/Genomes/M_tuberculosis/Interaction_data/interactions.txt"
			parsed_genes = input_parser(gene_filepath)
			add_to_biological_database(session_info, organism, database_table, parsed_genes)


	if option == '3':	
		print 'open'
		print len(get_SNPs_vcf('/Users/Admin/Work/Genomes/M_tuberculosis/S507/snps/507_High_Conf.vcf', '2000', 10000))
		print get_SNPs_vcf('/Users/Admin/Work/Genomes/M_tuberculosis/S507/snps/507_High_Conf.vcf', '2000', 10000)
	
	if option == '4':
		menu_pipeline = True
		print "Pipeline options:"
		print 'test bed'
		# From CDC1551
		test_genes = ['MT3205', 'MT3208']
		print test_genes
		
		print 'Correct retireval of test gene coordinates?'
		test_gene_full = get_db_subset("Jon", "CDC1551", "Gene", "Locus", test_genes[0])
		
		print 'Observed: ' + test_gene_full[0]['START'] + "-" + test_gene_full[0]['STOP']
		print 'Expected: 3483886-3484356'
		print 'Pass'
		
		print "Correct sequence retrieved?"
		
		gene_x_seq = get_seq_from_coords("Jon", "CDC1551", test_gene_full[0]['START'], test_gene_full[0]['STOP'])
		print gene_x_seq
		
		CDC1551_genes_H37Rv_homologues = get_homologues("Jon", unique_list(test_genes), "CDC1551", "H37Rv")
		print CDC1551_genes_H37Rv_homologues[1][0]
		
	if option == '9':
		# locational snp reporting
		print "ncRNA simmilarity"
		print conduct_blast_search("Jon", "H37Rv", "ACGACCCCCGCCAGGGGGGGAGGAGGCGAGGGTCGTCGTGCATCAGCCCCGGGGGGTCGGACTGATACACCCTCGGCTATGGCCGAGTAATGCT", "0.0001")
		print conduct_blast_search("Jon", "F11", "ACGACCCCCGCCAGGGGGGGAGGAGGCGAGGGTCGTCGTGCATCAGCCCCGGGGGGTCGGACTGATACACCCTCGGCTATGGCCGAGTAATGCT", "0.0001")
		print conduct_blast_search("Jon", "CDC1551", "ACGACCCCCGCCAGGGGGGGAGGAGGCGAGGGTCGTCGTGCATCAGCCCCGGGGGGTCGGACTGATACACCCTCGGCTATGGCCGAGTAATGCT", "0.0001")
		
		print "Target simmilarity"
		#print get_SNPs("Jon", "H37Rv", "Rv2672")
		
		print get_gene_item("Jon", "H37Rv", "Rv0204c", "START")
		print get_gene_item("Jon", "H37Rv", "Rv0204c", "STOP")
		print get_seq_from_coords("Jon", "H37Rv", get_gene_item("Jon", "H37Rv", "Rv0204c", "START"), get_gene_item("Jon", "H37Rv", "Rv0204c", "STOP"))
		print conduct_blast_search("Jon", "H37Rv", get_seq_from_coords("Jon", "H37Rv", get_gene_item("Jon", "H37Rv", "Rv0204c", "START"), get_gene_item("Jon", "H37Rv", "Rv0204c", "STOP")), "0.1")
		print conduct_blast_search("Jon", "F11", get_seq_from_coords("Jon", "H37Rv", get_gene_item("Jon", "H37Rv", "Rv0204c", "START"), get_gene_item("Jon", "H37Rv", "Rv0204c", "STOP")), "0.1")
		print conduct_blast_search("Jon", "CDC1551", get_seq_from_coords("Jon", "H37Rv", get_gene_item("Jon", "H37Rv", "Rv0204c", "START"), get_gene_item("Jon", "H37Rv", "Rv0204c", "STOP")), "0.1")
		
		# Get sequences for alignment
		print get_seq_from_coords("Jon", "F11", '242294', '243532')
		print get_seq_from_coords("Jon", "H37Rv", '241976', '243214')
		print get_seq_from_coords("Jon", "CDC1551", '242090', '243328')
		
		print "for kick blast"
		print conduct_blast_search("Jon", "H37Rv", "AGGCACACCGGTACACATGGGCAGACCCGGCGTGACTCTCGGGGGGCGTCTGACACCGCCTTCTGCGGGTCTTGCGCGGCCGGCCTTCACCCCGTCTTCCGGCACTTTCGATTGGTCACTAACCGGGCCTGC", "1")
		print conduct_blast_search("Jon", "F11", "AGGCACACCGGTACACATGGGCAGACCCGGCGTGACTCTCGGGGGGCGTCTGACACCGCCTTCTGCGGGTCTTGCGCGGCCGGCCTTCACCCCGTCTTCCGGCACTTTCGATTGGTCACTAACCGGGCCTGC", "1")
		print conduct_blast_search("Jon", "CDC1551", "AGGCACACCGGTACACATGGGCAGACCCGGCGTGACTCTCGGGGGGCGTCTGACACCGCCTTCTGCGGGTCTTGCGCGGCCGGCCTTCACCCCGTCTTCCGGCACTTTCGATTGGTCACTAACCGGGCCTGC", "1")


	if option == '5':
		print 'DB Status: '
		db_status_user(session_info)
		print ''
		this = get_available_organisms(session_info)
		for item in this:
			print item[0] + ' status:'
			db_status(session_info, item[0])
			print ''
		# Add a for each in user table feature
				
	if option == '8':
		print "ncRNA analysis pipeline and tools"
		ncRNA_path = '/Users/Admin/Work/Genomes/M_tuberculosis/Paida_ncRNA_Data/Final_Mtb_ncRNA/ncRNA-Table 1.csv'
		# Retrieve the ncRNA seq from the Paida's file.
		list_of_dicts  = []
		reader = csv.DictReader(open(ncRNA_path, 'r'), delimiter=',')
		for row in reader:
			list_of_dicts.append(row)
		print list_of_dicts
		for ncRNA in list_of_dicts:
			if len(ncRNA['NUMBER']) > 0:
				print ncRNA['CHROMOSOMAL LOCATION (kb)'], ncRNA['GENE LENGTH ']
				# Convert to ACTUAL sequence coordinates
				ncRNA_start = ncRNA['CHROMOSOMAL LOCATION (kb)'].replace(',', '')
				print ncRNA_start
				print ncRNA['NUMBER']
				print float(ncRNA_start)*1000
				ncRNA_start = int(float(ncRNA_start)*1000)
				
				# Now for the stop
				if float(ncRNA_start) > (float(ncRNA['GENE LENGTH '])*1000):
					ncRNA_end = float(ncRNA_start) + (float(ncRNA['GENE LENGTH ']))
					ncRNA_end = int(ncRNA_end)
					print 'greater'
				else:
					print 'less than'
					ncRNA_end = int(float(ncRNA['GENE LENGTH '])*1000)
					print 'len', ncRNA_start - ncRNA_end
				print ncRNA_start , '-' , ncRNA_end
				genes_there = get_annotations_from_locus_range("Jon", 'H37Rv', str(ncRNA_start), str(ncRNA_end))
				print genes_there
				ncRNA_seq = get_seq_from_coords("Jon", "H37Rv", str(ncRNA_start), str(ncRNA_end))
				print ncRNA_seq
		# Export fasta file using new co-ords
		for ncRNA in list_of_dicts:
			if len(ncRNA['START']) > 0:
				print ">ncRNA" + ncRNA['NUMBER']
				# print ncRNA['START'], ' ', ncRNA['STOP']
				print get_seq_from_coords("Jon", "H37Rv", str(ncRNA['START']), str(ncRNA['STOP']))
				print '\n'
				
	if option == 'Q' or option == 'q':
		main_loop = False
		print "Leaving Cell"
		
	print '\n'

# To Do
# VCF Parser and interpreter - load VCF into db.

# DB for available commands so the main loop can query and run them


