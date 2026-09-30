import os
import sys
import pandas as pd
import zipfile
import numpy as np
import gzip
from Bio.Seq import Seq
from scipy.cluster import hierarchy
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
plt.style.use('classic')
plt.rcParams['svg.fonttype'] = 'none' 
import matplotlib as mpl
mpl.rc('font',family='times new roman')
file_path_1 = '/home/willow/'
#file_path_1 = '/home/wkioncro/'
#file_path_2 = '/home/willow/'
file_path_2 = '/media/willow/MyPassport/'
#file_path_2 = '/vol/projects/wkioncro/'

main_folder = file_path_2 + 'Faber/2024-07-22_Suyen_Espinoza_24121PR_raw_FASTQ/'
main_folder = file_path_2 + 'Faber/2025-04-10_Suyen_Espinoza_25029PR/'
main_folder = file_path_2 + 'Faber/2025-07-28_Christopher_Birk_25062PR_raw_FASTQ_Pool2/'
main_folder = file_path_2 + 'Faber/2025-07-28_Christopher_Birk_25062PR_raw_FASTQ_Pool1/'

# After quantification, download featureCounts and BBDuk_L_R (only fastqc files)
# from cluster. Name BBDuk_L_R as BBDuk_L_R_fastqc

# Also this file has to be created by hand for each pool: Cdiff_seq_file_conversion_table.ods

# The following three files are output from running multiqc in each of the 
# three folders indicated, and must be converted to .tsv 
# (by just changing the file extension)
if False:
	os.system('multiqc ' + main_folder + 'quality_stats_by_lib')
	os.system('mv multiqc* ' + main_folder + 'quality_stats_by_lib/')
	os.system('multiqc ' + main_folder + 'BBDuk_L_R_fastqc')
	os.system('mv multiqc* ' + main_folder + 'BBDuk_L_R_fastqc/')
	os.system('multiqc ' + main_folder)
	os.system('mv multiqc* ' + main_folder + '/')

df = pd.read_table(main_folder + 'multiqc_data/featureCounts_assignment_plot.txt')
assigned_reads = {term:df['Assigned'][i] for i,term in enumerate(df['Sample'])}
unassiged_reads = {term:df['Unassigned: Unmapped'][i] for i,term in enumerate(df['Sample'])}
print(assigned_reads)
print(unassiged_reads)
df = pd.read_table(main_folder + 'BBDuk_L_R_fastqc/multiqc_data/fastqc_sequence_counts_plot.txt')
unique_trim_reads = {term:df['Unique Reads'][i] for i,term in enumerate(df['Sample'])}
duplicate_trim_reads = {term:df['Duplicate Reads'][i] for i,term in enumerate(df['Sample'])}
print(unique_trim_reads)
print(duplicate_trim_reads)
df = pd.read_table(main_folder + 'quality_stats_by_lib/multiqc_data/fastqc_sequence_counts_plot.txt')
unique_reads = {term:df['Unique Reads'][i] for i,term in enumerate(df['Sample'])}
duplicate_reads = {term:df['Duplicate Reads'][i] for i,term in enumerate(df['Sample'])}
print(unique_reads)
print(duplicate_reads)
#print(assigned_reads['vibrio04'] + unassiged_reads['vibrio04'])
#print(assigned_reads['vibrio04'])
#print(unassiged_reads['vibrio04'])
#print(unique_trim_reads['vibrio04'], duplicate_trim_reads['vibrio04'], unique_trim_reads['vibrio04'] + duplicate_trim_reads['vibrio04'])
#print(unique_reads['vibrio04'], duplicate_reads['vibrio04'], unique_reads['vibrio04'] + duplicate_reads['vibrio04'])

# This file has to be created by hand for each pool
df = pd.read_excel(main_folder + 'Cdiff_seq_file_conversion_table.ods')
File_types = [str(term) for term in list(df["File type"]) if not str(term) == 'nan']
File_names = [str(term) for term in list(df["File name"]) if not str(term) == 'nan']

pool_type = {file_name:File_types[i] for i,file_name in enumerate(File_names)}

df = {'Sample ID':[], 'Sample':[], 'Sequenced reads':[], 'Reads after trimming':[], 'Reads trimmed':[], 'Assigned reads':[], 'Reads not assigned':[], '% assigned reads':[], '% not assigned reads':[]}
for file_name in File_names:
	if file_name in unique_reads:
		reads_after_trim = 0
		if file_name in unique_trim_reads:
			reads_after_trim += unique_trim_reads[file_name] 
		if file_name in duplicate_trim_reads:
			reads_after_trim += duplicate_trim_reads[file_name]
		seq_reads = unique_reads[file_name] + duplicate_reads[file_name]
		unassigned = reads_after_trim - assigned_reads[file_name]
		df['Sample ID'].append(file_name)
		df['Sample'].append(pool_type[file_name])
		df['Sequenced reads'].append(seq_reads)
		df['Reads after trimming'].append(reads_after_trim)
		df['Reads trimmed'].append(seq_reads - reads_after_trim)
		df['Assigned reads'].append(assigned_reads[file_name])
		df['Reads not assigned'].append(unassigned)
		df['% assigned reads'].append(round(100*assigned_reads[file_name]/seq_reads, 1))
		df['% not assigned reads'].append(round(100*unassigned/seq_reads, 1))

df = pd.DataFrame(df)
df.to_excel('Cdiff_read_stats.xlsx') 

# This file should be the same for all pools
df = pd.read_excel('Cdiff_gene_info.xlsx')
gene_IDs = [str(term) for term in list(df["gene ID"]) if not str(term) == 'nan']
locus_tags = [str(term) for term in list(df["locus tag"]) if not str(term) == 'nan']
gene_names = [str(term) for term in list(df["Name"]) if not str(term) == 'nan']
gene_types = [str(term) for term in list(df["gene_biotype"]) if not str(term) == 'nan']
ID_to_name = {ID:gene_names[i] for i,ID in enumerate(gene_IDs)}
ID_to_locus = {ID:locus_tags[i] for i,ID in enumerate(gene_IDs)}
ID_types = {ID:gene_types[i] for i,ID in enumerate(gene_IDs)}
locus_to_name = {}
locus_to_types = {}

gene_types = sorted(list(set(gene_types)))
print(gene_types)

completed_quants = []
count_table = {}
genes = []
gene_type = {}
gene_type_genes = {}
total_counts = {}
total_genes = {}
missing_files = []
for term in pool_type:
	print(pool_type[term], term)
	exists = os.path.isfile(main_folder + 'featureCounts/' + term + '.count')
	if exists:
		completed_quants.append(pool_type[term] + ' ' + term)
		in_file = open(main_folder + 'featureCounts/' + term + '.count', 'r')
		lines = in_file.readlines()
		in_file.close()
		
		count_table[term] = {}
		gene_type[term] = {typ:0 for typ in gene_types}
		gene_type_genes[term] = {typ:0 for typ in gene_types}
		total_counts[term] = 0
		total_genes[term] = 0
		for line in lines:
			if not line[0] in ['#', 'G']:
				line_split = line.split('\t')
				#print(line_split)
				ID = line_split[0]
				locus = ID_to_locus[ID]
				if locus[-2] == '_':
					#print(locus)
					locus = locus[:-2]
					#print(locus)
					#print()
				count = int(line_split[6])
				
					
				#if locus in locus_to_name and not locus_to_name[locus] == ID_to_name[ID]:
				#	print('name', locus_to_name[locus], ID_to_name[ID])
				#	input(' .. ')
				
				
				#if 'asso_' in locus and ID_types[ID] in ['three_prime_UTR', 'five_prime_UTR']:
				#	locus = locus.replace('asso_', '')
				
				#if locus in locus_to_types:
				#	if locus_to_types[locus] in ['three_prime_UTR', 'five_prime_UTR'] and not ID_types[ID] in ['three_prime_UTR', 'five_prime_UTR']:
				#		locus_to_types[locus] = ID_types[ID]
				#		
				#else:		
				#	locus_to_types[locus] = ID_types[ID]
				#
				#if locus_to_types[locus] in ['three_prime_UTR', 'five_prime_UTR']:
				#	print(locus, locus_to_types[locus])
					
				
				#locus_to_name[locus] = ID_to_name[ID]
				#if locus in locus_to_types and not locus_to_types[locus] == ID_types[ID]:
				#	print('type', locus_to_types[locus], ID_types[ID])
				#	#input(' .. ')
				#print(locus, ID_types[ID], count)
				if not ID in count_table[term]:
					count_table[term][ID] = count
				else:
					#input(' .. ')
					count_table[term][ID] += count
				if not ID in genes:
					genes.append(ID)
				gene_type[term][ID_types[ID]] += count
				total_counts[term] += count
				gene_type_genes[term][ID_types[ID]] += int(count > 0)
				total_genes[term] += int(count > 0)
				
				
	else:
		missing_files.append(pool_type[term] + ' ' + term)



if True:
	
	sorted_gene_type_genes = np.transpose(sorted([[gene_type_genes[term]['CDS'], term] + [gene_type_genes[term][typ] for typ in gene_types if not typ == 'CDS'] + [total_genes[term]] for term in gene_type])[::-1])
	sorted_samples = sorted_gene_type_genes[1]
	
	out_table = {'Sample id':list(sorted_samples)}
	out_table['Sample type'] = [pool_type[term] for term in sorted_samples]
	out_table['CDS'] = list(sorted_gene_type_genes[0].astype(float))
	
	i = 2
	for typ in gene_types:
		if not typ == 'CDS':
			out_table[typ] = list(sorted_gene_type_genes[i].astype(float))
			i += 1
	
	out_table['total'] = list(sorted_gene_type_genes[i].astype(float))
	
	df = pd.DataFrame(out_table)
	df.to_excel("Biotype_gene_table.xlsx")  
	
	sorted_gene_type_counts = np.transpose(sorted([[gene_type[term]['CDS'], term] + [gene_type[term][typ] for typ in gene_types if not typ == 'CDS'] + [total_counts[term]] for term in gene_type])[::-1])
	sorted_samples = sorted_gene_type_counts[1]
	
	out_table = {'Sample id':list(sorted_samples)}
	out_table['Sample type'] = [pool_type[term] for term in sorted_samples]
	out_table['CDS'] = list(sorted_gene_type_counts[0].astype(int))
	
	i = 2
	for typ in gene_types:
		if not typ == 'CDS':
			out_table[typ] = list(sorted_gene_type_counts[i].astype(int))
			i += 1
	
	out_table['total'] = list(sorted_gene_type_counts[i].astype(int))
	
	df = pd.DataFrame(out_table)
	df.to_excel("Biotype_count_table.xlsx")
	
	
	colors = ['black', 'red', 'white', 'orange', 'darkblue', 'pink', 'darkgrey', 'darkred', 'purple', 'yellow', 'cyan', 'green']
	xticks = []
	plt.figure(figsize=(18, 4))
	x_val = 0
	#print(len(gene_types), gene_types)
	for i,term in enumerate(list(gene_type)):
		#print(term, gene_type[term])
		if not total_counts[term] == 0:
			last_val = 0
			handles = []
			
			for j,typ in enumerate(gene_types):
				#print(j, typ)
				handle = plt.bar(x_val, gene_type[term][typ]/total_counts[term], width = 1, bottom=last_val, color = colors[j], label = typ)
				last_val += gene_type[term][typ]/total_counts[term]
				handles.append(handle)
			xticks.append(pool_type[term] + ' ' + term)
			x_val += 1

	#print(handles)
	plt.legend(handles=handles, bbox_to_anchor=(1.01, 1), loc='upper left', fontsize = 10)
	plt.yticks(fontsize = 10)
	plt.ylabel('Biotype count fraction', fontsize = 10)
	plt.xticks([i for i in range(len(xticks))], [term for term in xticks], fontsize = 8)
	plt.xlim([-.5, len(xticks)-.5])
	plt.setp(plt.gca().xaxis.get_majorticklabels(), rotation=60, ha="right", rotation_mode="anchor")
	plt.tight_layout()
	plt.savefig('All_samples_biotypes_counts.png')
	plt.clf()


#colors = ['black', 'red', 'white', 'orange', 'darkblue', 'pink', 'darkgrey', 'darkred', 'purple', 'yellow', 'cyan', 'grey']
xticks = []
plt.figure(figsize=(18, 4))
x_val = 0
for i,term in enumerate(list(gene_type)):
	#print(term, gene_type[term])
	if not total_counts[term] == 0:
		last_val = 0
		handles = []
		
		for j,typ in enumerate(gene_types):
			handle = plt.bar(x_val, gene_type[term][typ]/total_counts[term], width = 1, bottom=last_val, color = colors[j], label = typ)
			last_val += gene_type[term][typ]/total_counts[term]
			handles.append(handle)
		xticks.append(pool_type[term] + ' ' + term)
		x_val += 1

print(handles)
plt.legend(handles=handles, bbox_to_anchor=(1.01, 1), loc='upper left', fontsize = 10)
plt.yticks(fontsize = 10)
plt.ylabel('Biotype count fraction', fontsize = 10)
plt.xticks([i for i in range(len(xticks))], [term for term in xticks], fontsize = 8)
plt.xlim([-.5, len(xticks)-.5])
plt.setp(plt.gca().xaxis.get_majorticklabels(), rotation=60, ha="right", rotation_mode="anchor")
plt.tight_layout()
plt.savefig('All_samples_biotypes_counts.png')
plt.clf()

#colors = ['black', 'red', 'white', 'orange', 'darkblue', 'pink', 'darkgrey', 'darkred', 'purple', 'yellow', 'cyan', 'grey']
xticks = []
plt.figure(figsize=(18, 4))
x_val = 0
print()
for i,term in enumerate(list(gene_type_genes)):
	if not total_genes[term] == 0:
		last_val = 0
		handles = []
		for j,typ in enumerate(gene_types):
			print(j, typ)
			handle = plt.bar(x_val, gene_type_genes[term][typ]/total_genes[term], width = 1, bottom=last_val, color = colors[j], label = typ)
			last_val += gene_type_genes[term][typ]/total_genes[term]
			handles.append(handle)
		xticks.append(pool_type[term] + ' ' + term)
		x_val += 1
print(handles)
plt.legend(handles=handles, bbox_to_anchor=(1.01, 1), loc='upper left', fontsize = 10)
plt.yticks(fontsize = 10)
plt.ylabel('Biotype gene fraction', fontsize = 10)
plt.xticks([i for i in range(len(xticks))], [term for term in xticks], fontsize = 8)
plt.xlim([-.5, len(xticks)-.5])
plt.setp(plt.gca().xaxis.get_majorticklabels(), rotation=60, ha="right", rotation_mode="anchor")
plt.tight_layout()
plt.savefig('All_samples_biotypes_genes.png')
plt.clf()

#print(gene_type_genes['vibrio56'], total_genes['vibrio56'])
print(missing_files)	
PC_mat = []
PC_labels = []
full_mat = []
full_labels = []
PC_ave = np.zeros(len(genes))
for term in count_table:
	vec = []
	for gene in genes:
		vec.append(np.log2(count_table[term][gene]+1))
	if pool_type[term] == 'PC':
		PC_labels.append(term)
		PC_mat.append(vec)
		PC_ave = PC_ave + np.array(2**np.array(vec)-1)/5
	full_mat.append(vec)
	full_labels.append(pool_type[term] + ' ' + term)
print(np.shape(full_mat))
if False:
	PC_ave = np.log2(PC_ave + 1)
	for i,term in enumerate(full_labels):
		if not 'PC' in term:
			scores = pearsonr(PC_ave, full_mat[i])
			plt.figure(figsize=(5, 5))
			plt.plot(2**np.array(PC_ave)-1, 2**np.array(full_mat[i])-1, marker = 'o', ms = 3, lw = 0)
			plt.loglog()
			plt.xlabel('PC ave.', fontsize = 14)
			plt.ylabel(term.split(' ')[0] + ' ('+term.split(' ')[1]+')', fontsize = 14)
			plt.title('Pearson r = '+str(round(scores.statistic,3)), fontsize = 14)
			plt.tight_layout()
			#.svg', dpi=350
			plt.savefig('PCave' + '_' + term.split(' ')[0] + '_' + term.split(' ')[1] + '_scatter.png')
			plt.clf()


gene_mat = []

for vec in full_mat:
	gene_mat.append(list((np.array(vec)>0).astype(int)))

Z = hierarchy.linkage(gene_mat, 'single')

plt.figure(figsize=(18, 6))

dn = hierarchy.dendrogram(Z, labels=full_labels, distance_sort='descending', leaf_rotation=90., leaf_font_size=8.)
hierarchy.set_link_color_palette(['m', 'c', 'y', 'k'])
#fig, axes = plt.subplots(1, 2, figsize=(8, 3))
hierarchy.set_link_color_palette(None)  # reset to default after use
plt.tight_layout()
plt.savefig('All_samples_dendrogram_genes.png')
plt.clf()
#exit()

Z = hierarchy.linkage(full_mat, 'single')

plt.figure(figsize=(18, 6))

dn = hierarchy.dendrogram(Z, labels=full_labels, distance_sort='descending', leaf_rotation=90., leaf_font_size=8.)
hierarchy.set_link_color_palette(['m', 'c', 'y', 'k'])
#fig, axes = plt.subplots(1, 2, figsize=(8, 3))
hierarchy.set_link_color_palette(None)  # reset to default after use
plt.tight_layout()
plt.savefig('All_samples_dendrogram_counts.png')
plt.clf()
if False:
	plt.figure(figsize=(5, 5))
	for i,term in enumerate(PC_labels[:-1]):
		for j in range(i+1, len(PC_labels)):
			scores = pearsonr(PC_mat[i], PC_mat[j])
			plt.plot(2**np.array(PC_mat[i])-1, 2**np.array(PC_mat[j])-1, marker = 'o', ms = 3, lw = 0)
			plt.loglog()
			plt.xlabel('PC ' + str(i+1) + ' ('+term+')', fontsize = 14)
			plt.ylabel('PC ' + str(j+1) + ' ('+PC_labels[j]+')', fontsize = 14)
			plt.title('Pearson r = '+str(round(scores.statistic,3)), fontsize = 14)
			plt.tight_layout()
			#.svg', dpi=350
			plt.savefig('PC' + str(i+1) + '_PC' + str(j+1) + '_scatter.png')
			plt.clf()

sample_ids = sorted([[sum([count_table[term][gene] for gene in genes]), term] for term in count_table])[::-1]
sorted_genes = sorted([[sum([count_table[term][gene] for term in count_table]), gene] for gene in genes])[::-1]


sample_ids = np.transpose(sample_ids)[1]
sorted_genes = np.transpose(sorted_genes)[1]

sample_types = [pool_type[term] for term in sample_ids]	
out_table = {' ':['', '', 'Sample id'] + list(sample_ids), '':['', '', 'Sample type'] + list(sample_types)}

sorted_CDS = [gene for gene in sorted_genes if ID_types[gene] in ['CDS']]
single_cell_ids = [term for term in sample_ids if 'SC' in pool_type[term]]
CDS_per_cell = {gene:[count_table[term][gene] for term in single_cell_ids] for gene in sorted_CDS}

xticks = []
logcount_mat = []
logcount_means = []
for i,gene in enumerate(sorted_CDS[:100]):
	print(gene, ID_to_name[gene], CDS_per_cell[gene])
	logcount_mat.append(np.log10(np.array(CDS_per_cell[gene]) + 1))
	logcount_means.append(np.mean(np.log10(np.array(CDS_per_cell[gene]) + 1)))
	xticks.append(gene + ' ' + ID_to_name[gene])

plt.violinplot(logcount_mat,
                  showmeans=False,
                  showmedians=True)

plt.plot(logcount_means, lw = 1.2)
plt.yticks(fontsize = 10)
plt.ylabel('log10(count + 1)', fontsize = 10)
plt.xticks([i for i in range(len(xticks))], [term for term in xticks], fontsize = 8)
plt.xlim([-.5, len(xticks)-.5])
plt.setp(plt.gca().xaxis.get_majorticklabels(), rotation=60, ha="right", rotation_mode="anchor")
plt.tight_layout()
plt.savefig('Single_cell_logcount_dists.png')
plt.clf()


for gene in sorted_genes:
	out_table[gene] = [ID_to_name[gene], ID_to_locus[gene], ID_types[gene]] + [count_table[term][gene] for term in sample_ids]

print(len(completed_quants))
df = pd.DataFrame(out_table)
print(np.shape(df))
df.to_excel("Count_table.xlsx")  
#completed_quants = sorted(completed_quants)
#for term in completed_quants:
#	print(term)
#/media/willow/My Passport1/Suyen/2024-12-18_Suyen_Espinoza_24207PR_raw_FASTQ/featureCounts




