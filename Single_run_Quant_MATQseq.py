#######################################################
# 
# This script takes a single run, runs adapter trimming,
# FastQC for quality control checks, and then 
# quantifies the reads based on bowtie2 and featureCounts.
# In this code, the following versions were tested:
#
#	BBDuk (BBMap version 39.01)
#	FastQC (v0.11.9)
#	Bowtie2 (2.3.5.1)
#	FeatureCounts (v2.0.8)
#	SAMtools (1.10)
#
#######################################################

import os
import sys
import pandas as pd
import zipfile
import numpy as np

########################################## User parameters #####################################################
lib_la = 'SINGLE'# Library layout
folder_path = 'Example_FASTQ_folder/'# Folder containing folder FASTQ_by_lib which has FastQ files to be quantified
run_ids = ['Cdiff3', 'Cdiff4']# Names of the FastQ, omitting .fq.gz
index_file = 'CP101905.2'# Name of index file
annotation_file = 'CP010905.2_complete_180723_updatedGP.gff3'# Name of annatation file, with type designator
feature_list = ['CDS', 'five_prime_UTR', 'three_prime_UTR', 'ncRNA', 'pseudogene', 'Riboswitch', 'ribozyme', 'RNase_P_RNA', 'rRNA', 'SRP_RNA', 'tmRNA', 'tRNA']# List of features to quantify against
gene_id_attribute = 'ID'# Feature identifier as found in the annotation file
num_threads = 4# Number of threads for featureCounts and bowtie2
################################################################################################################

# Retrieve run id
index = int(sys.argv[1])
run_id = run_ids[index]
print('Running on cell', run_id)

# First check if cell is already quantified. If so, then exit.
exists = os.path.isfile(folder_path +'featureCounts/' + run_id + '.count')
print(folder_path +'featureCounts/' + run_id + '.count', ' File exists: ', exists)
if exists:
	print(run_id + ' already quantified. Exiting...')
	exit()

# Trim left FastQ files
exists = os.path.isfile(folder_path + 'BBDuk_L/' + run_id + "_1"*int(lib_la == 'PAIRED') + ".fq.gz")
if not exists:
	os.system("mkdir " + folder_path + "BBDuk_L")
	if lib_la == 'PAIRED':
		os.system("bbduk.sh t=" + str(num_threads) + " in1=" + folder_path + "FASTQ_by_lib/" + run_id + "_1.fq.gz"+\
		 			         		" in2=" + folder_path + "FASTQ_by_lib/" + run_id + "_2.fq.gz"+\
		 			         		" out1=" + folder_path + "BBDuk_L/" + run_id + "_1.fq.gz"+\
		 			         		" out2=" + folder_path + "BBDuk_L/" + run_id + "_2.fq.gz ref=matqseq_primers.fa minlen=18 qtrim=rl trimq=20 ktrim=l k=17 mink=11 hdist=1 tpe tbo trimpolya=30")
		 
		
	elif lib_la == 'SINGLE':
		os.system("bbduk.sh t=" + str(num_threads) + " in=" + folder_path + "FASTQ_by_lib/" + run_id + ".fq.gz "+\
			"out=" + folder_path + "BBDuk_L/" + run_id + ".fq.gz ref=matqseq_primers.fa minlen=18 qtrim=rl trimq=20 ktrim=l k=17 mink=11 hdist=1 trimpolya=30")
else:
	print("Already left-trimmed FastQ files")

# Trim right FastQ files
exists = os.path.isfile(folder_path + 'BBDuk_L_R/' + run_id + "_1"*int(lib_la == 'PAIRED') + ".fq.gz")
if not exists:
	os.system("mkdir " + folder_path + "BBDuk_L_R")
	if lib_la == 'PAIRED':
		os.system("bbduk.sh t=" + str(num_threads) + " in1=" + folder_path + "BBDuk_L/" + run_id + "_1.fq.gz"+\
		 			         		" in2=" + folder_path + "BBDuk_L/" + run_id + "_2.fq.gz"+\
		 			         		" out1=" + folder_path + "BBDuk_L_R/" + run_id + "_1.fq.gz"+\
		 			         		" out2=" + folder_path + "BBDuk_L_R/" + run_id + "_2.fq.gz ref=nextera_and_primers.fa minlen=18 qtrim=rl trimq=20 ktrim=r k=17 mink=11 hdist=1 tpe tbo")
		 
	elif lib_la == 'SINGLE':
		os.system("bbduk.sh t=" + str(num_threads) + " in=" + folder_path + "BBDuk_L/" + run_id + ".fq.gz out=" + folder_path + "BBDuk_L_R/" + run_id + ".fq.gz ref=nextera_and_primers.fa minlen=18 qtrim=rl trimq=20 ktrim=r k=17 mink=11 hdist=1")
else:
	print("Already right-trimmed FastQ files")

# Check if FastQC summary results exist. If not, run FastQC
exists = os.path.isfile(folder_path + run_id + '_QC_results.txt')
if not exists:
	# First QC step: FastQC
	# These following three QC flags were chosen based on following the 
	# pipeline at https://github.com/SBRG/iModulonMiner/tree/main
	
	fastqc_fail_cols = ['Per base sequence quality', 'Per sequence quality scores', 'Per base N content']
	if lib_la == 'PAIRED':
		sample_pass = [[True], [True]]
		fastq_files = [folder_path + "BBDuk_L_R/" + run_id + "_1.fq.gz", folder_path + "BBDuk_L_R/" + run_id + "_2.fq.gz"]
		
	elif lib_la == 'SINGLE':
		sample_pass = [[True]]
		fastq_files = [folder_path + "BBDuk_L_R/" + run_id + ".fq.gz"]
	
	# Check if fastqc file is already generated.
	# If not, run FastQC (input type fastq, and use 4 threads)
	
	for i,fq_file in enumerate(fastq_files):
		result_file = fq_file.split('.')[0] + '_fastqc.zip'
		
		exists = os.path.isfile(result_file)
		if not exists:
			os.system('fastqc -f fastq -t 1 '+fq_file)
		
		exists = os.path.isfile(result_file)
		if not exists:
			print("FastQC failed. Exiting...")
			exit()
		
		# If FastQC ran successfully, read in the summary file
		# for each of the fastq files.
		
		with zipfile.ZipFile(result_file) as z:
			sum_file = z.open(fq_file.split('.')[0].split('/')[-1] + '_fastqc/summary.txt')
			for line in sum_file.readlines():
				line_mod = str(line).split('\t')[0].split('\\t')
				sample_pass[i].append([line_mod[0][2:], line_mod[1]])
				if line_mod[1] in fastqc_fail_cols:
					
					print(sample_pass[i][-1])
					sample_pass[i][-1].append('*')
					if not line_mod[0][2:] == 'PASS':
						sample_pass[i][0] = False
	
	# Only run quantification of all checks passed.
	
	pass_1st_QC = np.prod([int(term[0]) for term in sample_pass]) == 1
	if pass_1st_QC:
		print(run_id + ' passed QC check')
	else:
		print(run_id + ' did NOT pass QC check')
	
	# Write QC summary file. 
	# This file only includes the checks 
	# that FastQC includes.
	
	out_file = open(folder_path + run_id + '_QC_results.txt', 'w')
	out_file.write('FastQC\n' + run_id+'_1'*int(lib_la == 'PAIRED') + '\t' + str(sample_pass[0][0]) + '\n')
	for term in sample_pass[0][1:]:
		out_file.write('\t'.join(term) + '\n')
	out_file.write('\n')

	if lib_la == 'PAIRED':
		out_file.write(run_id+'_2\t' + str(sample_pass[1][0]) + '\n')
		for term in sample_pass[1][1:]:
			out_file.write('\t'.join(term) + '\n')
			
		out_file.write('\n')
		
	out_file.close()
else:

	# If the QC results file already exists,
	# skip running FastQC and simply check results.
	
	print(run_id + ' already checked for QC. ')
	in_file = open(folder_path + run_id + '_QC_results.txt')
	lines = in_file.readlines()
	in_file.close()
	sample_pass = []
	
	for line in lines:
		line_split = line.split()
		if len(line_split) == 2:
			if line_split[0].split('_')[0] == run_id:
				sample_pass.append(bool(line_split[1] == 'True'))
				
	pass_1st_QC = np.prod([int(term) for term in sample_pass]) == 1

########################################################
# Run quantification if QC checks passed.
########################################################

if pass_1st_QC:
	# and pass_2nd_QC
	print('QC check passed. Quantification will run.')
	
	os.system('mkdir ' + folder_path + 'bowtie2_aligned')
	os.system('mkdir ' + folder_path + 'featureCounts')
	
	#################################### Bowtie2 extra flags ############################################
	#
	# --local tells Bowtie2 to use local alignment rather than end-to-end alignment,
	# 	such that it will soft-clip bases on a read at the ends if those bases don't align well.
	#
	#################################### SAMtools extra flags ##########################################
	#
	# -b = output BAM format (binary version of SAM)
	#
	####################################################################################################
	
	if lib_la == 'PAIRED':
		
		os.system('bowtie2 -p ' + str(num_threads) + ' --local -x ' + 'bowtie2_index/' + index_file + '_index -1 ' +\
						folder_path + 'BBDuk_L_R/' + run_id + '_R1.fq.gz -2 ' +\
						folder_path + 'BBDuk_L_R/' + run_id + '_R2.fq.gz | samtools view -@ ' +\
						str(num_threads) + ' -b -h - > ' +\
						folder_path + 'bowtie2_aligned/' + run_id + '.bam')

		os.system('featureCounts -T ' + str(num_threads) + ' -a genomes/' + annotation_file +\
						' -t "' + ','.join(feature_list) + '" -g "' + gene_id_attribute +\
						'" -p --countReadPairs -o ' +\
						folder_path + 'featureCounts/' + run_id + '.count ' +\
						folder_path + 'bowtie2_aligned/' + run_id + '.bam')		
	
	elif lib_la == 'SINGLE':
		
		os.system('bowtie2 -p ' + str(num_threads) + ' --local -x ' + 'bowtie2_index/' + index_file + '_index -U ' +\
						folder_path + 'BBDuk_L_R/' + run_id + '.fq.gz | samtools view -@ ' +\
						str(num_threads) + ' -b -h - > ' +\
						folder_path + 'bowtie2_aligned/' + run_id + '.bam')
		
		os.system('featureCounts -T ' + str(num_threads) + ' -a genomes/' + annotation_file +\
						' -t "' + ','.join(feature_list) + '" -g "' + gene_id_attribute + '" -o ' +\
						folder_path + 'featureCounts/' + run_id + '.count ' +\
						folder_path + 'bowtie2_aligned/' + run_id + '.bam')
	
else:
	print('Quantification will not be run: QC checks did not pass.')




