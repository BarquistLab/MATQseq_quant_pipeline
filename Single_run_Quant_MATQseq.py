import os
import sys
import pandas as pd
import zipfile
import numpy as np

#file_path_1 = '/home/willow/'
file_path_1 = '/home/wkioncro/'
#file_path_2 = '/home/willow/'
#file_path_2 = '/media/willow/My\ Passport1/'
file_path_2 = '/vol/projects/wkioncro/'

# Retrieve run id

index = int(sys.argv[1])

lib_la = 'SINGLE'
#main_folder = '2024-08-02_Suyen_Espinoza_24121PR_raw_FASTQ'
#main_folder = '2024-12-18_Suyen_Espinoza_24207PR_raw_FASTQ'
#main_folder = '2024-07-22_Suyen_Espinoza_24121PR_raw_FASTQ'
#main_folder = '2025-04-10_Suyen_Espinoza_25029PR'
main_folder = '2025-07-28_Christopher_Birk_25062PR_raw_FASTQ/Pool2'
folder_path = file_path_2 + 'Single_cell_bacteria/CDiff/'+main_folder+'/FASTQ_by_lib/'
run_ids = [name.split('.')[0] for name in os.listdir(folder_path + 'seq_files') if not '_fastqc' in name]
run_id = run_ids[index]
print(run_ids)
print(run_id)
if True:
	# This code is for generating the initial read file fastqc
	# reports, and the is usually not necessary.
	fastq_files = [folder_path + "seq_files/" + run_id + ".fq.gz"]
	for i,fq_file in enumerate(fastq_files):
		result_file = fq_file.split('.')[0] + '_fastqc.zip'
		print(result_file)
		exists = os.path.isfile(result_file)
		if not exists:
			os.system('fastqc -f fastq -t 1 '+fq_file)
	#exit()

exists = os.path.isfile(file_path_2 +'Single_cell_bacteria/CDiff/'+main_folder+'/FASTQ_by_lib/featureCounts/' + run_id + '.count')
print(file_path_2 +'Single_cell_bacteria/CDiff/'+main_folder+'/FASTQ_by_lib/featureCounts/' + run_id + '.count', ' File exists: ', exists)
if exists:
	print(run_id + ' already quantified. Exiting...')
	exit()

# Trim left FastQ files
exists = os.path.isfile(folder_path + 'BBDuk_L/' + run_id + "_1"*int(lib_la == 'PAIRED') + ".fq.gz")
if not exists:
	os.system("mkdir " + folder_path + "BBDuk_L")
	if lib_la == 'PAIRED':
		os.system(file_path_1 + "bbmap/bbduk.sh -Xmx10g t=20 in1=" + folder_path + "seq_files/" + run_id + "_1.fq.gz"+\
		 			         		" in2=" + folder_path + "seq_files/" + run_id + "_2.fq.gz"+\
		 			         		" out1=" + folder_path + "BBDuk_L/" + run_id + "_1.fq.gz"+\
		 			         		" out2=" + folder_path + "BBDuk_L/" + run_id + "_2.fq.gz ref=" + folder_path + "matqseq_primers.fa minlen=18 qtrim=rl trimq=20 ktrim=l k=17 mink=11 hdist=1 tpe tbo trimpolya=30")
		 
		
	elif lib_la == 'SINGLE':
		#$BBDUK -Xmx10g t=20 in1=$dir/squid_seq_files/$r1 in2=$dir/squid_seq_files/$r2 out1=$dir/BBDuk_L/$r1 out2=$dir/BBDuk_L/$r2 \
		#ref=$dir/matqseq_primers.fa minlen=18 qtrim=rl trimq=20 ktrim=l k=17 mink=11 hdist=1 tpe tbo trimpolya=30
		os.system(file_path_1 + "bbmap/bbduk.sh -Xmx10g t=20 in=" + folder_path + "seq_files/" + run_id + ".fq.gz "+\
			"out=" + folder_path + "BBDuk_L/" + run_id + ".fq.gz ref=" + folder_path + "matqseq_primers.fa minlen=18 qtrim=rl trimq=20 ktrim=l k=17 mink=11 hdist=1 trimpolya=30")
else:
	print("Already left-trimmed FastQ files")

# Trim right FastQ files
exists = os.path.isfile(folder_path + 'BBDuk_L_R/' + run_id + "_1"*int(lib_la == 'PAIRED') + ".fq.gz")
if not exists:
	os.system("mkdir " + folder_path + "BBDuk_L_R")
	if lib_la == 'PAIRED':
		os.system(file_path_1 + "bbmap/bbduk.sh -Xmx10g t=20 in1=" + folder_path + "BBDuk_L/" + run_id + "_1.fq.gz"+\
		 			         		" in2=" + folder_path + "BBDuk_L/" + run_id + "_2.fq.gz"+\
		 			         		" out1=" + folder_path + "BBDuk_L_R/" + run_id + "_1.fq.gz"+\
		 			         		" out2=" + folder_path + "BBDuk_L_R/" + run_id + "_2.fq.gz ref=" + folder_path + "nextera_and_primers.fa minlen=18 qtrim=rl trimq=20 ktrim=r k=17 mink=11 hdist=1 tpe tbo")
		 
	elif lib_la == 'SINGLE':
		os.system(file_path_1 + "bbmap/bbduk.sh -Xmx10g t=20 in=" + folder_path + "BBDuk_L/" + run_id + ".fq.gz out=" + folder_path + "BBDuk_L_R/" + run_id + ".fq.gz ref=" + folder_path + "nextera_and_primers.fa minlen=18 qtrim=rl trimq=20 ktrim=r k=17 mink=11 hdist=1")
else:
	print("Already right-trimmed FastQ files")

folder_path = file_path_2 + 'Single_cell_bacteria/CDiff/'+main_folder
exists = os.path.isfile(folder_path + '/quality_stats_by_lib/' + run_id + '_QC_results.txt')
if not exists:
	# First QC step: FastQC
	os.system("mkdir " + folder_path + '/quality_stats_by_lib')
	fastqc_fail_cols = ['Per base sequence quality', 'Per sequence quality scores', 'Per base N content']
	if lib_la == 'PAIRED':
		sample_pass = [[True], [True]]
		fastq_files = [folder_path + "/FASTQ_by_lib/BBDuk_L_R/" + run_id + "_1.fq.gz", folder_path + "/FASTQ_by_lib/BBDuk_L_R/" + run_id + "_2.fq.gz"]
		
	elif lib_la == 'SINGLE':
		sample_pass = [[True]]
		fastq_files = [folder_path + "/FASTQ_by_lib/BBDuk_L_R/" + run_id + ".fq.gz"]
		
	for i,fq_file in enumerate(fastq_files):
		result_file = fq_file.split('.')[0] + '_fastqc.zip'
		print(result_file)
		exists = os.path.isfile(result_file)
		if not exists:
			os.system('fastqc -f fastq -t 1 '+fq_file)
		
		exists = os.path.isfile(result_file)
		if not exists:
			print("FastQC failed. Exiting...")
			exit()
		
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

	pass_1st_QC = np.prod([int(term[0]) for term in sample_pass]) == 1
	if pass_1st_QC:
		print(run_id + ' passed QC check')
	else:
		print(run_id + ' did NOT pass QC check')
	
	out_file = open(folder_path + '/quality_stats_by_lib/' + run_id + '_QC_results.txt', 'w')
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
	print(run_id + ' already checked for QC. ')
	in_file = open(folder_path + '/quality_stats_by_lib/' + run_id + '_QC_results.txt')
	lines = in_file.readlines()
	in_file.close()
	sample_pass = []
	#percent_unmap = []
	for line in lines:
		line_split = line.split()
		if len(line_split) == 2:
			
			if line_split[0].split('_')[0] == run_id:
				sample_pass.append(bool(line_split[1] == 'True'))
			#elif '%Hit_no_genomes:' in line:
			#	percent_unmap.append(float(line_split[1]))
				
	pass_1st_QC = np.prod([int(term) for term in sample_pass]) == 1
	#pass_2nd_QC = np.prod([int(float(term) <= 100) for term in percent_unmap]) == 1

# Run quantification
if pass_1st_QC:
	# and pass_2nd_QC
	print('QC check passed. Quantification will run.')
	#os.system(file_path_2 + 'salmon2 index -t ' + file_path_2 + 'Genomes/Ecoli/e_coli_k12_mg1655_U00096_mod.fa -i ' + file_path_2 + 'e_coli_k12_mg1655_index -k 29')
	#if lib_la == 'PAIRED':
	#	os.system(file_path_2 + 'salmon2 quant -i ' + file_path_2 + 'e_coli_k12_mg1655_index -l A --maxSoftclipFraction .8' +\
	#         				' -1 ' + file_path_2 + 'fastq/' + run_id + '/' + run_id + '_1.fastq.gz' +\
	#         				' -2 ' + file_path_2 + 'fastq/' + run_id + '/' + run_id + '_2.fastq.gz' +\
	#         				' -p 1 --validateMappings -o ' + file_path_2 + 'quants/' + run_id + '_quant')
	#elif lib_la == 'SINGLE':
	#	os.system(file_path_2 + 'salmon2 quant -i ' + file_path_2 + 'e_coli_k12_mg1655_index -l A --maxSoftclipFraction .8' +\
	#         				' -r ' + file_path_2 + 'fastq/' + run_id + '/' + run_id + '.fastq.gz' +\
	#         				' -p 1 --validateMappings -o ' + file_path_2 + 'quants/' + run_id + '_quant')
	#         
	os.system('mkdir ' + folder_path + '/FASTQ_by_lib/bowtie2_index')
	os.system('mkdir ' + folder_path + '/FASTQ_by_lib/bowtie2_aligned')
	os.system('mkdir ' + folder_path + '/FASTQ_by_lib/featureCounts')
	
	if not os.path.isfile(folder_path + '/FASTQ_by_lib/bowtie2_index/CP101905.2_index.1.bt2'):
		if index == 0:
			os.system('/vol/biotools/bin/bowtie2-build ' + folder_path + '/FASTQ_by_lib/genome/CP101905.2.fasta ' + folder_path + '/FASTQ_by_lib/bowtie2_index/CP101905.2_index')
		else:
			print('No index file and not first run. Exiting...')
			exit()
		
	if lib_la == 'PAIRED':
		#NOT UPDATED!!
		os.system('bowtie2 -p 20 --local -x ' + folder_path + '/FASTQ_by_lib/bowtie2_index/Vibrio_fischeri_ES114_index -N 1 -1 ' +\
						folder_path + '/FASTQ_by_lib/BBDuk_L_R/' + run_id + '_1.fq.gz -2 ' +\
						folder_path + '/FASTQ_by_lib/BBDuk_L_R/' + run_id + '_2.fq.gz | samtools view -@ 20 -bS -h - > ' +\
						folder_path + '/FASTQ_by_lib/bowtie2_aligned/' + run_id + '.bam')
		
		os.system(file_path_1 + '/anaconda3/envs/featurecounts/bin/featureCounts -T 10 -a ' + folder_path + '/FASTQ_by_lib/genome/Vibrio_fischeri_ES114.gff3 -t "gene" -p --countReadPairs -g "locus_tag" -o ' +\
						folder_path + '/FASTQ_by_lib/featureCounts/' + run_id + '.count ' +\
						folder_path + '/FASTQ_by_lib/bowtie2_aligned/' + run_id + '.bam')
						
	
	elif lib_la == 'SINGLE':
		#	bowtie2 -p 30 --local -x bowtie2_index_hiri/cdiff_index -U $file | samtools view -@ 10 -bS -h - > bowtie2_aligned_hiri/$bam_file			
		os.system('bowtie2 -p 30 --local -x ' + folder_path + '/FASTQ_by_lib/bowtie2_index/CP101905.2_index -U ' +\
						folder_path + '/FASTQ_by_lib/BBDuk_L_R/' + run_id + '.fq.gz | samtools view -@ 10 -bS -h - > ' +\
						folder_path + '/FASTQ_by_lib/bowtie2_aligned/' + run_id + '.bam')
		# featureCounts -T 20 -a genome/CP010905.2_complete_180723_updatedGP.gff3 \
		# -t "CDS,five_prime_UTR,three_prime_UTR,ncRNA,pseudogene,Riboswitch,ribozyme,RNase_P_RNA,rRNA,SRP_RNA,tmRNA,tRNA" -g "ID" -o featureCounts_HIRI/$count_file bowtie2_aligned_hiri/$bam_file
		os.system('featureCounts -T 20 -a ' + folder_path + '/FASTQ_by_lib/genome/CP010905.2_complete_180723_updatedGP.gff3'+\
						' -t "CDS,five_prime_UTR,three_prime_UTR,ncRNA,pseudogene,Riboswitch,ribozyme,RNase_P_RNA,rRNA,SRP_RNA,tmRNA,tRNA" -g "ID" -o ' +\
						folder_path + '/FASTQ_by_lib/featureCounts/' + run_id + '.count ' +\
						folder_path + '/FASTQ_by_lib/bowtie2_aligned/' + run_id + '.bam')
	
	
else:
	print('Quantification will not be run: QC checks did not pass.')




