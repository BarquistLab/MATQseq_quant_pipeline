import os
import pandas as pd


# Input/output locations relative to the pipeline's main directory.
anno_file = 'genomes/CP010905.2_complete_180723_updatedGP.gff3'
count_folder = 'Example_FASTQ_folder/FASTQ_by_lib/featureCounts'
output_folder = 'Example_FASTQ_folder'
gene_id_attribute = 'ID'# Feature identifier as found in the annotation file

# These are the feature types which will be searched for in the .count
# files and included in the final tables.
biotypes = ['CDS', 'five_prime_UTR', 'three_prime_UTR',
	'ncRNA', 'pseudogene', 'Riboswitch', 'ribozyme',
	'RNase_P_RNA', 'rRNA', 'SRP_RNA', 'tmRNA', 'tRNA']

# Read the annotation and store information indexed by the feature ID.
annotations = {}

with open(anno_file, 'r') as in_file:
	for line in in_file:
		if not line.startswith('#'):
			line_split = line.rstrip('\n').split('\t')
			if line_split[2] in biotypes:
				attributes = {}

				for item in line_split[8].split(';'):
					if '=' in item:
						key, value = item.split('=', 1)
						attributes[key] = value

				feature_id = attributes.get(gene_id_attribute)

				if feature_id is not None:
					# Collect locus_tag when available, otherwise fall back to associated_gene.
					locus_tag = attributes.get('locus_tag')

					if locus_tag is None:
						associated_gene = attributes.get('associated_gene')
						locus_tag = 'asso_' + associated_gene if associated_gene is not None else 'None'

					# Use the gene name when available, otherwise use the locus tag.
					gene_name = attributes.get('gene', locus_tag)

					annotations[feature_id] = {
						gene_id_attribute: feature_id,
						'Locus tags': locus_tag,
						'Name': gene_name,
						'Biotype': line_split[2],
						'Length': int(line_split[4]) - int(line_split[3]) + 1
					}


# Build the annotation portion of the count table.
df_count_table = {
	gene_id_attribute: [annotation[gene_id_attribute] for annotation in annotations.values()],
	'Locus tags': [annotation['Locus tags'] for annotation in annotations.values()],
	'Name': [annotation['Name'] for annotation in annotations.values()],
	'Biotype': [annotation['Biotype'] for annotation in annotations.values()],
	'Length': [annotation['Length'] for annotation in annotations.values()]
}


# Build the zero-count table.
# The '# features' column gives the total number of annotated features of each type.
df_zeros_table = {
	'Biotypes': biotypes,
	'# features': [
		sum(annotation['Biotype'] == biotype for annotation in annotations.values())
		for biotype in biotypes
	]
}


# Find all featureCounts output files.
count_files = sorted(file_name for file_name in os.listdir(count_folder) if file_name.endswith('.count'))


# Read each sample's featureCounts output and add its counts to the tables.
sample_counts = {}

for count_file in count_files:
	sample_name = os.path.splitext(count_file)[0]

	print(count_file)

	count_data = pd.read_csv(os.path.join(count_folder, count_file), sep='\t', comment='#')

	# The final column contains the counts for this sample.
	count_column = count_data.columns[-1]

	counts = dict(zip(count_data['Geneid'], count_data[count_column]))
	sample_counts[sample_name] = counts

	# Match featureCounts Geneids to the annotation IDs.
	df_count_table[sample_name] = [counts.get(feature_id, 0) for feature_id in annotations]

	# Count the number of zero-count features within each biotype.
	df_zeros_table[sample_name] = [
		sum(
			counts.get(feature_id, 0) == 0
			for feature_id, annotation in annotations.items()
			if annotation['Biotype'] == biotype
		)
		for biotype in biotypes
	]

	total_counts = sum(counts.get(feature_id, 0) for feature_id in annotations)

	print('Nonzero count features:', sum(count > 0 for count in df_count_table[sample_name]), '/', len(annotations))
	print('Total counts:', total_counts)
	print()


# Write the feature count table as a tab-separated file.
df = pd.DataFrame(df_count_table)
df.to_csv(os.path.join(output_folder, 'count_table.tsv'), sep='\t', index=False)


# Write the zero-count table as a tab-separated file.
df = pd.DataFrame(df_zeros_table)
df.to_csv(os.path.join(output_folder, 'zeros_table.tsv'), sep='\t', index=False)


# Calculate the percentage of total counts assigned to each biotype.
percentage_table = {'': ['% ' + biotype for biotype in biotypes]}

for sample_name, counts in sample_counts.items():
	total_counts = sum(counts.get(feature_id, 0) for feature_id in annotations)
	percentage_table[sample_name] = [100 * sum(counts.get(feature_id, 0) for feature_id, annotation in annotations.items() if annotation['Biotype'] == biotype) / total_counts for biotype in biotypes]

# Write the percentage table as a tab-separated file.
df = pd.DataFrame(percentage_table)
df.to_csv(os.path.join(output_folder, 'percentage_counts.tsv'), sep='\t', index=False)
