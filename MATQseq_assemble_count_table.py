import os
import pandas as pd

########################################## User parameters #####################################################
# Input/output locations relative to the pipeline's main directory.
anno_file = 'genomes/CP010905.2_complete_180723_updatedGP.gff3'
main_folder = 'Example_FASTQ_folder/'
count_folder = main_folder + 'featureCounts'
output_folder = main_folder

# Feature identifier used by featureCounts (-g).
gene_id_attribute = 'ID'

# Feature types to include in the output tables.
biotypes = ['CDS', 'five_prime_UTR', 'three_prime_UTR',
	'ncRNA', 'pseudogene', 'Riboswitch', 'ribozyme',
	'RNase_P_RNA', 'rRNA', 'SRP_RNA', 'tmRNA', 'tRNA']

# Sample types will be specified here from a user-defined .ods table.
# The keys must match the sample names derived from the .count filenames.
sample_types_file = main_folder + 'Cdiff_seq_file_conversion_table.ods'
sample_types_data = pd.read_excel(sample_types_file, engine='odf')
sample_types = dict(zip(sample_types_data['File name'], sample_types_data['File type']))
################################################################################################################


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
					# Use the locus tag when available, otherwise fall back to associated_gene.
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
					}


# Find all featureCounts output files.
count_files = sorted(file_name for file_name in os.listdir(count_folder) if file_name.endswith('.count'))

# Read the featureCounts output for each sample.
sample_counts = {}
missing_annotations = set()

for count_file in count_files:
	sample_name = os.path.splitext(count_file)[0]

	print(count_file)

	count_data = pd.read_csv(os.path.join(count_folder, count_file), sep='\t', comment='#')

	# The final column contains the counts for this sample.
	count_column = count_data.columns[-1]

	counts = dict(zip(count_data['Geneid'], count_data[count_column]))
	sample_counts[sample_name] = counts

	# Record any featureCounts IDs which are not present in the annotation.
	missing_annotations.update(feature_id for feature_id in counts if feature_id not in annotations)

	print('Nonzero count features:', sum(count > 0 for count in counts.values()), '/', len(counts))
	print('Total counts:', sum(counts.values()))
	print()


if missing_annotations:
	print('Warning:', len(missing_annotations), 'featureCounts IDs were not found in the selected annotation file.')

# Sort samples by their total number of assigned counts.
sample_names = sorted(sample_counts, key=lambda sample_name: sum(sample_counts[sample_name].values()), reverse=True)


# Sort features by their total number of assigned counts across all samples.
feature_ids = sorted(annotations, key=lambda feature_id: sum(sample_counts[sample_name].get(feature_id, 0) for sample_name in sample_names), reverse=True)


########################################## Count table #########################################################

# Reproduce the MATQ-seq count-table format:
#   row 1: gene/feature name
#   row 2: locus tag
#   row 3: biotype
#   remaining rows: samples
#
# Each feature is a column, and each sample contains its featureCounts value.

count_table = {' ': ['', '', 'Sample id'] + sample_names,
	      '': ['', '', 'Sample type'] + [sample_types.get(sample_name, '') for sample_name in sample_names]
}

for feature_id in feature_ids:
	count_table[feature_id] = [annotations[feature_id]['Name'], 
			       annotations[feature_id]['Locus tags'], 
			       annotations[feature_id]['Biotype']
			       ] + [sample_counts[sample_name].get(feature_id, 0) for sample_name in sample_names]

df = pd.DataFrame(count_table)
df.to_excel(os.path.join(output_folder, 'Count_table.xlsx'), index=False)


########################################## Biotype tables ######################################################

# Count the number of detected features of each biotype in each sample.
biotype_gene_table = {'Sample id': sample_names,
		'Sample type': [sample_types.get(sample_name, '') for sample_name in sample_names]
}

for biotype in biotypes:
	biotype_gene_table[biotype] = [sum(sample_counts[sample_name].get(feature_id, 0) > 0 for feature_id, annotation in annotations.items() if annotation['Biotype'] == biotype) for sample_name in sample_names]

biotype_gene_table['total'] = [sum(sample_counts[sample_name].get(feature_id, 0) > 0 for feature_id in annotations) for sample_name in sample_names]

df = pd.DataFrame(biotype_gene_table)
df.to_excel(os.path.join(output_folder, 'Biotype_gene_table.xlsx'), index=False)


# Count the number of assigned reads for each biotype in each sample.
biotype_count_table = {'Sample id': sample_names,
		 'Sample type': [sample_types.get(sample_name, '') for sample_name in sample_names]
}

for biotype in biotypes:
	biotype_count_table[biotype] = [sum(sample_counts[sample_name].get(feature_id, 0) for feature_id, annotation in annotations.items() if annotation['Biotype'] == biotype) for sample_name in sample_names]

biotype_count_table['total'] = [sum(sample_counts[sample_name].get(feature_id, 0) for feature_id in annotations) for sample_name in sample_names]

df = pd.DataFrame(biotype_count_table)
df.to_excel(os.path.join(output_folder, 'Biotype_count_table.xlsx'), index=False)
