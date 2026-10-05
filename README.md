# MATQ-seq Quantification Pipeline
A pipeline for the quantification of MATQ-seq data. This pipeline follows
an example based on C. difficile. The sequence and annotation files are not
provided in this tutorial for data protection reasons. To run this pipeline,
these files in addition to the FastQ files to be qauntified will have to be 
copied to the folders 'genomes' (for both sequence and annotation files) and 'seq_files', respectively.
This pipeline was tested based on the following software versions,

- Python 3.11.3
- BBDuk (BBMap 39.01)
- FastQC 0.11.9
- Bowtie 2 2.3.5.1
- featureCounts 2.0.8
- SAMtools 1.10

## 1 Read QC, read mapping, and quantification
### 1.1 Build reference genome index

Before mapping reads to features using 'bowtie2', it's first necessary to create an index file 
from the sequence file provided by the user,

```bash
bowtie2-build genomes/CP101905.2.fasta bowtie2_index/CP101905.2_index
```

### 1.2 Main quantification pipeline

The following Python script can now be run for quantification. At the top of the script is a set of 
parameters to specify the sequence, annotation, and FastQ file to be quantified. Additionally,
there is a list of features which the read mapping will be based on before quantifying, which 
the user will first need to fill out.

To decide which features to map reads against, the following 'awk' command can 
be used to list each feature and their frequencies in the user-provided annotation file:

```bash
awk -F '\t' '{print $3}' genomes/CP010905.2_complete_180723_updatedGP.gff3 | sort | uniq -c
```

Any additional software parameter choices which are not defined by the user at the top of the script
are described in the comments of the corresponding section of the code.

To run the pipeline for a single FastQ file, run

```bash
python Single_run_Quant_MATQseq.py 0
```

where the 0 indicates the index of the FastQ file in the list provided by the user.
This code will first run adapter trimming via 'bbduk', followed by 'fastqc' for quality control, then finally 'bowtie2' and 'featureCounts' for quantification.
This code will also check for already-complete steps, and skip ahead. Once complete, the FastQC report can be found next to the trimmed FastQ files in the 'BBDuk_L_R'
folder, and a summary of this report can be found in the 'FASTQ_by_lib' directory. See [Clostridioides difficile quantification and analysis](https://github.com/BarquistLab/C_diff_quantification_and_analysis)
for further details on how to read this report.

### 1.3 Iterative quantification

The use of the index at the end of the Python command can allow the user to then 
easily iterate over cells in the commandline, e.g. for 9 cells as

```bash
for i in {0..9}; do
	python Single_run_Quant_MATQseq.py $i
done
```

or submit the python script as an array job on a compute cluster.

## 2 Count table assembly

