[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2636733.svg)](https://doi.org/10.5281/zenodo.2636733)

# Mapping-based Genome Size Estimation (MGSE)

MGSE can harness the power of files generated in genome sequencing projects to predict the genome size. Required are the FASTA file containing a high continuity assembly and a BAM file with all available reads mapped to this assembly. The script construct_cov_file.py (https://doi.org/10.1186/s12864-018-5360-z) allows the generation of a COV file based on the (sorted) BAM file (also possible via MGSE directly). Next, this COV file can be used by MGSE to calculate the coverage in provided reference regions and to calculate the total number of mapped bases. Both values are subjected to the genome size estimation. Providing accurate reference regions is crucial for this genome size estimation. Different alternatives were evaluated and actual single copy BUSCOs (https://busco.ezlab.org/) appear to be the best choice. Running BUSCO prior to MGSE will generate all necessary files.


python MGSE.py \
--cov <FULL_PATH_TO_COVERAGE_FILE_OR_FOLDER> | --bam <FULL_PATH_TO_BAM_FILE> \
--out <FULL_PATH_TO_OUTPUT_DIRECTORY> \
--ref | --gff <FULL_PATH_TO_REF_GENE_FILE_OR_GFF3_FILE> | --busco <FULL_PATH_TO 'full_table_busco_run.tsv'> | --all <ALL_POS_USED_FOR_CALCULATION>
		
optional: \
--black <FULL_PATh_TO_FILE_WITH_CONTIG_NAMES_FOR_EXCLUSION> \
--gzip <ACITVATES_SEARCH_FOR_COMPRESSED_FILES> \
--bam_is_sorted PREVENTS_SORTING_OF_BAM_FILE \
--samtools <FULL_PATH_TO_SAMTOOLS> \
--bedtools <FULL_PATH_TO_BEDTOOLS> \
--name <NAME_OF_CURRENCT_ANALYSIS> \
--m <SAMTOOLS_MEMORY>[5000000000] \
--threads <SAMTOOLS_THREADS>[4]

				
WARNING: if --busco is used, the BUSCO GFF3 files need to be in the default folder relative to the provided TSV file

WARNING: MGSE requires absolute paths


Possible reference regions:

1) --ref = A very simple TAB-separated text file with information about chromosome, start, and end of regions which should be used as a reference set for the coverage calculation.

2) --gff = A GFF3 file with genes which should serve as reference regions.

3) --busco = This will extract the single copy BUSCOs from the provided TSV file and retrieves the corresponding annotation from GFF3 files generated while running BUSCO.

4) --all = All positions of the assembly will be included in the average coverage calculation.



python construct_cov_file.py \
--in <BAM_FILE> \
--out <OUTPUT_FILE> \
--bam_is_sorted PREVENTS_EXTRA_SORTING_OF_BAM_FILE \
--m <SAMTOOLS_MEMORY>[5000000000] \
--threads <SAMTOOLS_THREADS>[4]






## Reference:

Pucker B. Mapping-based genome size estimation. bioRxiv 607390; doi: https://doi.org/10.1101/607390

