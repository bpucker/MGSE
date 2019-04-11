# Mapping-based Genome Size Estimation (MGSE)

MGSE can harness the power of files generated in genome sequencing projects to predict the genome size. Required are the FASTA file containing a high continuity assembly and a BAM file with all available reads mapped to this assembly. The script construct_cov_file.py (https://doi.org/10.1186/s12864-018-5360-z) allows the generation of a COV file based on the (sorted) BAM file. Next, this COV file can be used by MGSE to calculate the coverage in provided reference regions and to calculate the total number of mapped bases. Both values are subjected to the genome size estimation. Providing accurate reference regions is crucial for this genome size estimation. Different alternatives were evaluated and actual single copy BUSCOs () appear to be the best choice. Running BUSCO prior to MGSE will generate all necessary files.


python MGSE.py \
--cov <FULL_PATH_TO_COVERAGE_FILE_OR_FOLDER> \
--out <FULL_PATH_TO_OUTPUT_DIRECTORY> \
--ref | --gff <FULL_PATH_TO_REF_GENE_FILE_OR_GFF3_FILE> | --busco <FULL_PATH_TO 'full_table_busco_run.tsv'>
		
optional: \
--black <FULL_PATh_TO_FILE_WITH_CONTIG_NAMES_FOR_EXCLUSION> \
--gzip <ACITVATES_SEARCH_FOR_COMPRESSED_FILES>

				
WARNING: if --busco is used, the BUSCO GFF3 files need to be in the default folder relative to the provided TSV file



python construct_cov_file.py \
--in <BAM_FILE> \
--out <OUTPUT_FILE> \
--bam_is_sorted <PREVENTS_EXTRA_SORTING_OF_BAM_FILE>


