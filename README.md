[![DOI](https://zenodo.org/badge/169685011.svg)](https://doi.org/10.5281/zenodo.2636732)

# Mapping-based Genome Size Estimation (MGSE)

MGSE can harness the power of files generated in genome sequencing projects to predict the genome size. Required are the FASTA file containing a high continuity assembly and a BAM file with all available reads mapped to this assembly. The script `construct_cov_file.py` (https://doi.org/10.1186/s12864-018-5360-z) allows the generation of a COV file based on the (sorted) BAM file (also possible via MGSE directly). Next, this COV file can be used by MGSE to calculate the coverage in provided reference regions and to calculate the total number of mapped bases. Both values are subjected to the genome size estimation. Providing accurate reference regions is crucial for this genome size estimation. Different alternatives were evaluated and actual single copy BUSCOs (https://busco.ezlab.org/) appear to be the best choice. Running BUSCO prior to MGSE will generate all necessary files.

<a href="https://doi.org/10.1101/607390 ">
<img alt="MGSE workflow (Pucker, 2021; doi:10.1101/607390)" src="https://raw.githubusercontent.com/bpucker/MGSE/df459458173ec0bcffba7ce61484e5765118c817/MGSE_concept.png" width="75%" height="75%">
</a>


```
Usage:
  python3 MGSE3.py [--cov <COV_FILE_OR_DIR> | --bam <BAM_FILE_OR_DIR> | --fasta <FASTA_FILE>
		   --fastq <FASTQ_FILE_OR_DIR>]
                   --out <DIR>
                   [--ref <TSV> | --gff <GFF> | --busco <FULL_TABLE.TSV> | --all]

Mandatory:
  Coverage data (choose one)
  [--cov               STR       Coverage file (COV) created by construct_cov_file.py or directory containing
                                 multiple coverage files]
  [--bam               STR       BAM file to automatically create the coverage file]
  [--fasta             STR       FASTA assembly file for read mapping
   --fastq             STR       Single FASTQ file or a directory containing FASTQ files]
  
  Output directory
  --out                STR       Output directory

  Reference regions to calculate average coverage (choose one)
  --ref                STR       TAB-separated file with sequence ID, start position, and end position
  --gff                STR       GFF3 file 
  --busco              STR       BUSCO annotation file (full_table.tsv)
  --all                          Use all positions of the assembly
		
Optional:
  --black              STR       Sequence ID list for exclusion
  --samtools           STR       Full path to samtools (if not in your $PATH)
  --bedtools           STR       Full path to bedtools (if not in your $PATH)
  --seqtech            STR       Sequencing technology of FASTQ (illumina|ont|pacbio)[illumina]
  --short_read_aligner STR       Full path to BWA MEM (if not in your $PATH)
  --long_read_aligner  STR       Full path to minimap2 (if not in your $PATH) 
  --name               STR       Prefix for output files []
  --feature            STR       Specific feature for analysis from GFF file (if other than 'gene')
  --m                  INT       Samtools sort memory [5000000000]
  --threads            INT       Samtools sort threads [4]
  --plot               BOOLEAN   Activate or deactivate generation of figures via matplotlib[FALSE]
  --black_list_factor  FLOAT     Black list factor for blacklisting of contigs with high coverage values [1.5]
  --ignore             BOOLEAN   Deactivate the black listing of contigs with high coverage values [FALSE]
  --gzip                         Search for files "*cov.gz" in --cov if this is a directory
  --gzip_fq                      Specify this flag if the FASTQ file(s) are compressed
  --bam_is_sorted                Do not sort BAM file
```
				
__WARNING:__
- MGSE requires absolute paths (at least use of absolute paths is recommended)
- Per default contigs with very high coverage values are put on a black list to prevent inflation of the genome size prediction    by plastome contigs (in plants). However, this function can be disabled via --ignore to estimate genome sizes with more        fragmented assemblies.


__Possible reference regions:__

1) `--ref` A very simple TAB-separated text file with information about chromosome, start, and end of regions which should be        used as a reference set for the coverage calculation.

2) `--gff` A GFF3 file with genes which should serve as reference regions.

3) `--busco` This will extract the single copy BUSCOs from the provided TSV file.

4) `--all` All positions of the assembly will be included in the average coverage calculation.

__Important pointers:__

1) '--fastq' If you have multiple FASTQ files,
    It is mandatory to put them in a folder and specify the folder path. While specifying the folder path ensure that it ends        with a backslash '/' since it is mandatory for the script to work correctly.

2) '--cov' Extension of uncompressed cov file should be .cov or .txt;
   Compressed cov file extension should be .cov.gz or.txt.gz


```
Usage:
  python2 construct_cov_file.py

  --in <BAM_FILE> --out <OUTPUT_FILE>

Mandatory:
  --in STR          Bam file
  --out STR         Output file

Optional:
  --bam_is_sorted   Do not sort BAM file
  --m INT           Samtools sort memory [5000000000]
  --threads INT     Samtools sort threads [4]
  --samtools STR    Full path to samtools (if not in your $PATH)
  --bedtools STR    Full path to bedtools (if not in your $PATH)
```


```
Usage:
  python3 busco2ploidy.py

  --in <BUSCO_TSV_FILE> --out <OUTPUT_DIR>

Mandatory:
--in  STR	BUSCO full_table.tsv
--out STR	Output directory




## Reference:
Natarajan S., Gehrke J., Pucker B. (2024). Mapping-based genome size estimation. bioRxiv 607390; doi: [10.1101/607390](https://doi.org/10.1101/607390).

