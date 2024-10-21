**Steps to test the MGSE3.py script:**

git clone https://github.com/bpucker/MGSE

cd Test_files

There are three possible test cases for the script.

```
Case 1) Testing with FASTA and FASTQ inputs:

python3 MGSE3.py [--fasta <FASTA> --fastq <FASTQ> --seqtech ont]
                 [--busco <BUSCO_TSV> | --gff <GFF3_FILE> | --all <ALL_POS>]
                 --out <DIR>

Case 2) Testing with bam file as input:

python3 MGSE3.py [--bam <BAM> --bam_is_sorted]
                 [--busco <BUSCO_TSV> | --gff <GFF3_FILE> | --all <ALL_POS>]
                 --out <DIR>

Case 3) Testing with coverage file as input:

python3 MGSE3.py [--cov <COV>]
                 [--busco <BUSCO_TSV> | --gff <GFF3_FILE> | --all <ALL_POS>]
                 --out <DIR>
```

**Test files details:**

FASTA     ->  Arabidopsis_thaliana.fa

FASTQ     ->  Arabidopsis_thaliana.fq

BAM       ->  Arabidopsis_thaliana.bam

COV       ->  Arabidopsis_thaliana.cov

BUSCO_TSV ->  Arabidopsis_thaliana_busco.tsv

GFF3_FILE ->  Arabidopsis_thaliana.gff


**More Information:** 

i) The fastq test file has long reads (Oxford Nanopore Technology sequencing reads)

ii) The .bam test file is already sorted.

iii) The test files have information corresponding to the first 10^5 base pairs on chromosome 1 of
    _Arabidopsis thaliana_ Col-0 accession. Hence you must get a genome size estimate of ~0.1 Mbp 
    when testing MGSE3.1.py with these test files.

**Note:**

For more details on the different parameters, please refer to the script README.


