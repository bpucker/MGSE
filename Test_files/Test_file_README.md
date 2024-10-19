Steps to test the MGSE3.1.py script with the test file are as follows:

git clone https://github.com/bpucker/MGSE

cd Test_files

You can find different input files and reference region files here. 

Please read the script README for the usage of these different files.

NOTE: 

i) The fastq test file has long reads; If testing with this fastq file input, use minimap2 as the aligner.

ii) The .bam test file is already sorted.

iii) The test files provided have content corresponding to the first 10^5 base pairs on chromosome 1 of
Arabidopsis thaliana Col-0 accession.

iv) You must get a genome size estimate of ~0.1 Mbp when testing MGSE3.1.py with these test files.


