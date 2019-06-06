### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.35 ###

__usage__ = """
					python construct_cov_file.py\n
					--in <BAM_FILE>
					--out <OUTPUT_FILE>
					
					--bam_is_sorted PREVENTS_EXTRA_SORTING_OF_BAM_FILE
					
					optional:
					--m <SAMTOOLS_MEMORY>[5000000000]
					--threads <SAMTOOLS_THREADS>[4]
					--samtools <FULL_PATH_TO_SAMTOOLS>
					--bedtools <FULL_PATH_TO_BEDTOOLS>
					
					feature requests and bug reports: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys

# --- end of imports --- #

def main( arguments ):
	
	bam_file = arguments[ arguments.index( '--in' )+1 ]
	output_file = arguments[ arguments.index( '--out' )+1 ]
	
	if '--samtools' in arguments:
		samtools = arguments[ arguments.index( '--samtools' )+1 ]
	else:
		samtools = "samtools"
	
	if '--bedtools' in arguments:
		bedtools = arguments[ arguments.index( '--bedtools' )+1 ]
	else:
		bedtools = "genomeCoverageBed"
	
	
	if '--bam_is_sorted' in arguments:
		sorted_bam_file = bam_file
	else:
		
		
		if '--m' in arguments:
			m = arguments[ arguments.index( '--m' )+1 ]
		else:
			m = "5000000000"
		
		if '--threads' in arguments:
			t = arguments[ arguments.index( '--threads' )+1 ]
		else:
			t = "4"
				
		print "sorting BAM file ..."
		sorted_bam_file = output_file + "_sorted.bam"
		cmd = samtools + " sort -m " + m + " --threads " + t + " " + bam_file + " > " + sorted_bam_file
		os.popen( cmd )
	
	# --- calculate read coverage depth per position --- #
	print "calculating coverage per position ...."
	cmd = bedtools + " -d -split -ibam " + sorted_bam_file + " > " + output_file
	os.popen( cmd )


if __name__ == '__main__':
	
	if '--in' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
	
	print "all done!"
