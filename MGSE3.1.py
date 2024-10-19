### Boas Pucker ###
### b.pucker@tu-bs.de ###
### Shakunthala Natarajan ###
### v3.1.0 ###

__usage__ = """
					python3 MGSE3.1.py
					--cov <FULL_PATH_TO_COVERAGE_FILE_OR_FOLDER>| --bam <FULL_PATH_TO_BAM_FILE> 
					| --fasta <FULL_PATH_TO_UNCOMPRESSED_FASTA_FILE> --fastq <FULL_PATH_TO_SINGLE_FASTQ_FILE_OR_FOLDER_OF_FASTQ_FILES_ENDING_WITH/> <PLACE MULTIPLE FASTQ_FILES_IN_FOLDER> --gzip_fq <SPECIFY_THIS_OPTION_IF_YOU_HAVE_COMPRESED_FASTQ_FILE(S)>
					--ngs <Illumina> or <ONT> or <PacBio>
					--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
					--ref | --gff <FULL_PATH_TO_REF_GENE_FILE_OR_GFF3_FILE> | --busco <FULL_PATH_TO 'full_table_busco_run.tsv'> | --all <ALL_POS_USED_FOR_CALCULATION>
					
					optional:
					--black <FULL_PATH_TO_FILE_WITH_CONTIG_NAMES_FOR_EXCLUSION>
					--gzip <ACITVATES_SEARCH_FOR_COMPRESSED_COVERAGE_FILES>
					--bam_is_sorted PREVENTS_SORTING_OF_BAM_FILE
					--samtools <FULL_PATH_TO_SAMTOOLS>
					--bedtools <FULL_PATH_TO_BEDTOOLS>
					--short_read_aligner <FULL_PATH_TO_BWA_MEM_SHORT_READ_ALIGNER>
					--long_read_aligner <FULL_PATH_TO_MINIMAP2_LONG_READ_ALIGNER>
					--name <NAME_OF_CURRENT_ANALYSIS>
					--feature <specify the specific feature (if other than gene) you want to choose from the GFF file
					--m <SAMTOOLS_MEMORY>[5000000000]
					--threads <SAMTOOLS_THREADS>[4]
					--plot <ACTIVATE_OR_DEACTIVATE_PLOTTING TRUE|FALSE>[FALSE]
					--black_list_factor <SPECIFY THE BLACKLIST FACTOR FOR BLACKLISTING CONTIGS><FLOAT_OR_INTEGER>
					--blackoff <ACTIVATE_OR_DEACTIVATE_PLOTTING TRUE|FALSE>[FALSE]
					
					WARNING: use of absolute paths is required
					WARNING: high coverage contigs are black listed by default. Use --blackoff to disable black listing.
					
					bug reports and feature requests: b.pucker@tu-bs.de					
					"""

import os, glob, sys, subprocess
try:	#import of re is optional to avoid unnessecary dependencies
	import re
except ImportError:
	pass
try:
	import matplotlib.pyplot as plt	#import of matplotlib is optional to avoid unnessecary dependencies
except ImportError:
	pass

try:
	import gzip
except ImportError:
	pass

# --- end of imports --- #

def get_mean( values ):
	"""! @brief calculate mean of given list of values """
	
	if values != []:
		return sum( values ) / float( len( values ) )
	else:
		return 0


def get_median( values ):
	"""! @brief calculate median of given list of values """
	values.sort()
	if len( values ) >= 1:
		if len( values ) % 2 == 0:
			return ( values[ int( ( len( values ) / 2 ) -1 ) ] + values[ int( len( values ) / 2 ) ] ) / 2.0
		else:
			return values[ int( len( values ) / 2 ) ]
	else:
		return 0


def load_coverage( cov_file ):
	"""! @brief load coverage from given file """
	
	if cov_file[-3:].lower() == "cov" or cov_file[-3:].lower() == "txt":	#uncompressed coverage file
		coverage = {}
		with open( cov_file, "r" ) as f:
			line = f.readline()
			prev_chr = line.split('\t')[0]
			cov = []
			while line:
				parts = line.strip().split('\t')
				if parts[0] != prev_chr:
					coverage.update( { prev_chr: cov } )
					prev_chr = parts[0]
					cov = []
				cov.append( float( parts[2] ) )
				line = f.readline()
			coverage.update( { prev_chr: cov } )
		return coverage
	
	else:	#compressed coverage file
		coverage = {}
		with gzip.open( cov_file, "rb" ) as f:
			line = f.readline()
			line = line.decode('utf-8')
			prev_chr = line.split('\t')[0]
			cov = []
			while line:
				parts = line.strip().split('\t')
				if parts[0] != prev_chr:
					coverage.update( { prev_chr: cov } )
					prev_chr = parts[0]
					cov = []
				cov.append( float( parts[2] ) )
				line = f.readline()
				line = line.decode('utf-8')
			coverage.update( { prev_chr: cov } )
		return coverage


def load_gene_positions( ref_gene_file ):
	"""! @brief load positions of reference genes """
	
	genes = []
	with open( ref_gene_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				genes.append( { 'chr': parts[0], 'start': int( parts[1] ), 'end': int( parts[2] ), 'comment': parts[3] } )
			except IndexError:
				genes.append( { 'chr': parts[0], 'start': int( parts[1] ), 'end': int( parts[2] ) } )
			line = f.readline()
	return genes



def get_avg_cov( coverage, genes, output_folder, plotting_status, black_user_specified, black_status, excl_use_spcf_seq, blf):
	"""! @brief calculate average sequencing depth based on read coverage in certain reference genes """

	values = []
	print(str(black_status))
	if black_status:  # calculation with reference genes when --blackoff not in arguments
		if black_user_specified:#removing genes on blacklisted contigs when --black in arguments
			for element in excl_use_spcf_seq:
				for k in range(len(genes) - 1, -1, -1):
					if genes[k]['chr'] == element:
						del genes[k]
		else:# removing genes on blacklisted contigs when --blackoff and --black not in arguments
			excluded_seq_ref_gene = []
			excluded_seq_ref_gene = identify_cov_outlier_contigs(coverage, blf)
			for elements in excluded_seq_ref_gene:
				for i in range(len(genes) - 1, -1, -1):
					if genes[i]['chr'] == elements:
						del genes[i]
	for gene in genes:
		tmp = coverage[ gene['chr'] ][ gene['start']:gene['end'] ]
		for val in tmp:
			values.append( val )
	values.sort()
	mean = get_mean( values )
	med = get_median( values )
	try:
		if plotting_status:
			fig, ax = plt.subplots()
			ax.boxplot( values )
			ax.set_yscale('log')
			ax.set_ylabel( "read coverage depth" )
			ax.set_title( "n (positions) = " + str( len( values ) ) )
			fig.savefig( output_folder + "avg_cov.png", dpi=600 )
			plt.close("all")
	except:
		pass
	
	return med, mean


def identify_cov_outlier_contigs( coverage, blf ):
	"""! @brief identify contigs with outlier coverage """
	
	avg_values = []
	vals_per_contig = []
	for key in coverage.keys():
		coverage[key].sort()
		avg_cov = get_median( coverage[ key ] )
		avg_values.append( avg_cov )
		vals_per_contig.append( { 'key': key, 'val': avg_cov } )

	avg_values.sort()
	total_avg_cov = get_median( avg_values )

	blacklist = []
	for entry in vals_per_contig:
		if entry['val'] > blf*total_avg_cov:
			blacklist.append( entry['key'] )
	return blacklist


def summarize_coverage( coverage, blacklist, black_status, blf):
	"""! @brief calculates cummulative coverage and total length of analysed sequence """
	
	if blacklist == []:
		if black_status:
			blacklist = identify_cov_outlier_contigs( coverage, blf )
			excluded_seq=blacklist.copy()
		else:
			excluded_seq=[]


	
	sys.stdout.write( "number of contigs on blacklist: " + str( len( blacklist ) ) + str(blacklist) + "\n" )
	sys.stdout.flush()
	
	total_cov = 0
	total_len = 0
	cm_total_cov = 0
	cm_total_len = 0
	
	for key in list( coverage.keys() ):
		value = coverage[ key ]
		if key not in blacklist:
			total_cov += sum( value )
			total_len += len( value )
		cm_total_cov += sum( value )
		cm_total_len += len( value )
	
	return total_cov, total_len, cm_total_cov, cm_total_len, excluded_seq


def construct_ref_gene_file( gff, ref_gene_file, feature_type ):
	"""! @brief generate reference region file based on provided GFF3 file """
	
	# --- load all feature from given GFF3 file --- #
	genes = []
	with open( gff, "r" ) as f:
		line  = f.readline()
		while line:
			if line[0] != "#":
				parts = line.strip().split('\t')
				if parts[2] == feature_type:
					try:
						ID = re.findall( "g\d+", parts[-1] )[0]
						genes.append( { 'chr': parts[0], 'start': parts[3], 'end': parts[4], 'ID': ID } )
					except IndexError:
						genes.append( { 'chr': parts[0], 'start': parts[3], 'end': parts[4], 'ID': "n/a" } )
			line = f.readline()
	# --- construct file with reference regions --- #
	with open( ref_gene_file, "w" ) as out:
		for gene in genes:
			out.write( "\t".join( map( str, [ gene['chr'], gene['start'], gene['end'], gene['ID'] ] ) ) + '\n' )


def construct_ref_gene_file_from_busco( ref_gene_file, busco_gene_positions ):
	"""! @brief construct ref gene file from BUSCO TSV file"""
	
	sys.stdout.write( "WARNING: NO BUSCO GFF FILES DETECTED - USING FULL GENE POSITIONS FROM TSV INSTEAD (THIS IS A NEW FEATURE)\n" )
	sys.stdout.flush()
	
	# --- construct file with reference regions --- #
	with open( ref_gene_file, "w" ) as out:
		for gene in busco_gene_positions:
			out.write( "\t".join( list( map( str, [ gene['chr'], gene['start'], gene['end'], gene['ID'] ] ) ) ) + '\n' )


def load_BUSCOs( busco_file ):
	"""! @brief load BUSCO gene IDs and BUSCO gene positions """
	
	genes, busco_gene_positions = {}, []
	with open( busco_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				try:
					if parts[1] == "Complete":
						if not ":" in parts[2]:	#exclude cases where column3 does not specify a real chromosome
							genes.update( { parts[0]: None } )
							busco_gene_positions.append( { 'chr': parts[2], 'start': min( [ int( parts[3] ), int( parts[4] ) ] ), 'end': max( [ int( parts[3] ), int( parts[4] ) ] ), 'ID': parts[0] } )
				except IndexError:
					pass
			line = f.readline()
	return genes, busco_gene_positions

def mapping( input_fasta, fastq_files, input_ngs, sam_file, bam_file, samtools, bwa, minimap2, t):
	"""! @brief generate BAM file from given FASTA and FASTQ files"""
	if input_ngs == 'Illumina':
		sys.stdout.write("You have given illumina reads. So, indexing started with BWA-MEM ...\n")
		sys.stdout.flush()
		cmd = bwa + " index " + input_fasta
		p = subprocess.Popen(args=cmd, shell=True)
		p.communicate()
		sys.stdout.write("Mapping started with BWA-MEM ...\n")
		sys.stdout.flush()
		if len(fastq_files) == 2:#paired end illumina reads
			cmd = bwa + " mem " + "-t " + t + " " + input_fasta + " " + fastq_files[0] + " " + fastq_files[1] + " > " + sam_file
			p = subprocess.Popen(args=cmd, shell=True)
			p.communicate()
		elif len(fastq_files) == 1:
			cmd = bwa + " mem " + "-t " + t + " " + input_fasta + " " + fastq_files[0] + " > " + sam_file
			p = subprocess.Popen(args=cmd, shell=True)
			p.communicate()
	elif input_ngs == 'ONT':
		index_file = input_fasta.split('/')[-1] + ".mmi"
		sys.stdout.write("You have given ONT reads. So, indexing started with minimap2 ...\n")
		sys.stdout.flush()
		cmd = minimap2 + " -x map-ont -d " + index_file + " " + input_fasta
		p = subprocess.Popen(args=cmd, shell=True)
		p.communicate()
		sys.stdout.write("Mapping started with minimap2 ...\n")
		sys.stdout.flush()
		cmd = minimap2 + " -ax map-ont " + index_file + " " + fastq_files[0] + " > " + sam_file
		p = subprocess.Popen(args=cmd, shell=True)
		p.communicate()
	elif input_ngs == 'PacBio':
		index_file = input_fasta.split('/')[-1] + ".mmi"
		sys.stdout.write("You have given PacBio reads. So, indexing started with minimap2 ...\n")
		sys.stdout.flush()
		cmd = minimap2 + " -x map-pb -d " + index_file + " " + input_fasta
		p = subprocess.Popen(args=cmd, shell=True)
		p.communicate()
		sys.stdout.write("Mapping started with minimap2 ...\n")
		sys.stdout.flush()
		cmd = minimap2 + " -ax map-pb " + index_file + " " + fastq_files[0] + " > " + sam_file
		p = subprocess.Popen(args=cmd, shell=True)
		p.communicate()
	else:
		sys.stdout.write("You have not specified your NGS technique correctly. Exiting ...\n")
		sys.stdout.flush()
		sys.exit(__usage__)
	sys.stdout.write("Converting SAM file to BAM file ...\n")
	sys.stdout.flush()
	cmd = samtools + " view -S -b --threads " + t + " " + sam_file + " > " + bam_file
	p = subprocess.Popen(args=cmd, shell=True)
	p.communicate()

def construct_cov_file( bam_file, sorted_bam_file, cov_file, m, t, samtools, bedtools, sorting ):
	"""! @brief generate COV file for given BAM file """
	
	if sorting:			
		sys.stdout.write( "sorting BAM file ...\n" )
		sys.stdout.flush()
		cmd = samtools + " sort -m " + m + " --threads " + t + " " + bam_file + " > " + sorted_bam_file
		p = subprocess.Popen( args= cmd, shell=True )
		p.communicate()
	
	# --- calculate read coverage depth per position --- #
	sys.stdout.write( "calculating coverage per position ....\n" )
	sys.stdout.flush()
	cmd = bedtools + " -d -split -ibam " + sorted_bam_file + " > " + cov_file
	p = subprocess.Popen( args= cmd, shell=True )
	p.communicate()


def main( arguments ):
	"""! @brief runs all parts of this script """
	
	prefix = arguments[ arguments.index( '--out' )+1 ]
	if prefix[-1] != "/":
		prefix += "/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	if '--name' in arguments:	#set name for current analysis
		prefix = prefix + arguments[ arguments.index( '--name' )+1 ]
	
	# --- if coverage file or folder with coverage files is provided --- #
	if '--cov' in arguments:
		input_cov = arguments[ arguments.index( '--cov' )+1 ]	#can be file or directory with files!
		if input_cov[-1] == "/":
			cov_files = sorted( glob.glob( input_cov + "*.cov" ) )
			if '--gzip' in arguments:
				cov_files += sorted( glob.glob( input_cov + "*.cov.gz" ) )
		else:
			cov_files = [ input_cov ]
	
	# --- if BAM file or folder of BAM files is provided --- #
	elif '--bam' in arguments:
		input_bam = arguments[ arguments.index( '--bam' )+1 ]	#can be file or directory with files!
		if input_bam[-1] == "/":
			bam_files = sorted( glob.glob( input_bam + "*.bam" ) )
		else:
			bam_files = [ input_bam ]
		
		
		if '--m' in arguments:
			m = arguments[ arguments.index( '--m' )+1 ]
		else:
			m = "5000000000"
		
		if '--threads' in arguments:
			t = arguments[ arguments.index( '--threads' )+1 ]
		else:
			t = "4"
		
		if '--samtools' in arguments:
			samtools = arguments[ arguments.index( '--samtools' )+1 ]
		else:
			samtools = "samtools"
		
		if '--bedtools' in arguments:
			bedtools = arguments[ arguments.index( '--bedtools' )+1 ]
		else:
			bedtools = "genomeCoverageBed"
		
		cov_files = []
		for bam_file in bam_files:
			cov_file = prefix + bam_file.split('/')[-1] + ".cov"
			if '--bam_is_sorted' in arguments:
				construct_cov_file( bam_file, bam_file, cov_file, m, t, samtools, bedtools, False )
			else:
				sorted_bam_file = prefix + bam_file.split('/')[-1] + "_sorted.bam"
				construct_cov_file( bam_file, sorted_bam_file, cov_file, m, t, samtools, bedtools, True )
			cov_files.append( cov_file )
	# --- if FASTA file and FASTQ file(s) are provided --- #
	else:
		input_fasta = arguments[arguments.index('--fasta') + 1]
		sam_file = prefix + input_fasta.split('/')[-1] + ".sam"
		bam_file = prefix + input_fasta.split('/')[-1] + ".bam"
		input_ngs = arguments[arguments.index('--ngs') + 1]
		input_fastq = arguments[arguments.index('--fastq') + 1]
		if input_fastq[-1] == "/":
			fastq_files = sorted(glob.glob(input_fastq + "*.fq"))
			fastq_files += sorted(glob.glob(input_fastq + "*.fastq"))
			if '--gzip_fq' in arguments:
				fastq_files += sorted(glob.glob(input_fastq + "*.fq.gz"))
				fastq_files += sorted(glob.glob(input_fastq + "*.fastq.gz"))
		else:
			fastq_files = [input_fastq]
		if '--samtools' in arguments:
			samtools = arguments[ arguments.index( '--samtools' )+1 ]
		else:
			samtools = "samtools"
		if '--short_read_aligner' in arguments:
			bwa = arguments[ arguments.index( '--short_read_aligner' )+1 ]
		else:
			bwa = "bwa"
		if '--long_read_aligner' in arguments:
			minimap2 = arguments[arguments.index('--long_read_aligner') + 1]
		else:
			minimap2 = "minimap2"

		if '--m' in arguments:
			m = arguments[arguments.index('--m') + 1]
		else:
			m = "5000000000"

		if '--threads' in arguments:
			t = arguments[arguments.index('--threads') + 1]
		else:
			t = "4"

		if '--bedtools' in arguments:
			bedtools = arguments[arguments.index('--bedtools') + 1]
		else:
			bedtools = "genomeCoverageBed"

		mapping(input_fasta, fastq_files, input_ngs, sam_file, bam_file, samtools, bwa, minimap2, t)

		cov_files = []
		cov_file = prefix + bam_file.split('/')[-1] + ".cov"
		if '--bam_is_sorted' in arguments:
			construct_cov_file(bam_file, bam_file, cov_file, m, t, samtools, bedtools, False)
		else:
			sorted_bam_file = prefix + bam_file.split('/')[-1] + "_sorted.bam"
			construct_cov_file(bam_file, sorted_bam_file, cov_file, m, t, samtools, bedtools, True)
		cov_files.append(cov_file)
	
	# ---- collect remaining MGSE options --- #
	if '--feature' in arguments:
		feature_type = arguments[ arguments.index( '--feature' )+1 ]
	else:
		feature_type = "gene"
	
	use_ref_genes = True
	if '--ref' in arguments:
		ref_gene_file = arguments[ arguments.index( '--ref' ) + 1 ]
	else:
		ref_gene_file = prefix + "ref_genes.txt"
		if '--gff' in arguments:
			gff = arguments[ arguments.index( '--gff' ) + 1 ]
			construct_ref_gene_file( gff, ref_gene_file, feature_type )
		elif '--busco' in arguments:
			busco_file = arguments[ arguments.index( '--busco' ) + 1 ]
			busco_genes, busco_gene_positions = load_BUSCOs( busco_file )
			construct_ref_gene_file_from_busco( ref_gene_file, busco_gene_positions )
		elif "--all" in arguments:
			use_ref_genes = False
		else:
			sys.exit( __usage__ )
	if '--black_list_factor' in arguments:
		blf = arguments[arguments.index('--black_list_factor') + 1]
		blf = float(blf)
	else:
		blf=1.5
	if '--black' in arguments:
		black_list_file = arguments[ arguments.index( '--black' )+1 ]
		blacklist = []
		black_status = True
		black_user_specified = True
		with open( black_list_file, "r" ) as f:
			line = f.readline()
			while line:
				blacklist.append( line.strip() )
				line = f.readline()
		excl_use_spcf_seq=blacklist.copy()
	else:
		excl_use_spcf_seq = []
		blacklist = []
		black_user_specified = False
		if '--blackoff' in arguments:	#disable black listing of contigs with high coverage values
			black_status = False
			excl_use_spcf_seq = []
		else:
			black_status = True
			excl_use_spcf_seq = []
	
	if '--plot' in arguments:
		plotting_status = arguments[ arguments.index( '--plot' )+1 ]
	else:
		plotting_status = False
	
	report_file = prefix + "report.txt"	#integrate date and time of run
	
	# --- run analyses --- #
	if use_ref_genes:
		genes = load_gene_positions( ref_gene_file )	#loading positions of reference genes
	
	with open( report_file, "w" ) as out: 
		for cov_file in cov_files:
			out.write( "processing: " + cov_file + "\n" )
			
			output_folder = prefix + cov_file.split('/')[-1].split('.')[0] + "/"
			if not os.path.exists( output_folder ):
				os.makedirs( output_folder )
			
			coverage = load_coverage( cov_file )
			total_cov, total_len, cm_total_cov, cm_total_len, excluded_seq = summarize_coverage(coverage, blacklist, black_status, blf)


			if not use_ref_genes:# calculation based on all positions
				if black_status:# calculation based on all positions when --blackoff is not given
					if '--black' not in arguments:# calculation based on all positions when --black not in arguments
						for seq in excluded_seq:
							coverage.pop(seq,None)
						values = list(sorted([x for chromosome in coverage.values() for x in chromosome]))
						avg_coverage_median = get_median(values)
						avg_coverage_mean = get_mean(values)
					else:# calculation based on all positions when --black in arguments
						for seqs in excl_use_spcf_seq:
							coverage.pop(seqs,None)
						values = list(sorted([x for chromosome in coverage.values() for x in chromosome]))
						avg_coverage_median = get_median(values)
						avg_coverage_mean = get_mean(values)
				else:# calculation based on all positions when --blackoff is given
					values = list(sorted([x for chromosome in coverage.values() for x in chromosome]))
					avg_coverage_median = get_median(values)
					avg_coverage_mean = get_mean(values)
			else:	#calculation based on reference genes
				avg_coverage_median, avg_coverage_mean = get_avg_cov( coverage, genes, output_folder, plotting_status, black_user_specified, black_status , excl_use_spcf_seq, blf )
			
			
			out.write( "average coverage in reference regions (median):\t" + str( avg_coverage_median ) + "\n" )
			out.write( "average coverage in reference regions (mean):\t" + str( avg_coverage_mean ) + "\n" )

			
			out.write( "total coverage (combined length of all mapped reads):\t" + str( total_cov ) + "\n" )
			out.write( "total sequence length:\t" + str( total_len ) + "\n" )
			
			try:
				est_genome_size_med = total_cov / ( avg_coverage_median*1000000.0 )
			except ZeroDivisionError:
				est_genome_size_med = 0
			try:
				est_genome_size_mean = total_cov / ( avg_coverage_mean*1000000.0 )
			except ZeroDivisionError:
				est_genome_size_mean = 0
			
			try:
				cm_est_genome_size_med = cm_total_cov / ( avg_coverage_median*1000000.0 )
			except ZeroDivisionError:
				cm_est_genome_size_med = 0
			try:
				cm_est_genome_size_mean = cm_total_cov / ( avg_coverage_mean*1000000.0 )
			except ZeroDivisionError:
				cm_est_genome_size_mean = 0
			
			out.write( "estimated genome size based on mean [Mbp]:\t" + str( est_genome_size_mean )  + "\n" )
			out.write( "estimated genome size based on median [Mbp]:\t" + str( est_genome_size_med ) + "\n" )


if '--cov' in sys.argv and '--out' in sys.argv and '--ref' in sys.argv:
	main( sys.argv )
elif '--cov' in sys.argv and '--out' in sys.argv and '--gff' in sys.argv:
	main( sys.argv )
elif '--cov' in sys.argv and '--out' in sys.argv and '--busco' in sys.argv:
	main( sys.argv )
elif '--cov' in sys.argv and '--out' in sys.argv and '--all' in sys.argv:
	main( sys.argv )
elif '--bam' in sys.argv and '--out' in sys.argv and '--ref' in sys.argv:
	main( sys.argv )
elif '--bam' in sys.argv and '--out' in sys.argv and '--gff' in sys.argv:
	main( sys.argv )
elif '--bam' in sys.argv and '--out' in sys.argv and '--busco' in sys.argv:
	main( sys.argv )
elif '--bam' in sys.argv and '--out' in sys.argv and '--all' in sys.argv:
	main( sys.argv )
elif '--fasta' in sys.argv and '--fastq' in sys.argv and '--ngs' in sys.argv and '--ref' in sys.argv:
	main( sys.argv )
elif '--fasta' in sys.argv and '--fastq' in sys.argv and '--ngs' in sys.argv and '--gff' in sys.argv:
	main( sys.argv )
elif '--fasta' in sys.argv and '--fastq' in sys.argv and '--ngs' in sys.argv and '--busco' in sys.argv:
	main( sys.argv )
elif '--fasta' in sys.argv and '--fastq' in sys.argv and '--ngs' in sys.argv and '--all' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
