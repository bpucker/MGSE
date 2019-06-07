### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.51 ###

__usage__ = """
					python MGSE.py
					--cov <FULL_PATH_TO_COVERAGE_FILE_OR_FOLDER>| --bam <FULL_PATH_TO_BAM_FILE>
					--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
					--ref | --gff <FULL_PATH_TO_REF_GENE_FILE_OR_GFF3_FILE> | --busco <FULL_PATH_TO 'full_table_busco_run.tsv'> | --all <ALL_POS_USED_FOR_CALCULATION>
					
					optional:
					--black <FULL_PATH_TO_FILE_WITH_CONTIG_NAMES_FOR_EXCLUSION>
					--gzip <ACITVATES_SEARCH_FOR_COMPRESSED_FILES>
					--bam_is_sorted PREVENTS_SORTING_OF_BAM_FILE
					--samtools <FULL_PATH_TO_SAMTOOLS>
					--bedtools <FULL_PATH_TO_BEDTOOLS>
					--name <NAME_OF_CURRENCT_ANALYSIS>
					--m <SAMTOOLS_MEMORY>[5000000000]
					--threads <SAMTOOLS_THREADS>[4]
										
					WARNING: if --busco is used, the BUSCO GFF3 files need to be in the default folder relative to the provided TSV file
					WARNING: use of absolute paths is required
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de					
					"""

import os, glob, sys
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
	
	if len( values ) >= 1:
		if len( values ) % 2 == 0:
			return ( values[ ( len( values ) / 2 ) -1 ] + values[ len( values ) / 2 ] ) / 2.0
		else:
			return values[ len( values ) / 2 ]
	else:
		return 0


def load_coverage( cov_file ):
	"""! @brief load coverage from given file """
	
	if cov_file[-3:].lower() == "cov":	#uncompressed coverage file
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


def get_avg_cov( coverage, genes, output_folder ):
	"""! @brief calculate average sequencing depth based on read coverage in certain reference genes """
	
	values = []
	for gene in genes:
		tmp = coverage[ gene['chr'] ][ gene['start']:gene['end'] ]
		for val in tmp:
			values.append( val )
	mean = get_mean( values )
	med = get_median( values )
	try:
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


def identify_cov_outlier_contigs( coverage ):
	"""! @brief identify contigs with outlier coverage """
	
	avg_values = []
	vals_per_contig = []
	for key in coverage.keys():
		avg_cov = get_median( coverage[ key ] )
		avg_values.append( avg_cov )
		vals_per_contig.append( { 'key': key, 'val': avg_cov } )
	total_avg_cov = get_median( avg_values )
	
	blacklist = []
	for entry in vals_per_contig:
		if entry['val'] > 1.5*total_avg_cov:
			blacklist.append( entry['key'] )
	return blacklist


def summarize_coverage( coverage, blacklist ):
	"""! @brief calculates cummulative coverage and total length of analysed sequence """
	
	if blacklist == []:
		blacklist = identify_cov_outlier_contigs( coverage )
	
	total_cov = 0
	total_len = 0
	cm_total_cov = 0
	cm_total_len = 0
	
	for key in coverage.keys():
		value = coverage[ key ]
		if key not in blacklist:
			total_cov += sum( value )
			total_len += len( value )
		cm_total_cov += sum( value )
		cm_total_len += len( value )
	
	return total_cov, total_len, cm_total_cov, cm_total_len


def construct_ref_gene_file( gff, ref_gene_file, feature_type="gene" ):
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


def construct_ref_gene_file_from_busco( gff_folder, busco_genes, ref_gene_file, feature_type ):
	"""! @brief construct ref gene file from BUSCO GFFs """
	
	gff_files = glob.glob( gff_folder + "*.gff" )
	genes = []
	for filename in gff_files:
		if filename.split('/')[-1].split('.')[0] in busco_genes:
			with open( filename, "r" ) as f:
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


def load_BUSCOs( busco_file ):
	"""! @brief load genes """
	
	genes = {}
	with open( busco_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				try:
					if parts[1] == "Complete":
						genes.update( { parts[0]: None } )
				except IndexError:
					pass
			line = f.readline()
	return genes


def construct_cov_file( bam_file, sorted_bam_file, cov_file, m, t, samtools, bedtools, sorting ):
	"""! @brief generate COV file for given BAM file """
	
	if sorting:			
		print "sorting BAM file ..."
		cmd = samtools + " sort -m " + m + " --threads " + t + " " + bam_file + " > " + sorted_bam_file
		os.popen( cmd )
	
	# --- calculate read coverage depth per position --- #
	print "calculating coverage per position ...."
	cmd = bedtools + " -d -split -ibam " + sorted_bam_file + " > " + cov_file
	os.popen( cmd )


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
	else:
		input_bam = arguments[ arguments.index( '--bam' )+1 ]	#can be file or directory with files!
		if input_bam[-1] == "/":
			bam_files = sorted( glob.glob( input_cov + "*.bam" ) )
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
			gff_folder = "/".join( busco_file.split('/')[:-1] ) + "/augustus_output/gffs/"
			busco_genes = load_BUSCOs( busco_file )
			construct_ref_gene_file_from_busco( gff_folder, busco_genes, ref_gene_file, feature_type )
		elif "--all" in arguments:
			use_ref_genes = False
		else:
			sys.exit( __usage__ )
	
	if '--black' in arguments:
		black_list_file = arguments[ arguments.index( '--black' )+1 ]
		blacklist = []
		with open( black_list_file, "r" ) as f:
			line = f.readline()
			while line:
				blacklist.append( line.strip() )
				line = f.readline()
	else:
		blacklist = []
	
	
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
			
			
			if not use_ref_genes:	#calculation based on all positions
				values = [ x for chromosome in coverage.values() for x in chromosome ]
				avg_coverage_median = get_median( values )
				avg_coverage_mean = get_mean( values )
			else:	#calculation based on reference genes
				avg_coverage_median, avg_coverage_mean = get_avg_cov( coverage, genes, output_folder )
			
			
			out.write( "average coverage in reference regions (median):\t" + str( avg_coverage_median ) + "\n" )
			out.write( "average coverage in reference regions (mean):\t" + str( avg_coverage_mean ) + "\n" )
			
			total_cov, total_len, cm_total_cov, cm_total_len = summarize_coverage( coverage, blacklist )
			
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
else:
	sys.exit( __usage__ )
