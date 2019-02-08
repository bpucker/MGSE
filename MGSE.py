### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python MGSE.py
					--cov <FULL_PATH_TO_COVERAGE_FILE_OR_FOLDER>
					--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
					--ref | --gff <FULL_PATH_TO_REF_GENE_FILE_OR_GFF3_FILE>
					
					optional:
					--black <FULL_PATh_TO_FILE_WITH_CONTIG_NAMES_FOR_EXCLUSION>
					
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
					except:
						genes.append( { 'chr': parts[0], 'start': parts[3], 'end': parts[4], 'ID': "n/a" } )
			line = f.readline()
	
	# --- construct file with reference regions --- #
	with open( ref_gene_file, "w" ) as out:
		for gene in genes:
			out.write( "\t".join( map( str, [ gene['chr'], gene['start'], gene['end'], gene['ID'] ] ) ) + '\n' )


def main( arguments ):
	"""! @brief runs all parts of this script """
	
	input_cov = arguments[ arguments.index( '--cov' )+1 ]	#can be file or directory with files!
	if input_cov[-1] == "/":
		cov_files = sorted( glob.glob( input_cov + "*.cov" ) )
	else:
		cov_files = [ input_cov ]
	
	prefix = arguments[ arguments.index( '--out' )+1 ]
	if prefix[-1] != "/":
		prefix += "/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	if '--ref' in arguments:
		ref_gene_file = arguments[ arguments.index( '--ref' ) + 1 ]
	else:
		ref_gene_file = prefix + "ref_genes.txt"
		gff = arguments[ arguments.index( '--gff' ) + 1 ]
		construct_ref_gene_file( gff, ref_gene_file )
	
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
	genes = load_gene_positions( ref_gene_file )	#loading positions of reference genes
	
	with open( report_file, "w" ) as out: 
		for cov_file in cov_files:
			out.write( "processing: " + cov_file + "\n" )
			
			output_folder = prefix + cov_file.split('/')[-1].split('.')[0] + "/"
			if not os.path.exists( output_folder ):
				os.makedirs( output_folder )
			
			coverage = load_coverage( cov_file )
			
			avg_coverage_median, avg_coverage_mean = get_avg_cov( coverage, genes, output_folder )
			out.write( "average coverage in reference regions (median):\t" + str( avg_coverage_median ) + "\n" )
			out.write( "average coverage in reference regions (mean):\t" + str( avg_coverage_mean ) + "\n" )
			
			total_cov, total_len, cm_total_cov, cm_total_len = summarize_coverage( coverage, blacklist )
			
			out.write( "total coverage (combined length of all mapped reads):\t" + str( total_cov ) + "\n" )
			out.write( "total sequence length:\t" + str( total_len ) + "\n" )
			
			est_genome_size_med = total_cov / ( avg_coverage_median*1000000.0 )
			est_genome_size_mean = total_cov / ( avg_coverage_mean*1000000.0 )
			
			cm_est_genome_size_med = cm_total_cov / ( avg_coverage_median*1000000.0 )
			cm_est_genome_size_mean = cm_total_cov / ( avg_coverage_mean*1000000.0 )
			
			out.write( "estimated genome size based on mean [Mbp]:\t" + str( est_genome_size_mean )  + "\n" )
			out.write( "estimated genome size based on median [Mbp]:\t" + str( est_genome_size_med ) + "\n" )
			
			out.write( "estimated genome size based on mean [Mbp] (+chloro +mito):\t" + str( cm_est_genome_size_mean ) + "\n" )
			out.write( "estimated genome size based on median [Mbp] (+chloro +mito):\t" + str( cm_est_genome_size_med ) + "\n" )


if '--cov' in sys.argv and '--out' in sys.argv and '--ref' in sys.argv:
	main( sys.argv )
elif '--cov' in sys.argv and '--out' in sys.argv and '--gff' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
