### Shakunthala Natarajan ###
### Boas Pucker ###
### v0.1 ###

__usage__ = """
	python3 busco2ploidy.py --in <FULL_PATH_TO_BUSCO_FILE> --out <OUTPUT_DIR>
			"""
######### imports #########
import matplotlib.pyplot as plt
from collections import Counter
import os
import sys
##### end of imports #####

# Function to parse the BUSCO file and extract gene occurrence counts
def parse_busco_file(file_path):
	gene_list = []
	with open(file_path, 'r') as f:
		for line in f:
			if not line.startswith('#'):
				columns = line.strip().split('\t')
				if columns:  # Ensure the line is not empty
					gene_list.append(columns[0])  # First column contains the gene name
	return gene_list

# Function to plot histogram of gene occurrence frequencies
def plot_gene_frequencies(gene_list,output_dir):
	gene_counts = Counter(gene_list)  # Count occurrences of each gene
	frequency_counts = Counter(gene_counts.values())  # Count occurrences of each frequency

	x_values = sorted(frequency_counts.keys())  # Unique occurrence frequencies
	y_values = [frequency_counts[x] for x in x_values]  # Count of genes occurring those many times

	plt.bar(x_values, y_values, color='blue', alpha=0.7)
	plt.xlabel("Occurrence Frequency")
	plt.ylabel("Number of BUSCO Genes")
	plt.title("Gene Occurrence Frequency Distribution")
	plt.xticks(x_values)  # Ensure all x-axis values are labeled
	plt.grid(axis='y', linestyle='--', alpha=0.7)
	plot_file=os.path.join(output_dir,'Gene_frequency_distribution.png')
	plt.savefig(plot_file)
	return frequency_counts  # Return frequency data for further calculations

# Function to compute Feff (Effective Copy Number Factor)
def calculate_Feff(frequency_counts,output_dir):
	total_genes = sum(frequency_counts.values())  # Total number of BUSCO genes
	Feff_numerator = sum(freq * count for freq, count in frequency_counts.items())  # Σ (C_i * N_i)
	Feff = Feff_numerator / total_genes if total_genes > 0 else 1  # Avoid division by zero
	return Feff

def main(arguments):
	busco_file = arguments[arguments.index('--in') + 1]
	output_dir = arguments[arguments.index('--out') + 1]
	gene_list = parse_busco_file(busco_file)
	frequency_counts = plot_gene_frequencies(gene_list,output_dir)
	Feff = calculate_Feff(frequency_counts,output_dir)
	print(f"Pseudo ploidy number: {Feff:.2f}")
	print("All done!")

if __name__ == '__main__':

	if '--in' in sys.argv and '--out' in sys.argv:
		main(sys.argv)
	else:
		sys.exit(__usage__)

