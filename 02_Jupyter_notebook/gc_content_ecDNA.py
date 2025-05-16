import csv
from Bio import SeqIO

def calculate_gc_content(sequence):
	gc_count = sequence.count('G') + sequence.count('C')
	return gc_count / len(sequence) * 100

def generate_gc_content_csv(fasta_file, window_sizes, output_prefix):
	for window_size in window_sizes:
		output_file = f"{output_prefix}_gc_content_{window_size}_bp.csv"
		with open(output_file, 'w', newline='') as csvfile:
			fieldnames = ['chromosome', 'start', 'end', 'GC_content']
			writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
			writer.writeheader()
			for record in SeqIO.parse(fasta_file, "fasta"):
				chromosome = record.id
				sequence = str(record.seq).upper()
				seq_len = len(sequence)
				for start in range(0, seq_len, window_size):
					end = min(start + window_size - 1, seq_len - 1)
					window_seq = sequence[start : end + 1]
					gc_content = calculate_gc_content(window_seq)
					writer.writerow({'chromosome': chromosome, 'start': start, 'end': end, 'GC_content': gc_content})

# specify the fasta file path and output prefix
fasta_file = "/Users/pfuderer/Documents/PhD/genomes/basic_ecDNAseq_indexfile_formatted.fasta"
window_sizes = [100, 1000, 2000, 5000, 10000, 20000, 50000]
output_prefix = "/Users/pfuderer/Documents/PhD/output/DNAscent/Jedrek_csv/20240529_G4Hunter/20240529"

# generate GC content CSV files for each window size
generate_gc_content_csv(fasta_file, window_sizes, output_prefix)
