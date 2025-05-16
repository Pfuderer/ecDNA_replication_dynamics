from Bio import SeqIO

# paths to the fasta file and the cytoband file
fasta_file_path = '/Users/pfuderer/Documents/PhD/genomes/basic_ecDNAseq_indexfile_formatted.fasta'
cytoband_file_path = '/Users/pfuderer/Documents/PhD/genomes/basic_ecDNAseq_indexfile_formatted_cytoband.txt'

# function to create a cytoband file from a fasta file
def create_cytoband_from_fasta(fasta_path, cytoband_path):
	with open(cytoband_path, 'w') as cytoband_file:
		for record in SeqIO.parse(fasta_path, "fasta"):
			chromosome_name = record.id
			sequence_length = len(record.seq)
			gieStain = 'gpos75'  # Use 'gpos75' as the stain
			# write to cytoband file
			cytoband_file.write(f'{chromosome_name}\t0\t{sequence_length}\t{chromosome_name}\t{gieStain}\n')

# call the function
create_cytoband_from_fasta(fasta_file_path, cytoband_file_path)
