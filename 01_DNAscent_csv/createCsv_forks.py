import os
import glob
import numpy as np
import pandas as pd
from Bio import SeqIO
from itertools import product, groupby
import gffutils
import csv
import time

# record running time
start_time = time.time()

#------------------ 00 ARGUMENTS ------------------

#whole fork or fork tip?

analogue_time = 12
date = '20240405_'

#------------------ 00 DATSETS ------------------

directories_HSR_hg38 = ['/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/2023_07_31_JJ_ONT_PromethION_2-8/20230425_1010_3H_PAO32479_809afcd8/dnascent_v312/',
'/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/2023_07_31_JJ_ONT_PromethION_2-8/20230510_1054_1G_PAO32432_0bc1c999/dnascent_v312/',
'/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/2023_07_31_JJ_ONT_PromethION_2-8/20230510_1054_3D_PAO33004_ced1b521/dnascent_v312/',
'/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/2023_07_31_JJ_ONT_PromethION_2-8/20230523_1042_3F_PAO33005_6c56f96b/dnascent_v312/']


directories_HSR_ecDNA = ['/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/2023_07_31_JJ_ONT_PromethION_2-8/20230425_1010_3H_PAO32479_809afcd8/dnascent_v312_ecDNA/',
'/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/2023_07_31_JJ_ONT_PromethION_2-8/20230510_1054_1G_PAO32432_0bc1c999/dnascent_v312_ecDNA/',
'/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/2023_07_31_JJ_ONT_PromethION_2-8/20230510_1054_3D_PAO33004_ced1b521/dnascent_v312_ecDNA/',
'/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/2023_07_31_JJ_ONT_PromethION_2-8/20230523_1042_3F_PAO33005_6c56f96b/dnascent_v312_ecDNA/']


# Colo320DM cell lines
directories_DM_hg38 = ['/home/plp27/rds/rds-pfuderer_3-28XcAxRnIRE/data/2023_03_21_JJ_ONT_PromethION1-8_d6c4705fpassfail/dnascent_v312/',
'/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/2023_07_31_JJ_ONT_PromethION_2-8/20230404_1021_2F_PAO32425_75987c73/dnascent_v312/',
'/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/2023_07_31_JJ_ONT_PromethION_2-8/20230523_1042_2E_PAO32230_cc7c73b8/dnascent_v312/',
'/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/2023_07_31_JJ_ONT_PromethION_2-8/20230523_1042_1G_PAO33157_73091287/dnascent_v312/']


directories_DM_ecDNA = ['/home/plp27/rds/rds-pfuderer_3-28XcAxRnIRE/data/2023_03_21_JJ_ONT_PromethION1-8_d6c4705fpassfail/dnascent_v312_ecDNA/',
'/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/2023_07_31_JJ_ONT_PromethION_2-8/20230404_1021_2F_PAO32425_75987c73/dnascent_v312_ecDNA/',
'/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/2023_07_31_JJ_ONT_PromethION_2-8/20230523_1042_2E_PAO32230_cc7c73b8/dnascent_v312_ecDNA/',
'/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/2023_07_31_JJ_ONT_PromethION_2-8/20230523_1042_1G_PAO33157_73091287/dnascent_v312_ecDNA/']

# Colo320DM cell lines with HU
directories_DM_hg38_HU = ['/home/plp27/rds/rds-pfuderer_4-qW4cLm3smyg/data/20240326_JJ_ONT_HU_Colo320_EdU_BrdU/TRAC-2-7311/20240110_1432_1G_PAQ70212_069eaaf1/dnascent_v312/',
'/home/plp27/rds/rds-pfuderer_4-qW4cLm3smyg/data/20240326_JJ_ONT_HU_Colo320_EdU_BrdU/TRAC-2-7312/20240110_1432_2G_PAQ61880_ba3b9440/dnascent_v312/',
'/home/plp27/rds/rds-pfuderer_4-qW4cLm3smyg/data/20240326_JJ_ONT_HU_Colo320_EdU_BrdU/TRAC-2-7313/20240110_1432_3G_PAQ89405_a452f409/dnascent_v312/'
]

directories_DM_ecDNA_HU = ['/home/plp27/rds/rds-pfuderer_4-qW4cLm3smyg/data/20240326_JJ_ONT_HU_Colo320_EdU_BrdU/TRAC-2-7311/20240110_1432_1G_PAQ70212_069eaaf1/dnascent_v312_ecDNA/',
'/home/plp27/rds/rds-pfuderer_4-qW4cLm3smyg/data/20240326_JJ_ONT_HU_Colo320_EdU_BrdU/TRAC-2-7312/20240110_1432_2G_PAQ61880_ba3b9440/dnascent_v312_ecDNA/',
'/home/plp27/rds/rds-pfuderer_4-qW4cLm3smyg/data/20240326_JJ_ONT_HU_Colo320_EdU_BrdU/TRAC-2-7313/20240110_1432_3G_PAQ89405_a452f409/dnascent_v312_ecDNA/'
]


cell_lines_dict = {'directories_HSR_hg38': directories_HSR_hg38, #all_dicts
					'directories_HSR_ecDNA': directories_HSR_ecDNA,
					'directories_DM_hg38': directories_DM_hg38,
					'directories_DM_ecDNA': directories_DM_ecDNA,
					'directories_DM_hg38_HU': directories_DM_hg38_HU,
					'directories_DM_ecDNA_HU': directories_DM_ecDNA_HU}

#------------------ 01 SEQUENCE CONTEXT ------------------

# 00: genomes
hg38 = '/home/plp27/rds/rds-pfuderer-G2yvliM05O4/genomes/GCF_000001405.40_GRCh38.p14_genomic.fna'
ecDNA = '/home/plp27/rds/rds-pfuderer-G2yvliM05O4/genomes/basic_ecDNAseq_indexfile_formatted.fasta'

# 01: create genome dictionary
def read_genome(fasta_file):
	genome_dict = {}
	for record in SeqIO.parse(fasta_file, 'fasta'):
		genome_dict[record.id] = str(record.seq).upper()
	return genome_dict

hg38_dict = read_genome(hg38)
ecDNA_dict = read_genome(ecDNA)

# 02: get sequence
def get_sequence(genome, chromosome, start, end):
	seq_string = genome[chromosome][start - 1: end]
	return seq_string

# 03: count number of homopolymers, di/trinucleotide repeats
def count_repeats(sequence, repeat):
	max_repeats = 0
	for item in repeat:
		count = 1
		while item * count in sequence:
			count += 1
		count -= 1
		if count > max_repeats:
			max_repeats = count
	return max_repeats

# 04: dinucleotide repeats
bases = ['A', 'C', 'G', 'T']
dinucleotides = [''.join(base) for base in product(bases, repeat=2)]
dinucleotides = [base for base in dinucleotides if not all(char == base[0] for char in base)]
dinucleotides.remove('CG')

# 05: trinucleotide repeats
trinucleotides = [''.join(base) for base in product(bases, repeat=3)]
trinucleotides = [base for base in trinucleotides if not all(char == base[0] for char in base)]

# how to apply repeat counting:
#homo_repeat = count_repeats(sequence, bases)
#di_repeat = count_repeats(sequence, dinucleotides)
#tri_repeat = count_repeats(sequence, trinucleotides)

#------------------ 02 OVERLAPPING GENES ------------------

# !no annotation file for ecDNA available

#------------------ 03 READ IN DATA ------------------

# make empty dataframe
#ToDo: gene orientation
column_names = ['chromosome', 'fork_lb', 'fork_ub', 'readID', 'read_lb', 'read_ub', 'strand', 'fork_length', 
	'a1_length', 'a2_length', 'a2_in_a1', 'a1_in_a1', 'a1_in_a2', 'a2_in_a2',
	'stall_score', 'fork_lb_dist', 'fork_ub_dist', 'fork_speed', 'fork_orientation',
	'As', 'Cs', 'Gs', 'Ts', 'CpGs', 'ATs', 'homo_repeat', 'di_repeat', 'tri_repeat',
	#'gene_name', 'gene_length', 'gene_orientation',
	'dataset', 'cell_line']

for cell_line, dataset in cell_lines_dict.items():
	print(f'Processing cell_line: {cell_line}')
	print(f'Processing dataset: {dataset}')
	rows = [column_names]
	for fn in dataset:
		print(f'File path {fn} and treatment {cell_line}')
		# choose correct genome dict
		if 'ecDNA' in fn:
			genome_dict = ecDNA_dict
		else:
			genome_dict = hg38_dict
		leftForkCtr = 0
		rightForkCtr = 0
		stallCtr = 0
		speedCtr = 0
		ignoreReads = []
		dataset_name = fn.split('/')[-3]
		# exclude origins
		file = glob.glob(os.path.join(fn, '*origins_DNAscent_forkSense.bed'))
		with open(file[0], 'r') as f:
			for line in f:	
				if line.startswith('#'):
					continue
				else:
					splitLine = line.rstrip().split()
					ignoreReads.append(splitLine[3])
		# exclude terminations
		file = glob.glob(os.path.join(fn, '*terminations_DNAscent_forkSense.bed'))
		with open(file[0], 'r') as f:
			for line in f:	
				if line.startswith('#'):
					continue
				else:
					splitLine = line.rstrip().split()
					ignoreReads.append(splitLine[3])
		# read in left forks
		file = glob.glob(os.path.join(fn, '*leftForks_DNAscent_forkSense_stressSignatures.bed'))
		with open(file[0], 'r') as f:
			for line in f:	
				if line.startswith('#'):
					continue
				else:
					splitLine = line.rstrip().split()
					#print(splitLine)
					chromosome = splitLine[0] #chromosome name
					fork_lb = int(splitLine[1]) #5’ boundary of the fork
					fork_ub = int(splitLine[2]) #3’ boundary of the fork
					readID = splitLine[3] # readID
					read_lb = int(splitLine[4]) #5’ boundary of the mapped read
					read_ub = int(splitLine[5]) #3’ boundary of the mapped read
					strand = splitLine[6] #strand mapped to (fwd or rev)
					fork_length = float(splitLine[7]) # fork track length [bp] - ft1
					a1_length = float(splitLine[8]) #first analogue segment length [bp] - ft2
					a2_length = float(splitLine[9]) #first analogue segment length [bp] - ft3
					a2_in_a1 = float(splitLine[10]) #second analogue segment length [bp] - ft4
					a1_in_a1 = float(splitLine[11]) #2nd in 1st - ft5
					a1_in_a2 = float(splitLine[12]) #1st in 1st - ft6
					a2_in_a2 = float(splitLine[13]) #2nd in 2nd - ft7
					stress = float(splitLine[14]) #stall score - ft8
					# add distance to end of mapped read info
					fork_lb_dist = fork_lb - read_lb
					fork_ub_dist = read_ub - fork_ub
					orientation = 'left'
					# ignore origin/termination readIDs
					if readID in ignoreReads:
						continue
					# for speed: only consider forks with 2kb distance to end of mapped read
					if fork_lb_dist >= 2000 and fork_ub_dist >= 2000:
						speed = ((fork_ub - fork_lb)/analogue_time)/1000 #/1000 for kb/min
						speedCtr += 1
					else:
						speed = 'NA'
					if stress >= 0:
						stallCtr += 1
					# 01: fork (tip) sequence context
					region_start = fork_lb #fork_lb - 2000
					region_end = fork_ub #fork_lb + 2000
					sequence = get_sequence(genome_dict, chromosome, region_start, region_end)
					# 02: normalise per region (if whole fork region_start = fork_lb, region_end = fork_ub)
					As = float(sequence.count('A'))/(region_end-region_start)
					Cs = float(sequence.count('C'))/(region_end-region_start)
					Gs = float(sequence.count('G'))/(region_end-region_start)
					Ts = float(sequence.count('T'))/(region_end-region_start)
					CpGs = float(sequence.count('CG'))/(region_end-region_start)
					ATs = float(sequence.count('AT'))/(region_end-region_start)
					homo_repeat = float(count_repeats(sequence, bases))/(region_end-region_start)
					di_repeat = float(count_repeats(sequence, dinucleotides))/(region_end-region_start)
					tri_repeat = float(count_repeats(sequence, trinucleotides))/(region_end-region_start)
					# 03: overlapping genes
					#gene = overlapping_gene_length(conversion(chromosome), fork_lb, fork_ub) #list with [name, length]
					# dataset & cell line
					row_data = [chromosome, fork_lb, fork_ub, readID, read_lb, read_ub, strand, fork_length,
						a1_length, a2_length, a2_in_a1, a1_in_a1, a1_in_a2, a2_in_a2,
						stress, fork_lb_dist, fork_ub_dist, speed, orientation,
						As, Cs, Gs, Ts, CpGs, ATs, homo_repeat, di_repeat, tri_repeat,
						dataset_name, cell_line]
				# add to list
				rows.append(row_data)
				leftForkCtr += 1
			print(f'Dataset: {fn} , left forks {leftForkCtr}')
		# read in right forks
		file = glob.glob(os.path.join(fn, '*rightForks_DNAscent_forkSense_stressSignatures.bed'))
		with open(file[0], 'r') as f:
			for line in f:	
				if line.startswith('#'):
					continue
				else:
					splitLine = line.rstrip().split()
					#print(splitLine)
					chromosome = splitLine[0] #chromosome name
					fork_lb = int(splitLine[1]) #5’ boundary of the fork
					fork_ub = int(splitLine[2]) #3’ boundary of the fork
					readID = splitLine[3] # readID
					read_lb = int(splitLine[4]) #5’ boundary of the mapped read
					read_ub = int(splitLine[5]) #3’ boundary of the mapped read
					strand = splitLine[6] #strand mapped to (fwd or rev)
					fork_length = float(splitLine[7]) # fork track length [bp] - ft1
					a1_length = float(splitLine[8]) #first analogue segment length [bp] - ft2
					a2_length = float(splitLine[9]) #first analogue segment length [bp] - ft3
					a2_in_a1 = float(splitLine[10]) #second analogue segment length [bp] - ft4
					a1_in_a1 = float(splitLine[11]) #2nd in 1st - ft5
					a1_in_a2 = float(splitLine[12]) #1st in 1st - ft6
					a2_in_a2 = float(splitLine[13]) #2nd in 2nd - ft7
					stress = float(splitLine[14]) #stall score - ft8
					# add distance to end of mapped read info
					fork_lb_dist = fork_lb - read_lb
					fork_ub_dist = read_ub - fork_ub
					orientation = 'right'
					# ignore origin/termination readIDs
					if readID in ignoreReads:
						continue
					# for speed: only consider forks with 2kb distance to end of mapped read
					if fork_lb_dist >= 2000 and fork_ub_dist >= 2000:
						speed = ((fork_ub - fork_lb)/analogue_time)/1000 #/1000 for kb/min
						speedCtr += 1
					else:
						speed = 'NA'
					if stress >= 0:
						stallCtr += 1
					# 01: fork (tip) sequence context
					region_start = fork_lb #fork_ub - 2000
					region_end = fork_ub #fork_ub + 2000
					sequence = get_sequence(genome_dict, chromosome, region_start, region_end)
					# 02: normalise per region (if whole fork region_start = fork_lb, region_end = fork_ub)
					As = float(sequence.count('A'))/(region_end-region_start)
					Cs = float(sequence.count('C'))/(region_end-region_start)
					Gs = float(sequence.count('G'))/(region_end-region_start)
					Ts = float(sequence.count('T'))/(region_end-region_start)
					CpGs = float(sequence.count('CG'))/(region_end-region_start)
					ATs = float(sequence.count('AT'))/(region_end-region_start)
					homo_repeat = float(count_repeats(sequence, bases))/(region_end-region_start)
					di_repeat = float(count_repeats(sequence, dinucleotides))/(region_end-region_start)
					tri_repeat = float(count_repeats(sequence, trinucleotides))/(region_end-region_start)
					# 03: overlapping genes
					#gene = overlapping_gene_length(conversion(chromosome), fork_lb, fork_ub) #list with [name, length]
					# dataset & cell line
					row_data = [chromosome, fork_lb, fork_ub, readID, read_lb, read_ub, strand, fork_length,
						a1_length, a2_length, a2_in_a1, a1_in_a1, a1_in_a2, a2_in_a2,
						stress, fork_lb_dist, fork_ub_dist, speed, orientation,
						As, Cs, Gs, Ts, CpGs, ATs, homo_repeat, di_repeat, tri_repeat,
						dataset_name, cell_line]
				# add to list
				rows.append(row_data)
				rightForkCtr += 1
			print(f'Dataset: {fn} , right forks {rightForkCtr}')
			print(f'Dataset: {fn} , speedCtr {speedCtr}')
			print(f'Dataset: {fn} , stallCtr {stallCtr}')
	# save
	outfile_name = '/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/plots/Jedrek_csv/' + '20240405/' + date + str(cell_line) + '.csv'
	with open(outfile_name, 'w', newline = '') as outfile:
		writer = csv.writer(outfile)
		for row in rows:
			writer.writerow(row)        
	# progress
	print(f'Done with {cell_line}; saved under {outfile_name}')

# record running time
end_time = time.time()

# end
print(f'Done. Execution time: {round((end_time - start_time) / 60, 1)} minutes')
