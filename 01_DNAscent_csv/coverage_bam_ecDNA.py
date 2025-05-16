import os
import pandas as pd
import pysam
import csv

#------------------ 00 SETUP ------------------

date = '202410017'

outpath = f'/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/plots/Jedrek_csv/{date}_bam_coverage'

if not os.path.exists(outpath):
	os.makedirs(outpath)

#------------------ 01 FUNCTION ------------------

# 01: extract data from bam file
def extract_bam_data(bam_files, output_csv, condition_name):
	with open(output_csv, mode='w', newline='') as csvfile:
		csv_writer = csv.writer(csvfile)
		# csv header with chromosome column added
		csv_writer.writerow(['chromosome', 'start_position', 'end_position', 'aligned_read_length', 'mapping_quality', 'dataset_name', 'readID'])
		# iterate through each bam file path
		for bam_file in bam_files:
			dataset_name = os.path.basename(os.path.dirname(os.path.dirname(bam_file)))
			bamfile = pysam.AlignmentFile(bam_file, 'rb')
			# iterate through each read in the bam file
			for read in bamfile.fetch():
				if read.is_unmapped:  # skip unmapped reads
					continue
				read_id = read.query_name
				chromosome = read.reference_name  # get the chromosome name
				start_pos = read.reference_start + 1  # 1-based start position
				end_pos = read.reference_end  # 1-based end position
				read_length = read.reference_length
				mapping_quality = read.mapping_quality
				# write to csv with chromosome
				csv_writer.writerow([chromosome, start_pos, end_pos, read_length, mapping_quality, dataset_name, read_id])
			# close bam file
			bamfile.close()

# 02: process bam files
def process_bam_files(directories, condition_name, output_csv):
	bam_files = []
	for directory in directories:
		bam_file = os.path.join(directory, 'alignments.sorted.bam')
		if os.path.exists(bam_file):
			bam_files.append(bam_file)
		else:
			print(f'bam file not found in directory: {directory}')
	# apply function 1
	extract_bam_data(bam_files, output_csv, condition_name)

#------------------ 02 INPUT ------------------

# HSR ecDNA directories
directories_HSR_ecDNA = [
	'/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/2023_07_31_JJ_ONT_PromethION_2-8/20230425_1010_3H_PAO32479_809afcd8/reads_minimap_ecDNA/',
	'/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/2023_07_31_JJ_ONT_PromethION_2-8/20230510_1054_1G_PAO32432_0bc1c999/reads_minimap_ecDNA/',
	'/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/2023_07_31_JJ_ONT_PromethION_2-8/20230510_1054_3D_PAO33004_ced1b521/reads_minimap_ecDNA/',
	'/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/2023_07_31_JJ_ONT_PromethION_2-8/20230523_1042_3F_PAO33005_6c56f96b/reads_minimap_ecDNA/'
]

# DM ecDNA directories
directories_DM_ecDNA = [
	'/home/plp27/rds/rds-pfuderer_3-28XcAxRnIRE/data/2023_03_21_JJ_ONT_PromethION1-8_d6c4705fpassfail/reads_minimap_ecDNA/',
	'/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/2023_07_31_JJ_ONT_PromethION_2-8/20230404_1021_2F_PAO32425_75987c73/reads_minimap_ecDNA/',
	'/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/2023_07_31_JJ_ONT_PromethION_2-8/20230523_1042_2E_PAO32230_cc7c73b8/reads_minimap_ecDNA/',
	'/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/2023_07_31_JJ_ONT_PromethION_2-8/20230523_1042_1G_PAO33157_73091287/reads_minimap_ecDNA/'
]

directories_DM_HU_ecDNA = [
	'/home/plp27/rds/rds-pfuderer_4-qW4cLm3smyg/data/20240326_JJ_ONT_HU_Colo320_EdU_BrdU/TRAC-2-7311/20240110_1432_1G_PAQ70212_069eaaf1/reads_minimap_ecDNA/',
	'/home/plp27/rds/rds-pfuderer_4-qW4cLm3smyg/data/20240326_JJ_ONT_HU_Colo320_EdU_BrdU/TRAC-2-7312/20240110_1432_2G_PAQ61880_ba3b9440/reads_minimap_ecDNA/',
	'/home/plp27/rds/rds-pfuderer_4-qW4cLm3smyg/data/20240326_JJ_ONT_HU_Colo320_EdU_BrdU/TRAC-2-7313/20240110_1432_3G_PAQ89405_a452f409/reads_minimap_ecDNA/'
]

#------------------ 03 RUN ------------------

# execute functions using HSR or DM ecDNA directories

# DM
process_bam_files(directories_DM_ecDNA, 'DM', f'{outpath}/{date}_DM_ecDNA_coverage_bam.csv')

# HSR
process_bam_files(directories_HSR_ecDNA, 'HSR', f'{outpath}/{date}_HSR_ecDNA_coverage_bam.csv')

# DM_HU
process_bam_files(directories_DM_HU_ecDNA, 'DM_HU', f'{outpath}/{date}_DM_HU_ecDNA_coverage_bam.csv')

#------------------ 04 FILTERED CSV ------------------

# read in the csv files
df_DM = pd.read_csv(f'{outpath}/{date}_DM_ecDNA_coverage_bam.csv')
df_HSR = pd.read_csv(f'{outpath}/{date}_HSR_ecDNA_coverage_bam.csv')
df_DM_HU = pd.read_csv(f'{outpath}/{date}_DM_HU_ecDNA_coverage_bam.csv')

# filter out reads with mapping quality less than 20 and length < 20kb
df_DM_filtered = df_DM[(df_DM['mapping_quality'] >= 20) & (df_DM['aligned_read_length'] >= 20000)]
df_HSR_filtered = df_HSR[(df_HSR['mapping_quality'] >= 20) & (df_HSR['aligned_read_length'] >= 20000)]
df_DM_HU_filtered = df_DM_HU[(df_DM_HU['mapping_quality'] >= 20) & (df_DM_HU['aligned_read_length'] >= 20000)]

# save filtered csv files
df_DM_filtered.to_csv(f'{outpath}/{date}_DM_ecDNA_coverage_bam_filtered.csv', index = False)
df_HSR_filtered.to_csv(f'{outpath}/{date}_HSR_ecDNA_coverage_bam_filtered.csv', index = False)
df_DM_HU_filtered.to_csv(f'{outpath}/{date}_DM_HU_ecDNA_coverage_bam_filtered.csv', index = False)

#------------------ END ------------------

print('Done.')
