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
date = '20240503_'

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

#------------------ 03 READ IN DATA ------------------

# make empty dataframe
#ToDo: gene orientation
column_names = ['chromosome', 'ori_lb', 'ori_ub', 'readID', 'read_lb', 'read_ub', 'strand',
	'dataset', 'cell_line']

for cell_line, dataset in cell_lines_dict.items():
	print(f'Processing cell_line: {cell_line}')
	print(f'Processing dataset: {dataset}')
	rows = [column_names]
	for fn in dataset:
		print(f'File path {fn} and treatment {cell_line}')
		dataset_name = fn.split('/')[-3]
		# read in origins
		oriCtr = 0
		file = glob.glob(os.path.join(fn, '*origins_DNAscent_forkSense.bed'))
		with open(file[0], 'r') as f:
			for line in f:	
				if line.startswith('#'):
					continue
				else:
					splitLine = line.rstrip().split()
					#print(splitLine)
					chromosome = splitLine[0] #chromosome name
					ori_lb = int(splitLine[1]) #5’ boundary of the fork
					ori_ub = int(splitLine[2]) #3’ boundary of the fork
					readID = splitLine[3] # readID
					read_lb = int(splitLine[4]) #5’ boundary of the mapped read
					read_ub = int(splitLine[5]) #3’ boundary of the mapped read
					strand = splitLine[6] #strand mapped to (fwd or rev)
					# dataset & cell line
					row_data = [chromosome, ori_lb, ori_ub, readID, read_lb, read_ub, strand,
						dataset_name, cell_line]
				# add to list
				rows.append(row_data)
				oriCtr += 1
			print(f'Dataset: {fn} , origins {oriCtr}')
	# save
	outfile_name = '/home/plp27/rds/rds-pfuderer_2-OvX291IjOjA/data/plots/Jedrek_csv/' + '20240503/' + date + str(cell_line) + '_origins.csv'
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
