# This script takes a tab-delimited txt file from StringTie
# with gene abundances and bins the expression (TPM) into bins of
# approximately 50,000 bp to cover the whole chromosome for chr1-5 of
# of the newly made assembly for Babesia duncani, taking the average TPM
# of the genes in that bin.
# Steven Abel, April 2022

from sys import argv

if len(argv) != 3:
	print("Run with 2 arguments: input_GeneAbundances.txt output_bins_avg_TPM.txt")
	exit(0)
else:
	pass

script, input_filename, output_filename = argv

current_file = open(input_filename)
write_page = open(output_filename, 'w')

chr1_binsize = (3133627/63); chr1_bin_number = 63
chr2_binsize = (1587141/32); chr2_bin_number = 32
chr3_binsize = (1428345/29); chr3_bin_number = 29
chr4_binsize = (1068143/21); chr4_bin_number = 21
chr5_binsize = (353949/7); chr5_bin_number = 7

# Create dictionaries to store the bin ranges, number of genes, and total TPM.
# Dict value list: (start of bin, end of bin, number of genes, total TPM, average TPM).
chr1_bins_genenumbers_expr = dict(); chr2_bins_genenumbers_expr = dict()
chr3_bins_genenumbers_expr = dict(); chr4_bins_genenumbers_expr = dict()
chr5_bins_genenumbers_expr = dict()
#chr1_bins_genenumbers_expr[0] = (0,chr1_binsize,0,0)
#chr1_bins_genenumbers_expr[1] = (chr1_binsize + 1,chr1_binsize * 2,0,0)
#chr1_bins_genenumbers_expr[2] = ((chr1_binsize * 2) + 1, chr1_binsize * 3,0,0)
for i in range(0,chr1_bin_number):
	chr1_bins_genenumbers_expr[i] = list()
	chr1_bins_genenumbers_expr[i].append((chr1_binsize * i) + 1)
	chr1_bins_genenumbers_expr[i].append(chr1_binsize * (i + 1))
	for j in range(0,3):
		chr1_bins_genenumbers_expr[i].append(0)
		
for i in range(0,chr2_bin_number):
	chr2_bins_genenumbers_expr[i] = list()
	chr2_bins_genenumbers_expr[i].append((chr2_binsize * i) + 1)
	chr2_bins_genenumbers_expr[i].append(chr2_binsize * (i + 1))
	for j in range(0,3):
		chr2_bins_genenumbers_expr[i].append(0)
		
for i in range(0,chr3_bin_number):
	chr3_bins_genenumbers_expr[i] = list()
	chr3_bins_genenumbers_expr[i].append((chr3_binsize * i) + 1)
	chr3_bins_genenumbers_expr[i].append(chr3_binsize * (i + 1))
	for j in range(0,3):
		chr3_bins_genenumbers_expr[i].append(0)
		
for i in range(0,chr4_bin_number):
	chr4_bins_genenumbers_expr[i] = list()
	chr4_bins_genenumbers_expr[i].append((chr4_binsize * i) + 1)
	chr4_bins_genenumbers_expr[i].append(chr4_binsize * (i + 1))
	for j in range(0,3):
		chr4_bins_genenumbers_expr[i].append(0)
		
for i in range(0,chr5_bin_number):
	chr5_bins_genenumbers_expr[i] = list()
	chr5_bins_genenumbers_expr[i].append((chr5_binsize * i) + 1)
	chr5_bins_genenumbers_expr[i].append(chr5_binsize * (i + 1))
	for j in range(0,3):
		chr5_bins_genenumbers_expr[i].append(0)

# Go through each line (gene) in input file. Check which bin it is in, and add
# the TPM of that gene into the total TPM for that bin, and increment the number
# of genes in that bin.
for line in current_file:
	line = line.rstrip()
	current_line = line.split('\t')
	if current_line[5] == "Mid":
		continue
	current_gene_ID = current_line[0]; current_gene_chr = current_line[2]
	current_gene_midcoord = float(current_line[5]); current_gene_TPM = float(current_line[9])
	current_gene_bin = 0
	
	if current_gene_chr == "chr1":
		for i in range(0,chr1_bin_number):
			if (chr1_bins_genenumbers_expr[i][0] < current_gene_midcoord < chr1_bins_genenumbers_expr[i][1]):
				current_gene_bin = i
				chr1_bins_genenumbers_expr[i][2] += 1
				chr1_bins_genenumbers_expr[i][3] += current_gene_TPM
	
	if current_gene_chr == "chr2":
		for i in range(0,chr2_bin_number):
			if (chr2_bins_genenumbers_expr[i][0] < current_gene_midcoord < chr2_bins_genenumbers_expr[i][1]):
				current_gene_bin = i	
				chr2_bins_genenumbers_expr[i][2] += 1
				chr2_bins_genenumbers_expr[i][3] += current_gene_TPM
				
	if current_gene_chr == "chr3":
		for i in range(0,chr3_bin_number):
			if (chr3_bins_genenumbers_expr[i][0] < current_gene_midcoord < chr3_bins_genenumbers_expr[i][1]):
				current_gene_bin = i
				chr3_bins_genenumbers_expr[i][2] += 1
				chr3_bins_genenumbers_expr[i][3] += current_gene_TPM
				
	if current_gene_chr == "chr4":
		for i in range(0,chr4_bin_number):
			if (chr4_bins_genenumbers_expr[i][0] < current_gene_midcoord < chr4_bins_genenumbers_expr[i][1]):
				current_gene_bin = i
				chr4_bins_genenumbers_expr[i][2] += 1
				chr4_bins_genenumbers_expr[i][3] += current_gene_TPM
				
	if current_gene_chr == "chr5":
		for i in range(0,chr5_bin_number):
			if (chr5_bins_genenumbers_expr[i][0] < current_gene_midcoord < chr5_bins_genenumbers_expr[i][1]):
				current_gene_bin = i
				chr5_bins_genenumbers_expr[i][2] += 1
				chr5_bins_genenumbers_expr[i][3] += current_gene_TPM
				
# Get the average gene expression within each bin, and 
# output the contents of the dicts to the output .txt file.
# NOTE: in case of 0 genes in bin, set average TPM to 0.
write_page.write("Chr1 bin number and average TPM\n")

for bin_number in chr1_bins_genenumbers_expr.keys():
	if chr1_bins_genenumbers_expr[bin_number][2] != 0:
		chr1_bins_genenumbers_expr[bin_number][4] = chr1_bins_genenumbers_expr[bin_number][3] / chr1_bins_genenumbers_expr[bin_number][2]
	else:
		chr1_bins_genenumbers_expr[bin_number][4] = 0
	write_page.write(f"{bin_number}\t{chr1_bins_genenumbers_expr[bin_number][4]}\n")
	
write_page.write("Chr2 bin number and average TPM\n")	
	
for bin_number in chr2_bins_genenumbers_expr.keys():
	if chr2_bins_genenumbers_expr[bin_number][2] != 0:
		chr2_bins_genenumbers_expr[bin_number][4] = chr2_bins_genenumbers_expr[bin_number][3] / chr2_bins_genenumbers_expr[bin_number][2]
	else:
		chr2_bins_genenumbers_expr[bin_number][4] = 0
	write_page.write(f"{bin_number}\t{chr2_bins_genenumbers_expr[bin_number][4]}\n")

write_page.write("Chr3 bin number and average TPM\n")
	
for bin_number in chr3_bins_genenumbers_expr.keys():
	if chr3_bins_genenumbers_expr[bin_number][2] != 0:
		chr3_bins_genenumbers_expr[bin_number][4] = chr3_bins_genenumbers_expr[bin_number][3] / chr3_bins_genenumbers_expr[bin_number][2]
	else:
		chr3_bins_genenumbers_expr[bin_number][4] = 0
	write_page.write(f"{bin_number}\t{chr3_bins_genenumbers_expr[bin_number][4]}\n")

write_page.write("Chr4 bin number and average TPM\n")	
	
for bin_number in chr4_bins_genenumbers_expr.keys():
	if chr4_bins_genenumbers_expr[bin_number][2] != 0:
		chr4_bins_genenumbers_expr[bin_number][4] = chr4_bins_genenumbers_expr[bin_number][3] / chr4_bins_genenumbers_expr[bin_number][2]
	else:
		chr4_bins_genenumbers_expr[bin_number][4] = 0
	write_page.write(f"{bin_number}\t{chr4_bins_genenumbers_expr[bin_number][4]}\n")

write_page.write("Chr5 bin number and average TPM\n")	
	
for bin_number in chr5_bins_genenumbers_expr.keys():
	if chr5_bins_genenumbers_expr[bin_number][2] != 0:
		chr5_bins_genenumbers_expr[bin_number][4] = chr5_bins_genenumbers_expr[bin_number][3] / chr5_bins_genenumbers_expr[bin_number][2]
	else:
		chr5_bins_genenumbers_expr[bin_number][4] = 0
	write_page.write(f"{bin_number}\t{chr5_bins_genenumbers_expr[bin_number][4]}\n")	
	
write_page.close()