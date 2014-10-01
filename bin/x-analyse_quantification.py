#!/usr/bin/python

import os
import sys
import shutil
import fileinput

infile=sys.argv[1]
output=sys.argv[2]
threshold=float(sys.argv[3])

def getMedian(numericValues):
	theValues = sorted(numericValues)
	if len(theValues) % 2 == 1:
		return theValues[(len(theValues)+1)/2-1]
	else:
		lower = theValues[len(theValues)/2-1]
		upper = theValues[len(theValues)/2]
	sum = float(lower) + float(upper)
	return (sum / 2)


if not os.path.isfile(infile):
	print "ERROR - %s does no exists\n" %infile
	exit

try:
	f = open(infile, 'r')
except ValueError:
	exit

cont = 0
rpkm_list = list()
reads_list = list()
length_list = list()
trid_list = list()
rpkm_sum = 0
rpkm_avg = 0
reads_sum = 0
reads_avg = 0
ngenes = 0
ngenes_sup = 0
tnuc_genes = 0
tlength_genes = 0
trpkm_genes = 0.0

data = f.readlines()

for lines in data:
	tokens = lines.split(';')
	for i in range(0,len(tokens)):
	#Canviar comes a punts

		if tokens[i].find('RPKM') != -1:
			rpkm = tokens[i].replace(' RPKM ', '')	
			rpkm = rpkm.replace(',', '.')
			rpkm_list.append(rpkm)
			ngenes += 1
			if float(rpkm) > threshold:
				ngenes_sup += 1
				trpkm_genes += float(rpkm)

		if tokens[i].find('reads') != -1:
			reads = tokens[i].replace(' reads ', '')
			reads = reads.replace(',', '.')
			reads_list.append(reads)

		if tokens[i].find('length') != -1:
			length = tokens[i].replace(' length ', '')
			length = length.replace(',', '.')
			length_list.append(length)
			tlength_genes += int(length)
		if tokens[i].find('transcript_id') != -1:
			trid_tokens = tokens[i].split(' ')
			trid = trid_tokens[1].replace('"', '')
			if trid in trid_list:
				nuc_tokens = tokens[i].split('\t')
				nuc = int(nuc_tokens[4]) - int(nuc_tokens[3])
				tnuc_genes += nuc
			else:
				trid_list.append(trid)
				nuc_tokens = tokens[i].split('\t')
				nuc = int(nuc_tokens[4]) - int(nuc_tokens[3])
				tnuc_genes += nuc
	rpkm_sum += float(rpkm_list[cont])	
	reads_sum += float(reads_list[cont])
	cont += 1

rpkm_avg = rpkm_sum / cont
reads_avg = reads_sum / cont
rpkm_med = getMedian(rpkm_list)
reads_med = getMedian(reads_list)
f.close()

f = open(output, 'w')
f.write('#;#Total;#Average;#Median\n')
f.write('Reads;' + str(reads_sum) + ';' + str(reads_avg) + ';' + str(reads_med) + '\n')
f.write('RPKM;' + str(rpkm_sum) + ';' + str(rpkm_avg) + ';' + str(rpkm_med) + '\n')
f.write('\n')
f.write('#;#Total;#Total Found;#Total Supported\n')
f.write('N Genes;' + str(ngenes) + ';' + str(ngenes_sup) +'\n')
f.write('N Nuc Genes;' + str(tlength_genes) + ';' + str(trpkm_genes) +'\n')
f.close()
