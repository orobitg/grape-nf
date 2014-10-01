#!/usr/bin/env python

import os
import sys
import shutil
import fileinput

annotation_file = sys.argv[1]

f = open(annotation_file, 'r')
lines = f.read()
f.close()
if lines.find('transcript_id') == -1:
	print "ERROR - GTF/GFF file: No transcript_id descrition. Flux Capacitor in Grape will crash. Correct annotation file" 
    
if lines.find('exon\t') == -1:

	tmp_file = annotation_file + '.tmp'
	fout = open(tmp_file,'w')
	for line in fileinput.input(annotation_file, inplace=1):
	    if 'cds' in line:
		line = line.replace('cds\t', 'exon\t')
	    elif 'CDS' in line:
		line = line.replace('cds\t', 'exon\t')
	    fout.write(line)
	fout.close()

	shutil.move(tmp_file, annotation_file)
