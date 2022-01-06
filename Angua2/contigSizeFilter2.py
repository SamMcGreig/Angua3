#! /usr/bin/python
#USAGE contigSizeFilter.py [input].fasta [sizeCutoff]
#a script for returning all sequences in a file longer than [sizeCutoff]
#with output file sorted by decreasing sequence length

import sys
from Bio import SeqIO

#parse arguments
inFile = sys.argv[1]
sizeCutoff = sys.argv[2]
outDir = sys.argv[3]

filename = inFile.split("/")
filename = filename[2]

outFile = outDir + 'sorted_' + sys.argv[2] + '_' + filename

#Get the lengths and ids, and sort on length         
len_and_ids = sorted((len(rec), rec.id) for rec in \
                     SeqIO.parse(inFile,"fasta"))
ids = reversed([id for (length, id) in len_and_ids])
del len_and_ids #free this memory
record_index = SeqIO.index(inFile, "fasta")
records = (record_index[id] for id in ids)

#filter sequences
keptRecords = []
for seq_record in records:
	if len(seq_record.seq) >= int(sizeCutoff):
		keptRecords.append(seq_record)
#output sequences to file
SeqIO.write(keptRecords, outFile, "fasta")
