import csv
import os
import glob
import re

reader_file = open("PB2/protein.faa",'r')
comp_file = open('PB2.fasta','w')

for line in reader_file:
        line = line.replace("\n","")
        if line[0] == '>':
                line = line.replace(">","\n>")
        comp_file.write(line)
comp_file.close()

comp_file = open('ncbi_dataset/data/protein.fasta','r')
pb2_file = open('ncbi_dataset/data/pb2.fasta','w')

for line in comp_file:
                if "polymerase" in line:
                        line = line.replace("]","]\n")
                        pb2_file.write(line)

pb2_file.close()
# \[organism=.+\]\n
# \[isolate=.+\]\n 
print("Done")
