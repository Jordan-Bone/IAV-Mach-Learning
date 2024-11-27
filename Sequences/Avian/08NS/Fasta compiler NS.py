import csv
import os
import glob
import re

reader_file = open("08NS/protein.faa",'r')
comp_file = open('08NS/protein.fasta','w')

for line in reader_file:
        line = line.replace("\n","")
        if line[0] == '>':
                line = line.replace(">","\n>")
        comp_file.write(line)
comp_file.close()

comp_file = open('08NS/protein.fasta','r')
ns1_file = open('ns1.fasta','w')
nep_file = open('nep.fasta','w')

for line in comp_file:
                if "nonstructural protein 2" in line or "nuclear" in line:
                        line = line.replace("]","]\n")
                        nep_file.write(line)
                elif "nonstructural protein 1" in line:
                        line = line.replace("]","]\n")
                        ns1_file.write(line)

nep_file.close()
ns1_file.close()
# \[organism=.+\]\n
# \[isolate=.+\]\n 
print("Done")
