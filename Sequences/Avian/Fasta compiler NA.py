import csv
import os
import glob
import re

reader_file = open("06NA/protein.faa",'r')
comp_file = open('06NA/protein.fasta','w')

for line in reader_file:
        line = line.replace("\n","")
        if line[0] == '>':
                line = line.replace(">","\n>")
        comp_file.write(line)
comp_file.close()

comp_file = open('06NA/protein.fasta','r')
neuraminidase_file = open('na.fasta','w')
leftovers_file = open('06NA/leftovers.fasta','w')

for line in comp_file:
                if "neuraminidase" in line or "Neuraminidase" in line or "NA" in line or "neuramindase" in line or "neuraminadase" in line:
                        line = line.replace("]","]\n")
                        neuraminidase_file.write(line)
                else:
                        line = line.replace("]","]\n")
                        leftovers_file.write(line)

leftovers_file.close()
neuraminidase_file.close()
# \[organism=.+\]\n
# \[isolate=.+\]\n 
print("Done")
