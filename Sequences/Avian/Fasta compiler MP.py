import csv
import os
import glob
import re

reader_file = open("07MP/protein.faa",'r')
comp_file = open('07MP/protein.fasta','w')

for line in reader_file:
        line = line.replace("\n","")
        if line[0] == '>':
                line = line.replace(">","\n>")
        comp_file.write(line)
comp_file.close()

comp_file = open('07MP/protein.fasta','r')
m1_file = open('m1.fasta','w')
m2_file = open('m2.fasta','w')
leftovers_file = open('07MP/leftovers.fasta','w')

for line in comp_file:
                if "matrix protein 1" in line.lower() or "M1" in line or ":1-2" in line:
                        line = line.replace("]","]\n")
                        m1_file.write(line)
                elif "matrix protein 2" in line.lower() or "M2" in line or "ion" in line or ":1-9" in line:
                        line = line.replace("]","]\n")
                        m2_file.write(line)
                else:
                        line = line.replace("]","]\n")
                        leftovers_file.write(line)

m1_file.close()
m2_file.close()
leftovers_file.close()
# \[organism=.+\]\n
# \[isolate=.+\]\n 
print("Done")
