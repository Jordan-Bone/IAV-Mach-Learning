import csv
import os
import glob
import re

reader_file = open("02PB1/protein.faa",'r')
comp_file = open('02PB1/protein.fasta','w')

for line in reader_file:
        line = line.replace("\n","")
        if line[0] == '>':
                line = line.replace(">","\n>")
        comp_file.write(line)
comp_file.close()

comp_file = open('02PB1/protein.fasta','r')
pb1_file = open('pb1.fasta','w')
pb1_f2_file = open('pb1-f2.fasta','w')

for line in comp_file:
                if "polymerase" in line and not "-F2" in line:
                        line = line.replace("]","]\n")
                        pb1_file.write(line)
                elif "-F2" in line:
                        line = line.replace("]","]\n")
                        pb1_f2_file.write(line)

pb1_file.close()
pb1_f2_file.close()
# \[organism=.+\]\n
# \[isolate=.+\]\n 
print("Done")
