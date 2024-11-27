import csv
import os
import glob
import re

reader_file = open("NP/protein.faa",'r')
comp_file = open('NP/protein.fasta','w')

for line in reader_file:
        line = line.replace("\n","")
        if line[0] == '>':
                line = line.replace(">","\n>")
        comp_file.write(line)
comp_file.close()

comp_file = open('NP/protein.fasta','r')
nucleocapsid_file = open('nucleocapsid.fasta','w')

for line in comp_file:
                if "capsid" in line:
                        line = line.replace("]","]\n")
                        nucleocapsid_file.write(line)

nucleocapsid_file.close()
# \[organism=.+\]\n
# \[isolate=.+\]\n 
print("Done")
