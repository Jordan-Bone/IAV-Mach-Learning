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
leftovers_file = open('08NS/leftovers.fasta','w')

for line in comp_file:
        if "nonstructural protein 2" in line or "nuclear" in line or "NS2" in line or (" 2 " in line and "non" in line.lower()):
                line = line.replace("]","]\n")
##                if " [isolate" in line:
##                        continue
                nep_file.write(line)
        elif "nonstructural protein 1" in line or "NS1" in line or (" 1 " in line and "non" in line.lower()) or ":1-2" in line:
                line = line.replace("]","]\n")
                ns1_file.write(line)
        else:
                line = line.replace("]","]\n")
                leftovers_file.write(line)

nep_file.close()
ns1_file.close()
leftovers_file.close()
# \[organism=.+\]\n
# \[isolate=.+\]\n 
print("Done")
