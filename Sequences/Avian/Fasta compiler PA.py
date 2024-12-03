import csv
import os
import glob
import re

reader_file = open("03PA/protein.faa",'r')
comp_file = open('03PA/protein.fasta','w')

for line in reader_file:
        line = line.replace("\n","")
        if line[0] == '>':
                line = line.replace(">","\n>")
        comp_file.write(line)
comp_file.close()

comp_file = open('03PA/protein.fasta','r')
pa_file = open('pa.fasta','w')
pa_x_file = open('pa-x.fasta','w')
leftovers_file = open('03PA/leftovers.fasta','w')

for line in comp_file:
                if ("polymerase" in line or "plymerase" in line or "PA" in line) and not "-X" in line:
                        line = line.replace("]","]\n")
                        pa_file.write(line)
                elif "-X" in line:
                        line = line.replace("]","]\n")
                        pa_x_file.write(line)
                else:
                        line = line.replace("]","]\n")
                        leftovers_file.write(line)


pa_file.close()
pa_x_file.close()
leftovers_file.close()
# \[organism=.+\]\n
# \[isolate=.+\]\n 
print("Done")
