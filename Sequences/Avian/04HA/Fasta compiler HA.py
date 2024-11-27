import csv
import os
import glob
import re

reader_file = open("04HA/protein.faa",'r')
comp_file = open('04HA/protein.fasta','w')

for line in reader_file:
        line = line.replace("\n","")
        if line[0] == '>':
                line = line.replace(">","\n>")
        comp_file.write(line)
comp_file.close()

comp_file = open('04HA/protein.fasta','r')
haemagglutinin_file = open('ha.fasta','w')
ha1_file = open('ha1.fasta','w')

for line in comp_file:
                if ("hemagglutinin" in line or "haemagglutinin" in line) and not "[polyprotein" in line:
                        line = line.replace("]","]\n")
                        haemagglutinin_file.write(line)
                elif "HA1" in line:
                        line = line.replace("]","]\n")
                        ha1_file.write(line)

haemagglutinin_file.close()
ha1_file.close()
# \[organism=.+\]\n
# \[isolate=.+\]\n 
print("Done")
