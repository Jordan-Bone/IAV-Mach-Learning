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
leftovers_file = open('04HA/leftovers.fasta','w')

for line in comp_file:
                if ("hemagglutinin" in line or "haemagglutinin" in line or "haemaglutinin" in line or "HA" in line or "hemagglutunin" in line or "hemaglutinin" in line) and not "[polyprotein" in line:
                        line = line.replace("]","]\n")
                        haemagglutinin_file.write(line)
                elif not ("sig" in line or "HA2" in line or "HA1" in line):
                        line = line.replace("]","]\n")
                        leftovers_file.write(line)

haemagglutinin_file.close()
ha1_file.close()
leftovers_file.close()
# \[organism=.+\]\n
# \[isolate=.+\]\n 
print("Done")
