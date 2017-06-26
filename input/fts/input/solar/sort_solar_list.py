# David R. Thompson
import os, sys

infile = open('solar_dc_20090102.101','r')
outfile = open('solar_dc_20090102.101.sorted','w')

linelist = []
for line in infile.readlines():
    freq = float(line[3:15])
    linelist.append((freq,line))
infile.close()

for freq,line in sorted(linelist):
    outfile.write(line)
outfile.close()
