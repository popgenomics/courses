#!/usr/bin/python
import fileinput
import sys

msfile = sys.argv[1] # name of the ms file
nA = int(sys.argv[2]) # number of diploids in spA
nB = int(sys.argv[3]) # number of diploids in spB

#msfile = sys.argv[1]
infile = open(msfile, "r")

res = ''
test = 0
cnt = -1
nInd = 0
#for line in fileinput.input():
for line in infile:
	line = line.strip()
	if 'segsites' in line:
		if test == 0:
			locusID=-1
		locusID += 1
		indID = -1
		segsites = int(line.split(':')[1])
		test = 1
	if test != 0:
		if segsites != 0:
			if 'segsites' not in line and "positions" not in line and line!="\n" and line!='//' and line!='':
				cnt +=1
				if cnt%2==0:
					cnt = 0
					nInd += 1
				line = line.replace('0', 'A')
				line = line.replace('1', 'G')
#				line += 'A'*(L[locusID]-segsites)
#				line += 'A'*(L-segsites)
				if nInd <= nA:
					species = 'spA'
				else:
					species = 'spB'
				individual = nInd
				allele = cnt + 1
				res += '>locus_{0}|{1}|{2}|allele{3}\n{4}\n'.format(locusID, species, individual, allele, line)
print(res)
