# load 'parse', a useful function from the library Biopython to parse fasta files, i.e, to read them
from Bio.SeqIO import parse

# load the library 'sys', to get arguments from the command line
import sys

def PCA(alignment):
	# function that prepares an input file for a PCA in R
	# alignment is a dictionary where entries are species: alignment[species]
	# alignment['spA']; alignment['spB']; etc ...

	# AA = 0
	# AG = GA = 1
	# GG = 2
	nSNPs = len(alignment['spA'][0])
	
	output = "species\tindividual"
	
	for SNP_i in range(nSNPs):
		output += '\tSNP_{snp}'.format(snp=SNP_i)
	output += '\n'
	
	# loop over individuals from spA
	cnt = 0
	for individual_i in range(0, len(alignment['spA']), 2):
		cnt += 1
		output += "spA\t{ind}".format(ind=cnt)
		# loop over positions
		for SNP_i in range(nSNPs):
			genotype = alignment['spA'][individual_i][SNP_i] + alignment['spA'][individual_i + 1][SNP_i]
			if genotype == "AA":
				recoded_genotype = 0
			else:
				if genotype == "GG":
					recoded_genotype = 2
				else:
					recoded_genotype = 1
			output += '\t{recoded_genotype}'.format(recoded_genotype=recoded_genotype)
		output += '\n'

	# loop over individuals from spB
	cnt = 0
	for individual_i in range(0, len(alignment['spB']), 2):
		cnt += 1
		output += "spB\t{ind}".format(ind=cnt)
		# loop over positions
		for SNP_i in range(nSNPs):
			genotype = alignment['spB'][individual_i][SNP_i] + alignment['spB'][individual_i + 1][SNP_i]
			if genotype == "AA":
				recoded_genotype = 0
			else:
				if genotype == "GG":
					recoded_genotype = 2
				else:
					recoded_genotype = 1
			output += '\t{recoded_genotype}'.format(recoded_genotype=recoded_genotype)
		output += '\n'
	return(output)


def pi(sequences):
	# function that computes pi for a list of sequences
	# number of individuals
	nIndividuals = len(sequences)
	
	# number of SNPs
	nSNPs = len(sequences[0])
	
	# computing pi
	pi = 0
	## loop over SNPs
	for position_i in range(nSNPs):
		list_of_alleles = []
		for individual_i in range(nIndividuals):
			list_of_alleles.append(sequences[individual_i][position_i])
		segregating_alleles = list(set(list_of_alleles))
		
		if len(segregating_alleles) > 1:
			## loop over alleles
			pi_tmp = 1
			for allele_i in segregating_alleles:
				pi_tmp -= (list_of_alleles.count(allele_i)/len(list_of_alleles))**2
			pi += pi_tmp
	return(pi)
			
			

# name of the fasta file
fileName = sys.argv[1]

# parse the fasta file
## the dictionary "alignment" will store the read sequences
alignment = {}

## parse : open and read the fasta file
infile = parse(fileName, 'fasta')
for seq in infile:
	sequence_name = seq.id
	species_name = sequence_name.split("|")[1]
	if species_name not in alignment:
		alignment[species_name] = []
	alignment[species_name].append(seq.seq)

# now, all sequences found in the fasta file are recorded in the dictionary "alignment"
pi_tot = pi(alignment['spA']+alignment['spB'])

pi_A = pi(alignment['spA'])
pi_B = pi(alignment['spB'])
pi_S = (pi_A + pi_B)/2.0

Fst = (pi_tot-pi_S)/pi_tot

print("piA\tpiB\tpiS\tpiT\tFst\n{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}".format(pi_A, pi_B, pi_S, pi_tot, Fst))


## PCA
PCA_input = PCA(alignment)
outfile = open("PCA.txt", "w")
outfile.write(PCA_input)
outfile.close()

