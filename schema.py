#! /usr/local/bin/python
"""Module for SCHEMA tools.

This module provides core functions used by all the SCHEMA tools.
 
    ******************************************************************
    Copyright (C) 2005  Allan Drummond, California Institute of Technology

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    *******************************************************************


SCHEMA was developed in the laboratory of Frances H. Arnold at the California Institute of Technology.

References:

Voigt, C. et al., "Protein building blocks preserved by recombination," Nature Structural Biology 9(7):553-558 (2002).
Meyer, M. et al., "Library analysis of SCHEMA-guided recombination," Protein Science 12:1686-1693 (2003).
Otey, C. et al., "Functional evolution and structural conservation in chimeric cytochromes P450: Calibrating a structure-guided approach," Chemistry & Biology 11:1-20 (2004)
Silberg, J. et al., "SCHEMA-guided protein recombination," Methods in Enzymology 388:35-42 (2004).
Endelman, J. et al., "Site-directed protein recombination as a shortest-path problem," Protein Engineering, Design & Selection 17(7):589-594 (2005).
"""

import sys, string, random
import pdb

DIGITS_LETTERS = string.digits + string.letters

def alignPDBResidues(residues, aligned_parent_protein, aligned_pdb_protein, library_parent_protein, chain_identifiers):
	"""Align the PDB residues to the aligned parent protein sequence, following the parent/PDB alignment."""
	residues = [res for res in residues if res.chain in chain_identifiers]
	new_residues = []
	j = 0
	# Prune/add to the residue list to make it match the library_parent_protein
	# Assumptions:
	#   1) The PDB sequence contains only amino acids, so that there are no explicit gaps.
	#   2) The aligned_parent_protein and library_parent_protein differ only in the number of gaps.

	# The PDB residues corresponding to sites in library_parent_protein will be accumulated
	# in the new_residues list.
	new_residues = []

	# The process of assembling is somewhat involved.  The indices i, j, and k
	# indicate the current position in the library parent, parent/PDB alignment and PDB
	# residue sequences, respectively.  There are four cases:
	# a) Library parent has a gap at the current position, and:
	#      1) PDB-aligned parent also has a gap, but the parent-aligned PDB does not
	#      2) Both have gaps
	# b) Library parent has an amino acid at the current position, and:
	#      3) PDB-aligned parent does not match the library parent
	#      4) PDB-aligned parent matches
	#           4.1) Parent-aligned PDB has a gap
	#           4.2) Parent-aligned PDB has a residue
	
	# The j index indicates positions in the aligned PDB protein, aligned_pdb_protein
	# and aligned_parent_protein
	j = 0
	# The k index indicates positions in the PDB residue sequence, residues
	k = 0
	for i in range(len(library_parent_protein)):
		li = library_parent_protein[i]
		# If the parent has a gap, 
		if li == '-':
			parj = aligned_parent_protein[j]
			pdbj = aligned_pdb_protein[j]
			if parj == '-' and pdbj != '-': # Parents have gaps, but PDB doesn't
				new_residues.append(residues[k])
				j += 1
				k += 1
			elif parj == pdbj == li: # All have gaps
				new_residues.append(pdb.Residue(None))
				j += 1
		else: # parent does not contain a gap at i
			# While the aligned parent doesn't match the library parent,
			# advance the indices.  This brings the sequences back into
			# register.
			pdbres = pdb.three_to_one_map[residues[k].residue]
			#print li, aligned_parent_protein[j], aligned_pdb_protein[j], pdbres
			while aligned_parent_protein[j] != li:
				pdbres = pdb.three_to_one_map[residues[k].residue]
				if aligned_pdb_protein[j] != '-':
					k += 1
				j += 1
				#print "*", li, aligned_parent_protein[j], aligned_pdb_protein[j], pdbres
			#print "end"
			# Now, if the aligned PDB has a gap, insert a gap
			# into the residues list, otherwise insert the residue.
			if aligned_pdb_protein[j] == '-':
				new_residues.append(pdb.Residue(None))
			else:
				pdbres = pdb.three_to_one_map[residues[k].residue]
				#print i, j, k
				#print li, aligned_parent_protein[j], aligned_pdb_protein[j], pdbres
				#print library_parent_protein[:i+1]
				#print aligned_parent_protein[:j+1]
				#print aligned_pdb_protein[:j+1]
				#print pdb.sequence(new_residues)
				#print pdb.sequence(residues)[:k+1]
				# Check to make sure we're adding the expected residue
				if aligned_pdb_protein[j] != pdbres:
					err_string = "Expected residue %s at aligned position %d based on alignment, but PDB had %s at corresponding residue %d.  Aborting..."
					raise ValueError, err_string % (aligned_pdb_protein[j], j+1, pdbres, k+1)
				new_residues.append(residues[k])
				k += 1
			j += 1
		if j >= len(aligned_pdb_protein) or k >= len(residues):
			break
	return new_residues

def getSCHEMAContactsWithCrossovers(contacts, parents, crossovers):
	"""Get contacts with correction for parental sequence identity
	and for fragments."""
	fragments = getFragments(crossovers, parents[0])
	filtered_contacts = []
	# First eliminate contacts between residues in the same fragment.  These
	# can never be broken by recombination.
	for (i, j, ri, rj) in contacts:
		same_fragment = False
		for (k,l) in fragments:
			if i>=k and i<l and j>=k and j<l:
				same_fragment = True
		# Now eliminate contacts where either residue is absolutely conserved.
		# These can never be broken by recombination.
		if not same_fragment:
			pi = [x[i] for x in parents]
			pj = [x[j] for x in parents]
			if not (pi.count(pi[0])==len(pi) or pj.count(pj[0])==len(pj)): # if neither i or j absolutely conserved
				# This is a bona-fide breakable contact.
				filtered_contacts.append((i,j,ri,rj))
	return filtered_contacts

def getSCHEMAContacts(contacts, parents):
	"""Get contacts with correction for parental sequence identity."""
	filtered_contacts = []
	# First eliminate contacts between residues in the same fragment.  These
	# can never be broken by recombination.
	for (i, j, ri, rj) in contacts:
		# Eliminate contacts where either residue is absolutely conserved.
		# These can never be broken by recombination.
		pi = [x[i] for x in parents]
		pj = [x[j] for x in parents]
		# if neither i or j absolutely conserved
		if not (pi.count(pi[0])==len(pi) or pj.count(pj[0])==len(pj)): 
			# This is a bona-fide breakable contact.
			filtered_contacts.append((i,j,ri,rj))
	return filtered_contacts

def getPDBContacts(residues, contact_distance):
	# Get set of residues within contact_distance angstroms in a PDB
	contacts = []
	# Loop through all residue pairs
	for i in range(0,len(residues)-1):
		for j in range(i+1, len(residues)):
			resi = residues[i]
			resj = residues[j]
			# If both residues are present (gaps == None)
			if resi and resj:
				# If residues are within contact_distance angstroms of each other, add as a contact
				contact = resi.isContact(resj, contact_distance)
				if contact:
					contacts.append((i, j, resi, resj))
	return contacts

def writeContactFile(contacts, outfile):
	outfile.write("# Fields are number, contacting residues i & j (in library coordinates), and residues i & j (in PDB coordinates)\n")
	outfile.write("# n	i	j	pdbi	pdbj\n")
	for k in range(len(contacts)):
		(i, j, resi, resj) = contacts[k]
		outfile.write("%d\t%d\t%d\t%s\t%s\n" % (k, i, j, resi.res_seq, resj.res_seq))

def readContactFile(f):
	contacts = []
	for line in f.readlines():
		if line[0] == '#':  # comment line
			continue
		flds = line.strip().split()
		contacts.append((int(flds[1]), int(flds[2]), flds[3], flds[4]))
	return contacts

def checkChimera(chimera_blocks, fragments, parents):
	"""Checks a chimera to see if it's valid."""
	# First see if there are the proper number of blocks
	if len(chimera_blocks) != len(fragments):
		return False
	# Then ensure that, within each block, the parent index
	# is no smaller than 1 or larger than the number of parents.
	for i in range(len(chimera_blocks)):
		if not (0 < int(chimera_blocks[i]) <= len(parents)):
			return False
	return True

def base(number, radix):
	"""Returns a string representation of the number in the base indicated by radix."""
	# Inverse function to int(str,radix) and long(str,radix)
	if not 2 <= radix <= 36:
		raise ValueError, "radix must be in 2..36"
	result = []
	addon = result.append
	if number < 0:
		number = -number
		addon('-')
	elif number == 0:
		addon('0')

	_divmod, _abc = divmod, DIGITS_LETTERS
	while number:
		number, rdigit = _divmod(number, radix)
		addon(_abc[rdigit])
	result.reverse()
	return ''.join(result)

def indexToFragment(index, fragments):
	"""Given an index into a full protein sequence, returns the zero-based index of the 
	fragment containing that position."""
	for frag_index in range(len(fragments)):
		(i,j) = fragments[frag_index]
		if index >= i and index < j:
			return frag_index
	return None

def getChimeraDisruption(chimera_blocks, contacts, fragments, parents):
	"""Takes a chimera block pattern, such as '11213312', and computes the SCHEMA
	disruption, the number of contacts broken by recombination."""
	parent_indices = [int(c)-1 for c in chimera_blocks]
	num_disruptions = 0
	for (i,j,ri,rj) in contacts:
		frag_i = indexToFragment(i, fragments)
		frag_j = indexToFragment(j, fragments)
		if parent_indices[frag_i] == parent_indices[frag_j]:
			# No disruption possible if both fragments come from the same parent.
			continue
		pair = (parents[parent_indices[frag_i]][i],	parents[parent_indices[frag_j]][j])
		# If pair doesn't exist in any parent, it's counted as disruptive
		if pair not in [(p[i], p[j]) for p in parents]:
			#print chimera_blocks, i+1, j+1
			num_disruptions += 1
	return num_disruptions

def getChimeraSequence(chimera_blocks, fragments, parents):
	"""Converts a chimera block pattern, such as '11213312', into a protein sequence
	by assembling fragments from the parents.  This sequence may then be used
	to compute mutational distances and so on."""
	parent_indices = [int(c)-1 for c in chimera_blocks]
	chimera = ''
	for i in range(len(parent_indices)):
		which_parent = parent_indices[i]
		(begin, end) = fragments[i]
		chimera += parents[which_parent][begin:end]
	return chimera

def getChimeraDistance(chimera, parent):
	"""Returns the number of mutations necessary to turn the chimera into the parent."""
	return len([(caa, paa) for (caa,paa) in zip(chimera, parent) if (caa != paa)])

def getChimeraShortestDistance(chimera_blocks, fragments, parents):
	"""Returns the minimum number of mutations necessary to turn the chimera into
	at least one of the parent sequences."""
	chimera = getChimeraSequence(chimera_blocks, fragments, parents)
	ms = [getChimeraDistance(chimera, p) for p in parents]
	return min(ms)

def readMultipleSequenceAlignmentFile(f):
	"""Reads a multiple sequence alignment file in ALN format."""
	parent_dict = {}
	keys = []
	for line in f.readlines():
		if line[0] == '#':  # comment line
			continue
		flds = line.strip().split()
		# Skip lines with too little data to be sequence info
		if len(flds) < 2:
			continue
		# Skip fields with incomprehensible data, such as alignment-score lines.
		if flds[0][0] not in DIGITS_LETTERS:
			continue

		key = flds[0]
		seq = flds[1]
		# Extract the parent sequence
		if parent_dict.has_key(key):
			parent_dict[key] += seq
		else:
			keys.append(key)
			parent_dict[key] = seq
	# Join the sequences and replace dots with dashes to indicate gaps.
	parent_list = [(key,''.join(parent_dict[key]).replace('.','-')) for key in keys]
	return parent_list

def readCrossoverFile(f):
	"""Read a file of crossover positions."""
	crossovers = [0]
	for line in f.readlines():
		line = line.strip()
		if line[0] == '#' or line == '':  # comment line or empty line
			continue
		try:
			crossovers = [int(x) for x in line.split()]
		except ValueError:
			print "Non-numeric data found in crossover file:\n", line.split()
	crossovers.sort()
	return crossovers

def mutationMatrix(fragments, parents):
	mut_dict = {}
	for f in range(len(fragments)):
		(i,j) = fragments[f]
		for p in range(len(parents)-1):
			for q in range(p+1,len(parents)):
				frag_p = parents[p][i:j]
				frag_q = parents[q][i:j]
				m = len([(a,b) for (a,b) in zip(frag_p,frag_q) if a != b])
				mut_dict[(f,p,q)] = m
				mut_dict[(f,q,p)] = m
	return mut_dict

def getChimeraShortestDistanceLookup(chimera_blocks, fragments, parents, mut_dict):
	ms = [0]*len(parents)
	for p in range(len(parents)):
		for f in range(len(fragments)):
			q = int(chimera_blocks[f])-1
			if q != p:
				key = (f,p,q)
				ms[p] += mut_dict[key]
	return min(ms)

def averageMutation(fragments, parents):
	avg_m = 0.0
	num_chimeras = 0
	p = len(parents)
	n = len(fragments)
	# Create a mutation matrix: how many mutations
	mut_dict = mutationMatrix(fragments, parents)			
	for i in range(p**n):
		# The next two lines turn i into a chimera block pattern 
		# (e.g., 0 -> '11111111', 1 -> '11111112', 2 -> '11111113'...)
		n2c = base(i,p)
		chimera_blocks = ''.join(['1']*(n-len(n2c))+['%d'%(int(x)+1,) for x in n2c])
		m = getChimeraShortestDistanceLookup(chimera_blocks, fragments, parents, mut_dict)
		avg_m += m
		num_chimeras += 1
	return avg_m/num_chimeras
	
def averageMutationSampled(fragments, parents, num_samples):
	# Estimate the average mutation level by taking num_samples samples.
	avg_m = 0.0
	num_chimeras = 0
	p = len(parents)
	n = len(fragments)
	# Create a mutation matrix: how many mutations
	mut_dict = mutationMatrix(fragments, parents)
	library_size = p**n
	if num_samples >= library_size: # Might as well be exact!
		return averageMutation(fragments, parents)
	
	for i in range(num_samples):
		rand_num = random.randint(0,library_size)
		# The next two lines turn i into a chimera block pattern 
		# (e.g., 0 -> '11111111', 1 -> '11111112', 2 -> '11111113'...)
		n2c = base(rand_num,p)
		chimera_blocks = ''.join(['1']*(n-len(n2c))+['%d'%(int(x)+1,) for x in n2c])
		m = getChimeraShortestDistanceLookup(chimera_blocks, fragments, parents, mut_dict)
		avg_m += m
		num_chimeras += 1
	return avg_m/num_chimeras
	
def averageEnergy(contacts, fragments, parents):
	avg_E = 0.0
	num_chimeras = 0
	p = len(parents)
	n = len(fragments)
	for i in range(len(parents)**len(fragments)):
		# The next two lines turn i into a chimera block pattern 
		# (e.g., 0 -> '11111111', 1 -> '11111112', 2 -> '11111113'...)
		n2c = base(i,p)
		chimera_blocks = ''.join(['1']*(n-len(n2c))+['%d'%(int(x)+1,) for x in n2c])
		E = getChimeraDisruption(chimera_blocks, contacts, fragments, parents)
		avg_E += E
		num_chimeras += 1
	return avg_E/num_chimeras

def getCrossoversFromFragments(fragments):
	"""Turns fragments, which are pairs of 0-based indices, into
	1-based crossover indices."""
	return [x+1 for (x,y) in fragments[1:]]

def getFragments(crossovers, parent):
	"""Turns crossover points, which are 1-based for readability, into (begin,end) 
	pairs of 0-based indices into the parent sequence."""
	xovers = crossovers[:]
	if 1 not in xovers:
		xovers = [1] + xovers
	xovers.append(len(parent)+1)
	fragments = [(xovers[i]-1, xovers[i+1]-1) for i in range(len(xovers)-1)]
	return fragments
	
def splitSequence(seq_length, min_length):
	"""Splits a sequence at a point randomly chosen from those which generate fragments
	no smaller than min_length."""
	return random.randint(min_length, seq_length-min_length)

def generateRandomCrossovers(seq_length, num_xovers, min_length):
	fragments = generateRandomFragments(seq_length, num_xovers+1, min_length)
	xovers = getCrossoversFromFragments(fragments)
	return xovers

def generateRandomFragments(seq_length, num_fragments, min_length):
	"""Generates a random partition that splits a sequence into num_fragments
	pieces, each of which is no shorter than min_length."""
	# We will assume that the inputs have been properly checked
	assert(num_fragments*min_length <= seq_length)

	# There's a trivial but very hard case which is that the only suitable
	# crossover pattern is one which divvies up the protein into fragments
	# of equal length.  Just return that pattern if it arises.
	if num_fragments*min_length == seq_length:
		fragments = []
		i = 0
		for f in range(num_fragments):
			fragments.append((i,i+min_length))
			i += min_length
		return fragments
	
	done = False
	while not done:
		fragments = [(0,seq_length)]
		tries = 0
		# It's possible that no suitable set of fragments will be found on any
		# given attempt, so we limit the number of tries at generating each
		# crossover point.
		max_tries = 4*num_fragments
		# Challenge: choose num_fragments-1 points randomly along the sequence
		# such that none of them lie within min_length of each other or of
		# the ends.
		crossovers = []
		while len(fragments) < num_fragments and not tries > max_tries:
			frag = random.choice(fragments)
			(start, end) = frag
			if end-start > 2*min_length:
				# There's enough room for a crossover
				crossover = start + splitSequence(end-start, min_length)
				# Now replace the whole piece with its two fragments
				fragments.remove(frag)
				fragments.append((start,crossover))
				fragments.append((crossover,end))
				if len(fragments) == num_fragments:
					done = True
			else:
				# Failed to find a good crossover pattern this time.  Try again.
				tries += 1

	fragments.sort()
	return fragments

def mean(x):
	"""Computes the average of a list of numbers."""
	if len(x) == 0:
		raise ValueError, "Cannot take average of zero-length list"
	return sum(x)/float(len(x))
