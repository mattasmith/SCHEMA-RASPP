#! /usr/local/bin/python

import sys, os, math, string, random, copy, time
import schema

HUGE_NUMBER = 1.0e10

# pse4lacalignment, min 8, 3 crossovers
# 24.667	10.222	8 16 24 
# 25.111	14.444	11 20 259 
# 25.556	15.333	11 20 257 
# 25.667	16.000	11 20 256 
# 27.000	17.111	11 22 255 
# 28.889	18.222	11 19 38 
# 29.556	20.444	11 19 43 
# 31.000	21.778	19 242 254 
# 30.778	22.222	20 242 251 
# 31.667	23.333	22 241 251 
# 32.556	24.000	22 240 251 
# 32.556	55.926	41 137 147 
# 34.556	56.222	41 137 149 
# 34.667	57.654	41 137 158 
# 35.333	58.444	41 137 162 
# 36.222	59.580	41 137 170 
# 37.222	60.741	41 125 170 
# 38.667	61.012	42 125 171 
# 39.333	62.000	41 125 176 
# 40.667	63.543	43 125 186 
# 41.889	64.148	43 119 192 

# 27.00 11.56   [10, 18, 26]
# 28.22 17.11   [12, 22, 36]
# 32.78 18.22   [14, 26, 38]
# 33.44 22.00   [22, 36, 49]
# 34.11 22.44   [22, 37, 51]
# 34.56 56.27   [44, 133, 148]
# 34.89 56.57   [44, 126, 143]
# 34.67 57.26   [44, 138, 154]
# 35.33 57.93   [44, 138, 161]
# 35.33 58.32   [44, 138, 162]
# 36.67 58.81   [44, 133, 162]
# 36.22 59.23   [44, 138, 167]
# 36.22 59.68   [44, 138, 171]
# 37.22 60.91   [44, 126, 171]
# 39.00 61.23   [49, 120, 171]

#4.444	24.605	41 54 66 
#5.444	25.062	40 54 70 
#7.444	26.086	35 47 66 
#8.333	27.136	21 41 70 

# 4.44  24.31   [42, 55, 67]
# 6.00  25.01   [40, 55, 67]
# 6.22  25.33   [38, 48, 67]
# 6.67  25.77   [38, 55, 71]
# 8.22  27.01   [22, 42, 67]

class Path:
	def __init__(self, indices, length):
		self.indices = indices
		self.length = length

	def frags(self,parent):
		return schema.getFragments(self.crossovers(), parent)
	
	def crossovers(self):
		return [c+1 for c in self.indices]
		
	def frag_lengths(self,parent):
		return [l-k for (k,l) in self.frags(parent)]

def make_2d_array(dim1, dim2, init_value = 0.0):
	arr = []
	for i in range(dim1):
		arr.append([init_value]*dim2)
	return arr

def min_index(arr):
	# return minimum and index
	cur_ind = 0
	cur_min = arr[cur_ind]
	for i in range(len(arr)):
		if arr[i] < cur_min:
			cur_ind = i
			cur_min = arr[i]
	return (cur_min, cur_ind)
	
def calc_average_energies_from_contacts(contacts, parents):
	num_residues = len(parents[0])
	num_parents = len(parents)
	avg_energies = []
	for (i,j,ri,rj) in contacts:
		sum = 0.0
		for p in range(num_parents-1):
			parp = parents[p]
			for q in range(p+1,num_parents):
				parq = parents[q]
				# If contact can be broken by these parents
				if (parp[i], parp[j]) != (parq[i], parq[j]):
					# Note that summing the disruption will not work in the
					# case of an arbitrary energy function, where we must
					# add the difference of energies between the disrupted
					# and undisrupted states.
					# Add two because we're only iterating over half the matrix
					sum += 2.0
		avg = sum/(num_parents**2)
		avg_energies.append((i,j,avg))
	# Now construct a matrix
	avg_energy_matrix = make_2d_array(num_residues, num_residues, 0.0)
	for (i,j,avg) in avg_energies:
		avg_energy_matrix[i][j] = avg
		avg_energy_matrix[j][i] = avg
	return avg_energy_matrix

def order_contacts(contacts):
	new_contacts = []
	for (i,j,ri,rj) in contacts:
		if i>j:
			new_contacts.append((j,i,rj,ri))
		else:
			new_contacts.append((i,j,ri,rj))
	new_contacts.sort()
	return new_contacts

def make_4d_energies(contacts, parents):
	ordered_contacts = order_contacts(contacts)
	energies = []
	for (i,j,ri,rj) in ordered_contacts:
		for p in range(len(parents)):
			parp = parents[p]
			for q in range(len(parents)):
				if p != q:
					parq = parents[q]
					pair = (parp[i], parq[j])
					if pair not in [(r[i], r[j]) for r in parents]:
						energies.append((i,j,p,q))
	return energies

def calc_average_energies(energies, parents):
	num_residues = len(parents[0])
	num_parents = len(parents)
	avg_energies = []
	#energy_dict = dict([(x,1) for x in energies])
	for i in range(num_residues-1):
		for j in range(i+1, num_residues):
			ij = [(ri,rj,p,q) for (ri,rj,p,q) in energies if ri==i and rj==j]
			nij = len(ij)
			nij_same = len([(ri,rj) for (ri,rj,p,q) in ij if ri==i and rj==j and p==q])
			#print i, j, nij, nij_same
			diff = nij - nij_same
			avg = float(diff)/(num_parents**2)
			avg_energies.append((i,j,avg))
	return avg_energies
	
def avg_energy_list_to_matrix(avg_energies, num_residues, init_value=0.0):
	# Now construct a matrix
	avg_energy_matrix = make_2d_array(num_residues, num_residues, init_value)
	for (i,j,avg) in avg_energies:
		avg_energy_matrix[i][j] = avg
		avg_energy_matrix[j][i] = avg
	return avg_energy_matrix

def calc_average_library_energy(avg_energy_matrix, num_residues, xover1, xover2):
	# Compute average energy of a library with crossovers at 
	# positions specified by crossovers.
	# The fragments are defined to end with xover1 and xover2
	lib_avg_energy = 0.0
	for i in range(xover1, xover2):
		for j in range(xover2, num_residues):
			lib_avg_energy += avg_energy_matrix[i][j]
	return lib_avg_energy

def calc_arc_lengths(avg_energies, parents):
	# Calculate the lengths, in average energies, from crossovers starting
	# at residue r1 to residue r2.
	# The fragments here begin at r1 and end at r2-1
	num_residues = len(parents[0])
	num_parents = len(parents)
	avg_energy_matrix = avg_energy_list_to_matrix(avg_energies, num_residues, 0.0)
	arc_lengths = make_2d_array(num_residues, num_residues, HUGE_NUMBER)
	for r1 in range(num_residues-1):
		for r2 in range(r1+1, num_residues):
			avg_energy = calc_average_library_energy(avg_energy_matrix, num_residues, r1, r2)
			arc_lengths[r1][r2] = avg_energy
	return arc_lengths

def add_list(x, y):
	return [a+b for (a,b) in zip(x,y)]

def collapse_parents(parents):
	"""Returns parental alignment with identical sites removed.

	RASPP treats fragment lengths in terms of numbers of non-identical
	mutations."""

	# The list of sites having perfect identity across parents
	identity_list = []
	num_parents = len(parents)
	# This procedure only makes sense with more than 1 parent
	if len(parents) < 2:
		return parents, []
	
	for i in range(len(parents[0])):
		sites = [p[i] for p in parents]
		c = sites.count(sites[0])
		if c == num_parents: # All parents have same residue
			# So store the position of the identity
			identity_list.append(i)
	# If there are no identities, bail.
	if len(identity_list) == 0:
		return parents, []
	
	# Now, if there are identities, create fragments (i,j)
	# that only include non-identities.
	frags = []
	cur = 0
	for i in identity_list:
		frags.append((cur,i))
		cur = i+1
	if cur < len(parents[0]):
		frags.append((cur,len(parents[0])))
	# Assemble the new parents by assembling the fragments
	new_parents = ['']*num_parents
	for i in range(num_parents):
		for (m,n) in frags:
			new_parents[i] += parents[i][m:n]
	return new_parents, identity_list

def collapse_matrix(matrix, identical_sites):
	dim = len(matrix)
	nsites = len(identical_sites)
	new_matrix = []
	for i in range(dim-nsites):
		new_matrix.append([0.0]*(dim-nsites))
	
	k = 0
	l = 0
	for i in range(len(matrix)):
		if not i in identical_sites:
			l = 0
			for j in range(len(matrix[i])):
				if not j in identical_sites:
					new_matrix[k][l] = matrix[i][j]
					l += 1
			k += 1
	return new_matrix

def translate_collapsed_indices(collapsed_crossovers, collapsed_sites):
	"""Converts crossover indices generated on parents whose identical sites have
	been removed back into indices relative to the full-length parents."""
	crossovers = collapsed_crossovers[:]
	crossovers.sort()
	collapsed_sites.sort()

	# For each collapsed site,
	for site in collapsed_sites:
		for i in range(len(crossovers)):
			if crossovers[i] > site:
				crossovers[i]+= 1
	return crossovers
		
def RASPP_SCHEMA(contacts, parents, num_crossovers, min_fragment_diversity):
	schema_contacts = schema.getSCHEMAContacts(contacts, parents)
	(collapsed_parents, identity_list) = collapse_parents(parents)
	energies = raspp.make_4d_energies(schema_contacts, parents)
	avg_energies = raspp.calc_average_energies(energies, parents)
	results = RASPP(avg_energies, parents, num_crossovers, min_fragment_diversity)
	
	for i in range(len(results)):
		(avg_E, collapsed_crossovers, l_min, l_max) = results[i]
		crossovers = translate_collapsed_indices(collapsed_crossovers, identity_list)
		results[i] = (avg_E, crossovers, l_min, l_max)
	return results

	
def RASPP(avg_energies, parents, num_crossovers, min_fragment_diversity):
	"""Find libraries with the lowest energy given constraints on fragment diversity."""
	
	# Collapse the parents to remove identical sites.  Necessary because RASPP's
	# definition of fragment "length" -- the number of changed residues in a fragment --
	# depends on the absence of universally conserved sites.
	(collapsed_parents, identical_sites) = collapse_parents(parents)
	if (num_crossovers+1)*min_fragment_diversity > len(collapsed_parents[0]):
		err_string = "%d crossovers with minimum fragment length of %d " % \
			  (num_crossovers, min_fragment_diversity) + \
			  "is impossible given parent non-identical sequence length of %d." % \
			  (len(collapsed_parents[0]),)
		raise ValueError, err_string
	

	results = []
	# Compute the arc lengths.	
	tstart = time.clock()
	arc_lengths = calc_arc_lengths(avg_energies, parents)
	ttot = time.clock()-tstart
	#print "# Arc lengths calculated in %1.2f sec" % ttot
	num_residues = len(collapsed_parents[0])
	
	# Now collapse the arc lengths.  It is more efficient to 
	# do this before calculating the arc lengths, but for user-friendliness
	# we do it now.
	arc_lengths = collapse_matrix(arc_lengths, identical_sites)
	
	# Iterate over all possible constraints on the fragment length
	min_l_min = min_fragment_diversity
	max_l_min = int(math.floor(num_residues/(num_crossovers+1.0)))+1
	for l_min in range(min_l_min, max_l_min+1):
		#print "l_min:", l_min, min_l_min, max_l_min
		min_l_max =int(math.ceil(num_residues/(num_crossovers+1.0)))
		max_l_max = num_residues-num_crossovers*l_min+1
		for l_max in range(min_l_max, max_l_max+1):
			#print "  l_max:", l_max, min_l_max, max_l_max+1
			# Find set of crossovers which minimize the average energy consistent with fragment-length constraints
			res = get_shortest_path(arc_lengths, collapsed_parents, num_crossovers, l_min, l_max)
			if res:
				results.append(res)
	# Convert results back into full-length parent indices
	for i in range(len(results)):
		(avg_E, collapsed_crossovers, l_min, l_max) = results[i]
		crossovers = translate_collapsed_indices(collapsed_crossovers, identical_sites)
		results[i] = (avg_E, crossovers, l_min, l_max)
	return results
	
def get_shortest_path(arc_lengths, parents, num_crossovers, l_min, l_max):
	"""Finds the set of crossovers which minimizes the average library energy.
	
	Exactly num_crossovers crossovers will be found, and they will satisfy the
	constraints that no resulting fragment be shorter than l_min or longer
	than l_max."""
	debug = False
	num_residues = len(parents[0])

	paths = []
	# First column:
	for i in range(num_residues):
		# The shortest path to each node in the first column is
		# just the arc length to that column.
		#
		# i represents the 0-based index of the beginning of
		# the new fragment.  E.g., i=5 means that the first
		# fragment begins at residue 0 and ends at residue 4
		# (length = 5), and the next fragment begins at res. 5.
		length = i
		if l_min <= length <= l_max:
			p = Path([i], arc_lengths[0][i])
			paths.append(p)

	# The shortest path to column k given the SPs to column k-1
	# is the one having min_i {SP^(k-1)_i + A(i,j)}.  (This path
	# ends at node j of column k.)

	# Goal: given SP to all nodes in col k-1, compute shortest path
	# to all nodes j in col k.
	# The paths 
	for k in range(1, num_crossovers):
		if debug:
			print "k=%d" % k
			for p in paths:
				print "%s\t%1.2f" % (p.frag_lengths(parents[0]), p.length)
		# The only allowable range for node j is where previous
		# fragments satisfy the length constraints.  There have
		# been k of them so far, so j >= k*l_min.
		min_j = k*l_min
		# There will be num_crossovers+1-k fragments after, and 
		# those fragments, if minimum length, cannot go past
		# the end of the protein, so j <= num_residues-(num_crossovers-k+1)*l_min.
		# For example, if there are 50 residues, this is the third
		# crossover of four, and l_min = 10, then j must be no more than 30.
		max_j = num_residues-(num_crossovers-k)*l_min+1
		#print "minj, maxj", min_j, max_j
		# At the beginning of this loop, paths contains the set of shortest
		# paths that are allowable and end nodes given by the final entry
		# in the .indices element of each Path object in shortest_paths.
		# At the end of this loop, shortest_paths will contain the set of shortest paths
		# that are allowable and end at nodes j.
		shortest_paths = []
		for j in range(min_j, max_j):
			best_path = None  #Path([], HUGE_NUMBER)		
			for p in paths:
				i = p.indices[-1]
				# Track the shortest path 
				frag_length = j-i
				if l_min <= frag_length <= l_max:
					# This is an allowable path
					# If it's shorter than the shortest path so far, save it
					if not best_path or p.length + arc_lengths[i][j] < best_path.length:
						best_path = Path(p.indices+[j], p.length + arc_lengths[i][j])
			# If at least one allowable path was found, we've got a shortest path to j.
			if best_path: #.length < HUGE_NUMBER:
				shortest_paths.append(best_path)
		# Now start all over again with this new set of shortest paths to column k.
		paths = shortest_paths

	if debug:
		print "final"
		for p in paths:
			print "%s\t%1.2f" % (p.frag_lengths(parents[0]), p.length)
		#for i in range(len(paths)):
		#	print "      %d\t%s\t%1.2f" % (i, paths[i].frag_lengths(parents[0]), paths[i].length)

	# Now filter out paths with illegal last-fragment lengths
	allowable_paths = []
	for p in paths:
		final_frag_length = num_residues - p.indices[-1]
		if l_min <= final_frag_length <= l_max:
			allowable_paths.append(p)
	paths = allowable_paths

	# Find the minimum path length
	lengths = [p.length for p in paths]
	if len(lengths) > 0:
		(avg_energy, i) = min_index(lengths)
		#frag_lengths = paths[i].frag_lengths(parents[0])
		#if min(frag_lengths) < l_min or max(frag_lengths)>l_max:
		#	print "Illegal lengths: min %d, max %d, actual %s" % (l_min, l_max, frag_lengths)
		return (avg_energy, paths[i].crossovers(), l_min, l_max)
	else:
		return None

def curve(results, parents, bin_width, max_samples=1e10):
	"""Compute a curve of average energy and average mutation, with the latter binned
	by bin_width.
	"""
	if len(results) < 1:
		# Nothing to do!
		return
	(e, xovers, lmin, lmax) = results[0]
	num_crossovers = len(xovers)
	#print "# No. of RASPP results:", len(results)
	# Because the RASPP curve does not involve the length min/max values, we can collapse
	# the set of RASPP results to only those with unique crossovers.  This can
	# greatly improve performance, since crossover patterns are often duplicated
	# across min/max values.
	unique_results = set([(avg_energy, tuple(crossovers)) for (avg_energy, crossovers, l_min, l_max) in results])
	#print "# No. of unique RASPP results:", len(unique_results)
	# Now compute the average mutation levels for these unique libraries.
	avg_E_ms = []
	for (avg_energy, crossovers) in unique_results:
		crossovers = list(crossovers)
		fragments = schema.getFragments(crossovers, parents[0])
		avg_m = schema.averageMutationSampled(fragments, parents, max_samples)
		avg_E_ms.append((avg_energy, avg_m, crossovers))
		
	ms = [m for (E, m, crossovers) in avg_E_ms]
	(min_m, max_m) = (min(ms), max(ms))
	num_bins = int((max_m-min_m)/bin_width) + 1

	# Assemble the RASPP curve.  If num_samples exceeds the library size,
	# then this curve is approximate.
	approx_curve = []
	for i in range(num_bins):
		approx_curve.append(None)
	
	for (E, m, crossovers) in avg_E_ms:
		bin_number = int((m - min_m)/bin_width)
		# If there's an existing value in this bin, check it
		if approx_curve[bin_number]:
			(E_old, m_old, crossovers_old) = approx_curve[bin_number]
			# If lower E in this bin, substitute it
			if E < E_old:
				approx_curve[bin_number] = (E, m, crossovers)
		else:  # Otherwise just add it
			approx_curve[bin_number] = (E, m, crossovers)

	# It may be that the approximate curve is exact.  If so, just return it.
	library_size = len(parents)**(num_crossovers+1)
	approximate = (library_size > max_samples)
	if not approximate:
		return [r for r in approx_curve if r]

	# If the curve IS approximate, we'll do a final pass so that
	# the bin values are correct.  Some libraries may still be
	# incorrectly binned because we've 
	# Compute the exact mutation numbers for the lowest-E libraries
	final_curve = []
	for r in approx_curve:
		if r:
			(E, approx_m, crossovers) = r
			fragments = schema.getFragments(crossovers, parents[0])
			true_avg_m = schema.averageMutation(fragments, parents)
			final_curve.append((E, true_avg_m, crossovers))
			#print "%1.2f\t%1.2f" % (true_avg_m, approx_m)
		
	return final_curve

