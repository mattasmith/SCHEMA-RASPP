#! /usr/local/bin/python
"""Script for randomly enumerating crossover points and evaluating average SCHEMA energy and mutation of the resulting libraries.

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

import sys, os, string, math, random, time
import pdb, schema, raspp

ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE = 'msa'
ARG_CONTACT_FILE = 'con'
ARG_OUTPUT_FILE = 'o'
ARG_MIN_FRAGMENT_SIZE = 'min'
ARG_NUM_LIBRARIES = 'libs'
ARG_MAX_CHIMERAS_PER_LIBRARY = 'chims'
ARG_NUM_CROSSOVERS = 'xo'
ARG_RANDOM_SEED = 'seed'
ARG_COMPARE = 'compare'  # Unused
ARG_HELP = 'help'

def parse_arguments(args):
	# Turn linear arguments into a dictionary of (option, [values,...]) pairs
	arg_dict = {}
	key = None
	for arg in args[1:]:
		if arg[0] == '-':
			key = arg[1:]
		else:
			if arg_dict.has_key(key):
				if arg_dict[key] is list:
					arg_dict[key] = arg_dict[key]+[arg]
				else:
					arg_dict[key] = [arg_dict[key],arg]
			else:
				arg_dict[key] = arg
	return arg_dict

def print_usage(args):
	print 'Usage: python', args[0].split(os.path.sep)[-1], " [options]"
	print "Options:\n", \
		'\t-%s <alignment file>\n' % ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE, \
		"\t-%s <contact file>\n" % ARG_CONTACT_FILE, \
		'\t-%s <# crossovers>\n' % ARG_NUM_CROSSOVERS, \
		'\t[-%s <# libraries to generate>]\n' % ARG_NUM_LIBRARIES, \
		'\t[-%s <random number seed>]\n' % ARG_RANDOM_SEED, \
		'\t[-%s <max. chimeras generated per library>]\n' % ARG_MAX_CHIMERAS_PER_LIBRARY, \
		'\t[-%s <min. fragment length>]\n' % ARG_MIN_FRAGMENT_SIZE, \
		'\t[-%s <output file>]' % ARG_OUTPUT_FILE


def confirm_arguments(arg_dict):
	# Are arguments okay?
	res = True
	arg_keys = arg_dict.keys()
	try:
		if len(arg_keys) == 0:
			res = False
			return
			
		if not ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE in arg_keys:
			print "  You must provide a library file (-%s <alignment file>)" % ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE
			res = False
		elif not os.path.isfile(arg_dict[ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE]):
			print "  Can't find library file %s" % arg_dict[ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE]
			res = False
			
		if not ARG_CONTACT_FILE in arg_keys:
			print "  You must provide a contact file (-%s <contact file>)" % ARG_CONTACT_FILE
			res = False
		elif not os.path.isfile(arg_dict[ARG_CONTACT_FILE]):
			print "  Can't find contact file %s" % arg_dict[ARG_CONTACT_FILE]
			res = False
			
		if not ARG_NUM_CROSSOVERS in arg_keys:
			print "  You must specify the number of crossovers (-%s <number of crossovers>)" % ARG_NUM_CROSSOVERS
			res = False
	except Exception, e:
		#print e
		res = False
	return res

def main(args):
	arg_dict = parse_arguments(args)
	if not confirm_arguments(arg_dict):
		if args[0].split(os.path.sep)[-1] == "schemarandom.py":
			print_usage(args)
		return

	# Flags and values
	print_E = False
	print_m = False
	
	# Inputs:
	#   The alignment/fragment file name.
	msa_file = arg_dict[ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE]

	# Read the alignment file to create a list of parents.
	# The parents will appear in the list in the order in which they appear in the file.
	parent_list = schema.readMultipleSequenceAlignmentFile(file(msa_file, 'r'))
	parents = [p for (k,p) in parent_list]
	
	# Get the contacts
	pdb_contacts = schema.readContactFile(file(arg_dict[ARG_CONTACT_FILE], 'r'))
	
	# Establish connection to output, either file or, if no output file is 
	# specified, to standard output.
	if arg_dict.has_key(ARG_OUTPUT_FILE):
		output_file = file(arg_dict[ARG_OUTPUT_FILE], 'w')
	else:
		output_file = sys.stdout

	# Get the number of libraries to evaluate.
	if arg_dict.has_key(ARG_NUM_LIBRARIES):
		num_libraries = int(arg_dict[ARG_NUM_LIBRARIES])
	else:
		num_libraries = int(1e3)

	# Get the minimum fragment size.
	if arg_dict.has_key(ARG_MIN_FRAGMENT_SIZE):
		min_length = int(arg_dict[ARG_MIN_FRAGMENT_SIZE])
	else:
		min_length = 4

	# Get the number of fragments -- one more than the number of crossovers.
	num_fragments = int(arg_dict[ARG_NUM_CROSSOVERS])+1
	
	
	num_parents = len(parents)
	library_size = num_parents**num_fragments

	if arg_dict.has_key(ARG_MAX_CHIMERAS_PER_LIBRARY):
		max_chimeras = min(library_size, int(arg_dict[ARG_MAX_CHIMERAS_PER_LIBRARY]))
	else:
		max_chimeras = library_size

	if arg_dict.has_key(ARG_RANDOM_SEED):
		random.seed(int(arg_dict[ARG_RANDOM_SEED]))
	

	# Make libraries consistent with RASPP
	(new_parents, identical_sites) = raspp.collapse_parents(parents)
	if len(new_parents[0]) < num_fragments*min_length:
		error_msg = "Minimum diversity length of %d is too large.\n%d " + \
					"fragments with diversity %d cannot be found in a " + \
					"sequence of length %d (with identities removed).  Aborting..."
		print error_msg % (min_length, num_fragments, min_length, len(parents[0]))
		return

	
	start_time = time.clock()

	output_file.write("# <E>\t<m>\tcrossover points\n")
	random_crossovers = []
	for libnum in range(num_libraries):
		crossovers = schema.generateRandomCrossovers(len(new_parents[0]), num_fragments-1, min_length)
		crossovers = raspp.translate_collapsed_indices(crossovers, identical_sites)
		random_crossovers.append(crossovers)
	for crossovers in random_crossovers:
		fragments = schema.getFragments(crossovers, parents[0])
		filtered_contacts = schema.getSCHEMAContactsWithCrossovers(pdb_contacts, parents, crossovers)
		all_chimeras = []
		if max_chimeras < library_size:
			# Assemble a random sample of chimeras, with replacement
			for n_chim in range(max_chimeras):
				chim_index = random.randint(0,library_size-1)
				n2c = schema.base(chim_index,num_parents)
				chimera_blocks = ''.join(['1']*(num_fragments-len(n2c))+['%d'%(int(x)+1,) for x in n2c])
				all_chimeras.append(chimera_blocks)
		else: # We'll be covering all chimeras in the library; might as well get a good sample.
			# The number of parents and fragments specifies all possible chimeras, regardless of
			# crossover point positions, so pre-generate all chimeras.
			max_chimeras = library_size
			for i in range(library_size):
				# The next two lines turn i into a chimera block pattern 
				# (e.g., 0 -> '11111111', 1 -> '11111112', 2 -> '11111113'...)
				n2c = schema.base(i,num_parents)
				chimera_blocks = ''.join(['1']*(num_fragments-len(n2c))+['%d'%(int(x)+1,) for x in n2c])
				all_chimeras.append(chimera_blocks)
			# Randomly assort the chimeras
			random.shuffle(all_chimeras)

		# Calculate average E and m for the library or subsample
		E_values = []
		m_values = []

		for chim_index in range(max_chimeras):
			chimera_blocks = all_chimeras[chim_index]
			E = schema.getChimeraDisruption(chimera_blocks, filtered_contacts, fragments, parents)
			m = schema.getChimeraShortestDistance(chimera_blocks, fragments, parents)
			E_values.append(E)
			m_values.append(m)
		average_E = schema.mean(E_values)
		average_m = schema.mean(m_values)
		xover_pat = '%d '*len(crossovers)
		xover_str = xover_pat % tuple(crossovers)
		output_file.write(('%1.4f\t%1.4f\t%s\n') % (average_E, average_m, xover_str))
		output_file.flush()
	total_time = time.clock()-start_time
	output_file.write('# Finished in %1.2f seconds (%d libraries, %d chimeras)\n' % (total_time, num_libraries, num_libraries*max_chimeras))
	if arg_dict.has_key(ARG_OUTPUT_FILE):
		output_file.close()

def main_wrapper():
	main(sys.argv)

main_wrapper()
