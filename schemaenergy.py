#! /usr/local/bin/python
"""Script for calculating SCHEMA energies.

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

import sys, string, os
import pdb, schema

ARG_PRINT_E = 'E'
ARG_PRINT_M = 'm'
ARG_PDB_ALIGNMENT_FILE = 'pdbal'
ARG_PARENT_INDEX = 'p'
ARG_CHIMERAS = 'chim'
ARG_CROSSOVER_FILE = 'xo'
ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE = 'msa'
ARG_CONTACT_FILE = 'con'
ARG_RANDOM_SEED = 'seed'
ARG_OUTPUT_FILE = 'o'
ARG_HELP = 'help'



def parse_arguments(args):
	# Turn linear arguments into a dictionary of (option, [values,...]) pairs
	arg_dict = {}
	key = None
	for arg in args[1:]:
		if arg[0] == '-':
			key = arg[1:]
			arg_dict[key] = None
		else:
			if arg_dict.has_key(key):
				if arg_dict[key]:
					if type(arg_dict[key]) is list:
						arg_dict[key] = arg_dict[key]+[arg]
					else:
						arg_dict[key] = [arg_dict[key],arg]
				else:
					arg_dict[key] = arg
			else:
				arg_dict[key] = arg
	return arg_dict

def print_usage(args):
	print 'Usage: python', args[0].split(os.path.sep)[-1], ' [options]'
	print 'Options:\n', \
		'\t-%s <alignment file>\n' % ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE, \
		'\t-%s <contact file>\n' % ARG_CONTACT_FILE, \
		'\t-%s <crossover file>\n' % ARG_CROSSOVER_FILE, \
		'\t[-%s <chimera list>]\n' % ARG_CHIMERAS, \
		'\t[-%s]\n' % ARG_PRINT_E, \
		'\t[-%s]\n' % ARG_PRINT_M, \
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
			print "  You must provide a library file (-%s <file>)" % ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE
			res = False
		elif not os.path.isfile(arg_dict[ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE]):
			print "  Can't find library file %s" % arg_dict[ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE]
			res = False
			
		if not ARG_CROSSOVER_FILE in arg_keys:
			print "  You must provide a crossover file (-%s <file>)" % ARG_CROSSOVER_FILE
			res = False
		elif not os.path.isfile(arg_dict[ARG_CROSSOVER_FILE]):
			print "  Can't find crossover file %s" % arg_dict[ARG_CROSSOVER_FILE]
			res = False
			
		if not ARG_CONTACT_FILE in arg_keys:
			print "  You must provide a contact file (-%s <file>)" % ARG_CONTACT_FILE
			res = False
		elif not os.path.isfile(arg_dict[ARG_CONTACT_FILE]):
			print "  Can't find contact file %s" % arg_dict[ARG_CONTACT_FILE]
			res = False
			
		if not (arg_dict.has_key(ARG_PRINT_E) or arg_dict.has_key(ARG_PRINT_M)):
			print "  No output specified; use -E to print disruption and/or -m to print mutation"
			res = False
	except Exception, e:
		raise e
		res = False
	return res

def outputEnergies(chimera_blocks, contacts, fragments, parents, output_file, output_string, print_E, print_m):
	if not schema.checkChimera(chimera_blocks, fragments, parents):
		output_file.write("# %s is not a valid chimera\n" % chimera_blocks)
		return
	output_vars = [chimera_blocks]
	E = None
	m = None
	if print_E:
		E = schema.getChimeraDisruption(chimera_blocks, contacts, fragments, parents)
		output_vars = output_vars + [E]
	if print_m:
		m = schema.getChimeraShortestDistance(chimera_blocks, fragments, parents)
		output_vars = output_vars + [m]
	#print output_vars
	output_file.write(output_string % tuple(output_vars))
	return (E,m)

def main(args):
	arg_dict = parse_arguments(args)
	if not confirm_arguments(arg_dict):
		if args[0].split(os.path.sep)[-1] == "schemaenergy.py":
			print_usage(args)
		return

	# Flags and values
	print_E = False
	print_m = False
	output_file = sys.stdout

	# Inputs:
	#   The alignment/fragment file name.
	msa_file = arg_dict[ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE]

	if arg_dict.has_key(ARG_PRINT_E):
		print_E = True
	if arg_dict.has_key(ARG_PRINT_M):
		print_m = True

	# Read the alignment file to create a list of parents.
	# The parents will appear in the list in the order in which they appear in the file.
	parent_list = schema.readMultipleSequenceAlignmentFile(file(msa_file, 'r'))
	parents = [p for (k,p) in parent_list]
	
	crossovers = schema.readCrossoverFile(file(arg_dict[ARG_CROSSOVER_FILE], 'r'))
	fragments = schema.getFragments(crossovers, parents[0])

	# Get the contacts
	pdb_contacts = schema.readContactFile(file(arg_dict[ARG_CONTACT_FILE], 'r'))
	contacts = schema.getSCHEMAContactsWithCrossovers(pdb_contacts, parents, crossovers)
	
	if arg_dict.has_key(ARG_OUTPUT_FILE):
		output_file = file(arg_dict[ARG_OUTPUT_FILE], 'w')

	# Now, what does the user want?
	output_string = '%s'
	output_file.write('# chimera')
	if print_E:
		output_string += '\t%d'
		output_file.write('\tE')
	if print_m:
		output_string += '\t%d'
		output_file.write('\tm')
	output_string += '\n'
	output_file.write('\n')
	
	if arg_dict.has_key(ARG_CHIMERAS): # Print values for chimeras
		chimeras = arg_dict[ARG_CHIMERAS]
		# Could be a) a chimera, b) a list of chimeras, or c) a file of chimeras.
		if type(chimeras) is list:
			# It's a list of chimeras
			for chimera_blocks in chimeras:
				outputEnergies(chimera_blocks, contacts, fragments, parents, output_file, output_string, print_E, print_m)
		elif os.path.isfile(chimeras):
			# It's a file of chimeras
			for line in file(chimeras,'r').readlines():
				chimera_blocks = line.strip()
				outputEnergies(chimera_blocks, contacts, fragments, parents, output_file, output_string, print_E, print_m)
		else:
			# It's a single chimera sequence
			chimera_blocks = chimeras
			outputEnergies(chimera_blocks, contacts, fragments, parents, output_file, output_string, print_E, print_m)
	else:
		# Enumerates all possible chimeras and their disruption and mutation values.
		p = len(parents)
		n = len(fragments)
		Es = []
		ms = []
		for i in xrange(len(parents)**len(fragments)):
			# The next two lines turn i into a chimera block pattern 
			# (e.g., 0 -> '11111111', 1 -> '11111112', 2 -> '11111113'...)
			n2c = schema.base(i,p)
			chimera_blocks = ''.join(['1']*(n-len(n2c))+['%d'%(int(x)+1,) for x in n2c])
			(E, m) = outputEnergies(chimera_blocks, contacts, fragments, parents, output_file, output_string, print_E, print_m)
			if (print_E):
				Es.append(E)
			if (print_m):
				ms.append(m)
		if (print_E):
			mean_str = "# Average disruption <E> = %1.4f\n" % schema.mean(Es)
			output_file.write(mean_str)
		if (print_m):
			mean_str = "# Average mutation <m> = %1.4f\n" % schema.mean(ms)
			output_file.write(mean_str)
	
	if arg_dict.has_key(ARG_OUTPUT_FILE):
		output_file.close()


def main_wrapper():
	main(sys.argv)

main_wrapper()
