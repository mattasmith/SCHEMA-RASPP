#! /usr/local/bin/python
"""Script for extracting sequences from a PDB file.

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

ARG_PDB_FILE = 'pdb'
ARG_OUTPUT_FILE = 'o'
ARG_HELP = 'help'
ARG_CHAINS = 'chains'

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
	print 'Usage: python', args[0].split(os.path.sep)[-1], " [options]"
	print 'Options:\n',\
		"\t-%s <PDB file>\n" % ARG_PDB_FILE, \
		"\t[-%s <PDB chain list, e.g. A B C>]\n" % ARG_CHAINS,\
		"\t[-%s <contacts output file>]" % ARG_OUTPUT_FILE

def confirm_arguments(arg_dict):
	# Are arguments okay?
	res = True
	arg_keys = arg_dict.keys()

	try:
		if len(arg_keys) == 0:
			res = False
			return
			
		if not ARG_PDB_FILE in arg_keys:
			print "  You must provide a PDB file (-%s <file>)" % ARG_PDB_FILE
			res = False
		elif not os.path.isfile(arg_dict[ARG_PDB_FILE]):
			print "  Can't find PDB file %s" % arg_dict[ARG_PDB_FILE]
			res = False
		
	except:
		res = False
	return res

def main(args):
	arg_dict = parse_arguments(args)
	if not confirm_arguments(arg_dict):
		if args[0].split(os.path.sep)[-1] == "pdbseq.py":
			print_usage(args)
		return

	# Flags and values
	
	# Inputs:
	#	The PDB file name.
	pdb_file = arg_dict[ARG_PDB_FILE]

	# The PDB chains
	# Many PDB files include multiple chains.  The chain_identifier list includes those
	# chains which correspond to the protein whose contacts are being evaluated.
	# Most often, chain 'A' (in the case of multiple chains) or chain ' ' (only one chain)
	# will be the appropriate choice.
	if arg_dict.has_key(ARG_CHAINS):
		chains = arg_dict[ARG_CHAINS]
		if type(chains) is list:
			chain_identifiers = chains + [' ']
		else:
			chain_identifiers = [chains, ' ']
	else:
		chain_identifiers = ['A',' ']
	
	# 	The file name for output.
	if arg_dict.has_key(ARG_OUTPUT_FILE):
		output_file = file(arg_dict[ARG_OUTPUT_FILE], 'w')
	else:
		output_file = sys.stdout
	
	# Read in the PDB file to create a list of residues.
	residues = pdb.File().read(file(pdb_file, 'r'))
	
	# Filter residues not in selected chains
	residue_seq = pdb.sequence(residues, chain_identifiers)
	if residue_seq == '':
		print "No residues found for chain(s) %s.  Aborting..." % chain_identifiers
		return
	
	# Print it
	output_file.write('# Residue sequence for chain(s) %s from PDB file %s\n%s' % \
		(chain_identifiers, pdb_file, residue_seq))
	
	if arg_dict.has_key(ARG_OUTPUT_FILE):
		output_file.close()

		
def main_wrapper():
	main(sys.argv)

main_wrapper()
