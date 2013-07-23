#! /usr/local/bin/python
"""Script for calculating SCHEMA contacts.

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
ARG_PDB_ALIGNMENT_FILE = 'pdbal'
ARG_PARENT_INDEX = 'p'
ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE = 'msa'
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
		"\t-%s <multiple sequence alignment file>\n" % ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE, \
		"\t[-%s <PDB-parent alignment file>]\n" % ARG_PDB_ALIGNMENT_FILE,\
		"\t[-%s <PDB chains list>]\n" % ARG_CHAINS,\
		"\t[-%s <contacts output file>]" % ARG_OUTPUT_FILE

def confirm_arguments(arg_dict):
	# Are arguments okay?
	res = True
	arg_keys = arg_dict.keys()

	try:
		if len(arg_keys) == 0:
			res = False
			return
			
		if not ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE in arg_keys:
			print "  You must provide an alignment file (-%s <file>)" % ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE
			res = False
		elif not os.path.isfile(arg_dict[ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE]):
			print "  Can't find library file %s" % arg_dict[ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE]
			res = False
			
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
		if args[0].split(os.path.sep)[-1] == "schemacontacts.py":
			print_usage(args)
		return

	# Flags and values
	
	# Inputs:
	#	The PDB file name.
	pdb_file = arg_dict[ARG_PDB_FILE]
	#   The alignment/fragment file name.
	msa_file = arg_dict[ARG_MULTIPLE_SEQUENCE_ALIGNMENT_FILE]
	#	The alignment between the reference parent (indicated by reference_parent_index)
	#	and the target protein sequence in the provided PDB file.  The amino acids in
	#	the aligned reference parent should correspond exactly to those in the 
	#	msa_file above.
	# If you don't provide a PDB alignment file, the program will assume that the ID of the PDB structure
	# contained in the HEADER field corresponds to one of the sequence IDs in the MSA.
	parent_pdb_alignment_file = None
	if arg_dict.has_key(ARG_PDB_ALIGNMENT_FILE):
		if not os.path.isfile(arg_dict[ARG_PDB_ALIGNMENT_FILE]):
			print "  Can't find PDB/parent alignment file %s" % arg_dict[ARG_PDB_ALIGNMENT_FILE]
			return 
		else:
			parent_pdb_alignment_file = arg_dict[ARG_PDB_ALIGNMENT_FILE]
	else:
		pdb_key = pdb.File().getIDCode(file(pdb_file,'r'))
		
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
	
	# Read the alignment file to create a list of parents.
	# The parents will appear in the list in the order in which they appear in the file.
	parent_list = schema.readMultipleSequenceAlignmentFile(file(msa_file, 'r'))
	parent_dict = dict(parent_list)

	# Generate the contacts
	# Read in the PDB file to create a list of residues.
	residues = pdb.File().read(file(pdb_file, 'r'))
	# Because the PDB file's residue sequence may differ from those of the parents, we
	# must align the PDB residues to one parent.
	if not parent_pdb_alignment_file:  # Just get PDB sequence from the multiple sequence alignment
		try:
			aligned_pdb = parent_dict[pdb_key]
			aligned_prot = parent_dict[pdb_key]
		except KeyError:
			print "Could not find sequence %s in the multiple sequence alignment file %s.  Aborting..." % (pdb_key, msa_file)
			return
	else: # Pull information from the parent/PDB alignment file.
		# Our objective is to find the sequence with the same key in both the parent MSA file and 
		# the parent/PDB alignment file.
		pdb_parent_seq_list = schema.readMultipleSequenceAlignmentFile(file(parent_pdb_alignment_file, 'r'))
		pdb_parent_seq_dict = dict(pdb_parent_seq_list)
	
		# Bail out if there are fewer than 2 sequences.
		if len(pdb_parent_seq_dict.keys()) < 2:
			print "Only found one uniquely named sequence in the PDB/parent alignment, %s.  Aborting..." % pdb_parent_seq_dict.keys()[0]
			return

		# Find the matching key
		pdb_key = None
		for k in parent_dict.keys():
			if pdb_parent_seq_dict.has_key(k):
				pdb_key = k

		# Bail out if no matching key is found
		if not pdb_key:
			print "Could not find parents %s in PDB/parent aligned sequences %s.  Aborting..." % (parent_dict.keys(),)
			return
		aligned_prot = pdb_parent_seq_dict[pdb_key]
		# Remove the sequence corresponding to the pdb_key, leaving only the parent sequence.
		del pdb_parent_seq_dict[pdb_key]
		# Take the first remaining sequence, which should be the parent sequence.
		aligned_pdb = pdb_parent_seq_dict.values()[0]

	# Check to make sure the parent sequence from both alignment files matches.
	if aligned_prot.replace('-','') != parent_dict[pdb_key].replace('-',''):
		print "The PDB-aligned parent and the named parent, %s, don't match!  Aborting..." % (pdb_key,)
		return
	# Check to ensure the aligned PDB sequence matches the residue sequence pulled directly from the PDB file.
	if aligned_pdb.replace('-','') != pdb.sequence(residues, chain_identifiers):
		print "The parent-aligned PDB sequence, %s, and the PDB file sequence, chain(s) %s in %s, don't match!  Aborting..." % (pdb_key, chain_identifiers, pdb_file)
		return
	#print aligned_prot
	#print aligned_pdb
	#print parent_dict[pdb_key]
	#print pdb.sequence(residues)
	
	# Align the residues with the parent protein.
	try:
		residues = schema.alignPDBResidues(residues, aligned_prot, aligned_pdb, parent_dict[pdb_key], chain_identifiers)
	except ValueError, ve:
		print ve
		return
	#print pdb.sequence(residues)
	#print parent_dict[pdb_key]
	
	# 	The contact file name for output.
	if arg_dict.has_key(ARG_OUTPUT_FILE):
		contact_file = file(arg_dict[ARG_OUTPUT_FILE], 'w')
	else:
		contact_file = sys.stdout
		

	# With an aligned set of residues and parents, we can now compute the SCHEMA contacts.
	# Note that for more than two parents, some of these contacts may only be broken by 
	# specific chimera patterns.
	contact_distance = 4.5  # Residues closer than this distance, in angstroms, are in contact.
	pdb_contacts = schema.getPDBContacts(residues, contact_distance)
	schema.writeContactFile(pdb_contacts, contact_file)
	
	if not contact_file == sys.stdout:
		contact_file.close()

		
def main_wrapper():
	main(sys.argv)

main_wrapper()
