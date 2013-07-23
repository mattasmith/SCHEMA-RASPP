#! /usr/local/bin/python

import sys, os, math, string

three_to_one_map =  {'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', \
			'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L', \
			'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', \
			'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y', '-':'-'}
amino_acids = three_to_one_map.values()

class Atom:
	"Reads a PDB file ATOM line and extracts fields"
	
	# Field definitions taken from http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html
	
	def __init__(self, line):
		self.line = line
		self.x = self.getX()
		self.y = self.getY()
		self.z = self.getZ()
		self.chain = self.getChainID()
		self.atom_name = self.getAtomName()
		self.alt_loc = self.getAltLoc()
		self.res_seq = self.getResidueSequenceNumber()
		self.icode = self.getInsertionCode()
		self.residue = self.getResidue()
		self.res_id = self.getResidueID()

	def getResidue(self):
		return (self.line[17:20]).strip()
	
	def getSerialNumber(self):
		return int(self.line[6:11])
		
	def getResidueSequenceNumber(self):
		return int(self.line[22:26])
	
	def getInsertionCode(self):
		return self.line[26]

	def getCoords(self):
		return self.getX(), self.getY(), self.getZ()

	def getX(self):
		return float(self.line[30:38])

	def getY(self):
		return float(self.line[38:46])

	def getZ(self):
		return float(self.line[46:54])
		
	def getChainID(self):
		return self.line[21]

	def getAtomName(self):
		return self.line[12:16].strip()

	def getAltLoc(self):
		return self.line[16]
	
	def getDistance(self, atom):
		return math.sqrt((self.x-atom.x)**2 + (self.y-atom.y)**2 + (self.z-atom.z)**2)
	
	def getResidueID(self):
		return ('%d'%self.getResidueSequenceNumber())+self.getInsertionCode()
		
	def __repr__(self):
		return '%s %d %s (%1.3f, %1.3f, %1.3f)' % (self.residue, self.res_seq, self.atom_name, self.x, self.y, self.z)

class Residue:
	"Reads a series of lines and extracts all atoms with the same \
	residue sequence ID"
	def __init__(self, lines):
		if lines:
			atom = Atom(lines[0])
			self.residue = atom.residue
			self.res_seq = atom.res_seq
			self.res_id = atom.res_id
			self.chain = atom.chain
			self.atoms = []
			for line in lines:
				atom = Atom(line)
				if atom.res_id == self.res_id:
					self.atoms.append(atom)
				else:
					# done reading
					break
		else:
			atom = None
			self.residue = '-'
			self.res_seq = -1
			self.chain = None
			self.atoms = []
				
	def getDistance(self, residue):
		smallest_distance = 1e6 # very large number
		for atom in self.atoms:
			for other_atom in residue.atoms:
				dist = atom.getDistance(other_atom)
				if dist < smallest_distance:
					smallest_distance = dist
		return smallest_distance

	def isContact(self, residue, min_dist=4.5):
		"""Returns True if residues have at least one pair of atoms within min_dist angstroms,
		otherwise returns False."""
		atoms_to_exclude = ['H', 'O'] # exclude hydrogens and backbone oxygens
		for atom in self.atoms[1:]: # exclude residue amino terminus
			if atom.atom_name in atoms_to_exclude:
				continue
			for other_atom in residue.atoms[1:]:
				if other_atom.atom_name in atoms_to_exclude:
					continue
				dist = atom.getDistance(other_atom)
				if dist <= min_dist:
					return True
		return False
	
	def __repr__(self):
		return '%s %d' % (self.residue, self.res_seq)

class File:
	def __init__(self):
		self.residues = []
		
	def read(self, f):
		lines = f.readlines()
		atom_lines = []
		for line in lines:
			if line[0:4] == 'ATOM':
				atom_lines.append(line)
				#print pdb.Atom(line)

		while atom_lines != []:
			residue = Residue(atom_lines)
			self.residues.append(residue)
			atom_lines = atom_lines[len(self.residues[-1].atoms):]
		return self.residues
	
	def getIDCode(self, f):
		for line in f.readlines():
			if line[0:6] == 'HEADER':
				return line[62:66]
		return None

def sequence(residues, chain_ids=['A',' ']):
	return ''.join([three_to_one_map[r.residue] for r in residues if r.chain in chain_ids])
	
def get(file_name, chain_ids=['A',' ']):
	residues = File().read(file(file_name,'r'))
	print sequence(residues, chain_ids)