#! /usr/local/bin/python

import sys, os, math, string, random

def random_contact_map(length, num_contacts):
	# Generate a random contact map
	residues = range(length)
	contacts = []
	pairs = []
	for i in range(num_contacts):
		# choose random contacts
		[i, j] = random.sample(residues, 2)
	
		if not (i,j) in pairs:
			if i < j:
				contacts.append((i,j,i,j))
			else:
				contacts.append((j,i,j,i))
			pairs.append((i,j))
			pairs.append((j,i))
	return contacts

aas = 'ACDEFGHIKLMNPQRSTVWY'
def random_protein(length):
	prot = []
	for i in range(length):
		prot.append(aas[random.randint(0,len(aas)-1)])
	return ''.join(prot)

def random_protein_with_identity(template_prot, seq_id):
	prot = []
	for i in range(len(template_prot)):
		if random.random() < seq_id:
			prot.append(template_prot[i])
		else:
			prot.append(aas[random.randint(0,len(aas)-1)])
	return ''.join(prot)

