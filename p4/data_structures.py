from sets import Set
import itertools

# Reads in a protein_nodes file and returns
# a dictionary mapping protein ids to the 
# protein name and its indication
def get_node_map(filename):
	file_input = open(filename, 'r')
	node_map = {}
	first_line_ignore = file_input.readline()
	line = file_input.readline()
	while line:
		parts = line.strip().split(",")
		protein_id = parts[0]
		protein_name = parts[1]
		indications = parts[2]
		node_map[protein_id] = {'name':protein_name, 'indication':indications}
		line = file_input.readline()
	file_input.close()
	return node_map

# Reads in a targets file and returns
# a dictionary mapping protein ids to the 
# set of drugs that bind them
def get_protein_map(targets_file):
	file_input = open(targets_file, 'r')
	protein_map = {}
	first_line_ignore = file_input.readline()
	line = file_input.readline()
	while line:
		parts = line.split(",")
		drug_id = parts[0]
		target_id = parts[1]
		if target_id in protein_map:
			protein_map[target_id].add(drug_id)
		else:
			protein_map[target_id] = Set([drug_id])
		line = file_input.readline()
	file_input.close()
	return protein_map

# Reads in a targets file and returns
# a dictionary mapping drug ids to the 
# set of targets that they bind
def get_target_map(targets_file):
	file_input = open(targets_file, 'r')
	drug_target_map = {}
	first_line_ignore = file_input.readline()
	line = file_input.readline()
	while line:
		parts = line.split(",")
		drug_id = parts[0]
		target_id = parts[1]
		if drug_id in drug_target_map:
			drug_target_map[drug_id].add(target_id)
		else:
			drug_target_map[drug_id] = Set([target_id])
		line = file_input.readline()
	file_input.close()
	return drug_target_map

# Reads in a drugs file and returns
# a dictionary mapping drug ids to the 
# 2D fingerprint that defines it, represented as a
# set of integers
def get_drug_map(drugs_file):
	file_input = open(drugs_file, 'r')
	drug_map = {}
	first_line_ignore = file_input.readline()
	line = file_input.readline()
	while line:
		parts = line.split(",")
		drug_id = parts[0]
		fingerprint = parts[2]
		fingerprint_set = Set(fingerprint.split())
		drug_map[drug_id] = fingerprint_set
		line = file_input.readline()
	file_input.close()
	return drug_map
