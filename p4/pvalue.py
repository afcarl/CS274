import sys
import argparse
import random
from data_structures import *
from tanimoto import get_tanimoto_score

# cutoff at which we accept tanimoto scores
THRESHOLD = 0.5

# Helper function to parse command line arguments
def get_args(args):
	parser = argparse.ArgumentParser(description='Run molecular dynamics simulation')
	parser.add_argument("drugs_csv", help='Drugs file')
	parser.add_argument("targets_csv", help='Targets file')
	parser.add_argument("proteinA", help='Protein A')
	parser.add_argument("proteinB", help='Protein B')
	parser.add_argument('--n', type=int, dest='n', help='Number of timesteps for program to iterate', default=100)

	options = vars(parser.parse_args())
	return options

# Computes the T_Summary value for a pair of proteins
# This is done by calculating the tanimoto score for all
# pairs of drugs across the two ligand sets, and summing
# those that are greater than the 0.5 threshold
def compute_t_summary(ligand_set_1, ligand_set_2, drug_print_map):
	total = 0.0
	for drug1 in ligand_set_1:
		for drug2 in ligand_set_2:
			result = get_tanimoto_score(drug_print_map[drug1], drug_print_map[drug2])
			if result > THRESHOLD:
				total += result
	return total

# Computes the p_value for a given pair of proteins.
# First we get the T_Summary for the proteins using their
# real ligand sets. Then we iteratively choose random sets 
# of ligands of the same size and compute the t_summary
# for each of those samples. If the randomly generated 
# t_summary is greater than our 'true' one then we add one
# to a counter. The p_value returned is the counter/num_iterations
def compute_p_value(protein_map, drug_print_map, proteinA, proteinB, num_samples):
	t_summary = compute_t_summary(protein_map[proteinA], protein_map[proteinB], drug_print_map)
	ligand_set_a_size = len(protein_map[proteinA])
	ligand_set_b_size = len(protein_map[proteinB])
	unique_ligands = set().union(*protein_map.values())
	counter = 0
	for i in range(num_samples):
		sample_set_a = random.sample(unique_ligands, ligand_set_a_size)
		sample_set_b = random.sample(unique_ligands, ligand_set_b_size)
		t_summary_check = compute_t_summary(sample_set_a, sample_set_b, drug_print_map)
		if t_summary_check >= t_summary:
			counter += 1
	return float(counter) / float(num_samples)

def main (args):
	options = get_args(args)
	drug_print_map = get_drug_map(options['drugs_csv'])
	protein_map = get_protein_map(options['targets_csv'])
	p_value = compute_p_value(protein_map, drug_print_map, options['proteinA'], options['proteinB'], options['n'])
	print p_value

if __name__ == '__main__':
	main (sys.argv)
