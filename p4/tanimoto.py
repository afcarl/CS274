import sys
from data_structures import *

# Computes the tanimoto score for two drugs, based
# on their 2D fingerprint
def get_tanimoto_score(print1, print2):
	intersect = len (print1.intersection(print2))
	union = len (print1.union(print2))
	return float(intersect) / float(union)

# Get all possible combinations of drug ids and compute the
# tanimoto score for each combination.
# Prints out that score and an indication of whether the drugs
# share a target
def compute_tanimoto_scores(drug_map, targets_map, output_file):
	output = open(output_file, 'w')
	combinations = list(itertools.combinations(sorted(drug_map.keys()), 2))
	sorted_combinations = sorted(combinations)
	for combo in sorted_combinations:
		drug1 = combo[0]
		drug2 = combo[1]
		score = get_tanimoto_score(drug_map[drug1], drug_map[drug2])
		targets1 = targets_map[drug1] if drug1 in targets_map else Set()
		targets2 = targets_map[drug2] if drug2 in targets_map else Set()
		isMatch = bool(targets1.intersection(targets2))

		to_write = ",".join([drug1, drug2, "{0:.6f}".format(score), str(int(isMatch))])
		output.write(to_write + "\n")
	output.close()

def main(drugs_file, targets_file, output_file):
	drug_map = get_drug_map(drugs_file)
	targets_map = get_target_map(targets_file)
	compute_tanimoto_scores(drug_map, targets_map, output_file)

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "Alas you clearly lack the skill to invoke my power.\n"
		print "----First argument: Drugs file, e.g. drugs.csv----\n"
		print "----Second argument: Targets file, e.g. targets.csv----\n"
		print "----Third argument: Output file, e.g. output.csv----\n"
		sys.exit()
	main(sys.argv[1], sys.argv[2], sys.argv[3])
