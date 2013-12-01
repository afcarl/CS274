import sys
from data_structures import *
from pvalue import compute_p_value

P_VALUE_THRESHOLD = 0.05

# Takes in a label (either 'name' or 'indication')
# and prints each node_id along with the corresponding
# label
def print_node_info(node_map, label):
	file_output = open(label + ".nodeAttr", 'w')
	file_output.write(label + "\n")
	for key, value in node_map.items():
		val_to_print = value[label]
		file_output.write(key + " = " + val_to_print + "\n")
	file_output.close()

# Prints edges in our network. For every pair of
# protein_ids in our data, we calculate the p_value
# between them. If the p_value is less than our
# threshold (i.e. significant) then we print the
# edge between those nodes
def print_edges(node_map, drug_map, protein_map):
	file_output = open('network.sif', 'w')
	combinations = list(itertools.combinations(sorted(node_map.keys()), 2))
	sorted_combinations = sorted(combinations)
	for combo in sorted_combinations:
		protein1 = combo[0]
		protein2 = combo[1]
		p_value = compute_p_value(protein_map, drug_map, protein1, protein2, 100)
		if p_value <= P_VALUE_THRESHOLD:
			file_output.write(protein1 + " edge " + protein2 + "\n")
	file_output.close()

def main(drugs_file, targets_file, nodes_file):
	drug_map = get_drug_map(drugs_file)
	protein_map = get_protein_map(targets_file)
	node_map = get_node_map(nodes_file)

	print_edges(node_map, drug_map, protein_map)
	print_node_info(node_map, 'name')
	print_node_info(node_map, 'indication')

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "Alas you clearly lack the skill to invoke my power.\n"
		print "----First argument: Drugs file, e.g. drugs.csv----\n"
		print "----Second argument: Targets file, e.g. targets.csv----\n"
		print "----Third argument: Protein nodes file, e.g. protein_nodes.csv----\n"
		sys.exit()
	main(sys.argv[1], sys.argv[2], sys.argv[3])
