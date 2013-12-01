Name: Rory MacQueen
SUNETID: macqueen
SUNETID#: 05524842

Instructions on how to run program:

Tanimoto:

python tanimoto.py [drugs.csv] [targets.csv] [outputfile.csv]

All arguments are required.

P-Value:

python pvalue.py [drugs.csv] [targets.csv] [protein_A_id] [protein_B_id] -n [NUM_TIMESTEPS]

The final argument (NUM_TIMESTEPS) is optional. If none is specified, program will default to a value of 100

Networkgen:

python networkgen.py [drugs.csv] [targets.csv] [protein_nodes.csv]

All arguments are required


Note that the file data_structures.py contains helper functions used by all the scripts. That file must reside in the same directory as the above scripts in order for those scripts to work
