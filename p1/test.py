import align
import os
import sys

def get_info(filename):
	info = {}
	s = set()
	file_input = open(filename, 'r')
	info['score'] = float(file_input.readline())
	file_input.readline()
	while True:
		line = file_input.readline()
		if line == None or line == "":
			break
		line2 = file_input.readline()
		s.add((line, line2))
		file_input.readline()
	info['matches'] = s
	file_input.close()
	return info


def equal_outputs(output1, output2):
	info1 = get_info(output1)
	info2 = get_info(output2)
	return info1 == info2

def main(folder):
	files = os.listdir(folder)
	for filename in files:
		filename = folder + "/" + filename
		basename = filename.split(".")[0]
		ext = filename.split(".")[-1]
		if ext == "input":
			align.main(filename, basename + ".test_output")
			if not equal_outputs(basename + ".output", basename + ".test_output"):
				print "WE HAVE A PROBLEM WITH: " + filename
			else:
				print filename + " : PASSED"

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print "Alas! You clearly lack the skill to invoke my power.\n"
		print "----First Argument: Folder of examples----"
		sys.exit()
	main(sys.argv[1])
