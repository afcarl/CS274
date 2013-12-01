# Author: Rory MacQueen <rorymacqueen@gmail.com>
# This is free software! I hereby give anyone permission
# to use, copy, modify, and redistribute this code in any way,
# but please include a citation to the author.
#
# This code implements the Smith-Waterman algorithm for aligning
# sequences. It handles both global and local alignment

import sys

# DELTA value used for evaluating equality of floating point numbers
DELTA = 0.001

# Helper function to evaluate equality of floating point numbers
def float_equality(floatA, floatB):
	return abs(floatA - floatB) < DELTA

# AlignmentCell:
#
# This class represents a cell in the Score matrices (M, Ix, Iy).
# Each cell has its own value, as well as an array of outgoing pointers,
# which is used to track the traceback of an alignment.
class AlignmentCell:

	def __init__ (self, value, pointers, owner, i_coord, j_coord):
		self._value = value
		self._pointers = pointers
		self._owner_matrix = owner
		self._i = i_coord
		self._j = j_coord

	def get_i_coord(self):
		return self._i

	def get_j_coord(self):
		return self._j

	def get_value(self):
		return self._value

	def set_value(self, new_value):
		self._value = new_value

	def get_pointers(self):
		return self._pointers

	def add_pointer(self, cell_to_point_to):
		self._pointers.append(cell_to_point_to)

	def get_owner_matrix(self):
		return self._owner_matrix

	# This method will be invoked when you call str(cell).
	# It will print out the relevant properties of a cell for debugging
	# purposes as well as for identifying the cells in a traceback.
	def __str__ (self):
		result = self._owner_matrix.get_name()
		result += ": | (" + str(self._i) + ", " + str(self._j) + ") " + str (self._value) + "|"
		result += " pointers: " + str([str(ptr._i) + "," + str(ptr._j) + " " + ptr._owner_matrix.get_name() for ptr in self._pointers])
		return result

# AlignmentMatrix:
#
# This class represents a score matrix. When initalized
# we pass in an sentinel value for the first column (col_init)
# and the first row (row_init). These sentinel values vary depending
# on which matrix we are building (e.g. M, Ix, Iy)
class AlignmentMatrix:

	def __init__ (self, num_rows, num_cols, name, row_init, col_init):
		self._num_rows = num_rows
		self._num_cols = num_cols
		self._name = name
		self._grid = [[AlignmentCell(0.0, [], self, i, j) for j in range(self._num_cols)] for i in range(self._num_rows)]

		# Initialize the boundaries to appropriate values
		self.set_cell_value(0, 0, float(0))
		for i in range(1, self._num_rows):
			self.set_cell_value(i, 0, row_init)
		for j in range(1, self._num_cols):
			self.set_cell_value(0, j, col_init)


	def get_cell(self, i, j):
		return self._grid[i][j]

	def set_cell_value(self, i, j, value):
		self.get_cell(i, j).set_value(value)

	def get_num_rows(self):
		return self._num_rows

	def get_num_cols(self):
		return self._num_cols

	def get_name(self):
		return self._name

	# Returns a list of candidate cells with which to start backtrace
	def get_backtrace_start_cells(self, isLocal):
		start_cells = []

		# If local, every cell is candidate for starting backtrace.
		# If global, then only the bottom roq and far-right column cells
		# are eligible. We leave bottom-right cell till the end so
		# that we don't double count it
		if isLocal:
			for i in range(0, self._num_rows):
				for j in range(0, self._num_cols):
					start_cells.append(self.get_cell(i, j))
		else:
			for i in range(1, self._num_rows - 1):
				start_cells.append(self.get_cell(i, self._num_cols - 1))
			for j in range(1, self._num_cols - 1):
				start_cells.append(self.get_cell(self._num_rows - 1, j))
			start_cells.append(self.get_cell(self._num_rows - 1, self._num_cols - 1))

		return start_cells

	# Helper function to print a matrix, for debugging
	def print_me(self):
		result = ""
		for i in range (len(self._grid)):
			for j in range(len(self._grid[0])):
				cell = self.get_cell(i, j)
				result += str(cell.get_value()) + " "
			result += "\n"
		print result

# AlignmentProblem:
#
# This class represents a sequence alignment task.
# It is initalized with the two sequences in question (sequenceA and sequenceB)
# as well as various other input information which is given to it.

class AlignmentProblem:

	def __init__ (self, sequenceA, sequenceB, isLocal, gap_penalties, match_matrix_map):
		self._sequenceA = sequenceA
		self._sequenceB = sequenceB
		self._isLocal = isLocal
		self._gap_penalties = gap_penalties
		self._gap_open_penalty_x = float(gap_penalties[0])
		self._gap_extension_penalty_x = float(gap_penalties[1])
		self._gap_open_penalty_y = float(gap_penalties[2])
		self._gap_extension_penalty_y = float(gap_penalties[3])
		self._match_matrix_map = match_matrix_map
		self._final_results = set()

		self._matrixM = AlignmentMatrix(len(sequenceA)+1, len(sequenceB)+1, "M", float("-inf"), float("-inf"))
		self._matrixIx = AlignmentMatrix(len(sequenceA)+1, len(sequenceB)+1, "Ix", float(0), float("-inf"))
		self._matrixIy = AlignmentMatrix(len(sequenceA)+1, len(sequenceB)+1, "Iy", float("-inf"), float(0))

	def get_results(self):
		return self._final_results

	def populate_matrices(self):
		# start at (1,1) since row 0 and column 0 are reserved
		# for sentinel values.
		num_rows = self._matrixM.get_num_rows()
		num_cols = self._matrixM.get_num_cols()
		for i in range(1, num_rows):
			for j in range(1, num_cols):
				self.populate_cell(i, j)

	def populate_cell(self, i, j):
		# We populate all three matrices at the same time
		self._populate_matrixM_cell(i, j)
		self._populate_matrixIx_cell(i, j)
		self._populate_matrixIy_cell(i, j)

	# Calculate winning cell out of contenders
	# Return winning cells and associated winning value.
	# Contenders must be an array of tuples, where each
	# tuple is (Cell, value)
	def _get_winning_cells(self, contenders):
		winner = max(contenders, key=lambda x: x[1])
		winners = []
		# We have to make sure to return all cells that have
		# a value matching that of the winning cell
		for contender in contenders:
			if float_equality(contender[1], winner[1]):
				winners.append(contender[0])
		return [winners, winner[1]]

	# These next three functions implement the core of the
	# Smith-Waterman algorithm. Each matrix is tracking each case
	# of the alignment: 1) a match between sequence A and sequence B.
	# 2) a character in sequence A aligned to a gap. 3) A character
	# in sequence B aligned to a gap.

	def _populate_matrixM_cell(self, i, j):
		# Need to subtract 1 from i, j here because
		# the first row and first column is reserved for
		# prepopulated, sentinel values
		char_match = self._sequenceA[i-1] + self._sequenceB[j-1]
		match_score = self._match_matrix_map[char_match]

		m_contender = self._matrixM.get_cell(i-1, j-1)
		ix_contender = self._matrixIx.get_cell(i-1, j-1)
		iy_contender = self._matrixIy.get_cell(i-1, j-1)

		winners, winning_value = self._get_winning_cells([
			[m_contender, m_contender.get_value() + match_score],
			[ix_contender, ix_contender.get_value() + match_score],
			[iy_contender, iy_contender.get_value() + match_score]])

		# Recall that in the local case, cells in a score matrix
		# must have values >= 0
		if winning_value < 0 and self._isLocal:
			winning_value = 0

		self._matrixM.set_cell_value(i, j, winning_value)
		for winner in winners:
			self._matrixM.get_cell(i, j).add_pointer(winner)

	def _populate_matrixIx_cell(self, i, j):
		m_contender = self._matrixM.get_cell(i-1, j)
		ix_contender = self._matrixIx.get_cell(i-1, j)

		winners, winning_value = self._get_winning_cells([
			[m_contender, m_contender.get_value() - self._gap_open_penalty_y],
			[ix_contender, ix_contender.get_value() - self._gap_extension_penalty_y]])

		if winning_value < 0 and self._isLocal:
			winning_value = 0

		self._matrixIx.set_cell_value(i, j, winning_value)
		for winner in winners:
			self._matrixIx.get_cell(i, j).add_pointer(winner)

	def _populate_matrixIy_cell(self, i, j):
		m_contender = self._matrixM.get_cell(i, j-1)
		iy_contender = self._matrixIy.get_cell(i, j-1)

		winners, winning_value = self._get_winning_cells([
			[m_contender, m_contender.get_value() - self._gap_open_penalty_x],
			[iy_contender, iy_contender.get_value() - self._gap_extension_penalty_x]])

		if winning_value < 0 and self._isLocal:
			winning_value = 0

		self._matrixIy.set_cell_value(i, j, winning_value)
		for winner in winners:
			self._matrixIy.get_cell(i, j).add_pointer(winner)

	# This method initiates the backtrace component of the algorithm.
	# We collect the start cells, and begin the recursive backtrace
	# from those cells.
	def perform_backtrace(self):
		start_cells = []
		start_cells += self._matrixM.get_backtrace_start_cells(self._isLocal) + self._matrixIx.get_backtrace_start_cells(self._isLocal) + self._matrixIy.get_backtrace_start_cells(self._isLocal)
		contenders = [[cell, cell.get_value()] for cell in start_cells]
		max_start_cells, max_value = self._get_winning_cells(contenders)
		self.max_match_value = max_value
		for cell in max_start_cells:
			self.recursive_backtrace(cell, "", "")

	# This method represents the backtrace through the matrices to reconstruct
	# an alignment. It is a recursive method. The base case depends on whether
	# this is a local or global alignment. For a global alignment, the base case
	# is when we reach a cell which has no outgoing pointers (i.e. a cell on the
	# edge of a score matrix). For a local alignment, the base case is when we hit
	# a cell which has a value of zero.

	def recursive_backtrace(self, cell, seqA_result, seqB_result):
		pointers = cell.get_pointers()
		matrix_name = cell.get_owner_matrix().get_name()

		# We can stop either if we have hit a cell with no more pointers
		# or if we are local and hit a cell with value 0
		if not pointers or (self._isLocal and float_equality(cell.get_value(), 0.0)):
			self._final_results.add((seqA_result, seqB_result))
			return

		# If we are in the M matrix, we align character i from sequence A with
		# character j from sequence B.
		# If we are in the Ix matrix, we align character i from sequence A with
		# a gap.
		# If we are in the Iy matrix, we align character j from sequence B with
		# a gap.

		if matrix_name == "M":
			seqA_result = self._sequenceA[cell.get_i_coord()-1] + seqA_result
			seqB_result = self._sequenceB[cell.get_j_coord()-1] + seqB_result
		if matrix_name == "Ix":
			seqA_result = self._sequenceA[cell.get_i_coord()-1] + seqA_result
			seqB_result = "_" + seqB_result
		if matrix_name == "Iy":
			seqA_result = "_" + seqA_result
			seqB_result = self._sequenceB[cell.get_j_coord()-1] + seqB_result
		
		# Invoke the recursive step on all cells pointed to by the current cell
		for next_cell in pointers:
			self.recursive_backtrace(next_cell, seqA_result, seqB_result)

# This function will actually print the results of the
# algorithm to an output file using the format
# specified in the handout

def print_results(matches, value, output_filename):
	file_output = open(output_filename, 'w')
	file_output.write(str(value))
	file_output.write("\n\n")
	for match in matches:
		file_output.write(match[0] + "\n")
		file_output.write(match[1] + "\n")
		file_output.write("\n")
	file_output.close()

# The main function handles reading in the
# file and setting up the alignment problem.
# It then calls methods to run the algorithm
# and print out the results.

def main(input_filename, output_filename):
	file_input = open(input_filename, 'r')
	sequenceA = file_input.readline().strip()
	sequenceB = file_input.readline().strip()
	isLocal = bool(int(file_input.readline()))
	gap_penalties = file_input.readline().split()
	size_of_A_alphabet = int(file_input.readline())
	A_alphabet = file_input.readline()
	size_of_B_alphabet = int(file_input.readline())
	B_alphabet = file_input.readline()

	# Note I represent the match matrix not as a matrix
	# per se but as simply a mapping for character pairs
	# to their match values. I found this to be an easy
	# and straighforward way of storing those matches.
	match_matrix_map = {}
	for i in range(0, size_of_A_alphabet * size_of_B_alphabet):
		entry = file_input.readline().split()
		key = entry[2] + entry[3]
		match_matrix_map[key] = float(entry[4])

	align_problem = AlignmentProblem(sequenceA, sequenceB, isLocal, gap_penalties, match_matrix_map)
	align_problem.populate_matrices()
	align_problem.perform_backtrace()
	print_results(align_problem.get_results(), align_problem.max_match_value, output_filename)

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print "Alas! You clearly lack the skill to invoke my power. The usage is:\n"
		print "----First Argument: Name of input file----"
		print "----Second Argument: Name of output file----"
		sys.exit()
	main(sys.argv[1], sys.argv[2])