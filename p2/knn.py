import sys
from random import shuffle
from math import sqrt

# Helper function to transpose matrix
def transpose_matrix(arr):
	return zip(*arr)

# Note, to write this function, I looked
# at a stackoverflow post on how to evenly
# split an array into chunks.
# http://stackoverflow.com/questions/2130016/splitting-a-list-of-arbitrary-size-into-only-roughly-n-equal-parts
def chunkIt(arr, num_chunks):
  average_size = len(arr) / float(num_chunks)
  result = []
  index = 0.0
  while index < len(arr):
    result.append(arr[int(index):int(index + average_size)])
    index += average_size

  return result

# Generates folds to use for cross-validation
def make_folds(matrix, num_folds):
	# Randomize the position of patients in the matrix
	shuffle(matrix)
	folds = chunkIt(matrix, num_folds)
	return folds

# Class to represent KNearestNeighbors problem
class KNearestNeighbors:
	def __init__(self, k, p):
		self._k = k
		self._p = p

	# Gets the distance between two points
	def _get_distance(self, pointA, pointB):
		if not len(pointA) == len(pointB):
			print "Trying to compare points of different dimension. Cannot do that"
			sys.exit()
		total_sum = 0.0
		# Computes Euclidean distance by getting sum of squares of distances
		for i in range(len(pointA)):
			valueA = pointA[i]
			valueB = pointB[i]
			total_sum += (valueB - valueA)**2
		return sqrt(total_sum)

	# Gets the k closest points to test_point out of trained_points
	def _get_closest_k_points(self, test_point, trained_points):
		points_with_distance = []
		# Gets distance to all points
		for point in trained_points:
			distance = self._get_distance(test_point[0], point[0])
			points_with_distance.append([point, distance])

		# Sort all points by distance and take the top k
		sorted_points = sorted(points_with_distance, key=lambda x: x[1])
		points_without_distances = [x[0] for x in sorted_points]
		return points_without_distances[0:self._k]

	# Gets prediction of class for a point by adding
	# up the votes of the k closest points
	def _get_prediction(self, closest_k):
		pos_votes = neg_votes = 0.0
		for point in closest_k:
			if point[1] == 1:
				pos_votes += 1
			elif point[1] == 0:
				neg_votes += 1
			else:
				print "ERROR: This guy isn't voting correctly"
				print point[1]
		pos_fraction = pos_votes / (pos_votes + neg_votes)
		# If fraction of positive votes is greater than threshold p,
		# then classify it as positive
		if pos_fraction >= self._p:
			return 1
		else:
			return 0

	# Main function for this class. Takes in positive and negative 'training' data.
	# Aggregates data into one training set.
	def train(self, pos_train_data, neg_train_data, pos_test_data, neg_test_data):
		# Combine folds of training data, and label each point appropriately
		labeled_points = []
		for data_point in pos_train_data:
			labeled_points.append([data_point, 1])
		for data_point in neg_train_data:
			labeled_points.append([data_point, 0])

		# Combine positive and negative test data with labels
		labeled_test_points = []
		for test_point in pos_test_data:
			labeled_test_points.append([test_point, 1])
		for test_point in neg_test_data:
			labeled_test_points.append([test_point, 0])

		# For each point, get closet k neighbours. Get prediction based
		# on those neighbours. Compare prediction to correct label
		# in order to calculate performance statistics
		tp=tn=fp=fn = 0.0
		for point in labeled_test_points:
			closest_k = self._get_closest_k_points(point, labeled_points)
			prediction = self._get_prediction(closest_k)
			correct_label = point[1]
			if prediction == 1 and correct_label == 1:
				tp += 1
			elif prediction == 1 and correct_label == 0:
				fp += 1
			elif prediction == 0 and correct_label == 1:
				fn += 1
			elif prediction == 0 and correct_label == 0:
				tn += 1

		return [tp, tn, fp, fn]

# Helper function to read data from file
def get_patients_from_file(filename):
	file_input = open(filename, 'r')
	line = file_input.readline()
	matrix = []
	while line:
		arr = line.split()
		# Convert all strs to floats here
		floats = [float(x) for x in arr]
		matrix.append(floats)
		line = file_input.readline()
	# Need to transpose matrix because we want rows to be samples
	# and columns to be features
	matrix = transpose_matrix(matrix)
	file_input.close()
	return matrix

# Helper function to print results
def print_results(sensitivity, specificity, accuracy, k, p, output_filename):
	file_output = open(output_filename, 'w')
	to_print = "k: " + str(k) + "\n" + "p: " + str(p) + "\n"
	to_print += "accuracy: " + "{0:.2f}".format(accuracy) + "\n"
	to_print += "sensitivity: " + "{0:.2f}".format(sensitivity) + "\n"
	to_print += "specificity: " + "{0:.2f}".format(specificity)
	file_output.write(to_print)
	print to_print
	file_output.close()

# Main function which reads data and runs algorithm
def main(pos_file, neg_file, k, p):
	all_patients = get_patients_from_file(pos_file)
	aml_patients = get_patients_from_file(neg_file)
	num_folds = 4
	all_folds = make_folds(all_patients, num_folds)
	aml_folds = make_folds(aml_patients, num_folds)

	total_tp = total_tn = total_fp = total_fn = 0.0
	for i in range(0, num_folds):
		# Splits training data to remove the one fold that we are going to test on
		to_train_all = all_folds[:i] + all_folds[i+1:]
		to_train_aml = aml_folds[:i] + aml_folds[i+1:]

		#aggregates the training folds together
		training_data_all = [item for sublist in to_train_all for item in sublist]
		training_data_aml = [item for sublist in to_train_aml for item in sublist]

		# Pick out the one fold we are going to test on
		to_test_all = all_folds[i]
		to_test_aml = aml_folds[i]

		# Run algorithm
		knn_problem = KNearestNeighbors(k, p) 
		tp, tn, fp, fn = knn_problem.train(training_data_all, training_data_aml, to_test_all, to_test_aml)

		total_tp += tp
		total_tn += tn
		total_fp += fp
		total_fn += fn

	# Calculate performance statistics to print out
	sensitivity = float(total_tp) / (total_tp + total_fn)
	specificity = float(total_tn) / (total_tn + total_fp)
	accuracy = float(total_tp + total_tn) / (total_tp + total_tn + total_fp + total_fn)

	print_results(sensitivity, specificity, accuracy, k, p, 'knn.out')


if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "Alas you clearly lack the skill to invoke my power.\n"
		print "----First argument: Positive training set file, e.g. ALL.dat----\n"
		print "----Second argument: Negative training set file, e.g. AML.dat----\n"
		print "----Third argument: Integer value for 'k', e.g. 5----\n"
		print "----Fourth argument: float value for 'p', e.g. 0.4----\n"
		sys.exit()
	main(sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4]))
