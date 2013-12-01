import sys
import random
from math import sqrt

# Helper function to compute mean of values
def mean(values):
	return sum(values) / len(values)

# This class represents a centroid in the KMeans algorithm
class Centroid:
	def __init__(self, vector, index):
		self.vector = vector

		# We maintain a set of old points and new points 
		# so that we can detect if we have converged by checking
		# to see if new_points differs from old_points
		self._new_points = []
		self._old_points = []
		self._index = index

	# Adds a point to this centroid.
	# Also returns boolean specifying if this point already existed
	# in this centroid. Necessary for checking convergence
	def add_point(self, point):
		already_had = point in self._old_points
		self._new_points.append(point)
		return already_had

	def get_points(self):
		return self._old_points

	def get_index(self):
		return self._index

	# This method computes the center of this centroid's
	# points
	def get_center_of_points(self):
		vector_value = zip(*self._old_points)[1]
		center = map(mean, zip(*vector_value))
		return center

	def update_points(self):
		self._old_points = self._new_points
		self._new_points = []

	# Computes the distance between a point and this centroid's center
	# based on Euclidean distance 
	def distance_to_me(self, point):
		if not len(point) == len(self.vector):
			print "Trying to compare points of different dimension. Cannot do that"
			sys.exit()
		total_sum = 0.0
		for i in range(len(point)):
			valueA = point[i]
			valueB = self.vector[i]
			total_sum += (valueB - valueA)**2
		return sqrt(total_sum)

	def __str__(self):
		string = ""
		string += "Centroid is centered at: " + str(self.vector)
		string += "\nPoints in centroid: \n"
		for point in self._old_points:
			string += str(point) + "\n"
		return string

# This class represents the KMeans algorithm
class KMeans:
	def __init__(self, k, expression_data, centroids):
		self._k = k
		# Store the index of each sample and each centroid
		self._expression_data = [[index, item] for index, item in enumerate(expression_data)]
		self._centroids = [Centroid(value, index) for index, value in enumerate(centroids)]

	# Gets the closest centroid to a point
	def _get_closest_centroid(self, point):
		closest_centroid = self._centroids[0]
		min_distance = float("inf")
		for centroid in self._centroids:
			distance = centroid.distance_to_me(point)
			if distance < min_distance:
				min_distance = distance
				closest_centroid = centroid
		return closest_centroid

	# Assigns points to centroids. Along the way, checks for 
	# convergence of the algorithm
	def _make_assignments_to_centroids(self):
		converged = True
		for point in self._expression_data:
			closest_centroid = self._get_closest_centroid(point[1])
			already_had_this_point = closest_centroid.add_point(point)
			if not already_had_this_point:
				converged = False
		return converged

	# Recomputes the center of each centroid and sets it
	def _recompute_centroids(self):
		for centroid in self._centroids:
			centroid.update_points()
			center = centroid.get_center_of_points()
			centroid.vector = center

	# Runs the algorithm for specified number of iterations
	# or until algorithm converges
	def run(self, num_iterations):
		for i in range(num_iterations):
			has_converged = self._make_assignments_to_centroids()
			if has_converged:
				print "Converged after " + str(i) + " iterations"
				break
			self._recompute_centroids()
	
	def get_results(self):
		result_arr = [[] for i in range(len(self._expression_data))]
		for centroid in self._centroids:
			for point in centroid.get_points():
				item_index = point[0]
				result_arr[item_index] = [item_index, centroid.get_index()]
		return result_arr

	@staticmethod
	def generate_random_centroids(data, num_centroids):
		return random.sample(data, num_centroids)

# Helper function to get the data from file
def get_data_from_file(filename):
	file_input = open(filename, 'r')
	line = file_input.readline()
	matrix = []
	while line:
		arr = line.split()
		# Convert all strs to floats here
		floats = [float(x) for x in arr]
		matrix.append(floats)
		line = file_input.readline()
	file_input.close()
	return matrix

def print_results(results, output_filename):
	file_output = open(output_filename, 'w')
	to_print = ""
	for item in results:
		to_print += str(item[0]+1) + "\t" + str(item[1]+1) + "\n"
	file_output.write(to_print)
	print to_print
	file_output.close()

def main(k, data_file, num_iterations, centroids_file):
	expression_data = get_data_from_file(data_file)
	if not centroids_file:
		centroids = KMeans.generate_random_centroids(expression_data, k)
	else:
		centroids = get_data_from_file(centroids_file)
	kmeans = KMeans(k, expression_data, centroids)
	kmeans.run(num_iterations)
	results = kmeans.get_results()
	print_results(results, "kmeans.out")

if __name__ == "__main__":
	centroids_file = None
	if len(sys.argv) == 5:
		centroids_file = sys.argv[-1]
		del sys.argv[-1]
	if len(sys.argv) != 4:
		print "Alas you clearly lack the skill to invoke my power.\n"
		print "----First argument: Integer value for 'k' number of centroids, e.g. 5----\n"
		print "----Second argument: Expression data file, e.g. expression.dat----\n"
		print "----Third argument: Integer value for num iterations allowed, e.g. 10----\n"
		print "----Fourth argument: Centroids file e.g. centroids.txt----\n"
		sys.exit()
	main(int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), centroids_file)
