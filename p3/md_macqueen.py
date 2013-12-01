import argparse
from math import sqrt
import sys

NUM_DIMENSIONS = 3

# Gets the distance between two points
def get_distance(pointA, pointB):
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

class AtomicInteraction:
	def __init__(self, atom1_id, atom2_id, distance, isBond):
		self.atom1_id = atom1_id
		self.atom2_id = atom2_id
		self.isBond = isBond
		self.ref_distance = distance
		self.force = 0

	def __str__(self):
		result = "Atom1 id: " + str(self.atom1_id) + "\t"
		result += "Atom2 id: " + str(self.atom2_id) + "\t"
		result += "isBond: " + str(self.isBond) + "\n"
		result += "ref distance: " + str(self.ref_distance) + "\n"
		return result

class Atom:
	def __init__(self, id_num, position_x, position_y, position_z, velocity_x, velocity_y, velocity_z, mass, bonded_atoms):
		# Use zero indexing
		self.id = id_num - 1
		self.position_vec = [position_x, position_y, position_z]
		self.velocity_vec = [velocity_x, velocity_y, velocity_z]
		self.acceleration_vec = [0, 0, 0]
		self.force_vec = [0, 0, 0]
		self.mass = mass
		# Use zero indexing
		self.bonded_atoms = [id_num - 1 for id_num in bonded_atoms]
		self.interactions = []

	def add_interaction(self, interaction):
		self.interactions.append(interaction)

	def get_bonded_atoms(self):
		return self.bonded_atoms

	def get_position(self):
		return self.position_vec

	def get_velocity(self):
		return self.velocity_vec

	def get_acceleration(self):
		return self.acceleration_vec

	def set_new_velocity(self, delta_t):
		self.velocity_vec = [self.velocity_vec[i] + 0.5*self.acceleration_vec[i]*delta_t for i in range(NUM_DIMENSIONS)]

	def set_new_position(self, delta_t):
		self.position_vec = [self.position_vec[i] + self.velocity_vec[i]*delta_t for i in range(NUM_DIMENSIONS)]

	def set_new_acceleration(self):
		self.acceleration_vec = [self.force_vec[i] / self.mass for i in range(NUM_DIMENSIONS)]

	def add_new_force(self, force, other_atom_pos):
		distance = get_distance(self.position_vec, other_atom_pos)
		self.force_vec = [self.force_vec[i] + force * ((other_atom_pos[i] - self.position_vec[i])/distance) for i in range(NUM_DIMENSIONS)]


	def __str__(self):
		result = "id: " + str(self.id) + "\n"
		result += "position: " + str(self.position_vec) + "\n"
		result += "velocity: " + str(self.velocity_vec) + "\n"
		result += "acceleration: " + str(self.acceleration_vec) + "\n"
		result += "force: " + str(self.force_vec) + "\n"
		result += "bonded atoms: " + str(self.bonded_atoms) + "\n"
		return result

class MolecularDynamics:
	def __init__(self, atoms, bonds, bond_spring_k, nonbond_spring_k, atomic_mass, delta_t):
		self.atoms = atoms
		self.bonds = bonds
		self.bond_spring_k = bond_spring_k
		self.nonbond_spring_k = nonbond_spring_k
		self.atomic_mass = atomic_mass
		self.delta_t = delta_t

	def set_new_velocities_and_positions(self):
		for atom in self.atoms:
			atom.set_new_velocity(self.delta_t)
			atom.set_new_position(self.delta_t)

	def set_new_velocities(self):
		for atom in self.atoms:
			atom.set_new_velocity(self.delta_t)

	def set_new_accelerations(self):
		for atom in self.atoms:
			atom.set_new_acceleration()

	def calculate_potential_energy(self):
		bond_potential_energy = 0.0
		non_bond_potential_energy = 0.0
		for interaction in self.bonds:
			atom1 = self.atoms[interaction.atom1_id]
			atom2 = self.atoms[interaction.atom2_id]
			distance = get_distance(atom1.get_position(), atom2.get_position())

			if interaction.isBond:
				force = self.bond_spring_k * (distance - interaction.ref_distance)
				potenial_energy = 0.5 * self.bond_spring_k * (distance - interaction.ref_distance)**2
				bond_potential_energy += potenial_energy

			else:
				force = self.nonbond_spring_k * (distance - interaction.ref_distance)
				potenial_energy = 0.5 * self.nonbond_spring_k * (distance - interaction.ref_distance)**2
				non_bond_potential_energy += potenial_energy
			interaction.force = force
		return [bond_potential_energy, non_bond_potential_energy]

	def calculate_kinetic_energy(self):
		kinetic_energy = 0.0
		for atom in self.atoms:
			total_velocity = sum([atom.get_velocity()[i] ** 2 for i in range(NUM_DIMENSIONS)])
			kinetic_energy += (0.5 * atom.mass * total_velocity)

		return kinetic_energy

	def calculate_forces(self):
		for atom in self.atoms:
			atom.force_vec = [0, 0, 0]
			for interaction in atom.interactions:
				otherAtomId = interaction.atom1_id if interaction.atom1_id != atom.id else interaction.atom2_id
				otherAtom = self.atoms[otherAtomId]
				atom.add_new_force(interaction.force, otherAtom.get_position())

	def run_simulation(self, num_iterations, erg_filename, rvc_filename, first_line):

		erg_output = open(erg_filename, 'w')
		header = "# step\tE_k\tE_b\tE_nB\tE_tot\n"
		erg_output.write(header)

		rvc_output = open(rvc_filename, 'w')
		rvc_output.write(first_line)
		self.write_out_atoms(rvc_output, 0, 0)

		inital_energy = lower_bound = upper_bound = 0
		for run in range(1, num_iterations+1):
			self.set_new_velocities_and_positions()
			bond_energy, non_bond_energy = self.calculate_potential_energy()
			self.calculate_forces()
			self.set_new_accelerations()
			self.set_new_velocities()
			kinetic_energy = self.calculate_kinetic_energy()
			total_energy = bond_energy + non_bond_energy + kinetic_energy
			
			# Record initial energy so we can check for instability
			# in subsequent runs
			if run == 1:
				inital_energy = total_energy
				lower_bound = inital_energy / 10.0
				upper_bound = inital_energy * 10.0

			if total_energy < lower_bound or total_energy > upper_bound:
				print "Oh dear! We are on run " + str(run) + " and this simulation has become most unstable!"
				print "Inital energy was: " + "{0:.1f}".format(inital_energy)
				print "Current energy is: " + "{0:.1f}".format(total_energy)
				print "I fear I must exit! Ta-ta!"
				sys.exit()

			if run % 10 == 0:
				to_write = str(run+1) + "\t" 
				to_write += "{0:.1f}".format(kinetic_energy) + "\t"
				to_write += "{0:.1f}".format(bond_energy) + "\t"
				to_write += "{0:.1f}".format(non_bond_energy) + "\t"
				to_write += "{0:.1f}".format(total_energy) + "\n"
				erg_output.write(to_write)

				self.write_out_atoms(rvc_output, run, total_energy)

		erg_output.close()
		rvc_output.close()

	def print_atom_dist(self):
		i = 17
		j = 96
		atom1 = self.atoms[i-1]
		atom2 = self.atoms[j-1]
		dist = get_distance(atom1.get_position(), atom2.get_position())
		print "{0:.4f}".format(dist)

	def write_out_atoms(self, output, step, energy):
		to_write = ""
		# Only add this header when not the inital print out
		if step > 0:
			to_write += "#At time step: " + str(step) + ", energy = " + "{0:.3f}".format(energy) + "kJ\n"
		for atom in self.atoms:
			position = atom.get_position()
			to_write += str(atom.id + 1) + "\t"
			to_write += "{0:.4f}".format(position[0]) + "\t"
			to_write += "{0:.4f}".format(position[1]) + "\t"
			to_write += "{0:.4f}".format(position[2]) + "\t"

			velocity = atom.get_velocity()

			to_write += "{0:.4f}".format(velocity[0]) + "\t"
			to_write += "{0:.4f}".format(velocity[1]) + "\t"
			to_write += "{0:.4f}".format(velocity[2]) + "\t"

			for bonded_atom_id in atom.bonded_atoms:
				to_write += str(bonded_atom_id + 1) + "\t"
			to_write = to_write.strip() + "\n"

		output.write(to_write)


def get_atoms_from_file(filename, atomic_mass):
	file_input = open(filename, 'r')

	# Skip the first line which is just a header for the file
	first_line = file_input.readline()
	out_first_line = "\t".join(first_line.split()) + "\n"
	line = file_input.readline()
	atoms = []
	while line:
		arr = line.split()
		atom = Atom(int(arr[0]), float(arr[1]), float(arr[2]), float(arr[3]), float(arr[4]), float(arr[5]), float(arr[6]), atomic_mass, [int(x) for x in arr[7:]])
		atoms.append(atom)
		line = file_input.readline()

	file_input.close()
	return [atoms, out_first_line]

def get_args(args):
	parser = argparse.ArgumentParser(description='Run molecular dynamics simulation')
	parser.add_argument("--if", dest='if', help='Input filename', required=True)
	parser.add_argument('--n', type=int, dest='n', help='Number of timesteps for program to iterate', default=1000)
	parser.add_argument('--kB', type=float, dest = 'kB', help='Bond spring constant', default=40000.0)
	parser.add_argument('--kN', type=float, dest = 'kN', help='Non-bond spring constant', default=400.0)
	parser.add_argument('--nbCutoff', type=float, dest='nbCutoff', help = 'Cutoff for nonbond interaction', default=0.50)
	parser.add_argument('--m', type=float, dest='m', help='Mass of atom', default= 12.0)
	parser.add_argument('--dt', type=float, dest='dt', help='Length of one timestep', default=0.001)
	parser.add_argument('--out', dest='out', help='Output filename', default="")

	# If out is empty string, make it the input filename
	options = vars(parser.parse_args())
	if not options['out']:
		# Remove .rvc
		options['out'] = options['if'].split(".")[-2]
	return options

def main (args):
	options = get_args(args)

	nbCutoff = options['nbCutoff']

	input_filename = options['if']
	atoms, first_line = get_atoms_from_file(input_filename, options['m'])
	number_of_atoms = len(atoms)
	interactions = []
	for index_i in range(number_of_atoms):
		atom1 = atoms[index_i]
		bonds = atom1.get_bonded_atoms()
		for index_j in range(index_i + 1, number_of_atoms):
			atom2 = atoms[index_j]
			distance = get_distance(atom1.get_position(), atom2.get_position())
			if index_j in bonds:
				interaction = AtomicInteraction(atom1.id, atom2.id, distance, True)
				interactions.append(interaction)
				atom1.add_interaction(interaction)
				atom2.add_interaction(interaction)
			else:
				if distance < nbCutoff:
					interaction = AtomicInteraction(atom1.id, atom2.id, distance, False)
					interactions.append(interaction)
					atom1.add_interaction(interaction)
					atom2.add_interaction(interaction)

	molecular_dynamics = MolecularDynamics(atoms, interactions, options['kB'], options['kN'], options['m'], options['dt'])
	erg_filename = options['out'] + "_out.erg"
	rvc_filename = options['out'] + "_out.rvc"
	molecular_dynamics.run_simulation(options['n'], erg_filename, rvc_filename, first_line)

if __name__ == '__main__':
	main (sys.argv)
