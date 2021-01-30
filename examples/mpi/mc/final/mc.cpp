#include <iostream>
#include <math.h>
#include <random>
#include <mpi.h>

// Initialize the random number generator with a pre-defined seed
//std::mt19937 mt(1);

// Initialize the random number generator with a random seed
std::random_device rd;
std::mt19937 mt(rd());

std::uniform_real_distribution<double> dist(0.0, 1.0);

// Generate an initial set of coordinates
int generate_initial_state(int num_particles, double box_length, double *coordinates) {

  // Fill coordinates with random values
  for (int i=0; i<3*num_particles; i++) {
    coordinates[i] = ( 0.5 - dist(mt) ) * box_length;
  }

  return 0;
}

// Evaluate the LJ potential for a given squared distance
double lennard_jones_potential(double rij2) {
  double sig_by_r2 = 1.0 / rij2;
  double sig_by_r6 = sig_by_r2*sig_by_r2*sig_by_r2;
  double sig_by_r12 = sig_by_r6*sig_by_r6;
  return 4.0 * ( sig_by_r12 - sig_by_r6 );
}

// Compute the standard tail energy correction for the LJ potential
double calculate_tail_correction( double box_length, double cutoff, int number_particles ) {
  double volume = box_length*box_length*box_length;
  double sig_by_cutoff1 = 1.0 / cutoff;
  double sig_by_cutoff3 = sig_by_cutoff1*sig_by_cutoff1*sig_by_cutoff1;
  double sig_by_cutoff9 = sig_by_cutoff3*sig_by_cutoff3*sig_by_cutoff3;

  double e_correction = sig_by_cutoff9 - 3.0 * sig_by_cutoff3;
  e_correction *= 8.0 / 9.0 * M_PI * number_particles / volume * number_particles;
  return e_correction;
}

// Compute the minimum image distance between two particles
double minimum_image_distance(double *r_i, double *r_j, double box_length) {
  double rij[3];
  rij[0] = r_i[0] - r_j[0] - box_length * round( (r_i[0] - r_j[0]) / box_length );
  rij[1] = r_i[1] - r_j[1] - box_length * round( (r_i[1] - r_j[1]) / box_length );
  rij[2] = r_i[2] - r_j[2] - box_length * round( (r_i[2] - r_j[2]) / box_length );

  double rij2 = ( rij[0] * rij[0] ) + ( rij[1] * rij[1] ) + ( rij[2] * rij[2] );
  return rij2;
}

// Compute the energy of a particle
double get_particle_energy(double *coordinates, int particle_count, double box_length, int i_particle, double cutoff2, MPI_Comm comm) {
  // Get information about the MPI communicator
  int my_rank, world_size;
  MPI_Comm_size(comm, &world_size);
  MPI_Comm_rank(comm, &my_rank);
  
  double e_total = 0.0;
  double *i_position = &coordinates[3*i_particle];

  for (int j_particle=my_rank; j_particle < particle_count; j_particle += world_size) {
    if ( i_particle != j_particle ) {
      double *j_position = &coordinates[3*j_particle];
      double rij2 = minimum_image_distance( i_position, j_position, box_length );
      if ( rij2 < cutoff2 ) {
	e_total += lennard_jones_potential(rij2);
      }
    }
  }

  // Sum the energy across all ranks
  double e_summed = 0.0;
  MPI_Allreduce(&e_total, &e_summed, 1, MPI_DOUBLE, MPI_SUM, comm);
  return e_summed;
}

// Calculate the total energy of the system
double calculate_total_pair_energy(double *coordinates, int particle_count, double box_length, double cutoff2) {
  double e_total = 0.0;
  for ( int i_particle=0; i_particle < particle_count; i_particle++ ) {
    for ( int j_particle=0; j_particle < i_particle; j_particle++ ) {
      double *r_i = &coordinates[3*i_particle];
      double *r_j = &coordinates[3*j_particle];
      double rij2 = minimum_image_distance( r_i, r_j, box_length );
      if ( rij2 < cutoff2 ) {
	e_total += lennard_jones_potential(rij2);
      }
    }
  }
  return e_total;
}

// Accept or reject a move based on the energy difference and system temperature
bool accept_or_reject( double delta_e, double beta ) {
  bool accept = false;
  if ( delta_e < 0.0 ) {
    accept = true;
  }
  else {
    double random_number = dist(mt);
    double p_acc = exp(-beta * delta_e);

      if ( random_number < p_acc ) {
	accept = true;
      }
      else {
	accept = false;
      }
  }
  return accept;
}

// Change the acceptance criteria to get the desired rate
double adjust_displacement( int n_trials, int n_accept, double max_displacement ) {
  double acc_rate = double(n_accept) / double(n_trials);
  double new_displacement = max_displacement;
  if ( acc_rate < 0.38 ) {
    new_displacement *= 0.8;
  }
  else if ( acc_rate > 0.42 ) {
    new_displacement *= 1.2;
  }
  return new_displacement;
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  double start_simulation_time = MPI_Wtime();
  double total_energy_time = 0.0;
  double total_decision_time = 0.0;

  int world_size, my_rank;
  MPI_Comm world_comm = MPI_COMM_WORLD;
  MPI_Comm_size(world_comm, &world_size);
  MPI_Comm_rank(world_comm, &my_rank);

  /******************
  * Parameter setup *
  ******************/

  double reduced_temperature = 0.9;
  double reduced_density = 0.9;
  int n_steps = 500000;
  int freq = 1000;
  int num_particles = 100;
  double simulation_cutoff = 3.0;
  double max_displacement = 0.1;
  bool tune_displacement = true;
  bool plot = true;

  double box_length = cbrt(num_particles / reduced_density);
  double beta = 1.0 / reduced_temperature;
  double simulation_cutoff2 = simulation_cutoff*simulation_cutoff;
  int n_trials = 0;
  int n_accept = 0;
  double *energy_array = new double[n_steps];

  double *coordinates = new double[3*num_particles];

  /*************************
  * Monte Carlo Simulation *
  *************************/
  if ( my_rank == 0 ) {
    generate_initial_state(num_particles, box_length, coordinates);
  }
  MPI_Bcast(coordinates, 3*num_particles, MPI_DOUBLE, 0, world_comm);

  double total_pair_energy = calculate_total_pair_energy(coordinates, num_particles, box_length, simulation_cutoff2);
  double tail_correction = calculate_tail_correction(box_length, simulation_cutoff, num_particles);

  // Beginning of main MC iterative loop
  n_trials = 0;
  for (int i_step=0; i_step<n_steps; i_step++) {

    int i_particle;
    double random_displacement[3];
    if ( my_rank == 0 ) {
      n_trials += 1;
      i_particle = floor( double(num_particles) * dist(mt) );
      for (int i=0; i<3; i++) {
        random_displacement[i] = ( ( 2.0 * dist(mt) ) - 1.0 ) * max_displacement;
      }
    }
    MPI_Bcast(&i_particle, 1, MPI_INT, 0, world_comm);
    MPI_Bcast(coordinates, 3*num_particles, MPI_DOUBLE, 0, world_comm);
    MPI_Bcast(random_displacement, 3, MPI_DOUBLE, 0, world_comm);


    // get the current energy of the test particle
    double start_energy_time = MPI_Wtime();
    double current_energy = get_particle_energy( coordinates, num_particles, box_length, i_particle, simulation_cutoff2, world_comm );
    total_energy_time += MPI_Wtime() - start_energy_time;

    // get the new coordinates of the test particle
    for (int i=0; i<3; i++) {
      coordinates[3*i_particle + i] += random_displacement[i];
      coordinates[3*i_particle + i] -= box_length * round(coordinates[3*i_particle + i] / box_length);
    }

    // get the new energy of the test particle
    start_energy_time = MPI_Wtime();
    double proposed_energy = get_particle_energy( coordinates, num_particles, box_length, i_particle, simulation_cutoff2, world_comm );
    total_energy_time += MPI_Wtime() - start_energy_time;

    if ( my_rank == 0 ) {
      // test whether to accept or reject this step
      double start_decision_time = MPI_Wtime();
      double delta_e = proposed_energy - current_energy;
      bool accept = accept_or_reject(delta_e, beta);
      if (accept) {
	total_pair_energy += delta_e;
	n_accept += 1;
      }
      else {
	// revert the position of the test particle
	for (int i=0; i<3; i++) {
	  coordinates[3*i_particle + i] -= random_displacement[i];
	  coordinates[3*i_particle + i] -= box_length * round(coordinates[3*i_particle + i] / box_length);
	}
      }

      double total_energy = (total_pair_energy + tail_correction) / double(num_particles);
      energy_array[i_step] = total_energy;

      if ( (i_step+1) % freq == 0 ) {
	if ( my_rank == 0 ) {
	  std::cout << i_step + 1 << "   " << energy_array[i_step] << std::endl;
	}

	if ( tune_displacement ) {
	  max_displacement = adjust_displacement(n_trials, n_accept, max_displacement);
	  n_trials = 0;
	  n_accept = 0;
	}
      }

      total_decision_time += MPI_Wtime() - start_decision_time;
    }
  }

  if ( my_rank == 0 ) {
    std::cout << "Total simulation time: " << MPI_Wtime() - start_simulation_time << std::endl;
    std::cout << "    Energy time:       " << total_energy_time << std::endl;
    std::cout << "    Decision time:     " << total_decision_time << std::endl;
  }
  
  delete [] energy_array;
  delete [] coordinates;
  MPI_Finalize();
  return 0;
}
