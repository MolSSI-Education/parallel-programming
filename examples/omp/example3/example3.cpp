#include <stdio.h>
#include <omp.h>
#include <math.h>

//number of atoms along each side of the square lattice structure
const int nside = 40;

//total number of atoms in the square
const int natoms = nside*nside;

//initial distance between neighboring particles in the lattice
const double datoms = 1.0;

//time between molecular dynamics steps
const double dt = 0.01;

//mass of the particles
const double mass = 1.0;

//number of molecular dynamics steps to run
const int niter = 1000;

//timings information
double force_zero_time = 0.0;
double force_calc_time = 0.0;
double vel_update_time = 0.0;
double coord_update_time = 0.0;

void calculate_velocities(double coords[][2], double velocities[][2], double* energy, double* potential)
{
  
  double start_loop;

  start_loop = omp_get_wtime();
  double forces[natoms][2];
  for (int i=0; i < natoms; i++) {
    forces[i][0] = 0.0;
    forces[i][1] = 0.0;
  }
  force_zero_time += omp_get_wtime() - start_loop;

  start_loop = omp_get_wtime();
  double v = 0.0;
  for (int i=0; i < natoms; i++) {
    for (int j=i+1; j < natoms; j++) {
      double dx = coords[j][0] - coords[i][0];
      double dy = coords[j][1] - coords[i][1];
      double dr2 = dx*dx + dy*dy;
      double dr = sqrt(dr2);
      double fx = (dx/dr)*(1.0/dr2);
      double fy = (dy/dr)*(1.0/dr2);
      forces[i][0] -= fx;
      forces[i][1] -= fy;
      forces[j][0] += fx;
      forces[j][1] += fy;
      v += 1.0/dr;
    }
  }
  *potential = v;
  force_calc_time += omp_get_wtime() - start_loop;

  //update the velocities
  start_loop = omp_get_wtime();
  double ke = 0.0;
  for (int i=0; i<natoms; i++) {
    ke += 0.5*mass*(velocities[i][0]*velocities[i][0] + velocities[i][1]*velocities[i][1]);
    velocities[i][0] += (dt/mass)*forces[i][0];
    velocities[i][1] += (dt/mass)*forces[i][1];
  }
  *energy = *potential + ke;
  vel_update_time += omp_get_wtime() - start_loop;
}

int main()
{

  double start_time = omp_get_wtime();

  //coordinates of all of the atoms in 2D space
  double coords[natoms][2];

  //velocity of all of the atoms in 2D space
  double velocities[natoms][2];

  //set the coordinates
#pragma omp parallel for
  for (int i=0; i< natoms; i++) {
    coords[i][0] = (i%nside)*datoms;
    coords[i][1] = (i/nside)*datoms;
  }

  //zero the velocities
  for (int i=0; i < natoms; i++) {
    velocities[i][0] = 0.0;
    velocities[i][1] = 0.0;
  }

  //main iterative loop
  double forces[natoms][2];
  double energy;
  double potential;
  double start_loop;
  for (int iter=0; iter < niter; iter++) {

    //determine the new velocities
    calculate_velocities(coords, velocities, &energy, &potential);

    //update the coordinates
    start_loop = omp_get_wtime();
    for (int i=0; i < natoms; i++) {
      coords[i][0] += velocities[i][0]*dt;
      coords[i][1] += velocities[i][1]*dt;
    }
    coord_update_time += omp_get_wtime() - start_loop;

    //output the energy at this iteration
    printf("Iteration: %i      Energy: %f      PE: %f\n",iter,energy,potential);
  }

    printf("\n");
    printf("Timings:\n");
    printf("   Force Zero:      %f\n", force_zero_time);
    printf("   Force Calc:      %f\n", force_calc_time);
    printf("   Velocity Update: %f\n", vel_update_time);
    printf("   Coords Update:   %f\n", coord_update_time);
    printf("   Total:           %f\n", omp_get_wtime()-start_time);

    return 0;
}
