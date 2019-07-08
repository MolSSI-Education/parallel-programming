#include <iostream>
#include <cmath>
#include <cstdio>
#include <list>
#include <vector>
#include <utility>
#include <algorithm>
#include <omp.h>
#include <mpi.h>

// Global data

const int natom = 8000;         // Number of atoms
const double sigma = 20;	// Particle radius
const double L = 1600;		// Box size
const double epsilon = 1.0;	// Binding energy
const double target_temp = 0.4; // Target temperature
const int nprint = 100;		// Print & temp scale every 100 steps
const int nneigh = 20;		// Recompute neighbor list every 20 steps
const int nstep=1000;		// Number of steps to take

const double r2cut_neigh = 35;  // Cut off for r-squared when computing neighbor list
const double r2cut_force = 30;  // Cut off for r-squared when computing forces
const double excess_vel = 1.6;  // Initialize with hot atoms

double time_force, time_neigh, time_total; // Timers

int nproc, me;
bool PRINT = false;

typedef std::pair<int,int> pairT;
typedef std::pair<double,double> xyT;
typedef std::vector<pairT> neighT;
typedef std::vector<xyT> coordT;
typedef std::vector<neighT> thrneighT;

// Enforce periodic boundary conditions
inline double periodic(double x, double L) {
    while (x>L) x-=L;
    while (x<0) x+=L;
    return x;
}

// Make the list of interacting atoms 
void neighbor_list(const coordT& coords, thrneighT& thrneigh) {
    double start = omp_get_wtime();
    //for (int ithr=0; ithr<thrneigh.size(); ithr++) {
#pragma omp parallel default(none) shared(thrneigh, coords, nproc, me)
    {
        int ithr = omp_get_thread_num();
	int nthread = omp_get_num_threads();
	neighT& neigh = thrneigh[ithr];
	neigh.clear();
	neigh.reserve(100*natom/nthread);
       	//for (int i=ithr; i<natom; i+=thrneigh.size()) {
	for (int i=(me*nthread)+ithr; i<natom; i+=nproc*nthread) {
	  double xi = coords[i].first;
	  double yi = coords[i].second;
	  for (int j=0; j<i; j++) {
            double xj = coords[j].first;
            double yj = coords[j].second;
            double dx = (xi-xj);
            double dy = (yi-yj);
            
            if (dx > (L/2)) dx = dx - L;
            else if (dx < (-L/2)) dx = dx + L;
            
            if (dy > (L/2)) dy = dy - L;
            else if (dy < (-L/2)) dy = dy + L;
            
            double r2 = (dx*dx+dy*dy)/(sigma*sigma);
            if (r2 < r2cut_neigh) {
	      neigh.push_back(pairT(i,j));
            }
	  }
	}
#pragma omp barrier
#pragma omp single
	{
	  //get the target number of pairs per thread
          int npair=0;
	  for (int i=0; i<nthread; i++) npair += thrneigh[i].size();
          npair = (npair-1)/nthread + 1;

	  for (int i=0; i<thrneigh.size()-1; i++) {
	    while(thrneigh[i].size() < npair) {
	      //steal pairs from the next thread
	      if (thrneigh[i+1].size() == 0) break;
	      thrneigh[i].push_back(thrneigh[i+1].back());
	      thrneigh[i+1].pop_back();
	    }
	    while(thrneigh.size() > npair) {
	      //donate pairs from the next thread
	      if (thrneigh[i].size() == 0) break;
	      thrneigh[i+1].push_back(thrneigh[i].back());
	      thrneigh[i].pop_back();
	    }
	  }
	}
    }
    time_neigh += omp_get_wtime() - start;
}

// Compute the forces, virial and potential energy
coordT forces(const thrneighT& thrneigh, const coordT& coords, double& virial, double& pe) {
    double start = omp_get_wtime();
    coordT f(natom,xyT(0.0,0.0));

    virial = pe = 0.0;

    // V(|ri-rj|) = epsilon*((sigma/r)^12 - 2*(sigma/r)^6)
    // dV/dxi = -12*epsilon*((sigma/r)^14 - (sigma/r)^8)*(xi-xj)/sigma**2
    // F[i][x] = -dV/dxi

    const double fac = epsilon*12.0/(sigma*sigma);
#pragma omp parallel default(none) shared(f, thrneigh, coords, virial, pe)
    {
      coordT f_thread(natom,xyT(0.0,0.0));
      double virial_thread = 0.0;
      double pe_thread = 0.0;
      int ithr = omp_get_thread_num();
      const neighT& neigh = thrneigh[ithr];
      for (neighT::const_iterator ij=neigh.begin(); ij!=neigh.end(); ++ij) {
        int i = ij->first;
        int j = ij->second;
        
        double xi = coords[i].first;
        double yi = coords[i].second;
        
        double xj = coords[j].first;
        double yj = coords[j].second;
        
        double dx = (xi-xj);
        double dy = (yi-yj);

        if (dx > (L/2)) dx = dx - L;
        else if (dx < (-L/2)) dx = dx + L;
        
        if (dy > (L/2)) dy = dy - L;
        else if (dy < (-L/2)) dy = dy + L;
        
        double r2 = (dx*dx + dy*dy)/(sigma*sigma);

        if (r2 < r2cut_force) {
	  double r6 = r2*r2*r2;
	  double r12 = r6*r6;
	  double vij = epsilon*(1.0/r12 - 2.0/r6);
	  double df = fac*(1.0/(r12*r2) - 1.0/(r6*r2));
	  double dfx = df*dx;
	  double dfy = df*dy;
            
	  f_thread[i].first += dfx;
	  f_thread[j].first -= dfx;
	  f_thread[i].second += dfy;
	  f_thread[j].second -= dfy;
            
	  pe_thread += vij;
	  virial_thread += dfx*dx + dfy*dy;
        }
      }
#pragma omp barrier
      //int ithr = omp_get_thread_num();
      int nthreads = omp_get_num_threads();
      int per_thread = natom / nthreads;
      for (int i=0; i<nthreads; i++) {
	int istart = ( (i+ithr) % nthreads ) * per_thread;
	int iend;
	if ( (i+ithr) % nthreads == nthreads - 1 ) {
	  iend = natom;
	  pe += pe_thread;
	  virial += virial_thread;
	}
	else {
	  iend = ( (i+ithr+1) % nthreads ) * per_thread;
	}
	for (int j=istart; j<iend; j++) {
	  f[j].first += f_thread[j].first;
	  f[j].second += f_thread[j].second;
	}
#pragma omp barrier
      }
    }

    // Reduce over MPI ranks
    coordT ftotal(natom, xyT(0.0,0.0));
    MPI_Allreduce(&f[0], &ftotal[0], 2*natom, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double virial_total, pe_total;
    MPI_Allreduce(&virial, &virial_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&pe, &pe_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    virial = virial_total; // To return values via the argument list
    pe = pe_total;

    time_force += omp_get_wtime() - start;
    return ftotal;
}

// Constrain a to be in [-b,b]
inline double restrict(double a, double b) {
    if (a > b) return b;
    else if (a < -b) return -b;
    else return a;
}

// Relax the initial random guess of atoms positions using steepest descent
void optimize(coordT& coords, thrneighT& thrneigh) {
    double dt = 0.1;
    double prev = 1e99;
    for (int step=0; step<600; step++) {
      if ((step%(3*nneigh)) == 0 || step<10) neighbor_list(coords, thrneigh);
        double virial,pe;
        coordT f = forces(thrneigh,coords,virial,pe);
        for (int i=0; i<natom; i++) {
            double x = coords[i].first;
            double y = coords[i].second;
            double fx= restrict(f[i].first, 2.0);
            double fy= restrict(f[i].second,2.0);

            coords[i] = xyT(periodic(x+dt*fx,L),periodic(y+dt*fy,L));
        }

        if ((step%50)==0 && PRINT) 
            std::cout << "optim: " <<  pe << " " <<  dt << std::endl;
        
        if (std::abs(pe-prev) < prev*0.01) break;
        prev = pe;
    }
}

// A simple pseudo-random number generator
double drand() {
    static const unsigned int a = 1664525;
    static const unsigned int c = 1013904223;
    static unsigned int x = 23111;
    static const double fac = 2.3283064365386963e-10;

    x = a*x + c;
    
    return fac*x;
}


// The main molecular dynamics program
void md() {
    const double dt=0.03;
    if (PRINT) std::cout << "Time step " << dt << std::endl;

    // initialize random coords and velocities
    coordT coords(natom, xyT(0.0,0.0));
    coordT v(natom, xyT(0.0,0.0));

    double vxmean=0.0, vymean= 0.0;

    double box = std::min(std::sqrt(1.0*natom)*sigma*1.25,L);
    for (int i=0; i<natom; i++) {
        double xi, yi;
        for (int attempt=0; attempt<5; attempt++) {
            xi = box*(drand()-0.5) + L/2.0;
            yi = box*(drand()-0.5) + L/2.0;
            double r2min = 1000000.0;
            for (int j=0; j<i; j++) {
                double xj = coords[j].first;
                double yj = coords[j].second;
                double dx = (xi-xj);
                double dy = (yi-yj);
                r2min = std::min(r2min,dx*dx+dy*dy);
            }
            if (r2min > 0.125*sigma*sigma) break;
        }
        coords[i]=xyT(xi,yi);
        double vx = (drand()-0.5)*std::sqrt(2.0*target_temp)*2.0*excess_vel;
        double vy = (drand()-0.5)*std::sqrt(2.0*target_temp)*2.0*excess_vel;
        vxmean += vx;
        vymean += vy;
        v[i] = xyT(vx,vy);
    }
    
    vxmean /= natom;
    vymean /= natom;
    
    // Adjust so the mean velocity is zero
    for (int i=0; i<natom; i++) {
        v[i].first -= vxmean;
        v[i].second -= vymean;
    }

    // Broadcast the coordinates and velocities
    if (MPI_Bcast((void*) &coords[0], natom*sizeof(xyT), MPI_BYTE, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
      throw("broadcast of coords failed");
    if (MPI_Bcast((void*) &v[0], natom*sizeof(xyT), MPI_BYTE, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
      throw("broadcast of vels failed");

    // Relax the initial random guess
    int nthreads;
#pragma omp parallel
    {
      nthreads = omp_get_num_threads();
    }
    thrneighT thrneigh(nthreads);
    optimize(coords, thrneigh);
    neighbor_list(coords, thrneigh);

    // make the initial forces
    double virial = 0.0;
    double temp = 0.0;

    double potential_energy;
    coordT f = forces(thrneigh,coords,virial,potential_energy);

    // Finally!  Do the molecular dynamics
    int step_stats = 0;
    for (int step=1; step<=nstep; step++) {

        // update the velocities to time t+dt/2, and the positions to time t+dt
        for (int i=0; i<natom; i++) {
            double vx = v[i].first, vy = v[i].second;
            double Fx = f[i].first, Fy = f[i].second;
            double x = coords[i].first, y = coords[i].second;
            vx += Fx*dt*0.5;
            vy += Fy*dt*0.5;
            x += vx*dt;
            y += vy*dt;
            v[i] = xyT(vx,vy);
            coords[i] = xyT(periodic(x,L),periodic(y,L));  // periodic boundary conditions
        }
        // make the forces at time t+dt
        if ((step%nneigh) == 0) {
	  neighbor_list(coords, thrneigh);
        }

        double virial_step;
        f = forces(thrneigh,coords,virial_step,potential_energy);
        virial += virial_step;

        // finish update of v to time t+dt
        double kinetic_energy = 0.0;
        for (int i=0; i<natom; i++) {
            double vx = v[i].first, vy = v[i].second;
            double Fx = f[i].first, Fy = f[i].second;
            vx += Fx*dt*0.5;
            vy += Fy*dt*0.5;
            v[i] = xyT(vx, vy);
            kinetic_energy += 0.5*(vx*vx+vy*vy);
        }

        temp += kinetic_energy/(natom - 1.0);
        step_stats += 1;

        if ((step%nprint) == 0) {

            if (step == nprint && PRINT) {
                printf("\n");
                printf("    time         ke            pe             e            T          P\n");
                printf("  -------    -----------   ------------  ------------    ------    ------\n");
            }

            temp = temp/step_stats;
            virial = 0.5*virial/step_stats;
            double pressure = (natom*temp + virial)/(L*L);
            double energy = kinetic_energy + potential_energy;

            // Temperature scaling for the first (nstep/3) steps ... indicated with a * on printing
            double vscale = std::sqrt(target_temp/temp);
            if (step>=(nstep/3)) vscale=1.0;
            const char scaling[2]={' ','*'};

            if (PRINT) printf("%9.2f   %12.5f   %12.5f  %12.5f %8.3f %12.8f %c\n",
                   step*dt, kinetic_energy, potential_energy, energy, temp, 
                   pressure, scaling[vscale!=1.0]);

            for (int i=0; i<natom; i++) {
                v[i].first *= vscale;
                v[i].second *= vscale;
            }

            temp = virial = 0.0;
            step_stats = 0;
        }
    }
    
}

int main(int argc, char** argv) {
    if (MPI_Init(&argc,&argv) != MPI_SUCCESS) 
      throw "MPI init failed";

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    PRINT = (me == 0);

    time_force = time_neigh = time_total = 0.0;
    double start = omp_get_wtime();

    md();

    time_total += omp_get_wtime() - start;

    if (PRINT) printf("times:  force=%.2fs  neigh=%.2fs  total=%.2fs\n", 
	   time_force, time_neigh, time_total);

    MPI_Finalize();
    return 0;
}
