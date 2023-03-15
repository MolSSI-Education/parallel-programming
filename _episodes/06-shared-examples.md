# OpenMP Hands-On

````{admonition} Overview
:class: overview

Questions:
- How can I use OpenMP to parallelize a code?

Objectives:
- Use OpenMP to implement shared-memory parallelization.
- Learn how to identify and fix race conditions.
- Optimize the performance of an OpenMP code.
````


## Example 1

From the top-level directory of the GitHub repository, do:

**$ cd examples/omp/example1**

If you open [example1.cpp](https://github.com/MolSSI-Education/parallel-programming/blob/gh-pages/examples/omp/example1/example1.cpp) with a text editor, you will see that it is a simple Hello World code.

``` cpp
#include <stdio.h>

int main()
{
  printf("Hello World\n");
  
  return 0;
}
```

Go ahead and build the code.

**$ cmake .  
$ make**

You could run the code just by typing `./example1` into the command line, but there is also a simple script that does the same thing.
Run the script now.

**$ ./run.sh**

> Hello World!

We will now turn this into a multi-threaded code.
To use OpenMP, you will need to add `#include <omp.h>` to the beginning of your code.

Then, you can add OpenMP parallelization through the use of a compiler directive, which is any line that begins with `#pragma` (short for "pragmatic information").
The directive `#pragma omp parallel` tells the compiler that the next block of code should be parallelized using OpenMP.
We are going to parallelize the line with the `printf`, so put curly brackets around that line in order to make it a distinct block of text, then put `#pragma omp parallel` before the block.
Your code should look like:

``` cpp
#include <stdio.h>
#include <omp.h>

int main()
{
#pragma omp parallel
  {
    printf("Hello World\n");
  }
  
  return 0;
}
```

Go ahead and rebuild the code.

**$ cmake .  
$ make**

Before running the code, we should also specify how many OpenMP threads we want to run on.
This is something that you decide at run-time, by setting a special system variable called `OMP_NUM_THREADS`.
Add a line to `run.sh` so that it sets `OMP_NUM_THREADS` to 4:

~~~
export OMP_NUM_THREADS=4
./example1
~~~

Now let's run the code:

**$ ./run.sh**

>Hello World!  
>Hello World!  
>Hello World!  
>Hello World!

The code prints `Hello World!` four times, but obviously this isn't because of any `for` loop.
Each of the four OpenMP threads that we specified with `OMP_NUM_THREADS` executes the `printf` call independently.

We will now edit the code so that each time `Hello World!` is printed we are told which of the threads is responsible.
Open `example1.cpp` and edit code to read:

``` cpp
#include <stdio.h>
#include <omp.h>

int main()
{
#pragma omp parallel
  {
    int thread_id = omp_get_thread_num();
    printf("Hello World! (%i)\n", thread_id);
  }
  
  return 0
}
```

Compile and run the code:

**$ cmake .  
$ make  
$ ./run.sh**

>Hello World! (2)  
>Hello World! (3)  
>Hello World! (1)  
>Hello World! (0)

Run the code a few more times, and you will find that the order in which the threads print their message isn't always the same.
This is because all of the threads run simultaneously, and depending on largely random factors some will finish slightly before others.
The order in which the threads finish determines the order in which the messages are printed.

Let's play around with these threads some more.
We will add another `printf` after the OpenMP-parallelized block, followed by a third `printf` inside a new OpenMP-parallelized block.

``` cpp
#include <stdio.h>
#include <omp.h>

int main()
{
#pragma omp parallel
  {
    int thread_id = omp_get_thread_num();
    printf("Hello World! (%i)\n", thread_id);
  }

  int thread_id = omp_get_thread_num();
  printf("This is another message! (%i)\n", thread_id);

#pragma omp parallel num_threads(3)
  {
    int thread_id = omp_get_thread_num();
    printf("Goodbye World! (%i)\n", thread_id);
  }
  
  return 0;
}
```

If you compile and run this code, you should get something like

>Hello World! (2)  
>Hello World! (0)  
>Hello World! (1)  
>Hello World! (3)  
>This is another message! (0)  
>Goodbye World! (0)  
>Goodbye World! (1)  
>Goodbye World! (2)

The following happens when the code is run:

 1. The code begins executing sequentially, with only the "master thread" running.
 2. When the master thread encounters an OpenMP-parallelized block, it `forks`, causing a number of threads equal to OMP_NUM_THREADS to begin executing.
 3. At the end of the OpenMP-parallelized block, the threads are all `joined` to the master thread, and sequential execution is resumed.
 4. The second `printf` is executed by only the master thread.
 5. The master thread again forks, but only creates a total of three threads. 
 6. All of the threads print a new line, and then are joined again.

https://en.wikipedia.org/wiki/Fork%E2%80%93join_model#/media/File:Fork_join.svg

## Example 2

We will now work on Example 2.

**$ cd ../example2**

In this directory is a code that does some simple math on some arrays, and then prints the average of the result:

**$ cmake.  
$ make  
$ ./run.sh**

>Average: 500000001.500000

We would like to improve the speed of this calculation through parallelization.
Before doing any sort of optimizations on a code, we should learn about which parts of the calculation are the most expensive - otherwise, we might waste our time parallelizing a part of the code that already doesn't take long to run, while ignoring much more costly parts of the code.
The most basic information we need is the time it takes to run each different part of the code.
To learn this, we will add some extra code to time each calculation.

Open `example2.cpp` using a text editor, and add the following to the beginning of `main`:

``` cpp
  double start_time = omp_get_wtime();
```

Now add the following to the end of main:

``` cpp
  printf("Total time: %f\n",omp_get_wtime()-start_time);
```

The code should now output the total about of time required to run:

>Average: 500000001.500000  
>Total time: 6.816225

Add similar timing statements around each of the `for` loops.
This should produce output similar to the following:

>Initialize a time: 1.402934  
>Initialize b time: 4.077723  
>Add arrays time: 0.900899  
>Average result time: 0.381398  
>  
>Average: 500000001.500000  
>Total time: 6.781963

Let's start by parallelizing the loop that adds arrays `a` and `b`.

``` cpp
  for (int i=0; i<N; i++) {
    a[i] = a[i] + b[i];
  }
```

Parallelize this loop in the following way:

``` cpp
#pragma omp parallel
{
  int i;
  int id;
  int nthreads;
  int istart;
  int iend;
  int Nthr;

  id = omp_get_thread_num();
  nthreads = omp_get_num_threads();

  Nthr = N / nthreads;
  istart = id * Nthr;
  iend = (id+1) * Nthr;
  if (id == nthreads-1) iend = N;

  for (i=istart; i<iend; i++) {
    a[i] = a[i] + b[i];
  }
}
```

>Add arrays time: 0.348737

A much easier way to parallelize this loop is to say:

``` cpp
#pragma	omp parallel for
  for (int i=0; i<N; i++) {
    a[i] = a[i] + b[i];
  }
```

>Add arrays time: 0.328928

Now we will parallelize the loop that averages the result.

``` cpp
#pragma omp parallel for
  for (int i=0; i<N; i++) {
    average += a[i];
  }
  average = average/double(N);
```

>Initialize a time: 1.397873  
>Initialize b time: 4.076062  
>Add arrays time: 0.342933  
>Average result time: 0.123558  
>  
>Average: 156250000.375000  
>Total time: 5.960761

This made the loop finish more quickly, but it also changed the final result.
What happened?
A strong hint comes from the fact that this result is the same thing we would get if we summed only the first 1/4 of the elements in array `a`.
The problem is that each thread works on its own copy of the `average` variable.
We want the loop to calculate the total value of `average`, summed across all threads.
This is known as a reduction operation, and we can tell the compiler that this is what we want by adding a `reduction` clause to the `pragma`:

``` cpp
#pragma omp parallel for reduction(+:average)
  for (int i=0; i<N; i++) {
    average += a[i];
  }
  average = average/double(N);
```

>Initialize a time: 1.383079  
>Initialize b time: 4.074402  
>Add arrays time: 0.346857  
>Average result time: 0.123527  
>  
>Average: 500000001.500000  
>Total time: 5.948352

We've made some nice improvements to a couple of the loops, but the real bottlenect happens when `a` and `b` are initialized.
We will now work on the loop that initializes `a`.
Start by adding an OpenMP pragma:

``` cpp
#pragma	omp parallel for
  for (int i=0; i<N; i++) {
    a[i] = 1.0;
  }
```

>Initialize a time: 0.391194  
>Initialize b time: 4.439494  
>Add arrays time: 0.296057  
>Average result time: 0.105363  
>  
>Average: 500000001.500000  
>Total time: 5.252643

Now do the same thing to b:

``` cpp
#pragma	omp parallel for
  for (int i=0; i<N; i++) {
    b[i] = 1.0 + double(i);
  }
```

>Initialize a time: 0.344843  
>Initialize b time: 0.350311  
>Add arrays time: 0.188084  
>Average result time: 0.085440  
>  
>Average: 500000001.500000  
>Total time: 0.987741

Notice that now even the last two loops are faster.
This is because of the principle of `first touch`.

## Example 3

Now we will look at Example 3.

**$ cd ../example3**

Compiling and running this code should produce something like the following:

**$ cmake .  
$ make  
$ ./run.sh**

>Iteration: 999      Energy: 92079.129718      PE: 16253.127101  
>  
>Timings:  
>   Force Zero:      0.000707  
>   Force Calc:      15.992926  
>   Velocity Update: 0.014128  
>   Coords Update:   0.001717  
>   Total:           16.014781

Clearly, the most expensive part of the code is the region where the forces are calculated, which consists of a double loop over all particles.

``` cpp
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
```

One approach to speeding this code up would be to add OpenMP parallelization to the inner loop, like this:

``` cpp
  for (int i=0; i < natoms; i++) {
#pragma omp parallel for reduction(+:v)
    for (int j=i+1; j < natoms; j++) {
       ...
    }
  }
```

>Iteration: 999      Energy: 286004.735108      PE: 18173.082058  
>  
>Timings:  
>   Force Zero:      0.001293  
>   Force Calc:      14.219656  
>   Velocity Update: 0.014499  
>   Coords Update:   0.005649  
>   Total:           14.246469

Unfortunately, this change causes the simulation to no longer produce the same results!

To understand why the results are different, look closely at how this loop updates the forces.
In the innermost loop is the line `forces[i][0] -= fx;`, which tells the computer "Copy the value of `forces[j][0]` from memory and store it in cache, then subtract `fx` from it, then replace the value of `forces[i][0]` in memory with the result.
When we run with multiple threads, different threads may be attempting to update 'forces[i][0]' at the same time.
If one thread copies the value of `forces[i][0]` into cache, then another thread updates it while the first thread is doing the subtration operation, the outcome will be that the first thread overwrites the update.
This scenario is known as a `race condition`, and happens any time two threads try to edit the same data at the same time.
If you run this calculation several times, you will get slightly different results each time depending on the order in which the threads "race" to update the forces.

One way to eliminate the race condition is to switch from looping over "unique pairs" to looping over "all pairs."
This way, we don't need to update `forces[i][0]` at all, only `forces[j][0]`.
The modified code looks like:

``` cpp
  double v = 0.0;
  for (int i=0; i < natoms; i++) {
#pragma omp parallel for reduction(+:v)
    for (int j=0; j < natoms; j++) {
      if (i != j ) {
        double dx = coords[j][0] - coords[i][0];
        double dy = coords[j][1] - coords[i][1];
        double dr2 = dx*dx + dy*dy;
        double dr = sqrt(dr2);
        double fx = (dx/dr)*(1.0/dr2);
        double fy = (dy/dr)*(1.0/dr2);
	//forces[i][0] -= fx;
	//forces[i][1] -= fy;
        forces[j][0] += fx;
        forces[j][1] += fy;
        v += 0.5/dr;
      } 
    }
  }
  *potential = v;
```

The main differences here are: (1) we changed the starting value of `j` from `i+1` to `0`, (2) we commented out the updates to `forces[i]`, (3) we multiple all contributions to `v` by `0.5` to avoid a double-counting error, and (4) we add an `if` check to avoid the `i == j` case.
This code produces the correct result while also being substantially faster.
   
>Iteration: 999      Energy: 92079.129718      PE: 16253.127101  
>  
>Timings:  
>   Force Zero:      0.000792  
>   Force Calc:      7.954821  
>   Velocity Update: 0.003412  
>   Coords Update:   0.005082  
>   Total:           7.969381  

That having been said, we can probably do somewhat better.
Every time the code enters an OpenMP-threaded region, the master thread must perform a `fork` (and later a `join`).
Because of this, there is an overhead cost every time the code encounters an OpenMP-parallelized region; the exact amount of overhead can vary, but it is usually on the order of 1 microsecond.
Think for a moment about how much time the code spends each time it enters our OpenMP region.
The total amount of time our code spends calculating the forces is 8.0 seconds, and the code enters the OpenMP region a number of times equal to the number of atoms times the number of iterations.
This works out to approximately 5 microseconds spent inside the OpenMP region each time it is called.
That isn't ideal - it is likely that a decent amount of that time is just `fork`/`join` overhead.

We can improve things by moving the parallelization to the outer loop over `i`.

``` cpp
#pragma omp parallel for reduction(+:v)
  for (int i=0; i < natoms; i++) {
    for (int j=0; j < natoms; j++) {
        double dx = coords[j][0] - coords[i][0];
        double dy = coords[j][1] - coords[i][1];
        double dr2 = dx*dx + dy*dy;
        double dr = sqrt(dr2);
        double fx = (dx/dr)*(1.0/dr2);
        double fy = (dy/dr)*(1.0/dr2);
        forces[i][0] -= fx;
        forces[i][1] -= fy;
	//forces[j][0] += fx;
	//forces[j][1] += fy;
        v += 0.5/dr;
    }
  }
```

>Iteration: 999      Energy: 92079.129718      PE: 16253.127101  
>  
>Timings:  
>   Force Zero:      0.000737  
>   Force Calc:      4.770037  
>   Velocity Update: 0.002868  
>   Coords Update:   0.005819  
>   Total:           4.784812  

That definitely helped - apparently the code was spending nearly half the time doing OpenMP overhead work.

## Example 4

Now we will work on a somewhat more realistic example of an MD code.

**$ cd ../example4  
$ cmake .  
$ make  
$ ./run.sh**

>Time step 0.03  
>optim: 1.84389e+17 0.1  
>optim: 11479.9 0.1  
>optim: -10384.3 0.1  
>optim: -14405.6 0.1  
>optim: -16274.2 0.1  
>optim: -17473.3 0.1  
>optim: -18386.3 0.1  
>optim: -19046.1 0.1  
>optim: -19516.6 0.1  
>optim: -19881.7 0.1  
>optim: -20184.1 0.1  
>optim: -20445.1 0.1  
>  
>    time         ke            pe             e            T          P  
>  -------    -----------   ------------  ------------    ------    ------  
>     3.00     2937.72894   -18197.02586  -15259.29692    0.371   0.04086690    
>     6.00     3426.85950   -18452.47277  -15025.61327    0.408   0.04039200    
>     9.00     3603.35315   -18693.26130  -15089.90815    0.441   0.03943282    
>    12.00     3615.37003   -19040.76407  -15425.39404    0.436   0.03848593    
>    15.00     3686.47202   -19111.90579  -15425.43377    0.453   0.03809054    
>    18.00     3702.65533   -19128.06154  -15425.40621    0.465   0.03778621    
>    21.00     3812.43640   -19237.88713  -15425.45072    0.472   0.03758076    
>    24.00     3849.80094   -19275.24791  -15425.44697    0.480   0.03741216    
>    27.00     3962.65623   -19388.13553  -15425.47930    0.492   0.03713681    
>    30.00     3973.40391   -19398.90927  -15425.50536    0.495   0.03703538    
>times:  force=14.85s  neigh=17.29s  total=32.24s

The code primarily consists of three sections: (1) The `neighbor_list` function, which generates a list of atom pairs that are close enough to be considered interacting, (2) a `forces` function, which calculates all contributions to the forces from the pairs in the neighbor list, and (3) the `md` function, which runs the calculation by calling `neighbor_list` and `forces` and then updating the atomic coordinates each timestep.


Note that the force evaluation is handled via an iterator over the neighbor list:

``` cpp
    for (neighT::const_iterator ij=neigh.begin(); ij!=neigh.end(); ++ij) {
       ...
    }
```

Loops with iterators are a little awkward to parallelize with OpenMP.
The approach we will take is to create a separate neighborlist for each thread.
First, define a `thrneighT` type, which will be a vector of neighborlists (add this near the beginning of `md.cc`, with the other typdefs).

``` cpp
typedef std::vector<neighT> thrneighT;
```

Then create a `thrneighT` which will contain all of the neighborlists for all of the threads.
Where appropriate, change `neighT` and `neigh` to `thrneighT` and `thrneigh`, respectively.
Here are the locations where these changes need to be made:

``` cpp
void neighbor_list(const coordT& coords, thrneighT& thrneigh) {
    double start = omp_get_wtime();
    for (int ithr=0; ithr<thrneigh.size(); ithr++) {
        neighT& neigh = thrneigh[ithr];
        neigh.clear();
        for (int i=ithr; i<natom; i+=thrneigh.size()) {
	  ...
	}
        ...
    }
...
coordT forces(const thrneighT& thrneigh, const coordT& coords, double& virial, double& pe) {
...
    for (int ithr=0; ithr<thrneigh.size(); ithr++) {
      const neighT& neigh = thrneigh[ithr];
      for (neighT::const_iterator ij=neigh.begin(); ij!=neigh.end(); ++ij) {
...
void optimize(coordT& coords, thrneighT& thrneigh) {
    double dt = 0.1;
    double prev = 1e99;
    for (int step=0; step<600; step++) {
      if ((step%(3*nneigh)) == 0 || step<10) neighbor_list(coords, thrneigh);
        double virial,pe;
        coordT f = forces(thrneigh,coords,virial,pe);
...
    // Relax the initial random guess
    int	nthreads;
#pragma	omp parallel
    {
      nthreads = omp_get_num_threads();
    }
    thrneighT thrneigh(nthreads);
    optimize(coords, thrneigh);
    neighbor_list(coords, thrneigh);
...
    coordT f = forces(thrneigh,coords,virial,potential_energy);
...
        // make the forces at time t+dt
        if ((step%nneigh) == 0) {
          neighbor_list(coords, thrneigh);
        }

        double virial_step;
        f = forces(thrneigh,coords,virial_step,potential_energy);
```

>    30.00     3973.40320   -19398.90856  -15425.50536    0.495   0.03703538    
>times:  force=14.50s  neigh=18.06s  total=32.67s

We are finally ready to work on parallelization of the `forces` subroutine.
Replacing the outer `for` loop in `forces` is easy enough:

``` cpp
#pragma	omp parallel default(none) shared(f, thrneigh, coords, virial, pe)
    {
       int ithr = omp_get_thread_num();
       const neighT& neigh = thrneigh[ithr];
       ...
    }
```

Unfortunately, there is one significant problem: the forces are incremented in a way that introduces a race condition.

``` cpp
          f[i].first += dfx;
          f[j].first -= dfx;
          f[i].second += dfy;
          f[j].second -= dfy;
```

Unlike the case with `Example 3`, it isn't straightforward to restructure the force evaluation in such a way that the forces on each atom are calculated by only one thread.
Instead, we will have each thread calculate some portion of the forces on all of the atoms, then perform a `reduction` operation on the forces.
The OpenMP `reduction` clause only works for for basic data types; it won't work for our defined type `coordT`.
As a result, we must perform the reduction manually.

To reduce the forces, first declare an array that will store the forces for a single thread.

``` cpp
#pragma	omp parallel default(none) shared(f, thrneigh, coords, virial, pe)
    {
      coordT f_thread(natom,xyT(0.0,0.0));
      double virial_thread = 0.0;
      double pe_thread = 0.0;
      ...
          f_thread[i].first += dfx;
          f_thread[j].first -= dfx;
          f_thread[i].second += dfy;
          f_thread[j].second -= dfy;

          pe_thread += vij;
          virial_thread += dfx*dx + dfy*dy;
        ...
	}
      }
    }
```

Then, `just before` the end of the OpenMP block, write the following:

``` cpp
#pragma omp critical
      {
        for (int i=0; i<natom; i++) {
          f[i].first += f_thread[i].first;
          f[i].second += f_thread[i].second;
        }
	pe += pe_thread;
	virial += virial_thread;
      }
```

Running on four threads, this gives improved performance for the forces:

>    30.00     3973.40314   -19398.90851  -15425.50536    0.495   0.03703538    
>times:  force=7.03s  neigh=18.99s  total=26.23s

Now let's work on the `neighbor_list` function.
Once again, we will convert the `for` loop over `thrneigh` into an OpenMP parallelized region.

``` cpp
#pragma omp parallel default(none) shared(thrneigh, coords)
    {
        int ithr = omp_get_thread_num();
	int nthread = omp_get_num_threads();
	...
    }
```

This improves the neighborlist creation time substantially:

>times:  force=3.67s  neigh=6.65s  total=10.49s

Running on 16 threads, we get this:

>times:  force=1.58s  neigh=2.26s  total=4.03s

One issue that limits our performance at high thread counts is the fact that we are not doing anything to `load balance` the neighborlist across threads.
We can improve the load balancing by adding some code to ensure that each thread has a neighborlist of similar length.
Add the following just before the end of the OpenMP parallelized block in `neighbor_list`:

``` cpp
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
```

Unfortunately, this really doesn't help our timings:

>times:  force=1.59s  neigh=2.37s  total=4.15s

The problem is that the added code does lots of operations that resize the neighborlists, which is a slow operation on a C++ `list` type.
To fix this, change the `neighT` type to be a vector:

````{tab-set-code}

```{code-block} c++
typedef std::vector<pairT> neighT;
```
````

Then, add the following to the `neighbor_list` function, just after the `neigh.clear();` line:

````{tab-set-code}

```{code-block} c++
neigh.reserve(100*natom/nthread);
```
````

These changes allow us to get our best timings yet:

>times:  force=1.55s  neigh=1.54s  total=3.27s

````{admonition} Key Points
:class: key

- It is extremely important to carefully avoid race conditions.
- Achieving good performance with OpenMP often requires larger code refactoring.
````
