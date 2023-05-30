---
html_meta:
  "description lang=en": "A detailed guide on how to use MPI in C++ for parallelizing code, handling errors, using non-blocking communication methods, and debugging parallelized code. Ideal for developers and programmers involved in high-performance computing. Three hands-on examples including parallelization of a Monte Carlo simulation code."
  "keywords": "Parallel Programming, MPI, Python, mpi4py, MPI-parallelized code, Message Passing Interface, MPI Standard, MPI operations, NumPy, MPI Timer, MPI Communicator, MPI Ranks, OpenMPI, MPICH, MS MPI, mpiexec, mpirun, monte carlo simulation"
  "property=og:locale": "en_US"
---

# MPI Hands-On - C++

````{admonition} Overview
:class: overview

Questions:
- How can I use MPI to parallelize a compiled code?

Objectives:
- Compile and run C++ codes that are parallelized using MPI.
- Use proper MPI error handling.
- Learn how to use non-blocking communication methods.
- Use a debugger with an parallelized code.
````


## 1. Example 1

### Writing Hello World

We'll start with the first example in [mpi/hello](https://github.com/MolSSI-Education/parallel-programming/tree/main/examples/mpi/hello), which is a simple Hello World code:

````{tab-set-code} 

```{code-block} cpp
#include <iostream>

int main(int argc, char **argv) {
    std::cout << "Hello World!" << std::endl;
    return 0;
}
```
````


Acquire a copy of the example files for this lesson, and then compile and run the example code:

````{tab-set-code} 

```{code-block} shell
$ git clone git@github.com:MolSSI-Education/parallel-programming.git
$ cd parallel-programming/examples/mpi/hello
$ mkdir build
$ cd build
$ cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_Fortran_COMPILER=mpifort ..
$ make
$ ./hello
```
````


````{tab-set-code} 

```{code-block} output
Hello World!
```
````


### Getting Started with MPI

Let's try running this code on multiple processes.
This is done using the `mpiexec` command.
Many environments also provide an `mpirun` command, which usually - but not always - works the same way.
Whenever possible, you should use `mpiexec` and not `mpirun`, in order to guarantee more consistent results.

> ## MPI - `mpiexec` vs `mpirun`
> MPI stands for **'message passing interface'** and is a message passing standard which is designed to work on a variety of parallel computing architectures. The [MPI standard](https://www.mpi-forum.org/docs/drafts/mpi-2018-draft-report.pdf) defines how syntax and semantics of a library of routines. There are a number of implementations of this standard including OpenMPI, MPICH, and MS MPI.
>
>  The primary difference between `mpiexec` and `mpirun` is that `mpiexec` is defined as part of the MPI standard, while `mpirun` is not.  Different implementations of MPI (i.e. OpenMPI, MPICH, MS MPI, etc.) are not guaranteed to implement `mpirun`, or might implement different options for `mpirun`.  Technically, the MPI standard doesn't actually require that MPI implementations implement `mpiexec` either, but the standard does at least describe guidelines for how `mpiexec` should work.  Because of this, `mpiexec` is generally the preferred command.
>

The general format for lanching a code on multiple processes is:

````{tab-set-code} 

```{code-block} shell
$ mpiexec -n <number_of_processes> <command_to_launch_code>
```
````


For example, to launch `hello` on 4 processes, do:

````{tab-set-code} 

```{code-block} shell
$ mpiexec -n 4 ./hello
```
````


````{tab-set-code} 

```{code-block} output
Hello World!
Hello World!
Hello World!
Hello World!
```
````


When you execute the above command, `mpiexec` launches 4 different instances of `./hello` simultaneously, which each print "Hello World!".

Typically, as long as you have at least 4 processors on the machine you are running on, each process will be launched on a different processor; however, certain environment variables and optional arguments to `mpiexec` can change this behavior.
Each process runs the code in `hello` independently of the others.

It might not be obvious yet, but the processes `mpiexec` launches aren't completely unaware of one another.
The `mpiexec` adds each of the processes to an MPI communicator, which enables each of the processes to send and receive information to one another via MPI.
The MPI communicator that spans all of the processes launched by `mpiexec` is called `MPI_COMM_WORLD`.

We can use the MPI library to get some information about the `MPI_COMM_WORLD` communicator and the processes within it.
Edit `hello.cpp` so that it reads as follows:

````{tab-set-code} 

```{code-block} cpp
#include <iostream>
#include <mpi.h>

int main(int argc, char **argv) {
  // Initialize MPI
  // This must always be called before any other MPI functions
  MPI_Init(&argc, &argv);

  // Get the number of processes in MPI_COMM_WORLD
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of this process in MPI_COMM_WORLD
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  // Print out information about MPI_COMM_WORLD
  std::cout << "World Size: " << world_size << "   Rank: " << my_rank << std::endl;

  // Finalize MPI
  // This must always be called after all other MPI functions
  MPI_Finalize();

  return 0;
}
```
````


Recompile the code:

````{tab-set-code} 

```{code-block} shell
$ make
```
````


In the above code we first include the MPI library header, `mpi.h`.
Then, we call `MPI_Init()`.
This function **must** be called before any other MPI functions, and is typically one of the first lines of an MPI-parallelized code.
Then, we call `MPI_Comm_size()` to get the number of processes in `MPI_COMM_WORLD`, which corresponds to the number of processes launched whenever `mpiexec` is executed at the command line.
Each of these processes is assigned a uniqe rank, which is an integer that ranges from `0` to `world_size - 1`.
The rank of a process allows it to be identified whenever processes communicate with one another.
For example, in some cases we might want rank 2 to send some information to rank 4, or we might want rank 0 to receive information from all of the other processes.
Calling `MPI_Comm_rank()` returns the rank of the process calling it within `MPI_COMM_WORLD`.

Go ahead and run the code now:

````{tab-set-code} 

```{code-block} shell
$ mpiexec -n 4 ./hello
```
````


````{tab-set-code} 

```{code-block} output
World Size: 4   Rank: 1
World Size: 4   Rank: 0
World Size: 4   Rank: 2
World Size: 4   Rank: 3
```
````


As you can see, the `MPI_Comm_size()` function outputs 4, which is the total number of ranks we told `mpiexec` to run with (through the `-n` argument).
Each of the processes is assigned a rank in the range of 0 to 3.

As you can see, the ranks don't necessarily print out their messages in order; whichever rank completes the `cout` first will print out its message first.
If you run the code again, the ranks are likely to print thier message in a different order:

````{tab-set-code} 

```{code-block} output
World Size: 4   Rank: 2
World Size: 4   Rank: 0
World Size: 4   Rank: 3
World Size: 4   Rank: 1
```
````


You can also try rerunning with a different value for the `-n` `mpiexec` argument.
For example:

````{tab-set-code} 

```{code-block} shell
$ mpiexec -n 2 ./hello
```
````


````{tab-set-code} 

```{code-block} output
World Size: 2   Rank: 0
World Size: 2   Rank: 1
```
````




### Error Handling with MPI

If an error forces an MPI program to exit, it should **never** just call `return` or `exit`.
This is because calling `return` or `exit` might terminate one of the MPI processes, but leave others running (but not doing anything productive) indefintely.
If you're not careful, you could waste massive amounts of computational resources running a failed calculation your thought had terminated.
Instead, MPI-parallelized codes should call `MPI_Abort()` when something goes wrong, as this function will ensure all MPI processes terminate.
The `MPI_Abort` function takes two arguments: the first is the communicator corresponding to the set of MPI processes to terminate (this should generally be `MPI_COMM_WORLD`), while the second is an error code that will be returned to the environment.

It is also useful to keep in mind that most MPI functions have a return value that indicates whether the function completed succesfully.
If this value is equal to `MPI_SUCCESS`, the function was executed successfully; otherwise, the function call failed.
By default, MPI functions automatically abort if they encounter an error, so you'll only ever get a return value of `MPI_SUCCESS`.
If you want to handle MPI errors yourself, you can call `MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN)`; if you do this, you must check the return value of every MPI function and call `MPI_Abort` if it is not equal to `MPI_SUCCESS`.
For example, when initializing MPI, you might do the following:

````{tab-set-code} 

```{code-block} cpp
if (MPI_Init(&argc,&argv) != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, 1);
```
````





## Example 2

### Basic Infrastructure

We will now do some work with the the example in [examples/mpi/average](https://github.com/MolSSI-Education/parallel-programming/tree/main/examples/mpi/average), which does some simple math.
Run the code now.

````{tab-set-code} 

```{code-block} shell
$ cd parallel-programming/examples/mpi/average
$ mkdir build
$ cd build
$ cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_Fortran_COMPILER=mpifort ..
$ make
$ ./average
```
````


````{tab-set-code} 

```{code-block} output
Average: 100000001.5
```
````


Let's learn something about which parts of this code account for most of the run time.
MPI provides a timer, `MPI_Wtime()`, which returns the current walltime.
We can use this function to determine how long each section of the code takes to run.

For example, to determine how much time is spent initializing array `a`, do the following:

````{tab-set-code} 

```{code-block} cpp
  // Initialize a
  double start_time = MPI_Wtime();
  double *a = new double[N];
  for (int i=0; i<N; i++) {
    a[i] = 1.0;
  }
  double end_time = MPI_Wtime();
  if (my_rank == 0 ) {
    std::cout << "Initialize a time: " << end_time - start_time << std::endl;
  }
```
````


As the above code indicates, we don't really want every rank to print the timings, since that could look messy in the output.
Instead, we have only rank 0 print this information.
Of course, this requires that we add a few lines near the top of the code to initialize MPI and query the rank of each process:

````{tab-set-code} 

```{code-block} cpp
  // Initialize MPI
  int world_size, my_rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
```
````


Also determine and print the timings of each of the other sections of the code: the intialization of array `b`, the addition of the two arrays, and the final averaging of the result.
Your code should look something like this:

````{tab-set-code} 

```{code-block} cpp
#include <iostream>
#include <mpi.h>

int main(int argc, char **argv) {
  // Initialize MPI
  int world_size, my_rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  int N = 200000000;

  // Initialize a
  double start_time = MPI_Wtime();
  double *a = new double[N];
  for (int i=0; i<N; i++) {
    a[i] = 1.0;
  }
  double end_time = MPI_Wtime();
  if (my_rank == 0 ) {
    std::cout << "Initialize a time: " << end_time - start_time << std::endl;
  }

  // Initialize b
  start_time = MPI_Wtime();
  double *b = new double[N];
  for (int i=0; i<N; i++) {
    b[i] = 1.0 + double(i);
  }
  end_time = MPI_Wtime();
  if (my_rank == 0 ) {
    std::cout << "Initialize b time: " << end_time - start_time << std::endl;
  }
  
  // Add the two arrays
  start_time = MPI_Wtime();
  for (int i=0; i<N; i++) {
    a[i] = a[i] + b[i];
  }
  end_time = MPI_Wtime();
  if (my_rank == 0 ) {
    std::cout << "Add arrays time: " << end_time - start_time << std::endl;
  }

  // Average the result
  start_time = MPI_Wtime();
  double average = 0.0;
  for (int i=0; i<N; i++) {
    average += a[i] / double(N);
  }
  end_time = MPI_Wtime();
  if (my_rank == 0 ) {
    std::cout << "Average result time: " << end_time - start_time << std::endl;
  }

  std::cout.precision(12);
  if (my_rank == 0 ) {
    std::cout << "Average: " << average << std::endl;
  }
  delete [] a;
  delete [] b;
  MPI_Finalize();
  return 0;
}
```
````


Now compile and run the code again:

````{tab-set-code} 

```{code-block} shell
$ make
$ ./average
```
````


````{tab-set-code} 

```{code-block} output
Initialize a time: 0.544075
Initialize b time: 0.624939
Add arrays time: 0.258915
Average result time: 0.266418
Average: 100000001.5
```
````



### Point-to-Point Communication

You can try running this on multiple ranks now:

````{tab-set-code} 

```{code-block} shell
$ mpiexec -n 4 ./average
```
````


````{tab-set-code} 

```{code-block} output
Initialize a time: 0.640894
Initialize b time: 0.893775
Add arrays time: 1.38309
Average result time: 0.330192
Average: 100000001.5
```
````


Running on multiple ranks doesn't help with the timings, because each rank is duplicating all of the same work.
In some ways, running on multiple ranks makes the timings worse, because all of the processes are forced to compete for the same computational resources.
Memory bandwidth in particular is likely a serious problem due to the extremely large arrays that must be accessed and manipulated by each process.
We want the ranks to cooperate on the problem, with each rank working on a different part of the calculation.
In this example, that means that different ranks will work on different parts of the arrays `a` and `b`, and then the results on each rank will be summed across all the ranks.

In this section, we will handle the details of the communication between processes using *point-to-point* communication.
Point-to-point communication involves cases in which a code explicitly instructs one specific process to send/recieve information to/from another specific process.
The primary functions associated with this approach are `MPI_Send()` and `MPI_Recv()`, which are involve the following arguments:

````{tab-set-code} 

```{code-block} cpp
int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,  MPI_Comm comm);
```
````


* `buf` --- pointer to the start of the buffer being sent
* `count` --- number of elements to send
* `datatype` --- MPI data type of each element
* `dest` --- rank of destination process
* `tag`  --- message tag
* `comm` --- the communicator to use

````{tab-set-code} 

```{code-block} cpp
int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
```
````


* `buf` --- pointer to the start of the buffer to receive the message
* `count` --- maximum number of elements the buffer can hold
* `datatype` --- MPI data type of each element
* `source` --- rank of source process ---  `MPI_ANY_SOURCE` matches any process
* `tag`  --- message tag (integer `>= 0`) --- `MPI_ANY_TAG` matches any tag
* `comm` --- the communicator to use
* `status` --- pointer to the structure in which to store status

We need to decide what parts of the arrays each of the ranks will work on; this is more generally known as a rank's workload.
Add the following code just before the initialization of array `a`:

````{tab-set-code} 

```{code-block} cpp
  // Determine the workload of each ran
  int workloads[world_size];
  for (int i=0; i<world_size; i++) {
    workloads[i] = N / world_size;
    if ( i < N % world_size ) workloads[i]++;
  }
  int my_start = 0;
  for (int i=0; i<my_rank; i++) {
    my_start +=	workloads[i];
  }
  int my_end = my_start	+ workloads[my_rank];
```
````


In the above code, `my_start` and `my_end` represent the range over which each rank will perform mathematical operations on the arrays.

We'll start by parallelizing the code that averages the result.
Update the range of the `for` loop in this part of the code to the following:

````{tab-set-code} 

```{code-block} cpp
  for (int i=my_start; i<my_end; i++) {
```
````


This will ensure that each rank is only calculating elements `my_start` through `my_end` of the sum.
We then need the ranks to communicate their individually calculated sums so that we can calculate the global sum.
To do this, add the following immediately after the end of the `for` loop:

````{tab-set-code} 

```{code-block} cpp
  if ( my_rank == 0 ) {
    for (int i=1; i<world_size; i++) {
      double partial_average;
      MPI_Status status;
      MPI_Recv( &partial_average, 1, MPI_DOUBLE, i, 77, MPI_COMM_WORLD, &status );
      average += partial_average;
    }
  }
  else {
    MPI_Send( &average, 1, MPI_DOUBLE, 0, 77, MPI_COMM_WORLD );
  }
```
````


The `MPI_DOUBLE` parameter tells MPI what type of information is being communicated by the `Send` and `Recv` calls.
In this case, we are sending a array of double precision numbers.
If you are communicating information of a different datatype, consult the following:

|**MPI data type**  |**C data type**     |
|:------------------|:-------------------|
|`MPI_BYTE`           |8 binary digits     |
|`MPI_CHAR`           |char                |
|`MPI_UNSIGNED_CHAR`  |unsigned char       |
|`MPI_SHORT`          |signed short int	 |	 
|`MPI_UNSIGNED_SHORT` |unsigned short int	 |	 
|`MPI_INT`            |signed int          |
|`MPI_UNSIGNED`       |unsigned int	 |	 
|`MPI_LONG`           |signed long int	 |	 
|`MPI_UNSIGNED_LONG`  |unsigned long int	 |	 
|`MPI_FLOAT`          |float               |
|`MPI_DOUBLE`         |double              |
|etc.               |                    |
|`MPI_PACKED`	    |define your own with|
|                   |[`MPI_Pack`](https://www.mpich.org/static/docs/v3.2/www3/MPI_Pack.html)/[`MPI_Unpack`](https://www.mpich.org/static/docs/v3.2/www3/MPI_Pack.html) |

Now compile and run the code again:

````{tab-set-code} 

```{code-block} shell
$ make
$ mpiexec -n 4 ./average
```
````


````{tab-set-code} 

```{code-block} output
Initialize a time: 0.63251
Initialize b time: 1.31379
Add arrays time: 1.89099
Average result time: 0.100575
Average: 100000001.5
```
````


You can see that the amount of time spent calculating the average has indeed gone down.

Parallelizing the part of the code that adds the two arrays is much easier.
All you need to do is update the range over which the `for` loop iterates:

````{tab-set-code} 

```{code-block} cpp
  for (int i=my_start; i<my_end; i++) {
```
````


Now compile and run the code again:

````{tab-set-code} 

```{code-block} shell
$ make
$ mpiexec -n 4 ./average
```
````


````{tab-set-code} 

```{code-block} output
Initialize a time: 0.636685
Initialize b time: 1.66542
Add arrays time: 0.466888
Average result time: 0.0871116
Average: 100000001.5
```
````


The array addition time has gone down nicely.
Surprisingly enough, the most expensive part of the calculation is now the initialization of the arrays `a` and `b`.
Updating the range over which those loops iterate speeds up those parts of the calation:

````{tab-set-code} 

```{code-block} cpp
  // Initialze a
  for (int i=my_start; i<my_end; i++) {
...
  // Initialize b
  for (int i=my_start; i<my_end; i++) {
```
````


````{tab-set-code} 

```{code-block} shell
$ make
$ ./average
```
````


````{tab-set-code} 

```{code-block} output
Initialize a time: 0.159471
Initialize b time: 0.183946
Add arrays time: 0.193497
Average result time: 0.0847806
Average: 100000001.5
```
````




### Reducing the Memory Footprint

The simulation is running much faster now thanks to the parallelization we have added.
If that's all we care about, we could stop working on the code now.
In reality, though, **time** is only one resource we should be concerned about.
Another resource that is often even more important is **memory**.
The changes we have made to the code make it run faster, but don't decrease its memory footprint in any way: each rank allocates arrays `a` and `b` with `N` double precision values.
That means that each rank allocates `2*N` double precision values; across all of our ranks, that corresponds to a total of `2*nproc*world_size` double precision values.
Running on more processors might decrease our run time, but it increases our memory footprint!

Of course, there isn't really a good reason for each rank to allocate the entire arrays of size `N`, because each rank will only ever use values within the range of `my_start` to `my_end`.
Let's modify the code so that each rank allocates `a` and `b` to a size of `workloads[my_rank]`.

Replace the initialization of `a` with:

````{tab-set-code} 

```{code-block} cpp
  double *a = new double[ workloads[my_rank] ];
  for (int i=0; i<workloads[my_rank]; i++) {
    a[i] = 1.0;
  }
```
````


Replace the initialization of `b` with:

````{tab-set-code} 

```{code-block} cpp
  double *b = new double[ workloads[my_rank] ];
  for (int i=0; i<workloads[my_rank]; i++) {
    b[i] = 1.0 + double(i + my_start);
  }
```
````


Replace the range of the loops that add and average the arrays to `for (int i=0; i<workloads[my_rank]; i++)`.

Now compile and run the code again:

````{tab-set-code} 

```{code-block} shell
$ make
$ ./average
```
````


````{tab-set-code} 

```{code-block} output
Initialize a time: 0.16013
Initialize b time: 0.176896
Add arrays time: 0.190774
Average result time: 0.0871552
Average: 100000001.5
```
````



### Collective Communication

Previously, we used **point-to-point communication** (i.e. `MPI_Send` and `MPI_Recv`) to sum the results across all ranks:

````{tab-set-code} 

```{code-block} cpp
  if ( my_rank == 0 ) {
    for (int i=1; i<world_size; i++) {
      double partial_average;
      MPI_Status status;
      MPI_Recv( &partial_average, 1, MPI_DOUBLE, i, 77, MPI_COMM_WORLD, &status );
      average += partial_average;
    }
  }
  else {
    MPI_Send( &average, 1, MPI_DOUBLE, 0, 77, MPI_COMM_WORLD );
  }
```
````


MPI provides many **collective communication** functions, which automate many processes that can be complicated to write out using only point-to-point communication.
One particularly useful collective communication function is `MPI_Reduce()`:

````{tab-set-code} 

```{code-block} cpp
  int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                 MPI_Op op, int root, MPI_Comm comm)
```
````


* `sendbuf` --- address of send buffer
* `recvbuf` --- address of receive buffer
* `count` --- number of elements in send buffer
* `datatype` --- MPI data type of each element
* `op` --- reduce operation
* `root` --- rank of root process
* `comm` --- the communicator to use

Possible values for `op` are:

|**Operation** | Description | Datatype|
|:-------------|:------------|:--------|
|`MPI_MAX`       |maximum      |integer,float|
|`MPI_MIN`       |minimum      |integer,float|
|`MPI_SUM`       |sum          |integer,float|
|`MPI_PROD`      |product      |integer,float|
|`MPI_LAND`      |logical AND  |integer|
|`MPI_BAND`      |bit-wise AND |integer,MPI_BYTE|
|`MPI_LOR`       |logical OR   |integer|
|`MPI_BOR`       |bit-wise OR  |integer,MPI_BYTE|
|`MPI_LXOR`      |logical XOR  |integer|
|`MPI_BXOR`      |bit-wise XOR |integer,MPI_BYTE|
|`MPI_MAXLOC`    |max value and location|float|
|`MPI_MINLOC`    |min value and location|float|

We will use the `MPI_Reduce()` function to sum a value across all ranks, without all of the point-to-point communication code we needed earlier.
Replace all of your point-to-point communication code above with:

````{tab-set-code} 

```{code-block} cpp
  double partial_average = average;
  MPI_Reduce(&partial_average, &average, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
```
````


Compiling and running with this change should produce the same results as before.

Note that in addition to enabling us to write simpler-looking code, collective communication operations tend to be faster than what we can achieve by trying to write our own communication operations using point-to-point calls.






## Example 3

Next, view [examples/mc](https://github.com/MolSSI-Education/parallel-programming/tree/main/examples/mpi/mc) which is a simple Monte-Carlo simulation.
Compile and run the code now.

````{tab-set-code} 

```{code-block} shell
$ cd parallel-programming/examples/mpi/mc
$ mkdir build
$ cd build
$ cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_Fortran_COMPILER=mpifort ..
$ make
$ ./mc
```
````



````{tab-set-code} 

```{code-block} output
...
497000   -6.28643
498000   -6.28989
499000   -5.96743
500000   -6.06861
Total simulation time: 2.59121
    Energy time:       2.47059
    Decision time:     0.0425119
```
````


As you can see, the code already has some timings, and the vast majority of time is spent in the calls to 'get_particle_energy()'.
That is where we will focus our parallelization efforts.

The function in question is:

````{tab-set-code} 

```{code-block} cpp
double get_particle_energy(double *coordinates, int particle_count, double box_length, int i_particle, double cutoff2) {
  double e_total = 0.0;
  double *i_position = &coordinates[3*i_particle];

  for (int j_particle=0; j_particle < particle_count; j_particle++) {
    if ( i_particle != j_particle ) {
      double *j_position = &coordinates[3*j_particle];
      double rij2 = minimum_image_distance( i_position, j_position, box_length );
      if ( rij2 < cutoff2 ) {
        e_total += lennard_jones_potential(rij2);
      }
    }
  }

  return e_total;
}
```
````


This looks like it should be fairly straightforward to parallelize: it consists of a single `for` loop which just sums the interaction energies of particle pairs.
To parallelize this loop, we need each rank to compute the interaction energies of a subset of these pairs, and then sum the energy across all ranks.

The `get_particle_energy()` function is going to need to know some basic information about the MPI communicator, so add the MPI communicator to its parameters:

````{tab-set-code} 

```{code-block} cpp
double get_particle_energy(double *coordinates, int particle_count, double box_length, int i_particle, double cutoff2, MPI_Comm comm) {
```
````


Now update the two times `get_particle_energy()` is called by `main`:

````{tab-set-code} 

```{code-block} cpp
    double current_energy = get_particle_energy( coordinates, num_particles, box_length, i_particle, simulation_cutoff2, world_comm );
    ...
    double proposed_energy = get_particle_energy( coordinates, num_particles, box_length, i_particle, simulation_cutoff2, world_comm );
```
````


Place the following at the beginning of `get_particle_energy()`:

````{tab-set-code} 

```{code-block} cpp
  // Get information about the MPI communicator
  int my_rank, world_size;
  MPI_Comm_size(comm, &world_size);
  MPI_Comm_rank(comm, &my_rank);
```
````


Change the `for` loop in `get_particle_energy` to the following:

````{tab-set-code} 

```{code-block} cpp
  for (int j_particle=my_rank; j_particle < particle_count; j_particle += world_size) {
```
````


The above code will cause each rank to iterate over particles with a stride of `world_size` and an initial offset of `my_rank`.
For example, if you run on 4 ranks, rank 0 will iterate over particles 0, 4, 8, 12, etc., while rank 1 will iterate over particles 1, 5, 9, 13, etc.

We then need to sum the energies across all ranks.
Replace the line `return e_total;` with the following:

````{tab-set-code} 

```{code-block} cpp
  // Sum the energy across all ranks
  double e_summed = 0.0;
  MPI_Reduce(&e_total, &e_summed, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  return e_summed;
```
````


Try to run it in parallel now:

````{tab-set-code} 

```{code-block} shell
$ make
$ mpiexec -n 4 ./mc
```
````


````{tab-set-code} 

```{code-block} output
16000   -2.19516e+19
17000   -2.19517e+19

===================================================================================
=   BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES
=   PID 81043 RUNNING AT Taylors-MacBook-Pro.local
=   EXIT CODE: 9
=   CLEANING UP REMAINING PROCESSES
=   YOU CAN IGNORE THE BELOW CLEANUP MESSAGES
===================================================================================
YOUR APPLICATION TERMINATED WITH THE EXIT STRING: Segmentation fault: 11 (signal 11)
This typically refers to a problem with your application.
Please see the FAQ page for debugging suggestions
```
````


That doesn't seem right at all.
What went wrong?

Our call to `MPI_Reduce` causes the energies to be summed onto rank 0, but none of the other ranks have the summed energies.
To have the energies reduced to all of the ranks, replace the `MPI_Reduce` call with a call to `MPI_Allreduce`:

````{tab-set-code} 

```{code-block} cpp
  MPI_Allreduce(&e_total, &e_summed, 1, MPI_DOUBLE, MPI_SUM, comm);
```
````


````{tab-set-code} 

```{code-block} shell
$ make
$ mpiexec -n 4 ./mc
```
````


````{tab-set-code} 

```{code-block} output
497000   -6.28644
498000   -6.2899
499000   -5.96744
500000   -6.06862
Total simulation time: 8.38658
    Energy time:       8.24201
    Decision time:     0.0563176
```
````


This is better, but we certainly aren't getting good timings.

Before we work on the timings problem, let's try another experiment.
Near the top of `mc.cpp` is some code that initializes the random number generator with a random seed.
Currently, the random number generator is being initialized with a fixed random seed of `1`:

````{tab-set-code} 

```{code-block} cpp
// Initialize the random number generator with a pre-defined seed
std::mt19937 mt(1);

// Initialize the random number generator with a random seed
//std::random_device rd;
//std::mt19937 mt(rd());
```
````


Let's try switching to using a random number to initialize the random number generator:

````{tab-set-code} 

```{code-block} cpp
// Initialize the random number generator with a pre-defined seed
//std::mt19937 mt(1);

// Initialize the random number generator with a random seed
std::random_device rd;
std::mt19937 mt(rd());
```
````


Now recompile and run again:

````{tab-set-code} 

```{code-block} shell
$ make
$ mpiexec -n 4 ./mc
```
````


````{tab-set-code} 

```{code-block} output
497000   -7.9279e+12
498000   -7.9279e+12
499000   -7.9279e+12
500000   -7.9279e+12
Total simulation time: 8.73895
    Energy time:       8.59179
    Decision time:     0.0549728
```
````


These are some crazy, unphysical numbers!
If you try running again with a single process (`mpiexec -n 1 ./mc`), you can confirm that the code gives much more reasonable energies when running on in serial.
The problem is that each iteration, the coordinates are updated by randomly displacing one of the particles.
Each rank randomly selects a particle to displace and the displacement vector.
Instead of contributing to the calculation of the same particle's interaction energies for the same nuclear configuration, each rank ends up calculating some of the interaction energies for different atoms and different coordinates.
This leads to utter chaos throughout the simulation.

You might be tempted to fix this by generating the random number generator seed on a single process, and then sending that information to the other processes, so that every process is using the same seed.
Although that might sound reasonable, it would still leave open the possibility that different processes could end up diverging over the course of a long simulation (remember, computers aren't infinitely accurate - slight descrepancies are guaranteed to happen, given enough time).

To fix this, we will have rank 0 be the only rank that randomly selects a particle or a displacement vector.
Rank 0 will then broadcast all necessary information to the other ranks, so that they keep in sync.

Replace the line where the coordinates are intially generated (where `generate_initial_state` is called) with this:

````{tab-set-code} 

```{code-block} cpp
  if ( my_rank == 0 ) {
    generate_initial_state(num_particles, box_length, coordinates);
  }
  MPI_Bcast(coordinates, 3*num_particles, MPI_DOUBLE, 0, world_comm);
```
````


At the beginning of the `for` loop in `main` you will see the following code:

````{tab-set-code} 

```{code-block} cpp
  // Beginning of main MC iterative loop
  n_trials = 0;
  for (int i_step=0; i_step<n_steps; i_step++) {
    n_trials += 1;
    int i_particle = floor( double(num_particles) * dist(mt) );
    double random_displacement[3];
    for (int i=0; i<3; i++) {
      random_displacement[i] = ( ( 2.0 * dist(mt) ) - 1.0 ) * max_displacement;
    }
```
````


Replace the above with the following:

````{tab-set-code} 

```{code-block} cpp
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
```
````


At the end of the `for` loop in `main` is the following code:

````{tab-set-code} 

```{code-block} cpp
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

      total_decision_time += MPI_Wtime() - start_decision_time;
```
````


Replace the above with the following, so that only rank `0` is executing it:

````{tab-set-code} 

```{code-block} cpp
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
```
````


Recompile and rerun the code.

````{tab-set-code} 

```{code-block} shell
$ make
$ mpiexec -n 4 ./mc
```
````


````{tab-set-code} 

```{code-block} output
497000   -6.06666
498000   -6.10058
499000   -5.98052
500000   -5.95301
Total simulation time: 16.7881
    Energy time:       15.9948
    Decision time:     0.0690625
```
````


This time the energies are much more consistent with what we expect; however, our timings are considerably worse than when we were only running on single process!
This is because the system we are running these calculations on is extremely small.
If you check the system parameters (under the `Parameter setup` comment in `mc.cpp`), you will see that this calculation only involves 100 particles.
The amount of work required to compute the energy of `100` Lennard-Jones particles is actually smaller than the amount of overhead associated with the extra MPI processes.
Let's make the simulation somewhat larger:

````{tab-set-code} 

```{code-block} cpp
  int n_steps = 100000;
  int freq = 1000;
  int num_particles = 10000;
```
````


Recompile and rerun the code on a single core.

````{tab-set-code} 

```{code-block} shell
$ make
$ ./mc
```
````


````{tab-set-code} 

```{code-block} output
97000   612.067
98000   609.113
99000   603.538
100000   599.461
Total simulation time: 41.0191
    Energy time:       39.9748
    Decision time:     0.011933
```
````


Finally, run the calculation in parallel.

````{tab-set-code} 

```{code-block} shell
$ mpiexec -n 4 ./mc
```
````


````{tab-set-code} 

```{code-block} output
97000   99.1126
98000   93.454
99000   91.1246
100000   87.397
Total simulation time: 22.2873
    Energy time:       15.4661
    Decision time:     0.0175401
```
````


Now we can clearly see a speedup when running in parallel.


````{admonition} Key Points
:class: key

- Where possible, use collective communication operations instead of point-to-point communication for improved efficiency and simplicity.
- Intelligent design choices can help you reduce the memory footprint required by MPI-parallelized codes
````
