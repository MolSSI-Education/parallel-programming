---
title: "MPI Hands-On - C++"
teaching: 0
exercises: 90
questions:
- "How can I use MPI to parallelize a compiled code?"
objectives:
- "Compile and run C++ codes that are parallelized using MPI."
- "Use proper MPI error handling."
- "Learn how to use non-blocking communication methods."
- "Use a debugger with an parallelized code."
keypoints:
- "Correct use of non-blocking can substantially reduce the time your code spends communicating between processes."
---

## 1. Example 1

### Writing Hello World

We'll start with the first example in [mpi/example1](https://github.com/MolSSI-Education/parallel-programming/tree/gh-pages/examples/mpi/example1), which is a simple Hello World code:

~~~
#include <iostream>

int main(int argc, char **argv) {
    std::cout << "Hello World!" << std::endl;
    return 0;
}
~~~
{: .language-cpp}

Acquire a copy of the example files for this lesson, and then compile and run the example code:

~~~
$ git clone git@github.com:MolSSI-Education/parallel-programming.git
$ cd parallel-programming/examples/mpi/hello
$ mkdir build
$ cd build
$ cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_Fortran_COMPILER=mpifort ..
$ make
$ ./hello
~~~
{: .language-bash}

~~~
Hello World!
~~~
{: .output}

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
{: .callout}

The general format for lanching a code on multiple processes is:

~~~
$ mpiexec -n <number_of_processes> <command_to_launch_code>
~~~
{: .language-bash}

For example, to launch `example1.py` on 4 processes, do:

~~~
$ mpiexec -n 4 ./hello
~~~
{: .language-bash}

~~~
Hello World!
Hello World!
Hello World!
Hello World!
~~~
{: .output}

When you execute the above command, `mpiexec` launches 4 different instances of `./hello` simultaneously, which each print "Hello World!".

Typically, as long as you have at least 4 processors on the machine you are running on, each process will be launched on a different processor; however, certain environment variables and optional arguments to `mpiexec` can change this behavior.
Each process runs the code in `example1.py` independently of the others.

It might not be obvious yet, but the processes `mpiexec` launches aren't completely unaware of one another.
The `mpiexec` adds each of the processes to an MPI communicator, which enables each of the processes to send and receive information to one another via MPI.
The MPI communicator that spans all of the processes launched by `mpiexec` is called `MPI_COMM_WORLD`.

We can use the MPI library to get some information about the `MPI_COMM_WORLD` communicator and the processes within it.
Edit `hello.cpp` so that it reads as follows:

~~~
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
~~~
{: .language-cpp}

Recompile the code:

~~~
$ make
~~~
{: .language-bash}

In the above code we first include the MPI library header, `mpi.h`.
Then, we call `MPI_Init()`.
This function **must** be called before any other MPI functions, and is typically one of the first lines of an MPI-parallelized code.
Then, we call `MPI_Comm_size()` to get the number of processes in `MPI_COMM_WORLD`, which corresponds to the number of processes launched whenever `mpiexec` is executed at the command line.
Each of these processes is assigned a uniqe rank, which is an integer that ranges from `0` to `world_size - 1`.
The rank of a process allows it to be identified whenever processes communicate with one another.
For example, in some cases we might want rank 2 to send some information to rank 4, or we might want rank 0 to receive information from all of the other processes.
Calling `MPI_Comm_rank()` returns the rank of the process calling it within `MPI_COMM_WORLD`.

Go ahead and run the code now:

~~~
$ mpiexec -n 4 ./hello
~~~
{: .language-bash}

~~~
World Size: 4   Rank: 1
World Size: 4   Rank: 0
World Size: 4   Rank: 2
World Size: 4   Rank: 3
~~~
{: .output}

As you can see, the `MPI_Comm_size()` function outputs 4, which is the total number of ranks we told `mpiexec` to run with (through the `-n` argument).
Each of the processes is assigned a rank in the range of 0 to 3.

As you can see, the ranks don't necessarily print out their messages in order; whichever rank completes the `cout` first will print out its message first.
If you run the code again, the ranks are likely to print thier message in a different order:

~~~
World Size: 4   Rank: 2
World Size: 4   Rank: 0
World Size: 4   Rank: 3
World Size: 4   Rank: 1
~~~
{: .output}

You can also try rerunning with a different value for the `-n` `mpiexec` argument.
For example:

~~~
$ mpiexec -n 2 ./hello
~~~
{: .language-bash}

~~~
World Size: 2   Rank: 0
World Size: 2   Rank: 1
~~~
{: .output}



### Error Handling with MPI

If an error forces an MPI program to exit, it should **never** just call `return` or `exit`.
This is because calling `return` or `exit` might terminate one of the MPI processes, but leave others running (but not doing anything productive) indefintely.
If you're not careful, you could waste massive amounts of computational resources running a failed calculation your thought had terminated.
Instead, MPI-parallelized codes should call `MPI_Abort()` when something goes wrong, as this function will ensure all MPI processes terminate.
The `MPI_Abort` function takes two arguments: the first is the communicator corresponding to the set of MPI processes to terminate (this should generally be `MPI_COMM_WORLD`), while the second is an error code that will be returned to the environment.

It is also useful to keep in mind that all MPI functions have a return value.
If this value is equal to `MPI_SUCCESS`, the function was executed successfully; otherwise, the function call failed.
By default, MPI functions automatically abort if they encounter an error, so you'll only ever get a return value of `MPI_SUCCESS`.
If you want to handle MPI errors yourself, you can call `MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN)`; if you do this, you must check the return value of every MPI function and call `MPI_Abort` if it is not equal to `MPI_SUCCESS`.
For example, when initializing MPI, you might do the following:

~~~
if (MPI_Init(&argc,&argv) != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, 1);
~~~
{: .language-cpp}




## Example 2

### Basic Infrastructure

We will now do some work with the the example in [examples/mpi/average](https://github.com/MolSSI-Education/parallel-programming/tree/gh-pages/examples/mpi/average), which does some simple math with NumPy arrays.
Run the code now.

~~~
$ cd parallel-programming/examples/mpi/average
$ mkdir build
$ cd build
$ cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_Fortran_COMPILER=mpifort ..
$ make
$ ./average
~~~
{: .language-bash}

~~~
Average: 100000001.5
~~~
{: .output}

Let's learn something about which parts of this code account for most of the run time.
MPI provides a timer, `MPI_Wtime()`, which returns the current walltime.
We can use this function to determine how long each section of the code takes to run.

For example, to determine how much time is spent initializing array `a`, do the following:

~~~
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
~~~
{: .language-cpp}

As the above code indicates, we don't really want every rank to print the timings, since that could look messy in the output.
Instead, we have only rank 0 print this information.
Of course, this requires that we add a few lines near the top of the code to initialize MPI and query the rank of each process:

~~~
  // Initialize MPI
  int world_size, my_rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
~~~
{: .language-cpp}

Also determine and print the timings of each of the other sections of the code: the intialization of array `b`, the addition of the two arrays, and the final averaging of the result.
Your code should look something like this:

~~~
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
~~~
{: .language-cpp}

Now compile and run the code again:

~~~
$ make
$ ./average
~~~
{: .language-bash}

~~~
Initialize a time: 0.544075
Initialize b time: 0.624939
Add arrays time: 0.258915
Average result time: 0.266418
Average: 100000001.5
~~~
{: .language-output}


### Point-to-Point Communication

You can try running this on multiple ranks now:

~~~
$ mpiexec -n 4 ./average
~~~
{: .language-bash}

~~~
Initialize a time: 0.640894
Initialize b time: 0.893775
Add arrays time: 1.38309
Average result time: 0.330192
Average: 100000001.5
~~~
{: .output}

Running on multiple ranks doesn't help with the timings, because each rank is duplicating all of the same work.
In some ways, running on multiple ranks makes the timings worse, because all of the processes are forced to compete for the same computational resources.
Memory bandwidth in particular is likely a serious problem due to the extremely large arrays that must be accessed and manipulated by each process.
We want the ranks to cooperate on the problem, with each rank working on a different part of the calculation.
In this example, that means that different ranks will work on different parts of the arrays `a` and `b`, and then the results on each rank will be summed across all the ranks.

In this section, we will handle the details of the communication between processes using *point-to-point* communication.
Point-to-point communication involves cases in which a code explicitly instructs one specific process to send/recieve information to/from another specific process.
The primary functions associated with this approach are `MPI_Send()` and `MPI_Recv()`, which are involve the following arguments:

~~~
int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,  MPI_Comm comm);
~~~
{: .language-cpp}

* `buf` --- pointer to the start of the buffer being sent
* `count` --- number of elements to send
* `datatype` --- MPI data type of each element
* `dest` --- rank of destination process
* `tag`  --- message tag
* `comm` --- the communicator to use

~~~
int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
~~~
{: .language-cpp}

* `buf` --- pointer to the start of the buffer to receive the message
* `count` --- maximum number of elements the buffer can hold
* `datatype` --- MPI data type of each element
* `source` --- rank of source process ---  `MPI_ANY_SOURCE` matches any process
* `tag`  --- message tag (integer `>= 0`) --- `MPI_ANY_TAG` matches any tag
* `comm` --- the communicator to use
* `status` --- pointer to the structure in which to store status

We need to decide what parts of the arrays each of the ranks will work on; this is more generally known as a rank's workload.
Add the following code just before the initialization of array `a`:

~~~
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
~~~
{: .language-cpp}

In the above code, `my_start` and `my_end` represent the range over which each rank will perform mathematical operations on the arrays.

We'll start by parallelizing the code that averages the result.
Update the range of the `for` loop in this part of the code to the following:

~~~
  for (int i=my_start; i<my_end; i++) {
~~~
{: .language-cpp}

This will ensure that each rank is only calculating elements `my_start` through `my_end` of the sum.
We then need the ranks to communicate their individually calculated sums so that we can calculate the global sum.
To do this, add the following immediately after the end of the `for` loop:

~~~
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
~~~
{: .language-cpp}

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

~~~
$ make
$ ./average
~~~
{: .language-bash}

~~~
Initialize a time: 0.63251
Initialize b time: 1.31379
Add arrays time: 1.89099
Average result time: 0.100575
Average: 100000001.5
~~~
{: .language-output}

You can see that the amount of time spent calculating the average has indeed gone down.

Parallelizing the part of the code that adds the two arrays is much easier.
All you need to do is update the range over which the `for` loop iterates:

~~~
  for (int i=my_start; i<my_end; i++) {
~~~
{: .language-cpp}

Now compile and run the code again:

~~~
$ make
$ ./average
~~~
{: .language-bash}

~~~
Initialize a time: 0.636685
Initialize b time: 1.66542
Add arrays time: 0.466888
Average result time: 0.0871116
Average: 100000001.5
~~~
{: .language-output}

The array addition time has gone down nicely.
Surprisingly enough, the most expensive part of the calculation is now the initialization of the arrays `a` and `b`.
Updating the range over which those loops iterate speeds up those parts of the calation:

~~~
  // Initialze a
  for (int i=my_start; i<my_end; i++) {
...
  // Initialize b
  for (int i=my_start; i<my_end; i++) {
~~~
{: .language-python}

~~~
$ make
$ ./average
~~~
{: .language-bash}

~~~
Initialize a time: 0.159471
Initialize b time: 0.183946
Add arrays time: 0.193497
Average result time: 0.0847806
Average: 100000001.5
~~~
{: .output}



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

~~~
  double *a = new double[ workloads[my_rank] ];
  for (int i=0; i<workloads[my_rank]; i++) {
    a[i] = 1.0;
  }
~~~
{: .language-cpp}

Replace the initialization of `b` with:

~~~
  double *b = new double[ workloads[my_rank] ];
  for (int i=0; i<workloads[my_rank]; i++) {
    b[i] = 1.0 + double(i + my_start);
  }
~~~
{: .language-cpp}

Replace the range of the loops that add and average the arrays to `for (int i=0; i<workloads[my_rank]; i++)`.

Now compile and run the code again:

~~~
$ make
$ ./average
~~~
{: .language-bash}

~~~
Initialize a time: 0.16013
Initialize b time: 0.176896
Add arrays time: 0.190774
Average result time: 0.0871552
Average: 100000001.5
~~~
{: .language-output}


### Collective Communication

Previously, we used **point-to-point communication** (i.e. `MPI_Send` and `MPI_Recv`) to sum the results across all ranks:

~~~
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
~~~
{: .language-cpp}

MPI provides many **collective communication** functions, which automate many processes that can be complicated to write out using only point-to-point communication.
One particularly useful collective communication function is `MPI_Reduce()`:

~~~
  int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                 MPI_Op op, int root, MPI_Comm comm)
~~~
{: .language-cpp}

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

~~~
  double partial_average = average;
  MPI_Reduce(&partial_average, &average, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
~~~
{: .language-cpp}

Compiling and running with this change should produce the same results as before.

Note that in addition to enabling us to write simpler-looking code, collective communication operations tend to be faster than what we can achieve by trying to write our own communication operations using point-to-point calls.

{% include links.md %}
