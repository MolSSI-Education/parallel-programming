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

## 1. Hello world 

### Writing hello world

Start from sequential version [`exercises/hello.cc`](https://github.com/wadejong/Summer-School-Materials/tree/master/MPI/exercises/hello.cc)
```c++
    #include <iostream>
    int main() {
        std::cout << "Hello" << std::endl;
        return 0;
    }
```
Build the sequential version with `make hello` or `icpc -o hello hello.cc`.

### Required elements of all MPI programs

* Include `mpi.h` --- older versions of some MPI implementations required it be the first header
* Initialize MPI --- by calling [`MPI_Init`](https://www.mpich.org/static/docs/v3.2/www3/MPI_Init.html), or [`MPI_Init_thread`](https://www.mpich.org/static/docs/v3.2/www3/MPI_Init_thread.html). Usually the first line of your main program will be similar to the following if you are not handling errors yourself (see just below)
```c++
    MPI_Init(&argc,&argv);
```
or this if you are handling errors
```c++
    if (MPI_Init(&argc,&argv) != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, 1);
```

* Finalize MPI --- by calling [`MPI_Finalize`](https://www.mpich.org/static/docs/v3.2/www3/MPI_Finalize.html) usually the penultimate line of your main program will be
```c++
    MPI_Finalize();
```
* Initializing MPI gives us access to the default communicator (`MPI_COMM_WORLD`)
    * An intra communicator encapsulates all information and resources needed for a group of processes to communicate with each other.  For our simple applications we will always being `MPI_COMM_WORLD` but for real applications you should be passing a communicator into all of your routines to enable reuse and interoperability.
    * We will look at communicators in more detail soon --- for now we just need to get the number of processes (by calling [`MPI_Comm_size`](https://www.mpich.org/static/docs/v3.2/www3/MPI_Comm_size.html)) and the rank (`0,1,2,...`) of the current process (by calling [`MPI_Comm_rank`](https://www.mpich.org/static/docs/v3.2/www3/MPI_Comm_rank.html)).
* Note how MPI wants access to the command line arguments (so we must modify the signature of `main`). You can use `NULL` instead of `argv` and `argc` but passing arguments to MPI is very useful.

### Error detection and exit
* MPI functions return `MPI_SUCCESS` on success or an error code (see documentation) on failure depending on how errors are handled.
* By default errors abort (i.e., the error handler is `MPI_ERRORS_ARE_FATAL`).  If you want MPI to return errors for you to handle, you can call [`MPI_Errhandler_set`](https://www.mpich.org/static/docs/v3.2/www3/MPI_Errhandler_set.html)`(MPI_COMM_WORLD, MPI_ERRORS_RETURN)`.
* To abort execution you cannot just `exit` or `return` because there's lots of clean up that needs to be done when running in parallel --- a poorly managed error can easily waste 1000s of hours of computer time.  You must call `MPI_Abort` to exit with an error code.

The new version ([`exercises/mpihello.cc`](https://github.com/wadejong/Summer-School-Materials/tree/master/MPI/exercises/hello.cc/mpihello.cc)) looks like this
```c++
    #include <mpi.h>
    #include <iostream>
    
    int main(int argc, char** argv) {
        MPI_Init(&argc,&argv);
        
        int nproc, rank;
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
        std::cout << "Hello from process " << rank << " of " << nproc << std::endl;
        
        MPI_Finalize();
        return 0;
    }
```

### Compiling MPI programs

Build your parallel version with `make mpihello` or `mpiicpc -o mpihello mpihello.cc`

* MPI provides wrappers for compiling and linking programs because there's a lot of machine specific stuff that must be done.
* For the summer school we are using the Intel C++ compiler and MPI library so we use the command `mpiicpc` (`mpiifort` for the Intel FORTRAN, `mpiicc` for Intel C) but more commonly you might be using the GNU stack (`mpicxx`, `mpicc`, `mpifort`, etc.) (see [here](https://software.intel.com/en-us/mpi-developer-reference-linux-compiler-commands) for more details).

### Running MPI programs

You can run the program sequentially just like any other program.

To run it in parallel we must use the use the `mpirun` command (this is system dependent and a common variant is `mpiexec`) since again there's lots of machine dependent stuff that needs doing.  At its simplest we must tell MPI how many processes we want to use --- here we use 4.
~~~
    mpirun -np 4 ./mpihello1
~~~

Your output might look like
~~~
    Hello from process 1 of 4
    Hello from process 0 of 4
    Hello from process 3 of 4
    Hello from process 2 of 4
~~~
but more likely will look something like
~~~
    Hello from process Hello from process 1 of 4
    Hello from process 2 of 4
    3 of 4
    Hello from process 0 of 4
~~~

We used four processes on the local machine (e.g., your laptop or the cluster login node).  More typically, and to avoid the ire of colleagues, we want to use multiple computers from the cluster.  You can manually provide to `mpirun` a hostfile that tells it which computers to use --- on most clusters this is rarely necessary since a batch system is used to
* time share the computers in the cluster
* queue jobs according to priority, resource needs, etc.


Here's an example batch job ([`exercises/mpihello.pbs`](https://github.com/wadejong/Summer-School-Materials/tree/master/MPI/exercises/mpihello.pbs)) for SeaWulf annotated so show what is going on:
~~~
    #!/bin/bash
    #PBS -l nodes=2:ppn=24,walltime=00:02:00
    #PBS -q molssi
    #PBS -N hello
    #PBS -j oe

    # Above says:
    # - job has 2 (dedicated) nodes with 24 processes per node with a 2 minute max runtime
    # - use the molssi queue
    # - name the job "hello"
    # - merge the standard output and error into one file

    # Output should appear in the file "<jobname>.o<jobnumber>" in the
    # directory from which you submitted the job

    # ================================================
    # If this is not in your .bashrc it needs to be here so that your job
    # uses the same compilers/libraries that you compiled with
    source /gpfs/projects/molssi/modules-intel

    # This will change to the directory from which you submitted the job
    # which we assume below is the one holding the executable
    cd $PBS_O_WORKDIR

    # Uncomment this if you want to see other PBS environment variables
    # env | grep PBS

    # Finally, run the executable using $PBS_NUM_NODES*$PBS_NUM_PPN
    # processes spread across all the nodes
    mpirun ./mpihello

    # You can run more things below
~~~
But I find the comments distracting, so here ([`exercises/mpihello.pbs`](https://github.com/wadejong/Summer-School-Materials/tree/master/MPI/exercises/mpihello_minimal.pbs)) is a minimal version.
~~~
     #!/bin/bash
     #PBS -l nodes=2:ppn=24,walltime=00:02:00
     #PBS -q molssi -N hello -j oe

     source /gpfs/projects/molssi/modules-intel
     cd $PBS_O_WORKDIR
     mpirun ./mpihello
~~~
You can copy and edit the file for your other jobs.  Note that other other systems running PBS will differ.

Submit the job from the directory holding your executable (or modify the batch script to use the full path to your executable)
~~~
    qsub mpihello.pbs
~~~

Useful PBS/Torque commands are
* `qstat` --- see all queued/running jobs
* `qstat -u <username>` --- to see just your jobs
* `qstat -f <jobid?>` --- to see detailed info about a job
* `qstat -Q` and `qstat -q` --- to see info about batch queues (for the summer school only `molssi` is available)
* `qdel <jobid>` --- to cancel a job

##  2. Sending and receiving messages --- point to point communication

### Essential elements
1. Process rank, message tag, MPI data type, communicator size
2. Blocking communication
3. One minimal set of six operations
3. Buffering and safe communication
4. Non-blocking communication
5. Implied weak synchronization
6. Other communication modes (synchronous send, buffered send)

A process is identified by its rank --- an integer `0,1,..,P-1` where `P` is the number of processes in the communicator (`P` is the size of the communicator).

Every message has a `tag` (integer) that you can use to uniquely identify the message.  This implements the CSP concept of message type --- MPI uses the term `tag` in order to distinguish from the datatype being sent (byte, integer, double, etc.).

In a conversation between a pair of processes, messages of the same `tag` are received in the order sent.

But messages from multiple processes can be interleaved with each other.


### Blocking communication

* When a blocking send function ([`MPI_Send`](https://www.mpich.org/static/docs/v3.2/www3/MPI_Send.html)) completes the buffer can immediately be reused without affecting the sent message.  Note that the receiving process may not necessarily have yet received the message.
* When a blocking recv function ([`MPI_Send`](https://www.mpich.org/static/docs/v3.2/www3/MPI_Recv.html)) completes the received message is fully available in the buffer.

```c++
    int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,  MPI_Comm comm);
```
* `buf` --- pointer to the start of the buffer being sent
* `count` --- number of elements to send
* `datatype` --- MPI data type of each element
* `dest` --- rank of destination process
* `tag`  --- message tag
* `comm` --- the communicator to use

```c++
    int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
```
* `buf` --- pointer to the start of the buffer to receive the message
* `count` --- maximum number of elements the buffer can hold
* `datatype` --- MPI data type of each element
* `source` --- rank of source process ---  `MPI_ANY_SOURCE` matches any process
* `tag`  --- message tag (integer `>= 0`) --- `MPI_ANY_TAG` matches any tag
* `comm` --- the communicator to use
* `status` --- pointer to the structure in which to store status

The actual source and tag of the received message can be accessed directly from the status.  Call [`MPI_Get_count`](https://www.mpich.org/static/docs/v3.2/www3/MPI_Get_count.html) to get the count.
```c++
    status.MPI_TAG;
    status.MPI_SOURCE;
    MPI_Get_count( &status, datatype, &count );
```

There's data types for everything, and you can define you own including non-contiguous data structures:

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


#### Exercise:

Write a program to send an integer (`=99`) from process 0 to process 1 and verify the value is correct.  If it is correct, then print "OK" and terminate correctly, otherwise abort.  Run your program.

Hint: start by copying `exercises/mpihello.cc`.

#### Exercise:

Write a program that has two processes exchanging a buffer of length `N` bytes.  Initialize the buffers to the process rank and verify the exchange happened correctly (i.e., the elements in the buffer received by process 1 should have the value 0).  Try `N=10**n` for `n=0,1,2,3,4,5,6,7,8` (i.e., go from small to very large messages).

There are several ways (both correct and incorrect) of writing this program.  You might try the simplest option first --- each process first sends its buffer and then receives its buffer.

Note that this is such a common operation that there is a special ([`MPI_Sendrecv`](https://www.mpich.org/static/docs/v3.2/www3/MPI_Sendrecv.html)) operation to make this less error prone, less verbose, and to enable optimizations.

#### Exercise:

Write a program to send an integer (`=99`) around a ring of processes (i.e., `0` sends to `1`, `1` sends to `2`, `2` sends to `3`, ..., `P-1` sends to `0`, with process 0 verifying the value is correct.  Your program should work for any number of processes greater than one.


### One minimal set of six operations

~~~
    MPI_Init
    MPI_Finalize
    MPI_Comm_size
    MPI_Comm_rank
    MPI_Send
    MPI_Recv
~~~

A timer is also useful --- [`MPI_Wtime`](https://www.mpich.org/static/docs/v3.2/www3/MPI_Wtime.html) returns a high precision wall clock (elapsed) time.  Note that clocks on each process are **not** synchronized.

### Non-blocking (asynchronous) communication

* When a non-blocking send function ([`MPI_Isend`](https://www.mpich.org/static/docs/v3.2/www3/MPI_Isend.html)) completes, the user must not modify the send buffer until the request is known to have completed (e.g., using ([`MPI_Test`](https://www.mpich.org/static/docs/v3.2/www3/MPI_Test.html)) or [`MPI_Wait`](https://www.mpich.org/static/docs/v3.2/www3/MPI_Wait.html)).

* When a non-blocking recv function ([`MPI_Irecv`](https://www.mpich.org/static/docs/v3.2/www3/MPI_Irecv.html)) completes, any message data is not completely available in the buffer until the request is known to have completed.

```c++
    int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,  MPI_Comm comm, MPI_Request *request);
    int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source,  int tag, MPI_Comm comm, MPI_Request *request);
```
* `request` --- a pointer to structure that will hold the information and status for the request

```c++
    int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status);
    int MPI_Wait(MPI_Request *request, MPI_Status *status);
    int MPI_Cancel(MPI_Request *request);
```
* `request` --- a pointer to the request being tested/waited-upon/cancelled
* `flag` --- a pointer to an int that will be non-zero (true) if the operation has completed
* `status` --- a pointer to the structure in which status will be stored if the operation has completed
* See also `MPI_Waitall`, `MPI_Waitany` and `MPI_Waitsome` for waiting on multiple requests


### Other P2P communication modes

*  Buffered send --- provide sender-side buffering to ensure a send always completes and to make memory-management more explicit
*  Synchronous send --- completes on the sender-side when the receive has also completed
*  Ready send --- if you know a matching receive has already been posted this enables optimizations, and this style of programming is explicitly safe from memory/buffer issues
*  One-sided operations --- remote memory access (RMA)


## 3. Global or collective operations

### Essential elements
1. Broadcast
2. Reduction
3. Implied global synchronization
4. Other global operations

In constrast to point-to-point operations that involve just two processes, global operations move data between **all** processes asscociated with a communicator with an implied **synchronization** between them.  All processes within a communicator are required to invoke the operation --- hence the alternative name "collective" operations.

Many chemistry, materials, and biophysics applications are written using only global operations to
* share/replicate information between all processes by broadcasting, and
* compute sums over (partial) results computed by each processes.

This approach does not necessarily scale to the very largest supercomputers, but can suffice for many needs.

We introduce broadcast and reduction, and then work through an example.

### Broadcast

([`MPI_Bcast`](https://www.mpich.org/static/docs/v3.2/www3/MPI_Bcast.html)) broadcasts a buffer of data from process rank `root` to all other processes.  Once the operation is complete within a process its buffer contains the same data as that of process `root`.
```c++
    int MPI_Bcast (void *buffer, int count, MPI_Datatype datatype, int root,  MPI_Comm comm)
```
* `root` --- the process that is broadcasting the data --- this **must** be the same in all processes

### Reduction

To combines values from all processes with a reduction operation either to just process `root` (([`MPI_Reduce`](https://www.mpich.org/static/docs/v3.2/www3/MPI_Reduce.html))) or distributing the result back to all processes (([`MPI_Allreduce`](https://www.mpich.org/static/docs/v3.2/www3/MPI_Allreduce.html))).

```c++
    int MPI_Reduce (const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);

    int MPI_Allreduce (const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm );
```
* `sendbuf` --- a pointer to the buffer that contains the local data to be reduced
* `recvbuf` --- a pointer to the buffer that will hold the result

There are many pre-defined reduction operation and you can also define your own

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

### Exercise

In [`exercises/pi_seq.cc`](https://github.com/wadejong/Summer-School-Materials/tree/master/MPI/exercises/pi_seq.cc) is a (now traditional) Monte Carlo program to compute the value of pi.  Make it run in parallel using broadcast and reduce.  Increase the number of points to demonstrate a speedup.

We will walk through the solution together since this is an important example.

### Exercise:

In [`exercises/trapezoid_seq.cc`](https://github.com/wadejong/Summer-School-Materials/tree/master/MPI/exercises/trapezoid_seq.cc) is a sequential program that uses the trapezoid rule to estimate the value of the integral

<img src="https://latex.codecogs.com/svg.latex?\Large&space;\int&#95;{-6}^{6}\exp(-x^2)\cos(3x)\,dx" title="Amdahl" />

It repeatedly increases the number of points by a factor two until the error is satisfactory.

Please make it run in parallel using MPI.  Initially just make the integration step run in parallel using all-reduce. Then, make only process 0 responsible for deciding if the error is satisfactory (and of course telling everyone else).

We will walk through the solution together since this is an important example.

### Other global operations

There are many other global operations --- barrier, gather, scatter, parallel prefix, etc.

There are also asynchronous variants --- which are interesting because they may permit overlap of multiple global operations, or overlap of work and communication, or reduce the impact of load imbalance.


### Another minimal set of six operations

~~~
    MPI_Init
    MPI_Finalize
    MPI_Comm_size
    MPI_Comm_rank
    MPI_Bcast
    MPI_Reduce
~~~
Or a total of eight if you include
~~~
    MPI_Send
    MPI_Recv
~~~
Or 9 if you want a timer
~~~
    MPI_Wtime
~~~



## 4. Debugging

There are some powerful visual parallel debuggers that understand MPI, but since these can be expensive we are often left just with GDB. There are variety of ways of using GDB to debug a parallel application:

* Intel MPI and MPICH provide an easy mechanism --- just add `-gdb` to your `mpirun` command.  At this point you are interacting with `gdb` attached to each of your processes.  By default your commands are sent to all processes and output is annoted process. You can control which process you are interacting with using the `z` command.  Some more info is in section 7 of [here](http://physics.princeton.edu/it/cluster/docs/mpich2-user-guide.pdf).

* Other MPI implementations have other mechanisms.  E.g., see here for [OpenMPI debugging](https://www.open-mpi.org/faq/?category=debugging).

* A more portable solution that assumes the MPI processes can create an X-window on your computer is `mpirun -np 2 xterm -e gdb executable` which creates an `xterm` for each process in your application.  This works great for a few processes, but does not scale and it can be complicated to get X-windows to work thru firewalls, etc.

## 5. Exercises

1. [easy] Skim through some of the other tutorials and documentation that have links provided above
2. [easy-medium] Write a program to benchmark the performance of reduce, all-reduce, broadcast as a function of both N and P.  Use N=1,2,4,8,...,1024*1024 doubles. And experiment with processes on the same node and on
different nodes (this means setting #nodes and #ppn correctly in the PBS file).
4. [easy] Parallelize Monte Carlo computation of pi starting from [`exercises/pi_seq.cc`](https://github.com/wadejong/Summer-School-Materials/tree/master/MPI/exercises/pi_seq.cc) using global operations
4. [easy] Work through the other various examples in the `exercises/` directory
5. [medium] Parallelize the recursively adaptive quadrature program [`exercises/recursive_seq.cc`](https://github.com/wadejong/Summer-School-Materials/tree/master/MPI/exercises/recursive_seq.cc)
6. [medium-hard] Write MPI versions of the example SCF, VMC, or MD codes in the main [chemistry examples directory](https://github.com/wadejong/Summer-School-Materials/tree/master/MPI/exercises).  This tree includes example programs for Hartree Fock, molecular dynamics (already seen in the OpenMP lecture), and variational quantum Monte Carlo.  Sequential, OpenMP, and MPI versions are provided, and the `README` in each directory gives more details.  There's lots of different approaches so don't take our parallel versions as definitive.
    * VMC is the easiest
    * MD is also easy to get started, but harder to get best performance
    * Hartree-Fock is actually not that hard but the complexity of the code and algorithm can obscure this
    * Look at the `README.md` files in the `mpi` sub-directories to get hints.

{% include links.md %}
