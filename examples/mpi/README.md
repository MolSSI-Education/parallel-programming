## Setting up your environment

All of the MPI examples will use [`mpi4py`](https://mpi4py.readthedocs.io/en/stable/) in Python.

**$ conda create --name sss_parallel
source activate sss_parallel
conda install -c conda-forge Python=3.7 mpi4py numpy**

If you are using `Windows`, you will also need to install [`MS-MPI`](https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi).
Click the link [`here`](https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi), download an MS-MPI install binary from the website, and use it to install MS-MPI onto your machine.

Let's confirm that everything is working.
First, check that you have a working version of the `mpiexec` command:

**$ mpiexec --version**

> mpiexec (OpenRTE) 3.1.4

The exact output isn't critical, as long as some version of `mpiexec` is found.

Now, check that the `mpi4py` package installed correctly:

**$ mpiexec -n 4 python -c "from mpi4py import MPI;print(MPI.COMM_WORLD.Get_size())" **

> 4\
4\
4\
4



## 1. Example 1

### Writing Hello World

We'll start with the first example in [`mpi/example1`](https://github.com/MolSSI-Education/parallel-programming/tree/gh-pages/examples/mpi/example1), which is a simple Hello World code:

```Python
if __name__ == "__main__":

    print("Hello World!")
```

Acquire a copy of the example files for this lesson, and then run MPI Example 1:

**$ git clone git@github.com:MolSSI-Education/parallel-programming.git
cd parallel-programming/examples/mpi/example1
python example1.py**

> Hello World!

### Getting Started with MPI

Let's try running this code on multiple processes.
This is done using the `mpiexec` command.
Many environments also provide an `mpirun` command, which usually - but not always - works the same way.
Whenever possible, you should use `mpiexec` and not `mpirun`, in order to guarantee more consistent results.
The general format for lanching a code on multiple processes is:

**$ mpiexec -n <number_of_processes> <command_to_launch_code>

For example, to launch `example1.py` on 4 processes, do:

**$ mpiexec -n 4 python example1.py**

> Hello World!\
Hello World!\
Hello World!\
Hello World!

When you execute the above command, `mpiexec` launches 4 different instances of `python example1.py` simultaneously, which each print "Hello World!".
Typically, as long as you have at least 4 processors on the machine you are running on, each process will be launched on a different processor; however, certain environment variables and optional arguments to `mpiexec` can change this behavior.
Each process runs the code in `example1.py` independently of the others.

It might not be obvious yet, but the processes `mpiexec` launches aren't completely unaware of one another.
The `mpiexec` adds each of the processes to an MPI communicator, which enables each of the processes to send and receive information to one another via MPI.
The MPI communicator that spans all of the processes launched by `mpiexec` is called `MPI.COMM_WORLD`.

In `mpi4py`, communicators are class objects, and we can query information about them through their class functions.
Edit `example1.py` so that it reads as follows:

```Python
from mpi4py import MPI

if __name__ == "__main__":

    world_comm = MPI.COMM_WORLD
    world_size = world_comm.Get_size()
    my_rank = world_comm.Get_rank()

    print("World Size: " + str(world_size) + "   " + "Rank: " + str(my_rank))
```

In the above code we first import `mpi4py`.
Then, we get the communicator that spans all of the processes, which is called `MPI.COMM_WORLD`.
The communicator's `Get_size()` function tells us the total number of processes within that communicator.
Each of these processes is assigned a uniqe rank, which is an integer that ranges from `0` to `world_size - 1`.
The rank of a process allows it to be identified whenever processes communicate with one another.
For example, in some cases we might want rank 2 to send some information to rank 4, or we might want rank 0 to receive information from all of the other processes.
Calling `world_comm.Get_rank()' returns the rank of the process that called it within `world_comm`.

Go ahead and run the code now:

**$ mpiexec -n 4 python example1.py**

> World Size: 4   Rank: 1\
World Size: 4   Rank: 0\
World Size: 4   Rank: 2\
World Size: 4   Rank: 3

As you can see, the `world_comm.Get_size()` function returns 4, which is the total number of ranks we told `mpiexec` to run with (through the `-n` argument).
Each of the processes is assigned a rank in the range of 0 to 3.

As you can see, the ranks don't necessarily print out their messages in order; whichever rank reaches the `print` function first will print out its message first.
If you run the code again, the ranks are likely to print thier message in a different order:

> World Size: 4   Rank: 2\
World Size: 4   Rank: 0\
World Size: 4   Rank: 3\
World Size: 4   Rank: 1

You can also try rerunning with a different value for the `-n` `mpiexec` argument.
For example:

**$ mpiexec -n 2 python example1.py**

>World Size: 2   Rank: 0\
World Size: 2   Rank: 1

## Example 2

### Basic Infrastructure

### Point-to-Point Communication

### Collective Communication


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

## Example 3

## Additional Considerations

### Asynchronous Communication

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

### MPI vs MPI4Py

### MPI Batch Scripts


Here's an example batch job ([`exercises/mpihello.pbs`](https://github.com/wadejong/Summer-School-Materials/blob/master/Tuesday-afternoon/exercises/mpihello.pbs)) for SeaWulf annotated so show what is going on:
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
    mpiexec ./mpihello

    # You can run more things below
~~~

Submit the job from the directory holding your executable (or modify the batch script to use the full path to your executable)
~~~
    qsub mpihello.pbs
~~~


### Weak Scaling and Strong Scaling

### Load Balancing

Data distribution can also be a challenge.  The finite memory of each node is one constraint.  Another is that all of the data needed by a task must be brought together for it to be executed.  A non-uniform data distribution can also lead to communication hot spots (processors that must send/recv a lot of data) or hot links (wires in the network that are heavily used to transmit data).  This last point highlights the role of network topology --- the communication pattern of your application is mapped onto the wiring pattern (toplogy) of your network.  A communication intensive application may be sensive to this mapping, and MPI provides some assistance for this (e.g., see [here](http://wgropp.cs.illinois.edu/courses/cs598-s16/lectures/lecture28.pdf)).  See also bisection bandwidth below.

The performance impact of poor load balance will become apparent at synchronization points (e.g., blocking global communications) where all processes must wait for the slowest one to catch up.

### Latency and Bandwidth

For point-to-point communication, the central concepts are latency (*L*, the time in seconds for a zero-length message) and bandwidth (*B*, speed of data transfer in bytes/second).  These enable us to model the time to send a message of *N* bytes as

<img src="https://latex.codecogs.com/svg.latex?\Large&space;T(N)=L+\frac{N}{B}" title="LB" />

For typical modern computers *L*=1-10us and *B*=1-100Gbytes/s.  It is hard to accurately measure the latency since on modern hardware the actual cost can depend upon what else is going on in the system and upon your communication pattern.  The bandwidth is a bit easier to measure by sending very large messages, but it can still depend on communication pattern and destination.

An important and easy to remember value is *N1/2*, which is the message length necessary to obtain 50% of peak bandwidth (i.e., *T(N)=2N/B*)

<img src="https://latex.codecogs.com/svg.latex?\Large&space;N&#95;{1/2}=LB" title="Nhalf" />

Inserting *L*=1us and *B*=10Gbyte/s, we obtain *N1/2*=10000bytes.

Excercise: Derive a similar a similar formula for the length necessary to acheive 90% peak bandwith.

Bisection bandwidth is another important concept especially if you are doing a lot of communication all at once.  The network connecting the computers (nodes) in your cluster probably does not directly connect all nodes with every other node --- for anything except a small cluster this would involve too many wires (*O(P^2)* wires, *P*=number of nodes).  Instead, networks use simpler topologies with just *O(P)* wires --- e.g., mesh, n-dimension torous, tree, fat tree, etc.  Imagine dividing this network in two halves so as to give the worst possible bandwidth connecting the halves, which you can derive by counting the number of wires that you cut. This is the bisection bandwidth.  If your application is doing a lot of communication that is not local but is uniform and random, your effective bandwidth is the bisection bandwidth divided by the number of processes.  This can be much smaller than the bandwidth obtained by a single message on a quiet machine.  For instance, in a square mesh of P processes, the bisection bandwidth is *O(sqrt(P))* and if all proceses are trying to communicate globally the average available bandwidth is *O(1/sqrt(P))*, which is clearly not scalable.

Thus, the communication intensity and the communication pattern of your application are both important.  Spreading communication over a larger period of time is a possible optimization.  Another is trying to communicate only to processes that are close in the sense of distance (wires traversed) on the actual network.

For global communication, the details are more complicated because a broadcast or reduction is executed on an MPI-implementation-specific tree of processes that is mapped to the underlying network topology.  However, for long messages an optimized implementation should be able to deliver similar bandwidth to that of the point-to-point communication, with a latency that grows roughly logarithmically with the number of MPI processes.

