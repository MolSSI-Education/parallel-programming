---
title: "MPI Hands-On - mpi4py"
teaching: 0
exercises: 90
questions:
- "How can I use MPI to parallelize a Python code?"
objectives:
- "Learn how to prepare an environment that includes mpi4py."
- "Learn the basics of writing an MPI-parallelized code."
- "Explore point-to-point and collective MPI operations"
keypoints:
- "Where possible, use collective communication operations instead of point-to-point communication for improved efficiency and simplicity."
- "Intelligent design choices can help you reduce the memory footprint required by MPI-parallelized codes"
---

## 1. Example 1

### Writing Hello World

We'll start with the first example in [`mpi/example1`](https://github.com/MolSSI-Education/parallel-programming/tree/gh-pages/examples/mpi/example1), which is a simple Hello World code:

~~~
if __name__ == "__main__":

    print("Hello World!")
~~~
{: .language-python}

Acquire a copy of the example files for this lesson, and then run MPI Example 1:

~~~
$ git clone git@github.com:MolSSI-Education/parallel-programming.git
$ cd parallel-programming/examples/mpi/example1
$ python example1.py
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
$ mpiexec -n 4 python example1.py
~~~
{: .language-bash}

~~~
Hello World!
Hello World!
Hello World!
Hello World!
~~~
{: .output}

When you execute the above command, `mpiexec` launches 4 different instances of `python example1.py` simultaneously, which each print "Hello World!".

Typically, as long as you have at least 4 processors on the machine you are running on, each process will be launched on a different processor; however, certain environment variables and optional arguments to `mpiexec` can change this behavior.
Each process runs the code in `example1.py` independently of the others.

It might not be obvious yet, but the processes `mpiexec` launches aren't completely unaware of one another.
The `mpiexec` adds each of the processes to an MPI communicator, which enables each of the processes to send and receive information to one another via MPI.
The MPI communicator that spans all of the processes launched by `mpiexec` is called `MPI.COMM_WORLD`.

In `mpi4py`, communicators are class objects, and we can query information about them through their class functions.
Edit `example1.py` so that it reads as follows:

~~~
from mpi4py import MPI

if __name__ == "__main__":

    world_comm = MPI.COMM_WORLD
    world_size = world_comm.Get_size()
    my_rank = world_comm.Get_rank()

    print("World Size: " + str(world_size) + "   " + "Rank: " + str(my_rank))
~~~
{: .language-python}

In the above code we first import `mpi4py`.
Then, we get the communicator that spans all of the processes, which is called `MPI.COMM_WORLD`.
The communicator's `Get_size()` function tells us the total number of processes within that communicator.
Each of these processes is assigned a uniqe rank, which is an integer that ranges from `0` to `world_size - 1`.
The rank of a process allows it to be identified whenever processes communicate with one another.
For example, in some cases we might want rank 2 to send some information to rank 4, or we might want rank 0 to receive information from all of the other processes.
Calling `world_comm.Get_rank()` returns the rank of the process that called it within `world_comm`.

Go ahead and run the code now:

~~~
$ mpiexec -n 4 python example1.py**
~~~
{: .language-bash}

~~~
World Size: 4   Rank: 1
World Size: 4   Rank: 0
World Size: 4   Rank: 2
World Size: 4   Rank: 3
~~~
{: .output}

As you can see, the `world_comm.Get_size()` function returns 4, which is the total number of ranks we told `mpiexec` to run with (through the `-n` argument).
Each of the processes is assigned a rank in the range of 0 to 3.

As you can see, the ranks don't necessarily print out their messages in order; whichever rank reaches the `print` function first will print out its message first.
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
$ mpiexec -n 2 python example1.py
~~~
{: .language-bash}

~~~
World Size: 2   Rank: 0
World Size: 2   Rank: 1
~~~
{: .output}

## Example 2

### Basic Infrastructure


Copy the code in this script in [`example2.py`](https://github.com/MolSSI-Education/parallel-programming/tree/gh-pages/examples/mpi/example2). This script does some simple math with NumPy arrays.
Run the code now.

~~~
$ python example2.py
~~~
{: .language-python}

~~~
Average: 5000001.5
~~~
{: .output}

Let's learn something about which parts of this code account for most of the run time.
MPI4Py provides a timer, `MPI.Wtime()`, which returns the current walltime.
We can use this function to determine how long each section of the code takes to run.

For example, to determine how much time is spent initializing array `a`, do the following:

~~~
    # initialize a
    start_time = MPI.Wtime()
    a = np.ones( N )
    end_time = MPI.Wtime()
    if my_rank == 0:
        print("Initialize a time: " + str(end_time-start_time))
~~~
{: .language-python}

As the above code indicates, we don't really want every rank to print the timings, since that could look messy in the output.
Instead, we have only rank 0 print this information.
Of course, this requires that we add a few lines near the top of the code to query the rank of each process:

~~~
    # get basic information about the MPI communicator
    world_comm = MPI.COMM_WORLD
    world_size = world_comm.Get_size()
    my_rank = world_comm.Get_rank()
~~~
{: .language-python}

Also determine and print the timings of each of the other sections of the code: the intialization of array `b`, the addition of the two arrays, and final averaging of the result.
Your code should look something like this:

~~~
import numpy as np

if __name__ == "__main__":

    # get basic information about the MPI communicator
    world_comm = MPI.COMM_WORLD
    world_size = world_comm.Get_size()
    my_rank = world_comm.Get_rank()

    N = 10000000

    # determine the workload of each rank
    workloads = [ N // world_size for i in range(world_size) ]
    for i in range( N % world_size ):
        workloads[i] += 1
    my_start = 0
    for i in range( my_rank ):
        my_start += workloads[i]
    my_end = my_start + workloads[my_rank]

    # initialize a
    start_time = MPI.Wtime()
    a = np.ones( N )
    end_time = MPI.Wtime()
    if my_rank == 0:
        print("Initialize a time: " + str(end_time-start_time))

    # initialize b
    start_time = MPI.Wtime()
    b = np.zeros( N )
    for i in range( N ):
        b[i] = 1.0 + i
    end_time = MPI.Wtime()
    if my_rank == 0:
        print("Initialize b time: " + str(end_time-start_time))

    # add the two arrays
    start_time = MPI.Wtime()
    for i in range( N ):
        a[i] = a[i] + b[i]
    end_time = MPI.Wtime()
    if my_rank == 0:
        print("Add arrays time: " + str(end_time-start_time))

    # average the result
    start_time = MPI.Wtime()
    sum = 0.0
    for i in range( N ):
        sum += a[i]
    average = sum / N
    end_time = MPI.Wtime()
    if my_rank == 0:
        print("Average result time: " + str(end_time-start_time))
        print("Average: " + str(average))
~~~
{: .language-python}

Now run the code again:

~~~
$ python example2.py
~~~
{: .language-bash}

~~~
Initialize a time: 0.03975701332092285
Initialize b time: 1.569957971572876
Add arrays time: 4.173098087310791
Average result time: 2.609341859817505
Average: 5000001.5
~~~
{: .language-output}


### Point-to-Point Communication

You can try running this on multiple ranks now:

~~~
$ mpiexec -n 4 python example2.py
~~~
{: .language-bash}

~~~
Initialize a time: 0.042365074157714844
Initialize b time: 1.9863519668579102
Add arrays time: 4.9583611488342285
Average result time: 2.9468209743499756
Average: 5000001.5
~~~
{: .output}

Running on multiple ranks doesn't help with the timings, because each rank is duplicating all of the same work.
We want the ranks to cooperate on the problem, with each rank working on a different part of the calculation.
In this example, that means that different ranks will work on different parts of the arrays `a` and `b`, and then the results on each rank will be summed across all the ranks.

We need to decide what parts of the arrays each of the ranks will work on; this is more generally known as a rank's workload.
Add the following code just before the initialization of array `a`:

~~~
    # determine the workload of each rank
    workloads = [ N // world_size for i in range(world_size) ]
    for i in range( N % world_size ):
        workloads[i] += 1
    my_start = 0
    for i in range( my_rank ):
        my_start += workloads[i]
    my_end = my_start + workloads[my_rank]
~~~
{: .language-python}

In the above code, `my_start` and `my_end` represent the range over which each rank will perform mathematical operations on the arrays.

We'll start by parallelizing the code that averages the result.
Update the range of the `for` loop to the following:

~~~
    for i in range( my_start, my_end ):
~~~
{: .language-python}

This will ensure that each rank is only calculating elements `my_start` through `my_end` of the sum.
We then need the ranks to communicate their individually calculated sums so that we can calculate the global sum.
To do this, replace the line `average = sum / N` with:

~~~
    if my_rank == 0:
        world_sum = sum
        for i in range( 1, world_size ):
      	    sum_np = np.empty( 1 )
            world_comm.Recv( [sum_np, MPI.DOUBLE], source=i, tag=77 )
            world_sum += sum_np[0]
        average = world_sum / N
    else:
        sum_np = np.array( [sum] )
        world_comm.Send( [sum_np, MPI.DOUBLE], dest=0, tag=77 )
~~~
{: .language-python}

The `MPI.DOUBLE` parameter tells MPI what type of information is being communicated by the `Send` and `Recv` calls.
In particular, we are sending a array of double precision numbers.
If you are communicating information of a different datatype, consult the following:

|**MPI4Py data type**  |**C data type**     |
|:------------------|:-------------------|
|`MPI.BYTE`           |8 binary digits     |
|`MPI.CHAR`           |char                |
|`MPI.UNSIGNED_CHAR`  |unsigned char       |
|`MPI.SHORT`          |signed short int	 |	 
|`MPI.UNSIGNED_SHORT` |unsigned short int	 |	 
|`MPI.INT`            |signed int          |
|`MPI.UNSIGNED`       |unsigned int	 |	 
|`MPI.LONG`           |signed long int	 |	 
|`MPI.UNSIGNED_LONG`  |unsigned long int	 |	 
|`MPI.FLOAT`          |float               |
|`MPI.DOUBLE`         |double              |

Now run the code again:

~~~
$ mpiexec -n 4 python example2.py
~~~
{: .language-bash}

~~~
Initialize a time: 0.04637002944946289\
Initialize b time: 1.9484930038452148\
Add arrays time: 4.914314031600952\
Average result time: 0.6889588832855225\
Average: 5000001.5
~~~
{: .output}

You can see that the amount of time spent calculating the average has indeed gone down.

Parallelizing the part of the code that adds the two arrays is much easier.
All you need to do is update the range over which the `for` loop iterations:

~~~
    for i in range( my_start, my_end ):
~~~
{: .language-python}

Now run the code again:

~~~
$ mpiexec -n 4 python example2.py
~~~
{: .language-bash}

~~~
Initialize a time: 0.04810309410095215
Initialize b time: 2.0196259021759033
Add arrays time: 1.2053139209747314
Average result time: 0.721329927444458
Average: 5000001.5
~~~
{: .output}

The array addition time has gone down nicely.
Surprisingly enough, the most expensive part of the calculation is now the initialization of array `b`.
Updating the range over which that loop iterates speeds up that part of the calation:

~~~
    for i in range( my_start, my_end ):
~~~
{: .language-python}

~~~
$ mpiexec -n 4 python example2.py
~~~
{: .language-bash}

~~~
Initialize a time: 0.04351997375488281\
Initialize b time: 0.503791093826294\
Add arrays time: 1.2048840522766113\
Average result time: 0.7626049518585205\
Average: 5000001.5
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
    a = np.ones( workloads[my_rank] )
~~~
{: .language-python}

Replace the initialization of `b` with:

~~~
    b = np.zeros( workloads[my_rank] )
    for i in range( workloads[my_rank] ):
        b[i] = 1.0 + ( i + my_start )
~~~
{: .language-python}

Replace the range of the loops that add and sum the arrays to `range( workloads[my_rank] )`.

Run the code again:

~~~
$ mpiexec -n 4 python example2.py
~~~
{: .language-bash}

~~~
Initialize a time: 0.009948015213012695\
Initialize b time: 0.5988950729370117\
Add arrays time: 1.2081310749053955\
Average result time: 0.7307591438293457\
Average: 5000001.5
~~~
{: .output}


### Collective Communication

Previously, we used **point-to-point communication** (i.e. `Send` and `Recv`) to sum the results across all ranks:

~~~
    if my_rank == 0:
        world_sum = sum
        for i in range( 1, world_size ):
            sum_np = np.empty( 1 )
            world_comm.Recv( [sum_np, MPI.DOUBLE], source=i, tag=77 )
            world_sum += sum_np[0]
        average = world_sum / N
    else:
        sum_np = np.array( [sum] )
        world_comm.Send( [sum_np, MPI.DOUBLE], dest=0, tag=77 )
~~~
{: .language-python}

MPI provides many **collective communication** functions, which automate many processes that can be complicated to write out using only point-to-point communication.
In particular, the `Reduce` function allows us to sum a value across all ranks, without all of the above code.
Replace the above with:

~~~
    sum = np.array( [sum] )
    world_sum = np.zeros( 1 )
    world_comm.Reduce( [sum, MPI.DOUBLE], [world_sum, MPI.DOUBLE], op = MPI.SUM, root = 0 )
    average = world_sum / N
~~~
{: .language-python}

The `op` argument lets us specify what operation should be performed on all of the data that is reduced.
Setting this argument to `MPI.SUM`, as we do above, causes all of the values to be summed onto the root process.
There are many other operations provided by MPI, as you can see here:

|**Operation** | Description | Datatype|
|:-------------|:------------|:--------|
|`MPI.MAX`       |maximum      |integer,float|
|`MPI.MIN`       |minimum      |integer,float|
|`MPI.SUM`       |sum          |integer,float|
|`MPI.PROD`      |product      |integer,float|
|`MPI.LAND`      |logical AND  |integer|
|`MPI.BAND`      |bit-wise AND |integer,MPI_BYTE|
|`MPI.LOR`       |logical OR   |integer|
|`MPI.BOR`       |bit-wise OR  |integer,MPI_BYTE|
|`MPI.LXOR`      |logical XOR  |integer|
|`MPI.BXOR`      |bit-wise XOR |integer,MPI_BYTE|
|`MPI.MAXLOC`    |max value and location|float|
|`MPI.MINLOC`    |min value and location|float|

Note that in addition to enabling us to write simpler-looking code, collective communication operations tend to be faster than what we can achieve by trying to write our own communication operations using point-to-point calls.


## Example 3

Next, view [`example3.py`](https://github.com/MolSSI-Education/parallel-programming/tree/gh-pages/examples/mpi/example3) which is a simple Monte-Carlo simulation.
Run the code now.

~~~
$ python example3.py
~~~
{: .language-bash}

~~~
1000 248.52688099543923
2000 10.588491394826892
3000 -0.9309007491547571
4000 -3.8247648102916196
5000 -4.715929587912762
6000 -5.362217832200815
7000 -5.570585267104749
8000 -5.649439720181915
9000 -5.65428738463388
10000 -5.73417919011543
Total simulation time: 21.389078855514526
    Energy time:       21.013432502746582
    Decision time:     0.09333038330078125
~~~
{: .output}

As you can see, the code already has some timings, and the vast majority of time is spent in the calls to 'get_particle_energy'.
That is where we will focus our parallelization efforts.

The function in question is:

~~~
def get_particle_energy(coordinates, box_length, i_particle, cutoff2):

    """
    This function computes the pairwise Lennard Jones energy of two particles in a periodic box.

    Parameters
    ----------
    r_i: list/array
        the potitional vection of the particle i
    r_j: list/array
        the potitional vection of the particle j
    box_length : float/int
        length of simulation box

    Return
    ------
    rij2: float
        the square of the shortest distance between the two particles and their images
    """


    e_total = 0.0

    i_position = coordinates[i_particle]

    particle_count = len(coordinates)

    for j_particle in range(particle_count):

        if i_particle != j_particle:

            j_position = coordinates[j_particle]

            rij2 = minimum_image_distance(i_position, j_position, box_length)

            if rij2 < cutoff2:
                e_pair = lennard_jones_potential(rij2)
                e_total += e_pair

    return e_total
~~~
{: .language-python}

This looks like it should be fairly straightforward to parallelize: it consists of a single `for` loop which just sums the interaction energies of particle pairs.
To parallelize this loop, we need each rank to compute the interaction energies of a subset of these pairs, and then sum the energy across all ranks.

The `get_particle_energy` function is going to need to know some basic information about the MPI communicator, so add the MPI communicator to its parameters:

~~~
def get_particle_energy(coordinates, box_length, i_particle, cutoff2, comm):
~~~
{: .language-python}

Now update the two times `get_particle_energy` is called by `main`:

~~~
        current_energy = get_particle_energy(coordinates, box_length, i_particle, simulation_cutoff2, world_comm)
        ...
        current_energy = get_particle_energy(coordinates, box_length, i_particle, simulation_cutoff2, world_comm)
~~~
{: .language-python}

Place the following at the beginning of `get_particle_energy`:

~~~
    # Get information about the MPI communicator
    my_rank = comm.Get_rank()
    world_size = comm.Get_size()
~~~
{: .language-python}

Change the `for` loop in `get_particle_energy` to the following:

~~~
    for j_particle in range(my_rank, particle_count, world_size):
~~~
{: .languge-python}

The above code will cause each rank to iterate over particles with a stride of `world_size` and an initial offset of `my_rank`.
For example, if you run on 4 ranks, rank 0 will iterate over particles 0, 4, 8, 12, etc., while rank 1 will iterate over particles 1, 5, 9, 13, etc.

We then need to sum the energies across all ranks.
Replace the line `return e_total` with the following:

~~~
    # Sum the energy across all ranks
    e_single = np.array( [e_total] )
    e_summed = np.zeros( 1 )
    comm.Reduce( [e_single, MPI.DOUBLE], [e_summed, MPI.DOUBLE], op = MPI.SUM, root = 0 )

    return e_summed[0]
~~~
{: .language-python}

Try to run it now:

~~~
$ mpiexec -n 4 python example3.py
~~~
{: .language-bash}

~~~
1000 -35480909996.566864
2000 -66252436255523.72
3000 -86936127660856.08
4000 -93141042416342.66
5000 -256171999678073.88
6000 -3.162015453630529e+21
7000 -3.1620181302289283e+21
8000 -3.162018130377518e+21
9000 -3.1620181324457333e+21
10000 -3.1620182854716e+21
Total simulation time: 31.748733043670654
    Energy time:       31.112581253051758
    Decision time:     0.21792912483215332
~~~
{: .output}

That doesn't seem right at all.
What went wrong?

Our call to `Reduce` causes the energies to be summed onto rank 0, but none of the other ranks have the summed energies.
To have the energies reduced to all of the ranks, replace the `Reduce` call with a call to `Allreduce`:

~~~
    comm.Allreduce( [e_single, MPI.DOUBLE], [e_summed, MPI.DOUBLE], op = MPI.SUM )
~~~
{: .language-python}

~~~
$ mpiexec -n 4 python example3.py
~~~
{: .language-bash}

~~~
1000 -5402881.246788438
2000 -5403807.559181325
3000 -5403898.801044374
4000 -5403916.261693102
5000 -5403921.433225453
6000 -5403923.534017933
7000 -5403924.646963553
8000 -5403925.292483066
9000 -5403925.63053995
10000 -5403926.272461226
Total simulation time: 43.26621890068054
    Energy time:       42.664116621017456
    Decision time:     0.16298675537109375
~~~
{: .output}

This still isn't consistent with our previous results.
The problem is that each iteration, the coordinates are updated by randomly displacing one of the particles.
Each rank randomly selects a particle to displace and the displacement vector.
Instead of contributing to the calculation of the same particle's interaction energies for the same nuclear configuration, each rank ends up calculating some of the interaction energies for different atoms and different coordinates.

To fix this, we will have rank 0 be the only rank that randomly selects a particle or a displacement vector.
Rank 0 will then broadcast all necessary information to the other ranks, so that they keep in sync.

Replace the line where the coordinates are intially generated (where `generate_initial_state` is called) with this:

~~~
    if my_rank == 0:
        coordinates = generate_initial_state(method=build_method, num_particles=num_particles, box_length=box_length)
    else:
        coordinates = np.empty([num_particles, 3])
    world_comm.Bcast( [coordinates, MPI.DOUBLE], root = 0 )
~~~
{: .language-python}

At the beginning of the `for` loop in `main` you will see the following code:

~~~
    for i_step in range(n_steps):

        n_trials += 1

        i_particle = np.random.randint(num_particles)

        random_displacement = (2.0 * np.random.rand(3) - 1.0) * max_displacement
~~~
{: .language-python}

Replace the above with the following:

~~~
    for i_step in range(n_steps):

    if my_rank == 0:
            n_trials += 1

            i_particle = np.random.randint(num_particles)
            i_particle_buf = np.array( [i_particle], 'i' )

            random_displacement = (2.0 * np.random.rand(3) - 1.0) * max_displacement
        else:
            i_particle_buf = np.empty( 1, 'i' )
            random_displacement = np.empty( 3 )
        world_comm.Bcast( [i_particle_buf, MPI.INT], root = 0 )
        i_particle = i_particle_buf[0]
        world_comm.Bcast( [random_displacement, MPI.DOUBLE], root = 0 )
        world_comm.Bcast( [coordinates, MPI.DOUBLE], root = 0 )
~~~
{: .language-python}

At the end of the `for` loop in `main` is the following code:

~~~
        start_decision_time = MPI.Wtime()

        delta_e = proposed_energy - current_energy

        accept = accept_or_reject(delta_e, beta)

        if accept:

            total_pair_energy += delta_e
            n_accept += 1
            coordinates[i_particle] += random_displacement

        total_energy = (total_pair_energy + tail_correction) / num_particles

        energy_array[i_step] = total_energy

	if np.mod(i_step + 1, freq) == 0:
            if my_rank == 0:
	       print(i_step + 1, energy_array[i_step])

            if tune_displacement:
                max_displacement, n_trials, n_accept = adjust_displacement(n_trials, n_accept, max_displacement)

        total_decision_time += MPI.Wtime() - start_decision_time
~~~
{: .language-python}

Replace the above with the following:

~~~
        if my_rank == 0:
            start_decision_time = MPI.Wtime()

            delta_e = proposed_energy - current_energy

            accept = accept_or_reject(delta_e, beta)

            if accept:

                total_pair_energy += delta_e
                n_accept += 1
                coordinates[i_particle] += random_displacement

            total_energy = (total_pair_energy + tail_correction) / num_particles

            energy_array[i_step] = total_energy

            if np.mod(i_step + 1, freq) == 0:

                print(i_step + 1, energy_array[i_step])

                if tune_displacement:
                    max_displacement, n_trials, n_accept = adjust_displacement(n_trials, n_accept, max_displacement)

            total_decision_time += MPI.Wtime() - start_decision_time
~~~
{: .language-python}

Try running the code again:

~~~
$ mpiexec -n 4 python example3.py
~~~
{: .language-bash}

~~~
1000 248.52688099525105
2000 10.588491394638726
3000 -0.9309007493429244
4000 -3.824764810479789
5000 -4.715929588100931
6000 -5.3622178323889855
7000 -5.570585267292914
8000 -5.649439720370088
9000 -5.65428738482205
10000 -5.734179190303595
Total simulation time: 13.671964883804321
    Energy time:       12.892877340316772
    Decision time:     0.15127253532409668
~~~
{: .output}

This time the results are much more consistent with what we expect.
