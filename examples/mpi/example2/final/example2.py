from mpi4py import MPI
import numpy as np

if __name__ == "__main__":

    # Get basic information about the MPI communicator
    world_comm = MPI.COMM_WORLD
    world_size = world_comm.Get_size()
    my_rank = world_comm.Get_rank()

    N = 10000000

    # determine the workload of each 
    workloads = [ N // world_size for i in range(world_size) ]
    for i in range( N % world_size ):
        workloads[i] += 1
    my_start = 0
    for i in range( my_rank ):
        my_start += workloads[i]
    my_end = my_start + workloads[my_rank]

    # initialize a
    start_time = MPI.Wtime()
    a = np.ones( workloads[my_rank] )
    end_time = MPI.Wtime()
    if my_rank == 0:
        print("Initialize a time: " + str(end_time-start_time))

    # initialize b
    start_time = MPI.Wtime()
    b = np.zeros( workloads[my_rank] )
    for i in range( my_start, my_end ):
        b[i - my_start] = 1.0 + i
    end_time = MPI.Wtime()
    if my_rank == 0:
        print("Initialize b time: " + str(end_time-start_time))

    # add the two arrays
    start_time = MPI.Wtime()
    for i in range( my_start, my_end ):
        a[i - my_start] = a[i - my_start] + b[i - my_start]
    end_time = MPI.Wtime()
    if my_rank == 0:
        print("Add arrays time: " + str(end_time-start_time))

    # sum the result
    start_time = MPI.Wtime()
    sum = 0.0
    for i in range( my_start, my_end ):
        sum += a[i - my_start]

    # sum the values across all ranks
    sum = np.array( [sum] )
    world_sum = np.zeros( 1 )
    world_comm.Reduce( [sum, MPI.DOUBLE], [world_sum, MPI.DOUBLE], op = MPI.SUM, root = 0 )
    if my_rank == 0:
        print("total: " + str(world_sum))

    # average the result
    average = world_sum / N
    end_time = MPI.Wtime()
    if my_rank == 0:
        print("Average result time: " + str(end_time-start_time))

    if my_rank == 0:
        print("Average: " + str(average))
