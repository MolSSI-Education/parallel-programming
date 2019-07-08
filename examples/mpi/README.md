
Hello world:
* `hello.cc` --- sequential hello world 
* `mpihello.cc` --- MPI parallel hello world
    * Init, Finalize, rank, size
* `mpihello.pbs` --- annotated PBS batch file to run MPI hello world
    * Writing a PBS batch script
* `mpihello_minimal.pbs` --- minimal PBS batch file to run MPI hello world

Point to point:
* `sendrecv99.cc` --- process 0 sends a value to process 1
    * Simple point-to-point communication
* `exchange1.cc` --- exchange of buffers between a pair of processes
    * Unsafe assumption of buffering
*  `exchange2.cc` --- exchange of buffers between a pair of processes
    * Various safe ways of exchanging buffers
* `ring.cc` --- passing a message around a ring of processes
    * Slightly more complex point-to-pint communication
* `ring-time.cc` --- passing a message around a ring of processes
    * Measures performance as a function of buffer size
* `client-server.cc` --- simple client-server model for a shared counter
    * Illustrates use of wild card receive and client-server model

Global/collective operations:
* `global.cc` --- collective operations
    * Illustrates basic use of reduce, all reduce, and broadcast
* `trapezoid_seq.cc` --- sequential version of trapezoid integration
* `trapezoid_par.cc` --- one possible parallel version of trapezoid integration
    * Illustrates basic use of broadcast+reduce, parallelization of a for loop, distributed logic
* `pi_seq.cc` --- sequential version of Monte Carlo calculation of PI
    * Starting point for exploring basic use of broadcast+reduce, and parallelization of a for loop
* `pi_par.cc` --- parallel version of  Monte Carlo calculation of PI
    * Illustrates use of reduce

Other:
* `recursive.cc` --- sequential version of adaptive quadrature for parallelization using P2P or global operations
* `nest_seq.cc` --- nest of loops to be parallelized
* `nest_par.cc` --- parallel nest
    * Illustrates use of counter to parallelize complex loop nest
