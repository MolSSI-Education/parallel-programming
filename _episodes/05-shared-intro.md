# Introduction to Shared-Memory Parallelization

````{admonition} Overview
:class: overview

Questions:
- What is shared-memory parallelization?

Objectives:
- Understand the advantages and limitations of shared-memory parallelization.
````


The shared-memory approach to parallelization involves running a calculation using a single process that spawns (i.e., starts running) multiple *threads*, each of which could be running on a different core.
These threads must all be run on the same node; if you need inter-node parallelization, you will need to implement distributed-memory parallelization.
As subsets of the same process, threads by default share access to all of the same memory (although some variables can be made private to individual threads), and you never need to worry about inter-process communication.

Although that might sound convenient, allowing multiple threads to access the same memory locations simultaneously can lead to **very** bizarre bugs.
For example, suppose you have an integer `n = 0` stored in RAM, and two threads try to add 1 to it at the same time.
It is not at all clear what the resulting value of `n` will be.
Back in episode 1, we mentioned that cores do not operate directly on data in RAM, but first fetch the values from RAM into their local cache.
If the two threads fetch the value of `n` from RAM at the same moment, they will both fetch the value `n = 0` and will subsequently both add one to it and produce a value of `n = 1`, which will be the final value in RAM.
But what if one of the threads is able to fetch the value of `n` from RAM, update it, and replace its value in RAM before the other thread can fetch the value of `n` from RAM?
In that case, the second thread will fetch the value `n = 1`, add one to it, and the final value of `n` will be 2.
The exact speed at which each thread accesses the location of `n` in RAM depends on numerous factors that are entirely uncontrollable for a programmer, including the physical location of each core, the algorithms used by the cache, the specific memory location of `n`, the detailed architecture of the cores, any microscopic defects in the cores, or even minor details regarding the physical environment of the cores.
You could literally get different results depending on whether it is a hot day or a cold day!
For these reasons, shared-memory parallelization algorithms tend to be far more challenging to write, debug, and maintain than one might assume.

Despite the subtle complexities, shared-memory parallelization has one significant advantage over distributed-memory parallelization: It tends to require considerably less RAM.
Whereas with distributed-memory parallelization the each process might have its own complete copy of all data, with shared-memory parallelization there is only a single copy of the data, regardless of the number of active threads.

One particularly popular approach to shared-memory parallelization is OpenMP, which we will explore during the hands-on tutorial in the next section.


````{admonition} Key Points
:class: key

- Shared-memory parallelization enables lower memory requirements than distributed-memory parallelization.
- Subtle bugs that are difficult to identify and fix are common when using shared-memory parallelization.
````
