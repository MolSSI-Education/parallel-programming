Parallel Programming
--------------------

This course by `The Molecular Sciences Software Institute <https://molssi.org/>`_ (MolSSI) teaches the fundamentals of parallel programming techniques, with an emphasis on MPI and OpenMP parallelization.
It assumes basic familiarity with Python programming, which is the subject of the MolSSI `Python Scripting for Computational Molecular Science <https://education.molssi.org/python_scripting_cms/>`_ lessons.
Episodes 4 and 6 assume basic familiarity with C++.
To see the full MolSSI's education mission statement, please see `here <http://molssi.org/>`_.
This lesson is under continual development, please report issues to the 
`workshop repository <https://github.com/MolSSI-Education/parallel-programming>`_. 

If you see a subject you would like to contribute to, submit a pull request!

.. admonition:: Prerequisites
   :class: attention

   Students should be familiar with opening the Terminal window and creating and navigating files in bash.

Workshop Lessons
================

Set-Up
#######
.. csv-table:: 
  :file: csv_tables/setup.csv
  :header-rows: 1

Introduction
############
.. csv-table:: 
  :file: csv_tables/intro.csv
  :header-rows: 1

Distributed-Memory Parallelization
##################################
.. csv-table:: 
  :file: csv_tables/distributed-memory.csv
  :header-rows: 1

Shared-Memory Parallelization
##################################
.. csv-table:: 
  :file: csv_tables/shared-memory.csv
  :header-rows: 1

.. toctree::
   :maxdepth: 2
   :hidden:
   :titlesonly:

   setup
   01-introduction
   00-distributed-memory
   00-shared-memory