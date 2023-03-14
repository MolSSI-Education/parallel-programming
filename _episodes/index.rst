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

Workshop lessons
================

.. list-table:: 
    :widths: 10 10 50 50
    :header-rows: 1
    :stub-columns: 1

    * - 
      - Lesson Title
      - Questions
      - Key Points
    * - 1
      - Set up
      - 
      - Download files required for the lesson
    * - 2
      - Introduction
      - What is parallelization and how does it work?
      - 
         * Understand the motivation for parallelizing code.
         * Understand how machine architecture affects the ability to parallelize code.
         * Be aware of the types of parallelization common in computational chemistry.
    * - Stub Row 3
      - Row 3
      - Column 2
      - Column 3 long. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nam sit amet mauris arcu.

.. toctree::
   :maxdepth: 2
   :hidden:
   :titlesonly:
  
   setup
   00-intro
   00-distributed-memory
   00-shared-memory