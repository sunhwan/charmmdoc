.. py:module:: repdstr

================================
The Parallel Distributed Replica
================================

By Paul Maragakis and Milan Hodoscek, 2005

Parallel distributed replica allows independent replicated systems
over specified number of processors. 

.. _repdstr_syntax:

::

   REPDstr NREP <int>

     Replicates the system <int> times

   REPDstr NREP <int> EXCHange [UNIT <int>] FREQuency <int> [ [ temp <real> ] temp <real> ... ]

        This is for replica exchange method (see details in c34test/rexc.inp)
   Currently it works so that when exchange occurs all the coordinates
   and velocities are exchanged, thus the lowest temperature is always on
   the first replica. The other method which exchange only the
   temperature will be implemented later.


      UNIT <int> - Optional keyword for exchange output file. Default for
                   int is 6. Must be opened after repd: each replica
                   writes to its own file. If no open statement all
                   replicas write to the same file with the default
                   fortran file name for this unit. Open before the repd
                   command is not very useful.

      FREQuency - when to exchange

      TEMPerature - the temperature of each system. It must be repeated NREP times

        Instead of repeated TEMP comands for equidistant temperature
      intervals one can use:

      STEM <starting-temperature> DTEM <temperature-increase>

 
.. _repdstr_io:

Once REPDstr command is activated the I/O capabilities of CHARMM
are expanded. In standard parallel mode CHARMM deals with I/O only on
the first process. The rest of processes get their data through
network or memory communication. So all I/O statements that are in the
script before REPDstr command are valid only on first process. In
distributed replica mode each replica needs its own and independent
I/O which is enabled after the REPDstr keyword in the input script.

As of May 2007 the following is working:

1. OPEN

   The command open read|write unit 1 card name somefile will open
   somefile_0 for replica 0, somefile_1 for replica 1, etc

2. READ/WRITE (what works? - coor, ...)

   writes to individual files one for each replica. The following is
   supported:
   
   - write title
   - read/write coor

3. STRE stream

   This will open stream_0 for replica 0, stream_1 for replica 1, etc
   It allows CHARMM to run different input files for each replica (or
   group of processors)

4. OUTU unit

   Will stream output to individual files as specified in the open
   command for particular unit. This command should precede STRE
   command if one wants both input and output files for each group of
   processors

5. IF ?MYREP .EQ. n THEN ....

   Works, too. Output only for the processor zero, unless OUTU is
   specified.

6. All the above works in parallel/parallel mode, ie each replica can
   be a parallel job in itself. The numeration of input and output
   files follows the processor numbers and not replica numbers.  NOTE:
   This might change because it maybe confusing? As of now it append
   _mynod to the file names, but _irepdstr would be a better
   choice. Also it would conform with old QM/MM RPATH using standard
   replica.

   The output is written only on a local process 0 for each replica,
   and similar is true also for streamn command.  The limitation is
   that the number of replicas must divide the number of processes
   allocated for parallel. Otherwise it bombs out with level -5.


.. _repdstr_examples:

.. note::

   If you are using mpich-1.2.X then you need to use -p4wd with the
   absolute path or -p4wd `pwd`


Example 1
---------

::

   read psf
   read coor
   nrep 4

This will replicate PSF and coordinates, so after nrep 4 there are
four independent runs with the same coordinates

Example 2
---------

::

   read psf
   nrep 4
   read coor name system.crd

This will replicate PSF but the coordinates will be read from 4
separate files: system.crd_0, system.crd_1, etc

Example 3
---------

::

   nrep 4
   stre inp

This will run for independent CHARMM jobs. Each inp_0, inp_1, inp_3,
and inp_4 can be different input files, with different PSFs,
parameters, etc

Example 4
---------

::

   open write unit 1 card name out

   nrep 4
   outu 1
   stre inp

The same as example 3 but now also output files out_0, out_1, ... will
be written. Note that OUTU must precede STREam command.



