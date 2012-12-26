.. py:module::domdec

====================
Domain Decomposition
====================

Domain decomposition (abreviated "domdec") is a method of parallelizing
Molecular Dynamics (MD) simulation. In Domain decomposition, the simulation box
is divided into Nx x Ny x Nz sub-boxes. Each CPU (or core) is assigned a home
sub-box. Each CPU is responsible for updating the coordinates of the atoms
residing in its home sub-box. The non-bonded forces are calculated for all
atom pairs in the home box plus around volume Rcut around the home box in
the positive x, y, and z direction. After the force calculation, the required
forces are communicated using MPI to the sub-boxes surrounding the home box.
The communication between the sub-boxes is implemented using the "Eighth shell"
method. For more details on the eighth shell method see:

* K. J. Bowers, R. O. Dror, D. E. Shaw, J. Comput. Phys., 221, 303-329 (2007).

* B. Hess, C. Kutzner, D. van der Spoel, E. Lindahl, J. Chem. Theory Comput., 4,
  435-446 (2008).

Special Notice:

The DOMDEC code for doing molecular dynamics
simulations with CHARMM is an evolving, highly scalable molecular
dynamics engine that has been released in the distribution version,
c37b1, without the usual year of the testing in the developmental
version because of the important speed up it can provide via
multiprocessor parallization relative to the previously available MD
codes in CHARMM. (For a discussion and benchmarks, see the News item
in www.charmm.org.) Although DOMDEC has been thoroughly tested, it is
likely that input script configurations that work fine in conventional
CHARMM may not function the same in DOMDEC, or additional bugs exist
that will be found as more people use the code. It is suggested,
therefore, that before doing long runs with DOMDEC, it be confirmed
that DOMDEC gives the correct results by comparing the results of a
test run with those obtained with one of the standard MD codes in
CHARMM. We note also that DOMDEC uses conventional periodic boundary
conditions rather than the image facility in CHARMM. We suggest that
DOMDEC is the method of choice for long standard dynamics NVT or NPT
runs with PME on up to 256 processors. Additional features will be
announced in updates to the documentation.

.. _domdec_description:

Description
-----------

In order to use the Particle-Mesh Ewald (PME) electrostatics, CHARMM must be
compiled with pref keywords COLFFT. Domdec splits
the CPUs given by the mpirun -command into direct and reciprocal CPUs. The
direct CPUs are responsible for bonded and non-bonded force calculation,
neighbor list search, etc. The reciprocal CPUs are responsible for calculating
the reciprocal part of the PME sum.


.. _domdec_limitation:

Current limitations of DOMDec
-----------------------------

- Must use COLFFT keyword in pref.dat in order to have correct Ewald electrostatics.
- In order to use SHAKE, must use "shake fast" and have FSSHK keyword in pref.dat
- DOMDec must be initialized before (fast) SHAKE is initialized.
- Must use Leapfrog Verlet integrator (in module dynamc.src):
  - LEAP, LANG, and CPT all work.
  - PERT and TSM do not work.
- Only supports orthogonal simulation boxes.
- Supports dmcons, rxncor, and absolute harmonic constraints (cons harm absolute).
- Only supports 3-atom solvent models (TIP3, SPC), e.g. no support for TIP4.
- Heavy atoms with hydrogen bonds must be immediately before the hydrogens in the
  topology.
- Dynamic Load Balancing is in beta phase for Constant Pressure simulations,
  if you have trouble, switch off Dynamic Load Balancing when performing constant
  pressure simulations.
- Minimization (mini -command) does not work with DOMDEC. The best way around this
  is to first do minimization and then turn on DOMDEC for dynamics.
- DOMDEC cannot be turned OFF within the script where the DOMDEC -command was
  given. For example you cannot run dynamics with DOMDEC, turn it off, and run
  dynamics without DOMDEC.
- IMAGe recentering command only operates on molecules that are a single DOMDEC
  group. For example, water molecules will be correctly recentered but most larger
  molecule are not.

Tips for improving performance:

- Saving trajectory (NSAVC, NSAVV) or restart file (ISVFRQ) during dynamics
  requires a all-to-all communication which, when done often, slows down the
  simulation. In most systems, setting NSAVC, NSAVV, and ISVFRQ to a value greater
  than 100 is adviced.
- Make sure the SSE instructions are being used. Look for
  "Using SSE version of non-bonded force loops" in the output.
- Use FFTW or MKL library.
- Try a different number of reciprocal cores to find the optimal value.

.. _domdec_syntax:

Syntax
------

Two ways to invoke domdec:

1. Add CHARMM script command

   ::

     DOMDec [NDIR NX NY NZ]
            [DLB {ON | OFF}]

   to ENERGY command.

   * NDIR NX NY NZ

     Defines the spatial division among processors/cores of the direct space
     calculation into NX x NY x NZ sub-boxes where each core gets a sub-box.
     In this way, NX x NY x NZ cores of the direct nonbond space calculation
     and the remaining cores are reserved for the reciprocal space calculation
     if PME is requested. For example, NDIR 2 2 2 will divide the simulation
     box into 2 x 2 x 2 (=8) sub-boxes. If the mpirun command asked for 12
     cores, then 4 cores would be reserved for simultaneous reciprocal-space
     calculation

     If NDIR is not defined, program will guess the values based on a simple
     algorithm where the number of reciprocal cores is set to 1/4th of the number
     of total cores.

   Since DOMDEC needs to have the whole system and energy parameters set,
   this command MUST be given after:
   - complete psf is constructed
   - cutoffs are specified via a previous nbonds, update, energy, or dynamics
     command
   - periodic system is specified with crystal.
   - Shake must not be initialized before domdec

   .. note::
      The user is responsible for running CHARMM with the correct number of
      CPUs. This means that the "mpirun" command must have enough CPUs to do the
      sub-box division.

   The sub-box sizes are limited by the cut-off and the number of
   sub-boxes as follows:

   BOXX    = system box size in X direction
   NX      = number of sub-boxes in X direction
   RCUT    = non-bonded cut-off + radius of the largest group
   BOXX/NX = sub-box size in X direction

   Then, for NX >= 2, the sub-box size in X direction must satisfy:
   BOXX/NX <= BOXX - 2*RCUT

   If your system violates this restriction, you can try reducing NX to 1 or by
   increasing NX.


.. _domdec_examples:

Examples of using Domdec
------------------------

1. Example with PME:

   PME electrostatics is used (CHARMM compiled with COLFFT):
   CHARMM is run with "mpirun -n 10" and DOMDec is initialized with:

   ::

     ENERGY DOMD NDIR 2 2 2

   This command divides the simulation box into 2 x 2 x 2 (=8) sub-boxes. The
   remaining 2 CPUs are assigned as reciprocal CPUs responsible for the PME
   reciprocal calculation.

2. Example without PME:

   No PME is used:
   CHARMM is run with "mpirun -n 8" and DOMDec is initialized with:

   ::

     ENERGY DOMD NDIR 2 2 2

   This command divides the simulation box into 2 x 2 x 2 (=8) sub-boxes, just
   like in the example 1. However, no reciprocal CPUs are assigned.

   DLB ON/OFF turns dynamic load balancing (DLB) on or off. Using DLB is adviced
   as it improves performance. By default DLB is ON.

3. Example 3

   ::

     ! System setup (psf) is done here

     energy eps 1.0 cutnb 11 cutim 11 ctofnb 9 ctonnb 7.5 vswi -
           Ewald kappa 0.320 pmEwald order 4 fftx 64 ffty 64 fftz 64 -
           domdec ndir 2 2 2 dlb on

     ! NOTE: SHAKE is initialized AFTER domdec command
     shake fast bonh tol 1.0e-8 para

     dynamics leap start timestep 0.002 nstep 100


.. _domdec_installing:

Installing CHARMM with DOMDEC enabled
-------------------------------------

The following install sequences should get you a working executable,
MPI libraries must be installed and mpif90 in your path. In the following,
use whichever architecture is correct for your machine, and whatever size
you choose. You can alter the executable limits at run time (using
dimension script command, see dimension.doc) regardless of how big or small
you compile CHARMM. You can use the following lines to compile CHARMM with
DOMDEC:

With FFTW-3.3 installed:

::

   $ install.com [host] large M fftw nolog +DOMDEC -CMPI

With MKL installed:

::

   $ install.com [host] large M mkl nolog +DOMDEC -CMPI

Without FFTW / MKL:

::

   $ install.com [host] large M nolog +DOMDEC +COLFFT -CMPI

Where [host] = em64t, osx, gnu

Compiling with PGI compiler (e.g. in kraken) with fftw:

::

   $ install.com gnu large M fftw nolog PGF95 +DOMDEC -CMPI NERSC

Notes on compiling with PGI:

PGI C compiler cannot be used to compile the SSE force kernels in
source/nbonds/enb_core_sse.c. If PGI compiler is used on this file, a warning
is issued compile time and the code resorts using the slower Fortran versions
of the force kernels. When compiling parallel CHARMM, "NERSC" flag used in
the above example switches the C compiler to gcc. When compiling in serial,
user has to switch manually to gcc compiler in order take advantage of the
SSE force kernels.

Notes on compiling with Pathscale:

FFTW3 include file fftw3.f03 does not compile correctly with Pathscale
compiler version 3.2.99. This is why FFTW and PATHSCALE pref keywords are
declared mutually exclusive.