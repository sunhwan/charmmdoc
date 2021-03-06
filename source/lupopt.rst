.. py:module:: lupopt

===========================
Low Energy Path OPTmization
===========================

This method optimizes a low energy path between a series of molecular
structures.  Energy minimization is done with constraints on center of mass
translation, rotation and orthogonality of step to path vector.

* Reference   : Choi, C. and Elber, R., J. Chem. Phys. 94:751  (1991)
* Source Code : rxncor/lupopt.src

By Krzysztof Kuczera, 12-Mar-1997, Lawrence, KS.


.. _lupopt_syntax:

Syntax for the LUPOpt Command
-----------------------------

::

   LUPOpt [NPATh integer] [UOUT integer] [INIT integer] -
          [EPSEner real] [MAXCycle integer] [STEP real] [IPVOpt integer] -
          [LPrint integer]

   [for 'INIT 2' this line should be followed directly by NPATH lines
    containing names of formatted CHARMM COOR files, no blank lines]

========  =======    ====================================================
Variable  Default    Meaning
========  =======    ====================================================
NPATH     MXPATH     Number of path points
UOUT      21         Unit number for output trajectory with optimized
                     path
INIT      1          Initialization mode:

                     * 1 - straight line in Cartesian space from 
                       MAIN to COMP coordinates
                     * 2 - read path from set of files, file names 
                       supplied below, 1 per line, no blank lines
EPSE      0.001      Structure will be classified as converged if 
                     energy change is lower than EPSE in one step
MAXC      100        Number of path optimization cycles. Each cycle
                     involves making one SD step for each of the
                     structures 2,3,...,NPATH-1
STEP      0.01       Length of optimization step in CHARMM units
IPVO      1          Path vector option

                     * 1 - standard option, path vector is I -> I+1
                     * 2 - symmetric option, path vector is I-1 -> I+1
LPRINT    1          Frequency of printing out path energies
========  =======    ====================================================


.. _lupopt_description:

Algorithm description
---------------------

Work is in Caretsian space, path optimized by constrained steepest
descent.  See C.Choi & R.Elber, J.Chem.Phys. 94:751-760 (1991).

1. An initial conformational path is read in 

   a) linear interpolation between MAIN and COMP coordinates
   b) series of structures in CHARMM COOR files (see LUPINI)
      structures I=1,2,...,NPATH
2. For each structure inside path I=2,3,...,NPATH-1

   a) Path vector is computed 
   b) A steepest descent step is taken, subject to
      rigid-body and path constraints (see LUPCNS)
   c) Step is accepted if energy decreases
   d) Convergence is checked by monitoring energy change
   
3. If procedure has converged along whole path, stop;
   otherwise return to step 2, possibly decreasing step.


Path initialization 
-------------------

Option 1 is simplest, but may lead to very poor initial guess, even
for buatne t->g!

Option 2 is more involved, but allows for more flexibility.
E.g. model structires along a straight line in dihedral space may
be generated from CHARMM or QUANTA, or a TRAVEL path may be input.

Path optimization
-----------------

For each structure gradient components along 7 LUP constraints are
removed (see LUPCNS) and a steepest descent (SD) step is performed.
The step length is formally STEP, due to CHARMM conventions it is
actually STEP*SQRT(NATOM), i.e. 

::

      Xi -> Xi - (STEP/GNORM)*Gi

where Xi -coordinate, Gi - gradient component,

::

      GNORM = SQRT[(Sum_i Gi**2)/Natom]

If the energy of the new structure is lower than before, the step is
accepted and new coords are stored on HEAP; if the energy is higher,
the step is rejected and coordinates are reset.
If steps are rejected for more than half of the structures in a
cycle, STEP is divided by 2.

Progress in path optimization is monitored by printing out energies
of all current structures every LPRI cycles.


Path constraints
----------------

The 7 constraints involve : 1-6 : rigid body translations and
rotations and 7: the path vector. See LUPCNS

As described by Ron Elber, the constraints are linear in Cartesian
coordinates, i.e. their Cartesian gradinets are 3*Natom dimensional
constant vectors. To simplify procedure the vectors are
orthonormalized. Gradient projections along these vectors are then
eliminated.

.. _lupopt_memory:

Memory Usage
------------

The INTEGER*4 arrays of heap pointers IBX,IBY,IBZ
are given a fixed dimension in LUPOPT:  MXPATH=200.
This can be changed by hand.

Memory could be a problem for large systems.

::

   The path structures are stored on the HEAP : NPATH*3*NATOM REAL*8
   Constraint vectors, also on heap             7*3*NATOM     REAL*8
   One working coordinate set                   3*NATOM       REAL*8

   Altogether:                              (NPATH+8)*3*NATOM REAL*8
   + neglible extras for pointers and coorio
