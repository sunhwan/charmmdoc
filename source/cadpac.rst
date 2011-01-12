.. py:module::cadpac

====================================================================================
Combined Quantum Mechanical and Molecular Mechanics Method Based on CADPAC in CHARMM
====================================================================================


by Paul Lyne (paul@tammy.harvard.edu)


.. _cadpac_description:

The CADPAC QM potential is initialized with the CADPac command.

::

   CADPac   [REMOve] [EXGRoup] (atom selection)

   REMOve:  Classical energies within QM atoms are removed.

   EXGRoup: QM/MM Electrostatics for link host groups removed.

The syntax of the CADPAC command in CHARMM follows closely that
of the GAMESS command.


.. _cadpac_usage:

Usage
-----

For complete information about CADPAC input see Chapter 1 in the CADPAC
distribution. 

A QM-MM job using CADPAC needs four input files.  The first is the
normal CHARMM input file containing the CADPac command. The second file is
the CADPAC input file specifying the basis set to be used and the Hamiltonian
that is needed. The third and fourth files are libfil.dat and modpot.dat
respectively. These are the library and model potential files that are 
supplied with CADPAC.


Cadpac Input File
-----------------

For the CADPAC input file the following cards must be present: TITLE, BASIS,
ATOMS, RUNTYP, START, FINISH.

=======  ========================================================================
TITLE    The keyword is always at the start of the input file and is followed
         by a one-line title on the next line of the input.

BASIS    This descirbes the basis set to be used for the QM region if a generic
         basis set is required.  Examples include STO3G,321G,631G,321G*,631G*.
         These are the most common.  Other basis sets are descibed in the CADPAC
         documentation. It is also possible to run a calculation using specific
         basis sets for individual atoms. If this feature is required then
         the BASIS keyword should be ommitted and the LIBRARY keyword is used
         for each atom in the QM region. For a more detailed description of the
         library command please refer to the official CADAPC documentation.
         All the basis sets that are supported by CADPAC are found in the files
         libfil.dat and modpot.dat.

ATOMS    This keyword is always required.

RUNTYP   For the purposes of QM-MM calculations this will either be ENERGY
         for a single point calculation or GRADIENT if the forces are also
         required.  For any minimization or dynamics calculations the GRADIENT
         keyword should be used.

START    This keyword is always required.

FINISH   This keyword is always required.
=======  ========================================================================

Hamiltonians
------------

The Hamiltonian is HF unless otherwise specified. The Hamiltonian can be
changed by inseerting the appropriate keyword after the RUNTYP key.

For example

======= =================================================
MP2     Performs an MP2 calculation
MP3     Performs an MP3 calculation
CI      Performs a Configuration Interaction calculation
        (please refer to the official CADAPC manual)
======= =================================================

For DFT calculations use the KOHNSHAM keyword:

=========================  ====================================================
KOHNSHAM LDA MEDIUM GRDWT  Performs an LDA calculation with a medium sized
                           grid for numerical quadrature.

KOHNSHAM BLYP LARGE GRDWT  Performs a non-local BLYP calculation with a large
                           sized grid
=========================  ====================================================

For other functionals see the official CADPAC manual.


CADPAC I/O
----------

CADPAC has hard wired units 1,2 and 3 for the libfil.dat, modpot.dat and
cadpac input file so avoid using these elsewhere in the CHARMM stream.
Other units that CADPAC commonly uses for the grid, integrals etc are
13,14,18,35,53,and 54.


Examples
--------

An example of a CADPAC input file to run with CHARMM:

::

   TITLE           ! Required
   this is a test  ! Put whatever you like on one line
   BASIS STO3G     ! Generic basis set to be used
   ATOMS           ! Required
   GRADIENT        ! Run type. Use this for optimizations
   START           ! Required
   FINISH          ! Required

The above input file tells CADPAC to use an STO-3G basis for the
atoms in the QM region. CADPAC will perform a gradient evaluation each
time that it is called by CHARMM. If you require just a single point 
calculation without gradients just use ENERGY instead of GRADIENT. The input
file above will perform a HF calculation. A DFT calculation is invoked as
follows:

::

   TITLE           ! Required
   this is a test  ! Put whatever you like on one line
   BASIS STO3G     ! Generic basis set to be used
   ATOMS           ! Required
   GRADIENT        ! Run type. Use this for optimizations
   KOHNSHAM LDA MEDIUM GRDWT 
   START           ! Required
   FINISH          ! Required

DF jobs are invoked by the KOHNSHAM card which takes the type of
functional and grid to be used as arguments. In this case an LDA functional is
used. Alternatives include BLYP, B3LYP. For details see the CADPAC 
distribution.

A sample shell script to run CHARMM with CADPAC is:

::

   #!/bin/tcsh -f
   # parameters:
   #               1       data file name
   #
   echo starting
   date
   echo $1
   set HOME= {where CADPAC data files are}
   # data set and output in home directory
   #
   set data=$HOME/$1.inp
   set output=$HOME/$1.out2
   # make a temporary directory to hold the workfiles
   cd /tmp
   mkdir $1
   cd $1
   #   basis set library file assigned to fort.1
   #   pseudopotential library on fort.2
   #   the CADPAC input file is copied to UNIT 3
   cp $HOME/$1.str fort.3
   cp $HOME/$1.par .
   cp  $HOME/libfil.dat fort.1
   #cp  $HOME/modpot.dat fort.2
   #
   #   run the program
   charmm.exe  < $data
   rm -r ../$1

An example file can be found in test/c25test/cwat.inp.  This
input file also uses cwat.str and the sample run script runcwat.


.. _cadpac_status:

Status
------

CADPAC/CHARMM interface status (February 1997)

- CADPAC, GAMESS and QUANTUM keywords cannot coexist in pref.dat

- CADPAC recognizes atoms by their masses as specified in the 
  RTF file

- The program runs on ALPHA, SGI, C90, IBMRS, HPUX platforms.

- There are references to a parallel version in the code. This has not been
  fully tested yet and so won't be included until a future release.
