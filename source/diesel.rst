.. py:module:: diesel

============================================================================================
Combined Quantum Mechanical and Molecular Mechanics Method Based on DIESEL(GAMESS) in CHARMM
============================================================================================

by Milan Hodoscek
(milan@helix.nih.gov,milan@cmm.ki.si)

Multi reference CI program DIESEL is connected to CHARMM
program in a QM/MM method. To obtain the integrals for input to
DIESEL program it is run from the GAMEss command. 

.. _diesel_description:

Syntax
------

The DIESEL QM potential is initialized with the GAMEss command.

::

   GAMEss   DIESel <int> <int> ... / for the rest of options see :doc:`gamess` /

In order to run DIESEL the standard GAMEss command must be used
with the added DIESel keyword. The integer numbers after this keyword
represent which energy is used in the CHARMM code for further
processing.

DIESEL is the program to perform multi reference CI calculations.


.. _diesel_usage:

Usage
-----

In order to run DIESEL with CHARMM one has provide separate
input files for GAMESS (see :doc:`gamess`) and for DIESEL. The
information provided by GAMESS for DIESEL is the file which contains
MO one and two electron integrals. In order to obtain such integrals
one must specify RUNTYPE=MCSF in GAMESS control input. The following is
an example for water with 3-21g basis set:

::

 $CONTRL COORD=UNIQUE NOSYM=1 ICHARG=0 SCFTYP=mcscf $END
 $SYSTEM MEMORY=1400000 TIMLIM=100000 $END
 $BASIS  GBASIS=N21 NGAUSS=3  $END
 $DET    NCORE=1 NACT=6 NELS=8 $END
 $MCSCF  maxit=1 micit=1 $end
 $DATA


 $END

This produces the necessary moints file which can be read by DIESEL
input. This should be improved because we really don't need to run
MCSCF, especially because it is very memory demanding. But it is
currently the cheapest way to produce integrals over MO basis, without
extensive modification of GAMESS. ???

Typical DIESEL script would be for example:

::

   #! /bin/sh

   export DIESEL_EXE_DIR=/software/qc/diesel/1.14pre/Binaries/Intel

   cat <<! >diesel.in
   MOIntegralFileFormat             = GAMESSC1
   MOLCASRootDir                    = `pwd`
   TempDir                          = `pwd`

   NumberOfElectrons                = 8
   Multiplicities                   = { 1 }
   IrReps                           = { 0 }
   Roots                            = { 1 2 3 4  }

   SelectionThresholds              = { 1 1e-3 1e-5 }
   MaxDavidsonIters                 = 40
   MaxHamiltonStorageMem            = 100MB

   !

   $DIESEL_EXE_DIR/diesel <diesel.in >diesel.out 2>diesel.prot.out

For complete information about DIESEL input see its own user's
guide. In order to provide the correct input files for GAMESS and
DIESEL one has to specify the following ENVIronment commands in CHARMM
input script.

::

   envi input   "h2o.str"
   envi output  "h2o.gms"
   envi punch   "test.dat"
   envi aoints  "test.f8"
   envi moints  "test.f9"
   envi dictnry "test.f10"
   envi work15  "test.f15"
   envi dasort  "test.f20"
   envi dafl30  "test.f30"
   envi jkfile  "test.f23"
   envi casints "test.f13"
   envi civectr "test.f12"

   ! For DIESEL
   envi dieselscript "CI.job"
   envi dieselout    "dies.out"

DIESEL provides many energies (for various multiplicities and roots)
in the same run. In order to get one of them for further processing
(minimization for example; but be careful: no derivatives are
available so one has to do very costly numerical calculation) the two
integers after DIESEL keyword are for multiplicity and root,
respectively.

.. note::
   For complete example look at test/c28test/dieseltst.inp]

.. _diesel_installation:

Installation
------------

To obtain the program write to Michael Hanrath
(michael.hanrath@uni-koeln.de). Since the program is written in C++ it was
not practical to put it under CHARMM tree and compile them
together. CHARMM only knows how to execute the script which runs
DIESEL. DIESEL itsel is also just a driver for other programs. The
only input it needs from GAMESS is the one and two electron binary
file. 



.. _diesel_status:

DIESEL/CHARMM interface status (February 2001)
----------------------------------------------

- no derivatives in DIESEL
- C1 symetry 

Problems to be solved:

- avoid running MCSCF
- cleanup the file name mess.

.. _diesel_implementation:

Implementation
--------------

C++ is not very practical to compile with fortran programs so
CHARMM/GAMESS and DIESEL are completely separated. When the integrals
are transformed to MO basis by GAMESS, CHARMM calls system routine to
run shell script for DIESEL. In it one has to specify the path to
DIESEL distribution.

