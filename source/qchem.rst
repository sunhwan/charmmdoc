.. py:module:: qchem

====================================================================================
Combined Quantum Mechanical and Molecular Mechanics Method Based on Q-Chem in CHARMM
====================================================================================

H. Lee Woodcock
(hlwood@nih.gov)

based on the GAMESS(US) interface from Milan Hodoscek
(milan@par10.mgsl.dcrt.nih.gov,milan@kihp6.ki.si)
and
the GAMESS(UK) interface from Paul Sherwood
(p.sherwood@dl.ac.uk)

Ab initio program Q-Chem is connected to CHARMM program in a
QM/MM method. This method is based on the interface to the GAMESS (US
version), the latter being an extension of the QUANTUM code which is
described in J. Comp. Chem., Vol. 11, No. 6, 700-733 (1990).

The QM/MM interface between Q-Chem and CHARMM is described in the
following work which should be cited when used...

* H. Lee Woodcock, M. Hodosceck, A. T. B. Gilbert, P. M. W. Gill, H. F. Schaefer,
  B. R. Brooks; Interfacing CHARMM and Q-Chem to perform QM/MM and QM/MM reaction
  pathway calculations. J. Comp. Chem.; 2007; 28 (9); 1485-1502.

.. _qchem_description:

The Q-Chem QM potential is initialized with the QCHEM command
-------------------------------------------------------------

::

  QCHEm    [REMOve] [EXGRoup] [NOGUess] [BLURred [RECAll INT]] [COORdinates]
           [QCLJ] [PARAllel [INT]] [SHESsian] [RHESsian] [CHARge] [MICRo] [SAVE]
           [RESTart] [RESEt] (atom selection)

  REMOve:  Classical energies within QM atoms are removed.

  EXGRoup: QM/MM Electrostatics for link host groups removed.

  NOGUess: Obtains initial orbital guess from previous calculation.
           Default is to recalculate initial orbitals each time.

  BLURred: MM charges are treated as a gaussian function (equivalent to ECP)
           width of the gaussian function is specified by default in WMAIN
           array (usually by SCALar command). The value for charge is taken
           from PSF. Some values of WMAIN have special meaning:

           WMAIN.LE.   0.0 treat this atom as point charge in the QM/MM potential
           WMAIN.GE.5000.0 treat this atom as an infinitely diffuse Gaussian

  RECAll:  Use the RECAll array (as specified in scalar.doc) to set BLUR
           widths instead of the main WMAIN array. This is necessary when
           using Gaussian BLUR MM charges with the Replica Path or NEB
           methods as these make use of the WMAIN array. See QM-MM_DGMM.inp
           in the test directory for an example.

  COORdinates: This keyword will activate CHARMM to obtain an updated geometry
           from Q-Chem as the calculation proceeds. This can be particularly
           useful as Q-Chem can perform QM optimiztions using delocalized
           internal coordinates in the presence of a fix field of point charges.
           This can significantly speed QM/MM minimizations and can be used in
           an iterative approach. Note: to use this the JOBTYPE in the Q-Chem
           control file should be set to OPT (i.e. JOBTYPE = OPT).

  QCLJ:    Activates Q-Chem to use CHARMM's Lennard-Jones parameters when
           performing QM calculations in a fixed field of point charges. This
           can be particularlly useful as the QM region can be overly attracted
           to bare point charges.

  PARAllel: Allows the user to specify how many processors they wany the Q-Chem
           calculation to utilize. Previously, Q-Chem would use the same number
           of processors as CHARMM was using, however, in most cases the Q-Chem
           calculation will be much more expensive so having 1 CHARMM process
           and 4 Q-Chem processes is more efficient. Note: This currently only
           works with parallel versions of CHARMM although it can be extended
           to work with serial versions.

  SHESsian: Save Hessian computed via the QM/MM Normal Mode Procedure. The
           Hessian will be saved as an ascii file named: hessian.dat. Typically
           VIBRAN recomputes the Hessian each time it is needed; for QM or QM/MM
           calculations this is inefficient and thus saving the Hessian becomes
           very important. Particularly, this is used when employing the VSA
           method (see vibran.doc).

  RHESsian: Read a previously saved Hessian (hessian.dat) from a file (see SHES).

  CHARge:  Read QM charges from an file (charges.dat). A charges.dat file is
           created by Q-Chem by seting the REM keyword QMMM_CHARGES = TRUE.
           This file contains the Muliken charges for the QM region in the
           same order that is specified in the PSF file.

  MICRo:   Turns on the QM/MM Micro-iteration scheme of Kastner et al. J. Chem.
           Theory Comput., 3 (3), 1064-1072, 2007. Currently, this is best used
           in conjunction with a loop where CHARMM alternates between MM and
           QM/MM micro and macro cycles. Also, the CHARge keyword should be
           used to set the charges on QM region during the MM cycles. This is a
           new feature that requires further testing so be careful!

  RESEt:   Resets all QM/MM options to their initial defaults. This is needed
           for the QM/MM Micro-iteration approach to alternate between MM and
           QM/MM stages. Note: after using this a new "QCHEM" command must be
           issued!

  SAVE:    Activates CHARMM to save the converged SCF orbitals from a given
           energy calculation. Note: This ideally should be called once using
           a specific Q-Chem control file that contains specialized SCF
           convergence options. This option shoudl then be followed by a new
           "QCHEm" call that specifies "RESTart" and performs the actual
           QM/MM minimization. See below for example...

  RESTart: This tells CHARMM to restart a QM/MM calculation using previously
           saved orbitals. The "QCHEm" command that uses this as an option
           should be proceeded by a "QCHEm SAVE" command to preform an
           initial calculation saving the orbitals. See below for example...

           Example:

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            envi qchemcnt  "qcnt1.inp"    ! File that contains special SCF
            envi qcheminp  "q1.inp"       ! convergence options
            envi qchemexe  "qchem"
            envi qchemout  "q1.out"
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           QCHEm SAVE REMOve SELEct RESId 1 SHOW END
           ENERgy

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           envi qchemcnt  "qcnt2.inp"     ! Regular Q-Chem control file
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           QCHEm RESTart NOGUess REMOve SELEct RESId 1 SHOW END
           MINI ABNR NSTEp 10 NPRInt 1

  =======================================================================

           The atoms in selection will be treated as QM atoms.
           Link atom may be added between an QM and MM atoms with the
           following command:

  =======================================================================

  ADDLinkatom  link-atom-name  QM-atom-spec  MM-atom-spec

        link-atom-name ::= a four character descriptor starting with QQ.

        atom-spec::= {residue-number atom-name}
                     { segid  resid atom-name }
                     { BYNUm  atom-number     }

When using link atoms to break a bond between QM and MM
regions bond and angle parameters have to be added to parameter file
or better use READ PARAm APPEnd command.

If define is used for selection of QM region put it after all
ADDLink commands so the numbers of atoms in the selections are not
changed. Link atoms are always selected as QM atoms.


.. _qchem_usage:

Usage
-----

CHARMM input scripts are the same as before except the addition of ENVIronment
commands and the QCHEm command itself. Q-Chem commands are in a separate file
call qchem.inp, (or with an alternative name indicated by the "QCHEMCNT"
environment variable). The Q-Chem input file has the same structure as it
would have for a normal Q-Chem run, except that the specification of the
geometry, in the molecule section, is omitted. Note: the charge and
multiplicity are still included in the molecule section.

Names of the files for Q-Chem are specefied with environment
variables as follows. These four ENVIronment variables must be set!

::

     use ENVIronment command inside CHARMM

     ENVI qchemcnt  "qchem.inp"
     ENVI qcheminp  "q1.inp"
     ENVI qchemexe  "qchem"
     ENVI qchemout  "qchem.out"

or use the following for (t)csh

::

     setenv qchemcnt qchem.inp
     setenv qcheminp q1.inp
     setenv qchemexe qchem
     setenv qchemout qchem.out

or use the following for ksh,sh,bash

::

     export qchemcnt=qchem.inp
     export qcheminp=q1.inp
     export qchemexe=qchem
     export qchemout=qchem.out

1. The QCHEMCNT variable specifies the main Q-Chem input file which contains
   the $rem section, $molecule section (without geometry), $comment section,
   ect..,

2. The QCHEMINP variable is the final input file that will get passed to
   Q-Chem. CHARMM actually writes this file and adds the correct geometry and
   any external/point charges (e.g. MM atoms) to an $external_charges section.

3. The QCHEMEXE is the location of the qchem script. Specify the entire path
   unless $QC/bin is included in your default path.

4. The QCHEMOUT file specifies the Q-Chem output file. This file get
   overwritten for each optimization/time step. In the future, there will be a
   mechanism to save old output files.

Q-Chem input file parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following $rem variables must be specified in the QCHEMCNT file in order
to perform CHARMM QM/MM or pure QM calculations.

::

  qm_mm                 true
  jobtype               force
  symmetry              off
  sym_ignore            true
  print_input           false
  qmmm_print            true

1. qm_mm = true: Turns QM/MM on in Q-Chem

2. jobtype = force: Needed to do QM/MM optimizations. Set to "SP" if QM/MM
   energy is desired.

3. symmetry = off: Turn off symmetry

4. sym_ignore = true: Prevents Q-Chem from reorienting molecule

5. print_input = false: Use this if you have a large molecule and do not want
   1000s of atoms echoed back to the output file.

6. qmmm_print = true: Reduces some of the print out during QM/MM calculations.
   This prevents external charges from being printed out if
   there are more than 50 of them.

Sample QCHEMCNT file (qchem.inp):
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

  $comment
  Input file comes from CHARMM
  $end

  $rem
  exchange              HF
  basis                 6-31G*
  qm_mm                 true
  jobtype               force
  symmetry              off
  sym_ignore            true
  print_input           false
  qmmm_print            true
  $end

  $molecule
  0 1
  $end

The above is for 6-31G calculation of any neutral molecule.

[NOTE: For another example look at test/cquantumtest/alanine_qchem.inp]

.. _qchem_installation:

Installation
------------

One of the main benefits of using Q-Chem to do QM/MM calculations with CHARMM
is the ease of which you can get up and running jobs. All you have to do is
compile CHARMM in the following way....

::

  install.com <machine-type> <CHARMM size> QC <other CHARMM options>

This will compile the serial version of CHARMM to run with a serial version of
Q-Chem. To compile a parallel version of CHARMM to run with a parallel or
serial version of Q-Chem you could use the following script....

::

  #!/bin/csh
  # Compile Parallel CHARMM with Q-Chem support

  # USE STANDARD MPI (i.e. MPICH)
  setenv MPI /base/mpi/directory
  setenv MPI_LIB $MPI/lib
  setenv MPI_LIB $MPI/include

  # SET THE PATH TO MPIF77
  set path=($MPI/bin $path)

  install.com <machine-type> <CHARMM size> M QC MPICH <other CHARMM options>

.. _qchem_status:

Q-Chem/CHARMM interface status (July 2007)
------------------------------------------

- Parallel version is fully functional

- Replica/Path and Nudged Elastic Band Methods function in a highly parallel
  and parallel/parallel fashion.

- I/O including standard input and output are separated for
  Q-Chem.

- All CHARMM testcases are still OK when CHARMM is compiled
  with Q-Chem inside.

- QCHEM, GAMESS, GAMESSUK, CADPAC and QUANTUM keywords cannot coexist in
  pref.dat

- Q-Chem recognizes atoms by their masses as specified in the
  RTF file

- Delocalized Gaussian Blurred MM charges have been implemented for both
  energies and analytic gradients

- Full (i.e. no restraints/constraints) QM/MM 2nd derivatives (i.e. Hessians)
  are available.

.. _qchem_functionality:

Functionality
-------------

1. QM/MM optimizations (analytic gradients) using Q-Chem can be performed
   using the following methods.

   - HF*     (RHF, UHF, ROHF)
   - DFT*    (RHF, UHF, ROHF)
   - RIMP2   (RHF, UHF)
   - MOS-MP2 (RHF, UHF)
   - SOS-MP2 (RHF, UHF)
   - SCS-MP2 (RHF, UHF)
   - MP2     (RHF, UHF)
   - CCSD    (RHF, UHF)

   (*) Analytic derivatives run in parallel.

2. QM/MM single point energies using Q-Chem can be performed using the
   following methods (in addition to the above).

   - Local MP2 (RHF, UHF)
   - CCSD(T)   (RHF, UHF)

3. Additional analytic derivative and energy point methods will be made
   available in future releases. To request support for methods please contact
   H. Lee Woodcock (hlwoodr_at_nih_dot_gov) and/or post request to the CHARMM
   forums.

.. _qchem_rpath:

RPath
-----

1. Additional ENVIronment variable: To do QM/MM Replica/Path or Nudged Elastic
   Band calculations with CHARMM and Q-Chem you must define one extra variable.

   ::

      ENVI QCHEMPWD "/path/to/working/rpath/directory"

2. After defining this above ENVIronment variable all that is left to do is
   add the "rpath" keyword to the QCHEm call. For example...

   ::

      QCHEm RPATh REMOve select qm_region end

   This will create nrep directories in /path/to/working/rpath/directory and each
   point of the pathway will be computed in a different directory.

   Note: you must be running a parallel version of CHARMM with the same number of
   processors as you have replicas (i.e. pathway points).

.. _qchem_pert:

PERT
----

To run ab initio QM/MM free energy perturbation you need to specify additional
environment variables in the QM/MM setup...

1. sainp: state A control file (same as QCHEMCNT; specific for state A)
2. sbinp: state B control file (same as QCHEMCNT; specific for state B)
3. stateainp: auto generated Q-Chem input file for state A
4. statebinp: auto generated Q-Chem input file for state B
5. stateaout: specify Q-Chem output for state A QM calculation
6. statebout: specify Q-Chem output for state B QM calculation

Example...

::

   envi qchemexe  "qchem"               ! Command to call quantum program
   envi qchemcnt  "data/qchem_pert.inp" ! Non Pert Control file
   envi qcheminp  "q1.inp"              ! Non Pert Quantum input file
   envi qchemout  "qchem.out"           ! Non Pert Quantum output file
   envi sainp     "data/s0.inp"         ! State 0 control file
   envi sbinp     "data/s1.inp"         ! State 1 control file
   envi stateainp "state0.inp"          ! State 0 quantum input file
   envi statebinp "state1.inp"          ! State 1 quantum input file
   envi stateaout "state0.out"          ! State 0 quantum output file
   envi statebout "state1.out"          ! State 1 quantum output file

See test/cquantumtest/qmmm_pert.inp for a complete example.

Please see :doc:`pert` for a complete description of running free energy
perturbation in CHARMM.

.. _qchem_normal_mode_analysis:

Normal Mode Analysis
--------------------

To run full QM/MM Normal Mode Analysis (i.e. QM/MM 2nd derivatives, Hessians)
you need to run QM/MM with the VIBRan module (see :doc:`vibran`) of CHARMM. To
perform this calculation just run QCHEm has usual...

Example:

::

  QCHEm REMOve SELEct RESId 1 SHOW END

Then invoke the VIBRan module...

::

  VIBRan
  DIAG
  END

In addition, you must add the following line to the QCHEMCNT file (the file that
controls the the REM variables passed to Q-Chem).

::

  QMMM_FULL_HESSIAN     TRUE

Please see the "QM-MM_Normal_Modes.inp" testcase in the test directory for the
full example.

Currently, this only works with standard point charge QM/MM models (i.e. not
Gaussian blurred charges), but this will be extended in the future.

