.. py:module:: subst

====================================
Command Line Substitution Parameters
====================================

The following are substitution parameters available within CHARMM;

General
-------

=============   =========================================================================================
:sub:`PI`       Pi, 3.141592653589793
:sub:`KBLZ`     The Boltzmann factor (0.001987191)
:sub:`CCELEC`   1/(4 PI epsilon) in AKMA units (332.0716)
:sub:`SPEEDL`   Speed of light
:sub:`CNVFRQ`   Conversion from root(Kcals/mol/AMU) to frequencies in CM-1
:sub:`TIMFAC`   Conversion from AKMA time to picoseconds
=============   =========================================================================================

Control and system variables
----------------------------

==============   ========================================================================================
:sub:`BOMLEV`    The error termination level (-5 to 5)
:sub:`WRNLEV`    The warning print level (-5 to 10)
:sub:`PRNLEV`    The standard print level (-1 to 15)
:sub:`IOLEV`     The I/O level (-1 to 1)
:sub:`IOSTAT`    The status of most recent OPEN command (-1=failed,1=OK)
:sub:`TIMER`
:sub:`FASTER`
:sub:`LFAST`
:sub:`OLMACH`
:sub:`OUTU`
:sub:`FLUSH`
:sub:`FNBL`
:sub:`NBFACT`
:sub:`LMACH`
:sub:`MYNODE`    Current node number (0 to NUMNODE-1)
:sub:`NUMNODE`   The number of nodes (distributed memory)
:sub:`NCPU`      The number of CPUs (shared memory use)
:sub:`SYSSTAT`
:sub:`CPUTIME`   CPU time used at last call to TIMRE (TIMER NOW, TIMER DIFF)
:sub:`ELATIME`   Elapsed time  at last call to TIMRE (TIMER NOW, TIMER DIFF)
:sub:`PID`       Process ID of the CHARMM process (or the rank 0 process when
                 MPI is used).
==============   ========================================================================================

PSF counts
----------

==============   ========================================================================================
:sub:`NSEG`      Number of segments
:sub:`NRES`      Number of residues
:sub:`NATOM`     Number of atoms
:sub:`NGRP`      Number of groups
:sub:`NBOND`     Number of bonds
:sub:`NTHETA`    Number of angles
:sub:`NPHI`      Number of dihedrals
:sub:`NIMPHI`    Number of improper dihedrals
:sub:`NACC`      Number of acceptors
:sub:`NDON`      Number of donors
:sub:`NNB`       Number of explicit nonbond exclusions
:sub:`CGTOT`     Total system charge
:sub:`MASST`     Total system mass
:sub:`NATI`      Total number of image plus primary atoms
==============   ========================================================================================

Parameter counts
----------------

==============   ========================================================================================
:sub:`NATC`      Number of atom types
:sub:`NCB`       Number of bond parameters
:sub:`NCT`       Number of angle parameters
:sub:`NCSB`      Number of stretch-bend parameters
:sub:`NCP`       Number of diheral parameters
:sub:`NCI`       Number of improper dihedral parameters
:sub:`NCOOP`     Number of out-of-plane parameters
:sub:`NCH`       Number of hydrogen bond parameters
:sub:`NCN`       Number of vdw parameter pairs
:sub:`NCQ`       Number of bond charge increments
==============   ========================================================================================

Other counts
------------

==============   ========================================================================================
:sub:`NCSP`      Number of restrained dihedral (CONS DIHE command).
:sub:`NTRA`      Number of image transformations
:sub:`TOTK`      Number of Ewald K vectors (not PME)
:sub:`NIC`       Number of Internal Coordinate entries in the IC table
==============   ========================================================================================

Dimension Limits
----------------

==============   ========================================================================================
:sub:`MAXA`      Number of atoms
:sub:`MAXATC`    Number of atom types
:sub:`MAXB`      Number of bonds
:sub:`MAXIMP`    Number of improper dihedrals
:sub:`MAXNB`     Number of explicit nonbond exclusions
:sub:`MAXP`      Number of dihedrals
:sub:`MAXPAD`    Number of donors and acceptors
:sub:`MAXRES`    Number of residues
:sub:`MAXSEG`    Number of segments
:sub:`MAXT`      Number of angles
:sub:`MAXCB`     Number of bond parameters
:sub:`MAXCH`     Number of hydrogen bond parameters
:sub:`MAXCI`     Number of improper dihedral parameters
:sub:`MAXCN`     Number of vdw pair parameters
:sub:`MAXCP`     Number of dihedral parameters
:sub:`MAXCT`     Number of angle parameters
:sub:`MAXCSP`    Number of restrained dihedrals
==============   ========================================================================================

Coordinate manipulation parameters
----------------------------------

===============   =======================================================================================
:sub:`XAXI`       vector of defined axis (set by the COOR { AXIS | ORIE RMS | LSQP | HELIx } command).
:sub:`YAXI`
:sub:`ZAXI`
:sub:`RAXI`       length of vector (often 1.0)
:sub:`XCEN`       origin/center of axis vector
:sub:`YCEN`
:sub:`ZCEN`
:sub:`XMIN`       Extreme values (COOR STAT command)
:sub:`YMIN`
:sub:`ZMIN`
:sub:`WMIN`
:sub:`XMAX`
:sub:`YMAX`
:sub:`ZMAX`
:sub:`WMAX`
:sub:`XAVE`       Average values (COOR STAT command).
:sub:`YAVE`
:sub:`ZAVE`
:sub:`WAVE`
:sub:`MASS`       mass of selected atoms
:sub:`RMS`        Root mean squared difference between two structures.
:sub:`XMOV`       displacement of atoms from best fit (COOR ORIE command).
:sub:`YMOV`
:sub:`ZMOV`
:sub:`THET`       Angle of rotation from best fit (degrees)
:sub:`SHIFT`      Translation of best fit move projected on rotation axis.
:sub:`AREA`       Requested surface area (COOR SURF command).
:sub:`VOLUME`     Requested volume (COOR VOLUme command).
:sub:`NVAC`       Number of vacuum points
:sub:`NOCC`       Number of occupied points
:sub:`NSEL`       Number of selected points
:sub:`FREEVOL`    Total free volume
:sub:`MIND`       Minimum distance (COOR MIND command).
:sub:`NPAIR`      Number of pairs (COOR DIST command).
:sub:`NCONTACT`   Number of contacts (COOR DMAT command).
:sub:`RGYR`       Radius of gyration (COOR RGYR command).
:sub:`XCM`        Center of mass (COOR RGYR command).
:sub:`YCM`
:sub:`ZCM`
:sub:`XDIP`       Dipole moment  (COOR DIPOle command)
:sub:`YDIP`
:sub:`ZDIP`
:sub:`RDIP`       Dipole magnitude
:sub:`CHARGE`     Charge of selected atoms
:sub:`NHBOND`     total number of hydrogen bonds (COOR HBONd command).
:sub:`AVNOHB`     Average number of hydrogen bonds
:sub:`AVHBLF`     Average hydrogen bond life
:sub:`MINDA1`     First  atom of minimum distance atom pair
:sub:`MINDA2`     Second atom of minimum distance atom pair
:sub:`NHYDRR`     Average number of solvent molecule -solute contacts (COOR ANALys)
:sub:`NHYDAR`     Average number of solvent atom  - solute contacts (COOR ANALys)
:sub:`NHYDAA`     Average number of solvent atom -solute atom contacts (COOR ANALys)
:sub:`ENTROPY`    Configurational entropy estimate (COOR COVA THERmo)
:sub:`SROT`       Rotational entropy    (COOR INERtia ENTRopy)
:sub:`STRA`       Translational entropy (COOR INERtia ENTRopy)
:sub:`SVIB`       Vibrational entropy   (VIBRAN DIAG ENTRopy )
:sub:`SSUM`       SROT+STRA+SVIB
:sub:`NALPHA`     Number of residues in alpha helix
:sub:`ALPHA`      Fraction of residues in alpha helix
:sub:`NBETA`      Number of residues in beta strands
:sub:`BETA`       Fraction of residues in beta strands
:sub:`PHASE`      Sugar pucker phase
:sub:`AMP`        Sugar pucker amplitude
===============   =======================================================================================


SCALar STATistics command substitution parameters
-------------------------------------------------

==============   ========================================================================================
:sub:`SMIN`      Minimum value
:sub:`SMAX`      Maximum value
:sub:`SAVE`      Average value
:sub:`SVAR`      Variance about average
:sub:`SWEI`      Total weight used in the averaging
:sub:`STOT`      Total of selected atoms
:sub:`NSEL`      Number of selected atoms
==============   ========================================================================================


Quick command substitution parameters
-------------------------------------

==============   ========================================================================================
:sub:`XVAL`      X position of group of atoms
:sub:`YVAL`      X position of group of atoms
:sub:`ZVAL`      X position of group of atoms
:sub:`DIST`      Distance between two atom analysis
:sub:`THET`      Angle for three atom analysis
:sub:`PHI`       Dihedral for four atom analysis
==============   ========================================================================================


Shape analysis
--------------

==============   ========================================================================================
:sub:`SFIT`
:sub:`THET`
:sub:`XAXI`
:sub:`YAXI`
:sub:`ZAXI`
:sub:`RAXI`
==============   ========================================================================================


Saddle point calculation (TRAVel)
---------------------------------

==============   ========================================================================================
:sub:`SADE`      Saddle point energy
:sub:`SADI`      Saddle point index
:sub:`SADO`      Saddle point order
==============   ========================================================================================


Energy calculation results
--------------------------

==============   ========================================================================================
:sub:`XCM`       Center of mass (from MMFP energy term calcuation)
:sub:`YCM`
:sub:`ZCM`
:sub:`XCM2`      Spatial extent
:sub:`YCM2`
:sub:`ZCM2`
:sub:`RGEO`      average distance from reference
:sub:`ENPB`      electrostatic free energy of solvation (from PBEQ)
:sub:`RMAX`      maximum distance to origin for the SSBP energy term
==============   ========================================================================================


Minimization results
--------------------

================   ======================================================================================
:sub:`MINCONVRG`
:sub:`MINECALLS`
:sub:`MINGRMS`
:sub:`MINSTEPS`
================   ======================================================================================


PERT results
------------

==============   ========================================================================================
:sub:`TPDEL`     Thermodynamic Perturbation energy change
:sub:`TPTOT`     Thermodynamic Perturbation total energy
:sub:`TIDEL`     Thermodynamic Integration energy change
:sub:`TITOT`     Thermodynamic Integration total energy
:sub:`SLDEL`     Slow Growth energy change
:sub:`SLTOT`     Slow Growth total energy
:sub:`DFLC`      DIFFLC, ie the fluctuation about the average energy difference
:sub:`AVKE`      Average kinetic energy
:sub:`KEFL`      Fluctuation in kinetic energy
==============   ========================================================================================


Atom selection parameters
-------------------------

==============   ========================================================================================
:sub:`NSEL`      Number of selected atoms from the most recent atom selection.
:sub:`SELATOM`   Atom number of first selected atom
:sub:`SELCHEM`   Chemical type of first selected atom
:sub:`SELIRES`   Residue number of first selected atom
:sub:`SELISEG`   Segment number of first selected atom
:sub:`SELRESI`   Resid of first selected atom
:sub:`SELRESN`   Residue type of first selected atom
:sub:`SELSEGI`   Segid of first selected atom
:sub:`SELTYPE`   Atom name of first selected atom
==============   ========================================================================================


Crystal parameters
------------------

===============   =======================================================================================
:sub:`XTLA`       Unit cell dimensions
:sub:`XTLB`
:sub:`XTLC`
:sub:`XTLALPHA`   Unit cell angles
:sub:`XTLBETA`
:sub:`XTLGAMMA`
:sub:`XTLXDIM`    Number of crystal degrees of freedom (cube=1,triclinic=6,..)
===============   =======================================================================================


Data from most recently read (or current) trajectory file
---------------------------------------------------------

==============   ========================================================================================
:sub:`NFILE`     Number of frames in the trajectory file
:sub:`START`     Step number for the first frame
:sub:`SKIP`      Frequency at which frames were saved
                 (NSTEP=NFILE*SKIP when not using restart files)
:sub:`NSTEP`     Total number of steps in the simulation
:sub:`NDEGF`     Number of degrees of freedom in the simulation
                 (Can be use to get the temperature with velocity files).
:sub:`DELTA`     The dynamics step length (in picoseconds).
:sub:`NTOT`      Total number of frames actually read
==============   ========================================================================================


Nonbond list counts
-------------------

==============   ========================================================================================
:sub:`NNBA`      Number of atom  pairs (main list)
:sub:`NNBG`      Number of group pairs (main list)
:sub:`NNBI`      Number of crystal atom pairs (Phonons only)
:sub:`NRXA`      Number of atom  exclusions due to replicas
:sub:`NRXG`      Number of group exclusions due to replicas
==============   ========================================================================================


Correlation Function Results
----------------------------

==============   ========================================================================================
:sub:`AVER`      Series average (CORREL's SHOW command)
:sub:`FLUC`      Series fluctuation
:sub:`P2`        P2 average
:sub:`P2R3`
:sub:`P2RA`
:sub:`R3R`
:sub:`R3S`
:sub:`P0`        Polynomial best fit components (MANTime POLY command).
:sub:`P1`
:sub:`P2`
:sub:`P3`
:sub:`P4`
:sub:`P5`
:sub:`P6`
:sub:`P7`
:sub:`CFNORM`    Multiplicative normalization factor for correlation function
==============   ========================================================================================


Generalized Born Solvation Calculation Results
----------------------------------------------

==============   ========================================================================================
:sub:`GBAL`      Generalized Born alpha values
:sub:`GBAT`      Generalized Born atomic solvation energy contributions values
:sub:`SIGX`      Partial contribution to GB force, X-direction (diagnostic)
:sub:`SIGY`      Partial contribution to GB force, Y-direction (diagnostic)
:sub:`SIGZ`      Partial contribution to GB force, Z-direction (diagnostic)
:sub:`T_GB`      Partial contribution to GB force, 3-body term (diagnostic)
==============   ========================================================================================


Vibrational analysis of thermodynamic properties
------------------------------------------------

==============   ========================================================================================
:sub:`FTOT`      Vibrational free energy.
:sub:`STOT`      Vibrational entropy.
:sub:`HTOT`      Vibrational enthalpy.
:sub:`CTOT`      Vibrational heat capacity.
:sub:`ZTOT`      Zero point correction energy.
:sub:`FCTO`      Classical vibrational free energy.
:sub:`ETOT`      Total harmonic limit classical free energy
                 (to compare with free energy perturbation simulations).
:sub:`TRAC`      Trace of the Hessian for selected atoms
:sub:`SROT`      Rotational entropy    (COOR INERtia ENTRopy)
:sub:`STRA`      Translational entropy (COOR INERtia ENTRopy)
:sub:`SVIB`      Vibrational entropy   (VIBRAN DIAG ENTRopy )
:sub:`SSUM`      SROT+STRA+SVIB
==============   ========================================================================================

Multi-Site lambda-dynamics trajectory analysis
----------------------------------------------

==============   ========================================================================================
:sub:`TMIN`      Minimum number of transitions for any site in the system
:sub:`TMAX`      Maximum number of transitions for any site in the system
:sub:`FPL`       Fraction of the snapshots which represent full Physical Ligands
:sub:`POP#`      Population for the substituent associated with indicated BLOCK number at the low
                 threshold value (e.g. the ?pop2 contains the population for substituent in BLOCK 2
                 given CUTLO threshold)
:sub:`DDG#_#`    Relative free energy between the first and second substituents listed at the low
                 threshold value (e.g. ?ddg2_5 is the relative free energy between the substituents
                 associated with BLOCKS 2 and 5).
==============   ========================================================================================

Miscellaneous
-------------

===============  ========================================================================================
:sub:`VIOL`      Total violation for all NOE restraints (NOE WRITe/PRINt ANAL)
:sub:`DRSH`      the DRSH value in subroutine PSHEL (undocumented)
                 (also undocumented in fcm/mmfp.fcm in violation of coding stds.)
:sub:`DCOEFF`    The diffusion constant (COOR ANALysis SOLVent command).
:sub:`TIME`      simulation time(ps) for current frame in trajectory reading
:sub:`STEP`      Step number for current frame in trajectory reading
:sub:`PQRES`     Final value of target function in RMSDyn 2D-projection
:sub:`WHAMFE`    total free energy of WHAM
:sub:`SSBPLRC`   long-range free energy correction for SSBP
:sub:`SSBPLRCS`  standard deviation of SSBP long-range correction
===============  ========================================================================================



pref.dat Keywords
-----------------

The following keywords have substitutions of 0 or 1 for set or not-set,
respectively.

::

   ACE       FOURD     MOLVIB  	OLDDYN    PNOE     RXNCOR
   ADUMB     GAMESS    MPI     	PARAFULL  POLAR    SCALAR
   ASPENER   GENBORN   MTS     	PARALLEL  PRIMSH   SHAPES
   BLOCK     GENCOMM   MULTCAN 	PARASCAL  PVM      SINGLE
   CFF       GENETIC   NIH     	PARVECT   PVMC     SOCKET
   CMPI      IMCUBES   NOCORREL	PATHINT   QHMCM    SOFTVDW
   CRAYVEC   LATTICE   NOIMAGES	PBEQ      QTSM     TNPACK
   DIMB      LDM       NOMISC  	PBOUND    QUANTA   TRAVEL
   DMCONS    LONEPAIR  NOST2   	PBOUNDC   QUANTUM  VECTOR
   DOCK      MC        NOVIBRAN	PERT      REPLICA  WCA
   EISPACK   MCSS      NO_BYCC 	PM1       RGYCONS
   FMA       MMFF      NO_BYCU 	PMEPLSM   RISM

Example: In this case GENBORN keyword was in the pref.dat so
:sub:`?genborn` is substituted with 1 and the IF test will
evaluate as true.

::

 	  if ?genborn .eq. 1 then goto dogenborn

	   CHARMM>    if ?genborn .eq. 1 then goto dogenborn
 	   RDCMND substituted energy or value "?GENBORN" to "1"
 	   Comparing "1" and "1".
 	   IF test evaluated as true.  Performing command


See :doc:`energy` for the energy related substitution parameters.
