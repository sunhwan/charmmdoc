.. py:module:: energy

===============================================
Energy Manipulations: Minimization and Dynamics
===============================================

The main purpose of CHARMM is the evaluation and manipulation of
the potential energy of a macromolecular system. In order to compute
the energy, several conditions must be met. There are also several
support commands which directly relate to energy evaluation.

.. index:: energy; description
.. _energy_description:

Syntax for Energy Commands
--------------------------

There are two direct energy evaluation commands. One is parsed
through the minimization parser and the other involves a direct call
to GETE.  See :doc:`minimiz` and
:ref:`interface <usage_gete>`.  In addition to getting the energy,
the forces are also obtained.


The ENERgy command
^^^^^^^^^^^^^^^^^^

The energy command is processed through the minimization parser

::

   ENERgy [ nonbond-spec ] [ hbond-spec ] [ image-spec ] [ print-spec ] [ COMP ]
          [ domdec-spec ] [ COMP ] [  INBFrq 0    ] [  IHBFrq 0  ]
          [  IMGFrq 0  ] [NOUPdate] [ openmm-spec ]

   hbond-spec        *note Hbonds:(chmdoc/hbonds.doc).
   nonbond-spec      *note Nbonds:(chmdoc/nbonds.doc).
   image-spec        *note Images:(chmdoc/images.doc)Update.
   domdec-spec       *note Domdec:(chmdoc/domdec.doc).
   openmm-spec       *note OpenMM:(chmdoc/openmm.doc).

If the COMP keyword is specified, then the comparison coordinate
set is used, but this disables the use of the fast routines. The keyword
NOUPdate turns off all update routines, and thus requires all lists
to be present already.


The GETE command
^^^^^^^^^^^^^^^^

The command is a direct call to GETE routine.

::

  GETE  [ COMP ] [ PRINt [ UNIT int ] ] [ openmm-spec ]
                 [ NOPRint            ]

For this command to work, all list must be set up. This is best done
through the UPDAte command. The COMP keyword will cause the comparison
coordinate set to be used. The PRINt keyword will result in a subsequent
call to PRINTE in order to print the energy. If the PRINt keyword is not
specified, then NO indication that the energy has been called will be given.


The UPDAte command
^^^^^^^^^^^^^^^^^^

The command sets up required lists for GETE.

::

   UPDAte [ nonbond-spec ] [ hbond-spec ] [ image-spec ] [ COMP ]
          [  INBFrq 0    ] [  IHBFrq 0  ] [  IMGfrq 0  ]
          [  EXSG {list-of-segment-names} | EXOF ]

The update command will set up the codes lists and also create a
nonbond list (unless INBFrq is 0) and a new hbond list (unless IHBFrq is 0).
If the COMP keyword is specified, then the comparison coordinates will be
used in setting up the nonbond and hbond lists.

EXSG keword with optional following list of segment names allows to
exclude some nonbonded interactions (ELEC & VDW). If list of names is empty
ALL INTERsegment nonbonded interactions will be excluded. If list is not
empty all INTER and INTRA segment nonbonded interactions for listed
segments will be ecluded. EXOF turns off this option.
H-bond energies (HBON) are not affected at the moment (Dec 3, 1991).


.. index:: energy; skipe
.. _energy_skipe:

Skipping selected energy terms
------------------------------

There is a facility to skip any desired energy terms during
energy evaluation. For each energy term there is associated a logical
flag determining whether that energy term is to be computed.

Specifications are processed sequentially. The default operation
is INCLude which implies that subsequent energy term are to be removed
from the energy calculation.

.. note::

   EXCLude implies that the
   energy term is to be computed.

If for some reason, the list presented here is out of date, the
data in SKIPE(energy.src) and in ENER.FCM of the source should be
consulted.

Syntax:

::

                   [ INCLude ]
                   [ EXCLude ]
   SKIPe  repeat(  [   ALL   ]  )
                   [   NONE  ]
                   [   item  ]

   item::=
             [ BOND ]   [ ANGL ]  [ UREY ]   [ DIHE ]
             [ IMPR ]   [ VDW  ]  [ ELEC ]   [ HBON ]
             [ USER ]   [ HARM ]  [ CDIH ]   [ CIC  ]
             [ CDRO ]   [ NOE  ]  [ SBOU ]   [ IMNB ]
             [ IMEL ]   [ IMHB ]  [ XTLV ]   [ XTLE ]
             [ EXTE ]   [ RXNF ]  [ ST2  ]   [ IMST ]
             [ TSM  ]   [ QMEL ]  [ QMVDW]   [ ASP  ]
             [ EHARM]   [ GEO  ]  [ MDIP ]   [ STRB ]
             [ VATT ]   [ VREP ]  [ IMVREP ] [IMVATT]
             [ OOPL ]   [ CMAP ]  [ EPOL ]   [ CPUC ]

description:

  =======  ===================================================================
  BOND     bond energy
  ANGL     angle energy
  UREY     Urey-Bradley energy term
  DIHE     dihedral energy
  IMPR     improper dihedral energy
  VDW      van der Waal energy
  ELEC     electrostatic energy
  HBON     hydrogen bond energy
  USER     user supplied energy (USERLINK)
  HARM     harmonic positional constraint energy
  CDIH     constrained dihedral energy
  CPUC   - constrained puckering energy
  CIC      internal coordinate constraint energy
  CDRO     quartic droplet potential energy
  NOE      NOE general distance restraints
  SBOU     solvent boundary energy
  IMNB     image van der Waal energy
  IMEL     image electrostatic energy
  IMHB     image hydrogen bond energy
  XTLV     crystal van der Waal energy
  XTLE     crystal electrostatic energy
  EXTE     extended electrostatic energy
  RXNF     reaction field energy
  ST2      ST2 water-water energy
  IMST     image ST2 water-water energy
  TSM      TMS free energy term.
  QMEL     energy for the quantum mechanical atoms and their
           electrostatic interactions with the MM atoms using the AM1
           or MNDO semi-empirical approximations
  QMVDW    van der Waals energy between the quantum mechanical and
           molecular mechanical atoms
  ASP      solvation free energy term based on Wesson and Eisenberg
           surface area method
  EHARM    second harmonic restraint term (for implicit Euler integration)
  GEO      Mean-Field-Potential energy
  MDIP     MDIPole mean fields constraints
  STRB     strech-bend interaction (MMFF)
  VATT     VdW attraction (MMFF)
  VREP     VdW repulsion (MMFF)
  IMVREP   image VdW repulsion (MMFF)
  IMVATT   image VdW attraction (MMFF)
  OOPL     out-of-plane (MMFF)
  CMAP     2D dihedral cross term energy correction map
  EPOL     polarization energy computed from PIPF (see pipf.doc)
  =======  ===================================================================

Examples:

::

  SKIP ALL EXCL BOND - do just bond energy
  SKIP EXCL ALL      - return flags to default state
  SKIP ELEC VDW      - throw out electrostatics and van der Waals energy


.. index:: energy; interaction
.. _energy_interaction:

Interaction energies and forces
-------------------------------

The INTEraction command computes the energy and forces
between any two selections of atoms.

::

   INTEraction [ COMP ] [ NOPRint ] 2x(atom-selection) [UNIT int]

If only one atom selection is given, then a self energy will be computed.
This routine is quite efficient and may be used within a CHARMM loop
without too much overhead, though there are some restrictions.
The COMP keyword causes the comparison coordinates to be used.
The NOPRint keyword will prevent the results from being printed.

This routine works in the same manner as the GETE command in that
all of the lists (CODES, nonbond, and Hbond) must be specified before
invoking this command. One difference is that SHAKE will not be respected
with this command (i.e. if the coordinates don't satisfy the constraints,
neither will the energy).

The following energy terms may be computed by this routine
(unless suppressed with the SKIP command);

===============   =========================================================
Bond              Energy defined by the two atoms involved.
Angles            Energy allocated to the central atom (auto energy only).
Dihedral          Energy defined between central two atoms
Improper          Energy defined by first atom (auto energy only)
van der Waal      ATOM option only. Energy defined by two atoms involved.
Electrostatic     ATOM option only. Energy defined by two atoms involved.
Hbond             Energy defined by heavy atom donor and acceptor atom.
Harmonic cons     Energy allocated to central atom (auto energy only).
Dihedral cons     Energy defined by central two atoms.
User energy       Atom selections may be passed to USERE in the selection
                  common (DEFIne command).
                  Fill forces and energies as desired.
===============   =========================================================

All other energy terms will be zeroed. For terms listed "auto energy only",
the corresponding atom must be present in both atom selections.
For the remaining terms, one atom of the pair must be present in each
of the atom selections. The energy division matches the method used in
the analysis facility.

This command will not work with the selection of images atoms,
or the selection of ST2 waters. All energy terms not listed above will
not be computed. The nonbond list must be generated with the ATOM and VATOM
options. [T.Lazaridis, July 1999: Now INTE can work with the GROUP option]

The individual energy terms are stored in the energy common
and are available in commands and titles via the "?energy-term"
substitution.

The forces for all kept energy terms will be returned in
the force arrays. Note, that it is possible for atoms to have a force
that were not selected in either selection specification. This may
happen for angle or dihedral terms on the first and last atoms. It may
also happen in a similar manner for improper dihedrals, hydrogen bonding
terms, and dihedral constraints.


.. index:: energy; eten
.. _energy_eten:

The 10-12 van der Waals potential
---------------------------------

The ETEN command is used to switch between the use of a 6-12 van der
Waals potential (default), and a 10-12 potential.

::

   ETEN {ON}
        {OFF}

Setting the flag "ETEN" to "ON or OFF" switches the van der Waals to
a modified Lennard-Jones function containing an attractive r^-10 term and
repulsive r^-12 and r^-6 terms. This was introduced to support simulation of
the Go models built by the webserver at
http://mmtsb.scripps.edu/webservices/gomodel.html

When the 10-12 potential is turned on, all energy evaluations will be
carried out using this potential, including minimizations, normal mode
analysis, etc. Issuing the ETEN command with any keyword other than ON will
turn off the 10-12 potential, reverting to the 6-12 potential.

The 10-12 potential energy may be turned off without reverting to the
6-12 potential using the SKIPE command with the VDW item, since this potential
replaces the VDW energy.

This option does not support CFF, MMFF, IMAGE, GRAPE, ewald, multi-
body dynamics, and fast vector. It also does not does not support van der
Waals shifting, force switching, or switching, as well as soft core van der
Waals.

This option fully supports BLOCK.


.. index:: energy; fast
.. _energy_fast:

FASTER
------

::

   FASTer {integer}
          {OFF    }
          {ON     }
          {DEFAult}
          {SCALar } ! for testing only
          {VECTor } ! for testing only
          {CRAYvec } ! Use parallel code designed for a CRAY
          {PARVec  } ! Use parallel/vector code best SMP machines and Convex

Instead of using an integer value, FASTer command can be issued
with one of the following keywords. The FASTer keyword or integer defines
which versions of the energy routines to be used.

+------+----------------+--------------------+-----------------------------------------------------+
|      |  Keyword       | Equivalent integer | Description                                         |
+======+================+====================+=====================================================+
|FASTer|  OFF           |    -1              | Always use slow routines                            |
+------+----------------+--------------------+-----------------------------------------------------+
|      |  DEFAult       |    0               | Use fast routine if possible, no error              |
|      |                |                    | if cannot (default)                                 |
+------+----------------+--------------------+-----------------------------------------------------+
|      |  ON            |    1               | Use best optimized routine for the current machine  |
|      |                |                    | (Error message if cannot)                           |
+------+----------------+--------------------+-----------------------------------------------------+
|      |  SCALar        |    2               | Use fast scalar routine (Error message if cannot)   |
+------+----------------+--------------------+-----------------------------------------------------+

There exist a general and a fast version of the internal
energy routines (bond, angle, dihedral, and improper dihedral).  The
is also a fast version of nonbond energy evaluation (roughly 30-50%
faster).  These routines were designed for long minimization or
dynamics calculations.

To request the FAST routine, the FASTer command should be used
with a positive integer or an appropriate  keyword.  A negative
integer will disable the fast energy routines.  If the fast routines
are requested and it is not possible to use the fast routines, a
warning will be issued, and the general routines will be used in their
place.

The fast routines are more efficient in several ways;

1) arrays are included in common files rather than passed
2) second derivatives have been removed
3) analysis and print options have been removed

The restrictions are that;
1) the MAIN coordinate set must be used in the energy evaluations
2) second derivatives may not be requested
3) The PSF, parameter, and codes arrays must be used (from the common files)
4) a limited set of nonbond options must be used.

The current nonbond options supported by the fast nonbond routine
are as follows.

::

         ATOM [CDIE] [SHIFt  ]  VATOM [VSHIft  ]
              [RDIE] [SWITch ]        [VSWItch ]
                     [FSWItch]        [VFSWitch]
                     [FSHIft ]

        GROUP [CDIE] [SWITch ]  VGROUP [VSWItch ]
              [RDIE] [FSWItch]


.. _energy_needs:

Requirements before energy manipulations can take place
-------------------------------------------------------

Before the energy of a system can be evaluated and manipulated,
a number of data structures must be present.

First, a PSF must be present.

Second, a parameter set must be present. It must contain all
parameters which are required by the PSF being used.

Third, coordinates must be defined for every atom in the system.
An undefined coordinate has a particular value, and if two coordinates
have the same value, division by zero will occur in the evaluation of
the energy. If the positions of hydrogens are required, the hydrogen
bond generation routine, see :doc:`hbonds`, must be
called before the energy is evaluated.

Fourth, provisions must be made for having a hydrogen bond list
and a non-bonded interaction list. Having non-zero frequencies for
updating this lists is one way, one can also read these lists in, see
:ref:`read <io_read>`, or generate them with separate
commands, see :doc:`HBgen <hbonds>`, or
:doc:`NBgen <nbonds>`.


.. _energy_optional:

Optional actions you can take to modify the energy manipulations
----------------------------------------------------------------

There exist several commands which can modify the way the
potential energy is calculated or can affect the way energy
manipulations are performed.

The Constraint command, see :doc:`cons`, can
be used to constraints of various kinds. First, it can be used to set
flags for particular atoms which will prevent them from being moved
during minimization or dynamics. Second, it can be used to add
positional constraint term to the potential energy. This term will be
harmonic about some reference position. The user is free to set the
force constant. Third, the user can place a harmonic constraint on the
value of particular torsion angles in an attempt to force the geometry
of a molecule. Other constraints are also available.

The SHAKe command, see :ref:`SHAKE <cons_shake>`, is
used to set constraints on bond lengths and also bond angles during
dynamics. It is very valuable in that it permits a larger step size to
be used during dynamics. This is vital for dynamics where hydrogens
are explicitly represented as the low mass and high force constant of
bonds involving hydrogen require a ridiculously small step size.

The user interface commands can be used to modify the
calculation of the potential and to add another term to the potential
energy.  See :ref:`interface <usage_modify>` for details.


.. _energy_substitution:

Substitution
------------

The following command line substitution values may be included in
any command or title.  To get the total energy, the syntax;

::

      ...... ?TOTE .....

should be used.

Energy related properties:

=========  =============================================================
 ?TOTE     total energy
 ?TOTK     total kinetic energy
 ?ENER     total potential energy
 ?TEMP     temperature (from KE)
 ?GRMS     rms gradient
 ?BPRE     boundary pressure applied
 ?VTOT     total verlet energy (no HFC)
 ?VKIN     total verlet kinetic energy (no HFC)
 ?EHFC     high frequency correction energy
 ?EHYS     slow growth hysteresis energy correction
 ?VOLU     the volume of the primitive unit cell
           = A.(B x C)/XNSYMM. Defined only if images are present,
           or unless specified with the VOLUme keyword.
 ?PRSE     the pressure calculated from the external virial.
 ?PRSI     the pressure calculated from the internal virial.
 ?VIRE     the external virial.
 ?VIRI     the internal virial.
 ?VIRK     the virial "kinetic energy".
=========  =============================================================

Energy term names:

=========  =============================================================
 ?BOND     bond (1-2) energy
 ?ANGL     angle (1-3) energy
 ?UREY     additional 1-3 urey bradley energy
 ?DIHE     dihedral 1-4 energy
 ?IMPR     improper planar of chiral energy
 ?CMAP     2D dihedral cross term energy correction map
 ?STRB     Strech-Bend coupling energy (MMFF)
 ?OOPL     Out-off-plane energy (MMFF)
 ?VDW      van der waal energy
 ?ELEC     electrostatic energy
 ?HBON     hydrogen bonding energy
 ?USER     user supplied energy term
 ?HARM     harmonic positional restraint energy
 ?CDIH     dihedral restraint energy
 ?CPUC     puckering restraint energy
 ?CIC      internal coordinate restraint energy
 ?CDRO     droplet restraint energy (approx const press)
 ?NOE      general distance restraint energy (for NOE)
 ?SBOU     solvent boundary lookup table energy
 ?IMNB     primary-image van der waal energy
 ?IMEL     primary-image electrostatic energy
 ?IMHB     primary-image hydrogen bond energy
 ?EXTE     extended electrostatic energy
 ?EWKS     Ewald k-space sum energy term
 ?EWSE     Ewald self energy term
 ?RXNF     reaction field electrostatic energy
 ?ST2      ST2 water-water energy
 ?IMST     primary-image ST2 water-water energy
 ?TSM      TMS free energy term
 ?QMEL     Quantum (QM) energy with QM/MM electrostatics
 ?QMVD     Quantum (QM/MM) van der Waal term
 ?ASP      Atomic solvation parameter (surface) energy
 ?EHAR     Restraint term for Implicit Euler integration
 ?GEO      Mean-Field-Potential energy term
 ?MDIP     Dipole Mean-Field-Potential energy term
 ?PRMS     Replica/Path RMS deviation energy
 ?PANG     Replica/Path RMS angle deviation energy
 ?SSBP     ???????  (undocumented)
 ?BK4D     4-D energy
 ?SHEL     ???????  (undocumented)
 ?RESD     Restrained Distance energy
 ?SHAP     Shape restraint energy
 ?PULL     Pulling force energy
 ?POLA     Polarizable water energy
 ?DMC      Distance map restraint energy
 ?RGY      Radius of Gyration restraint energy
 ?EWEX     Ewald exclusion correction energy
 ?EWQC     Ewald total charge correction energy
 ?EWUT     Ewald utility energy term (for misc. corrections)
=========  =============================================================

Energy Pressure/Virial Terms:

=========  =============================================================
 ?VEXX      External Virial
 ?VEXY
 ?VEXZ
 ?VEYX
 ?VEYY
 ?VEYZ
 ?VEZX
 ?VEZY
 ?VEZZ
 ?VIXX      Internal Virial
 ?VIXY
 ?VIXZ
 ?VIYX
 ?VIYY
 ?VIYZ
 ?VIZX
 ?VIZY
 ?VIZZ
 ?PEXX      External Pressure
 ?PEXY
 ?PEXZ
 ?PEYX
 ?PEYY
 ?PEYZ
 ?PEZX
 ?PEZY
 ?PEZZ
 ?PIXX      Internal Pressure
 ?PIXY
 ?PIXZ
 ?PIYX
 ?PIYY
 ?PIYZ
 ?PIZX
 ?PIZY
 ?PIZZ
=========  =============================================================

Examples:

1. Save the structure with a lower NOE restraint energy.

   ::

      READ COOR CARD      UNIT 1  ! Read the first structure
      READ COOR CARD COMP UNIT 2  ! Read the second structure
      ENERGY                      ! Compute energy of first structure
      SET 1 ?NOE                  ! save the NOE energy value
      ENERGY COMP                 ! Compute the energy of the second structure
      IF ?NOE LT @1  COOR COPY    ! replace first structure if second has
                                  ! a lower energy.

2. Write some energy values when saving coordinates

   ::

      ....
      COOR ORIE RMS MASS
      ENERGY
      OPEN WRITE CARD UNIT 22 NAME RESULT.CRD
      WRITE COOR CARD UNIT 22
      * Final coordinates
      * energy=?ENER and electrostatic energy=?ELEC
      * mass weighted rms deviation from xray structure is ?RMS
      *


.. _energy_running_average:

Running Energy Averages (ESTATS)
--------------------------------

The ESTATS command is a basic statistical facility that allows the
calculation and manipulation of the mean and variance of the potential
energy and its components over a number of potential energy calculations,
without the need for writing out trajectories or coordinate files--i.e.
the calculations are done "on the fly." ESTATS can be used in dynamics runs
or in other sampling procedures that result in serial calls to the ENERGY
subroutine.  The facility will collect data points at specified sampling
intervals along a collection run for a specified step length and calculate
the running statistics.  An initial portion of the collection run may be skipped
(e.g. for eliminating the equilibration period from the statistics during
dynamics). Results may be written either to standard output or to a file.
The facility will, if requested, "variable-ize" the calculated averages, i.e.
allow assignment of the values to CHARMM script variables.  The facility can
also write the individual potential energy values to a file.

Syntax:

::

   ESTAts   [LENGTH <integer>] [SKIP <integer>] [IPRF <integer>]
            [NPRI <integer>] [IUNW <integer>] [NEPR <integer>]
            [IUPE <integer>]
            [UPLM <real>] [LOLM <real>] [FRPI] [VARI]

            [BOND] [ANGLe] [UREY-Bradley] [DIHEdral] [IMPRoper]
            [VDWaals] [ELECtrostatics] [HBONding] [USER]
            [SBOUnd]  [ASP]

========= ==================================================================
LENGth    length of trajectory (number of total energy calculations)
          from which sampling is to take place  (default 0).
SKIP      specifies a length of energy data points (calls to ENERGY)
          after which the data collection is to begin (default 0).
IPRFreq   specifies the frequency with which data points will be
          collected (i.e. every IPRFrequency energy calculations).
          [BOND], [ANGLE], etc.
          the energy term keywords specify which components of the potential
          energy are to have their statistics calculated.  HBONding is
          the hydrogen bonding energy; USER is the user-defined energy;
          SBOUnd is the solvent boundary potential; ASP is the implicit
          solvation energy (e.g. from eef1).  Statistics on the total
          potential energy are always calculated.
IUNWrite  fortran unit onto which statistics are to be written
          (default is no printing)
NPRInt    period for writing the energy statistics to standard output
          (default is no printing)
IUPE      fortran unit onto which the potential energies are to be
          written (default is no printing)
NEPRint   period for writing potential energies
UPLM      limit above which an energy value will be discarded from the
          statistics (default 99999999).
LOLM      limit below which an energy value will be discarded from the
          statistics (default -99999999).
FPRI      keyword specifying the final statistics are to be written to
          standard output at the end of the collection
STOP      stops the data collection and prints current statistics
VARI      keyword specifying that values of averages and variances will
          be assignable to CHARMM script variables.
          The values can be accessed as follows:
          * ?AENE,?VENE  mean and variance for potential energy
          * ?ABON,?VBON  mean and variance for bonds
          * ?AANG,?VANG  for angles
          * ?AURE,?VURE  for Urey-Bradley terms
          * ?ADIH,?VDIH  for dihedral terms
          * ?AIMP,?VIMP  for improper dihedral terms
          * ?AVDW,?VVDW  for van der Waals
          * ?AELE,?VELE  for electrostatics
          * ?AHBO,?VHBO  for hydrogen bond terms
          * ?AUSE,?VUSE  for user energy
          * ?ASBO,?VSBO  for solvent boundary potential (sbound)
          * ?AASP,?VASP  for solvation term
========= ==================================================================

Note that ALL component energy terms for which statistics are being calculated
must be in the proper range (>LOLM and <UPLM) in order for a given data point
to be included.  Discarded data points will result in statistics that are
based on less than LENGTH data points.  The number of discarded data points
will be printed to standard output with the final statistics at the end of
the collection.

EXAMPLE:

::

   ESTAts LENGTH 1000000 SKIP 100000 IPRFreq 5 NPRINT -1 FPRInt -
   VDW ELEC BOND ANGL IMPR SBOU DIHE -
   IUNWrite 11 NUPRint 10000 NEPR 1000 IUPE 10

This specifies that the statistics are to be done on 180,000 data points
(1,000,000 - 100,000)/5.  Statistics will be done on the specified
energy terms in addition to the potential energy.  Statistics will be
written every 10,000 steps to unit 11 and the potential energies
will be written every 1000 steps to unit 10.  No printing will be
done to standard output (NPRINT -1) except for the final statistics
(FPRInt).

This statistics file is written out according to the following format:

::

   EPOT             10    -1912.29336237      226.63620520
   BOND             10      212.58550818       91.35922427
   ANGL             10      299.99516787       95.65303762
   UREY             10       39.09373234       18.64669506

The first column indicates the energy term. The second column indicates
the number of data points included in the calculations (number of values
over which statistics are taken). The third column gives the average and
the fourth column gives the fluctuation (standard deviation).

Differences between ESTATS and statistics calculated in standard
dynamics:

1) In ESTATS, the denominator in the standard deviation formula is
   (N-1)^(1/2), where N is the number of data points.  In dynamics,
   the denominator is N^(1/2).
2) In dynamics, the initial energy terms are considered step "0" and
   not included in the statistics; hence for a direct comparison, it
   is necessary to specify SKIP 1 in the ESTATS command and increase
   the LENGTH of the collection by 1.


.. _energy_multe:

Multiple Energy Evaluations for Related Conformations
-----------------------------------------------------

The purpose of the MLTE command is to make repeated energy evaluations
for related conformations through rapid reading of pre-computed coordinates.
The MLTE command is invoked from the command line after a valid UPDAte
command has specified the energy function desired.

::

  MLTE [ ATNI atom_i ] [ ATNJ atom_j ] [ FILI filename_i ] [ FILJ filename_j ]
         [ OUTF filename_output ]
         repeat ( [ETER eterm] )

Cartesian coordinates for a collection of coordinates for two sets of atoms
are stored in "write coor dumb" format in files whose names are indicated by
FILI and FILJ, respectively.  The number of the first atom of group_i is given
by ATNI and that of group_j is given by ATNJ.  The energy for all
conformational pairings (each conformation of group_i with each conformation
of group_j) will be computed and printed in the file indicated by OUTF.  The
energy terms to be included is indicated with repeated invocations of the
ETER command, with energy term names listed in the substitution section.
