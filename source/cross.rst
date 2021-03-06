.. py:module:: cross

=================================================
Reactive Molecular Dynamics with Surface Crossing
=================================================

by  David R. Nutt
and Jonas Danielsson (jonas.b.danielsson@gmail.com)
and Markus Meuwly (m.meuwly@unibas.ch)

Questions and comments regarding RMD should be directed to

* Stephan Lutz (stephan.lutz@unibas.ch)
* Markus Meuwly (m.meuwly@unibas.ch)

Reference:

* D. R. Nutt, M. Meuwly, Biophysical Journal 90 (4), 1191-1201 (2006).

The Reactive Molecular Dynamics (RMD) method allows to simulate
dynamics using multiple potential energy surfaces provided by the user, such
that the dynamics always takes place on the lowest surface. Crossings are
detected automatically and occur by a smooth switching centered in time around
which the crossing was detected. The current implementation assumes
that the long-range interactions are handled with the SHIFT command
for electrostatics and the SWITCH for the Lennard-Jones potential. To include
RMD in the compilation, the flag RMD should be included in pref.dat

.. _cross_syntax:

Syntax
------

Syntax for the Reactive Molecular Dynamics commands

::

  RXMD   UPAR integer [ XTIM integer ] [ UNIT integer ] [ FREQ integer ] -
                      [ BCHK real ] [ BVDW ]

======= ====== ==================================================================
UPAR    Int    Unit containing the parameters for the two surfaces.
               This unit must be open (FORMatted) for reading when RXMD
               is called.

XTIM     10    Sets the number of timesteps during which the smooth
               switching between the two surfaces takes place around a
               crossing.

UNIT     -1    Unit to which a coordinate dump (in PDB format) will occur
               at the crossing point.

FREQ      1    Frequency (in time-steps) at which RXMD print out the energy
               difference between the currently active surface and the surface
               which comes closest to a crossing.

BCHK     -1    Cutoff distance for energy comparison checks for surfaces
               including a dissociated bond on the current surface which is
               present on another one. If the distance between the bonding
               atoms is larger than this value the energy between those
               surfaces is not compared at all. Inactive for negative values.

BVDW           Activates special Lennard-Jones potential for the Van-der
               Waals interactions between dissociated bonds. This interaction
               replaces the standard Van-der Waals interaction coming from
               the ATOM section in the surface crossing parameter file or
               the standard CHARMM Van-der Waals interaction as defined by
               the CHARMM parameter file. It is only considered for the 1-2
               (possible bonding atom partners) and 1-3 interactions (related
               angles which must be listed in the ANGL section). The
               description of individual Van-der Waals parameters (epsilon &
               sigma) for every broken bond is described below.
======= ====== ==================================================================


.. _cross_description:

Description of the Reactive Molecular Dynamics command
------------------------------------------------------

To invoke RMD, the CROS command should be given before the
DYNAmics command. At the point where RXMD is called, it is important
that the crossing parameter file is opened on the unit UPAR as FORMATTED.
It is also recommended to proceed the RXMD command with an UPDAte so
that bonded parameter lists are up-to-date.

Energy terms present in both the PSF and the crossing parameter file, will
have the term removed and replaced with that defined in the parameter file.
In this sense, the CROSS module can be used to introduce specific parameters
in the reactive center, or replace a harmonic bond with an anharmonic one, etc.

It is possible to do minimizations on the potential energy surfaces defined in
the crossing parameter file via the RXMD command. During the minimization, the
switching function is turned off so that the minimal surface is followed in
each step.

A typical input sequence for a RXMD simulation is as follows

::

  OPEN UNIT 9 READ FORMATTED NAME cross.dat

  UPDATE
  RXMD UPAR 9 XTIM 10 UNIT 14

  OPEN UNIT 14 WRITe FORMatted NAME crossing.pdb

  OPEN UNIT 11 READ FORMatted NAME old.res
  OPEN UNIT 12 WRITe FORMatted NAME new.res
  OPEN UNIT 13 WRITe UNFORMatted NAME  cross_traj.dcd
  DYNAmics LEAPfrog VERLET RESTart -
  ...

In the first step, RXMD will detect which potential energy surface is
lower and initiate the dynamics on that. The output from a
RXMD simulation is the same as that from a standard dynamics run, with
the addition of the timepoints of the detected crossings
(and the structure at the crossing dumped in UNIT in PDB format). If IPRINT
is > 5, additional detailed information about added and removed energy
terms and the energy difference between the potential energy surfaces are
printed. The crossing procedure uses an energy term called RXMD that can be
bypassed by the SKIPE command. The removal of other energy terms with SKIPE
will also remove the corresponding energy terms in the user-defined PES in a
consistent way.

For the Electrostatics, any method can be used, but the energy difference
between the two surfaces is always calculated with the SHIFTed cut-off
method for the electrostatics, and will therefore not be consistent if
any other method used for the total energy and forces.

The use of the RXMD module together with any procedure that modifies the
PSF is for the moment not supported, since the method assumes the PSF to be
kept as it is defined at the moment the RXMD command is given.

A detailed list of nonbonding interactions which are removed (exclusions) or
introduced (inclusions) relatively to the preloaded PSF bond definitions for
every defined ARMD surface is printed if PRNLEV is set to 6 before of the RXMD
command is called. It is recommended to check this list and compare it to the
PSF definition when setting up a new ARMD system to avoid possible
misinterpretations of the ARMD algorithm.


.. _cross_extra_parameter_file:

Extra Parameter Files
---------------------

The user has to provide an extra parameter file to describe the two
potential energy surfaces involved in the crossing. The format should look
like:

::

  SURF nr_of_surfaces (int)
  surface1 surface2 delta1
  surface1 surface3 delta2
  surface1 surface4 ...
  ...
  [ BART
  surface1 surface2 btol1
  surface1 surface3 btol2
  surface1 surface4 ...
  ... ]
  ATOM nr_of_atoms (int)
  atom1 q1_1 epsilon1_1 sigma1_1 q1_2 epsilon1_2 sigma1_2 (int,real*6) ...
  atom2 q2_1 ...
  ...
  BOND nr_of_harmonic_bonds (int)
  atom1a atom1b k1_1 r1_1 k1_2 r1_2 (int*2,real*4) ...
  atom2a atom2b k2_1 ...
  ...
  MORS nr_of_morse_bonds (int)
  atom1a atom1b d1_1 r1_1 b1_1 d1_2 r1_2 b1_2 (int*2,real*6) ...
  atom2a atom2b d2_1 r2_1 ...
  ...
  ANGL nr_of_angles (int)
  atom1a atom1b atom1c k1_1 t1_1 k1_2 t1_2 (int*3,real*4) ...
  atom2a ...
  ...
  DIHE nr_of_dihedrals (int)
  atom1a atom1b atom1c atom1d k1_1 m1_1 p1_1 k1_2 m1_2 p1_2 (int*4,real,int,
  real*2,int,real) ...
  ...
  [ BVDW nr_of_lj-parameters (int)
  atom1a atom1b totsig1 toteps1 rep-exp1 att-exp1
  atom2a atom2b totsig2 toteps2 rep-exp2 att-exp2
  ... ]

Symbols:

============= ==================================================================
delta         potential energy shift between any surface > 1 and surface 1
btol          energy tolerance for surface switching (Default: 0.0001 kcal/mol)
q             partial charge
sigma,epsilon vdw parameters
k             force constant
r             equilibrium bond length
d             dissociation energy
b             beta parameter in Morse potential
t             equilibrium angle
m             dihedral multiplicity
p             phase angle
totsig,toteps combined vdw parameters for atom_a and atom_b
============= ==================================================================

``_y`` means that this parameter should be used on surface y, so q1_2 means the
partial charge of atom 1 on surface 2, and m3_1 means the dihedral
multiplicity of dihedral 3 on surface 1.

All forcefield terms that differ between any of the states or from the PSF
should be defined. If a term is absent in one of the states
the corresponding force constant (BOND, ANGL, DIHE) or dissociation energy
(MORS) should be given a negative value. For dihedrals, a multiplicity of 0
indicates an improper dihedral.

Note that all blocks must be present (expect of BART which is optional),
even if no new energy term of that kind is defined. For example, even if no
new dihedrals are defined, the file should still has a line reading
'DIHE 0'. There should be no comments or empty lines in this file.

The BART section defines an additional list of energy thresholds that allows
events where the potential energies come close (but not quite cross) to
induce surface switching. If this section is missing every btol defaults to
10e-4 kcal/mol.

The usage the optional BVDW section is described in more detail below.

An example of the parameter input file is given below, for a NO molecule
with and without a bond to a heme group. (94 and 95 is the NO ligand, 21 the
heme iron, 22-25 the pyrrole nitrogen, and 15 the nitrogen of the histidine
coordinating on the opposite side)

::

    SURF 2
    1 2 -25.0
    BART
    1 2 0.001
    ATOM 2
    94 -0.063 -0.200 2.000  0.021 -0.200 1.850
    95  0.063 -0.159 2.050 -0.021 -0.120 1.700
    BOND 6
    94   95   1147.5 1.151 826.5 1.141
    21   22    270.2 1.958 270.2 2.100
    21   23    270.2 1.958 270.2 2.100
    21   24    270.2 1.958 270.2 2.100
    21   25    270.2 1.958 270.2 2.100
    15   21     65.0 2.200  65.0 2.100
    MORS 1
    21 94 -1.000 0.000 0.000 30.0 1.740 3.200
    ANGL 14
    21 94 95 -1.000 0.000 35.0 134.0
    22 21 94 -1.000 0.000 50.0  90.0
    23 21 94 -1.000 0.000 50.0  90.0
    24 21 94 -1.000 0.000 50.0  90.0
    25 21 94 -1.000 0.000 50.0  90.0
    15 21 94 -1.000 0.000 50.0 180.0
    22 21 23  14.39  90.0 80.0  90.0
    23 21 24  14.39  90.0 80.0  90.0
    24 21 25  14.39  90.0 80.0  90.0
    25 21 22  14.39  90.0 80.0  90.0
    15 21 22  50.0   90.0 50.0 107.0
    15 21 23  50.0   90.0 50.0 107.0
    15 21 24  50.0   90.0 50.0 107.0
    15 21 25  50.0   90.0 50.0 107.0
    DIHE 0


.. _cross_lj_treatment:

Special L-J treatment
---------------------

To prevent clashes between different tertiary structure elements, the
standard CHARMM forcefield generally sets the VdW parameters epsilon and
sigma to unnaturally large values. To obtain reasonable transition
barriers for a bond formation reaction, these parameters need to be
scaled down for the specific atom pairs and preferentially also between
atom pairs in a 1-3 distance defining an angle (ANGL) over the bond which
is described by this atom pair.

The BVDW option gives the user the option to define individual parameters
for epsilon and sigma of the VdW interaction between an atom pair
describing a bond or its 1-3 interactions which need to be
listed in the BOND, MORS or ANGL section of the parameter file. Furthermore
the corresponding bond or angle interaction must be deactivated on at least
one of the defined surfaces. Additionally, the exponents of the
Lennard-Jones 12-6 potential as defined in the CHARMM force field must be
redefined for each of this specific VdW interactions which allows for the
application other functional forms. If BVDW was requested, epsilon, sigma
and the repulsive and attractive Lennard-Jones exponents for any broken
bond or angle defined by two atoms are read from the BVDW section in the
RMD parameter file.

An adapted example for the parameter input file making use of the BVDW
option is given below. The Morse potential describing a bond between atoms
21 and 94 on surface 2 is deactivated on surface 1. For the VdW
interactions acting in this state between the two atoms and the atom pairs
located in a 1-3 distance around them (identified by force constants
k of -1.0 in the ANGL section) specific parameters are added after the
DIHE section.

::

    SURF 2
    1 2 -25.0
    ATOM 2
    94 -0.063 -0.200 2.000  0.021 -0.200 1.850
    95  0.063 -0.159 2.050 -0.021 -0.120 1.700
    BOND 6
    94   95   1147.5 1.151 826.5 1.141
    21   22    270.2 1.958 270.2 2.100
    21   23    270.2 1.958 270.2 2.100
    21   24    270.2 1.958 270.2 2.100
    21   25    270.2 1.958 270.2 2.100
    15   21     65.0 2.200  65.0 2.100
    MORS 1
    21 94 -0.250 1.500 0.000 30.0 1.740 3.200
    ANGL 14
    21 94 95 -1.000 0.000 35.0 134.0
    22 21 94 -1.000 0.000 50.0  90.0
    23 21 94 -1.000 0.000 50.0  90.0
    24 21 94 -1.000 0.000 50.0  90.0
    25 21 94 -1.000 0.000 50.0  90.0
    15 21 94 -1.000 0.000 50.0 180.0
    22 21 23  14.39  90.0 80.0  90.0
    23 21 24  14.39  90.0 80.0  90.0
    24 21 25  14.39  90.0 80.0  90.0
    25 21 22  14.39  90.0 80.0  90.0
    15 21 22  50.0   90.0 50.0 107.0
    15 21 23  50.0   90.0 50.0 107.0
    15 21 24  50.0   90.0 50.0 107.0
    15 21 25  50.0   90.0 50.0 107.0
    DIHE 0
    BVDW 7
    21 94 0.316 2.75 12 6
    21 95 0.245 3.02 12 6
    22 94 0.316 2.75 12 6
    23 94 0.316 2.75 12 6
    24 94 0.316 2.75 12 6
    25 94 0.316 2.75 12 6
    15 94 0.316 2.75 12 6

