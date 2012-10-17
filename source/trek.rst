.. py:module:: trek

========================================================
TReK: a program for Trajectory REfinement and Kinematics
========================================================

.. note::

   Current Version 2.10 , July 5-2003

   Please report problems or send suggestions to :

     | Stefan Fischer
     | Tel. +(49)6221-548879
     | e-mail: stefan.fischer@iwr.uni-heidelberg.de

   Check for application examples and bugfixes under:
   http://www.iwr.uni-heidelberg.de/groups/biocomp/fischer


.. note::

   - For the busy user: please read at least the sections marked "!!!"
     before starting to use TReK and CPR.

   - This interface to TReK replaces the TRAVel module of CHARMM, which is
     no longer supported.  An effort was made to reproduce the "look&feel"
     of TRAVEL, whose input-scripts should work with TReK.


TReK is a collection of tools to find and smooth the minimum-energy path
between two known structures representing the reactant and the product states
of a reaction (here, the word "reaction" is used to design both a chemical
reaction and/or a conformational transition). The energy can be an empirical
energy (MM), a pure quantum (QM) potential,
or a combined QM/MM energy.

The path is described by a series of structures which form a chain, or
trajectory in conformational space, whose two end-points are the reactant
and product structures.  Watched in sequence, these molecular structures
constitute a "movie" of the essential motions along the reaction, hence
the appellation "kinematics".

The main tool is a new implementation of the CPR (Conjugate Peak Refinement)
algorithm (S.Fischer and M.Karplus, Chem. Phys. Letters 194, p.252, 1992),
which refines an initial guess of the path (for ex. the linear interpolation
between end-points) and finds the trajectory that follows the valleys
of the high-dimensional energy surface, as well as all the exact saddle-points
along that path.

Other tools are the SCM (Synchronous Chain Minimization) algorithm
or the SDP (Steepest Descent Path) method, which both allow to smooth
the shape of a path whose saddle-points are already known.

Other tools allow to analyze or manipulate reaction-paths.


.. _trek_syntax:

Syntax of the TReK program
--------------------------

Keywords in [...] are optional. Default values are given in (...).
Choose one from list :  {...|...|...} or [...|...|...]

Main command (entering and leaving TReK)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

   TREK  [MAXPoints int (100) ]   [ { XAXI | YAXI | ZAXI } [ ROTAtion ] ]
                                  [ SCALe [{ COMP | WMAIn }] ]

   { END | QUIT }

Subcommands (within  TReK program)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

   [VERBose int (2) ]
   [DISPlay int (80)]
   [CHROno { RESEt | PRINt } ]

   TRAJectory READ [UNIT int (40)] [REFErence] [ori]
     file1.crd
     file2.crd
        .
        .
     fileN.crd
   DONE

   TRAJectory READ  NAME file.dcd [UNIT int (40)] [REFErence] [ori]
                                  [BEGIn int (1)] [SKIP int (1)] [STOP int (0)]

   TRAJect WRITe NAME file.dcd  RESTART file.cpr  [UNIT int (40)]

   TRAJectory ANALyze [SCAN [STEP real]]

   TRAJectory { INCRease | DECRease [PEAK] [MINI] [SADD] }
              [ NGRId int | STEPsize real ] [ori]

   SADDle int [ int [ int ...]]

   CPR  [SADDle] [cpr_exit] [saddledef] [restart]
        [scanning] [relaxspeed] [oscillation] [miscrefine] [linextrem] [ori]

   SCM  [NCYCle int (1)]    [PROJtol real (.01)] [ANGLe real (90.)]
        [MINUpdate int (5)] [MAXUpdate int (30)] [linextrem] [ori]

   CROS [MINDist real (e-4)] [MINStep int (50)] [FIRStep real] [ANGL real (20)]

   SDP  [SAVDistance real (.05)] [MINGrad real (e-3)]
        [NREActant int (maxpoints/2)] [NPROduct int (maxpoints/2)]
        [MODE int (4)] [ANGLe real (30.)] [MINCycle int (50)] [linextrem]

   COPY [COMP] { SADDle | ORDEr int | INDEx int }


Detailed keywords
^^^^^^^^^^^^^^^^^

::

   cpr_exit    :== exit-criteria for CPR :

                   [NCYCle int (1)]  [NECAlls int (999999)]  [HIGHsad]

   saddledef   :== saddle-point criteria :

                   [SADGrad real (.05)] [SADMini int (sqrt(N))]

   restart     :== read a TReK restart-file :

                   RESTart filename [SCAN] [UNIT int (40)]

   scanning    :== interpolation step-size for path-scanning :

                   [INTErpol int (3)] [NGRId int (5) | STEPsiz real]
                   [STPLowest real (0)]  [SCAN]

   relaxspeed  :== stringency of path-relaxation :

                   [TOL1proj real (1.)] [TOL2proj real (3.)] [LINMax int (2000)]

   oscillation :== oscillation detection and prevention :

                   [LOOPreduc int (sqrt(N))] [FRAMe int (10)]
                   [TOLOscill real (.15)] [PROJincr real (2.)] [MAXOscil int (4)]

   miscrefine  :== miscellaneous :

                   [DELTa real (-e-7)]  [DISPlay int (0)]
                   [NTANgent int (3)]  [REMOvemod int (0)]

   linextrem   :== one-dimensional line-extremization :

                   [BRAKetstep real (.1)] [BRKScale real (2.)] [LXEVal int (8)]
                   [BRKMagnif real (5.)] [FIRStep real (.05)] [EXITmode int (3)]
                   [TOLMax real (e-4)]  [TOLGrad real (.05)] [ATOMax real (0.75)]
                   [TOLStep real (e-10)] [TOLEne real (e-7)]

   ori         :== coordinate orientation :

                   [ ORIEnt | NOORient ] , default is NOORient.


CHARMM command variables set by the CPR command
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

   ?SADI       :== index of the highest fully refined saddle-point
   ?SADO       :== order along the path of that saddle-point
   ?SADE       :== energy of that saddle-point


Note about file-names
^^^^^^^^^^^^^^^^^^^^^

::

   Examples of allowed names are

      ~/subdir/file.crd
     ../subdir/file.crd
      ./subdir/file.crd
        subdir/file.crd

In case the file is already open, it must be closed before entering TReK,
otherwise the "TRAJ READ" command will fail.

.. _trek_maincmd:

TReK Main Command
-----------------

Before invoking TReK for the first time, the following must have been done
in CHARMM :

- The molecular system was GENErated.
- One structure was read into the main CHARMM coordinate-set
  (for ex. the thoroughly minimized structure of the reactant conformer).
- Optional: the desired Images were set-up, Image Centering turned off.
- Optional: the desired atoms were fixed with ``CONS FIX SELECT ... END``.
- Optional: the desired QM regions and link-atoms were setup.
- The FASTER mode was set to a value compatible with the non-bond settings.
- **AND at least one ENERgy call was made, with the above settings being the same as used for the minimization of the path end-points.**

Once TReK has been launched, a typical session would involve :

- Reading the structures of the initial path, and if applicable, declare
  already refined saddle-points (see :ref:`typic1 <trek_input>`).
- Refining the path with CPR or SCM commands (see respective sections).
- Saving the resulting path (see :ref:`typic2 <trek_output>`).
- Optional: Printing the energy profile along the current path
  (see section "TRAJECTORY ANALYZE").
- Optional: Copying a specified path-point to the CHARMM coordinate-set,
  for further analysis (see :ref:`typic3 <trek_output>`).
- Leaving TReK with QUIT or END.


END !!!
^^^^^^^

Exits the TReK program, back into the CHARMM command processor, without
terminating the TReK session.
All TReK data-structures are maintained in memory, allowing to re-enter TReK
to continue refinement. In that case, none of the keywords of the TREK main
command may be used (such as MAXPoints).


QUIT !!!
^^^^^^^^

Terminates the TReK session, erasing all its data structures and freeing all
its memory.  It allows to re-start TReK as if it had never been called before.


Optional keywords on the TReK command-line are :

MAXPoints int !!!
^^^^^^^^^^^^^^^^^

When entering the TReK program from CHARMM, it is possible to specify the
maximum number of structures that will make up a path during a
given CPR or SDP refinement.  If this value is reached during refinement and
even more path-points are needed, then the refinement will stop gracefully,
allowing to save the path with a TRAJECTORY WRITE command.
Later, restart TReK with a larger value for MAXPOINTS to continue the
refinement.

From experience of transitions in proteins, about 20-30 points
are needed for each saddle-point. For ex. if less than 5 saddle-points are
expected along the path, set MAXPoints to ~150.  Large conformational
transitions in proteins (ex. hemoglobin T to R) can have more than 100
saddle-points.


Using TReK with IMAGES
^^^^^^^^^^^^^^^^^^^^^^

TReK supports IMAGES in CHARMM.  Before starting TReK,

- the image centering MUST be turned off :         IMAGE FIXED SELE all END
- the image updating should be set to automatic :  NBOND IMGFRQ -1

::

   [ { XAXI | YAXI | ZAXI } [ ROTAtion ] ] keywords :

When the only involved symmetry operation is a
translation along a single axis and/or a rotation along a single axis, then
this axis and/or that rotation must be declared explicitly upon starting TReK.
Currently, only the three main axis are supported.

If there are translations in more than one dimension (e.g. a membrane),
then these keywords are not necessary and should NOT be used.

Examples:

1) A periodic DNA helix along the y-axis requires : TREK YAXI ROT
2) An dimeric protein with C2 symmetry (one monomer per image) : TREK ZAXI ROT
3) A periodic chain along the x-axis, without rotational symmetry : TREK XAXI
4) A lipid-bilayer with periodic boundaries : TREK
5) A crystal : TREK


SCALE keyword
^^^^^^^^^^^^^

It is possible to find the minimum-energy path on a "stretched" energy surface,
where each atomic coordinate is multiplied by a scaling factor greater or
equal to 1.  This allows for example to express the path in terms
of the mass-weighted coordinates.

Coordinate-scaling is activated by issuing the SCALe keyword on the TReK
command.  By default, the scaling-factors are taken for each atom-coordinate
from the X,Y,Z values stored in the main coordinate-set.  Adding the COMP
keyword will result in these values being taken from the comparison
coordinate set.
Alternatively, if the WMAIn keyword is added, the scaling factors are taken
from the WMAIN vector and the X,Y and Z dimensions of each atom are given
the same scaling-factor.  For ex. using    WMAIN = sqrt(MASS)
will yield the mass-weighted coordinates.

Use the SCALAR command, see :doc:`scalar`, to set the
appropriate values of the X,Y,Z coordinate-vectors or of the WMAIN vector,
before starting TReK.
The scaling-factors should not be modified after starting TReK.

All RMS(distances) printed from within TReK will use the scaled coordinates,
for ex. the reaction-coordinate LAMBDA printed by the "TRAJ ANALyze" command.
However, with COPY or "TRAJ WRITE" the coordinates are re-scaled back to
normal.

When coordinate-scaling is activated, re-orientation is automatically disabled
during the CPR refinement (as if issuing the NOORient keyword), but if desired,
orienting can still be used during the reading of the initial structures by
"TRAJ READ".


.. _trek_usage:

Usage
-----

Many useful suggestions for setting-up a path refinement for a large transition
in a protein are given in the article entitled
"Automated computation of low-energy pathways for complex rearrangements in
proteins: application to the conformational switch of RasP21", by
F. Noe, F. Ille, J.C. Smith & S. Fischer, Proteins 59, p.534-544 (2005).


Using TReK with QM/MM
^^^^^^^^^^^^^^^^^^^^^

TReK can be used to find paths and saddle-points on combined
quantum/empirical energy potentials.  For this to work well, some default
settings must be modified with the CPR command (see section on
"Numerical energy potentials").


Preparing the path end-points !!!
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Make sure that the internal coordinates of atoms that are not significantly
involved in the reaction are the same for the reactant and the product.

For example, a phenyl side-chain has two rotationally equivalent orientations,
in which the delta1 and epsilon1 atoms are exchanged with the delta2 and
epsilon2 atoms.

When the reactant and the product coordinate sets are from
different origins (for example from different X-ray structures) or the result
from different calculations (for example separated by a long dynamics run
or an HBUILD), then for every group or side-chain that has equivalent
(for ex. Phe, Tyr, Asp, Glu, -NH2, -NH3, -CH3 and -CH2- )
or nearly equivalent orientations (for ex. Val, Leu, Lys, Arg),
it is worth checking that the corresponding atoms are named consistently
in the reactant and product coordinates. Otherwise, CPR will proceed with
refining all the transition paths associated with these changes in equivalent
positions (rotating the Phe ring in the example above, or exchanging the H's
of a -CH2- group), which is very CPU time consuming and will make it
difficult to analyze the reaction one is interested in.

.. warning::

   In general, it is more efficient to use end-point structures that are
   very well energy-minimized, for ex. with RMS(gradient) < 0.0001 .

This is best achieved by minimizing them with the Newton-Raphson
(ABNR) method in CHARMM.

For molecules with many local minima (ex. proteins), it
is advisable to run some quenched molecular dynamics (MD): take structures at
periodic intervals along a MD trajectory and minimize them. Use the ones with
lowest energy as end-states of the reaction. An alternative is to do annealing
MD (i.e. slowly decrease the temperature).

Using end-points which are poorly minimized or that are stuck in high-energy
local minima will result in a path whose energy-profile "hangs" between the
end-states, without clear activation barriers.


Fixing atoms !!!
^^^^^^^^^^^^^^^^

When the reaction involves only a sub-region of the molecular system
(e.g. catalysis in a protein), then fixing the atoms that do not really
participate in the reaction speeds up the path-calculation significantly:

If some atoms have been fixed, then the fixed atoms of the reactant,
of the product and of all initial intermediates MUST have
exactly the same Cartesian coordinates.
If they differ, then TReK will set the fixed atom coordinates of all
structures to those of the reactant and give a warning !

To find out how to fix as many atoms as possible without biasing the
transition, progressively add more distant atoms to the moving region and
recompute the path.  When the energy-barrier does not change significantly
with more atoms, the final set of moving atoms defines the molecular region
that really participates in the reaction.

.. _trek_input:

Input
-----

The TReK program has its own I/O commands.

TRAJECTORY READ !!!
^^^^^^^^^^^^^^^^^^^

The purpose of this command is to read the initial path to be refined.
In its simplest form it consists of the two end-points of the reaction.
If desired, one or several intermediate points can also be read.
The points are read in sequence, as they would appear along the reaction:

1) Reactant structure.
2) Intermediate(s).
3) Product structure.

This sequence of structures defines the initial-guess for the path, which
is build by interpolating linearly in Cartesian coordinates from one
structure to the next.
If the two end-points are identical (for ex. when studying the 360 degree
rotation of a side-chain), then at least 2 additional intermediates must be
provided (obviously, since it is impossible to define a linear interpolation
path between identical structures) !  More generally, any two structures
adjacent along the initial path are not allowed to be identical.

Intermediate(s) should also be provided when direct linear interpolation
between the two end-states would result in initial path structures which
have extreme energies. For ex. the 180 degree flip of a phenyl side-chain.
In that example, an intermediate could be constructed by rigidly rotating
the Phe ring by 90 degrees. Alternatively, build the initial path by combined
interpolation in Cartesian and Internal coordinates as described in the next
section (see :ref:`ringrotate <trek_initialpath>`).

Varying the initial path allows to search for different minimum energy paths:
Using different intermediates will direct the refinement toward different and
sometimes better paths. This is analogous to single structure optimization,
where a structure of global minimum energy is searched by starting many
minimizations from different initial conformers.

The command is also used to read a partially refined path, whose refinement
is to be continued (i.e. the initial path is the partially refined path).

There are two ways to read coordinate-sets :

  1) From a series of formated CHARMM coordinate files.

AND/OR

  2) From an unformatted CHARMM dynamics trajectory file.

The "TRAJECTORY READ" command can be issued several times consecutively (in
either of these two forms). The successive path-points will be appended to
constitute the input-path used by CPR, SCM, CROS or SDP.

When reading from an unformatted dynamics trajectory file (for example the
output of an earlier refinement), it is possible to restrict the input
to subsection of the trajectory file with the BEGIN, SKIP, STOP keywords
(see examples below in :ref:`trek_cprcmd_exampl`).

Once a path refinement has been started (for ex. by issuing the CPR command),
no more path-points can be read and the "TRAJECTORY READ" command is disabled.


Reference structure for orienting !!!
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If there are no fixed atoms (see above), then it is advisable to remove the
rotational and translational degrees of freedom. This is done by aligning
the path-points onto a reference structure (so as to minimize the
RMS-difference in coordinates between each path-point and this reference
structure, as would the CHARMM command  "COOR ORIEnt RMS").
This is achieved by issuing the "ORIEnt" keyword on the "TRAJECTORY READ"
command line.
This is generally recommended for the path-points of the initial guess-path.
It is not recommended for the path-points of a partially refined path,
whose refinement is being continued.  The default setting is "NOORient".

By default, the first structure that is read is used as the reactant structure
(i.e. 1st path-point) AND is also used as the reference structure for
orienting.

Alternatively (recommended procedure), the reference structure can be
read explicitly, by specifying the "REFERence" keyword along with
the "TRAJ READ" command.
In that case, the first structure to be read becomes the reference point,
but does NOT belong to the path. This is useful when the refinement of a very
long path is split into several refinements of path subsections, each of
which should align path-points onto a unique reference structure.

When saving the path to a file (TRAJ WRITE NAME ...), the reference point
is NOT written as part of that path.

The reference point can only be specified once, before other path-points
are read. If desired, the reference point can be centered and oriented
according to its principal axis (in the same way as would "COOR ORIEnt"),
by specifying the "ORIEnt" keyword with the "TRAJ READ REFER" command.

Path-points that are subsequently added or modified by a path-algorithm
(CPR, SCM) remain essentially aligned with the other path-points.
They can be re-aligned onto the reference point by specifying the "ORIEnt"
keyword on the command-line of the respective algorithm, although this is
not recommended in most cases.

If some atoms have been fixed, or if there are crystal images, then all
coordinate re-orienting is automatically disabled (the reference point can
still be read, but it will be ignored).


Declaring saddle-points !!!
^^^^^^^^^^^^^^^^^^^^^^^^^^^

After the initial-path has been read, if some path-points had already been
refined to saddle-points in previous CPR runs, it is necessary to declare
these path-points as saddle-points before using commands like "CPR", "SCM"
or "TRAJ DECR".

The purpose of this is :
- Before CPR, to avoid that CPR must re-refine those points all over again (very CPU consuming).
- Before "SCM", to leave the saddle-points unchanged.
- Before "TRAJ DECR", to prevent that saddle-points get removed.

There are two ways to declare path-points as saddle-points :

1) The "SADDLE" command.

   For example if the saddle-points are in sequential position n1, n2, ...
   along the path (do not confuse the nX with the point-indices IDX of the
   previous CPR-run!), then the command to issue is

      ``SADDLE n1 n2 ...``

   The command can be re-issued on successive lines, if necessary.
   Ignore the warning printed after each "SADDLE" command
   (see :ref:`Bugs <trek_bugs>`).

2) An alternative to using the SADDLE command, is to read the restart-file
   corresponding to the current path (see :ref:`CPRcmd <trek_cprcmd>`).
   This is done on the CPR command-line.

   - If preparing to run a CPR calculation, then the restart-file should
     only be read at the moment of issuing the first CPR command.
   - If preparing to use another tool (e.g. SCM, etc.), then the CPR command
     can be called with "NCYCLE 0", so that no CPR-cycles are performed,
     but the command is only used to read the restart-file :

          ``CPR NCYCLE 0  RESTART filename``

     This flags the saddle-points according to the LinMin array of the
     restart-file : all points which are listed with     LinMin < 0
     are considered as saddle-points (no checking is done to verify that
     these points really qualify as saddle-points).

     .. note::

        If a restart-file is used, do not declare saddle-points twice
        by listing them again with the "SADDLE" command.


.. _trek_initialpath:

Create initial path by interpolating sidechains in IC
-----------------------------------------------------

This section describes a method for generating an initial path for CPR
when the two end-structures of the protein are so different that the
default interpolation in Cartesian coordinates gives poor results.
The method combines Cartesian and IC interpolation with sidechain shrinking.
It has been described in detail in:   Proteins 59, p.534-544 (2005).

The tool to generate this initial path consists of a PERL-script (called
vector.pl) and a CHARMM-script that uses standard commands.
These scripts are located in the CHARMM-tree subdirectory ``./support/trek/initial_path/``

Summary
^^^^^^^

Interpolation in Internal Coordinates (IC) describes well the torsional
transitions between rotameric states of the sidechains. However, when
IC-interpolation is applied to the backbone-torsion angles, it can lead to
severe disruption of the backbone fold. In contrast, Cartesian interpolation
approximately preserves the backbone fold in most cases (if the protein
remains compact during the transition), but often leads to extreme
deformation of the sidechains.
Therefore, the two interpolation methods are combined here:
First, the backbone atoms are interpolated in Cartesian coordinates
so as to preserve the backbone fold, and then the sidechain atoms are
built onto the interpolated backbone, using internal coordinates
that are interpolated between the values of the end-states.

The initial path build by this "combined interpolation" can be further
improved, so as to avoid unrealistic events such as bond-crossing
or ring-penetration. These are effectively avoided by shrinking all sidechains
in each intermediate point of the initial path.
This is achieved by reducing all bond lengths of the sidechains
(for ex. to half the original size) before building them onto the backbone.
At first glance, this approach would seem to be unphysical and
would result in energetically unfavorable paths. However, the minimization
process applied to the path-points during the CPR computation rapidly
restores the shrunken sidechains to their normal size, while undesirable
events no longer occur in the resulting path.


vector.pl
^^^^^^^^^

For both the combined interpolation and the shrinking, the PERL script
"vector.pl" is called from the CHARMM-script "comb_interpol.inp".

vector.pl reads two IC tables from CHARMM and calculates the shortest change
between each pair of torsion-angles in these IC-tables, for ex. -120 --> +120 :
delta(angle) = 120 (and not 240).
It writes out the resulting delta(angles) in an IC table, which can then
be read into CHARMM. The script-parameters are displayed when executing it
without any parameter or with the parameter -h.


comb_interpol.inp
^^^^^^^^^^^^^^^^^

Is the main CHARMM input-file.
Calls

   - user_def.str
   - vector.pl
   - mini_interpol.str  (optional)
   - make_traj.str      (optional)


user_def.str
^^^^^^^^^^^^

Here, the user must define all necessary settings.  For. example some
recommended values are :

::

   SET steps =  20    ! Put 20 intermediate structures in the initial path.
   SET shrink = 0.5   ! Halves the size of the sidechains.


mini_interpol.str
^^^^^^^^^^^^^^^^^

Optional, use only if CPR cannot deal with un-minimized initial intermediates.
If requested, it is called at the end of comb_interpol.inp.
Does a quick and constrained minimization of each frame.


make_traj.inp
^^^^^^^^^^^^^

Creates a CHARMM trajectory from the interpolated frames for viewing in VMD.
An example of the resulting initial path (4 intermediates, shrinking to 50%)
is available under ./support/trek/initial_path/frames/


Remarks
^^^^^^^

If there are no fixed atoms, then the product state is "COOR ORIENTed" upon
the reactant state.

When doing CPR on the so-obtained initial path, to avoid that many of
the initial intermediate path-points are deleted right away by CPR,
use a larger value of INTERpol for a ~500 CPR-cycles, for ex.:

   ``CPR NCYCLE 500 INTER 10 ...``


.. _trek_cprcmd:

Conjugate Peak Refinement (CPR)
-------------------------------

Description of the algorithm
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Provided with an initial guess for the transition path between a reactant
and a product state, CPR refines the minimum-energy path and accurately
determines all saddle-points of the high-dimensional energy surface along
that path, in systems with up to several thousands of degrees of freedom.
CPR uses the configuration energy and its first derivative.

Ref.:  S.Fischer and M.Karplus, Chem. Phys. Letters 194, p.252-261, (1992).

The method is suited for the study of large-scale conformational change in
macromolecules using empirical energy functions, or for simulating
enzymatic reactions with the help of combined quantum/classical (QM/MM)
force-fields.  Because non-relevant saddle-points are abundant on the
energy surface of even small molecules, CPR's reliability in finding the
transition-state that connects the desired reactant and product makes it
useful also for the study of reactions with pure QM potentials.

An essential feature of the way CPR finds a path is that it does not drive
(steer) the reaction along some pre-defined reaction coordinate, which
would have been pre-defined as a function of the degrees of freedom.
Rather, CPR self-consistently optimizes a chain of points describing the
path, without applying any constraints that would affect the shape or the
length of that path.  The primary focus of the algorithm is on accurately
determining the saddle-points along the path.  Path-points which are not
saddle-points are optimized only to the extent that the energy along the
path segments connecting these points decreases monotonously from the
saddle-point to its adjacent local minima.  The purpose of points between
saddle-points is to ensure the continuity of the path, not to find the
bottom of the energy valley, which can be obtained more efficiently
afterwards with the SCM method or by steepest descent on each side of the
saddle-point(s).  The number of path-points is not fixed during the
refinement, but is allowed to increase and/or decrease.  This enables CPR
to accommodate any degree of complexity of the underlying energy surface,
including the presence of multiple saddle-points along the path.

The method is a heuristic procedure, where in each CPR cycle:

- the highest local energy maximum along the path (called the "peak") is
  determined by stepping along the path and computing the energy at each
  step (this is referred to as "scanning" the path).
- the path is modified by either improving, removing or inserting one path-
  point, so that the new path avoids the peak.

Improved as well as newly inserted points are optimized by a controlled
conjugate-gradient minimization, which prevents the point from falling
into local minima along the path and which converges to the saddle-point
if the path is crossing a saddle-region of the energy surface.
The path is fully refined when, after a number of such CPR cycles,
the only remaining energy-peaks along the path are the exact saddle-points.

The settings of the algorithm are independent of molecular size and the
nature of the reaction.


Before invoking CPR for the first time !!!
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following must have been done within TReK :

- Reading the coordinate sets (at least two) of the initial path.
  The reactant and the product conformations must be very well minimized,
  with a RMS(gradient) < 0.001 .
  The initial path-points can be read either from a series of
  CHARMM coordinate files (ASCII format) or from a binary CHARMM dynamics
  trajectoryry file (see :ref:`cpr1invok <trek_input>`).

- Flag those path-points that are already known saddle-points with the SADDLE
  command, or use a CPR restart-file, so that they do not get refined again
  (see :ref: `cpr2invok <trek_input>`).

- Optional:  Setting the desired amount of information printed out during
  refinement with the VERBOSE command.  Also, CHARMM warnings about
  distorted angles, etc., which are abundant in the early phases of a path
  refinement, can be suppressed with "WRNLev 0".


CPR command !!!
^^^^^^^^^^^^^^^

This command starts the path-refinement.  The number of CPR-cycles is
specified with NCYCLE.  It can also be limited by specifying an
approximate number of energy-calls NECAlls, which will not be exceeded.
NECAlls allows better control than NCYCLE to determine
the amount of CPU-time before exit, because the number of energy-calls
made during one CPR-cycle can vary a lot from one cycle to another.

Tentative saddle-points are improved until their RMS-gradient reaches
the value specified with the SADGrad keyword, and this for a successive
number of line-minimizations equal to  SADMini.

In the initial phase of a refinement, only a moderate amount of CPU-time
is spent on refining individual saddle-points. This is achieved by setting
the default value of SADMini to SQRT(n), (where n is the number of moving
atoms) instead to 3n-1.
This allows to refine the overall path and the approximate saddle-points.
If desired, the exact saddle-points can be ultra-optimized in subsequent
CPR runs, by issuing the CPR command with the SADDLE keyword (see below),
the latter being more CPU-time consuming.

If the switch HIGHSAD is used, then the refinement will stop as soon as the
highest tentative saddle-point has been found.
This allows to get a quick estimate of the overall barrier height along the
path, for example to decide whether it is worthwhile to proceed with
refining of the rest of the path.

The default values of all CPR settings (including the values of SADGrad &
SADMini) were optimized for the use with analytical energy and
gradient functions, and are independent of the number of atoms or the
type of reaction.  So in general, they do not need to be modified.
Note that all values are remembered from one CPR call to the next, as well
as when temporarily exiting TReK with the END command.

Upon exit, CPR prints a summary of the energy and location of any
saddle-points found so far. It also summarizes the location of some of
the deepest local minima along the path, if any are present.
If all energy-peaks along the path are saddle-points (according to the SADGrad
and SADMini criteria, see above), then the printed output will contain the
words "FULLY refined".  If the printout says "NOT fully refined", continue
with more CPR refinement-cycles.

The CPR command can be issued several times in the same TReK session,
for example to allow for frequent saving of the path.  Once a path reaches
the status of "FULLY refined", the remaining CPR commands are ignored,
costing no extra CPU-time.

The energy-profile along the current path, and the location of the
saddle-points can be obtained with the "TRAJ ANALysis" command,
or by looking at the content of the TReK restart-file
(see :ref:`eneprof <trek_output>`).


Saving the path !!!
^^^^^^^^^^^^^^^^^^^

After the path has undergone a number of refinement cycles, save
the whole trajectory to a standard CHARMM dynamics file (binary format).
This is done with the TRAJECTORY WRITE command (see :ref:`cprsav <trek_output>`).
This saving should be done frequently, since a given CPR-cycle
can be unpredictably long (unless the NECAlls keyword is used).


Ultra optimizing the saddle-points !!!
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once an overall path has been refined with the default settings for
SADGrd and SADMini, it is possible to continue the refinement of the
saddle-points, to make sure that they are true first-order
saddle-points.

The definition of a first-order saddle-point:

- it has a vanishing gradient, and
- it has exactly one negative eigenvalue of its 2nd derivative matrix (this
  can be checked with the VIBRAN module of CHARMM).

CPR flags a path-point as a first-order saddle-point when:

1) its energy is a local maximum along the path, and
2) it is the result of 3n-1 successive line-minimizations (n = number of atoms)
   which all gave a gradient smaller than SADGrad.

.. note::

   the final RMS-gradients at the saddle-point will be much smaller than
   SADGrad, because of the larger number (SADMini) of line-minimizations
   that will be performed.

The CPU-time required to do the 3n-1 line-minimizations is large.
In practice, the criterion  SADMini = 3n-1  is necessary only for small
molecules (n < 300 atoms).  For larger proteins, it usually suffices that
SADMini ~ 1000, which brings the RMS(gradient) in the range < 10^-6 ,
where the limits of numerical accuracy of the computer are reached anyway.

In addition to SADGrad and SADMini, a number of other CPR settings
have slightly different optimal values for ultra-optimizing the saddle-points.
For convenience, they are all modified simultaneously by adding the
"SADDLE" switch to the CPR command (see example below).
The only effect of this switch is to modify the default settings for the
following keywords :

::

   [SADGrad e-3] [SADMini 1000] [TOL1proj .5] [TOL2proj .5]
   [NTANgent 6]  [ATOMax 0]     [DELTa e-7]   [NOORient]  [LOOPreduc 0]

Specifying another value for any these keywords in addition to
the SADDLE switch has precedence over these values.
For ex. "SADDle SADMini 500" will use SADMini = 500, not 1000.

.. note::

   Do not confuse this CPR switch with the "SADDLE" command used
   to declare saddle-points (see :ref:`cprsad <trek_input>`).


Reading the restart-file !!!
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A CPR restart-file is written each time a path is saved
(see :ref:`cprrestart <trek_output>`).  It can then be read the first time the
CPR command is given during a TReK session.

The restart-file serves two purposes :

1) It allows to flag already refined saddle-points, so that they don't get
   refined again at the beginning of each CPR-run.
   Refined saddle-points satisfy to the criterion  LinMin(i) = SADMini ,
   where LinMin(i) counts how many successive times the RMS-gradient at the
   path-point was less than SADGrd.  To distinguish flagged saddle-points in
   the restart-file (or in the printout of "TRAJ ANAL"), they are given
   a negative sign, i.e.  LinMin(i) = -SADMini (see :ref:`restartinfo <trek_output>`).

   Path-points that have been declared tentatively as saddle-points are also
   listed with a negative LinMin(i) value, but ``|LinMin| < SADMini``.  If such
   points are to be left untouched in the next CPR-run, set their LinMin
   to -999999 in the restart-file (not really worth the trouble, since
   usually it takes only 1 CPR-cycle to flag them as saddle-point again).

2) It allows to save CPU-time upon initialization of CPR, by taking
   the information about the energy along the path from the restart-file,
   rather than having to rescan the whole path. However, when the relative
   energy change along a segment is less than 10^-6, this segment will be
   rescanned anyway to protect against numerical accuracy problems when
   writing/reading trajectory files.

   For each path-point, the energy and RMS(gradient) is listed, one line
   per point, on the first half of the line.
   Listed on the second half of the line is the scanning information for
   the segment preceding the path-point, such as the segment-length LenSeg(i).

It is imperative that the content of the restart-file correspond to the initial
path that was read with "TRAJ READ NAME" in TReK, otherwise unpredictable
behavior will occur.

If some segments in the restart-file do not have the correct
scanning information, for example when two paths and their restart-files
have been spliced, the user can force the re-scanning of the joining segments
by setting their StepSiz(i) field to any value < 0 .

The re-scanning of the whole path can also be forced by using the SCAN keyword,
in which case only the LinMin(i) information of the restart-file will be used
to declare previously refined saddle-points.  For example, this is required
when a path that underwent a number of SCM-cycles is refined further with CPR.

.. note::

   - If any of the numeric fields contain "*****", TReK dies ungracefully.

   - Atom-coordinates in CHARMM-trajectory are writen and read with only
     32-bit precision (i.e. REAL*4), whereas they are stored in memory in
     64-bit precision (i.e. REAL*8). Just like an MD-run, the progress of
     a path refinement by CPR is very sensitive to small changes in the
     initial conditions. Therefore, stopping/restarting a refinement will
     result in somewhat different paths than an uninterrupted refinement.

.. _trek_cprcmd_exampl:

Example CPR input-files !!!
^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

   * Beginning a CPR refinement.
   *
   ! Generate the system :
   STRE gene.str
   ! Fix atoms, if desired :
   CONS FIX SELECT nomov END
   ! Read a coordinate-set and call the energy :
   OPEN UNIT 13 READ FORM NAME  coord_file1.crd
   READ COOR CARDS UNIT 13
   CLOSE UNIT 13
   ! This energy-call is mandatory (use desired non-bond settings) !!! :
   ENER IHBFR 0 ...

   ! Ready to start TReK :
   TREK
      ! Read the reference-coordinates onto which the path will be oriented
      ! (it will not be part of the path) :
      TRAJ READ REFErence   !  NOORient (optional)
        coord_file0.crd
      DONE

      ! Reading the initial-guess path (orient only structures that are
      ! put into the initial path):
      TRAJ READ ORIEnt
        coord_file1.crd     ! this is the reactant structure.
        coord_file2.crd
        coord_file3.crd
      DONE
      ! Adding some points from a dynamics trajectory :
      TRAJ READ NAME  md.dcd  BEGIN 20 SKIP 50 ORIENT
      ! Adding more points :
      TRAJ READ ORIENT
        coord_file4.crd
        coord_file5.crd     ! this is the product structure.
      DONE

      ! Optional: check the energy-profile of the initial path:
      TRAJ ANAL

      ! Starting refinement :
      CPR NCYCLE 100
      ! Writing the path to a trajectory file :
      TRAJ WRITE NAME  path.dcd  RESTART path.cpr

      ! Continue refining, while frequently saving the improved path :
      CPR NCYCLE 100
      TRAJ WRITE NAME  path.dcd  RESTART path.cpr
            .
            .
            .
      CPR NCYCLE 100
      TRAJ WRITE NAME  path.dcd  RESTART path.cpr

      ! Optional: check the energy-profile of the final path:
      TRAJ ANAL
      ! Copy the 3rd point along the path to the main-coordinate set (for ex.):
      COPY ORDER 3
      ! Returning to CHARMM :
   QUIT

   ! Continuing the CPR refinement.
   ! Re-initialize TReK with more space for path-points :
   TREK  MAXPoints 200
      TRAJ READ REFErence  NOORient
        coord_file0.crd
      DONE
      ! Reading the previously saved CPR trajectory (do not orient !!!) :
      TRAJ READ NAME  path.dcd

      ! Re-starting CPR, repeat :
      CPR NECAlls 10000   RESTART path.cpr
      TRAJ WRITE NAME  path.dcd  RESTART path.cpr
            .
            .
            .

      ! Refining 1st-order saddle-points, repeat :
      CPR NECAlls 10000 SADDLE
      TRAJ WRITE NAME  path2.dcd  RESTART path2.cpr
            .
            .
            .

      ! Copying the saddle (if already found) to the main-coordinate set :
      COPY SADDLE
      ! Temporarily returning to CHARMM :
   END

   ! Analyzing the saddle-point:
   IC FILL
   PRINT IC

   RETURN



Numerical and quantum energy potentials !!!
-------------------------------------------
If the energy and/or its derivative is computed numerically (for. ex. quantum
potential or finite difference gradient), some CPR settings (TOLMax,
TOLStep, TOLEne) must be relaxed, otherwise the lack of precision in the
energy/gradient would perturb the convergence of CPR.

``DELTA`` should be increased, so that it is larger than the one-dimensional
finite-difference step used to compute the derivatives (see below for a
description of DELTA).

Example of practical values when using SCF quantum-potentials :

::

   CPR DELTA -1.0e-4 TOLMax 0.1 TOLStep 1.0e-5 TOLEne 1.0e-5

If the SADDle keyword is added (i.e. 1st-order saddle-points are to
be optimized), relax the saddle-point convergence criterion. A recommended
value is SADGrad 0.01 .

Because high-energy regions of conformational space are explored in the
early phase of a CPR refinement, it often happens that quantum-methods fail
to reach SCF convergence.  This can be aleviated by introducing intermediates
into the initial path, so that the linear interpolations between path-points
do not generate structures that are too distorted.  These intermediates do
not need to very good, since they will be removed automatically by CPR later
during the refinment.
Also, it sometimes helps to do the few first CPR cycles (until enough "good"
intermediates are added by CPR itself) with orthogonal (rather than congugate)
line-minimizations.  This is achieved by using the "DELTA 0" setting.
For. example :

::

   TREK
      CPR NECALLS 100  -
          DELTA 0.0 TOL1proj 3.0  FIRStep 0.001 BRKMagn 2.0   -
          EXIT -3  TOLMax 0.15 TOLGrad 0.1 TOLStep 0.001 TOLEne 0.01

      TRAJ WRITE NAME  ugly.dcd
   QUIT

Do not forget to "QUIT" TReK after using these settings (not "END" !), so that
all CPR-settings are returned to their efficient values.


Scanning of the path
^^^^^^^^^^^^^^^^^^^^

CPR searches for energy-maxima along the path (peaks) by "stepping"
along the path and evaluating the energy at every step.
In TReK, the step-size varies along the path. This is a major change
relative to the implementation of CPR in TRAVEL, where the step-size was
constant along the path. In TReK, the step-size is a fraction 1/m
of the length of each path-segment (which varies along the path).

The keyword INTERpol specifies the value of m (m>1).  The default
value is m = 3.  A larger value is more robust and results in a path with
more points, requiring more computations.  m = 2 saves CPU-time (useful
for QM potentials) and is sufficient in most cases.

Boundaries for the step-size can be enforced, with the lowest value specified
with the STPLow keyword (default is 0), and the largest value with the
keyword STEPsize (see below).


Summary of main CPR settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

   NCYCLE = maximum number of CPR-cycles to be performed.

   NECAlls = maximum number of energy-calls to be done.

   SADGrd  = RMS-Gradient in successive line-minimizations leading to a saddle.
             Don't set lower than 10^-3, since the final RMS-gradient will be
   	  much lower anyway due to the many (SADMini) conjugate
             line-minimizations.

   SADMini = Number of successive line-minimizations with
             RMS(gradient) < SADGrd, that qualify a path-point to be
             considered as a saddle-point.  Synonymous to SADCyc.

   DISPlay = Nb. of path-points for which the energy & gradient are printed
             each CPR-cycle (0 by default).  Only useful for small molecules.


Internal CPR-settings  (usually do not need to be tampered with)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

  DELTA    = Finite difference step [in Angstroem] along a path-segment
             for obtaining the 1-dimensional 2nd-derivative used in the
             construction of directions conjugate to that path segment.

             - DELTA = 0 : the set of directions is constructed orthogonally
               to the path segment, not conjugate.

             - DELTA < 0 : the finite difference step is equal to |DELTA|.
               This is the default. When the gradient vector is computed by
	       finite-difference or with QM-energy, make sure that |DELTA| is
	       larger than the finite-difference step used for 1st derivatives
	       or the numerical accuracy of the QM method.

	     - DELTA > 0 : the finite difference step is equal to the distance
               of the scanning point on the segment to the maximized point
               on the segment, so that the value of DELTA is only used as step
               when the line-maximization has failed to move the scanning
               point.  This is the default when ultra-optimizing saddle-points.

  STEPsiz  = Upper limit on the size (in Angstrom) of the steps taken
             when searching for energy peaks along the path.
             By default equals the RMS-distance between reactant and product,
             divided by NGRId+1  (NGRId = 6 by default).
             Of course, if the reactant and the product are the same structure
             (for example in the case of a full group rotation), then NGRId
             makes no sense and the STEPsiz value can be specified explicitly.

  NTANGENT = Number of energy probes on path tangent in search for local max.

  TOL1PROJ = Gradient projection tolerance when refining a path-point.
  TOL2PROJ = Gradient projection tolerance when adding a path-point.

  LINMax   = Maximum number of line minimizations per CPR-cycle.

  FRAME    = Frame-length (in number of cycles) for oscillation-detection
             (maximum allowed value is 20).

  TOLOSCIL = Oscillation tolerance ratio for energy or gradient.
             - TOLOSCIL = 0 : oscillation will not be detected.

  PROJINCR = Increase factor in TOL1prj and TOL2prj when oscillation occurs.

  MAXOSCIL = Maximum tolerated number of oscillation-occurrences.

  REMOVEMOD= m.  If m> 0 , do m line-minimizations for adjacent minimum
             upon point removal. If m <= 0 , add minimum only to avoid looping.


Line-extremization settings.  (usually do not need to be tampered with)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These settings are common to the CPR and SCM commands.  Their default
values have been optimized for the use with analytical energy and derivative
functions, assuming 64-bit precision real numbers (i.e. REAL*8).

If non-analytical energy and/or derivative functions are used, the following
values should be made less stringent (see "Numerical energy potentials").

::

  TOLMAX   = Tolerated gradient, 1-dim. maximizations.
  TOLGRAD  = Tolerated gradient, 1-dim. minimizations.
  TOLSTEP  = Smallest 1-dim. extremization step, in Angstrom.
  TOLENE   = Smallest fractional energy change.

Internal line-extremization settings (do not tamper) :

::

  BRAKETST = Maximal 1-dim. bracketing step, in Angstrom.
  FIRSTEP  = First bracketing step.
  BRKSCALE = Dynamic bracketing Scaling factor.
  BRKMAGNI = Bracket magnification-limit factor.
  EXITMOD  = Exit-Mode: exit with respect to TOLGRA, or also TOLStep & TOLEne.
  LXEVAL   = Max. number of energy-evaluations during each line-extremization.
  ATOMax   = Max. tolerated motion of any atom during a line-maximization.
             When the initial guess for a path is very poor or some of its
             points are very strained, some atoms can experience excessive
             forces, so that they get displaced too much during the
             line-minimizations (some even get pushed out of the structure !).
             To avoid this, ATOMax is set to 0.75 Angstrom during the
             early phases of CPR or SCM refinements.
             If ATOMax = 0, then no limit is imposed on atomic motions.


Obsolete CPR features
^^^^^^^^^^^^^^^^^^^^^

In the previous implementation of CPR called TRAVel, the step-size was not
variable along the path, but was assigned a constant value specified with
the STEPsiz keyword.

For the purpose of backward compatibility, this still can be done by setting
INTERpol < 1, although this is NOT ADVISABLE.

Reduction of STEPsiz (division by SQRT(2) ) then occurs every LOOPred
occurrence of algorithm-looping (unless LOOPred is set to 0).
For ex. if LOOP=4, then the reduced STEPsiz will be 1/2 its original value
after 8 occurrences of looping.  When INTErpol > 1, the default value of
LOOPREd is 0 .


.. _trek_output:

Output
------

TRAJECTORY WRITE !!!
^^^^^^^^^^^^^^^^^^^^

This command writes the path to a file in the format of a standard CHARMM
dynamics trajectory (binary format). This file can then be read with the
"TRAJECTORY READ NAME" command (see :ref:`output <trek_input>`) in order to continue
the refinement, or be treated like any normal CHARMM trajectory, for ex. for
visualization with programs like VMD.


Content of the restart-file
^^^^^^^^^^^^^^^^^^^^^^^^^^^

When writing a path, a CPR restart-file is also written (in ASCII format).
This file contains information about the path, if writing takes place right
after a CPR command. If writing after other commands (for ex. "TRAJ INCR" or
SCM), information such as the energies written to the restart-file are not
pertinent, only the LinMin column (used to flag saddle-points) is meaningful.

The restart-file can be read on the CPR command-line in later TReK sessions
(not with the "TRAJ READ .." command), for ex. when continuing the path
refinement (see :ref:`readrestart <trek_cprcmd>`).

Writing the restart-file does not cost CPU-time, since the data is
taken from the internal data-structures of TReK.  Because the lines of this
file are 160 characters long, it is easier to view it with the "less" utility,
using the -S key:  "less -S restart.cpr" .

The header of the restart-file looks as follows:

::

     * CPR RESTART-FILE (FORMAT 1.1),  CORRESPONDING TO THE PATH IN FILE :
     *                name_of_file_saved_with_traj_write
     *                        title lines ...
     *                            ...
     *
                       THE CPR PARAMETERS WERE :
 SADGRD    SADMini   STEPSZ  REDLUP   DELTAS  LORIEN  PRTOL1   PRTOL2    INTERP
 5.0000E-2      4  9.1980E-2    0   1.000E-07    T    1.000     3.000       3

The title-lines start with an "*", can have any content and appear in unlimited
number.  By default, the second title-line contains the name of the
corresponding trajectory-file, as specified when it was saved.

In addition to the title, the header shows the value that some CPR settings
had in last CPR run.  These header-lines do not start with a "*" and are
mandatory.

After the header, the path is described one path-segment per line.
The first part of the line lists information meaningful to the user.
It is the same as printed by the "TRAJ ANAL" command (described below) :

::

  N   IDX  Lambda  X(N)-Xref  Energy  Grad  LinMin  Curvat


The second part of the line lists internal CPR data. This information is used
by the CPR algorithm when continuing a refinement (see :ref:`scaninfo <trek_cprcmd>`):

::

  StepSiz   LenSeg   StpMax   MaxEnergy   NxMn  PvMn  Up  Dwn

This information describes the path-segment from point (N-1) to N,
i.e. the segment PRECEDING the point listed on that line (the values for
the first path-point are dummy-values) :

::

   StepSiz = The step-size last used to scan the segment.
             Setting it < 0 will force the rescan of that one segment.
   LenSeg  = Length of the segment.
   StpMax  = Number of steps to the peak on the segment, counting from point N-1.
   MaxEner = Energy of the peak on the segment (if any).
   NxMn    = Number of steps to the first local minimum, counting from point N-1.
   PvMn    = Number of steps to the last  local minimum, counting from point N-1.
   Up      = .TRUE. if the energy ends with an increase on that segment.
   Dwn     = .TRUE. if the energy starts with a decrease on that segment.


TRAJECTORY ANALYZE !!!
^^^^^^^^^^^^^^^^^^^^^^

This command re-calculates and print for every path-point the
following values in column-format.
To get an energy profile along the path, plot "Energy" versus "Lambda".
Saddle-points are identified by a negative value in the column "LinMin".

::

   Column:
   .......
   N          = The order of the point along the path.  It is the number that
                is referred to by the SADDLE command, or with BEGIN and STOP keywords.

   Idx        = The identifier number given to each point in the last CPR
                calculation.  It has no meaning for the next runs.  It should
                not be used to identify saddle-points with the SADDLE command.
                Idx is given a negative sign when the structure has not been
                modified (for ex. by CPR) since the last time it was read
                with the "TRAJ READ" command.

                 N
   Lambda     = SUM{ |X(i) - X(i-1)| } / SQRT(3n)  ,
                i=2

                where X(i) is the coordinate-set of path-point i (n atoms).
                This is the total path-length up to path-point N.
                It is the intrinsic reaction-coordinate.  Quantities such as
                the energy, etc., should be plotted as a function of Lambda,
                not N.  Units are Angstrom (see *note lambda: Bugs).

   X(N)-Xref  = The RMS-difference between the coordinates of the point N and
                the reference-structure. Units are Angstrom
                (see *note rmsdiff: Bugs).

   Energy     = The energy (in kcal/mol)

   Grad       = The RMS-gradient of the energy (in kcal/(mol*angstrom)).

   LinMin     = A negative number indicates a tentative saddle-point :

                -SADMini  the path-point satisfies the SADGrd & SADMini criteria
                         for being a saddle-point (see *note analsad: CPRcmd).
                         The saddle-point is sure to be first-order
                         if SADMini was set to 3*n-1 (n = number of moving atoms),
                         and SADGrd was set < .05.

                -M       the path-point has satisfied the condition
                         RMS(gradient)<SADGrd in M successive line-minimizations,
                         M < SADMini, yet it is already flagged as a saddle-point
                         because the numerical machine precision has been reached,
                         so that all requested SADMini line-minimizations could
                         not be performed. Generally is a true saddle-point
                         nevertheless.

                -999999  the point has been declared initially as a saddle-point
                         with the SADDLE command.

                 0       a normal path-point.

                +M       the path-point has satisfied the condition
                         RMS(gradient)<SADGrd in M successive line-minimizations,
                         M < SADMini.  Not yet a saddle-point, but might become
                         one at a later stage.

                +999999  a path-point declared as a saddle-point with the SADDLE
                         command, but another (lower) local energy max. is still
                         present within the path-segment starting at this point.

   Curvat     = The path curvature, as measured by the angle (in degrees) of the
                two path segments joining at this point and divided by half the
                sum of the two adjacent segment-lengths (in Angstrom).
   	     Ex.: zero means the path is linear.

   Grad/Proj  = The ratio of the gradient to its projection onto the plane
                orthogonal to the path at this point. Always > 1 . The larger,
                the better the path follows the gradient (can be improved by
                running SCM for small molecules).

"Curvat" and "Grad/Proj" are not very meaningful for large molecules.

With the SCAN switch, the analysis is extended to points lying along the
path-segments, built by linear interpolation between the explicit path-points.
The size of the interpolation-steps is specified with the STEP keyword
(default is the STEPsize value last used with CPR).


::

   COPY [COMP]

One structure of the path can be copied to the CHARMM main (or comparison)
coordinate set, where it will be available for further manipulation after
exiting the TReK program.

The point can be selected according to its order number N along the path
(ORDER N) or according to its current index (INDEX idx).

It is also possible to select the currently highest saddle-point along the
path by using the SADDLE switch, if such a saddle-point(s) has already
been refined or declared.


.. _trek_tools:

Miscellaneous tools
-------------------

TRAJECTORY INCREASE  or  TRAJECTORY DECREASE
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This command allows to increase or to decrease the density of points along
the path.

If increasing, new points will be inserted by interpolation
between the existing points, so that the distance between adjacent points
corresponds to the specified STEPsize.

Conversely, if decreasing, points will be removed if they are closer
than STEPsize to the point preceding them along the path. This allows to
generate a path with nearly equidistant points, which is required for SCM.

"TRAJ DECREAse" will not delete previously refined or declared saddle-points
(unless the "SADD" keyword is given), it will preserve path-points which are
local minima along the path (unless the "MINI" keyword is given)
and it will not delete path-points if this introduces energy peaks
within the new path-segments (unless the "PEAK" keyword is given).
This allows to use "TRAJ DECR" to remove all unneeded points from a path,
for example before smoothing the path with SCM.

To avoid removing the saddle-points of an old path, declare them either with
the command "SADDLE" or by reading the CPR restart-file corresponding to that
path (see :ref:`decr <trek_input>`).

A sample input-script :

::

  TREK MAXPP 200
    ! Reading the refined CPR-path :
    TRAJ READ NAME   fully_refined.dcd

    ! Flagging the saddle-points previously found by CPR :
    CPR NCYCLE 0  RESTART fully_refined.cpr

    ! Remove as many points as possible (i.e. large STEP values), but
    ! preserve minima & saddle-points and don't introduce energy peaks:
    TRAJ DECREAse STEP 999.9

    ! Get a path with points separated approx. by 1/10 of the
    ! distance between reactant and product, do not preserve
    ! minima or saddle-points and allow new energy peaks:
    TRAJ DECREAse NGRID 10  MINI SADDLE PEAK
    TRAJ WRITE NAME ngrid10.traj  RESTART ngrid10.cpr


    TRAJ WRITE NAME shorter.traj  RESTART shorter.cpr
  QUIT


Miscellaneous
^^^^^^^^^^^^^

From inside TReK, it is possible to use all the CHARMM commands handled
by the SUBROUTINE MISCOM(), such as STREAM, GOTO, LABEL, SET, INCR, DECR, etc.
See :doc:`miscom` for more details.


VERBOSE int
^^^^^^^^^^^
Sets the amount of details being printed out. Values 0-2 are useful for
production runs (default is 2), values 3-6 should only be used for debugging.

The print-level variable PRNLEV is set to 3 when entering TReK,
but this can be overridden by issuing the PRNLevel command.
This is useful to prevent abundant warnings when the structure is very
distorted, which is frequent in the early phase of a path refinement.
Setting "WRNL 0" is also useful in this respect.


.. _trek_scmcmd:

Synchronous Chain Minimization
------------------------------

SCM is an algorithm that smooths the shape of an existing path and
brings it closer to the bottom of the adiabatic valley.
Unlike CPR, SCM does not change the number of points along the path.

SCM
^^^

Before issuing this command, an initial path must be available.
The density of path-points along the initial path must be high enough
to give a relatively continuous description of the shape of the path.
Also, the path-points should be approx. equidistant along the path.

The best use of SCM is to start from a path resulting from a CPR refinement.
Alternatively, the initial path could be generated by a straight
line-interpolation (see the "TRAJ INCRease" command) between the reactant
and product, but this is very inefficient computationally.

SCM works by synchronously minimizing all its path-points,
under the constraint that the points move within hyper-planes
orthogonal to the path. These planes are updated no more often than every
MINUpdate cycles of conjugate minimization (per point), but no later than
after MAXUpdate of such cycles.

The number of conjugate line-minimizations done between plane-updates
on a given point is controlled by the value taken by the angle between the
two path segments joining at that point. This angle varies from 0.0 degrees
(when the path is linear) to 180.0 degrees (when the path is reversing its
direction).  The difference of SCM with other locally updated plane (LUP)
implementations, is that successive conjugate line-minimizations are
continued only as long as this angle is decreasing.  When this angle is
increasing, minimization of that point is stopped once the angle exceeds
a specified value ANGLe, although at least MINUpdate line-minimizations
are done no matter what the behavior of this angle.
When the setting of  ANGLE=180 , SCM becomes equivalent to the LUP method
by Ron Elber.  In that case, and upon extensive minimization, the path
may converge towards a discontinuous state and some path-points may fall
into the reactant and/or the product states.

SCM has converged when the projection of the gradient onto the path is
less than PROJTOL at every point.
This can be quite time consuming, so it is recommended to force an exit
from SCM by specifying the maximum number of global cycles NCYCLE.
In every SCM-cycle, there will be between MINUPDATE and MAXUPDATE
line-minimizations done at each point, so that NCYCLE should be kept small
to allow for periodic saving of the path (see :ref:`scmsav <trek_output>`).

Path-points that are already saddle-points should be declared
(see :ref:`scmdecl <trek_input>`, command "SADDLE"), so that they are kept unchanged
during the SCM refinement.

See :ref:`extrm <trek_cprcmd>` for a description of the line-extremization keywords.


Example input-file
^^^^^^^^^^^^^^^^^^

::

   * Smoothing a refined path with SCM.
   *
   ! Set the desired segment-length in the "equidistant" path,
   ! for ex. 1/30 of the total path-length :
   set finalstep 0.02
   ! Generate the system :
   STRE gene.str
   ! Fix atoms, if desired
   CONS FIX SELECT nomov END
   ! Read any coordinate set and call the energy :
   OPEN UNIT 13 READ FORM NAME  system.crd
   READ COOR CARDS UNIT 13
   CLOSE UNIT 13
   ! This energy-call is mandatory (use desired non-bond settings) !!! :
   ENER IHBFR 0
   ! Start TReK :
   TREK
      ! Read the coordinate-set upon which all path-points will be oriented :
      TRAJ READ REFErence  NOORient
        coord_file0.crd
      DONE
      ! Reading the path previously refined with CPR :
      TRAJ READ NAME  path.dcd

      ! Dummy CPR-command (since NCYCle=0), only to get the saddle-point(s)
      ! declared from the restart-file :
      CPR NCYCLE 0  RESTART saddle.cpr

      ! Trick to make points along the path nearly equidistant (increasing then
      ! decreasing the density of points), so that SCM works better :
      CALC incrstep = @finalstep / 10
      TRAJ INCREASE STEP @incrstep
      ! Keep the saddle-points :
      TRAJ DECREASE STEP @finalstep MINI PEAK

      ! Do a few SCM-cycles:
      SCM NCYCLE 100  ORIEnt
      ! Writing out the path to a trajectory file :
      TRAJ WRITE NAME  path.dcd  RESTART path.cpr

      ! Repeat and save :
      SCM NCYCLE 100  ORIEnt
      TRAJ WRITE NAME  path.dcd  RESTART path.cpr
            .
            .
            .
      ! Get the final energy-profile (after SCM runs, the restart-file holds
      ! invalid energies, only its saddle-point flag LinMin(i) is valid) :
      TRAJ ANAL
   QUIT
   RETURN


Combining SCM and CPR
^^^^^^^^^^^^^^^^^^^^^

The SCM and CPR commands can be given in alternation within a same TReK
session.  This allows to get a path that is both smooth in shape and
has no energy peaks along the path-segments.

It is recommended to run CPR on the path resulting from SCM-runs, because
the SCM-runs may have introduced new energy peaks along the path, which
will be removed (or refined into saddle-points) by CPR.

The restart-file written after the SCM runs still indicates
which points are saddle-points, but it no longer has the proper energy
scanning information.  Therefore, the "SCAN" keyword must be added to
the first CPR command, to insure that the path gets searched for energy-peaks.

CPR and SCM runs can be alternated in this way until the last CPR-run does
not modify the path anymore. This then means that the previous SCM-run did
not introduce any new peaks, so that the resulting path is completed.

.. _trek_sdpcmd:

Steepest Descent Path
---------------------

SDP
^^^

The SDP command provides a carefully controlled descent along the energy
valley, down from a 1st-order saddle-point (obtained from CPR, for example).
This is a rather slow procedure, practical only for small systems.
For large systems, in it recommended to smooth the existing CPR-path with
the SCM method (see :ref:`smooth <trek_scmcmd>`).

Four modes of descent are available :

::

   Mode 1 :   The size of the step down the gradient is reduced until the angle
              between the path and the gradient is less than the value specified
              by ANGLE. This mode is the truest to the definition of the adiabatic
              path, but it also is very slow. It is recommended only for very
              small molecules.

   Mode 2 :   The step is taken along the gradient until the new gradient is
              orthogonal to it (= strict steepest descent).

   Mode 3 :   The step along the gradient is taken as large as possible, as long
              as the energy keeps decreasing (= loose steepest descent).

   Mode 4 :   The steps are those of a strict conjugate-gradient descent.
              While not following the exact bottom of the adiabatic valley
              from step to step, the average path obtained by saving path-points
              after several steps is not significantly worse than the path
              obtained by saving after several steps in Mode 1.
              This mode is the fastest and the one recommended for large
              molecules. If a very accurate adiabatic path is demanded, the
              resulting path can be further improved with the SCM method.
              This mode is the default.

The input to SDP is a chain of points which straddles the saddle.  This chain
can be created with the "CROSS" command, for example (see next).  Up to
NREACTANT and NPRODUCT points are then added to the chain, on the reactant and
product side respectively. The total number of points can not exceed MAXP,
though, so when using SDP, it is advisable to start TReK with a large enough
value for MAXP. A new point is added to the path every time the sum of the
steps taken since the last addition reaches SAVDISTANCE. When the path enters
a region where the energy-gradient is less than specified with MINGRAD,
SDP stops extending the path on that side of the saddle-point, provided that
it already did at least MINCYCLE steps (to allow moving away from the saddle
region, where the gradient is also vanishingly small).
See :ref:`linextr <trek_cprcmd>` for a description of the line-extremization keywords.


CROSsing
^^^^^^^^

This command is used to prepare a minimal path-chain for input to SDP from
a path refined with CPR.  A fully refined saddle-point and two points lying
on each side of it are the required input to CROS. The two surrounding points
serve as initial guess for the saddle-point crossing direction. The output of
CROSS is also a chain with three points, the second of which is the unmodified
saddle-point, while the two surrounding points are located close to the saddle
and near the bottom of the adiabatic valley.  This is simply achieved by
stepping down along the gradient in Mode 1 (see :ref:`descmod <trek_sdpcmd>` for a
description of the descent modes), after having taken the first step
on each side of the saddle-point in direction of the two initial surrounding
path-points.

The keywords of CROSS are :

::

   MINStep int (50)   = The minimum number of steps that are taken, before the new
                        point adjacent to the saddle (on either side) is saved.
   MINDist real (e-4) = The minimum distance moved from the saddle, before the new
                        point adjacent to the saddle (on either side) is saved.
   FIRStep real       = The size of the first step away from the saddle in
                        direction of the initial adjacent path-point.
                        By default, this size is 1/20th of the distance from
                        the saddle to the respective initial point.
   ANGL real (20)     = The angle between the stepping direction and the gradient,
                        which is a setting of the steepest-descent in Mode 1.

By reiterating the CROS command, these points are gradually improved and
yield the reactive mode at the top of the barrier.
CROS is NOT designed to refine a saddle point, which instead must be
provided to it.

If the path provided to CROS has more than three points, it is assumed that
this path is the result of a CPR refinement during the same TReK session and
all points are deleted, except for the highest saddle-point and its two
direct neighbors (provided that CPR had already fully refined the saddle-point,
otherwise the CROS command is ignored).


Sample input-file
^^^^^^^^^^^^^^^^^

::

   * Getting the steepest-descent path.
   *
   ! Generate the system :
   STRE gene.str
   ! Fix atoms, if desired
   CONS FIX SELECT nomov END
   ! Read any coordinate set and call the energy :
   OPEN UNIT 13 READ FORM NAME  system.crd
   READ COOR CARDS UNIT 13
   CLOSE UNIT 13
   ! This energy-call is mandatory (use desired non-bond settings) !!! :
   ENER IHBFR 0
   ! Start TReK with more space to put new path-points :
   TREK  MAXPoints 200
      ! Reading the saddle-point and two neighbors from a fully refined CPR path :
      TRAJ READ NAME  cpr_path.dcd   BEGIN 11 STOP 13
      ! Analyzing the path :
      TRAJ ANAL
      ! Improving the barrier crossing mode :
      CROSsmode
      ! Repeat to improve further :
      CROSsmode
      ! Verify that the 2 surrounding points are close to the saddle-point :
      TRAJ ANAL
      ! Specify the distance between saved points (in Ang.) in steepest descent:
      SDP  SAVDist 0.05
      ! Writing out the path to a trajectory file :
      TRAJ WRITE NAME  path.dcd
      ! Smooth and further refine the path :
      SCM  NCYCle 5
      TRAJ WRITE NAME  path.dcd
      ! Analyzing the path :
      TRAJ ANAL
      ! Returning to CHARMM :
   QUIT
   RETURN

.. _trek_bugs:

TReK known bugs
---------------

Definition of the RMS-difference in coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Normally (and in CHARMM), the RMS-difference between the coordinates of two
structures X1 and X2 is defined as follows (n = number of atoms) :

::

   RMS = | X1 - X2 | / SQRT(n) ,

                             3n
   where | X1 - X2 | = SQRT( SUM{ x1(i) - x2(i))^2 } )
                             i=1

But in TReK, all input and output of any value involving a distance in
conformational space (for ex. LAMBDA or the curvature in "TRAJ ANAL"),
the coordinate RMS-difference is computed as :

::

   RMS'= | X1 - X2 | / SQRT(3*n)       , in other words

   RMS = SQRT(3)*RMS'


Prematurely exceeding MAXP
^^^^^^^^^^^^^^^^^^^^^^^^^^

In some molecular systems, it has happened (very rarely) that nearly identical
points are added repeatedly one next to the other at the same location along
the path, so that the maximum allowed number of path-points NMAXP is reached
and the path-refinement stops.
If this happens, simply continue the refinement after excising
these extraneous points as well as a few additional path-points on each side
along the path.


SADDLE command
^^^^^^^^^^^^^^

A warning is printed after each "SADDLE" command:  Ignore it !  The command
is accepted anyway.  After the CPR command has been issued, path-points
that have been accepted as saddle-points are listed explicitly.


