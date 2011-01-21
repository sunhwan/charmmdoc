.. py:module:: fourd

================================================
4 Dimension dynamics: Description and Discussion
================================================

The energy embedding technique entails placing a molecule into a
higher spatial dimension {Crippen,G.M. & Havel,T.F. (1990)
J.Chem.Inf.Comput.Sci. Vol 30, 222-227}.  The possibility of surmounting
energy barriers with these added degrees of freedom may lead to lower
energy minima.  Here, this is accomplished by molecular dynamics in four
dimensions.  Specifically, another cartesian coordinates was added
to the usual X, Y, and Z coordinates in the LEAPfrog VERLet algorithm.

To employ 4D energy embedding, the energy function and force field
in CHARMM was modified to include fourth dimension coordinates.  An 
additional harmonic energy function has been included to control the 
extent to which a molecule is embedded.  This is quantitatively done by
altering the value of its force constant, initially given by the parameter
K4DI.

The 4D energy embedding procedure can be broken down into three
parts: 4D coordinate generation, relaxation, and back projection.  Fourth
dimensional coordinates can be generated in several ways.  An energy, E4FILL,
in the Fourth dimension can be specified with random coordinates generated
as to sum up to the 4D harmonic energy that a user specifies (i.e. E4FILL 50.0
will give coordinates such that the total sums approximately 50.0 Kcal).  
This method may seem a bit abrupt since a molecule is suddently "thrown" 
into a higher dimension, hence,  molecular dynamics can be used to 
allow a molecule to more slowly obtain fourth dimension coordinates.  
This is done by specifying an initial 4D temperature, FSTT4, with subsequent 
velocities assigned accordingly.  Finally, both these methods may be applied 
simultaneously.  Relaxation involves allowing the molecule to explore the 
potential energy surface and is essentially equilibration.  Alternatively, 
minimization in 4D can be done with the steepest descent algorithm followed 
by 4D dynamics.  Now all that remains is to project this structure back into 
three dimensions.  This last step is thus termed the back projection and is 
achieved by increasing the fourth dimensional force constant linearly
from its initial value of K4DI to MULTK4*K4DI step-wise over the period INC4D
to DEC4D.  This results in a stronger force, confining the 4th dimension 
coordinates to smaller values (i.e. eventually back to 3D).

A problem inherent in the final step of 4D energy embedding is that
"sometimes all projections lead to a bad final conformation" {Crippen,G.M &
Havel,T.F.(1990)J.Chem.Inf.Comput.Sci.Vol 30,222-227}.  Thus, the structure
is rotated into its principal axis of intertia (center of mass) both before
and after its back projection.  When this step is applied the message

::

   ROTATION APPLIED TO PRINCIPAL AXES
   
will appear.  Dynamc4.src is essentially dynamc.src in 4 dimensions.  Note
that even though qeuler still exists in dynamc4.src it has not yet been 
tested.  Also, the usual shake algorithm will only be applied to 
3-dimensional space.

.. _fourd_syntax:

Syntax for the Dynamics Command
-------------------------------

::

   DYNAmics { [LEAPfrog] } VER4 {STRT   } {[TIMEstp real]} [NSTEp integer] -
            { [LANGevin] }      {STARt  } {[FIL4dimension]} 
                                {RESTart} {[SKBOnd]} {[SKANgle]} {[SKDIhedral]}
                                          {[SKVDerWaals]} {[SKELectrostatics]}

            four dimension-spec nonbond-spec hbond-spec frequency-spec -
            unit-spec temperature-spec options-spec 

   hbond-spec::=           updated as in normal LEAPfrog VERLet.
   nonbond-spec::=         updated in 4 dimensions.

   four dimension-spec::=  [K4DInitial real] [INC4Dforce integer]
                            [DEC4Dforce integer] [MULTK4di real]
                             [E4FILLcoordinates real]

   frequency-spec::=       [INBFrq integer] [IEQFrq integer] [IHBFrq integer]
                            [IHTFrq integer] [IPRFrq integer] [NPRInt integer]
                             [NSAVC integer]  [NSAVV integer]  [NTRFrq integer]
                              [ILBFrq integer]  [ISVFRQ integer]
                               [IEQ4 integer] [IHT4 integer] 

   unit-spec::=            [IUNCrd integer] [IUNRea integer] [IUNVel integer]
                            [IUNWri integer] [KUNIt integer]  [CRAShu integer]
                             [BACKup integer]

   temperature-spec::=     [FINAlt real] [FIRStt real] [TEMInc real]
                            [TSTRuc real] [TWINDH real] [TWINDL real]
                             [FNLT4 real] [FSTT4 real] [TIN4 real]
                              [TWH4 real] [TWL4 real]

   options-spec::=         [IASOrs integer] [IASVel integer] [ICHEcw integer]
                            [ISCAle integer] [ISCVel integer] [ISEEd integer]
                             [SCALe real] [NDEGg integer] [RBUFfer real]
                              [AVERage] [ECHEck real] [TOL real]
                               [ICH4 integer]


.. _fourd_description:

Options common 4D dynamics & minimization
-----------------------------------------

The following table describes the keywords which apply to only four
dimension dynamics & minimization.  The remaining parameters are described in 
:doc:`dynamc` and :doc:`minimiz`.

::

   FOURdimensions [INC4d int] [DEC4d int] [K4DI real] [MULTK4 real]  - 
                    [ SKBO ] [ SKAN ] [ SKDI ] [ SKVD ] [ SKEL ] [ SKCO ]  -
                      [FIL4 [E4FILL real ] ]  [ SHAKe ]

=======  =======  ==============================================================
Keyword  Default  Purpose
=======  =======  ==============================================================
INC4D    NSTEP     The step number (specifically, the time in a
                   dynamics run) at which the back projection from 
                   4 to 3 dimensions will begin.  Note the default
                   value of NSTEP will result in no back projection.

DEC4D    NSTEP     The step number at which the back projection from
                   4 to 3 dimensions will end.

K4DI     50.0     The initial force constant for the 4th dimensional
                  harmonic energy term.

MULTK4   1.0      The factor by which K4DI will increase linearly from 
                  INC4D to DEC4D.

FSTT4    FIRSTT    The initial temperature, in the 4th dimension, at which the 
                   velocities have to be assigned to begin the dynamics run.
                   If an equal amount of kinetic energy is needed in all 4
                   dimensions, the default value should be used.  This is
                   because the velocities are all assigned independently in
                   accordance to the initial temperature.

FNLT4    FINALT    The desired final (equilibrium) temperature, in the 4th
                   dimension, for the system.  A final temperature of zero
                   degrees is recommended during a back projection (from 
                   INC4D to DEC4D).

IEQ4     IEQFRQ    The step frequency for assigning or scaling the 4th
                   dimension velocities to FNLT4 temperature during the
                   equilibration stage of the dynamics run.

IHT4     IHTFRQ    The step frequency for heating the molecule in the 4th
                   dimension, in increments of TIN4 degrees in the heating
                   portion of a dynamcis run.

TIN4     TEMINC    The temperature increment to be given to the system every
                   IHT4 steps.  Important in the 4th dimension heating stage.
 
TWH4     TWINDH    The temperature deviation from FNLT4 to be allowed on the
                   high temperature side.  Used only during 4th dimension
                   equilibration.

TWL4     TWINDL    The temperature deviation from FNLT4 to be allowed on the
                   low temperature side.  Used only during 4th dimension
                   equilibration.            

ICH4     ICHECW    The option for checking to see if the average 4th
                   dimension temperature of the system lies within the 
                   allotted temperature window (between FNLT4+TWH4 and
                   FNLT4-TWL4) every IEQ4 steps. 

FIL4               The flag to fill the 4th dimension coordinates.  The 
                   harmonic energy potential of these coordinates will sum
                   to E4FILL.  If not present (recommended), the 4th  
                   dimension coordinates are set to zero and the system will
                   'go into the 4th dimension' as a result of their 
                   initial velocities.

E4FILL    0.0      The total harmonic potential energy from which the initial
                   4th dimension coordinates will be calculated.  Only used  
                   when the flag FIL4 is present.    

SKBO               Flag to skip 4th dimension bond energies (i.e.only
                   compute bond energies in 3 dimensions).

SKAN               Flag to skip 4th dimension angle energies. 

SKDI               Flag to skip 4th dimension proper dihedral energies.

SKVD               Flag to skip 4th dimension Van der Waals energies. 

SKEL               Flag to skip 4th dimension electrostatic energies. 

SKCO               Flag to skip 4th dimension restraint (so restraining Forces 
                   are calculated in 3D only).

SHAKe              Command to place all 4D W's into same W every iteration 
                   (NOTE:energy not conserved).  The 4D forces are not normally 
                   mass weighted, but if SHA4 is used then they are.  Maybe it 
                   should be a 4D option in the future.
=======  =======  ==============================================================

Other Commands:

::

   CONS FIX4 ...     Used in analogy to the FIX command to FIX 4th D coordinates
                     with CONS (meaning one can FIX something in 3D only).


   SCALar FDEQ (0.0) The equilibrium value(s) that the 4th D function will use as
                     the center of the harmonic.  Used for restraining the
                     4th D to non zero values (i.e. forcing a system into
                     the 4th Di).  It should be set with the SCALAR
                     option for individual atoms (if one wants to set different
                     atoms into different 4th D coordinate minima).
                     (1/2)*K4d*W**2, where W=FDIM(I)-FDEQ(I)

   SCALar FDIM (0.0) The coordinate(s) (in analogy to X,Y, & Z) of the 4th D.
                     It should be set with the SCALAR option for individual atoms
                     (if one wants to set different atoms into different 4th D 
                     coordinates).


.. _fourd_recommended:

Recommended CHARMM input for 4d dynamics
----------------------------------------

1) Beginning with a 3d structure and no 4d coordinates, a structure is
   equilibrated in 4d and then back projected (forced back) to 3d.

   ::
   
      DYNAMCS LEAP VER4 START K4DI 50.0 NSTEP 20000 -
       TIMESTEP .001 FSTT4 300.0 FNLT4 300.0 CUTBN 8.0 -
       IHTFRQ 0 IEQFRQ 100 IEQ4 100 NPRINT 10 -
       IUNREA -1 IUNWRI 16 -
       IHBFRQ 25 FIRSTT 1000.0 FINALT 1000.0 TEMINC 0.0 TIN4 0.0

      DYNAMCS LEAP VER4 RESTART NPRE 0 NSTEP 15000 - 
       K4DI 50.0 INC4D 0 DEC4D 15000 MULTK4 10.0 -
       TIMESTEP .001 FSTT4 300.0 FNLT4 300.0 CUTBN 8.0 -
       IHTFRQ 0 IEQFRQ 100 IEQ4 100 NPRINT 10 -
       IUNREA 16 IUNWRI 17 -
       IHBFRQ 25 FIRSTT 1000.0 FINALT 100.0 TEMINC 3.0 TIN4 1.0

2) Beginning with a 4d structure with 10.0 Kcal initially in the 4th 
   dimension.

   ::
   
      DYNAMCS LEAP VER4 START K4DI 50.0 NSTEP 20000 -
       FIL4 E4FILL 10.0 -
       TIMESTEP .001 FSTT4 300.0 FNLT4 300.0 CUTBN 8.0 -
       IHTFRQ 0 IEQFRQ 100 IEQ4 100 NPRINT 10 -
       IUNREA -1 IUNWRI 16 -
       IHBFRQ 25 FIRSTT 1000.0 FINALT 1000.0 TEMINC 0.0 TIN4 0.0

3) Fixing the 4th D coordinates of some bulk solvent and setting the
   solute coordinates "out" in 4D space and along with its equilibrium
   value.  Following this the energy is determined..
   
   ::
   
      CONS FIX4 SELE SEGID BULK END
      SCALAR FDIM SET 10.0 SELE SEGID SOLV END
      FOUR K4DI 50.0 SKBO SKAN SKDI SKCO
      ENERGY

