.. py:module:: mbond

==============================
Multi-body Dynamics:  Overview
==============================

In multi-body dynamics, aggregates of atoms are gathered into
"bodies".  For a dynamics run, the system comprises one or more bodies
and zero or more atoms which are not part of any body.  By gathering
the atoms in this way, the total number of variables in the system is
considerably reduced which is expected to significantly improve the
computational performance.  Furthermore because such a simulation aims
to reproduce the characteristic (i.e. low-frequency) motion of the
system, relatively long time steps are possible.  The final advantage
of this scheme is that bond-lengths may be explicitly constrained
(between bodies and in the atomistic regions) in a computationally
efficient manner.

For detail description of MBO(N)D method,  refer to the paper of 
"MBO(N)D: A multibody method for long-time molecular dynamics simulations"
in Journal of Computational Chemistry, Vol 21, 1 (2000).

There are two steps needed before starting a dynamics run

1. Identifying the bodies and the atomistic regions and
	
2. Generating (or loading pre-calculated) modes for each of the
   bodies. 

All of the standard CHARMM output and analysis mechanisms work with
multi-body dynamics.  One new file format is used (to store computed
mode shapes).

The MBOND command is used for setting up the system and controlling
its activation. It can be used either as a single line command (mainly for
control and status reports) or as an opening for a command block (for
setting up the system substructuring and mode assignments). 

All the single line commands can be also used from within the command block. 
The single line commands are:

::

   MBOND   { [ON    ] }
   	{ [OFF   ] }
   	{ [CLEAr ] }
   	{ [STATus { [SHORt               ]} ]}
   		  { [LONG  [BODY char*4] ]}
   		  { [VLONg [BODY char*4] ]}
   	{ [ELECtrostatics]}
   	{ [ATOM   { [SCBOnd real] } ]}
                     { [NOBOnd     ] }
   	{ [STIFF ] }

The command block has the following syntax:

::

   MBOND [NBODY integer]
   	SUBStructure
   	.
   	.
   	END
   	MODES
   	.
   	.
   	.
   	END
   END

NBODY is used to specify the MAXIMUM number of bodies to be supported.
It only needs to be specified once per CHARMM run.  To minimize memory
usage, it is useful to keep this number small.  If you change the
value of NBODY within the same CHARMM run, then all previously defined
structures are erased.  NBODY must be a positive number.

Within the MBOND command block all of CHARMM's miscellaneous commands
can be used. This allows for command-line control: STREam files, GOTO,
LABEl, IF etc.    


.. _mbond_substructuring:

Sub structuring 
---------------

To identify collections of atoms to be grouped into bodies, use the
``SUBSTRUCTURING sub-comands`` of MBOND.  A body is identified using the
SELECTION commands of CHARMM.  Generally bodies will be contiguous
atoms, though there is no restriction that they be so (e.g. you could
define a set of water molecules as a body).

Bodies can have any name of four characters, except 'ALL '.

By default, the bonds between bodies and all bonds in the atomistic
regions will be constrained.  To change this behavior, use the HINGE
subcommand (currently only BOND and OFF are available).

It is possible to specify that all atoms in the system should be
modeled as separate bodies i.e. do an atomistic simulation using the
multi-body formulation.  This is the function of the PARTICLE
sub-command.

The syntax for the SUBSTRUCTURE sub-command is:

::

	MBOND  [NBODY integer]
 
	    SUBStructure

	   	BODY   char*4  SELEct (...) END    [SHOW]

                MAKEBODY SELE (...) END            [SHOW]

        	HINGE { [BOND]} [ANGLe] [DIHEdral] [OFF]

		[ PARTICLE ]

    	    END
	END

These options are documented on the next screen.


.. _mbond_subdesc:

Syntax Description
------------------

* ``SUBStructure``

  Begins the substructuring commands. When you exit (END
  subcommand), the "topology" information is passed to MBOND. 


* ``BODY char*4  SELEct (....) END  [SHOW]``

  Define a BODY using the CHARMM regular SELECT facility. The body
  does not have to be continuous in terms of atomic numbers.
  The "body_name" can be any character*4 string.
  
  SHOW flashes out information about the body just defined.

* ``MAKEBODY SELEct (....) END [SHOW]``

  Define bodies using the CHARMM regular SELECT facility. 
  i.e. 
  ``Makebody sele resname tip3 end show``
  put each residue name "TIP3" to be in each separated body.

* ``HINGE { [BOND ]} [ANGLe] [DIHEdral] [OFF]``

  Define the type of constraints between bodies. Currently only
  BOND constraints are available. 

  options:
  
  ======== ====================================================
  BOND     constrain all bonds 
           (i.e. except those within a body)  (default)
  ANGLe    constrain hinge angles (currently not available)
  DIHEdral constrain hinge dihedrals (currently not available)
  OFF      NO constraints.
  ======== ====================================================
  
* ``PARTICLE``

  Each atom in the system is to be considered as a separate
  body.  The MBOND formulation will be used for dynamics
  simulations.  This is sometimes useful for direct comparison
  with CHARMM atomistic simulations.

.. _mbond_modes:

Modes
-----

If you simply select sets of atoms as bodies, then these will be
modeled as rigid bodies by default.  To identify a body as flexible,
a set of modes can either be generated or loaded for it.  A new file
format (see ) has been created with all necessary information for a
multi-body dynamics simulation.  In addition, modes can be sorted
(either by frequency or by delocalization).

The syntax for the MODE sub-command is:

::

   MBOND  [NBODY integer]

    	MODEs { [KEEP  ] } [NMODe integer]
   	      { [DELEte] }
                 {[HARMonic] [FCON] [RCUT] [VDWR]} Optional-selection (see below)

   		LOAD     char*4 {[NMODe integer]             } UNIT integer
   			        {[MODE  integer THRU integer]}   [SHOW]  [NOVEC]

   		UNLOad   [char*4 ]
   			 [ALL ]

   		USEM     {[char*4]} 	{[NMODe integer]             }
   			 {[ALL]}       	{[MODE  integer THRU integer]}

   		GENErate char*4 [NMODe integer] [DELO] [SHOW]

   			MINI minimiz_spec

   			DIAG  	{[VACUum]}
   			  	{FIXEd}

   			SORT	{FREQ}		[SHOW]
   				{DELOcalization}

   			WRITE {[NMODe integer]		  } UNIT integer
   			      {[MODE integer THRU integer]}
   			      {[ALL]}	
   		END

   		SORT	{[char*4]} 	{[FREQ]}
   			{[ALL]}		{[DELOcalization]}

   		PED {[char*4]} [Mode-Spec] [Magnitude Spec] [TOL real]

   	    END

   ***Optional-selection { [MULTiple] [BOND] [ANGLE] [ALL ] [VDWF] }

These options are documented on the next screen.

.. _mbond_modedesc:

Mode Syntax
-----------

::

   MODES { [KEEP  ] } [NMODe integer]
         { [DELEte] }
          { [HARMomic] [FCON] [RCUT] [VDWR]} Optional-selection (see below)

   ***Optional-selection { [MULTiple] [BOND] [ANGLE] [ALL ] [VDWF] 

Invokes the flexible-body mode-definition module. All multi-atom
bodies (defined in the substructuring module) are assumed "rigid"
unless explicitly defined as "flexible" in the MODES module (by LOAD,
USEM or GENErate).  

By applying "HARMonic" command after "MODES", the
simple harmonic potential, instead of the regular CHARMM potential,
will be used to generate modes for the body. (For more detail of the
harmonic potential mode, see :ref:`harmonic`

options:

======== ==========================================================
KEEP     (the default) keeps the mode information (eigenvectors,
         frequencies, nominal coordinates) in CHARMM memory after
         passing it to MBOND.  
DELETE   removes this information after interfacing with MBOND.
NMODE    sets a default value for the NMODE keyword in the LOAD
         command. Default is 20.  
HARMonic invokes the harmonic potential mode generation for the
         flexible-body.
======== ==========================================================

Following commands works only if "HARM" command is defined.

======== ==============================================================
FCON     defines a single force constant of the harmonic potential
         to be used for mode-generation. (No need to be defined when
         MULT command is on.)
RCUT     defines the cut-off distance for pairing atoms. Only distance
         criteria is used to define a pair of atoms.
VDWR     sets an initial cut-off distance to be a sum of van der Waals
         radius of interacting atoms.
         
         Final Cut-off distance = RCUT + VDWR-atom1 + VDWR-atom2
======== ==============================================================

**Optional-selection Commands for Harmonic potential modes**

======== ===================================================================
BOND
ANGLE    
ALL      These options overide the distance criteria for defining
         pair of atoms. BOND command forces to make bonded atoms as
         a interaction pair. ANGLE command forces to connect
         atoms which have a common atom center. ALL command will invoke
         both BOND and ANGLE command. The remaining of atoms are
         paired by the distance criteria.

MULTIple forces to use the different force constant based on the
         nature of interactions. Each harmonic force constant is
         obtained from the CHARMM parameters, such as bond, angle,
         Lennard-Jones (VDW) parameters.

VDWF     In MULT selection, the force constant (K) between paired atoms,
         that are connected by the distance criteria, is defined by the
         second derivative of the Lennard-Jones potential.

         :math:`K = \frac{d^2 V}{d r^2}` at where r is the given distance.

         IF VDWF is on, r is the potential minima in the previous equation.
======== ===================================================================

:: 

   LOAD  char*4 {NMODe integer             }  UNIT integer 
         {MODE  integer THRU integer}    [SHOW]  [NOVEC]

Load the mode-information of the body  with "body_name" character*4 
from the formatted file opened in UNIT (UNFOrmatted option not
available).

The LOAD command reads the body-based nominal coordinates, the
"modal stiffness" matrix (F^-1 K F), and either the first NMODE
eigenvectors or the eigenvectors in the range "MODE ... THRU ...".

If this body is already defined as "flexible" a call to LOAD
will overwrite the previous mode-information (frees the
currently allocated space and reallocates a new one).

options:

===== ==========================================================
SHOW  echoes the information read from the file.

NOVEC a flag that allows for reading the modal-stiffness
      matrix and properties without reading the eigenvectors.
===== ==========================================================
 
*An appropriate unit has to be opened before using the LOAD
command.*

::

   UNLOad   {[char*4]}
      {[ALL]}

Frees the space where the mode information for the body  with
"body_name" character*4 is stored, effectively turning it to a
"rigid" body.   

Can be used for changing a "flexible" body into a "rigid" body, and for
explicitly "cleaning up" before changing the character of an already
defined "flexible" body.

UNLOAD is identical to: "LOAD char*4 NMODE 0", and is not
necessary when LOADing or GENErating new modes for a "flexible" body.

UNLOAD ALL removes from memory all mode information for all bodies.

::

   USEM     {[char*4]}     {NMODE no_of_modes     }
      {[ALL]}     {MODE ifirst THRU ilast}

Use the first NMODE modes or MODE/THRU modes in the dynamics
calculation. The default is all the loaded modes. The number
of modes used cannot exceed the number    of modes loaded by
LOAD (or GENErate).  

The USEM command does not free memory so that the loaded modes can be
reached again if needed. 

USEM ALL ... will use the specified mode numbers for all bodies.

::

   GENErate char*4 {[NMODe integer]} [DELO] [SHOW]
         {ALL}

         MINI minimiz_spec

         DIAG  {[VACUum]}
            {FIXE}

         SORT  {FREQ}      [SHOW]
            {DELOcalization}


         WRITE {[NMODe integer]       } UNIT integer
               {[MODE integer THRU integer]}
         END
   
Generate calculates modes for the specified body.  At most,
NMODes (default 20) are generated.  ALL modes for that body
may be determined if sufficient memory is available.  Any
previously loaded (or generated) modes for that body are cleared.

The DELOcalization is a measure of how localized a particular
mode is.  It is the ratio of the fourth moment of the
displacement along the mode vector and the second-moment
squared, i.e., 

Sum_over_atoms (r^4/(r^2)^2).

In very localized modes this ratio is close to 1,
in very delocalized modes it is close to 0.

.. note::
   NOTE, at very high frequencies multiple occurrences of
   localized modes (e.g., C-H stretching), will have a small 
   "DELOcalization" value.

Selecting this keyword computes the DELOcalization of each
generated mode and stores it as the second property of that
mode (the frequency is the first). 

All existing MINImization options are suppored (by the
existing CHARMM minimization module).  All atoms which do not
belong to the     specified body are fixed for the duration of the 
GENErate command (all atoms are restored to their initial
coordinates on exit).

DIAG supports two different methods.  In both, it computes 
the Hessian for atoms in the selected body, and diagonalizes
it. The 6 translational/rotational modes are discarded. 
FIXED environment includes the interactions between atoms in
this body and other atoms in the system when computing the 
Hessian.  VACUUM (the default) ignores these external
interactions .

By default, generated modes are SORTed by frequency (most
negative to most positive).  You can change the order by 
SORTing according to the delocalization factor (and re-sorting
by frequency). 

WRITe saves the range of selected modes to the specified UNIT
in the format which LOAD uses.  They are saved in the current
order (i.e. sorted by frequency or by delocalization).

All of CHARMMs control and command parsing options are
supported within the generate block.  In particular, this
allows for different stream files to be created and used to 
generate modes for each of the bodies in a system.

::

   SORT  {[char*4]}  {[FREQ]}
      {[ALL]}     {[DELOcalization]}

Modes can be sorted at other times other than GENEration.  For
example if modes are computed using some other program, then
loaded into CHARMM, they can be still be SORTed.  The possible
keys are FREQuency and DELOcalization.  If the DELOcalization
has not already been computed, it will be automatically.

::

   PED   {[char*4]} [Mode-Spec] [Magnitude Spec] [TOL real]

For a given magnitude specification, it computes the expectation
value for the energy contribution change for each internal
coordinate term (bond, angle, dihedral, and improper dihedral)
and prints that term if the fluctuation is greater than the
tolerance (default TOL 0.0001).

[Mode spec] is the usual choice of either NMODES N (i.e. the first N
modes) or MODE M THRU N.  See :ref:`Normal Modes <vibran>`
for a description of [Magnitude Spec].

Modes must have been loaded for the specified body (either by LOAD or
GENERATE).  Right now, only VACUUM modes are supported (i.e.
contributions from  atoms outside the body are NOT included).


.. _mbond_other:

Other
-----

::

   MBOND   { [ON    ] }
      { [OFF   ] }
      { [CLEAr ] }
      { [STATus { [SHORt               ]} ]}
           { [LONG  [BODY char*4] ]}
           { [VLONg [BODY char*4] ]}
      { [ELECtrostatics]}
      { [ATOM   { [SCBOnd real ]} ]}
                     { [NOBOnd      ]}
      { [STIFfness]}

All these single line commands can be also used from within the MBOND
command block (but not from within the SUBS and MODE sub-blocks

====== =================================================================
ON     Activate MBOND and take the substructuring into
       account when calculating energy etc. This is the
       default after setting up the substructuring (in the
       MBOND command block). 

OFF    Un-activate MBOND (all data structures are kept intact).

CLEAR  Un-activate MBOND and remove all related data
       structures from memory.  Also resets ELEC, ATOM etc.

STATus Printout MBOND status report (including number of
       bodies, size, mode information, etc). For the LONG and
       VLONG options one can specify a specific BODY name,
       the default is all bodies.  You can also specify ALL.

       option:
       * SHORt    - give global status
       * LONG (default)  - adds body information
       * VLONg       - adds mode information

ELEC   Use MBOND's multipole expansion algorithm to calculate
       electrostatic interactions (skip CHARMM's electrostatic
       energy evaluations). The default is to use CHARMM's
       electrostatics (!).  CURRENTLY DISABLED.

ATOM   Intra-body forces can be computed either explicitly or
       by employing a linear force approximation (using the
       STIFFNESS MATRIX).  The default is the full atomistic
       calculation specified by the ATOM option.  

       To compensate for Cartesian curvelinear effects the SCBOND 
       option may be invoked. This will scale all the intra-body 
       bond-energy terms by the factor <real>. The scaling factor 
       should be in the range [0.0,1.0] (where 0.0 is identical to 
       NOBOnd, and 1.0 is identical to the regulat MBOND ATOM option). 
       Inputs outside of this range are reset to 1.0).
       The NOBOND keyword will skip intra-body bond-energy terms 
       altogether (equivalent to SCBOnd 0.0).

STIFf  Use a linear force approximation to compute the
       intra-body forces (using the STIFFNESS MATRIX).
       This is the opposite of the ATOM option.
====== =================================================================

.. _mbond_dynamics:

Dynamics
========

Multi body dynamics are invoked by using the option MBOND on the
regular CHARMM DYNAmics command.  This enables a number of keywords
particular to MBOND and disables some which are not.  In general,
however, all the standard DYNAmics keywords are supported unless
specifically mentioned here.

Bodies must have been defined prior to invoking multi-body dynamics
(see MBOND SUBStructure for details).  If no modes have been loaded
(or generated) for the bodies, they are assumed to be rigid.

Control is passed from CHARMM to the MBO(N)D package until the end of
the dynamics simulation.  During the run, CHARMM is utilized to
evaluate the forces.  It also implements some of the controlling logic
for heating, equilibration, thermostats etc.

Only velocity scaling is supported for heating and equilibration
protocols.  However the initial velocities may be assigned using any
of the CHARMM options (including using a RESTART file).  The thermostat
which is available is a simple Berendsen method. 

The general protocol is to heat a system using an atomistic simulation,
then equilibrate it and generate modes.  When the changeover is made
to a multi-body simulation, the system needs to be equilibrated once
more before production dynamics can begin.

Multiple Time Scales are supported in MBOND.  There is a special MTS
keyword within the MBOND command to enable/control Multiple Time Scale
dynamics. See :ref:`mbond_mts` for details.
 
The following CHARMM features are NOT supported in multi-body
simulations: SHAKE, CONSTANT PRESSURE, NOSE.  In addition,
the integrators the MBOND supports are Lobatto, Velocity Verlet,
Velocity Central Difference and 4th Order Runge Kutta.  There is no
LEAPFROG or 4D Verlet.

These MBOND dynamics options are documented on the next screen.

.. _mbond_dyndesc:

DynDesc
=======

::

   mbond-spec::=  [[ {LOBAtto}   ]           ]
         [[ VVER     ] [NCYC integer] [VTOL real]  ]
         [[ VCD      ]           ]

         [[ RK ]                 ]

         [ MBPRLev integer ] [ VECFreq integer ] [ TFRQ integer ]
           [ IUVEctor integer] [ ATOM ]

The following four integrator choices are available for multi-body dynamics

======== =======================================
LOBATTO  The default integrator.
VVER     Iterated velocity-verlet integration.
VCD      Verlet Central difference
RK       4th Order Runge-Kutta
======== =======================================

======== ======== ================================================================
Keyword  Default  Description
======== ======== ================================================================
NCYC     1        The maximum number of iterations in any verlet style
                  integration step  (default 1)

VTOL     1E-10    Convergence criterion for verlet style integrators
                  in multi-body dynamics (default 1E-10)

ATOM              If this is not specified (or has not been
                  previously), the linear force approximation is used
                  within bodies.  If specified, the full force
                  computation is performed.

MBPRlev  0        Print level for multi-body dynamics.

                  ======== ===== ================================================
                  MBPRlev   -1   No printout from within MBOND, only the
                                 data passed to CHARMM will be printed out
                                 (every NPRInt timesteps).
                  MBPRlev    0   Basic printout from within MBOND.
                  MBPRlev    1   Print regular (un-substructured) CHARMM energy
                  MBPRlev    3   Full printout from within MBOND.
                  MBPRlev    9   Extended (debugging) printout from within MBOND.
                  ======== ===== ================================================

VECFreq  NPRINT   Printout frequency for MBOND vector information
                  (MBPRLEV at least 3).

IUVEctor 6        Unit to write MBOND debugging information.

TFRQ     NPRINT   Frequency for writing out the temperatures.
======== ======== ================================================================


.. _mbond_mts:

Multiple Time Scale (MTS)
=========================

Multiple Time Scale dynamics are a natural complement to body based
dynamics.  MBOND currently supports segregation of the computation of
non-bonded forces into different timescales.  By default, the number
of bins and the relative ratios are computed automatically (although
these may also be specified exactly).

Multiple time step dynamics in MBOND is controlled by the MTS keyword
to the MBOND command.  This must be invoked BEFORE dynamics is begun.
Once you have selected options in this way they remain in effect for
all subsequent dynamics runs. 

::

   MBOND
      MTS  
          MASS [CALIbration f]
          DIST
         LINEar  A B
         CUTOn R
                   SLFG BIN #1 RSCUT #2 BUFF #3
                         ! #1 - bin number
                         ! #2 - Cut off distance in Angstrom
                         ! #3 - Buffer region in Angstrom
          END
          RATIos K L M..... ! 1 < K < L < M < ...
                       ! K,L,M specify the different
                                 ! update frequencies
          MAXStep t
               ADD ## SELE (...) END  ! ## - bin number
          CLEAr
      END   

The MTS keyword signifies that subsequent dynamics runs are to utilize
multiple time steps.  This keyword sets up the necessary data
structures to support multiple time steps.  Dynamics is only initiated
with the DYNA MBOND command.

The number of stages for multiple time steps can either be specified
by the user (implicitly through the RATIos keyword) or computed
automatically.  If it is not provided, it will always be computed
automatically.  In addition, the optimal frequencies for updating the
various stages will be automatically computed if not provided.

The MASS keyword specifies that the "interaction mass" is used to segregate
the non-bonded force.  The interaction mass depends on the
substructuring employed: for an atom in an atomistic region, it is the
mass of that atom; for an atom in a rigid body, it is the combined
mass of all atoms in that body; for an atom in a flexible body, it
is currently the mass of the atom.

When examining atom pairs in the construction of the non-bonded list,
the quantity SQRT((M1*M2)/M1+M2)) [where Mi is the interaction mass of
the i'th atom; which is the mass of the body if the atom belongs to a rigid
body, or is the mass of the atom itself if it belongs to a flexible body
or is a particle] is used to determine at what rate the non-bonded forces
for this pair will be computed.  This ensures that interaction involving
atoms in bodies will, in general, be computed at a slower rate.

It is possible to calibrate the interpretation of this quantity (use
the CALIbrate keyword).  The default value is 0.5 and typically values
are close to the fastest required update frequency (bearing in mind the
timestep which is specified later in the dynamics run).

If you use the DISTANCE keyword (recommended), the interaction mass
is multiplied by a DISTANCE factor C(r) (where r is the distance
between the two atoms).  The form of the function C is specified by
the sub-options

::

   CUTON r0 C(r) = 1    r < r0
   LINEAR A B      C(r) = A*(r - r0) + B   r >= r0

Usually r0 ~ 4A.  A must be positive and B should be 1.
The DISTance keyword implies the MASS keyword.

In addition, you can use the capability to select distance-dependent 
binning (Similar to the CHARMM MTS - refer to :doc:`mts`) by the keyword, SLFG.
User can select only a distance-depenent binning, instead of mass-weighted 
one described previously. By using this method, nonbond forces are 
calculated by using switching functions in order to have smooth transitions 
of forces. This method currently works with atom-based cut-off only.

Also by using the keyword, SLFG, you can invoke the capability to distribute
only internal motions (bond, angle, dihedral, etc.) into the fastest bin.
If user uses "SLFG BIN 1 RSCUT 0.0 BUFF 0.0", user distribute all internal
force calculations to the 1st bin. This can help to get the better stability 
for hinge-off simulation.

If you wish to explicitly define both the number of bins and the
different non-bond computation rates, use the RATIos keyword.  Up to
four (integer) values may be provided.  The fastest rate is always 1
and the values provided must be monotonically increasing and greater
than one.  The number of values provided determines the number of bins
i.e. if one value is provided, the non-bonded interactions will be
separated into two stages; if 3 values are provided, 4 stages will be
used.

MAXStep is used to limit the biggest timesetep that is taken when
automatic ratios are computed.  This is useful if you are using a
large value for the timestep in the DYNAMICS command and the automatic
computation may be producing huge net timesteps for the slowest update
rate. 

The ADD keyword is used to specify interaction, which involve selected atoms,
into the user specified bins in MTS. i.e. ADD 1 SELEct RESNAME TIP3 END -
This means that all interaction involving "tip3" residues are place into 1st
bin.

CLEAR is used to disable body based multiple time steps.  All ratios
(explicitly provided or automatically generated) are cleared.


NOTES:
1) A valid substructuring must have been defined before the MTS keyword
   can be used.

2) Update frequencies in dynamics (non-bonded list generation, reporting
   etc) are modified to be multiples of the lowest common multiple of the
   MTS frequencies.

3) The MBOND Multiple Time Step only segregates non-bonded interactions.
   Other force computations (e.g. dihedrals etc) are all computed at the
   fastest rate. (We will be adding other keywords extend this).

4) There's not much benefit to using the MASS keyword by itself.

5) In theory, the interaction mass of atoms in a flexible body can be
   larger than the mass of the individual atom.  But in practice, the
   benefits of doing this are not very great.

6) If you already have a large base timestep (15 or 20 fs), it is best to
   just specify the exact ratios rather than computing them
   automatically.

7) When specifying exact ratios, bear in mind that the LCM will determine
   the non-bond update frequency.  So 
   
   ``RATIo 2 3 4 6        LCM = 12``
   
   may be a better choice than 
   
   ``RATIo 2 3 4 5        LCM = 60``


EXAMPLES

::

   MBOND
      SUB
        (substructuring commands)
      END
      MTS            ! Enable Multiple Time Step Dynamics
         MASS CALIB .707      ! Use interactio mass with a factor of sqrt(2)
      END            ! The number of bins and the frequencies
   END               ! have not been specified and will be 
               ! computed automatically

   MBOND
      SUB
        (substructuring commands)
      END
      MTS            ! Enable Multiple Time Step Dynamics
        CLEAR        ! Eliminate previously assigned parameters
        DIST
      LINEAR 2
      CUTON 4
      END
        RATIo 2 4 6     ! 4 stages of force computation
      END


   MBOND
      MTS
        RATIO 2  5
        DIST
           SLFG BIN 1 RSCUT 5.0  BUFF 1.0 ! Bin #1 has the interactions 
                                          ! between 0 and 5A
           SLFG BIN 2 RSCUT 8.0  BUFF 1.0 ! Bin #2 has the all interactions 
                                          ! between 5 and 8.0A 
                                          ! Rest of interactions are in Bin #3
        END
      END
   END


   MBOND
      MTS 
        MASS CALI  0.5
        RATIO  3  6
        DIST                        ! 3 stages of force computation
           Linear 1.0               ! but all interactions involving TIP3
           CUTON 4                  ! residues are place into 1st bin.
        END
        ADD 1 SELEct RESNAME TIP3 END 
      END
   END

.. _mbond_langevin:

Langevin
========

Langevin Dynamics is invoked for multibody dynamics by specifying the
LANGEVIN keyword in the DYNAMICS MBOND command.  The usual CHARMM
control parameters, ISEED and TBATH, are read (the defaults are 314159
and FINALT (which itself defaults to 298.0)).

Friction coefficients for the forces are specified separately using
the MBOND LANGEVIN keyword.  The available options are

::

   MBOND
      LANGEVIN
        BODY {Name} COEF real [real, real,... real]
        BODY ALL COEF real
        ATOM ALL COEF real [real, real]
        ATOM SELEct (standard selection) END COEF real [real, real]
        LIST
        CLEAR
      END
   END

Currently, a different coefficient can be specified for EVERY degree
of freedom that a body or particle has.  So for each body, there are 3
translational, 3 rotational and 1 for each mode USEd for that body.
To simplify things, a single coefficient may be specified
(implementation detail - this is replicated 6+nModes times) for a
body.  If you specify more than one, you MUST specify exactly
6+nModes. 

Likewise, each particle has 3 degrees of freedom and either 1 or 3
coefficients may be specified.  You can specify ALL particles or
select them using any CHARMM selection sequence.  Likewise ALL bodies
can be done at once or you can deal with one at a time.

1) It is OK to specify values for everything using ALL and then
   overwrite the specific values for a specific set of particles or a
   specific body.

2) As bodies may have different numbers of modes loaded, currently you
   can only specify a single coefficient when using ALL bodies.

Specification of the coefficients obviously requires that
substructuring has been defined AND that the appropriate mode
specifications made.  If you subsequently change any of this, then
there is no way to track how the coefficients map to the new scheme
and they are ALL flagged as invalid.

Some examples

::

   MBOND
     SUBST
        STREAM blah.str   ! Defines substructuring
     END
     MODES
        OPEN READ UNIT 16 modes.mod
        LOAD bdy1 NMODe 4 UNIT 16
     END
     LANGEVIN
        BODY ALL COEF 0.8
        BODY bdy1 COEF 0.5 0.5 0.5 0.6 0.6 0.6 1.1 1.2 1.1 0.9
        ATOM ALL COEF 1.0
     END
   END

   DYNAMICS MBOND  -
      LANGEVI ISEED 3250 TBATH 200 -
      LOBATTO NCYC 3 -
   etc


.. _mbond_output:

The output of a multi-body dynamics simulation includes almost all of
the same information generated during an atomistic simulation,
together with additional information pertinent to the body based
approach.

The standard dynamics files (trajectory, velocity, energy and restart
files) are all supported during multi-body dynamics (there are no
changes in format).  Energy files have an additional line of BODY
specific information per entry.  Trajectories can be analyzed as usual
and runs restarted if need be.

NOTE that currently RESTART fits the current coordinates and
velocities, the available modes/modal amplitudes and the coordinates
at which the modes were generated.  Care should be taken to check that
this fitting process has not significantly modified the trjaectory.

The energy reporting, controlled by NPRINT, has an additional line of
information each step, which specifies the following terms:

* Deformation Energy
* Electrostatic energy as computed by MBOND's HBMA
* Momentum
* Aggregate Body Temperature
* Temperature in the atomistic regions

The overall temperature and that of the body and atomistic regions can
be separately controlled by the TFRQ option in the DYNAmics command.
In addition, this will report the individual kinetic energy of each
body.

MBPRLEV controls the reporting of the following, 

::

   MBPRLEV = 0
      Modal and Kinetic Energy for each body
      Linear and Angular Momentum for each body

   MBPRLEV = 1 - The above plus
      Convergence rate during initial fitting

   MBPRLEV = 2 - The above plus
      All initial input data for each body and particle

   MBPRLEV = 9 - The above plus
      All sizes for the multibody dynamics arrays
      Complete description of the initial data
      Modal Amplitudes, body velocities and angular velocities

.. _mbond_harmonic:

Body-Based Mode Generation by Harmonic Potential
================================================

1) Introduction

   Generating reasonable body-based modes is a crucial part of MBO(N)D. Those
   modes are used for modeling flexible bodies. In recent years, there have 
   been several research papers in which people described a simple approach
   to generate reasonable thermal fluctuations of proteins. In the approach,
   a single parameter harmonic potential was used instead of the more complicated
   force field. Based on their observations, low-frequency mode-vectors 
   generated by this simple potential was in good agreement with those calculated
   from the full potential.

   From the MBO(N)D point of view, this approach for generating body-based 
   modes is attractive. The instantaneous structure can be treated as an
   equilibrium structure in this approach so that no artificial energy minimization
   is necessary. Other approach, using the full CHARMM poteintail, requires 
   small steps of energy minimization before generating modes because 
   inappropriate modes produced from slight distorted peptide-plains. In the same
   vein, no negative frequency eigenvectors are generated. In the regular 
   MBO(N)D mode generation approach, several of the body-based modes are usually
   of negative-frequency. However, harmonic potential approach will eliminate 
   those problems altogether.

   References:
   
   1. M.M. Tirion, Phys. Rev. Letts., 77, 1905 (1996)
   2. I. Bahar, A. R. Atilgan, and B. Erman, Folding & Design, 2, 173 (1997)


2) Command Description

   *Generation of body-based modes by using the harmonic potential 
   instead of the regular CHARMM potential.*

   By applying "HARMonic" command after "MODES", the simple harmonic 
   potential, instead of the regular CHARMM potential, will be used to 
   generate modes for the body. This "HARM" command works only when
   the body-based modes are generating. If you want to read the existed
   modes, this command will not do any effects.

   The syntax for the MODE sub-command for using the harmonic potential, is:

   ::
   
      MBOND  [NBODY integer]

         MODEs { [HARMomic] [FCON] [RCUT] [VDWR]} Optional-selection (see below)

               - Regular mode generation command, such as USE, SORT, WRITE, etc.
                 (Check interface.note)

              END

      END

      ***Optional-selection { [MULTiple] [BOND] [ANGLE] [ALL ] [VDWF]}

      ----------------------------------------------------------  

      MODES { [HARMomic] [FCON] [RCUT] [VDWR] } Optional-selection 

      options:

              HARMonic invokes the harmonic potential mode generation for the
                       flexible-body.
 
              FCON    defines a single force constant of the harmonic potential
                      to be used for mode-generation. (No need to be defined when
                      MULT command is on.)

              RCUT    defines the cut-off distance for pairing atoms. Only distance
                      criteria is used to define a pair of atoms. 

              VDWR    sets an initial cut-off distance to be a sum of van der Waals 
                      radius of interacting atoms.
                      Final Cut-off distance = RCUT + VDWR-atom1 + VDWR-atom2

              ### Optional-selection Commands for Harmonic potential modes ###
  
              BOND
              ANGLE
              ALL     These options overide the distance criteria for defineing 
                      pair of atoms. BOND command forces to make bonded atoms as
                      a interaction pair. ANGLE command forces to connect 
                      atoms which have a common atom center. ALL command will involke
                      both BOND and ANGLE command. The remaining of atoms are 
                      paired by the distance criteria.

              MULTIple forces to use the different force constant based on the
                       nature of interactions. Each harmonic force constant is
                       obtained from the CHARMM parameters, such as bond, angle,
                       Lennard-Jones (VDW) parameters.

              VDWF    In MULT selection, the force constant (K) between paired atoms,
                      that are connected by the distance criteria, is defined by the 
                      second derivative of the Lennard-Jones potential.
                
                           K = d^2 V /dr^2  at where r is the given distance.  (1) 

                      IF VDWF is on, 
                           K = d^2 V/ d r^2 at where r is the potential minima  (2)


3) Examples

   ::
   
      ----- A simple single parameter approach --------

       (1)
           MODEs HARM VDWR RCUT 2.0  FCON 100.0
             Open write ...
             generate  ....
             Diag ...
             write ...
           END
   
             -  A single parameter Harmonic potential approach.
             -  Connecting pairs are defined by the distance cutoff.
                When The distance between two atoms are  less than 
                  2.0 + vdw radius of atom 1 + vdw radius of atom2,
                those atoms are connected by a single spring.

       (2) 
 
          MODEs HARM BOND RCUT 5.0  FCON 100.0
             .....
          END
    
             - A single parameter harmonic potential approach
             - First, atoms which are  covalent bonded to each other 
               are connected by a spring with force constant = 100kcal/mole*A^2
               Then, rest of atoms are examined by the distance criteria.
               Here, two atoms, not covalent bonded, are connected by a spring
               if the distance between those atoms are less than 5A.

       (3) 

          MODEs HARM ALL RCUT 2.0  FCON 50.0
            ....
          END

            - A single parameter harmonic potential approach
            - First, atoms which are covalent bonded to each other or are connected
              by a common atoms (ANGLE) are connect by a spring with 
              force constant = 50Kcal/mole*A^2. Then, rest of atoms are checked
              by the distance. If the distance between two atoms are less than
              2A, those atoms, which are not already paired by the first check
              (BOND and ANGLE)  are connected by a spring.


      -------- Multipe harmonic potential approach  -------

       (4)
 
          MODEs HARM MULTI VDWF RCUT 5.0  
             .....
          END

            - Multiple parameter harmonic potential approach
            - Atom pairs are solely created by checking the distance between
              atoms. Then, each pair's force constant is calculated by the
              second derivative of Lennard-Jones potential. "VDWF" command
              involk to use Eq. (2), defined above, for the calculations of 
              the force constants.

       (5) 
  
          MODEs HARM MULTI ALL VDWF RCUT 5.0
              ....
          END

            - Multiple parameter harmonic potential approach
            - First, atoms which are covalent bonded to each other or 
              are connected by a common atoms (ANGLE) are connect by a spring
              with force constants which are calculated from the CHARMM
              parameters (BOND and angle parameters) Then, the rest of atoms
              are checked by the distance criteria (<5.0A). IF the distance
              is less than 5A, the atoms are connected by a spring with a
              force constant calcualted by Eq. (2).

       (6) 
          MODEs HARAM
            .....
          END

            - This simple command will be treated as the same command as
              (5).

If you have any question, please contact Masa (watanabe@moldyn.com)
