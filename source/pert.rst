.. py:module:: pert

#####################################
Free Energy Perturbation Calculations
#####################################

The :chm:`PERTurbe` command allows the scaling of energy between PSFs for use in
energy analysis, comparisons, slow growth free energy simulations,
widowing free energy simulation, and for slow growth homology modeling.
This is a rather flexible implementation of free energy perturbation
that allows connectivity to change.  Also, three energy restraint
terms (harmonic, dihedral and NOE) and the SKIP command flags are subject
to change which allows a flexible way in which to compute free energy
differences between different conformations.  This code in implemented
in a non-intrusive manner and works with all minimizers and
integrators.  :chm:`SHAKE` can now be applied to bonds which are mutated as
well; an appropriate constraint corrections is calculated
automatically in these cases.


.. _pert_syntax:

Syntax of Free Energy Perturbation Commands
-------------------------------------------

::

   PERTurb  [OFF] [INBFrq int nonbond-specs] [RESEt] [MMFP] soft-core-spec
            atom-selection [INBFrq  0] [CHEM]               ]

   soft-core-spec::= [SCL0 SCR0 real] [SCL1 SCR1 real]

The :chm:`PERT OFF` command disables the free energy routines and the current
(lambda=1) PSF is used for subsequent commands.

When :chm:`OFF` is not specified, this command saves the current PSF as the lambda=0
state.  The atom-selection indicated which atoms have changed.  This is to
make the calculation run more efficiently.  If only a small percent of the
atoms have changed, this doubles the performance.  The nonbond specs are
included to make sure the nonbond exclusion lists are properly setup.  This
then allows the connectivity to change during the simulation.  :chm:`INBFrq` should
not be set to zero here unless the exclusion lists have already been setup
in a previous command.

Adding the :chm:`MMFP` keyword signals that the user intends to
change :doc:`MMFP <mmfp>` restraints as part of an alchemical
transformation. If it is omitted, the MMF-potentials are used as constant 
energy terms as was the case in CHARMM versions up to c30a2(x).

The ``soft-core-spec`` specifies whether the shift-based soft-core potential
applies to the the simulations.  The soft-core only applies on the
repulsive part (as defined by the WCA separation) of the Lennard-Jones
potential.  To get meaningful free energy, the soft-core potential
only is to be used when the electrostatic and WCA attractive
interaction between the selected atoms and the rest of the system are
turned off by using the :chm:`scalar CHARGE` and :chm:`scalar WCAD`.  There are two sets
of parameters that correspond to the energy evaluation
in the initial and final state.  :chm:`SCL0` specifies that soft-core
potential is used for the initial state.  :chm:`SCR0` is the
corresponding parameter that controls the strength of the soft-core
potential for the initial state.  :chm:`SCL1` and :chm:`SCR0` are the corresponding
values for the final state.  When :chm:`SCR*` equals to zero, all WCA repulsive
interaction of the selected atoms is off.  When :chm:`SCR*` equals to one,
full WCA repulsive interaction of the selected atoms is on.
A typical usage of the soft-core potential is to compute the free energy
from non-interaction state to repulsively interaction state for
a given selection of atoms.  The simulation can be conducted with
a series of Widom insertion/deletion stages.  The following values
of the :chm:`SCR*` are typical to perform the simulation in stages.

====== ==== ====
Window SCR0 SCR1
====== ==== ====
1      0    0.2
2      0.2  0.3
3      0.3  0.4
4      0.4  0.5
5      0.5  0.6
6      0.6  0.7
7      0.7  0.8
8      0.8  0.9
9      0.9  1.0
====== ==== ====

Because of the number of stages, the repulsive potential change in
each stage is so small that a single Widom insertion or deletion
(the sampling is done at the two end points, no intermediate windows
required) is sufficient for the free energy in each stages.
Notice this implementation of soft-core potential only works
with FEP or WHAM (exponential average).  Once the selected atoms
interact with the rest of the system with full repulsion, it is
straightforward to compute free energy contributions from
WCA attraction and electrostatics.  For details,
see Deng Y. and Roux B., J. Phys. Chem. B, 108 (42) 16567-16576

.. note::
   The WCA separation is implemented for:
   
   * group based slow (fast off) vdW routine,
   * atom based slow routine  with distance switch cutoff (fast off),
   * group based fast (fast gene) vdW routine and
   * atom based fast (fast gene) vdW energy with distance switch cutoff.

:: 

   ENERgy   ...   [ RESET ]  [ free-energy-step-spec ]
   DYNAmics ...              [ PUNIt integer         ]   [WHAM integer] 
   MINImize ...
   
     [ RESET ]     ! Resets all all accumulation data and counters.
                    (automatic for the first PERT or after a PERT OFF command)
   
   free-energy-step-spec::=
    [PWINdow [LAMBda real] ] [PSTArt int] [PSTOp int] [LSTArt real] [LSTOp real] -
    [PSLOwgrowth           ]   
   
            [PINCrement int]  [PEQUilibrate int]  [LAVErage]  [LINCrement real]
   
   
     [PWINdow   ]  ! specifies the windowing algorithm (default)
     [PSLOwgrowth] ! specifies the slow growth algorithm
     
     [LAMBda real] ! specifies the lambda value for windowing methods or for
                     energy or minimization calculations.
   
     [PSTArt int]  ! starting dynamics step number for accumulation (default 1)
     [PSTOp  int]  ! ending dynamics step number for accumulation (default 0)
   
     [LSTArt real] ! specifies the starting lambda value (default 0.0)
     [LSTOp  real] ! specifies the final lambda value (default 1.0)
   
     [PINCrement int] ! specifies number of steps to next window (auto mode).
     [PEQUil int]  ! specifies number of steps for equilibration (auto mode).
     [LAVErage]    ! specifies that lambda = (LSTART+LSTOP)/2    (auto mode).
     [LINCrement]  ! Specifies the lambda increment between windows (auto mode).
   
     [PSSP]        ! use soft core potentials for interactions in reac.
                   ! and product list.  This option is remembered. With
                   ! the PSSP keyword, two parameters, ALAM and DLAM can
                   ! be set.
   
     [ALAM real]   ! Separation parameter for elec. interaction (defaults to 5A^2)
   
     [DLAM real]   ! Separation parameter for LJ interaction (defaults to 5A^2)
   
     [NOPSsp]      ! Turn off use of soft core interactions.
   
     [ END   ]     ! Turns off the free energy code
  

The :chm:`PSTArt` and :chm:`PSTOp` values are relative to the number of dynamics steps
since :chm:`PERT` command was first enabled, or if a :chm:`PERT RESET` command is used
(regardless of whether the DYNAmics commands is run more than once or
whether the dynamic run involved the use of restart files).

By specifying the auto mode parameters (:chm:`PINCrement`, :chm:`PEQUilibrate`, :chm:`LINCrement`),
a new window will start at the conclusion of the current window with modified
parameters.  Also, in auto mode, the run will terminate at the end of a window
where the :chm:`LSTOP` value is 1.0 or 0.0.

The new commands :chm:`PSSP`/:chm:`NOPSsp` and the optional parameters :chm:`ALAM` and :chm:`DLAM`
control the interactions between soft core potentials and PERT.  After
you specify PSSP in an energy, minimization or Dynamics command, soft
core LJ and elec. interactions are used in all reactant and product
nonbonded list terms. The separation parameters for elec. and LJ
interactions can be set with the ALAM and DLAM options, the default of
5A^2 should be reasonable.  The option is remembered, i.e., after the
first invocation of PSSP all further calls to EPERT use soft core
interactions. To turn this off, use the NOPSsp keyword.  Use of
softcore interactions is also turned off after a PERT RESEt or a PERT
OFF command.  For example (assuming that PERT has been already set
up):

::

   ENER LAMB 0.5 ! calc. system energy for L=0.5
   ENER LAMB 0.5 PSSP ! calc. system energy for L=0.5 using soft core potentials
   ENER LAMB 0.6 ! calc. system energy for L=0.6 using soft core potentials
                 ! since the PSSP option is remembered
   ENER LAMB 0.6 NOPS ! calc. system energy for L=0.6, use of soft 
                 ! core pots. turned off.

The PERT/PSSP code supports only thermodynamic integration and slow-growth.
Please ignore all results starting with TP> when using PSSP.

The PUNIt option allows the free-energy-step-spec to be specified 
more than once and acts as a scheduler for a particular simulation.
The format of the PUNIt file is;

::

   * title
   *
   repeat-lines-of(free-energy-step-specs)
   END

The end is optional and terminates the free energy run.
Lack of an END (i.e. and end-of-file or blank lines) will put PERT into
auto mode which will continue until LSTOP becomes 1.0 or 0.0 (based on
the LINCrement value).

The WHAM option allows to write a formatted file to post-process with 
the Weighted Histogram Analysis Method (see below).


.. _pert_description:

Description of PERT Commands
----------------------------

The PERTurb command copies and saves the current PSF and restraint
data for harmonic, NOE and dihedral restraints to the initial (lambda=0)
saved state.  The SKIP command flags are also saved to allow linear scaling
of entire energy terms.

The structure may then be modified or perturbed with patches,
SCALar commands, with the DELEte command, or by generating or reading
a new PSF.

The Basic mode of operation is;

::
   
   ....
   PERTurb        ! Define the lambda=0 state.
   PATCH ....     ! Define the lambda=1 state.
   DYNAmics ....  ! Run MD on intermediate surface...
   ....

The PSF in use when dynamics or energy minimization is invoked
becomes the final (lambda=1) state.  The actual energy computed
is a linear combination of these two endpoints.

The PATCH command may be replaced with any other command that
modifies the PSF.  Some examples which modify the PSF;

::

   SCALAR CHARGE SET -0.55 SELE ATOM A 1 O* SHOW END ! change a charge

   DELETE ATOM SELE ALL END
   READ SEQUENCE ....
   GENERATE ...    ! generate a new different PSF.

   DELETE CONNECTIVITY ....  ! modify the PSF by changing the connectivity.
                   ! see a word of warning below!

   SCALAR TYPE SET 14 SELE ATOM A 1 O SHOW END ! change the vdw atom type 

   OPEN ....     ! Read a new PSF
   READ PSF ... 

It is not required that the PSF be modified.  If one wants to
carry out coordinate perturbation only, it is sufficient to modify
the harmonic restraints, the NOE restraints, or the dihedral restraints.
In this way, it is possible to calculate the free energy differences
between different conformers.  (However, this option should not be
used with simultaneous change of SHAKE constraints)

Note that with this implementation, because two PSFs are used,
that the connectivity may change.  The use of 1-4 interactions and
nonbond exclusions is fully supported.  This allows this method to be
used for examining changes that involve bond changes, such as cystine
bridge formation. [Note added by Stefan Boresch (stefan@mdy.univie.ac.at)
Changes in connectivity that involve bond breaking or forming are
highly problematic and may not converge.  This is explained in
detail in Boresch & Karplus, J. Phys. Chem. B 1999, 103, 103-136.
The flexibility made possible by the implementation of PERT puts
the responsibility of what can be done and what not on the user!]

The advanced mode of operation is;

::

   ....
   PERTurb
   PATCH ....
   DYNAmics ....
   ....
   PERTurb
   PATCH ....
   DYNAmics ....
   ....
   PERTurb
   PATCH ....
   DYNAmics ....
   ....

In this way, several changes can be affected in a single CHARMM run.
For example, the first patch might be the removal of charges, and the
second patch could correspond to a change in atom size, and the third
patch could simply consist of modifying dihedral restraints so as to
affect a conformational change.  The free energy differences and
fluctuations will be calculated for each window as well as the total
for all previous windows.

.. _pert_restrictions:

RESTRICTIONS:
-------------

The number of atoms in both sets must match!  If the system of
interest has different numbers of atoms, then dummy atoms must be
used.  The mapping of atoms between the first and last structure is
one to one. 

The following CHARMM features are not currently supported for use
with free energy perturbation;

:: 

   INTEraction_energy

These commands will continue to work, but will only use the final
(lambda=1) structure.  Most other energy related CHARMM features
are supported.  

The IMAGE/CRYSTAL facility has been supported now for some time;
however, IMAGE/CRYSTAL needs to be set up *after* the PERT command!!!

The following CHARMM energy related features cannot be modified with
the PERT command (e.g. cannot be part of what is changing, and are
only determined by the final state).

* HBON - hydrogen bond energy
* ST2  - ST2 energy
* CIC  - internal coordinate constraint energy
* CDRO - quartic droplet potential energy
* USER - user supplied energy (USERLINK)
* RXNF - Reaction field energy
* IMNB - image van der Waal energy
* IMEL - image electrostatic energy
* IMHB - image hydrogen bond energy
* IMST - image ST2 energy
* SBOU - solvent boundary energy
* UREY - Urey Bradley energy terms
* XTLV - Crystal vdw terms
* XTLE - Crystal electrostatics

Extended electrostatics is implemented within PERT and can be used with 
the following CHARMM commands:

::

   NBONDS GROUP SWITCH CDIE VDW VSWITCH EXTEND GRAD QUAD -
   CTONNB 12.0 CTOFNB 12.0 CUTNB 12.0 WMIN 1.2 WRNMXD 1.2 EPS 1.0

.. note::
   The ctonnb, ctofnb and cutnb values should be the same when 
   implementing extended electrostatics in PERT to prevent problems with 
   mixing of usage of switching functions and extended electrostatics


.. _pert_references:

References
----------

* M Mezei and D.L. Beveridge, in Annals of the NYAS, "Free Energy
  Simulations" 482 (1986)
* T. P. Straatsma Ph.D. Thesis "Free Energy Evaluation by
  Molecular Dynamics Simulations"
* Kollman, P. A.; et al.  J. Am. Chem. Soc. 1987, 109, 1607.
* Kollman, P. A.; et al.  J. Am. Chem. Soc. 1987, 109, 6283.
* Kollman, P. A.; et al.  J. Chem. Phys. 1989, 91, 7831.
* Bell, C. D.; Harvey, S. C.,   J. Phys. Chem. 1986, 90, 6595.
* van Gunsteren, W.F. et al. in: Computer Simulation of Biomolecular
  Systems: Theoretical and Experimental Applications, vol. 2, eds. van
  Gunsteren W.F. and Weiner P.K. (Escom, Leiden, 1994), p. 349


.. _pert_examples:

Examples
--------

The input file
^^^^^^^^^^^^^^

:: 

   * A SIMPLE TEST RUN FOR PERT
   *
   bomlev -1
   OPEN READ FILE UNIT 1 NAME ~/c22pt/toph19.mod
   READ      RTF   UNIT  1
   OPEN READ FILE UNIT 2 NAME ~/c22pt/param19.mod
   READ      PARAMETER  UNIT  2
   READ      SEQUENCE  CARD
   *  FIRST SEQUENCE FOR SECOND DERIVATIVE TEST
   *
       2
   AMN CBX
   GENERATE  A
   GENERATE  B DUPLICATE A
   
   OPEN UNIT 3 READ CARD NAME perttest.crd
   READ COOR CARD UNIT 3
   
   ! modify the charge for the lambda=0 state
   SCALAR CHARGE SET -0.55 SELE ATOM A 1 O* SHOW END
   
   ! minimize initial state so initial forces will be small.
   MINI ABNR NSTEP 100 CTOFNB 12.0 CUTNB 14.0
   
   PERT  ! save all PSF data for the lambda=0 state
   
   ! modify the charge again for the lambda=1 state
   SCALAR CHARGE SET -0.15 SELE ATOM A 1 O* SHOW END
   
   ! carry out free energy run from first to final state
   
   OPEN READ CARD UNIT 88 NAME perttest.punit
   DYNA VERLET STRT  NSTEP 12000 TIMESTEP 0.001 -
       IPRFRQ 100 IHTFRQ 0 IEQFRQ 100 NTRFRQ 2000  -
       IUNCRD 50  ISEED 314159  -
       NPRINT 100 NSAVC 0 NSAVV 0 INBFRQ 25 IHBFRQ 0 -
       CTOFNB 12.0 CUTNB 14.0 -
       FIRSTT 300.0 FINALT 300.0 TEMINC 100.0   -
       IASORS 0 IASVEL 1 ISCVEL 0 ICHECW 1 TWINDH 20.0 TWINDL -20.0 -
       PUNIT 88
   
   PERT OFF
   energy ! just a check at lamda=1
   STOP

The punit file
^^^^^^^^^^^^^^

::

   * PUNIT FILE FOR SIMPLE TEST CASE
   * use window method with 2000 steps of equilibration
   * and 8000 steps of analysis for each of 10 evenly spaces
   * windows.
   *
     LSTART 0.0  LAMBDA 0.0  LSTOP 0.05  PSTART  12000  PSTOP  20000  PWIND
     LSTART 0.05 LAMBDA 0.1  LSTOP 0.15  PSTART  22000  PSTOP  30000  PWIND
     LSTART 0.15 LAMBDA 0.2  LSTOP 0.25  PSTART  32000  PSTOP  40000  PWIND
     LSTART 0.25 LAMBDA 0.3  LSTOP 0.35  PSTART  42000  PSTOP  50000  PWIND
     LSTART 0.35 LAMBDA 0.4  LSTOP 0.45  PSTART  52000  PSTOP  60000  PWIND
     LSTART 0.45 LAMBDA 0.5  LSTOP 0.55  PSTART  62000  PSTOP  70000  PWIND
     LSTART 0.55 LAMBDA 0.6  LSTOP 0.65  PSTART  72000  PSTOP  80000  PWIND
     LSTART 0.65 LAMBDA 0.7  LSTOP 0.75  PSTART  82000  PSTOP  90000  PWIND
     LSTART 0.75 LAMBDA 0.8  LSTOP 0.85  PSTART  92000  PSTOP 100000  PWIND
     LSTART 0.85 LAMBDA 0.9  LSTOP 0.95  PSTART 102000  PSTOP 110000  PWIND
     LSTART 0.95 LAMBDA 1.0  LSTOP 1.0   PSTART 112000  PSTOP 120000  PWIND
     END
   
Or equivalently using auto mode:

::
   
   * PUNIT FILE FOR SIMPLE TEST CASE
   * use window method with 2000 steps of equilibration
   * and 8000 steps of analysis for each of 10 evenly spaces
   * windows.
   *
    LSTART 0.0  LAMBDA 0.0  LSTOP 0.05  PSTART  12000  PSTOP  20000  PWIND PEQUIL 2000 PINCR 10000 LINCR 0.1


Or also equivalently as:

:: 

   * PUNIT FILE FOR SIMPLE TEST CASE
   * use window method with 2000 steps of equilibration
   * and 8000 steps of analysis for each of 10 evenly spaces
   * windows.
   *
     LSTART 0.0  LAMBDA 0.0  LSTOP 0.05  PSTART  12000  PSTOP  20000  PWIND
     LSTART 0.05 LAMBDA 0.1  LSTOP 0.15  PEQUIL   2000  PINCR  10000  PWIND
     LSTART 0.15 LAMBDA 0.2  LSTOP 0.25  PEQUIL   2000  PINCR  10000  PWIND
     LSTART 0.25 LAMBDA 0.3  LSTOP 0.35  PEQUIL   2000  PINCR  10000  PWIND
     LSTART 0.35 LAMBDA 0.4  LSTOP 0.45  PEQUIL   2000  PINCR  10000  PWIND
     LSTART 0.45 LAMBDA 0.5  LSTOP 0.55  PEQUIL   2000  PINCR  10000  PWIND
     LSTART 0.55 LAMBDA 0.6  LSTOP 0.65  PEQUIL   2000  PINCR  10000  PWIND
     LSTART 0.65 LAMBDA 0.7  LSTOP 0.75  PEQUIL   2000  PINCR  10000  PWIND
     LSTART 0.75 LAMBDA 0.8  LSTOP 0.85  PEQUIL   2000  PINCR  10000  PWIND
     LSTART 0.85 LAMBDA 0.9  LSTOP 0.95  PEQUIL   2000  PINCR  10000  PWIND
     LSTART 0.95 LAMBDA 1.0  LSTOP 1.0   PEQUIL   2000  PINCR  10000  PWIND
     END


An annotated output example
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This output is a short excerpt from perttest.out where the input
lines start with "|" and output lines start with ".....|".

The CHARMM command is:

::

   | CHARMM>    DYNA VERLET STRT  NSTEP 12000 TIMESTEP 0.001 -
   | CHARMM>        IPRFRQ 100 IHTFRQ 0 IEQFRQ 100 NTRFRQ 2000  -
   | CHARMM>        IUNCRD 50  ISEED 314159  -
   | CHARMM>        NPRINT 100 NSAVC 0 NSAVV 0 INBFRQ 25 IHBFRQ 0 -
   | CHARMM>        CTOFNB 12.0 CUTNB 14.0 -
   | CHARMM>        FIRSTT 300.0 FINALT 300.0 TEMINC 100.0   -
   | CHARMM>        IASORS 0 IASVEL 1 ISCVEL 0 ICHECW 1 TWINDH 20.0 TWINDL -20.0 -
   | CHARMM>        PUNIT 88

and the relevent punit data is in perttest.punit:

:: 
   
   |* PUNIT FILE FOR SIMPLE TEST CASE
   |* use window method with 2000 steps of equilibration
   |* and 8000 steps of analysis for each of 10 evenly spaces
   |* windows.
   |*
   | LSTART 0.0  LAMBDA 0.0  LSTOP 0.05  PSTART  1200 -
   |   PSTOP  2000  PWIND PEQUIL 200 PINCR 1000 LINCR 0.1

The output starting at line 1618 (in the middle of the dynamics command) is:

::

   .....| PERTURBATION> Free energy perturbation results:
   
This indicates that a "window" was just completed.

::

   .....| PERTURBATION> results, LSTART=    0.050000  LSTOP=    0.150000  LLAST=    0.100000 Number of steps used=       800
   
This says that the window started at lambda=0.05 and ended at 0.15.
The window "center" was at lambda=0.1
A total of 800 steps was used for collecting averages and fluctuations.

::

   .....| PERTURBATION> result: EXPAVE=0.155456E+01 EXPFLC=0.225213E+00 DIFAVE=   -0.256530 DIFFLC=    0.088960

The values:

* EXPAVE is the time average of :math:`e^{(E_f(t)-E_i(t)-E_f(0)+E_i(0))/kT}`
* EXPFLC is the fluctuation of this value about its average
* DIFAVE is the time average of :math:`E_f(t)-E_i(t)-E_f(0)+E_i(0)`
* DIFFLC is the fluctuation of this value about its average

.. note::
   this value should not be much larger than kT for a
   good window schedule.  If this value is too large,
   then smaller window lambda steps should be used.
   The value here (0.09) indicates that a much larger
   window would have been OK.  ef(t) is the energy at
   lambda=LSTOP, ei(t) is the energy at lambda=LSTART
   ef(0) is the initial energy at LSTOP, ei(0) is the
   initial energy at LSTART

::

   .....| PERTURBATION> TP Windowing result, EPRTOT=    1.392400  EFORWARD=    0.914282 EPREF=    1.177303

This is the old format. In the new format the values EPRTOT,EFORWARD,
EBACKWARD,EPREFF,EPREFB are reported where the forward energy is from LLAST
to LSTOP and the backward from LLAST to LSTART.  Separating the window into
two halves (double wide sampling) improves the accuracy of the the TP method.

::

   .....| PERTURBATION> TI Windowing result, EPRTOT=    1.400218  EFORWARD=    0.920773 EPREF=    1.177303

* EPRTOT is the total energy for this window and all previous (since a PERT RESET).
* EFORWARD is the energy for this current window.  
* EPREF is the initial energy difference (ef(0)-ei(0))

::

   .....| PERTRES> LSTART=     0.05000 LSTOP=     0.15000 EPRTOT=     1.40022  EFORWARD=     0.92077 EPREF=     1.17730 DIFAVE=    -0.25653 DIFFLC=     0.08896

This is the same data on a one line format.  To use this, grep (search)
for "PERTRES".

::

   .....| PERTURBATION> Averages for the last      800  steps:
   .....|PAVE DYN: Step         Time      TOTEner        TOTKe       ENERgy  TEMPerature
   .....|PAVE PROP:             GRMS      VEREner        VERKe       EHFCor        VIRKe
   .....|PAVE EXTERN:        VDWaals         ELEC       HBONds                      USER
   .....|PAVE PRESS:            VIRE      VIRI           PRESSE      PRESSI    VOLUme
   .....| ----------       ---------    ---------    ---------    ---------    ---------
   .....|PAVE>      800      0.00000      0.92077     -0.00039      0.92077      0.00000
                                                      <delF*v>       <delE>
   .....|PAVE PROP>          0.00000      0.00000      0.00000      0.00000      0.00000
   .....|PAVE EXTERN>        0.00000      0.92077      0.00000                   0.00000
   .....|PAVE PRESS>         0.00000      0.00000      0.00000      0.00000      0.00000

This is the average values for this window (TI results).

.. note::

   the <delF*v> is a "correction" term that shows how well the window is
   equilibrated.  This value should be close to zero and much smaller than
   its fluctuation.  If it is not, then the assumptions required for a
   free energy calculation using windowing are not met.
   This can occur if there is a "snapping" event with releases energy in
   an irreversible manner, or if the system is not at equilibrium.
   For a slow growth "window", this is a correction term which should be
   multiplied by the estimated delay of the equilibration at the current
   step and then added to the total.
   For example, if the configuration distribution lags behind the energy
   potential by 100 steps, this value should be scaled by 100.
   When forcing a "continuous" change through slow growth, there tends to
   be a delay since the structure does not respond instantly to the
   potential.

::

   .....| ----------       ---------    ---------    ---------    ---------    ---------
   .....| PERTURBATION> Fluctuations for the last      800  steps:
   .....|PFLC>      800      0.00000      0.08911      0.00985      0.08896      0.00000
   .....|PFLC PROP>          0.00000      0.00000      0.00000      0.00000      0.00000
   .....|PFLC EXTERN>        0.00000      0.08896      0.00000                   0.00000
   .....|PFLC PRESS>         0.00000      0.00000      0.00000      0.00000      0.00000
   
This is the fluctuation data for the current window (TI results).

::
   
   .....| ----------       ---------    ---------    ---------    ---------    ---------
   .....| PERTURBATION> TOTALS since last reset:
   .....|PTOT>      800      0.00000      1.40022     -0.00038      1.40022      0.00000
   .....|PTOT PROP>          0.00000      0.00000      0.00000      0.00000      0.00000
   .....|PTOT EXTERN>        0.00000      1.40022      0.00000                   0.00000
   .....|PTOT PRESS>         0.00000      0.00000      0.00000      0.00000      0.00000
   
This is the total for this window and all previous (since the last PERT RESET).

::

   .....| ----------       ---------    ---------    ---------    ---------    ---------
   .....|
   .....| PERTURBATION> EOF on punit file: PERT in auto mode.
   
This indicates that no more data was found on the PUNIT file.

:: 

   .....| PERTURBATION> Free energy perturbation calculation continues.
   
Now we start a new window.

::

   .....| PERTURBATION> PSTART=        3200  PSTOP=        4000
   
We will run equilibration until step 3200, then we will collect data until step 4000.

::

   .....| PERTURBATION> LSTART=    0.150000  LSTOP=    0.250000  LAMBDA=    0.200000
   
This indicates the boundaries of the new window.

::

   .....| PERTURBATION> Windowing will be used.
   
We are using windowing (fixed lambda values), instead of slow growth.


.. _pert_constraints:

Constraints
-----------

If SHAKE is applied to bond terms which are changed as the
result of an alchemical mutation a constraint correction is calculated
where required in slow-growth mode and TI in windowing mode.  The
exponential formula in windowing mode is not supported.  The user has
to beware of subtle problems regarding a possible "moment of inertia"
term that may be or may be not included in this correction (S. Boresch
& M. Karplus, to be published) In order for the constraint correction
to work properly, attention has to be given to the following points:

1. SHAKE must not be applied to angle terms
2. the PARA option has to be used (it's not clear, how to support
   reference coordinates in the context of an alchemical mutation)
3. the SHAKE command has to issued after the PERT command.  (This
   is similar to setting up IMAGEs in connection with PERT).  A
   typical input will look something like

   :: 

 	    PERT SELE ... END
	    !change psf; after ALL changes have been made
   	  SHAKe BOND PARA
   	  DYNA ... ! carry out MD simulations etc.
	    PERT OFF

4. One should not mix situations where a constraint correction for
   SHAKE is necessary with the use of harmonic, NOE and dihedral
   restraints to calculate conformational free energy differences.

Items (2) and (3) can lead to error messages in situations where there
is actually no problem, e.g. you just want to apply SHAKe to your
solvent which is not affected by the mutation, so you specify SHAKe
before PERT and "bomb".  Nevertheless, I thought better safe than
sorry and if wanted one can override the warnings with a BOMLEV -3.
Item (4) is simply untested.


.. _pert_wham:

WHAM
----

The post-processing Weighted Histogram Analysis Method (WHAM) can 
be used to help reaching better statistical convergence on free energy 
perturbations calculations.  The approach represents a generalization 
of the histogram method developed by Ferrenberg and Swendsen (1989).  
The central idea, which goes back to the maximum overlap method 
developed by Bennet (1976) to estimate free energy differences,
consists in constructing an optimal estimate of the unbiased distribution 
function as a weighted sum over the data extracted from all the 
simulations and determining the functional form of the weight factors that 
minimizes the statistical error.  The WHAM approach can be used to calculate 
the PMF along coordinates (Kollman, 1992; Boczko, 1993; Brooks, 1993; 
Roux, 1995) and can also be used to post-process free energy perturbation
calculations in which no PMF is desired. 

Assuming that the simulations are performed (as in PERT) with a potential
function given by a linear switching from E_0  to E_1, that is, 

.. math::

   E_k = (\lambda_{k-1})*E_1 + (\lambda_k)*E_0

The free energy constants F_k corresponding to any :math:`\lambda_k` window are 
given by,

.. math::

   e^{-F_k/kBT} = \sum_i \left( \sum_t \frac{ \mathrm{Top}_i(t) }{ \mathrm{Bot}_i(t) } \right)

where 

.. math::

   \mathrm{Top}_i(t) = e^{-E_k[X_i(t)]/kBT}
   
and

.. math::

   \mathrm{Bot}_i(t) = \sum_j  { Ntime(j) * e^{+F_j/kBT-E_j[x_i(t)]/kBT} }

where :math:`E_j[x_i(t)] = (\lambda_j-1)*E_1[x_i(t)] + (\lambda_j)*E_0[x_i(t)]`
is the potential energy function of the j-th window evaluated with 
the configuration taken from the i-th simulation.  The WHAM equations 
for the F_k must be solved iteratively.  

The syntax of the command is:

::

   WHAM MAXWindow <integer>  MAXTime <integer> unit <integer> -
        tolerance <real>  nstep <integer> [guess] -
        ioffset <integer>   nskip <integer> {lambda <real> lambda <real> ...}

where

* MAXWindow is the total number of windows
* MAXTime   is the total number of time-step configuration for each window
* unit      is the unit of the formatted file containing all the information
* tolerance is the tolerance on the F_k to reach convergence
* nstep     is the maximum number of iterations on the WHAM equations
* guess     to flag that an initial guess is provided for the F_k. Those
  are read directly from the input stream with one line per
  window in the format [ window <integer>   F() <real> ] 
* nskip     use only every nskip data point to reach convergence (faster)
* lambda    give any value of lambda for which you want the free energy (a list)
* ioffset   reference energy level at window number "ioffset".

The file containing the information can be written by PERT during dynamics
if the keyword WHAM is used in the dynamics command (see above).  In 
principle, the WHAM could also be used with non-linear perturbations, 
but then the code in PERT would have to evaluate several energies since 
those could not be obtained by lambda interpolations.  

The total free energy of is stored in 'WHAMFE' substitution viriable.

Some references on WHAM:

* S. Kumar, D. Bouzida, R.H. Swendsen, P.A. Kollman, and J.M. Rosenberg.
* J. Comp. Chem. 13, 1011--1021 (1992).
* E.M. Boczko and C.L. Brooks III.  J. Phys. Chem. 97, 4509--4513 (1993).
* A.M. Ferrenberg and R.H. Swendsen.  Phys. Rev. Lett. 63, 1195--1198 (1989).
* C.M. Bennet.  J. Comp. Phys. 22, 245--268 (1976).
* C.L. Brooks III and L. Nilsson.  J. Am. Chem. Soc. 115, 11034--11035 (1993).
* B. Roux, Comp. Phys. Comm. 91, 275-282 (1995).


.. _pert_pssp

PERT/PSSP
---------

Some details concerning the implementation of the PERT/PSSP code,
including present limitations:

Introduction
^^^^^^^^^^^^

The PERT free energy module of CHARMM is based on a linear dependence
on the coupling parameter.  While simplifying implementation, this
approach is prone to van der Waals endpoint problems.  One widely used
method to overcome the van der Waals endpoint problem is the use of
soft core Lennard Jones and electrostatic interactions for those
energy terms that cause problems.  This capability has been added,
following Zacharias, Straatsma and McCammon, J. Chem. Phys. 1994, 100,
9025.

Outline of the implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the following L denotes the coupling parameter lambda. Subscripts
_i and _f denote initial and final state respectively. The variables
de (ALAM) and dv (DLAM) can be set by the user; reasonable defaults (5
A^2) are used.

The functional form of the soft core routines in combination with PERT
is as follows:

::

   U_LJ(L) = U_LJ,0 + 

                        A_f                      B_f
       L     * (-------------------  -  --------------------) +
                 (r^2 + dv*(1-L))^6      (r^2 + dv*(1-L))^3

                        A_i                      B_i
       (1-L) * (-------------------  -  --------------------) +
                 (r^2 + dv*L)^6             (r^2 + dv*L)^3


   U_ELEC(L) = U_ELEC,0 + 

                      qi_f*qj_f       
       L     * (----------------------) +
                 sqrt(r^2 + de*(1-L)) 


                      qi_i*qj_i       
       (1-L) * (----------------------)
                   sqrt(r^2 + de*L) 


Of course, the effects of tapering functions (SHIF, SWIT etc.)  have
to be taken into account properly; this is particularly important for
the calculation of the forces and of the derivative dU/dL

In principle, soft core potentials are only required for atoms that
'vanish' at one of the endpoints (i.e., dummy atoms). In this
implementation, a simpler approach was used: When PERT is activated,
CHARMM uses three nonbonded lists (six with IMAGE/PERT): (1) one for
the "environment" part (the part of the system that remains
unchanged), (2) one for the "reactant" part (interactions between
initial state atoms themselves and initial state atoms and the
environment), and (3) one for the "product" part (interactions between
final state atoms themselves and final state atoms and the
environment). All reactant and product list interactions are calculated
using soft core potentials.  Since (see equations above) at the
endpoints the soft core expressions reduce to normal interactions, use
of the soft core potentials is equivalent to a modified
path, but the overall result of the free energy simulation is unchanged.
(Note: the effect on free energy components has not been explored
systematically!)

Obviously, using soft core potentials breaks the standard scheme how
PERT calculates dU(L)/dL since instead of

	   U(L) = U_0 + (1-L)*U_i + L*U_f  ! standard PERT

we now have

	   U(L) = U_0 + (1-L)*U_i(L) + L*U_f(L)  ! PERT/PSSP

with corresponding differences for dU/dL.  One sees that the
standard PERT scheme gives approximately "half" of the required
derivative, but we still need the terms (1-L)*[dU_i(L)/dL] and 
L*[dU_f(L)/dL].  The modified energy routines I use do these
additional calculations.

Summing up, modifications of the code were necessary in the
following three areas:

(a) Additions to the parser: PSSP/NOPSsp keywords, ALAM, DLAM
parameters; initializations and resets

(b) Modifications to subroutine EPERT itself, making (i) sure that the
correct energy subroutines are called if PSSP is active, and (ii) that
the additional contributions to dU/dL are temporarily stored and
added to the LJ and elec free energy contributions calculated in the
standard PERT way. To achieve (i), subroutines FASTST,
EVDW (enbonda.src) and EGROUP (enbondg.src) were modified slightly
as well.

(c) Special purpose nonbonded energy routines (based on the
standard slow energy routines) have been added to the file epert.src.
Currently only a subset of nonbonded options is supported (see below)

All new nonlocal variables (no arrays are needed!) have been added
in pert.fcm (epert.fcm is unchanged)

The outline given here together with the comments in the code
should make the inner workings of the PERT/PSSP code clear.

To quickly 'grep' for all changes, seach for lines(comments) starting
with Cpssp

For a description of how to use the new functionality (activated by
the PSSP keyword), see the modified PERT documentation and the new
test cases.


Comments and present limitations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To the best of my knowledge, all reported uses of soft core potentials
in free energy simulations have been based on thermodynamic
integration (TI), not the "exponential formula" ("thermodynamic
perturbation", "free energy perturbation (FEP)").  The present
implementation also supports only TI.  While the output claims to give
values obtained with the exponential formula (TP> lines in output
file), the reported values ARE W R O N G if soft core potentials are
used. This is similar to all cases when the constraint correction is
needed, which also only works with TI !!! At present, it is not clear
whether the exponential formula can be supported easily.

Only a limited subset of nonbonded options is supported at presented.
Nonsupported options are hopefully caught and should make CHARMM die.
Nevertheless, check against the following list:

At present, the following limitations apply: 

* Only constant dielectric (CDIE)
* For group based cutoffs (GROU / VGROU), the following nonbonded options
  are supported at present:

  ::
  
                VDW (LJ)                 ELEC
                -----------------------------
                  VSWI                   SWIT

* For atom based cutoffs (ATOM / VATO), the following nonbonded options
  are supported at present:
  
  ::

                VDW (LJ)                 ELEC
                -----------------------------
                  VSWI                   SHIF
                                         FSWI
                                         EWAL (tradional or PME)

Finally, note that there is no support for parallel architectures.

Outlook / TODO
^^^^^^^^^^^^^^

* Support parallel architectures (someone else needs to do this, since
  I have no hardware) -- this should probably postponed until a merge
  with the standard slow energy routines

* All PERT/PSSP specific modifications could easily be put behind a
  separate compilation flag (e.g., ##IF PSSP) if this were 
  desired.

* Support additional tapering functions: While it doesn't seem
  necessary to support all combinations of nonbonded options in
  CHARMM, support for FSHIft and RDIE is planned.

* Use the slow energy routines instead of special purpose routines in
  epert.src.  Further, maybe a merge with the existing soft core
  routines in CHARMM (intended primarily for docking) is possible.
  This would lead to a unified, general and flexible soft core
  facility. Also, this would cure (maybe?) (most of) the parallel code
  incompatibilities...

* Consider support of the exponential formula; if this cannot be
  done easily, remove the output lines to avoid confusion.


.. _pert_patch

PATCHING DUMMY SIDECHAINS FOR FREE ENERGY PERTURBATIONS
-------------------------------------------------------

The command MKPRES has been introduced to write a PATCH for adding a dummy
sidechain onto a backbone with the goal of performing free energy calculations.
The command generates the list of needed dihedral angles and non-bonded
exclusions.   Only the cross internal energy terms between the dummy sidechain
and the bakbone are introduced, the rest is generated from the normal generate
command. Basically, it is ok to use such a mixed topology (single vs dual)
by branching at the carbon CB.   With this treatment, all the bonds and angles
are kept, some dihedrals may need to be turned off for statistical consistency
of the reference state (this requires some thinking by the user, sorry...).
It can be shown that free energy differences calculated with these end-points
are correct (even though the individual free energies are themselves different
than those with ideal gas of free particle in which all internal energy terms
are turned off).  Proline is not supported by this method.  Glycine might 
be ok, but remain vigilant.

One can introduce dummy atoms, which retain all the covalent interactions,
in a "hybrid residue", in such a way that the influence of the bonded
interactions with dummy atoms do not influence the final free energy change.
The simulations thus can be done using a transformation protocol in which 
all covalent bond contributions are maintained invariant throughout the
calculations; only the nonbonded interactions are varied.  

The hybrid method is a scheme which retains some features of both the single
and dual topology techniques.  The dummy atoms, which are covalently linked
to the protein in question, have no nonbonded interactions at one  or the
other of the two end point reference states. The potential energy function
describing the transformation is constructed such that all internal energy
terms are invariant with respect to the thermodynamic coupling parameter
lambda.  This simulation procedure therefore has similarities with both
the "single topology" and "dual topology" methods.

The coupling of the dummy atoms to the real atoms cancels out exactly from
the calculated free energy differences.  The equivalence holds as long as
the coupling between the dummy atoms and the rest of the system satisfies
certain conditions.  First, there can be only one bond between the dummy atoms
of a mutated residue and the real atoms in the  rest of the system, because
multiple bonds would add spurious coupling between the real atoms. Second,
to avoid spurious coupling between the dummy atoms and the rest of the system,
there cannot be multiple bond angles and dihedral torsion angles between 
the dummy atoms of a transformed residue and more than two real atoms in
the rest of the system.

The theoretical arguments explaining the approach have been elaborated in
the following references:

* Boresch, S.; Karplus, M. J. Phys. Chem. A 1999, 103, 103-118.
* Boresch, S.; Karplus, M. J. Phys. Chem. A 1999, 103, 119-136.
* Shobana S.,B. Roux, and Olaf S. Andersen, J. Phys. Chem. B 2000, 104, 5179-5190

Syntax:

::
  
   MKPRES {PatchName} unit <integer> atom_selection atom_selection atom_selection atom_selection 

The four atom_selections correspond to the following pieces of the psf:

* first:  is the invariant backbone of state 0
* second: is the mutated sidechain for state 0
* third:  is the invariant backbone of state 1
  (normally this should correspond identically to the first selection)
* fourth: is the mutated sidechain for state 1

Here is an example for liking a dummy valine to an alanine:

:: 
   
   set Residue1 = ALA
   set Residue2 = VAL
   
   !first state
   read sequence card
   *  residue1
   *
     3
   ALA @Residue1 ALA 
   generate SEG1 setup 
   
   ! second state (must have complete second segment to generate internal
   ! energy terms of second state)
   read sequence card
   *  residue2
   *
     3
   ALA @Residue2 ALA
   generate SEG2 setup 
   
   define BACK select type CA .or. type HA* .or. type N .or. type HN .or. -
          type C .or. type O .or. type HT* .or. type OT* .or. type CB show end
   
   open write card unit 10 name mkpres.rtf
   write title unit 10
   ** Patch for alchemical mutation of ALA to VAL
   **
   *
   
   MKPRES @PatchName unit 10 -
          select segid SEG1 .and. back end -
          select segid SEG1 .and. resid 2 .and. (.not. back ) end -
          select segid SEG2 .and. back end -
          select segid SEG2 .and. resid 2 .and. (.not. back ) end 
   
Two patches are written in unit 10.  The first one is for real alanin/dummy
valine while the second patch turns the system into dummy alanine/real valine.
Please check your patch before lengthy calculations!


.. _pert_mmfp:

Details about the use of MFF-potentials in PERT
-----------------------------------------------

The MMFP keyword makes it possible to modify (add/remove)
MMF-potentials during an alchemical PERT mutation.  If the keyword was
specified when activating PERT, the original call to the MMFP routines
in the constant-energy section of EPERT (epert.src) is omitted.
Instead, two new calls to the MMFP routines, one in the lambda=0
section and one in the lambda=1 section, are used.  It is an optional
command, so the original use of the MMFP potentials in PERT (as
constant energy terms, independent from lambda) is still available
(and the default behaviour).  The possibility of lambda-dependent
MMF-potentials should facilitate free energy simulations where
auxiliary restraints are needed, such as calculations of absolute
binding free energies.

.. note::
   It is *important* to choose a MAXGEO value
   (see *note geo:(chmdoc/mmfp.doc)) appropriate for
   the sum of all GEO restraints that will be set up *before* and *after*
   the call to PERT. 

The testcase 'pert-mmfp.inp' demonstrates the new functionality.


.. _pert_chemical:

Chemical
--------

Experimental support of "chemical paths" in PERT free energy simulations
(i.e., simulations in which a solute is decoupled from solvent or in
which a ligand is decoupled from the protein). CHARMM has to
be compiled with the CHEMPERT keyword in pref.dat for this functionality
to be available.

The CHEM keyword on the PERT commandline changes how energies and
forces are computed when PERT is active. (Note that there are no
provisions to turn this mode off other than specifying PERT OFF!)

Its intended use is as follows: Suppose you have a solute in SEGI SOLU
and water in SEGI WAT. Now specify

::

   fast off ! important, otherwise SCALAR RSCA has no effect!!

   PERT sele segi SOLU end CHEM

   scalar char set 0. sele segi SOLU end ! turn off solute charges
   scalar rsca set 0. sele segi SOLU end ! turn off solute LJ interactions

Then lambda=0 corresponds to a fully interacting system and lambda=1
corresponds to water plus the "solute" in the gas phase. However, within 
the SOLU segment (the solute) *all* interactions, in particular the
intramolecular nonbonded terms, are computed normally.

If the CHEM keyword were missing, then at lambda=1 no solute
intramolecular nonbonded interactions would be needed, requiring a
separate calculation correcting for this.

The code (should) work(s) with most nonbonded options, including EWALD
and PSSP soft-cores (note that these impose there own limitations
which nonbonded options can be used). The EWALD support is admittedly
a wild hack. Quite generally I do recommend to test correct
functionality by applying tests similarly to the new testcases
chempert1.inp and chempert2.inp adapted to one's particular system one
wishes to study.


.. _pert_lrcorrection:

Long-range correction
---------------------

Since c29 CHARMM supports an isotropic Lennard-Jones long range
correction (lrc nonbonded option, ##LRVDW keyword in pref.dat). So far,
there is no support for this correction in connection with PERT.
In principle, this support has now been added; both the energy, as
well as the virial correction are supported. HOWEVER, because
of the way the LRC interfaces with the energy routines, the correction
is NOT computed correctly if LJ interactions are turned off by 
scalar RSCA set 0. Thus, in the presumably most interesting
case, simulations where a solute or a ligand is decoupled from
the rest of the system (cf. CHEMical PERT above), the LRCorrection
cannot be computed in the most obvious manner. It is hoped that
full support can be added, but this will entail a rewrite of how
the LRC is computed in general.

::

   PERT-BLOCK Lambda Dynamics

New code is added in current version so that PERT and BLOCK 
lambda-dynamcs can be used together. Keyword QLDM of lambda-dynamcs
will trigger the mixture of PERT/BLOCK. No additional setup is required
as long as the traditional PERT and BLOCK setup is kept. FAST OFF is
recommended for this setup. The acdemic usage of this method will be
published shortly.  -- H. Li and W. Yang
