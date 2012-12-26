.. py:module:: block

=====
BLOCK
=====

The commands described in this section are used to partition a
molecular system into "blocks" and allow for the use of coefficients
that scale the interaction energies (and forces) between these blocks.
This has a number of applications, and specific commands to carry out
free energy simulations with a component analysis scheme have been
implemented. The lambda-dynamics, an alternative way of performing
free energy calculations and screening binding molecules, has also been
implemented.  Subcommands related to :chm:`BLOCK` will be described here.  To
see how to output the results of a dynamics run, please see :doc:`dynamc`
documentation (keywords are :chm:`IUNLDM`, :chm:`NSAVL`, and :chm:`LDTITLE`).
Please refer to :doc:`pdetail` for detailed description of the lambda
dynamics and its implementation.

:chm:`BLOCK` was recently modified so that it works with the :chm:`IMAGE`
module of CHARMM.  As some changes to the documentation were necessary
anyways, it was decided to also improve the existing documentation.
The :ref:`Syntax and Function <block_syntax>` section below are relatively unchanged;
the added documentation is in the :ref:`Hints <block_hints>` section (READ IT if you are using
BLOCK for the first time!).  Comments/suggestions to boresch@tammy.harvard.edu.

BLOCK was modified so that it works with the Ewald (simple and PME)
method of CHARMM. The Syntax and Function of BLOCK module are unchanged.

.. _block_syntax:

Syntax of BLOCK commands
------------------------

::

   BLOCk [int]

   Subcommands:

   miscellaneous-command-spec    !   see *note miscom:(chmdoc/miscom.doc).

   CALL int atom-selection

   LAMBda real

   COEFficient int int real -
         [BOND real] [ANGL real] [DIHEdral real] [ELEC real] -
         [VDW real] [VDWA real] [VDWR real]

   NOFOrce

   FORCe

   FREE_energy_evaluation  [OLDLambda real] [NEWLambda real] -
        FIRSt int [NUNIT int] [BEGIn int] [STOP int] [SKIP int] -
        [TEMPerature real] [CONTinuous int] [IHBF int] [INBF int]
        [IMGF int]

   INITialize

   CLEAr

   Energy_AVeraGe  [OLDLambda real] [NEWLambda real] -
        FIRSt int [NUNIT int] [BEGIn int] [STOP int] [SKIP int] -
        [CONTinuous int] [IHBF int] [INBF int] [IMGF int]

   COMPonent_analysis  DELL real NDEL int [TEMPerature real] -
        FIRSt int [NUNIT int] [BEGIn int] [STOP int] [SKIP int]
        [IHBF int] [INBF int] [IMGF] int

   AVERage {DISTance int int}
           {STRUcture}
        [PERT] [TEMPerature real] [OLDLambda real] [NEWLambda real] -
        FIRSt int [NUNIT int] [BEGIn int] [STOP int] [SKIP int]

   LDINitialize int real real real real [real]

   RMBOnd RMANgle

   LDMAtrix

   LDBI int

   LDBV int int int int real real int

   LDRStart

   LDWRite IUNL int NSAVL int

   RMLAmbda {internal_energy_spec}
              internal_energy_spec ::== BOND THETa|ANGLe PHI|DIHEd IMPHi|IMPR
   SAVE

   UNSAve

   QLDM  [THETa]

   QLMC  [MCTEmperature real] [FREQ int] [MCSTep int] [MAX real]

   MCIN  int {real .... real}

   MCDI  real

   MCRS

   MCLEar

   MSLD int_1 int_2 ... int_nblocks { FNEXponential [real] }
                                 { FNSIn }
                                 { F2Exponential }
                                 { F2Sin }
                 ! note: int_1 must be 0, block 1 = environment = Site 0

   MSMAtrix

   LANG [TEMP real]

   RSTP int real
   " Dual-topology Softcore"
   [PSSP]        ! use soft core potentials for interactions in between
                 ! blocks.  This option is remembered. With
                 ! the PSSP keyword, two parameters, ALAM and DLAM can
                 ! be set.

   [ALAM real]   ! Separation parameter for elec. interaction (defaults to 5A^2)

   [DLAM real]   ! Separation parameter for LJ interaction (defaults to 5A^2)

   [NOPSsp]      ! Turn off use of soft core interactions.
   " -- H. Li and W. Yang

   MCFRee EXFReq int FINI real FFIN real FLAT real

   MCLAmd int LAMD0 real LAMD1 real ....... LAMD[int-1] real

   HYBH real                     ! HYBrid Hamiltonan module (HYBH).

   OUTH int                      ! HYBH

   TSTH real [<update-spec>]     ! HYBH

   PRIN                          ! HYBH

   PRDH                          ! HYBH

   CLHH                          ! HYBH

   END


.. _block_function:

1) ``BLOCk [int]`` enters the block facility.  The optional integer is
   only read when the block structure is initialized (usually the first
   call to block of a run) to specify the number of blocks for space
   allocation.  If not specified, the default of three is assumed.

2) ``END`` exits the block facility.  The assignment of blocks, the
   coefficient weighting of the energy function, the force/noforce
   option, etc.  remain in place.  For the terms of the energy function
   that are supported, each call to :doc:`ENERGY <energy>` (either directly or through
   :doc:`MINIMIZE <minimiz>`, :doc:`DYNAMICS <dynamc>`, etc.  commands) results in an energy and force
   weighted as specified.  The matrix of interaction coefficients is
   printed upon exiting.

3) ``CALL`` removes the atoms specified by ``atom-selection`` from their
   current block and assigns them to the block number specified by the
   integer.  Initially all atoms are assigned to block 1.  If atoms are
   removed from any block other than block 1, a warning message is
   issued.  If blocks are assigned such that some energy terms (theta,
   phi, or imphi) are interactions between more than two blocks, a
   warning is issued when the ``END`` command is encountered.  (Take such
   warnings seriously; this is a severe error and indicates that
   something is wrong.  However, the problem might be not the ``CALL``
   statement (or the atom selection) itself; quite possibly your hybrid
   molecule was generated improperly)

4) ``LAMBda`` sets the value of lambda to "real".  This command is only
   valid when there are three blocks active.  Otherwise multiple COEF
   commands may be used to set the interaction coefficients manually.

   ::

      LAMBda x

   is equivalent to (let y=1.0-x)

   ::

      COEF 1 1 1.0
      COEF 1 2 y
      COEF 1 3 x
      COEF 2 2 y
      COEF 2 3 0.0
      COEF 3 3 x

5) ``COEF`` sets the interaction coefficient between two blocks (represented
   by the integers) to a value (the real number).  When the block facility
   is invoked, all of the atoms are initially assigned to block 1 and all
   interaction coefficients are set to one.  The required real value
   (first specified) scales all energy terms expect those specific terms
   which are named with alternative corresponding scale factors.

   The name :chm:`VDWA` and :chm:`VDWR` correspond to the Attractive and Repulsive
   terms in the Lennard-Jones potential respectively. That is they allow
   one to independently scale the attractive (:math:`r^6`) and repulsive terms (:math:`r^{12}`)
   independently.

6) ``NOFOrce`` specifies that in subsequent energy calculations, the
   forces are not required.  This is economical when using the
   post-processing commands (:chm:`FREE`, :chm:`EAVG`, :chm:`COMP`).
   Forces may be turned back on with the :chm:`FORCe` command; this is
   necessary before running minimizations and dynamics if there was a
   prior :chm:`NOFO` command.

7) ``FREE`` calculates a free energy change using simple exponential
   averaging, i.e. the "exponential formula".  If the old and new lambdas
   (:chm:`OLDL`, :chm:`NEWL`) are specified (can only be done when three blocks are
   active), the perturbation energy is calculated for these values (i.e.
   :chm:`FREE` gives you the free energy difference between NEWLambda and
   OLDLambda via perturbation from the lambda value at which your
   trajectory was calculated.  If not, the current coefficient matrix is
   used (:chm:`FREE` should be used with three blocks, and the use of :chm:`OLDL` and
   :chm:`NEWL` is recommended).  :chm:`FIRSt_unit`, :chm:`NUNIt`, :chm:`BEGIn`, :chm:`STOP`,
   and :chm:`SKIP` specify the trajectory/ies that is/are to be read (for a further
   description see the :chm:`TRAJ` command elsewhere in the CHARMM
   documentation).  :chm:`TEMPerature` defaults to 300 K and gives the
   temperature value to be used in :math:`k_B T`. :chm:`CONTinuous` specifies the
   interval for writing cumulative free energies.  A negative value
   causes binned (rather than cumulative average) values to be written.
   Be careful to make sure that you use correct non-bonded lists (see the
   hints section!)

8) ``INITialize`` is called automatically when the BLOCK facility is
   first entered and may also be called manually at some other point.
   All atoms are assigned to block one and all interaction coefficients
   are set to their initial value.

9) ``CLEAr`` removes all traces of the use of the BLOCK facility.  The
   next command should generally be :chm:`END`, and then CHARMM will operate
   as if BLOCK had not ever been called.

10) ``EAVG`` The average value of the potential energy during a simulation
    can be calculated with the :chm:`EAVG` (Energy_AVeraGe) command.  The parsing
    is very much like the FREE command above.  The most frequent use of
    this command is to calculate the average value of dV/dlambda during
    the course of a simulation for use in thermodynamic integration.
    :chm:`CONTinuous` specifies the interval for writing cumulative free
    energies.  A negative value causes binned (rather than cumulative
    average) values to be written.  Be careful to make sure that you use
    correct non-bonded lists (see the hints section!)  The command accepts
    the :chm:`OLDL` / :chm:`NEWL` option, similarly to :chm:`FREE`, but for
    :chm:`EAVG` it is recommended to set up the interaction matrix (using
    :chm:`COEF` commands) yourself -- see the hints section.

11) ``[COMP]`` The :chm:`COMP` module is essentially a modified version of the
    EAVG module which aside from calculating :math:`\left< dU/dl \right>  = \left< U_1 - U_0 \right>`
    at a given value of :math:`\lambda_i` will also give you expectation values of
    this quantity at :math:`\lambda_{i\pm1}`, :math:`\lambda_{i\pm2}`, etc. based on perturbation theory.
    :chm:`COMP` requires 4 blocks.  Put the usual WT (reactant) in block 2 and
    MUT (product) in block 3.  Put the portion of the environment whose
    contribution to the free energy change is desired into block 4 (this
    can be everything else, or just a subset) (Note that the same can be
    achieved easily with the :chm:`EAVG` command) You have to set up your own
    coefficient matrix.  Much of the parsing is like the :chm:`EAVG` command.
    :chm:`CONT` is not supported.  Two special subcommands (required) are DELL
    and :chm:`NDEL`.  The normal output of COMP is :math:`\left< U_1 - U_0 \right>`
    evaluated at the lambda of the simulation.  However, :chm:`COMP` also evaluates the same
    ensemble averages perturbed to ``lambda = lambda +/-
    {0,1,2,...NDEL}*DELL``.  This (sometimes) helps the quadrature in
    thermodynamic integration.  Note that :chm:`NDEL` must be at least 1, and
    :chm:`DELL` should not be zero.  (You have to specify these values; the
    default values will lead to an invalid input, i.e. you bomb...) Be
    careful to make sure that you use correct non-bonded lists (see the
    :ref:`hints <block_hints>` section!)  A word of warning: If your initial ensemble average
    (at the lambda of the simulation) is not well converged, then your
    perturbed values are most likely random numbers.  The approach taken
    by :chm:`COMP` is theoretically sound, but it should only be applied if
    convergence has been established!  The output format of :chm:`COMP` is
    somewhat messy: :chm:`COMP` first prints :math:`\left< dU/dl \right>  = \left< U_1 - U_0 \right>`
    at lambda =

    ::

		 lambda - NDEL*DELL
                 lambda - (NDEL-1)*DELL
                 ...
		 lambda
		 lambda + DELL
                 ...
                 lambda + NDEL*DELL;

    then it prints an average (integral) value over these results.  The
    meaning of this last value is unclear to me.  In earlier versions of
    this documentation, :chm:`COMP` has been recommended over :chm:`EAVG`.  In my
    experience the opposite is true.  There is little :chm:`COMP` can do which
    you can't do with EAVG (aside from obtaining expectation values for
    :math:`\left< dU/dl \right>`).  (Maybe the unclear output of the :chm:`COMP`
    module is the main reason why I don't like it).

12) ``[AVER]`` The :chm:`AVERage` command is used to extract ensemble average
    structural properties from a dynamics simulation.  Features in this
    implementation allow averages taken over ensembles that are perturbed
    from that which the simulation corresponds to.  This is particularly
    useful for calculating the average structure expected at lambda=0.0
    from a simulation run at lambda=0.1, for example.  One may calculate
    average structures ``[STRUcture]`` and average distances ``[DISTance int
    int; where the two integers are the atom numbers between which the
    average distance is requested]``, currently.  The :chm:`PERT` keyword indicates
    that a perturbed ensemble from the dynamics trajectory is desired,
    with :chm:`TEMPerature` giving the temperature to use in the exponential for
    the perturbation (defaults to 300 K), OLDLambda and NEWLambda are the
    lambdas for which the simulation was run and for which the ensemble is
    requested, respectively (only valid if three blocks are active; if
    these are not specified, the perturbation energy is calculated with
    the current coefficient matrix), and the remaining keywords are used
    to specify the trajectory.  NOTE: TO THE BEST OF MY KNOWLEDGE THIS
    COMMAND HAS NOT BE MAINTAINED (so you are on your own if you use it!)

13) ``LDINitialize`` specifies input parameters for running lambda
    dynamics.  It sets up the value of ``lambda**2``, the velocity of
    the lambda, its mass and reference free energy (or biasing potential).
    E.g, the following input lines set up
    parameters for the third lambda with ``[lambda(3)]**2 = 0.4``,
    ``lambdaV(3) = 0.0``, ``lambdaM(3) = 20.0``, and ``lambdaF(3)=5.0`` (note that ``lambdaF(1)``
    should always be set to zero).

    ::

      LDIN 3   0.4   0.0   20.0   5.0

    For more details, see :ref:`block_hints_lambda_dynamics`.

14) LDMAtrix will automatically map the input lambda**2 values onto the
    coefficient matrix of the interaction energies (and forces) between
    blocks.

15) LDBI provides an option on applying biasing potentials on lambda
    variables. The integer value specifies the total number of biasing
    potentials to be used. E.g,

    ::

      LDBI 3

    will include total of 3 biasing potentials in the simulation.

16) LDBV sets up the specific form of the biasing potentials. At the
    moment, the functional form is of power law and allows three different
    classes (for details see "the actual simulations"). The input format is

    ::

      LDBV INDEX  I   J  CLASS  REF  CFORCE NPOWER

    e.g.

    ::

      LDBV   2    2   3    3    0.0   50.0   4

    will assign the second biasing potential acting between lambda(2) and
    lambda(3). The potential form belongs to the third class with reference
    value of zero, the force constant of 50 kcal/mol and the power of four.

17) LDRStart is used to restart the lambda dynamics runs.

18) LDWRite specifies the FORTRAN output unit No. and the frequency
    for writing lambda histogram by assigning an integer to IUNL and an
    integer to NSAVL. (IUNL and NSAVL can be reset in DYNAmic command,
    see :doc:`dynamc`)

19) RMBOnd and RMANgle are used when no scaling of bond and angle energy
    terms is desired.

20) RMLA is used when no scaling of bond, angle, proper torsion, and
    improper torsion terms are desired. This option always works with block module.
    The keywords: "RMBOnd" and "RMANgle" work only in lambda-dynamics.

    COEF command can work in the same way when lambda-dynamics or hybrid-MC/MD are
    not used.

    e.g.

    "RMLA BOND" = "COEF real BOND 1.0"

21) SAVE saves the decomposed-energy file for post processing in the TSM
    module. This command gives a choice for free energy calculation with
    block module to get free energy without saving the trajectory file.
    The condition and the name for the decomposed-energy file can be defined
    in the dynamics module. (see dynamic.doc, keyword: IBLC, NBLC)

22) UNSAve removes the traces of the use of SAVE command shown above.

23) QLDM turns on lambda-dynamics option. LDIN command also turns on
the lambda-dynamics option only when QLMC turns off.

24) QLMC turns on hybrid-MC/MD option. If QLMC option is on, LDIN commands
    do not activate the QLDM option.

    In this version, we do not re-assign the velocity of the atoms when
    chemical variables (lambda) are changed by MC method. Therefore, the kinetic
    terms suddenly change into the different phase space. The stochastic dynamics
    may diminish such artificial effects and help to reach the canonical ensemble.
    QLMC and QLDM are exclusive and latest choice is active.
    QLMC command should specify conditions for hybrid-MC/MD.

    e.g.

    ::

      QLMC MCTEmperature 300.0 FREQ 10 MCST 5 MAX 0.9

    IN the above example, the temperature used for sampling the chemical space
    by MC method is 300.0 [Kelvin]; MC sampling works every 10 molecular dynamics
    steps (using for sampling of the atomic space); in one MC sampling, 5 trials
    are examined; the scale factor (lambda^2) for the selected ligand is assigned
    to 0.9 and the rest of ligands (L-1) have the scale factor 0.1/(L-1).
    Different temperature can be defined in the lambda-dynamics and hybrid MC-MD
    for atomic variables and chemical variables.

25) MCIN allows the intermediate states in which only two ligands have non-zero
    lambda values in hybrid-MC/MD method.

    e.g. (Three ligands system)

    ::

       MCIN 5 0.0 0.25 0.5 0.75 1.0

    5 means that each ligand may have one these five scalings:

    ::

       0.0, 0.25, 0.5, 0.75, and 1.0.

    In this condition, CHARMM recognizes the following chemical states:

    +----------+-------------------+
    |          |   (SCALE FACTOR)  |
    |          +------+-----+------+
    | STATE NO.| LIG_A|LIG_B|LIG_C |
    +----------+------+-----+------+
    |   1      | 1.0  |0.0  |0.0   |
    +----------+------+-----+------+
    |   2      | 0.0  |1.0  |0.0   |
    +----------+------+-----+------+
    |   3      | 0.0  |0.0  |1.0   |
    +----------+------+-----+------+
    |   4      | 0.25 |0.75 |0.0   |
    +----------+------+-----+------+
    |   5      | 0.75 |0.25 |0.0   |
    +----------+------+-----+------+
    |   6      | 0.25 |0.0  |0.75  |
    +----------+------+-----+------+
    |   7      | 0.75 |0.0  |0.25  |
    +----------+------+-----+------+
    |   8      | 0.0  |0.25 |0.75  |
    +----------+------+-----+------+
    |   9      | 0.0  |0.75 |0.25  |
    +----------+------+-----+------+
    |  10      | 0.5  |0.5  |0.0   |
    +----------+------+-----+------+
    |  11      | 0.5  |0.0  |0.5   |
    +----------+------+-----+------+
    |  12      | 0.0  |0.5  |0.5   |
    +----------+------+-----+------+

26) MCDI (increment) specifies the step size to move in lambda chemical
    movement. It allows intermediate states in which more than two ligands
    can have non-zero lambda values in hybrid-MC/MD method. "MCDI" requires the
    uniform interval for the definitions of the intermediate states.
    Step size must satisfy:

    ::

      Stepsize = 1.0/integer.

    Example: Three ligands system

    ::

      MCDI 0.25   ! 0.25 shows the step size to move in lambda chemical movement.

    In this condition, CHARMM recognizes next chemical states.

    +----------+-------------------+
    |          |   (SCALE FACTOR)  |
    |          +------+-----+------+
    | STATE NO.| LIG_A|LIG_B|LIG_C |
    +----------+------+-----+------+
    |   1      | 1.0  |0.0  |0.0   |
    +----------+------+-----+------+
    |   2      | 0.0  |1.0  |0.0   |
    +----------+------+-----+------+
    |   3      | 0.0  |0.0  |1.0   |
    +----------+------+-----+------+
    |   4      | 0.25 |0.75 |0.0   |
    +----------+------+-----+------+
    |   5      | 0.75 |0.25 |0.0   |
    +----------+------+-----+------+
    |   6      | 0.25 |0.0  |0.75  |
    +----------+------+-----+------+
    |   7      | 0.75 |0.0  |0.25  |
    +----------+------+-----+------+
    |   8      | 0.0  |0.25 |0.75  |
    +----------+------+-----+------+
    |   9      | 0.0  |0.75 |0.25  |
    +----------+------+-----+------+
    |  10      | 0.5  |0.5  |0.0   |
    +----------+------+-----+------+
    |  11      | 0.5  |0.0  |0.5   |
    +----------+------+-----+------+
    |  12      | 0.0  |0.5  |0.5   |
    +----------+------+-----+------+
    |  13*     | 0.25 |0.25 |0.5   |
    +----------+------+-----+------+
    |  14*     | 0.25 |0.5  |0.25  |
    +----------+------+-----+------+
    |  15*     | 0.5  |0.25 |0.25  |
    +----------+------+-----+------+


    It is possible for MCDI to produce a state in which  three ligands take
    non-zero lambda values as shown with the asterisk (states 13, 14 and 15).
    "MCDI" seems to be more general, but "MCIN" allows non-uniform
    intervals. Thus, small step sizes can be assigned near end points.

27) MCRS ignores the force for lambda coming from the restraining potential
    in lambda-dynamics. It also ignores the restraining potential energy when
    chemical space is sampled by MC method. CMC/MD (Chemical Monte Carlo &
    molecular dynamics) method can be carried out by combining this command
    with QLMC.

28) MCLEar removes the traces of the use of QLMC command shown above.
    BLOCK CLEAr command also removes the all traces of the use of QLMC.
    MCLEar removes the traces of QLMC, while BLOCK CLEar removes all traces of the
    BLOCK module.

29) LANG turns on the interaction between lambda variable and langevin
    heatbath.  In general, weak interaction between lambda variables and atoms
    produced large deviations from the target temperature. Different temperatures
    for lambda and atoms make nonequilibrium states and gave incorrect free
    energies.  Therefore, we recommend that LANG turn on in any lambda-dynamics
    simulations.  LEAP FROG integration method is required when using the LANG
    option.

30) RSTP adds the restraining potential for the unbound states ligands
    in lambda-dynamics and hybrid-MC/MD method to keep the physical low energy
    states.  The type of the restraining potential used with RSTP is;

    ::

      R = alpha *(1 - lambda^2)*  ( V - F )
       i                    i        i   i

    It disappears when this ligands is in bound state (lambda=1).

    e.g.

    ::

      REST 3 0.3

    3 means the type of the restraining potential; 0.3 shows the alpha value.

    There are three types for the restraining potential.

    * Type 1 Both environmental atoms and the ligands feel the restraining potential.
      Umbrella sampling technique is used to remove the bias effect coming
      from the restraining potential.

    * Type 2 The fixed average structure of the environmental atoms are assigned into
      Block 2. The restraining potential was calculated Ri is defined as a
      function of the fixed environmental atoms and the ligands.
      When the system is flexible and the difference between the real
      coordinates of the environmental atoms and fixed average coordinates
      are considerably large, the convergence tends to slow.

    * Type 3 When the environmental atoms form the specific structure and vibrated
      around the minimum, the fixed average structure of the environmental
      atoms are similar to those of the real time coordinates.
      Therefore, the force coming from the restraining potential can be
      approximated zero as an average.  If such a condition is satisfied,
      the environmental atoms can be ignored the force coming from the
      restraining potential and the ligands only feel the restraining
      potential.This approximation may have a problem when we handle the
      unstructured system like gas or liquid.

    The utility program, post_ldm_mcmd.exe is prepared for calculating the free
    energy differenes both without or with the restraining potential in
    lambda-dynamics or hybrid-MC/MD method.

    This program is saved in "support/post_analysis".

31) MCFRee EXFReq int FINI real FFIN real FLAT real is the main subcommand for
    the definition of simulated scaling simulations. Here, EXFReq int is to set up
    the frequency for Monte Carlo acceptance and rejection of the lambda space
    move. FINI real is to set up the initial modification factor, usually as
    2.71828 following the original Wang-Landau algorithm. FFIN real is to set up
    the cutoff value for the final modification factor. FLAT real is to set up
    the cutoff value for each cycle of flatness judgment.

    Reference: Li, H., Fajer, M., and Yang, W. 2007. Simulated scaling method for
    efficient localized conformational sampling and simultaneous alchemical free
    energy simulation: A general method for MM, QM, and QM/MM simulations.
    J. Chem. Phys. 126:024106.

32) MCLAmd int LAMD0 real LAMD1 real ...... LAMD[int-1] real is an additional
    facility for the flexible usage of the simulated scaling method. Here, [int]
    is to define the number of lambda values. LAMD0 is the first lambda value,
    LAMD1 is the second one, ...., LAMD[int-1] is the last one.

33) HYBH , HYBrid_Hamiltonian module. Implementation of the truncation
    scheme described in "Ensemble Variance in Free Energy Calculations by
    Thermodynamic Integration: Theory, Optimal "Alchemical" Path, and
    Practical Solutions", A.Blondel (2004) J.Comp.Chem 25, 985-993.
    Details on the method should be sought therein. In brief, the
    implementation is based on dual topology (although single topology
    could be used under some conditions), the bonded terms (bond, angle
    and Urey-Bradly) are kept unchanged, dihedral and impropers are
    modified according to simple quadratic scheme (w_product=(3.l+1).l/4),
    and electrostatic and van der Waals are treated together with a
    truncation scheme reminiscent of soft-core vdw to minimize the
    numerical fluctuations of the integrant (hence Optimal "Alchemical"
    Path). Ewald sums and correction terms associated appeared soft
    enough to be treated according to linear scaling of the charges,
    allowing direct analytical calculation of dEwald/dl. A benefit of the
    method, in addition to the fact that the integrant has limited
    numerical fluctuations, is that it also produce a linear evolution of
    the integrant along lambda (or l) in regular cases.

    The implementation attempts to supports most of non-bonds, image
    and Ewald sums options and warnings are made. Slow routines are not
    currently supported. However, it is advised to test the results when
    new combination of options are used. CMAP is not currently supported.

    Associated commands are called from within the BLOCk module and are:

    * ``HBYH real``: Switchs the module on and sets the lambda parameter.
      Due to the theoretical properties of the method, evenly spaced
      values should be sufficient (eg. l=(2i-1)/20). The product part
      (bloc 3) is weighted according to l as explained above, and the
      reactant part (bloc 2) is weighted according to (1-l) as explained
      above.

    * ``OUTH int``: Sets the output unit for the dE/dl terms.

    * ``TSTH real [<update-spec>]``: Sets dl and tests the derivatives (dE/dl)
      by finite differences (E(l+dl)-E(l-dl))/2/dl. None zero components
      of the energy are printed.

    * ``PRIN``: Prints dE/dl with the usual ENERGY printing format.

    * ``PRDH``: Writes dE/dlambda components to outh unit. Replaces the
      automatic writting performed during dynamics, for example, when
      re-reading a trajectory for post-processing.

      The current form of the output is formatted, two line per dynamic step.

      ::

        R l dDIHEr dIMDIHEr dVDWr dELECr dEWKSUMr dEWSELFr d(EWEXCL+EWQCOR+EWUTIL)r
        P l dDIHEp dIMDIHEp dVDWp dELECp dEWKSUMp dEWSELFp d(EWEXCL+EWQCOR+EWUTIL)p
        Format: (a1,1x,f6.4,7(1x,1pg24.16e2))

    * ``CLHH``: Clears the data structure for truncation scheme and switchs off
      the module without changing the rest of the block setup. Note, the
      BLOCk/CLEAr command also switchs off the module.

    No analysis routine is currently supplied as careful convergence
    analysis should be undertaken. It is advised that additions of the
    terms be made at least in real*8 format as truncation errors might
    be significant otherwise.

    Testcases c35test/block_hybh.inp & block_hybh_ew.inp are provided.

34) MSLD invokes Multi-Site lambda-dynamics. The integers which follow
    the keyword indicate the "Site" to which atoms within each block are
    assigned. The first block must be assigned to Site 0 (the "environment"
    atoms). Currently, QLDM THETA must be specified prior to invoking MSLD.

    Several different functional forms of lambda have been implemented. The
    default functional form is FNEX 5.5. (Note: these functions are for
    lambdas associated with all blocks except for block 1--ie. the environment
    atoms at site 0.)

    i) n-block normalized exponential: FNEX [c]

       ::

         num(Site_a,sub_i) = exp(c*sin(theta(Site_a,sub_i))


         lam(Site_a,sub_i) =    num(Site_a,sub_i)
                             -------------------------
                               ----
                               \
                               /    num(Site_a,sub_j)
                               ----
                               all j

    ii) n-block normalized sin: FNSI

        ::

          num(Site_a,sub_i) = sin(theta(Site_a,sub_i))^2

          lam(Site_a,sub_i) =    num(Site_a,sub_i)
                              -----------------------------
                                ----
                                \
                                /    num(Site_a,sub_j)
                                ----
                               all j

    iii) 2-block exponential: F2EX  (based on the logistic function)

         ::

           lam(Site_a,sub_1) = exp(theta(Site_a)) / [ 1.0 + exp(theta(Site_a)) ]

           lam(Site_a,sub_2) = 1.0 / [ 1.0 + exp(theta(Site_a)) ]

    iv) 2-block sin: F2SI   (based on constant pH-MD and theta-dynamics)

        ::

          lam(Site_a,sub_1) = sin(theta(Site_a))^2

          lam(Site_a,sub_2) = 1.0 - sin(theta(Site_a))^2

    The MSMA keyword is the Multi-Site lambda-dynamics equivalent to the
    LDMAtrix command and will automatically map the input lambda values
    onto the coefficient matrix of the interaction energies (and forces)
    between blocks.

    Assuming that groups of atoms have already been defined to correspond to
    "site1sub1" etc., here is an example of a Multi-Site lambda-dynamics
    setup in an input file.

    ::

      BLOCK 7
          Call 2 sele site1sub1 end
          Call 3 sele site1sub2 end
          Call 4 sele site2sub1 end
          Call 5 sele site2sub2 end
          Call 6 sele site2sub3 end
          Call 7 sele site2sub4 end
          qldm theta
          lang temp 310.0
          ldin 1 1.0  0.0  12.0  0.0 5.0
          ldin 2 0.50 0.0  12.0  0.0 5.0
          ldin 3 0.50 0.0  12.0  3.2 5.0
          ldin 4 0.30 0.0  12.0  0.0 5.0
          ldin 5 0.40 0.0  12.0 -0.5 5.0
          ldin 6 0.15 0.0  12.0  8.5 5.0
          ldin 7 0.15 0.0  12.0 15.1 5.0
          rmla bond thet
          msld 0 1 1 2 2 2 2 fnex 5.5
          msma
      END

    After this setup, minimizations and dynamics can be invoked as usual. MSLD
    is currently only compatible with the default dynamics routine (leapfrog
    Verlet) and can be used with Langevin dynamics (LANG) using the LEAP
    integrator.

    Analysis of the generated lambda trajectories can be performed using
    options in the trajectory command for multiple blocks at one or two Sites
    (see TRAJ LAMB in dynamc.doc). For hybrid molecules that have multiple
    blocks at more than two Sites, we suggest running the TRAJ LAMB command
    with the "print" option to write out lambda and theta values at each step.

    Currently, Multi-Site lambda-dynamics is compatible with LDBI and
    LDBV. However, the LDBV defined biases are not yet taken into
    account in the TRAJ analysis routine.


.. _block_hints:

HINTS
-----

A warning is in order: the BLOCK module is quite user-unfriendly, AND
the user (=you) has to know what he/she is doing, otherwise you won't
get anywhere!  (Of course, this could be a blessing in disguise) There
are two applications for BLOCK: (i) Mere use as an energy partitioning
facility, which may, e.g., very helpful as an alternative to the
INTEraction energy command and (ii) use in free energy simulations.
The focus here is on free energy applications.  The following paragraphs
assume that you are familiar with the theory of free energy difference
simulations (e.g. Brooks et al. Advances in Chem. Physics, Vol. LXXI,
1988, chapter V); the emphasis here is to show how a rough tool as
BLOCK can be used to implement the theory in a program and (of course)
how to use it.

Using BLOCK in order to calculate a free energy difference consists
out of two rather dissimilar parts (as far as practical problems are
concerned): (i) Run your system at various values of lambda and save
trajectories. (ii) Postprocess the trajectories with the FREE or the
EAVG command (possibly COMP), use the quantities which these modules
give you to calculate the free energy difference.

(i) The actual simulations
^^^^^^^^^^^^^^^^^^^^^^^^^^

It's probably easiest to use a concrete example, and the free energy
difference between ethane and methanol in aqueous solution is used for
that purpose.  BLOCK is a so-called dual topology method (D. Pearlman,
JPC 1994, 98, 1487) i.e. one has to duplicate any atom that is
different with respect to any of its parameters.  In the
ethane/methanol case this means that you have to run with a solute
which looks something like

::

            H1
               \             /H4
               \  C1E ---- C2-H5
            H2 = {   }       \H6
               /  C1M --- OG
               /            \HG1
            H3


(and there is water.)

Conceptually, this system is divided into three regions:

-	environment: water, H1, H2, H3 (the region where nothing changes)
-	reactant: C1E, C2, H4, H5, H6 (ethane half)
-	product: C1M, OG, HG1 (methanol half),

where of course the role of reactant and product is interchangeable.

The steps involved to start running dynamics are as follows:

(1) set up the hybrid (generate psf).  In principle straightforward,
    but there is a practical pitfall: The autogenerate angles and
    dihedrals option(s) may produce artificial dihedrals/angles between
    the two/three parts of the system, e.g. you don't want angles
    H1-C1E-OG etc. or dihedrals H3-C1M-C2-H4 etc.  Also, make sure to
    specify nonbonded exclusions between the reactant and product part,
    otherwise you'll get endless distance warnings and may even bomb if
    two atom positions coincide.

(2) Place the hybrid into water (stochastic or periodic boundary
    conditions -- yes, IMAGE is now supported) as usual

(3) Partition the system, i.e. enter BLOCK
    The following script fragment will do the trick:

    ::

   	block 3
   	call 2 sele <reactant> end
   	call 3 sele <product> end
   	end

    (reactant and product have to be defined according to your system).
    BLOCK 3 initializes the block module with 3 blocks, all atoms are in
    block 1.  The two CALL commands bring the reactant and product part of
    the system into block 2 and 3 respectively.

(4) Run the necessary MD simulations.  Let's assume that you decide to
    use the following values of lambda, lambda = 0.1, 0.3, 0.5, 0.7, 0.9.
    You want to start your simulation at lambda = 0.1 and you have already
    partitioned your system as shown in (3).  (This information is kept
    within the same script between calls to block, but it is not saved in
    restart files or the psf, i.e. you have to repeat this step (as well
    as step (3)) in every input file).  Enter block again, e.g.

    ::

   	block
   	lamb 0.1
   	end

    From now on interactions between the 3 blocks will be scaled according
    to the following matrix (lambda = l = 0.1 ==> 1-l = 0.9):

    =====  === === ====
    block   1   2   3
    =====  === === ====
    1      1.0 1-l  l
    2      1-l 1-l  0.0
    3      l   0.0   l
    =====  === === ====

Please note that BLOCK will first calculate an interaction, then check
to which block the two atoms belong and scale the energy (and forces)
appropriately.  Therefore, if the distance between 2 atoms is zero
(e.g. in the ethane/methanol example I would define C1M and C1E on top
of each other!) then you need non-bonded exclusions, otherwise you
encounter a division by 0 error!

The LAMB command is a shortcut for the following sequence of COEF
commands, the following code fragment should be self-explanatory:

::

	block
	coef 1 1 1.0
	coef 1 2 0.9
	coef 1 3 0.1
	coef 2 2 0.9
	coef 2 3 0.0
	coef 3 3 0.1
	end

BLOCK only accepts and uses symmetric matrices, i.e. it doesn't
matter whether you specify COEF 1 2 or COEF 2 1.

Whenever you now call the energy routines, the energies/forces
returned from them will be scaled according to the matrix you have set
up.  Minimizers and Dynamics can be used as always.  So you are ready
to run dynamics, and for arguments sake say that you run at every
value of lambda 10,000 steps equilibration and 20,000 steps production
(i.e. you save coordinates to trajectories) You don't need to save
every step, every 5th to 20th step is probably more than enough.  (If
you saved every step you'd obtain highly correlated data, i.e. you
have larger trajectories, but you won't gain anything in terms of
convergence.)

(ii) Post-processing -- how to obtain a free energy difference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At this point in our example, you would have five trajectories
corresponding to lambda = 0.1, 0.3, ..., 0.9 The BLOCK module now has
to be used to obtain the average quantities you need for either the
exponential formula (FREE) or thermodynamic integration (EAVG,COMP)
from the trajectories you generated in step (i)

(1) At this point, some issues regarding the non-bonded list have to
    be considered.  No special considerations were necessary while running
    dynamics (aside from having some non-bonded exclusions where
    necessary); you just set up list updates as usual.  During
    post-processing there are two considerations: (a) efficiency -- you
    just want to calculate the necessary subset of interactions (otherwise
    your post-processing run will take about as much time as the
    simulation itself), and (b) proper list-updating.

    (a) Efficiency: In none of the post-processing routines do you need
        the interactions between particles that belong to the environment;
        therefore you should avoid calculating them.  This can be done easily
        by specifying

        ::

	       cons fix sele <environment> end

        Note that this is not necessary, but it will reduce the CPU time
        necessary from hours to minutes (and results are identical!)  However,
        if you had atoms belonging to reactant or product or both FIXed during
        the simulations in step (i), you MUST NOT FIX them now; otherwise
        you'll omit contributions.

    (b) List updating: While the efficiency considerations in principle
        are optional, you have to follow one of the two strategies below
        otherwise you'll get erroneous results.  If you used IMAGE, you have
        to use the second protocol!  Originally, the BLOCK post-processing
        commands would not do any list updating.  This meant that you had to
        have a nonbonded list which would include all possible interactions
        before starting post-processing -- don't forget that you post-process
        over, e.g., 20 ps and particles will move quite far.  You can easily
        create such a nonbonded list by specifying a CUTNB value of, e.g. 99.
        or 999. Ang (surely, all possible interactions will be included).  A
        CHARMM script looks approximately as follows:

        ::

         	!set up system (psf, initial coordinates)
         	block
         	!partition system
         	end
         	cons fix sele <environment> end
         ==>	energy cutnb 99. <all other options as during dynamics>
         	!open trajectories
         	block
         	!postprocessing
         	end

        In this case, do not use the inbf, ihbf and imgf options of the
        post-processing commands, they will default to 0, i.e. no update.
        This approach, however, CANNOT work with IMAGES!  Proper use of IMAGEs
        requires that the minimum image convention is checked periodically,
        i.e. particles have to be repartitioned between primary and image
        region.  As the BLOCK post-processing commands now understand INBF,
        IHBF and IMGF, this doesn't pose a problem.  However, the automated
        update is not supported (if you specify a negative value, you'll get a
        mild warning and the system will default to +1), and I recommend that
        you use 1 for all frequencies (don't forget, the frames in your
        trajectory are several steps apart, i.e. in general an update may be
        necessary)  The above scheme now looks like:

        ::

         	!set up system (psf, initial coordinates)
         	block
         	!partition system
         	end
         	cons fix sele <environment> end
         	! set up images if needed
         ==>	energy <all options, incl. CUTNB,  as during dynamics>
         	!open trajectories
         	block
         	eavg <other options> inbf 1 ihbf ? (imgf 1)
         	end

        Unless you have explicit hbond terms, ihbf can of course be 0!
        (Please note that there may or may not be problems with CRYSTAL, see
        Limitations section)

(2) The actual post-processing commands.  In the following I'll
    explain how to set things up for FREE, EAVG and COMP (as well as why).
    To speed up things further, you'll also want to specify the NOFOrce
    option at some point.

    FREE: This module allows you to calculate the necessary ensemble
    average for the exponential formula.  Using our example, you can for
    example estimate the free energy difference between l=0.1 (a value at
    which you ran a trajectory) and l=0.0, or, based on your l=0.1
    trajectory the free energy difference between l=0.0 and 0.2 (double
    wide sampling), i.e.

    ::

      A(0.0)-A(0.1) = -k_B*T*ln <exp[-(U(l=0.0)-U(l=0.1))/kT]>_(l=0.1)

    or

    ::

      A(0.2)-A(0.0) = -k_B*T*ln <exp[-(U(l=0.2)-U(l=0.0))/kT]>_(l=0.1)

    You should set up your system with 3 blocks and the usual environment,
    reactant and product partitions.  Before entering block to issue the
    free command, you have to open the trajectory/ies.

    ::

   	! all the stuff shown above for non-bond lists
   	open file unit 10 read name dat01.trj
   	block
   	free oldl 0.1 newl 0.0 first 10 nunit 1 [temp 300. -
   		inbf 1 imgf 1]
   	end

    or, for double wide sampling, the free line would be replaced by

    ::

   	free oldl 0.0 newl 0.2 first 10 nunit 1 [temp 300. -
   		inbf 1 imgf 1]

    Here dat01.trj is the trajectory which contains your 20 ps of dynamics
    at lambda = 0.1.  Based on the oldl/newl values (which correspond to
    A(newl) - A(oldl)), FREE generates the appropriate interaction matrix,
    which it prints; I recommend that you try to understand why it
    generates this matrix! FIRST is the unit number of the first
    trajectory file (10 in our example), NUNIT is the number of
    trajectories (1 in our example).  These (and the other options
    regarding the trajectories work as in any other post-processing
    command in CHARMM, see e.g. the TRAJ command) The update frequencies
    are optional depending on how you decided to handle your non-bonded
    updates.  temp defaults to T=300 K, cf. equations above.

    If you specify CONT +n, you'll get a cumulative average every n steps;
    in this case the last value equals the final result; if you specify CONT
    -n, you'll get the average over every n frames, plus of course the
    final result at the end.

    Note that trajectories are not rewound after use; i.e. before any
    subsequent FREE (or EAVG,COMP) command you have to rewind (or reopen)
    them!

    Once you have all the free energy pieces you need, you simply add them
    up to obtain the free energy difference (beware of sign errors
    depending on how you defined oldl/newl)

    EAVG: The main use of this module lies in obtaining the required
    ensemble averages for thermodynamic integration.  The most significant
    difference to EAVG is that you have to specify your own interactions
    matrix.  BLOCK uses linear coupling in lambda in the potential energy
    function, i.e.

    ::

	   V(l) = V0 + (1-l)*V_reac + l*V_prod,

    where V0 contains all the intra-environment terms, V_reac are the
    intra-reactant and reactant-env. interactions, and V_prod are the
    intra-product and product-env. interactions, respectively.  The
    quantity of interest in TI is dV/dl; for the above potential energy
    function we have

    ::

	   dV/dl = V_prod - V_reac

    It's very easy to obtain this quantity from EAVG.  Use 3 blocks,
    partition the system as before.

    ::

   	! all the stuff shown above for non-bond lists
   	open file unit 10 read name dat01.trj
   	block
   	coef 1 1  0.
   	coef 1 2 -1.
   	coef 2 2 -1.
   	coef 1 3  1.
   	coef 2 3  0.
   	coef 3 3  1.
   	eavg first 10 nunit 1 [inbf 1 imgf 1 cont +-n]
   	end

    You will calculate the average interaction energy over all the frames
    in the trajectory according to the following (symmetric) matrix

    ::

	     0.0
       -1.0  -1.0
        1.0   0.0  1.0;

    i.e. it's easy to see that the above script will give you <V_prod -
    V_reac>_(l=0.1).  If you post-process the other trajectories (l=0.3,
    0.5, ..,0.9) in an analogous fashion, you just have to approximate the
    TI integral by the trapezoidal formula (for basic Newton Cotes
    formulae (open and closed) see, e.g., Numerical Recipes), i.e. in this
    case you would have

    ::

	   dA = 0.2 * (dV(0.1)+dV(0.3)+...+dV(0.9)),

    where dV(0.1) = <V_prod - V_reac>_(l=0.1), etc.

    The above is an example of the basic use of EAVG.  You automatically
    get the formal components according to interaction type.  Cont +-n
    works similarly to the FREE case.  If you wanted to exclude the
    intramolecular contributions from ethane and methanol you could set up
    a slightly different coefficient matrix, i.e.

    ::

   	coef 1 1  0.
   	coef 1 2 -1.
   	coef 2 2  0.
   	coef 1 3  1.
   	coef 2 3  0.
   	coef 3 3  0.

    and you'll get only the solute-solvent contributions.  You can use
    more blocks (m > 3) to extract only a subset of interactions, e.g.

    ::

   	block 1: environment - x
   	block 2: reactant
   	block 3: product
   	block 4: x,

    where x is the region of interest, e.g. a specific sidechain in a
    protein (but not the one that is mutated!)

    Using EAVG with an appropriate coefficient matrix, e.g.

    ::

   	coef 1 1  0.
   	coef 1 2  0.
   	coef 1 3  0.
   	coef 1 4  0.
   	coef 2 2  0.
   	coef 2 3  0.
   	coef 2 4 -1.
   	coef 3 3  0.
   	coef 3 4  1.
   	coef 4 4  0.

    will give you (after integration over lambda) the free energy
    contribution of the interaction of sidechain x with the mutation site.
    Note that such formal free energy components may be (strongly)
    path-dependent.  These last two examples have hopefully provided a
    flavor of what can be done with the EAVG module.

    COMP: This module is also used for thermodynamic integration.  It
    always operates with four (and only four) blocks, just as the advanced
    example last given for EAVG, so it facilitates COMPonent analysis.
    Here I want to focus on the second unique aspect of COMP, it's
    capability to extrapolate additional datapoints, and so I consider in
    the framework of our ethane/methanol example the "special" case where
    I want the total free energy difference (as before in EAVG).  In order
    to do this, the system needs to be partitioned as follows

    ::

   	block 1: --
   	block 2: reactant
   	block 3: product
   	block 4: environment

    Whereas EAVG gave us <V_prod - V_reac>_l only for those lambda values
    at which we had actually done the simulations, COMP gives us
    additional values via perturbation (see Bruce Tidor's thesis).  Using

    ::

   	! all the stuff shown above for non-bond lists
   	open file unit 10 read name dat01.trj
   	block
   	coef 1 1  0.
   	coef 1 2  0.
   	coef 1 3  0.
   	coef 1 4  0.
   	coef 2 2 -1.
   	coef 2 3  0.
   	coef 2 4 -1.
   	coef 3 3  1.
   	coef 3 4  1.
   	coef 4 4  0.
   	comp first 10 nunit 1 [inbf 1 imgf 1] dell 0.06667 ndel 1
   	end

    will now give us <V_prod - V_reac>_l at l=0.03334, l=0.1 and
    l=0.16667.  If we use the same script on the other trajectories, we
    have 15 instead of 5 datapoints for the integration, i.e. we can
    obtain dA as

    ::

	   dA = 0.06667 * (dV'(0.03334)+dV(0.1)+...+dV'(0.96667)),

    where dV(0.1) = <V_prod - V_reac>_(l=0.1), etc. and the ' indicates
    that this is a perturbed quantity.   In principle, this
    should give a better numerical integration; however, in practice
    everything depends on how well your actual data (l=0.1, 0.3, ...,0.9)
    are converged.

    There is no check whether your ndel/dell combination is meaningful;
    and you cannot run COMP without using the perturbation feature, i.e.
    NDEL should be set to at least 1 (valid values are 1 through 5).  The
    defaults (if you don't specify ndel/dell) lead to an invalid input
    (This should be fixed...)

.. _block_hints_lambda_dynamics:

(iii) Lambda-dynamics simulations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In an efforts to make the transition from using previous subcommands
to running the lambda dynamics as smoothly as possible, we purposely
parallel new syntax after the COEF subcommand.  There are
total of eight new keywords for setting up new dynamics.  They are
classified according to their functionality.

(a) LDINitialize and LDMAtrix

    These two keywords are basic commands for starting the lambda
    dynamics run.  The correct use of them is tied together with the BLOCK
    and CALL commands.  Using the same example as the one given in "the
    actual simulations", the input script fragment will be as following:

    ::

        block 3
        call 2 sele <reactant> end
        call 3 sele <product> end
        LDIN 1    1.0    0.0    20.0    0.0
        LDIN 2    0.9    0.0    20.0	0.0
        LDIN 3    0.1    0.0    20.0	0.0
        LDMA
	     end

    Here, the LDINitialize command models after the COEF command with
    the format

    ::

        LDIN  INDEX   LAMBDA**2   LAMBDAV   LAMBDAM   LAMBDAF

    Several comments are in order.  First, notice that [lambda(1)**2]
    = 1.0 and [lambda(2)]**2 + [lambda(3)]**2 = 1.0.  They are quite
    similar to the inputs of COEF subcommand.  However, since one
    index instead of a pair is required here,  only diagonal elements
    of the interaction coefficient matrix are specified.  To fill up
    the matrix, LDMA is provided to finish the job automatically.
    In general, if there is total of N blocks, the first one is
    by default assumed to be the region where nothing changes.
    Therefore, [lambda(1)**2] = 1.0 is always true. The condition

    .. math::
       :label: 1

       \sum_{i=2}^N \lambda(i)^2 = 1.0


    has to be satisfied for the partion of the system Hamiltonian.
    Due to some technical reasons in our implementation (details
    see :doc:`pdetail`), we have used [lambda(i)**2] instead of lambda(i)
    in our partion of the system Hamiltonian.  Next, to make sure the above
    condition is met at any given simulation step, we have also enforced a
    condition containing velocities of the lambda variables

    .. math::
       :label: 2

       \sum_{i=2}^N \lambda(i)*\lambda_V(i) = 0.0

    We used lambdaV(i) = 0.0 in the above script just to simplify the
    input.  As far as the mass parameter lambdaM is concerned, the minimum
    requirement is that the value of mass has to be chosen such that the
    time step (or frequency) of lambda variables is consistent with that
    used for spatial coordinates x, y, z.  Since the lambda variable is
    introduced into the system by using extended Lagrangian,
    considerations gone into the similar quantities, such as the
    adjustable parameter Q in a Nose thermostat are applicable to the
    choice of lambdaM.  Some crude estimation can be made by examining
    the derivative of the system Hamiltonian with respect to the
    lambda, the curvature (simple harmonic approximation) or energy
    difference between two end-point states (0 and 1).  Our experience
    has indicated that a conservative choice of the mass, i.e. a little
    bit heavier mass than that of the crude estimate, serves us well
    so far.

    The biasing potential LAMBDAF has two functions: (1) In the screening
    calculations LAMBDAF corresponds to the free energy difference of the
    ligands in the unbound state. Such calculations can identify ligands
    with favorable binding free energy and a ranking of the ligands can be
    obtained from the probability of each ligand in the lambda=1 state;
    (2) In precise free energy calculations, LAMBDAF corresponds to the best
    estimate of free energy from previous calculations. Therefore the
    estimate of free energy can be improved iteratively.


(b) LDBI and LDBV

    In order to provide better control over simulation efficiency and
    sampling space, an option of applying biasing (or umbrella)
    potentials is furnished.  LDBI specifies how many biasing
    potentials will be applied and LDBV supplies all the details.
    The general input format is

    ::

       LDBV INDEX  I   J  CLASS  REF  CFORCE NPOWER


    Let us look at the following script

    ::

       block
       LDBI   3
       LDBV   1    2   2    1    0.2   40.0   2
       LDBV   2    3   3    2    0.6   50.0   2
       LDBV   3    2   3    3    0.0   20.0   2
       end

    It states that there is total of 3 biasing potentials. The first one
    (INDEX = 1) is acting on lambda(2) itself (I = J = 2), the second one
    on lambda(3) and the third one is coupling lambda(2) and lambda(3)
    together.  At the moment, five different classes of functional forms
    are supported:

    CLASS 1:

    .. math::

       V = \begin{cases}
          \mathrm{CFORCE} \cdot (\lambda - \mathrm{REF})^\mathrm{NPOWER} , & \text{ if } \lambda < \mathrm{REF}  \\
          0, & \text{ otherwise }
       \end{cases}


    CLASS 2:

    .. math::
       V = \begin{cases}
          \mathrm{CFORCE} \cdot (\lambda - \mathrm{REF})^\mathrm{NPOWER} , & \text{ if } \lambda > \mathrm{REF}  \\
          0, & \text{ otherwise }
       \end{cases}

    CLASS 3:

    .. math::
       V = \mathrm{CSFORCE} \cdot [\lambda(I) - \lambda(J)]^\mathrm{NPOWER}

    CLASS 4:

    ::
               __
              |   CFORCE*(1.0 - ((lambda - REF)**2)/REF**2) if lambda < REF
          V  =|
              |   0                                         otherwise
              |__


    CLASS 5:

    ::

          V  =    CFORCE*lambda(I)

    .. note::
       the CLASS 5 biasing potential is the same as invoking the
       biasing potential LAMBDAF in LDIN (except these biases will not
       currently be taken into account in the TRAJ analysis routines).


(c) LDRStart

    LDRStart is used only if for some reason, e.g. execution of EXIT command,
    the logical variable QLDM for the lambda dynamics has been set to false.
    In this case, to restart the dynamics, LDRStart can be used to reset
    QLDM = .TRUE..  However, if LDIN is also being used in restarting the
    dynamics, it will automatically reset QLDM. Therefore, LDRS does not
    need to be called in this case.


(d) LDERite

    LDWRite provides specifications for writing out lambda dynamics, i.e.
    the histogram of the lambda variables, the biasing potential etc.  The
    integer variable IUNLdm is the FORTRAN unit on which the output data
    (unformatted) are to be saved. The value of the integer NSAVL sets step
    frequency for writing lambda histograms.  IUNLdm is defaulted to -1 and
    NSAVL is defaulted to 0.  Both IUNLdm and NSAVl can be reset in DYNAmics
    command (Please refer to :doc:`dynamc` for details).

    the following script will set IUNLdm with unit No. 8 and NSAVL equal to 10:

    ::

      LDWRite IUNL 8 NSAVL 10


(e) RMBOnd and RMANgle

    Since each energy term is scaled by lambda, RMBOnd and RMANgle can prevent
    bond breaking caused by such scaling during dynamic simulations. Alternatively
    one can fix bonds (and angles) using SHAKE. But is is not always possible.

(f) MSLD

    Multi-Site lambda-dynamics is a generalized version of the original
    lambda-dynamics. Greater numerical stability of the simulations
    is acheived with the MSLD definitions of lambda which implicitly
    satisfy the constraints a) that each lambda value varies between 0 and 1 and
    b) that the lambda values for a given Site sum to 1 (see the functional
    forms listed above). Any system set up for the original lambda-dynamics
    (i.e. that has multiple blocks at only one Site) can be run using MSLD.
    In this case, the system would be set up in BLOCK as before, but the LDMA
    command would be replaced by the MSLD commands.

    For example, the original lambda-dynamics, using the theta-dynamics option
    (qldm test) setup would be:

    ::

        BLOCK 4
            Call 2 sele site1sub1 end
            Call 3 sele site1sub2 end
            Call 4 sele site1sub3 end
            qldm theta
            lang temp 310.0
            ldin 1 1.0  0.0  12.0  0.0 5.0
            ldin 2 0.50 0.0  12.0  0.0 5.0
            ldin 3 0.20 0.0  12.0  3.2 5.0
            ldin 4 0.30 0.0  12.0 -1.0 5.0
            rmla bond thet
            ldma                    ! use for original lambda-dynamics
        END

    and the MSLD setup would be:

    ::

        BLOCK 4
            Call 2 sele site1sub1 end
            Call 3 sele site1sub2 end
            Call 4 sele site1sub3 end
            qldm theta              ! required for MSLD
            lang temp 310.0
            ldin 1 1.0  0.0  12.0  0.0 5.0
            ldin 2 0.50 0.0  12.0  0.0 5.0
            ldin 3 0.20 0.0  12.0  3.2 5.0
            ldin 4 0.30 0.0  12.0 -1.0 5.0
            rmla bond thet
            msld 0 1 1 1 fnex 5.5   ! use for MSLD
            msma                    ! use for MSLD
        END

    Lambda trajectory files written by MSLD can be analyzed by TRAJ LAMB
    commands.  The header contains all the information required to process
    the trajectory (e.g. number of blocks, which blocks are assigned to
    which site etc.). The lambda trajectory files are specified in the
    DYNAMICS commands using keywords:

    ::

       IUNLDM unit ! where unit corresponds to the unit number of the
                     lambda trajectory file
       NSAVL freq  ! where freq corresponds to the frequency of writing
                     the lambda values

    The TRAJ LAMB command will process the lambda trajectory file and print
    out statistics related to individual sites ("single-site" statistics):

    * the population of each block (population = the number
      of snapshots in which each block(i) has lambda(i) = 1, or more
      specifically, the number of snapshots in which each block(i) has
      lambda(i) > threshold).

    * the number of transitions at each Site (i.e. the number of times
      the identity of the block with lambda(i) > threshold changes).

    * and the relative free energies for each pair of blocks at each Site.
      (without and with the correction for the fixed lambda biased invoked in
      the LDIN command)

    Output is provided for two threshold values (default 0.8 and 0.9) for
    approximating lambda(i) = 1 to provide an estimate of the sensitivity of
    the results to the specific threshold used:

    ::

             lambda(i) = 1, if lambda(i) > threshold

    For systems with more than one site (i.e. sites at which multiple blocks
    are modeled), a complete physical ligand is present at a given snapshot
    when there is a block with lambda > threshold at each Site. For a given
    system, there are a total of N(site_1) x N(site_2) x ... N(site_n)
    possible ligands where N(i) is the number of blocks at Site i.

    For systems with two sites, in addition to the general "single-site"
    statistics, the TRAJ LAMB command will account for all combinations of
    the blocks and print out "multi-site" statistics:

    * the populations of each "ligand" for two thresholds (population =
      the number of snapshots in which each "ligand" exists, i.e. the
      combination of blocks corresponding to the ligand each have lambda = 1)

    * the number of transitions between these ligands * the relative free
      energies of each pair of ligands (without and with the correction for
      the fixed lambda biased invoked in the LDIN command)

    For systems with more than two sites, it is recommended that you use
    the TRAJ LAMB PRINT command to print out the lambda values for each
    snapshot and perform the population analysis and compute the relative
    free energies yourself.

    See TRAJ LAMB in dynamcs.doc for a complete list of options.
    E.g.:

    * To read header information only:

      ::

        open unit 24 read file name scratch/msld_prod.lmd
        traj lamb query unit 24
        close unit 24

    * To process the trajectory file and print out lambda values at each
      timestep:

      ::

        open unit 24 read file name scratch/msld_prod.lmd
        traj lamb print first 24 nunit 1
        close unit 24

    While the trajectory is being processed, the following internal variables are
    stored:

    ========== =======================================================================
    ``TMIN``   Minimum number of transitions for any site in the system
    ``TMAX``   Maximum number of transitions for any site in the system
    ``FPL``    Fraction of the snapshots which represent full Physical Ligands
    ``POP#``   Population for the substituent associated with indicated BLOCK number
               at the low threshold value (e.g. the ?pop2 contains the population for
               substituent in BLOCK 2 given CUTLO threshold)
    ``DDG#_#`` Relative free energy between the first and second substituents
               listed at the low threshold value (e.g. ?ddg2_5 is the relative free
               energy between the substituents associated with BLOCKS 2 and 5).
    ========== =======================================================================

    If for any reason you wish to suppress the storage of internal variables
    (for example, if you have many substituents in your system and alreadyt
    many internal variables have been stored such that processing the MSLD
    trajectory gives a fatal error indicative of too many variables) then
    include the keyword "nosub" in the trajectory command, i.e.:

    ::

       open unit 24 read file name scratch/msld_prod.lmd
       traj lamb print first 24 nunit 1 nosub
       close unit 24


.. _block_limitation:

Limitation
----------

(1) Please be advised (again) that the AVERage command is unsupported,
    and I would not be surprised if it does not work (anymore).  Unless
    someone who understands this module better than I do maintains it, I
    recommend that we remove it.

(2) BLOCK now coexists with IMAGE "peacefully" and essentially
    transperantly to the user.  It works correctly for the case of a
    periodic water-box (cf. the block3.inp testcase).  I would, however,
    check carefully whether things really work before I would use it on
    something fancier like infinite alpha helices.  Similarly, it is not
    clear to me whether things work with the CRYSTAL facility.  If one
    modifies block3 as to use CRYSTAL instead of IMAGE things (seem to)
    work.  On the other hand, I know that I didn't support XTLFRQ in the
    post-processing routines as I don't understand its meaning.  I'll fix
    things if someone is willing to help me with the bits and pieces I
    don't understand.

(3) Bond and bond angle terms (including Urey-Bradleys).  Be advised
    that if you run a simulation at lambda = 0 or lambda = 1 you may
    effectively remove bond (and bond angle terms) as they get scaled by
    zero.  In other words, you would have ghost particles that can move
    freely through your systems, and this leads to all sorts of nasty
    side-effects.  Furthermore, this approach is not sound theoretically
    (S. Boresch & M. Karplus, unpublished).  So in general, avoid running
    at lambda = 0 and 1.  If you have your bonds constrained you're safe
    as the constraint will keep things together (that won't take care of
    angles however!)  In order to avoid artifacts from noisy, diverging
    bond and bond angle contributions throw them out during
    post-processing, e.g. by using the SKIP BOND ANGL UREY command before
    starting block post-processing.  If you want to see what can go wrong,
    look at the block2 test-case...

Dual Topology Soft Core Potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The new commands PSSP/NOPSsp and the optional parameters ALAM and
DLAM control the interactions between soft core potentials and BLOCK,
which is essentially the same as the PSSP command in the PERT soft
core (see pert.doc). After you specify PSSP inside BLOCK, soft core
LJ and electrostatic interactions will be used inside block interactions.
For the atom based NBOND command (NBOND ATOM), the block coefficents
(lambda) of VDW and ELEC can be defined as the different values. For
the group based case (NBOND GROUP), they share the same lambda value
currently. The separation parameters for elec. and LJ interactions can
be set with the ALAM and DLAM options, the default of 5A^2 should be
reasonable. The option is memorized, i.e., after the first invocation
of PSSP, all further calls of EVDW will use soft core interactions.
To turn this off, please use the NOPSsp keyword inside BLOCK/END pair.
So far, FAST OFF is recommended." -- New by H. Li and W. Yang

Adaptive Integration (ADIN) Method for Hybrid MD/MC Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to overcome the trapped distribution at certain lambda value in the
chemical space hybrid MD/MC simulation, adaptive integration method was
implemented. In this method, the biasing free energy potential is derived
by linearly integrating the ensemble average of energy derivatives at various
lambda values. By adaptive integration method, free energy difference between
two end states can be quickly computed. It is noted that this technique works
well when free energy has linear relationship with lambda value. It can crash
when there is severe end point singularity problem. Its general efficiency is
lower than the simulated scaling method, which does not suffer from end point
singularity problems. - by Lianqing Zheng and Wei Yang

Theta-dynamics
^^^^^^^^^^^^^^

This is an alternative method for the original lambda-dynamics. Lambda**2 is
replaced by sin(theta)**2 and (1-lambda**2) by cos(theta)**2. Theta, instead
of lambda, now is the variable for propagation. This implementation can avoid
the artifacts brought in by the constant external works in the Lagarangian
Multiplier boundary treatment. In the theta-dynamics, history dependent
approaches can work very nicely with no danger of being trapped at the end
points. - by Lianqing Zheng and Wei Yang

Multi-Site lambda-dynamics (MSLD)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is a more generalized lambda-dynamics method that allows
multiple substituents on multiple Sites on a common framework to be
evaluated simultaneously. Different functional forms of lambda have been
implemented which inherently satisfy the constraints that each lambda
should vary between 0 and 1 and the sum of the lambda values at a given
Site must equal 1. This strategy reduces the need to use Lagrangian
Multipliers and renormalization schemes and, for most systems, the timestep
can be increased in dynamics to 2 fs when SHAKE is invoked.
- by Jennifer L. Knight and Charles L. Brooks III

.. _block_examples:

Examples
--------

Here is an example of independently scaling the attractive
and repulsive terms in the Lennard-Jones interaction:

::

  ! scale the interaction parameters
  block 2
  call 2 sele segid heli end
  coeff 1  1 0.0 ! turn off the interactions between atoms in set 1
  coeff 1  2 1.0 vdwa 0 vdwr 1.0 ! scaling ratio to scale interactions
                                 ! between protein and other atoms
  coeff 2  2 1.0 ! leave interactions within the protein unchanged
  end

In this example we turn off the attractive term (vdwa) in the LJ interaction
and have only hard-core repulsion.

