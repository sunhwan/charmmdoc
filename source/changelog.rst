
=============================
CHARMM Developer's Change Log
=============================

Entries in each node are recorded by CHARMM developers to indicate new
and modified features of CHARMM during the development cycle, i.e., the 
alpha version period.

::

   ------------------------------------------------------
    CHARMM22.0.b  Release           April     22, 1991
    CHARMM22.0.b1 Release           September 30, 1991
    CHARMM22      Release           January    1, 1992
         c22g1    Release           February  15, 1992
         c22g2    Release           July       7, 1992
         c22g3    Release           November   3, 1992
         c22g4    Release           March      1, 1993
         c22g5    Release           August     1, 1993

    CHARMM23.0
         c23a1    Developmental     August    15, 1992
         c23a2    Developmental     October   25, 1992
         c23f     Developmental     March      1, 1993
         c23f1    Developmental     March     15, 1993
         c23f2    Developmental     August    15, 1993
         c23f3    Release           February   1, 1994
         c23f4    Release           August    15, 1994
         c23f5    Release           March     15, 1995

    CHARMM24.0
         c24a1    Developmental     February  15, 1994
         c24x1    Evaluation        February  15, 1994
         c24a2    Developmental     August    15, 1994
         c24a3    Developmental     March     15, 1995
         c24b1    Release           August    15, 1995
         c24b2    Release           February  15, 1996
         c24g1    Release           August    15, 1996
         c24g2    Release           February  15, 1997
 
    CHARMM 25.0
         c25a0    Developmental     August    15, 1995
         c25a1    Developmental     February  15, 1996
         c25a2    Developmental     August    15, 1996
         c25a3    Developmental     February  15, 1997
         c25b1    Release           August    15, 1997
   ------------------------------------------------------

 
.. _changelog_c21-c22:


Summary of Modifications of Developmental CHARMM21 to CHARMM22
--------------------------------------------------------------

* Linear pressure ramping added to CPT code (see :doc:`pressure`)
* Frequency based crystal update is now supported
* Relevent new keyword is IXTFrq (see :doc:`image`)
* Constant Pressure and Temperature (CPT) dynamics (See :doc:`pressure`)

  * TRICLINIC unit cell is now supported.
  
* Miscellaneous commands:

  * UPPEr and LOWEr keywords added (see :doc:`miscom`)
  
* Minimization: new keyword (FMEM) for ABNER minimizer (see :doc:`minimiz`)
* Internal coordinates (see :doc:`intcor`)

  * New commands:
    
    * IC SAVE
    * IC RESTore
    * IC RANDom [iseed]

  * Internal coordinates converted to double precision.
  
* Coordinate Manipulation (See :doc:`miscom` and CORMAN.DOC)

  * New inline command varaibles added:
    :sub:`THETa`, :sub:`XMOVe`, :sub:`YMOVe`, :sub:`ZMOVe`, :sub:`RMS`
    
  * New CORMAN commands added:
    
    * COOR HELIx
    * COOR PUCKer
    * COOR COVAraince
    * COOR SEARch ... RBUFF ...
    
* Energy, Angles
  
  * Urey-Bradley 1-3 terms have been added as an option.
  * Format of parameter file affected.  (See :doc:`io`)
  * Energy analysis code added (ANALysis ON command). (See :doc:`analys`)
  
* NOE distance restraints (See :doc:`cons`)
  
  * Overhaulled to become a general distance restraint term.
  * Commands syntax overhaulled as well.

* PSF common structure modified
  
  * Unused PSF arrays removed.  All size limits increased.
  * Binary file format changed to INTEGER*4 and REAL*8
  * PSF numbers added to ?variable list (See :doc:`miscom`).

* Output redirecting implemented. (See :doc:`miscom`)
  
  * OUTU replaces all writes to unit 6.
  
* ATLIM modified to allow a limit of several days.
  
  * PASMID has been changed to an integer which points the
  * current day.  See :doc:`miscom`
  
* Free energy perturbation commands added. (See :doc:`pert`)
  
  * Several new commands and features have been modified
    to allow free energy perturbation simulations to be performed.

* Partition function and classical free energy codee added to the vibrational
  analysis code. (See :doc:`vibran`)
  
  * Atom selection added for EDIT commands.
  * Atom selection added for WRITE SECOnd-derivatives CARD command.

* New time series commands and options (See :doc:`correl`)

  ::
  
      ENTER PUCKer
      ENTER HELIx
      ENTER RMS
      ENTER ENERgy
      ENTER RMS [MASS] atom-selection
      ENTER ATOM CROSsproduct
      ENTER FLUC CROSsproduct
      ENTER VECT CROSsproduct
      ENTER HBOND
      ENTER MODE
      ENTER RMS [MASS] [ORIEnt]
            ...
      TRAJ ... atom-selection
      MANTIME SQUARE (vectors now allowed)
      MANTIME ABS    (vectors now allowed)
      MANTIME ACOS
      
* Off-by-one error removed in time series data (time series now do not start
  at time zero, but at time DELTA*SKIP).

* Langevin dynamics modified.
  
  * An improved algorithm has been incorporated which gives a more accurate
    integration at low gamma values as well as the proper brownian dynamics
    limiting values in the large gamma limit (and is more efficient).
  * The gaussian random generator has been replaced to give a much more
    accurate distribution and uses only one random number call per atom
    by using an error function lookup table.

* Miscellaneous commands added. (See :doc:`miscom`)
  
  * DIVIde, EXONent, RANDom, and SHOW
  
* New miscellaneous variables added.
  
  * :sub:`RAND`

* Precision and index limits improved.
  
  * The entire program (except for the graphics section) has been
    converted to REAL*8 and INTEGER*4 from REAL*4 and INTEGER*2.

* Constant Pressure and Temperature (CPT) dynamics added. (See :doc:`pressure`)
  
  * Pressure analysis code added.
  * NTRFRQ usage modified so that it works for IMAGES and CRYSTAL.

* Heuristic nonbond update feature added. (See :doc:`nbonds`)
* New (consistent) energy print format with search line indicators.
* Graphics subsection added for workstations.
* New GRADient option added for most minimization methods for
  searching for saddle points.
* FAST option is now the default.  It is no longer necessary to have the
  command "FAST 1" in order to use the efficient energy routines.
* Constrained reference now only set for selected atoms for the CONS HARMonic
  command (the old method limited versatility). (See :doc:`cons`)
* Parallelization for shared memory multi-processor machines has been 
  implemented. Functionality for the fast energy routines has been increased.
  The vector/parallel routines will now to no electrostatics and novdw
  as well as simple cut-offs.
* SPECIfy  command. Controls various options such as I/O buffer flushing
  maximum number of processors to be used and whether to use the fast
  nonbond list generator.
* ``SYSTem "unix bourne shell commands"`` This command permits the user to issue
  Unix shell commands from the program. The command string must be enclosed
  in double quotes to prevent the CHARMm parser from converting the string
  to uppercase.
* SHAKE FAST This command specifies the use of the new vector/parallel SHAKE
* Deleted Features:
  
  * The old VAX analysis facility has been removed.
  * Sigma van der Waal switching and shifting options has been removed.
  * BARRI command removed.

.. _changelog_c20-c22:

Major Enhancements and Developments in CHARMM22
-----------------------------------------------

As CHARMM20 is not clearly defined, it is not straightforward to sort
out major differences between the current version of CHARMM
(CHARMM22.0) and a previous version (CHARMM20 or CHARMm21).
The VAX version CHARMM on HUCHE1 turns out to be a "developmental"
version towards CHARMM21 and contains the crystal facility, BLOCK, etc.
The following is prepared by comparing the developmental VAX version
CHARMM21 source code and that of CHARMM22.0.

Obsolete Modules Deleted from CHARMM20
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

(1) GRAMPS
    It is supported only in the VAX version CHARMM20.
    TH:[MK.PROT.SOURCE.VAX]GRAMPS.FLX contains an interactive routine that
    writes several files for the command language interpreter for
    producing computer graphics on the Evans & Sutherland
    Multi-Picture-System called GRAMPS.  This obsolete feature is no
    longer supported in CHARMM22.

(2) PARAmeter Optimization
    PARMOP is not incorporated in the VAX version CHARMM20 either except
    at the point of command parsing.  It seems that the feature has never
    been included in the central version.

New Features in CHARMM22
^^^^^^^^^^^^^^^^^^^^^^^^

(1) BLOCK

    The developmental CHARMM21 VAX version supports some BLOCK commands.
    The BLOCK commands are used to partition the molecular system into
    blocks and allows for the use of coefficients that scale the
    interaction energies between the blocks.  Specific commands to carry
    out free energy simulations with a component analysis scheme have been
    implemented.

(2) CRYStal

    The CRYStal commands are used to build a crystal with any space group
    symmetry, to optimize its lattice parameters and molecular coordinates
    and to carry out a vibrational analysis.  The CRYSTAL program is
    incorporated into the IMAGE module.  The VAX developmental version has
    a separate CRYSTL module.

(3) COOR COVAri

    The new COORdinate subcommand COVAriance is added.  It computes
    covariances of the spatial atom displacements of a dynamics trajectory
    for selected pairs of atoms.

(4) CORR HELIx / CORR PUCKer

    The New CORRelation commands HELIx and PUCKer are introduced.  The
    HELIx command computes time series of the helical axis orientation and
    PUCKer computes that of the sugar pucker phase and amplitude.

(5) DRAW, GRAP

    The new module GRAPHICS provides CHARMM the capability of displaying
    molecular structures when run on a graphics workstation.  (Currently
    works only on Apollo machines.)

(6) HBTRim

    The HBTRim command deletes hydrogen bonds that have an energy of
    interaction that is higher than the specified cutoff.  This command is
    used to reduce a list of hydrogen bonds to that of important hydrogen
    bonds.

(7) MOLVIB

    MOLVIB is a general purpose vibrational analysis program, suitable for
    small to medium sized molecules (less than 50 atoms).  It performs
    canonic force field calculations (KANO), crystal normal mode analysis
    for k=0 (CRYS) and other vibrational analyses in internal coordinates
    or in Cartesian coordinates.  Details are documented in :doc:`molvib`.

(8) PERT

    The PERTurbe command allows the scaling between PSFs for use in energy
    analysis, comparisons, slow growth free energy simulations, and
    widowing free energy simulations.  This is a rather flexible
    implementation of free energy perturbation that allows connectivity to
    change.  Also, three energy restraint terms (harmonic, dihedral and
    NOE) are subject to change which allows a flexible way in which to
    compute free energy differences between different conformations.

(9) QUANTUM

    Quantum mechanical and molecular mechanical combined force field
    method is implemented by employing the semi-empirical SCF method of
    the MOPAC program.  This module has not been tested nor documented.
    The code does not confirm CHARMM coding standards.  The future of the
    code is not certain at the time of the current release.

(10) RMSD

     The new RMSDyn routine is a modified CORMAN routine by William D.
     Laidig, which computes the RMS difference between two trajectory files
     and make a matrix of results.  
     
(11) RXNCOR

     The RXNCor command is used for defining a reaction coordinate for any
     molecule based on its structure and impose an umbrella potential along
     that reaction coordinate  (i.e., to run activated dynamics along this
     coordinate) in order to trace out the free energy profile during the
     structural change along the coordinate.

(12) SOLANA

     The solvent analysis facility computes solvent averaged properties,
     e.g., the solvent velocity autocorrelation function, mean-square
     displacement function, solvent-solvent radial distribution functions,
     solvent-reference site radial distribution function, and the solvent -
     reference site deformable boundary force.

(13) TRAJ

     The new TRAJectory command is used to merges or to break up a dynamics
     coordinate or velocity trajectory into different numbers of units.

(14) TSM

     The Thermodynamics Simulation Method module performs the free energy
     simulation.

(15) Urey-Bradley Energy Term

     Urey-Bradley 1-3 terms have been added.  The developmental CHARMM21
     also includes U-B terms.

(16) Update

     Two new non-bonded neighbour list updating schemes are introduced; one
     has something to do with an automated updating procedure and the other
     with the list generation algorithm.
     
     When INBFRQ is set to -1 (which is the default), heuristic testing
     is performed every time ENERGY is called and a list update is done if
     necessary.
     
     A new routine NBNDGC (nbndgc.src), a modification of NBONDG, is
     introduced.  NBNDGC is based on a cubical grid searching algorithm and
     generates the nonbonded list in linear time, as opposed to quadratic.
     On the Convex C220, which is a vector machine, it is faster than
     NBONDG for any system larger than a few hundred atoms.

(17) Integrator

     The leap-frog integrator has been implemented.  While the "old" Verlet
     integrator is still available via the DYNA VERLet command (and is the
     default), the new integrator can be accessed by DYNA LEAP.  The velocity
     Verlet integrator is also added in CHARMM. This new velocity Verlet 
     integrator can be called by DYNA VVER.
     
(18) Constant Pressure & Temperature Dynamics (DYNCPT)

     The constant pressure/temperature dynamics algorithm is implemented
     following the paper by Berendsen et al. (J. Chem. Phys. (1984) 81(8)
     p.3684).


Modification of CHARMM20 to CHARMM22
------------------------------------

(1) ANALysis

    The VAX version analysis facility is replaced by an energy
    contribution array (ECONT).  All evaluated energy terms are
    partitioned into each atomic contribution and collected in the array,
    which is accessible through the SCALAR command.

(2) XRAY

    The XRAY command of CHARMM20 is replaced by the READ XRAY command in
    CHARMM22.  In CHARMM22, all I/O functions are parsed in mainio.src.
    The subroutine XRAY is changed to RDXRAY, which generates a card file
    compatible with Richard Feldmann's XRAY display program.

(3) NOE

    NOE constraint has been overhauled.  It now handles general distance
    restraint terms.

(4) MISCOM

    The miscellaneous command parser (miscom.src in CHARMM22) is modified.
    
    (1) The SKIPE command is parsed in MISCOM.
    (2) New command parameter (@x) handling commands are added: DIVIde,
        EXPOnentiate, GET, MULTiply and SHOW.
    (3) The RANDOM command is added to set random number specifications.
    (4) The STOP command is parsed in MISCOM.
    (5) The QUICk (or Q) command is added to carry out a quick coordinate
        analysis.

(5) HANDLE

    The subroutine HANDLE is improved to accept command line arguments
    given with the CHARMM command issued to an operating system.  It works
    on most UNIX, UNICOS and VAX/VMS versions.

(6) Command Parameters

    In CHARMM20, we have ten command parameters @n, where n is a single
    digit, 0 through 9.  It is expanded to support any single
    alpha-numeric character so that one can use upto 36 command
    parameters (0-9, a-z).

(7) Dynamic Memory Allocation

    Most of UNIX versions now support VEHEAP.  VEHEAP was originally
    implemented by employing VAX/VMS system calls.  It expands the HEAP
    common block when more HEAP space is needed.  In UNIX versions, we use
    the UNIX system library routine malloc(), if available (the
    availability depends on the machine), to perform the same function.  

(8) File Format / Compatibility

    All binary files except dynamics trajectory are written in double
    precision format and not compatible with old versions.  For PSF,
    topology, parameter, etc. one should use CARD format to transfer
    previous version files to CHARMM22.  Trajectory files are written in
    single precision and compatible with all CHARMM versions and QUANTA.
    Old version dynamics restart files are not compatible with CHARMM22.

(9) Random Number Generator

    All random number routines are implemented in double precision (64-bit
    words).  Box-Muller algorithm is used for generating a Gaussian random
    deviat.  A machine specific random number routine (RANV of CONVEX
    VECLIB) is used in a CMU version.


.. _changelog_c22-c23:

Major Enhancements and Developments in CHARMM23
-----------------------------------------------

As an on-going project, CHARMM development has been carried out with
CHARMM version 23 series.  CHARMM development entails two objectives.
First, we maintain an integrated macromolecular science package
running on a wide range of computing devices.  Second, we incorporate
and exploit molecular simulation methodologies at the frontier of
current research.

In order to establish the first objective, we maintain all source
and support files under CVS (Concurrent Versions System) control.  The
ROOT repository is tammy.harvard.edu:/prog/chmgr/CVS.  CHARMM23 is
stored in /prog/chmgr/CVS/c23a.  A particular version is retrieved
with the version name as the rivision tag (e.g., c23f3).

Since we branched out from the CHARMM22 release version c22g2, we
have made two alpha versions  and four FORTRAN versions.

::

     c23a1    Developmental     August    15, 1992
     c23a2    Developmental     October   25, 1992
     c23f     Developmental     March      1, 1993
     c23f1    Developmental     March     15, 1993
     c23f2    Developmental     August    15, 1993
     c23f3    Release           February   1, 1994

c23f3 is the current release version.  As the "f" in c23f stands for
FORTRAN version, we converted FLECS source into FORTRAN.  The
conversion task had been completed as of c23f2.  Now CHARMM is written
in full FORTRAN except several machine dependant codes written in C.
The universal languages (C and FORTRAN) make it easier to port to new
machines in a broad range of architectural designs and to incorporate
new methodologies into a research version of CHARMM.

During the c23 development cycle, we have added and tested several
new features as described below.  We have also ported c23 to new
machines and supported c23f versions on the following platforms.


Platforms Supported
^^^^^^^^^^^^^^^^^^^

   ===========   ======================================
   RREFX key     Platforms
   ===========   ======================================
   ALLIANT       Alliant
   ALPHA         DEC alpha workstation
   APOLLO        HP-Apollo, both AEGIS and UNIX
   ARDENT        Stardent
   CONVEX        Convex Computer
   CRAY          Cray Research Inc.
   DEC           DEC ULTRIX
   HPUX          Hewlett-Packard series 700
   IBM           IBM-3090 running AIX
   IBMMVS        IBM's MVS platform
   IBMRS         IBM RS/6000
   IBMVM         IBM's VM platform
   IRIS          Silicon Graphics
   MACINTOSH     Apple Macintosh computers (system 7)
   SUN           Sun Microsystems
   VAX           Digital Equipment Corp. VAX VMS
   ===========   ======================================

New Features in CHARMM23
^^^^^^^^^^^^^^^^^^^^^^^^

(1) Cray Fast Code (Douglas J. Tobias)

    Vector/parallel code for energy calculation, shake, and nonbonded list
    generation on the Cray was implemented.  Dynamic heap and stack
    allocation on the Cray was added.

(2) PARALLEL (Bernard R. Brooks)

    General code for support of CHARMM on MIMD machines is completed.
    This includes control of the I/O levels for all file I/O.  For
    parallel machines or workstation clusters, only node zero performs I/O
    and it broadcasts are to other nodes.
    
    All compuationally intensive code exercised in MD is now fully
    parallel which includes: DYNAMC, ENERGY (and most subsections), SHAKE,
    PRSSRE, DYNLNG, IMAGES,...  Almost all comutationally intensive code
    in the first order minimizers is fully parallel.  Other usage of the
    energy routines are parallel (such as the energy time series in CORREL).

(3) Dynamics Integrator

    1. Leap-Frog Integrator (Bernard R. Brooks)
    
       Berendsen's method was modified so that it would work for very
       small systems and for very weak coupling constants.  Now it is
       possible to use SHAKE with CPT and get correct pressures and
       temperatures.  Another change is to calculate the change in potential
       energy due to the constant pressure algorithm.  The energy lost due to
       the changes in box size is now added to the kinetic energy during the
       constant temperature procedure.   This allows the constant presure
       code to nearly conserve energy and allows the constant temperature
       code to be used with weak coupling times.  This correction was made
       when we found that water box simulations with the Berendsen's method
       were running about 10 degrees too cold when both temperature and
       pressure coupling times of 1ps were used.  Now the correct target
       temperature is achieved, even in the limit of very weak couplings.

    2. EULER Dynamics Integrator (Bernard R. Brooks)
   
       The incorporation of of the Langevin/Implicit Euler dynamics
       integrator has been achieved.  The effect is to remove the energy in
       the high frequency degrees of freedom which eliminates the noise in
       free energy studies where bonds are being modified.  To support the
       Implicit Euler integration, a Truncated Newton Minimizer has been
       added.  This minimizer may be used directly using the MINI TN command.
       The minimizer is not yet fully implemented (it works, but is not as
       efficient as it will be), but it is already very competitive relative
       to existing minimization methods.  MINI TN does not work with SHAKE. 
       This code has been developed by Tamar Schlick at NYU.  It has been
       integrated within CHARMM with some modifications.

    3. EHFC: High Freequency Correction (Bernard R. Brooks)
    
       The leap-frog dynamics integrator has been modified to have an
       improved high frequency correction (HFC) term.  With the old term,
       energy was conserved within a harmonic degree of freedom, but total
       energy would drift as energy exchanged between high and low frequency
       degrees of freedom.  The new code avoids this problem.  The total
       energy and kinetic energy that is printed in the first line of
       dynamics energy printout has reverted to the standard Verlet energies,
       and these match the output of the old integrator.  The HFC terms
       (total energy, and kinetic energy) are now printed on the second line.
       The fluctuation of the HFC total energy is usually an order of
       magnitude smaller than that of the total energy.  The HCF total energy
       is a good indicator of problems with NVE dynamics because small
       changes in total energy are not lost in the noise of high frequency
       oscillations.

    4. Velocity Verlet Integrator  (Masa Watanabe)
    
       Velocity Verlet method has been implemented.  Two integrator
       (Verlet and Leap-frog) methods presented in CHARMM have their own
       flavors, but Verlet method handles velocities rather awkward and may
       introduce some numerical imprecision.  On the other hand, the
       Leap-frog integrator minimizes loss of precision on a computer, but it
       does not handle the velocities in a satisfactory manner.  Velocity
       Verlet integrator can store positions, velocities, and accelerations
       all at the same time and minimizes round-off error.

    5. Nose-Hoover Constant Temperature Method (Masa Watanabe)
    
       The constant temperature method has been implemented based on
       S. Nose, JCP 81, 511 (1984) and W.G. Hoover, Phy. Rev. A 31, 1695 (1985).
       This is an another type of constant temperature method, but an
       equilibration time in the vicinity of the desired temperature is
       faster than other routines which are available in CHARMM.  Also
       multi-temperature controls are also developed in order to equilibrate
       the system faster and keep the system in the desired temperature well.
       This method works with Verlet and Velocity Verlet integrators.

    6. Multiple Time-Scaled Method (Masa Watanabe)
    
       Tuckerman et al proposed a reversible RESPA algorithm recently
       (Tuckerman, Berne, Martyna, JCP 97, 1990 (1992)).  Previous MTS
       methods have the disadvantages of loosing accuracy due to the
       approximation of holding the slow variables fixed while integrating
       the equations for the fast variables.  But in this reversible RESPA
       equations of motions are derived from Liouville operators and Trotter
       theorem.  The method gives more accurate dynamics than previous
       methods.  In this implementation, one can specify up to three
       different time steps in dynamic simulation run.

(4) RISM (Reference Interaction Site Model) (Georgios Archontis)

    The RISM module allows the user to calculate the site-site radial
    distribution functions g(r) and pair correlation functions c(r) for a
    multi-component molecular liquid.  These functions can then be used to
    determine quantities such as the potential of mean force or the cavity
    interaction term between two solute molecules into a solvent, and the
    excess chemical potential of solvation of a solute into a solvent.  The
    change in the solvent g(r) upon solvation can be determined and this
    allows for the decomposition of the excess chemical potential into the
    energy and entropy of solvation.

(5) MMFP (Miscellaneous Mean Field Potential) (Benoit Roux)

    The MMFP Commands are primarily used for setting up special
    restraining potentials on some or all of the atoms.  The key word MMFP
    is used to enter the MMFP environement.  In the MMFP environment, all
    miscellaneous commands (label, goto, if, etc...), and string
    substitutions (with @1, @2, etc...) are supported.  The key word END
    returns to the main parser. The restraining potentials are used in all
    energy calculations, unless SKIP is used.  The subcommand RESET clears
    the potential.  This module is still under development and only the
    subcommand GEO is released.  The subcommand GEO (standing for
    geometrical) is used to setup various restraining potential
    (spherical, planar or cylindrical restraints) on some or all atoms.
    The selection specification should be at the end of the command.  The
    default atom selection includes all atoms.  Future subcommands will
    include continuum electrostatic reaction field and solvent mean field
    potentials. Expected date of release is Spring 1994.

(6) NMR Analysis (Benoit Roux)

    The NMR commands may be used to obtain a set of time series for a
    number of NMR properties from a trajectory.  Among the possible
    properties are relaxation rates due to dipole-dipole fluctuations (T1,
    T2, NOE, ROE), chemical shift anisotropy and Deuterium order
    parameters for oriented samples.

(7) REPLICA (Leo Caves)

    Tool to support LES and MCSS calculations.  Performs replication
    of arbitrary regions of PSF.  Data structure interfaces to non-bond
    list generation routines, to perform appropriate exclusions.  In
    association with BLOCK can provide appropriate energy/force
    normalizations for various classes of methods employing replicas.
    
    Introduced REPLICA and REPDEB preprocessor directives.  Code for
    cray multi-tasking list generation routine used inference and has not
    been tested.  Convex parallel code works fine.  Added miscellaneous
    parameters to report number of atom/group pairs from non-bonded
    routines: ?NNBA, ?NNBG, ?NNBI for atom/group/images respectively.  For
    replica-based exclusions from the list there are ?NRXA and ?NRXG for
    atom and group exclusions.

(8) Clustr code integrated into CORREL (Charles L. Brooks III)

    The CLUSTER command clusters time series data obtained within the
    CORREL facility.  The data are grouped into sets with similar time
    series values, using euclidean distance as the dissimilarity measure
    between different time frames of a set of time series.  It is useful,
    for example, for grouping together similar conformations or energy
    levels.

(9) GRAPHICS (Richard M. Venable)

    Graphics code converted to FORTRAN and overhauled.  Versions that
    work with Xwindows and GL are in progress.  A new preflx keyword,
    NODISPLAY, builds a version which produces HPGL, PLUTO FDAT, and
    LIGHT.atm files without requiring any screen display capabilities.
    The SG (IRIS) code incorporation is relatively untested.  Postscript
    file output similar to HPGL (but much nicer looking, hopefully) is
    also implemented.

Major Modifications
^^^^^^^^^^^^^^^^^^^

(1) Command Line Handling

    1. Extension of Command Line Parameter Handling (Leo Caves)
    
       A command line parameter token can now be a string rather than
       just one of the single characters 0-9 and A(a)-Z(z).  For substitution,
       a token is indicated by the use of the @ character as before.  The
       token is end-delimited by any non-alphanumeric character.  In the case
       that the token is not found in the parameter table, a check is made to
       see if the first character of the token is itself a token in the
       parameter table. If this single character token is in the table, the
       corresponding value is substituted -- this is the necessary scheme to
       allow backwards compatibilty with the old parameter substitution,
       which allowed parameters embedded in strings.  For unambiguous token
       detection, "protect" the token with brackets {} --- this allows for
       the use of non alphanumerics in tokens such as -, _.

    2. New Parsing Options (Bernard R. Brooks)
    
       The IF command will be expanded to allow commands such as:

       ::
       
            IF ?ENER .GT. ?VDW  THEN GOTO label
            
            or
            
            IF ?NSEL .LT. 8 THEN GOTO label

    3. MSCNUM (Bernard R. Brooks)

       New code for flexible miscellaneous command substitutions has been
       fully incoporated.  Additional types were needed to make this code more
       flexible.  Three types are supported, REAL(\*8), INTEGER, CHARACTER.
       There are three subroutines which can be called; integer (SETMSI),
       character (SETMSC), and real (SETMSR) to specify a command substitution
       variable.  Now it is possible for ?NATOM to return an integer, ?RSM to
       return a real number, and ?SEGID to return the segment identifier of the
       first selected atom.

(2) QUANTUM

    Quantum mechanical and molecular mechanical combined force field
    method was implemented by employing the semi-empirical SCF method of
    the MOPAC program in the CHARMM version 22.  The QUANTUM code has been
    modified extensively to meet CHARMM standards.
    
    There were several problems with the quantum code that have been
    fixed.  The van der Waal group nonbond list was missing due to an
    improper interpretation of the group-group exclusion list in CHARMM
    (It's a two state list, not a 3 state as in the atom-atom exclusion
    list).  All vdw interactions between QM and MM group where any QM atom
    had an exclusion or a 1-4 interaction with any MM atom were not
    computed.  This caused major problems in certain situations where
    there was a strong electrostatic attraction with no compensating vdw
    interaction.
    
    New code to add link and place link atoms has been written.

(3) Frequency Based Crystal Update (Ryszard Czerminski)

    The modification allowes for automated, frequency based, crystal
    update.  New variable (IXTFRQ) is introduced which controls frequency
    of the crystal update.

(4) Ability to Linearly Increase/Decrease Pressure (Ryszard Czerminski)

    The goal was to allow for linear increase (decrease) of the
    pressure during single dynamic run.  New variables/keywords were
    introduced (PIXX - initial value of XX component of pressure tensor,
    PFXX - final value etc... for other components).

(5) Atom Selection

    1. Atom Parse (Bernard R. Brooks)
    
       A new atom name parsing subroutine has been developed.  This makes
       the code simpler and facilitates further advancements in atom
       parsing.  One new feature allows an atom selection to be used to
       select a series of atoms.  This is very useful in CORREL for
       specifying clusters of atoms for analysis.  When the atom selection
       feature is used to specify 4 atoms of a dihedral, the first 4 selected
       atoms will be chosen.

    2. New Tokens (Bernard R. Brooks)
    
       * new operator; ``.BYGROUP. <factor>``
       * new token; ``IGROup  <int1> : <int2>``

       have been added to allow the selection of atoms based on electrostatic
       groupings.
       
       Several keynames have been added to allow the query of the
       characterstics of selected atoms;

       ::
       
          ?SELATOM  - number of first atom selected
          ?SELIRES  - number of first residue selected
          ?SELISEG  - number of first segment selected

          ?SELTYPE  - name of first atom selected
          ?SELRESI  - resid of first residue selected
          ?SELSEGI  - segid of first residue selected
          ?SELRESN  - residue type of first atom selected
          ?SELCHEM  - chemical type of first atom selected

       These new keywords are in addition to the existing keyword;
       
       ::
       
          ?NSEL    - Number of atoms selected

(6) Correlation

    1. New MANTim Options in CORREL (Bernard R. Brooks)
    
       A histogram option to time series manipulation has been developed.
       This is executed by the command;

       ::
       
         MANTime time-series-name HISTogram min-value max-value num-steps

       The selected time series is replaced with a histogram which contains
       the probability of finding the time series within a given value range.
       Also, new options (RATIo and KMULt) added to the CORREL MANTIME command.

    2. Dihedral Time Series in CORREL (Bernard R. Brooks)
    
       Fixed problems with the diheral code in correl to account for
       torsional timeseries.  The correct fluctuation is now determined.
       The extra processing has been removed from the SHOW command because
       the data may no longer be valid for this processing when MANTIME
       commands are present in a script.  A new command option "MANTime
       CONTinuous-dihedral" has been added to allow a dihedral timeseries to
       be unfolded to a continuous function. 

    3. Extension of Solanal ANALysis command (Arnaud Blondel)
    
       A command -CROSs- was added to allow a cross analysis on two
       selected subsets of atoms.  For the moment the exclusion of the couple
       of atoms belonging to the same SEGId is not implemented.  The keyword 
       CROSs cannot be selected with the following options: WATer, SITE,
       IKIRkg, ISDIst, IFDBf.  IVAC, IMSD and IFMIn have not been tested with
       CROSs.

(7) SCALAR Command Enhancement (Bernard R. Brooks)

    The ASP arrays (IGNOre, ASPV and VDWS) are now accessible.  There
    is a sort option for the SHOW command.  There is a new MASS keyword
    for the STATistics and AVERage commands
    
    A new SCALAR READ option has been added.  It allows values to be
    entered from a file.  The use is:

    ::
    
      OPEN READ CARD UNIT 12 NAME file.dat
      SCALar WMAIn READ 12 SELE ... END

    which will read selected entries to the weighting array.


(8) SURFACE (Bernard R. Brooks)

    New analytic surface area code and energy terms for ASP (Atomic
    Solvation Parameters) energy and forces have been fully integrated
    (and parallelized for multi-machines).  This has been achieved by the
    incorporation and adaptation of the code from Wesson and Eisenberg.
    The default for the COOR SURFace command is now the analytic surface
    area.  The anaylitic answer is less expensive and more accurate.  The
    older Lee and Richard's algorithm may still be invoked by specifying a
    nonzero RPRObe value.  The maximum number of contacts that a sphere
    may have has been increased from 15 to 35.


(9) QAUGMENT (Bernard R. Brooks)

    It is desirable for a patch to be able to augment the charge of an
    atom.  The current code could only set a charge.  The new code can add
    or subtract a value from the charge.  This is done by using a patch
    charge value near 100.0.   For example, a charge of 100.15 will add
    0.15 to the current charge. A charge value of -101.0 will subtract 1.0
    from the current charge.  Charge values less than -90.0 or larger than
    90.0 are no longer allowed for generate or patch without charge
    augment.  It allows more flexible patches to be developed where the
    prior charge on modified atoms need not be known.

(10) COORdinate Commands

     1. VACUUM_OP: COOR SEARCH Subcommand (Bernard R. Brooks)
     
        The ability to manipulate pixel bitmaps generated from the COOR SEARCH
        command has been developed. The new syntax for the COOR SEARCH command is;

        ::
        
            COOR SEARch {PRINt [UNIT int]} {            } {[VACUum]} {[RESEt]} [SAVE]
                        {[NOPRint]       } {[RCUT  real]} { FILLed } { AND   }
                                           {[RBUFf real]} { HOLES  } { OR    }
                                                                     { XOR   }

        The new keywords are;
        
           ===== ===============================================================        
           SAVE  save the resultant bitmap for subsequent operations
           AND   logical AND the new bitmap with the previously saved map
           OR    logical OR  the new bitmap with the previously saved map
           XOR   logical XOR the new bitmap with the previously saved map
           HOLES search for holes (vacuum points surrounded by filled points)
           ===== ===============================================================

     2. New COOR DIST command (Bernard R. Brooks)
     
        The COOR DISTance command has been overhauled and has additional
        features.  One such feature is the ability to get g(r) plots from
        trajectory files using atom selections.  It has several other
        features.  The new syntax is:

        ::
        
            COOR DISTance

                {  WEIGhting vector-spec               atom-selection           }
                {                                                               }
                { [UNIT int] [CUT real] [ENERGy [CLOSe]] 2X(atom-selection) -   }

                        { [Nonbonds] } { [NO14exclusions] } { [NOEXclusions] }  -
                        { NONOnbonds } {    14EXclusions  } {    EXCLusions  }

                     [TRIAngle]   [ HISTogram HMIN real HMAX real HNUM integer  -
                                     [HSAVe] [HPRInt] [HNORm real] [HDENsity real] ]


(11) JOIN/RENUMBER Command (Bernard R. Brooks)

     A "JOIN segid RENUMBER" feature is added in the JOIN command.
     This allows resid's to be made sequential within a single segment.

(12) PREFX.SRC overhauled. (Bernard R. Brooks)

     The PREFX program has been overhauled.  The new code has the
     following features: 

     - It allows "!" comments at the end of valid FORTRAN statements.
     - Conversion to single precision is performed ONLY if the SINGLE
       keyword is present.
     - It allows the use of identifier comments in ## statements.
       For example:
       
       ::
       
          ##IF PERT (pertprint)
          ...
          ##ELSE (pertprint)
          ...
          ##ENDIF (pertprint)

     This makes the code easier to read and allows ##ENDIF statements to be
     uniquely identified.  A fatal error is flagged if the identifiers do
     not match.

.. _changelog_c23-c24:


Major Enhancements and Developments in CHARMM24
-----------------------------------------------

During the C24 development cycle, February 15, 1994 to February 15, 1996,
we made two bugfix-updates in the c23 releases and three alpha versions
and one beta version in the c24 development line.  c24x1 is the MMFF
implementation in CHARMM developed at the Molecular Simulations Inc.

::

        CHARMM23.0
             c23f4    Release           August    15, 1994
             c23f5    Release           March     15, 1995

        CHARMM24.0
             c24a1    Developmental     February  15, 1994
             c24x1    Evaluation        February  15, 1994
             c24a2    Developmental     August    15, 1994
             c24a3    Developmental     March     15, 1995
             c24b1    Release           August    15, 1995

Only bugfixes are incorporated into CHARMM23 and all new developments
and enhancements have been carried out with the CHARMM24 developmental
versions.  All modifications are thoroughly recorded in the
ChangeLog.c24 file and the following is the summary of new features
and major enhancements in CHARMM 24.

New Features in CHARMM24
^^^^^^^^^^^^^^^^^^^^^^^^

(1) New Ports and Parallel Versions

    1. Enhancement to Parallel Code (Bernard R. Brooks and Milan Hodoscek)

       There has been continued development of the parallel code for
       CHARMM.  This includes new features run in parallel, new machine types
       supported, new parallelization methods, and code made to run more
       efficiently.  Due to conflict in routine names with library routines,
       the subroutines: WRITEC and READC had to be renamed.
       
       Initial code to allow the use of the Terra parallel computer has
       been added.  Added preflx keyword SGIMP for multiprocessor SG machines
       using PVM massage passing library.  The difference between PVM and
       (SGIMP, PVM) is that all the processes are spawned on one host and
       some communication parameters are not supported on MP machines. It can
       be used on a single processor SG for testing purpose. Use PVM only on
       a cluster of any type of workstation. 

    2. Convex Exemplar SPP-100 and generic PVM Ports (Charles L. Brooks, III and Stephen H. Fleischman)
    
       A port of CHARMM version 24a2 to general PVM based parallelism
       using existing parallel code as well as a port to the Convex parallel
       machine are included.


    3. Cray T3D Port (Charles L. Brooks, III and Barry C. Bolding)

       A port of CHARMM version 24a2 to the Cray T3D parallel computer using
       existing parallel code is included.

    4. Port of parallel CHARMM to Convex Exemplar SPP-1000 and generic MPI (Charles L. Brooks, III and Stephen H. Fleischman)

       A port of CHARMM version 24a3 to general MPI based parallelism
       using existing parallel code as well as a port to the Convex parallel
       machine are included.


    5. Thinking Machine's CM5 Port (Robert Nagle)

       Previous communication scheme was based on a simple send and
       receive model.  By using TMC's active message layer, communication
       bandwith can be increased by anywhere from 50% to 5X.

    6. OS/2 Port (Stefan Boresch)

       CHARMM (c23f4 and c24a3) has been ported to the OS/2 operating
       system, version 2.x and higher.  The Watcom Fortran compiler (v. 9.5,
       patch-level (c)) has been used.  A new pre-processor keyword, OS2, has
       been introduced, and all OS/2 related changes hide behind the OS2
       keyword.  There is currently no install script.  Please contact me
       if you want to build an OS/2 version of CHARMM (boresch@tammy.harvard.edu).

(2) Fast Multipole Code for Electrostatic interactions (Robert Nagle)

    This is an initial implementation of a fast multipole method,
    based on John Board's work.  A new non-bond option (FMA) has been added.
    This replaces cut-off parameters with a no cut-off hierarchical
    technique.  The advantages of this method are that you can control the
    error and that it is amenable to parallelization.  FMA is an O(N)
    technique but the constant is large and so FMA will, in general, be
    slower for systems of less that 5000 atoms, for the same accuracy.
    
    Two options, LEVEL and TERMS, govern how many hierarchical levels
    are used and how many terms are retained in the expansion, respectively.
    In the method, each box at every level is subdivided into 8 sub-boxes
    - you should select LEVEL so that the boxes at the lowest (i.e.
    finest) level contain 10-20 atoms on average: 3 or 4 will be typical
    choices.  You then select TERMS to control the accuracy that you
    require: 4 will often suffice but I would generally recommend 6 or
    even 8.  See the references in :doc:`fma` for a detailed description of
    the error bounds.
    
    NOFMA is the nonbond option which turns off the multipole method.
    Compilation of FMA is controlled by the flag, FMA, in pref.dat.
    
    FAST ON is required for this initial implementation.  This
    implementation is not yet parallelized.

(3) Energy Embedding by the Addition of a Higher Spatial Dimension (Elan Z. Eisenmesser / Carol Post)

    The energy embedding technique entails placing a molecule into a
    higher spatial dimension [Crippen, G. M. & Havel, T. F. (1990) J.
    Chem. Inf. Comput. Sci. Vol 30, 222-227].  The possibility of
    surmounting energy barriers with these added degrees of freedom may
    lead to lower energy minima.
    
    With the recent success of using four dimensions in the GROMOS
    force field [Van Schaik, R. C., Berendsen, H. J. C., Torda, A. E., &
    van Gunsteren, W. F. (1993) J. Mol. Biol. Vol 234, 751-762], creating
    a similar option in CHARMM should also prove advantageous.
    Specifically, another cartesian coordinate was added to the usual X,
    Y, and Z coordinates and was appropriately named FDIM for Fourth
    DIMension.  This implementation has led to alterations in some
    existing code along with the addition of several algorithms.

(4) DIMB (Diagonalization In a Mixed Basis) Method (David Perahia, Liliane Mouawad, Herman van Vlijmen)

    The DIMB (Diagonalization In a Mixed Basis) method (see L. Mouawad
    and D. Perahia (1993), Biopolymers, 33, 599) is an iterative method to
    calculate the N lowest normal modes of molecules.  It is especially
    targeted to do large molecules, since it does not require the full
    Hessian to be stored in memory or on disk.  In short, the method
    does repetitive reduced-basis diagonalizations in bases that consist
    partially of the approximate eigenvectors, and partially of Cartesian
    coordinates.  Eigenvectors are saved to file during the process.  Before
    that is done, a new basis is again created, which consists of the
    approximate eigenvectors at that point + the residual vectors (Lanczos
    vectors).  This accelerates the convergence.  A very good property of
    this method is that the final eigenvectors are as accurate as the user
    wants them to be, so the results are no different from a full-blown
    diagonalization.
    
    Because the method is iterative, it takes longer to converge than
    a regular diagonalization.  Sizewise it can handle almost anything on
    a moderately sized computer.  David Perahia calculated a few dozen modes
    of Hemoglobin (~600 residues = ~6000 atoms = ~18000 d.o.f.) on a
    SGI workstation with 90 Mb memory.  I have done several calculations
    on 900 residue systems.  The actual time to reach convergence depends
    on the available memory, the desired accuracy, and the number of
    requested normal modes.
    
    One other area where the method saves memory is in the storage of the
    original Hessian.  Since this matrix is usually sparse for large systems,
    a compressed Hessian is set up, which contains all non-zero elements.
    
    In addition, I added the option to used this compressed Hessian in the
    reduced-basis diagonalization option of VIBRAN.  Before, the same size
    limits applied to full diagonalizations and reduced-basis diagonalizations.
    This should not be: people usually want to do reduced-basis calculations
    because the molecule is too big for the Hessian to be stored in memory.
    The option VIBRAn REDUce CMPAct will fill the compact Hessian and 
    form the reduced-basis Hessian from this compact Hessian.  Overall, this
    is a big saving on memory space.

(5) Arithmetic Expression Interpreter (Benoit Roux)

    An interpretor of arithmetic expression has been added to the
    CHARMM command parser.  It is called at the level of the miscellaneous
    command handling using simply by the word CALC (for calculator).
    It can be used to evaluate algebraic numerical expression.  The command
    supports all mathematical numerical expression with arbitrary number
    of nesting of recursive parentheses, e.g.,

    ::
    
       exp[1.0-cos(2*(log(2*pi))**2)/0.5]

    The parsing is actually very crude since the expression is translated
    back and forth between character string and a real variable to handle
    the logic (there is no real subroutine recursion).


(6) TNPACK Update (Tamar Schlick, Phillipe Derreumaux and Eric Barth)

    The truncated-Newton minimization package TNPACK, developed by
    T. Schlick and A. Fogelson, has been incorporated into CHARMM and
    adopted for biomolecular energy minimization.  TNPACK is based on the
    preconditioned linear conjugate-gradient technique for solving the
    Newton equations.  The structure of the problem --- sparsity of the
    Hessian --- is exploited for preconditioning.
    
    Thorough experience with the new version of TNPACK in CHARMM has
    been described in a paper now in press in the Journal of Computational
    Chemistry: Applications are reported for a series of molecular systems
    including Alanine Dipeptide (N-Methyl-Alanyl-Acetamide), a dimer of
    N-Methyl-Acetamide, Deca-Alanine, Mellitin (26 residues), Avian
    Pancreatic Polypeptide (36 residues), Rubredoxin (52 residues), Bovine
    Pancreatic Trypsin Inhibitor (58 residues), a dimer of Insulin (99
    residues), and Lysozyme (130 residues).  Through comparisons among the
    minimization algorithms available in CHARMM, we find that TNPACK
    performs significantly better than ABNR in terms of CPU time when
    curvature information is calculated by a finite-difference of
    gradients (the "numeric" option of TNPACK).  The CPU gain is 50% or
    more (speedup factors of 1.5 to 2.5) for the largest molecular systems
    tested and even greater for smaller systems (CPU factors of 1 to 4 for
    small systems and 1 to 5 for medium systems).  With the analytic
    option, TNPACK converges more rapidly than ABNR for small and medium
    systems (up to 400 atoms) as well as large molecules that have
    reasonably good starting conformations; for large systems that are
    poorly relaxed (i.e., the initial Brookhaven Protein Data Bank
    structures are poor approximations to the minimum), TNPACK performs
    similarly to ABNR.
    
    TNPACK uses curvature information to escape from undesired
    configurational regions and to ensure the identification of true local
    minima.  It converges rapidly once a convex region is reached and
    achieves very low final gradient norms, such as of order 10E-8, with
    little additional work.  Even greater overall CPU gains are expected
    for large-scale minimization problems by making the architectures of
    CHARMM and TNPACK more compatible with respect to the
    second-derivative calculations.
    
    This work should be the focus of future developments.  Such work
    involves sparse storage of the Hessian, efficient sparse
    Hessian/vector multiplications, and separation of the gradient and
    Hessian calculations.

(7) X-window graphics extensively modified (Richard M. Venable)

    Several new features have been added to the X-window version of CHARMM
    graphics.  This code has also been tested on a wider variety of
    hardware platforms (for example: SGI).
    Changes include: double-buffering, clipping, StaticColor, symbol fonts,
    window title, modified colormap calls, and a misc.  Bug fixes in the
    labeling of the X axis.  A NODISPLAY compile option has been added to
    the X windows version of CHARMM graphics in which only derivative
    files are produced.  The GRAPhics NOWIndow option can be used to
    generate the same effect at run time.


(8) Minimum Image Periodic Boundary Code (Charles L. Brooks, III, William A. Shirley and Stephen H. Fleischman)

    Simple minimum periodic boundary conditions are added for cubic,
    truncated octahedra and rhomboidal (dodecahedra) periodicities which
    augments the image facility and enhances parallel scaling on scalar
    parallel machines as well as significantly reducing the memory
    requirements.  This code is developed and fully tested for the
    simulation cells described above when the cell edgelength is the same
    in all dimensions.  The (trivial) extension to non-identical cell
    sides will be added.  However, it is critical to see reasonable
    performance on all scalar parallel platforms where simulations using
    images are currently employed that this enhancement be added now.

(9) GAMESS Code (Bernard R. Brooks and Milan Hodoscek)

    The CHARMM-GAMMES interface is under development.  The interface
    part is completed and testing is in progress.

Major Enhancements in CHARMM24
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

(1) New Dihedral / Improper Dihedral Energy Routines (Arnaud Blondel)

    The previous energy routines used the derivatives d(cos(phi))/dr
    to calculate the forces and the second derivatives.  This choice
    introduced an artificial singularity at sin(phi)=0.
    
    The new routines use the derivative d(phi)/dr and thus have no
    singularities.  This removes the tests to avoid numerical overflow or
    the switch functions in the vector improper routines.
    
    The new dihedral routines now support cases where planar conformation
    is not an extremum.  Thus a value other than 0 or 180 can be specified
    in the dihedral parameters.  The dihedral constraints can also use the
    dihedral functional form using the key word PERIod and giving a
    non-zero number.

(2) Extended Pressure System, Langevin Piston Code (Bernard R. Brooks, Scott E. Feller and Yuhong Zhang)

    The constant pressure code has been overhauled.  The old method
    based on Berendsen's method has been replaced with a Langevin Piston
    Method.  When no friction is applied, this method becomes the standard
    method based on Nose and Klein (adapted from Andersen).  At the limit
    of infinite friction with no random force, this reverts to the
    Berendsen method.
    
    The unit cell information has been added to the trajectory file
    format.  This implementation required an update to the image and
    crystal code which cleaned up some ancient problems.  Options for
    including the surface tension (gamma-Area) term is also completed and
    tested.  This has been developed for the accurate simulation of
    interfacial systems.

(3) Anisotropic Harmonic Restraints (Bernard R. Brooks)

    The global scale factors: "XSCAle", "YSCAle", and "ZSCAle" have
    been added to the "CONS HARM" command.  This allows using the CONS
    HARM to enforce a planar or linear restraint.  This feature is also
    useful for use in conjunction with our COORPLAS program (for generating
    3-D coordinates from plastic models).

(4) New RESDistance Facility (Bernard R. Brooks)

    A new facility, RESD, has been created to allow general distance
    restraints based on a linear combination of distances.  This is useful
    for searching reaction pathways.

(5) New READ PARAm APPEnd Option (Bernard R. Brooks)

    An append option has been added to the READ PARAM CARD command.
    This allows just a few parameters to be modified without editing an
    entire parameter file.  A modification to the binary parameter file
    format was necessary.  Old binary files may not be appended, but they
    are still supported.

(6) New READ PSF APPEnd Option (Bernard R. Brooks)

    An append option has been added to the READ PSF command.  This allows PSFs
    to be easily merged to make a larger PSF.  No modification to the binary
    parameter file format was necessary.  This option works with both FILE
    and CARD options.

(7) Best Fit Option to CORREL TRAJectory Command (Bernard R. Brooks)

    The TRAJectory command in correl now accepts an ORIENt keyword with an
    optional [MASS] qualifier in conjunction with a second atom selection
    that will best fit selected atoms with respect to the rms deviation
    from the reference structure (in the comparison coordinate set).  This
    operation is done prior to the determination of any time series value.
    This operation will not affect any time series value that is based
    only on relative distances and angles.


(8) QM/MM Exclude Group Option (Bernard R. Brooks)

    An option EXGRoup has been added which causes all atoms in the group
    of the link atom host to be excluded from the QM/MM electrostatic
    interaction terms.  Code for specifying the charge of link atoms and
    their placement has also been added.

(9) Enhancements to the Ewald Code (Bernard R. Brooks,  Scott E. Feller and Steve Bogusz)

    The EWALD electrostatic option now runs efficiently for parallel
    architectures.  Also, the maximum K-space values can be specified
    independently for each direction.  Several bugs were fixed.
    Additional ways to compute ERFC() were added, including a lookup
    table.

(10) MMFP/SSBP Upgrade (Benoit Roux and Dmitrii Beglov)

     The Miscellaneous Mean-Field Potentials (MMFP) has been upgraded.
     The  spherical solvent boundary potential (SSBP) has also been
     incorporated into EPERT.  A new "membrane-like" planar potential
     has been introduced using Gaussians to provide a smooth free energy
     function based on hydropathy profile of individual amino acids
     and solvent exposure.  This is useful to orient membrane proteins.
     A new primary shell of hydration has been added to the MMFP facility
     to provide one layer of solvent around a flexible polypeptide.
     For more information, see Beglov & Roux, Biopolymers 35: 171-178 (1995).
     
     A solvent boundary potential for the simulation of water at
     constant pressure is also added to the Miscellaneous Mean Field
     Potential module.  The boundary potential is an approximation but
     follows from a rigorous statistical mechanical treatment of the
     boundary.  In light of the difficulties raised by the previous
     treatments, a different route was chosen to formulate and develop the
     solvent boundary potential for computer simulations of a finite
     representation of an infinite bulk system.  The present theoretical
     formulation is based on a separation of the multidimensional
     solute-solvent configurational integral in terms of n "inner" solvent
     molecules nearest to an arbitrary solute, and the remaining "outer"
     bulk solvent molecules.
    
     This formulation, which differs significantly from previous
     treatments, provides further insight into the statistical mechanical
     basis of the solvent boundary potential and is helpful in constructing
     useful approximations for computer simulations in dense liquids.
     An approximation to the solvent boundary potential is constructed for
     simulations of bulk water at constant pressure, including the
     influence of van der Waals (done with RISM) and electrostatic
     interactions (done with a Kirkwood-like multipole expansion).
     The approach has been tested with success on several typical systems
     (water, ions, n-butane and alanine dipeptide).

(11) Upgrade of the NMR module (Benoit Roux)

     The NMR module is upgraded to have better output style.  The old
     version used the value of PRNLEV to choose the printed quantities.
     Since this was a non-standard style in CHARMM, a series of logical
     flags have been included in the command calls to print some chosen
     quantities.  In addition, the chemical shift anisotropy (CSA, used in
     solid state NMR of membrane proteins in oriented samples) has been
     redefined in term of a zmatrix to prevent confusion.  The deuterium
     quadrupolar splittings (DQS) command is also upgraded.  A bug in a
     call to NORMAL was fixed.

(12) New Options to CORREL (Lennart Nilsson)

     Two new MANTime options have been added to CORREL: CROS and DOTP.
     CROSsprod name  Q(T) = Q(T) x Q2(T) produces the 3D crossproduct of
     the two 3D vectors formed by the selected and named timeseries and
     DOTProd name Q(T) = x-comp of Q(T)= Q(T) . Q2(T) gives x-comp of Q2(T)
     angle in degrees between the two vectors.

(13) The COOR HBONd Command (Lennart Nilsson)

     An option for the analysis of H-bond patterns from trajectories
     has been added to corman.

     ::
     
        COORdinates  HBONd 2X(atom-selection) [CUT <real>] [CUTA <real>] 
                 [IUNIt <int>]  [BRIDge <resnam>]
                 [FIRSt int] [NUNIts int] [NSKIp int] [BEGIn int] [STOP int]

     The HBONd command analyses a trajectory for hydrogen bonding
     patterns.  For each acceptor/donor in the first selection the average
     number and average lifetime of hydrogen bonds to any atom in the
     second selection is calculated.  A hydrogen bond is assumed to exist
     when two candidate atoms are closer than the value specified by CUT
     (default 2.4A, (reasonable criterion, DeLoof et al. (1992) JACS 114,
     4028), and if a value for CUTAngle is given the angle formed by D-H..A
     is greater than this CUTAngle (in degrees, 180 is a linear H-bond);
     the default is to allow all angles.  The current implementation
     assumes that hbonding hydrogens are present in the PSF and also uses
     ACCEptor and DONOr information from the PSF to determine what pairs
     are possible.
     
     If output is wanted to a separate file the IUNIt option can be
     used.  If the BRIDge option is used the routine calculates average
     number and lifetime of bridges formed between all pairs of atoms in
     the two selections; a bridge is counted a residue of the type
     specified with the BRIDge <resnam>  hydrogen bonds (using same
     criteria as for direct hbonding) to at least one atom in each
     selection.  The typical use of this would be to find water bridges.
     Here again, results are presented for each atom in the first selection.
     
     In order not to find hbonds between bonded atoms UPDATE is
     called, which requires coordinates to be present when invoking this
     module.  Since this is done just to get the non-bond exclusion lists,
     the cut-offs are set to very small values, and could influence
     subsequent energy evaluations if the non-bond cutoffs are not then
     respecified.

(14) NORESET Option for SHAKE (Lennart Nilsson)

     The NORESET option is added to allow multiple shake commands.
     It is useful to be able to define shake on bonds, bonh or so on
     several different sets of atoms, with different shake options.  The
     NORESET keyword to shake command allows this by not zeroing counter.

(15) Trajectory Reading (Lennart Nilsson)

     READCV is modified to read coordinates at multiples of skip FROM
     the actual first coordinate set in a trajectory file.

(16) Make BLOCK work with IMAGE/CRYSTAL and vice versa (Stefan Boresch)

     In order to make BLOCK work / coexist with the IMAGE module two
     things had to be changed: (1) A memory allocation problem in the BLOCK
     datastructure and (2) the post-processing modules needed to be
     overhauled to allow for nonbonded list updates while reading frames
     from the trajectory.
     
     Ad (1), memory allocation: BLOCK uses two data-structures, one
     containing the interaction matrix between blocks, and one containing
     the block number for each atom (IBLCKP).  This array was allocated so
     far as INTEG4(NATOM) on the heap.  However, when IMAGE atoms are
     present, the energy routines attempt to find out to which block an
     IMAGE atom belongs.  This at one point or the other causes a memory
     access violation.  The solution consists out of two parts.  (i) The
     IBLCKP data-structure is now allocated as INTEG4(MAXAIM) on the heap;
     therefore there is always enough space provided.  (ii) The entries for
     the IMAGE atoms have to be initialized, and this has to be done at
     EVERY image update.  However, similar things are already done for a
     number of other quantities like masses, vdW params, charges etc.  All
     this is done among a number of other things in subroutine MKIMAT in
     upimag.src, where I have added an appropriate statement.
     
     Ad (2), changes to post-processing routines: Real/Image atoms
     leave/enter the simulation box/system dynamically.  Therefore, the
     nonbonded/image interaction lists have to be updated during
     post-processing.  The hooks were already in the program, subroutine
     BLUPLST.  The real changes hide in this routine, most changes in
     BLFREE, BLEAVG and BLCOMP are either cosmetic or ensure proper
     printout.  Post-processing routines FREE, EAVG and COMP will actually
     print IMAGE terms if present.  The routine BLUPLST is a sibling of
     routine updeci in heurist.src.  The heuristic update scheme itself is
     removed, as I feel that one should update the lists at every frame.
     Also, the CRYSTAL specific section of UPDECI is not present in BLUPLST
     as I don't understand it.  Therefore, care should be exercised when
     using BLOCK with CRYSTAL!  Negative values of INBFRQ/IMGFRQ are
     trapped, in this case they are set to 1; Printout from the update /
     list generation routines is suppressed by temporarily raising the
     PRNLEV to 1.
     
     The BLOCK documentation (:doc:`block`) has been revised and reflects
     these modifications.  A new testcase block3.inp has been added to
     test/c24test.

(17) Constraint correction for PERT (Stefan Boresch)

     The current version of PERT cannot handle situations where SHAKE
     is applied to bonds which change in length due to an alchemical
     mutation
     as SHAKE and PERT do not "communicate".  Furthermore, in such cases a
     constraint correction has to be computed and added to the free energy
     difference.  Two steps are required to fix this problem:

     (1) The constraint list needs to be updated as a function of the
         coupling parameter lambda.
     (2) The constraint correction has to be calculated.

     Only thermodynamic integration (both for slow-growth and
     windowing)
     is supported; the exponential formula will give nonsense results.  (If
     someone wants to fix this, please look at Pearlman/Kollman, JCP 1991,
     94, 4532 and Severance et al. J. Comput. Chem. 1995, 16, 311.)
     
     The method to calculate the constraint corrections is based on
     extracting the respective Lagrangian multipliers from the SHAKe
     routine; this approach is briefly described in van Gunsteren et al.
     Computer Simulation of Biomolecular Systems: Theoretical and
     Experimental Applications; ESCOM: Leiden 1994; Vol. 2, pp 315-348.
     The approach fully includes inertial contributions, it is left to the
     user to account for those correctly in the context of the problem.
     
     The new code is mostly transparent and does not really require
     additional documentation.  However, some information is added to
     :doc:`pert`.  A new testcase pert2.inp is also added to test/c24test.

(18) Non-Cubic Crystal Building Problem Fix (Wonpil Im and Ryszard Czerminski)

     The crystal build facility uses the symmetrized rotated shape matrix
     XTLABC obtained from lattice parameters.  However, it does not apply
     the same rotation to the unit cell moiety, which may result in bad
     contacts in non-cubic crystals.  The problem is fixed by calling the
     subroutine ROTXTL.  Some tests for the rotation are added by Ryszard.

.. _changelog_c24-c25:

Major Enhancements and Developments in CHARMM25
-----------------------------------------------

During the C25 development cycle, August 15, 1995 to August 15, 1997,
we made three bugfix-updates in the c24 releases and three alpha versions
and one beta version in the c25 development line.

::

        CHARMM24.0
             c24b2    Release           February  15, 1996
             c24g1    Release           August    15, 1996
             c24g2    Release           February  15, 1997

        CHARMM25.0
             c25a0    Developmental     August    15, 1995
             c25a1    Developmental     February  15, 1996
             c25a2    Developmental     August    15, 1996
             c25a3    Developmental     February  15, 1997
             c25b1    Release           August    15, 1997

Only bugfixes are incorporated into CHARMM24 and all new developments
and enhancements have been carried out with the CHARMM25 developmental
versions.  All modifications are thoroughly recorded in the
ChangeLog.c25 file and the following is the summary of new features
and major enhancements in CHARMM 25.


New Features in CHARMM25
^^^^^^^^^^^^^^^^^^^^^^^^

(1) Merck Molecular Force Field (MMFF) (Thomas A. Halgren, Ryszard Czerminski, Jay L. Banks,
    Bernard R. Brooks, and Youngdo Won)

    Merck Molecular Force Field (MMFF) developed by Tom Halgren at
    Merck has been implemented in CHARMM.  Ryszard introduced MMFF into
    c23f2, which made the c24x1 (February 15, 1994) version for evaluation.
    As CHARMM was evolved through the c24 development project, Jay
    incorporated MMFF into c24b1 in a less intrusive manner.  Bernie and
    other developers reviewed c24b1/MMFF and suggested some corrections.
    Youngdo took the Jay's code and Bernie's suggestions and made the
    checkin code of MMFF.  MMFF is documented in doc/:doc:`mmff`.

(2) CADPAC (Paul Lyne)

    An interface is added to allow CHARMM to run with CADPAC6.0 when
    performing QM-MM calculations. CADPAC6.0 can perform HF, MP2, MP3 and DF
    calculations. 

(3) Particle Mesh Ewald Code (Bernard R. Brooks)

    The Particle Mesh Ewald (PME) method has been implemented.  This
    code is based on code sent by Tom Darden at NIEHS/NIH.  It has been
    modified so as to conform with CHARMM coding standards.  This version
    is much faster than the standard Ewald code and accuracy does not appear
    to be a problem when reasonable options are used.  This code uses the
    new "smooth" algorithm.  See :doc:`ewald` for more details.
    
    The code is now running in parallel and the following features are
    supported:
    
    - PERT (free energy calculation) with PME (including pressures)
    - Assymetric units with CRYSTAL (when NOPEr>0 in CRYStal BUILd command)
      is now supported with PME.
    - Total charge (Qtot<>0) energy and pressure correction term has been added.
    - Accurate pressures for the triclinic (and all other) cases (and for PERT)
    - Ewald energy components have been separated and can be turned off
      with the SKIP command ('EWKS','EWSE','EWEX','EWQC','EWUT').
      (k-space,self term,exclusion,total Q correction,utility)

(4) External Force to Selected Atoms (Lennart Nilsson)

    A new command has been added which calculates a new energy term
    corresponding to a static or periodically varying external force on an atom
    selection.

(5) Distance matrix and radius of gyration restraints (Charles L. Brooks, III, Felix B. Sheinerman and Erik Boczko)

    New restraint energy terms added to permit restraint of system based
    on its radius of gyration and/or the value of a reaction coordiante what
    describes the degree of nativeness based on the number of native side
    chain contacts.  New Keywords are RGYCONS and DMCONS.

(6) HTML Doc Files (Rick Venable and Charles L. Brooks, III)

    Added html documentation files and developed/modified doc2html.com
    originally developed at NIH.  All relevant files added to support/htmldoc

Major Enhancements in CHARMM25
------------------------------

(1) PARALLEL CODE reorganized and extended (Milan Hodoscek, Charles L. Brooks)

    The parallel code has been updated and organized into three parts:
    paral1.src, paral2.src and paral3.src.  The new code is faster and there
    has been a significant additions to support other platforms.  We now
    support about 15 platforms including ALPHAMP, T3D, T3E, Terra,
    Global-Works-Server and others.

(2) Linux Port (Milan Hodoscek)

    The Linux port is done with the GNU Fortran compiler, version 0.5.18.
    For now, all Linux related changes are under the GNU keyword.

(3) Nonbond Energy Code Overhaul with Semi-Automatic Code Expansion (Bernard R. Brooks)

    The program PREFLX (PREFX) has been overhauled to allow semi-automatic
    code expansion in the moving of inner loop if-tests to the outside of
    do-loops.  The nonbond energy routines are cleaned and organized to 
    utilized the semi-automatic expansion.  Obsolete ZTBL code is removed.

(4) Ewald code (Bernard R. Brooks)

    Memory needs of the EWALD electrostatic option have been reduced,
    and multiple parallel options are now supported.  Pressure code has been
    fixed as well.
    
    The calling sequence to ENBOND was modified so that a flag (QEWEX)
    can be sent indicating whether the nonbond exclusion correction should be
    performed for the Ewald calculation.  This corrects several problems
    (such as Ewald with MTS and Ewald with PERT) and this simplifies some code
    relative to the handling of the exclusion lists.  Also there were several
    changes to EPERT so that the Ewald method will report a correct internal
    virial (for pressure).  The Ewald method was enabled for the GROUP option
    so that group lists can be used.  This reduces the amount of memory and the
    time needed to handle the nonbond lists (good for limited memory parallel
    machines).
    
    A version of EWALD was developed for MMFF.  The usual MMFF electrostatic
    term: qq/(r+d)  is split into two terms:  qq/r -  qq*d/(r*(r+d))  The first
    term is handled by the Ewald method in the usual manner (real-space and
    k-space parts) and the second term is truncated at the cutoff distance
    using a switching function (from CTONNB to CTOFNB).  Since the second term
    is quite small at the cutoff distance, the use of a switching function
    should not introduce significant artificial forces.

(5) Restrained Distance Code Enhancement (Bernard R. Brooks)

    The restained distance method has been extended to allow the use of
    a one sided function (positive or negative).  It also allows a non-unit
    exponent for the individual distance terms.  The code is now much more
    general in its ability to define distnace based retraints based on multiple
    distances.

(6) COOR DIPOLE (Bernard R. Brooks)

    A COOR DIPOle command has been added.  This command computes tha charge
    and dipole (multipoles) for selected atoms.

(7) READ PSF APPEnd (Bernard R. Brooks)

    The READ PSF APPEnd command option has been modified so that it does
    not initialize the coordinates of existing atoms.  Only the new appended
    atoms will have undefined coordinates.

(8) Replica within Images (Bernard R. Brooks)

    The replica code has been enhanced so that it workes with images and the
    crystal facility.

(9) The crystal facility has been extended (Bernard R. Brooks)

    The follwing new features have been added:
    
    - The "DODE" has been renamed "OCTA" (for truncated OCTAhedron).
      (pressure bug fixed for OCTA)
    - A new type "RHDO" has been added (for RHombic DOdecahedron).
    - The CRYSTAL BUILD command is now much faster and more accurate.
      The use of the double atom search has been limited.
    - The documentation has been updated to give detailed information
      regarding crystal types.
    - The WRITe/PRINt IMAGE command is no longer iterative (in accord
      with the existing documentation).

(10) Overhaul of Harmonic restraints (Bernard R. Brooks)

     The CONS HARM command has been overhauled and extended.  The new syntax
     has three different types of harmonic restraints:

     ::
     
         CONStraint HARMonic { [ABSOlute]  absolute-specs    }  force-const-spec
                             {  BESTfit    coordinate-spec   }
                             {  RELAtive  2nd-atom-selection }
                             {  CLEAr                        }

     The ABSOlute is the old method.  The BESTfit causes the reference set to be
     logically bestfit rotated/translated before computing the restraint energy.
     The RELAtive allows two portions of one PSF to be restrained to the same
     internal geometry by the bestfit least squares rotation (no reference
     coordinates used).

     Some features and changes:
     
     - Multiple restraints (same or different types) are allowed.
     - HARMonic restraint I/O is no longer supported.
     - The old command syntax still functions (no rewrite of scripts required).
     - The READ/PRINt/WRITe CONS commands now have a "PSF 0" option for PERT.
     - PERT supports all of these restraint types.
     
     Restriction:
     
     - Each atom may participate in AT MOST one harmonic restraint term.

(11) Enhancement to REPLICA/PATH (Bernard R. Brooks)

     The REPLICA/PATH method has been extended to allow for bestfit translation
     and/or rotations between adjacent replicas before computing the restraint
     energies.  Getting the forces right was the hard part.   This allows entire
     molecules to be replicated (or sections with significant freedom).

     ::
     
         RPATh  [ KRMS real ] [ KANGle real ] [ COSMax real ] [MASS] [WEIGht]
                      [ KMAXrms real ] [RMAXrms real ] [ ROTAtions ] [ TRANslations ]


(12) CHARMM/GAMESS enhanced (Milan Hodoscek)

     The version of GAMESS has been updated to the March-97 version
     from Ameslab.  Also, QM/MM gaussian blur of MM charges has been
     implemented as an option.

(13) A few small changes to MMFP and NMR (Benoit Roux)

     A few small changes to some MMFP subroutine have been made.  The
     main thing is a second atom select for the SSBP command that allows
     the present of atoms outside the boundary radius.  This could be
     useful when the boundary is used only for an active site.   The
     relaxation time due to the chemical shift anisotropy addition has
     been added to NMR.

(14) COOR DMAT (Charles L. Brooks, III)

     The dist keyword has been removed from the covariance command and
     a new analysis command has been added under the coor subsyntax.  This
     command is accessed with the command COOR DMAT and provides some
     general tools for the calculation, manipulation and storage/extraction
     of distance matrix based properties.  This routine has some overlap
     with the new distance command introduced by Bernie Brooks but also
     provides significant complementarity in extending the range of
     properties computed.



