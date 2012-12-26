.. py::module:perturb

=====================================================
Perturbation: Thermodynamic Perturbation Calculations
=====================================================

See also:

* Details: How to run perturbation calculations. :doc:`pdetail`
* Implementation: How it is implemented. Programming details. :doc:`pimplem`
* CFTI: Conformational Energy/Free Energy Calculation (Krzysztof Kuczera) :doc:`cfti`

.. _perturb_syntax:

Syntax for the Perturbation Command
===================================

Syntax of the set up of the perturbation command.

::
   
   [SYNTAX TSM]

   TSM

   		Chemical Perturbation Parameters:

   1.  REACtant atom_selection_list | NONE

   2.  PRODuct atom_selection_list   | NONE

   3.  LAMBda <real> [ POWEr <int> ]

   4.  SLOW TEMP <real> LFROm <real> LTO <real> [ POWEr <int> ]


   5.  DONT {REACtant} {internal_energy_spec} [SUBTract]
            {PRODuct} {internal_energy_spec}

   6.  GLUE {CM FORCe <real> MIN <real>} [SUBR] [SUBP]
            {ATOMs FORCE <real> MIN <real> atom_spec atom_spec

   7.  NOKE {REAC}
            {PROD}

   8.  SAVE UNIT <integer> [FREQ <integer>]

   9.  COLO atom_spec PCHArge <real> [RCHArge <real>]
            atom_spec ::= segid resnum type

   10. PIGGyback PIGGy atom_spec BACK atom_spec
            atom_spec ::= segid resnum type

   11. UMBRella 4x( atom_spec) VACTual <real>
              atom_spec: segid resnum type


   	Internal Coordinate (IC) Perturbation Parameters:

   12. FIX  {ic-spec} [TOLI <real>]

   13. MAXI <integer>

   14. MOVE {ic-spec} 2x{atom-selection} BY <real>

   15. SAVIc [ICUNit <integer>] [ICFReq <integer>] [NWINdows <integer>]

   16. END

   	internal_energy_spec ::== BOND THETa|ANGLe PHI|DIHEd IMPHi|IMPR

   	ic-spec ::=	{[DISTance] 2x{atom-spec} }
   			{[BOND] 2x{atom-spec}     }
   	   	 	{[ANGLe] 3x{atom-spec}    }
   		 	{[THETa] 3x{atom-spec}    }
   		 	{[DIHEdral] 4x{atom-spec} }
   		 	{[PHI] 4x{atom-spec}      }

    	atom_spec ::= segid resid type

   	atom-selection ::= see (*Note Select: (SELECT).)


   ***** Note: must have non-bonded exclusions between reactant and product 
         atoms in rtf.

   -----------------------------------------------------------------------

   TSM CLEAr
         Clears heap data structures used in perturbation setup, cancels
         constraints and perturbations, and resets logical flags.


.. _perturb_description:

Explanation of the Perturbation Setup
=====================================

Currently  the  perturbation  setup  is initiated by invoking the
command  TSM with nothing else on the command line.  This is followed by
a  number  of  other  commands,  listed below, and terminated with an END
command.  Two types of thermodynamics perturbations are available: chemical
perturbation and internal coordinate perturbation.  Each is discussed
separately below.


.. _perturb_chempert:

Chemical Perturbation
---------------------

For chemical perturbations,  a minimum of three commands are necessary 
besides TSM and END:  REAC - to specify the reactant atom list; PROD - to 
specify the product atom list;  LAMBda or SLOW to specify lambda for windowing
or the slow growth technique.

::

   1.  REACtant atom_selection_list | NONE
 
           Specifies the reactant atom list (see *Note details: (pdetail).).
   The atom selection list uses the standard CHARMM selection command syntax
   (see  *Note  Select: (SELECT).).   Subsequent invocations of this command
   clears the selections of any earlier invocation.
 
   2.  PRODuct atom_selection_list | NONE
 
           Specifies the product list (see above).
 
   3.  LAMBda <real> [ POWEr <int> ]
 
           The hybrid Hamiltonian is defined, in this implementation, as
 
           H(lambda) = ( (1 - lambda)**N )V(reac) + (lambda**N)V(prod).
 
   This  command  specifies lambda and N.  It also indicates that the window
   method is to be used (see *Note details: (pdetail).).
 
   4.  SLOW TEMP <real> LFROm <real> LTO <real> POWEr <int>
 
           This   command    specifies    that   the   "slow   growth"  (see
   *Note  details: (pdetail).)  method  be used. LFROm and LTO indicates the
   limits  of  integration.   POWEr  has  the  same  meaning in the previous
   command.
 
   5. DONT {REACtant} {internal_energy_spec} [SUBTract]
           {PRODuct} {internal_energy_spec}
           internal_energy_spec :== BOND THETa|ANGLe PHI|DIHEd IMPHi|IMPR
 
           This command indicates that the specified internal energy term(s)
   for  the reactant or product atoms is (are) to be ignored as perturbation
   interactions.   That  means  that  the  specified  interactions  are  not
   factored  by  lambda**N  or  (1 -lambda)**N  and do not contribute to the
   value of V(reac) or V(product).  The interaction is, however, computed in
   full and treated as part of H(env) (see *Note details: (pdetail).).  More
   than one internal energy type can be specified at a time but reactant and
   product must be specified in separate commands.  The optional sub-command
   SUBTract  causes  specified  perturbation forces on the environment atoms
   (see   *Note   details:  (pdetail).) to  be  subtracted.   This  is  very
   non-Newtonian  and  was  included  as an early (and largely unsuccessful)
   attempt  to  generate  configurations  at  lambda  =  0  and  1  for  the
   non-existent  group.   See  the  discussion  of the endpoint problem, (in
   *Note  details:  (pdetail).),  also  known  as  the  lambda  goes to zero
   catastrophe  (Beveridge,  1987).   Anyway,  what  this does is remove the
   forces  due  to  terms  specified in the DONT option from the environment
   atoms  involved.   In  addition,  the  terms  are  not lambda factored or
   included  in  H(lambda)  or in V(reac) or V(pert). The forces are left on
   the   perturbed  atoms  with  the  hope  that  this  would produce usable
   configurations.  It  does  not.   Rather,  the  atoms  drag along and bad
   non-bond  contacts  result  when  lambda is 0 or 1 (the original intended
   use).
 
           We  generally  use  the DONT option for both reactant and product
   for  the  bond  stretching  and bending terms.  These terms are generally
   uncoupled  from  the  interactions  of interest and it appears that their
   exclusion,  even  with  the  resulting non-physical Hamiltonian, does not
   significantly   affect   the  relative  free  energies  of  solvation  or
   drug/enzyme  binding.   The  same cannot be said for torsions which other
   implementations leave out.
 
           The  educated  advice  is use the DONT REAC (and PROD) BOND THETA
   and  donot!!!  use  the  SUBTract option.  Note that we have had problems
   with  the  fraying  of  bonds  to hydrogen near the lambda endpoints when
   the DONT BONDs options were NOT used.
 
           Note  that the REAC and PROD must be issued before the respective
   DONT  command.  Subsequent invocations of a DONT REAC/PROD command clears
   the applicable flags first.
 
   6. GLUE {CM FORCe <real> MIN <real>} [SUBR] [SUBP]} |
           {ATOMs FORCE <real> MIN <real> atom_spec atom_spec }
           atom_spec ::= segid resnum type
 
           Here  we have another failed attempt.  The GLUE ATOM command sets
   an  harmonic  force  between  a reactant and a product atom.  One of each
   must  be  given.   The  GLUE CM indicates that the centers of mass of the
   reactant  and  product  atoms  are to be connected by the harmonic force.
   FORCe  is  the  force  constant  in the same units as the bond stretching
   force  constant kcal/mol/A**2.  MIN is the minimum distance in angstroms.
   Unfortunately, using MIN set to zero causes problems with SHAKE (how does
   a  floating  point  zero  divide  check  error  grab  you, Buckaroo?). We
   intended  to use this on systems where there were no environment atoms to
   keep  the groups together.  Shake would be used since the "GLUE" force is
   unphysical.   However,  the  aforementioned error made use of this option
   undesirable.  As it turns out, in solvated systems our concerns about the
   two  groups  flopping  around,  with  attendant  problems  with  sampling
   convergence, was unfounded.  Best not to use this option.
 
   7.  NOKE {REAC}
            {PROD}
 
           This  specifies that the kinetic energy should not include either
   contributions  form  the REACtant or PRODuct atoms.  REACtant and PRODuct
   must  be  selected  in  separate commands.  When the number of degrees of
   freedom are calculated in the subroutine DCNTRL, the ones due to REACtant
   or PRODuct (depending on which command(s) is/are issued) are not counted.
   When  the  kinetic  energy  is calculated in the subroutine DYNAMC, these
   degrees  of  freedom  are  ignored. Same for the temperature calculation.
   This  option  should be used for the non-existent  atoms  at  lambda  = 0
   (product  atoms) or 1 (reactant atoms) since the atoms donot exist (hence
   they  are  termed  non-existent  buckaroo) they should not be expected to
   contribute their (3/2)kT per degree of freedom to the kinetic energy.
 
           Normally,  the kinetic energy contributions from the reactant and
   product atoms are factored by (1 - lambda)**N and lambda**N respectively.
 
   8. SAVE UNIT <integer> [FREQ <integer>]
 
           This  command  determines where the output generated in DYNAMC is
   to be sent (UNIT command) and the frequency of output (FREQ command).  We
   generally  use  a  frequency  of  one  (output on every step).  Note that
   before  dynamics  are  run  the file must be opened for formatted writing
   with  the CHARMM OPEN command.  Binary output is not currently supported.
 
   9. COLO atom_spec PCHArge <real> [RCHArge <real>]
           atom_spec ::= segid resnum type
 
           Sometimes  the van der Waal characteristics and the "identity" of
   of  an atom remains the same while the charge interactions with this atom
   atom  is  perturbed from a reactant value to a product value.  An example
   is  the  methanol  ->  ethane  mutation  where  the methanol OH atoms are
   full-fledged   reactant   atoms   and   one  ethane  methyl  group  is  a
   full-fledged  product atom.  The common methyl group is treated as a COLO
   atom.  Interactions appropriately factored by lambda are calculated using
   a product charge, specified by the PCHArge command and a reactant charge,
   either the charge in the residue topology file or that specified with the
   RCHArge  command.   The  COLO  command is issued once for each COLO atom.
   The colo atoms cannot be in either the reactant or product lists.  If any
   of this is unclear see *Note details: (pdetail). An explanation about how
   this whole thing works is given there.
 
   10. PIGGyback {PIGGy} atom_spec {BACK} atom_spec
                 {REACtant}        {PRODuct}
            atom_spec ::= segid resnum type
 
           This  command  allows  you to convert isolated atom into another.
   It  is intended for mutations such as Br- -> Cl- and Ar -> Ne.  The atoms
   cannot  be  bonded  to  anything  else  (a  wrndie error will be issued).
   Currently,  only  one  reactant/product atom pair can be used.  The PIGGy
   atom  must  be  a REACtant atom and the BACK atom must be a PRODuct atom,
   otherwise  an  error  will  be  flagged.  The REAC and PROD commands must
   already have been issued.  Note the synonyms for PIGGy and BACK.
 
           When  this  option is  in  effect the forces on the back atom are
   added  to  those  on  the  PIGGy  atom at each step of the dynamics.  The
   coordinates  of  the  BACK atom are made equal to those of the PIGGy atom
   at each time step.
 
   11. UMBRella 4x( atom_spec) VACTual <real>
              atom_spec ::= segid resnum type
 
           This  command specifies an umbrella sampling correction to all of
   the  averages in post-processing. The four atom specifications define the
   the  dihedral angle involved.  The command is repeated for each dihedral.
   If  there are multiple dihedral angles through the axis of two atoms, all
   all  of  them  should  be  specified.  It  is  assumed that the surrogate
   potential  term  is  in  the  parameter  file  for the particular type of
   torsion.   VACTual is the coefficient for the real potential.  Currently,
   only  the  three-fold  term  is  supported.  A further limitation is that
   although  you  can  specify particular dihedral angles for this treatment
   all  torsions  with  that  type  will  use  the modified potential in the
   parameter  file.   This  part  of  the program is slated for modification
   as   soon  as  possible.   For  an  explanation  of  the  terms  and  how
   the  umbrella correction works, see *Note details: (pdetail).
 
.. _perturb_icpert:

Internal Coordinate (IC) Perturbation
-------------------------------------

To setup an ic perturbation, you need to  1) specify the internal 
coordinate(s) to be constrained during the perturbation (FIX), 2) specify 
which atoms will move during the perturbation and which atoms will remain 
fixed (MOVE), and 3) indicate how and where the perturbation data will be 
saved (SAVI).

::

   12. FIX  {ic-spec} [TOLI <real>]

   	ic-spec ::=	{[DISTance] 2x{atom-spec} }
   			{[BOND] 2x{atom-spec}     }
   	   	 	{[ANGLe] 3x{atom-spec}    }
   		 	{[THETa] 3x{atom-spec}    }
   		 	{[DIHEdral] 4x{atom-spec} }
   		 	{[PHI] 4x{atom-spec}      }

           The FIX command defines an internal coordinate to be constrained: 
   DIST or BOND specify distance constraints, ANGL or THET bond angle 
   constraints, and DIHE or PHI dihedral angle constraints.  TOLI sets the 
   tolerance (the allowed deviation in Angstroms (distance constraints) or 
   degrees (angle constraints)) to be used for the specified constraint in 
   the constraint resetting procedure (see *Note details: (pdetail).).  The 
   default is 10E-10 (Angstroms or degrees).  It is very important to note 
   that the reference value of the constraint is set to the value of the 
   internal coordinate at the instant the command is issued.  One or more 
   i.c. constraints can be specified per simulation.

   13. MAXI <integer>

           MAXI sets the maximum number of iterations to be used in the 
   iterative i.c. constraint resetting procedure (see *Note details: 
   (pdetail).).  The default is 500.

   14. MOVE {ic-spec} BY <real> INTE {atom-selection}

           The MOVE commands specify the internal coordinates to be perturbed 
   and define the atoms to be moved by the perturbation.  Several MOVE 
   commands may be used to set up a perturbation consisting of changes in 
   several internal coordinates.  Since i.c. perturbations are really only 
   useful in conjunction with i.c. constraints, for each MOVE command there 
   should be a corresponding FIX command with the same ic-spec.  Following 
   the BY keyword is a real number which is the amount that the internal 
   coordinate will be changed by the perturbation (in Angstroms for distances 
   and degrees for angles).

           The INTE selection part of the MOVE command defines the solute 
   partition, that is, the atoms to be moved by the perturbation.  The 
   perturbation is applied to all of the atoms specified in the selection 
   using displacements determined from moving the internal coordinate BY 
   value.  However, some atoms may not be moved even though they are included 
   in the solute partition with the INTE selection because zero displacements 
   will be computed for them (e.g. if they lie on the rotation axis, like the 
   central atom in a perturbed angle, or either of the two central atoms in a 
   perturbed dihedral angle).  If a double selection is given (e.g. INTE SELE 
   atom-selection END SELE atom-selection END), then the two selected groups 
   of atoms are considered separate sections of the solute partition.  In 
   that case, to accomplish the overall perturbation, each section is moved 
   half of the BY value.

   	The INTE selection also specifies which contributions are to be 
   included in the perturbation interaction energies.  The calculation of the 
   perturbation interaction energies is based on the interaction energy 
   calculation which is done when the CHARMM INTEre command is issued (see 
   *Note interaction: (energy). for more details).  Thus, the perturbation 
   interaction energies may contain the following energy contributions: bond, 
   bond angle, dihedral, improper dihedral, van der Waals, electrostatic, 
   hydrogen bond, harmonic positional constraint, and harmonic dihedral 
   constraint.  In addition to these contributions, which are the usual 
   CHARMM interaction energy terms, the perturbation interaction energies may 
   also contain the following image contributions: van der Waals, 
   electrostatic, and hydrogen bond.  (Note that the CHARMM INTEre command is 
   parsed by the main CHARMM command parser.  It should not be confused with 
   the INTE part of the MOVE cammand which is parsed by the TSM command 
   parser.  We apologize for any confusion which may result from the use of 
   the INTE keyword in the TSM command.  It seemed appropriate since it 
   indicates a similar interaction energy calculation.)

   	To explain how the interaction energy calculation works, we define 
   two "selection groups".  The first selection group contains all of the 
   atoms in the system.  The second selection group contains all of the atoms 
   included in the INTE selection.  The rules which are used to determine 
   which contributions are included in the interaction energies are as 
   follows: a bond is included if the two atoms defining the bond are in 
   different selection groups;  a bond angle if the central atom is in both 
   selection groups; a dihedral angle (intrinsic torsion or harmonic 
   constraint) if the two central atoms are in different selection groups; an 
   improper dihedral angle if the first atom is in both selection groups; a 
   nonbonded interaction (van der Waals and electrostatic) if both atoms are 
   in different selection groups and the interaction is in the nonbonded 
   list; a hydrogen bond if the donor and acceptor are in different selection 
   groups; and finally, a harmonic positional constraint is included if the 
   atom is in both selection groups.

   	The user should decide carefully which interaction energy contribu-
   tions she wants to have included before running the perturbation simula-
   tion.  Then she must appropriately design the INTE selection.  For 
   example, suppose she wants to compute the free energy as a function of the 
   dihedral angle in an extended-atom (four atom) model for butane.  A change 
   in the dihedral angle only changes the position of the methyl group.  The 
   user might therefore select only a terminal methyl group (e.g. segment 
   BUTA, residue 1, atom C4) using the INTE command:

   INTE SELE (ATOM BUTA 1 C4) END

   With this selection, the intrinsic torsional potential would not be 
   included in the interaction energies since the CHARMM interaction energy 
   routine only computes torsional terms if the central two atoms of the 
   torsion are in different selection groups.  Of course, this is not 
   generally a problem, since the missing term could be simply added to the 
   thermodynamics after processing the interaction energies.  If the user 
   preferred to have the intrinsic torsional contribution included in the 
   interaction energies, she would add the methylene group (atom C3) to the 
   INTE selection, e.g. she could use the following selection in place of the 
   one above:

   INTE SELE ((ATOM BUTA 1 C3) .or. (ATOM BUTA 1 C4)) END

   With this specification, the C3 atom is included in the solute partition.  
   However, its position is not changed by the perturbation since it lies on 
   the axis about which the solute atoms are rotated in the dihedral angle 
   perturbation.  Now the two central atoms, C2 and C3, are included in 
   different selection groups, so the intrinsic torsional contribution is 
   included in the interaction energies.
   	There is a subtle point that must be considered when the perturba-
   tion consists of moving more than one internal coordinate.  As an example, 
   suppose we want to perturb both the dihedral angles, which we call phi and 
   psi, in an extended atom (five-atom) model for pentane.  Further suppose 
   that, in the double-wide sampling, we want the perturbation in one 
   direction to increase both phi and psi, and the perturbation in the other 
   direction to decrease them.  We might try the following MOVE commands 
   (e.g. +/- 5 degree perturbations of each dihedral angle; C2 and C4 are 
   selected so the intrinsic torsion terms are included in the interaction 
   energies):

   MOVE DIHE PENT 1 C1 PENT 1 C2 PENT 1 C3 PENT 1 C4 BY 5.0 -
     INTE SELE ((ATOM PENT 1 C1) .OR. (ATOM PENT 1 C2)) END
   MOVE DIHE PENT 1 C2 PENT 1 C3 PENT 1 C4 PENT 1 C5 BY 5.0 -
     INTE SELE ((ATOM PENT 1 C4) .OR. (ATOM PENT 1 C5)) END

   However, in the algorithm which changes bond and dihedral angles, a 
   perturbation in the forward direction corresponds to a counterclockwise 
   rotation of the atoms to be moved around the bond vector (e.g. C2�C3 or 
   C3�C4 in pentane).  Thus, with the above MOVE commands, the forward 
   perturbation decreases phi while it increases psi.  That is not what we 
   wanted.  To fix the problem, we simply reverse the sign of the BY value in 
   one of the MOVE commands:

   MOVE DIHE PENT 1 C1 PENT 1 C2 PENT 1 C3 PENT 1 C4 BY 5.0 -
     INTE SELE ((ATOM PENT 1 C1) .OR. (ATOM PENT 1 C2)) END
   MOVE DIHE PENT 1 C2 PENT 1 C3 PENT 1 C4 PENT 1 C5 BY �5.0 -
     INTE SELE ((ATOM PENT 1 C4) .OR. (ATOM PENT 1 C5)) END

   Now the forward perturbation decreases both phi and psi and the reverse 
   perturbation decreases them.  The user should carefully consider how the 
   atoms will be moved when choosing the signs of the BY value when more than 
   one internal coordinate is perturbed.


   15. SAVIc [ICUNit <integer>] [ICFReq <integer>] [NWINdows <integer>]
             [SUPP]
          (for "on-the-fly" free energy and average energy calculations:)
             [RUNA] [PEVEry <integer>] [RUNIt <integer>] [RPRInt <integer>]
   	  [TEMP <real>]
           The SAVI command specifies how and where the perturbation data 
   (which consists primarily of internal coordinate values and interaction 
   energies) will be saved during the simulation.  The integer following the 
   ICUN keyword is the number of the fortran unit to which the perturbation 
   data is written.  The perturbation file should be opened for formatted 
   writing on this unit using the CHARMM OPEN command before the dynamics 
   command is issued.  The integer following ICFR is the frequency (in 
   dynamics steps) with which the data is written to the file.  The default 
   ICFR value is 10.  If the frequency is zero (e.g. if the SAVI command is 
   not issued) or ICFR 0 is indicated, then a level 0 warning is issued since 
   there is no need to do perturbations if the data is not going to be saved.  
   The integer, m, following the NWIN keyword indicates the number of 
   "double-wide" subintervals that the BY value, dx, will be divided into.  
   Thus, the 2m perturbations, dxi = i*dxm, where dxm = dx/m and i = -m,-
   m+1,..., -1, 1,...m-1,m, are all carried out during the simulation, 
   yielding 2m free energy differences.  For example, if the BY value is 1.0 
   and the NWIN value is 2, perturbations which change the internal 
   coordinate by -1.0,-0.5,0.5, 1.0 are carried out.  The default NWIN value 
   is 1.  The SUPP keyword suppresses printing of the internal coordinate
   values to the output file.
          Running or "on-the-fly" free energy changes and average energy 
   changes can be calculated for internal coordinate perturbations through
   the invocation of the "RUNAverage" keyword. The routine will include all
   data points (i.e. every ICFR'th step in the trajectory) in these cal-
   culations.  The results will be written to the specified file
   in RUNIt every RPRInt sampled data points (default is no writing out of
   averages).  Hence for ICFR of 5 and RPRInt of 10, the averages will be 
   calculated every 5 steps and written out every 50 (RPRInt*ICFR) steps.
   (See testcase for examples.) TEMP specifies the temperature in degrees
   Kelvin at which the free energy is to be calculated (default 300).
   PEVEry specifies the period (in ICFR number of steps) for writing out
   the usual tsm output file containing the energies and the internal 
   coordinates (default 1--file is written every ICFR steps).

   The "on-the-fly" output file is formatted as follows:
   RUNAV>              5  168.5417     1.84955880     2.02536215
   RUNAV>              5  168.0417     0.84931404     0.90105846
   RUNAV>              5  167.0417    -0.69647324    -0.62795693
   RUNAV>              5  166.5417    -1.21666615    -0.91659629
   RUNAVI>   1  167.54174244
   RUNAVI>   2  155.47953779
   RUNAV>             10  170.7559     1.62569412     1.84119225
   RUNAV>             10  170.2559     0.77839900     0.83143112
   RUNAV>             10  169.2559    -0.67517686    -0.62242291
   RUNAV>             10  168.7559    -1.20371332    -0.99498136
   RUNAVI>   1  169.75592284
   RUNAVI>   2  155.45919138

   For the "RUNAV>" lines, the 2nd column gives the number of data points
   included in the averages. The second line gives the average value of
   the internal coordinate after a particular perturbation. (For perturba-
   tions involving multiple internal coordinates, the value of only the
   only the first internal coordinate specified in the input file is given).
   The third and fourth columns give the free energy change and average
   energy change, respectively, for the given perturbation.
   For the "RUNAVI>" lines, the first column gives the number of the
   internal coordinate (in its order of appearance in input file) and
   the second column gives the average value for that coordinate over
   the unperturbed trajectory.  All coordinates involved in the 
   perturbations (i.e. specified by the MOVE command) are listed.

   16.  END
 
           Terminates  the  perturbation  setup.   At this point the program
   does additional error checking and prints out values of some parameters.
 
   -------------------------------------------------------------------------
 
   TSM CLEAr
 
           A  separate  command  ( NOT!!  a setup command ) to clear logical
   flags and release HEAP memory allocated for perturbation data structures.
   It  is  not  necessary  to use this command unless you have more than one
   dynamics  run  in  a  single  job  and  want  to  reset  or  turn off the
   perturbation.    Definitely  invoke  this  command  before  entering  the
   perturbation setup a second time.
 

.. _perturb_post_processing:

Post-Processing of Perturbation Output
======================================
 
Syntax for Post-Processing Commands
-----------------------------------
 
::

    1.   TSM POST [PSTAck <int>] [PLOT] [TI] [NODEriv] [COMPonents] [ENDPoints]
 
                  [IC] [MAXP <integer>] [MAXW <integer>] [SURF] [MAXS <integer>]
                       [NODEriv] [INTE]
 
    2.   PROCess FIRSt <int> [NUNIt <int>] BINSize <int>
                 [CTEM] [TEMP <real>] [DELTa <real>] 
                 [BEGin <integer>] [STOP <integer>] [SKIP <int>] [NMAX <int>]
                 LAMBda <real> [ONE] [ZERO] [UMBRella] [EAVG] 
 
    3.   END
 

.. _perturb_ppost:

Description of Post-Processing Commands
---------------------------------------

Post-processing of perturbation data is initiated by the following 
command:

::
 
   TSM POST [PSTAck <int>] [PLOT] [TI] [COMPonents] [ENDPoints]
 
            [IC] [MAXP <integer>] [MAXW <integer>] [SURF] [MAXS <integer>]
                 [INTE]

            [NODEriv]


Summary of Parameters:
----------------------

1. Chemical Perturbation Post-processing Parameters

   ::
   
         PSTAck:        Array size for plotting x,y points. Default 100.
                        Needed for thermodynamic integration and/or
                        plotting.
         PLOT:          Create PLT2 output.
         TI:            Thermodynamic integration:
                        Delta A = int 0 to 1 <dE/dLambda>.
         COMP:          Do Vdw, Elec and Intern component analysis.
         ENDP:          Calculate TI integral to endpoints, i.e., full 
                        (0,1).


2. Internal Coordinate (IC) Perturbation Post-processing Parameters

   ::
   
         IC:            Specifies that ic perturbation output will be 
                        processed.

         MAXP:          Maximum number of ic perturbations.
         MAXW:          Maximum number of ic perturbation windows (NWIN in
	 	 	SAVI command).
         SURF:          Generate thermodynamic surfaces from ic 
                        perturbation data.
         MAXS:          Maximum number of points in surface.
         NODEriv:       Only calculate the free energy.
         INTE:          Calculate average interaction energies from ic 	
			perturbation data.

3. NODEriv Subcommand

   ::
   
         NODEriv:       Only calculate the free energy.


Chemical Perturbation Post-processing
-------------------------------------

There  are  two methods of calculating the relative free energies
and  relative  temperature derivative properties: the perturbation method
and the Thermodynamic Integration method (see :doc:`Detail <pdetail>`.).
By  default  the  perturbation method is used.  The optional parameter TI
specifies  the  thermodynamic  integration  technique.  PSTAck determines
the  size  of  arrays  to  be  allocated from the stack.  For the default
perturbation  method  this  allocates  space for plotting values.  The TI
method   requires  arrays  for  actually  calculating  the  thermodynamic
properties.   The  parameter  PLOT  indicates  that output is created for
PLT2  data  files.   They must be edited out of the output file.  NODEriv
indicates  that  only the free energy is to be calculated and not Delta E
or Delta S.


Internal Coordinate Perturbation Post-processing
------------------------------------------------

The MAXP, MAXW, and MAXS parameters are used to allocate memory for 
the processing of the perturbation data.  The integer following the 
keyword MAXP is the maximum number of perturbed internal coordinates in 
the data files to be processed (e.g. the number of MOVE commands in the 
perturbation dynamics input files).  The default MAXP value is 1.  The 
integer following MAXW is the maximum number of subintervals in each 
window (e.g. the NWIN value on the SAVI command command line in the 
dynamics input files).  The default MAXW value is also 1.  If the SURF 
keyword is present on the TSM POST command line, then thermodynamic 
surfaces will be constructed using the thermodynamic differences computed 
from the perturbation data.  We say more about this below.  The integer 
following MAXS is the maximum number of points in a thermodynamic surface.  
If all of the N data files to be processed using PROC commands (see below) 
have the same number of subintervals, m, the MAXS value is equal to mN + 
1.  The default MAXS value is 100.

In addition to the free energy differences, the internal energy and 
differences are computed by default using finite-difference 
temperature derivatives.  However, the calculation of the derivative 
properties may be turned off using the NODE keyword.  If the NODE keyword 
is present on the TSM POST IC command line, the internal energy and 
entropy differences are not computed.  If the INTE keyword is present, the 
average interaction energies of the solute partition (the atoms specified 
by the INTE selection in the MOVE command) with the bath partition (the 
remaining atoms), as well as the average total energies in the unperturbed 
and perturbed systems are computed and printed in the output file from the 
post-processing run.  By default the average interaction energies and 
average total energies are not computed.

.. _perturb_pproc:

The TSM POST command is followed by one or more PROC commands which 
specify the processing of perturbation data files, and terminated with the 
END command.  The syntax of the PROC command is as follows:

::

   PROC FIRSt int [NUNIts int] BINSize int [CTEMp] [TEMP real] [DELTa real]

   LAMBda <real> [ONE] [ZERO] [UMBRella] [EAVG] [SKIP <int>]
          [NMAX <int>]
 
   [BEGIn int] [STOP int]


Summary of Parameters:
----------------------

========= =================================================================
FIRSt     Fortran unit number of first file.

NUNIT     Number  of  fortran i/o units.  One can 'tack' on trajectories
          from separate files in the manner used for trajectory commands
          in  correl.   The  files  must be opened for read access prior
          initiating the post processing  with the TSM POST command, and
          the units must be numbered consecutively, starting with FIRSt.
          This is due to the fact that the post-processor command reader
          does  not  handle  MISCOM commands (*Note Misc: (MISC).). This
          probably  should  be  corrected.   Only  formatted  files  are
          handled currently, remember to open the files as formatted. We
          also note that it  does not matter which order  multiple files
          in a  given window are processed.  The post-processing program
          does not check to see if the dynamics steps are contiguous. It 
          only  checks  to see that  all of the  files in a given window 
          have the same header, as they should.  Default = 1.

BINSize   The number of data points per bin for error calculation.

CTEMp     A   flag  to  indicate  that  average  temperature  is  to  be
          calculated.  Because  this is being calculated while the other
          averages  are being accumulated the temperature is not used in
          calculating the thermodynamic properties.   To use the average
          temperature in the thermodynamic calculations, the user has to
          process the data twice,  manually specifying the average temp-
          erature from the first processing run as the TEMP value in the
          second run. By default the average temperature is not computed.

TEMP      The temperature for calculating properties.  Default = 298.

DELTa     The  temperature  increment  for finite difference derivatives
          (calculate  delta  E  and  delta  S).   No  meaning  if  TI is
          specified.   A  level 0 warning is issued.   The default is 10
          degrees.
========= =================================================================

The following parameters are only used when processing chemical 
perturbation data.

========= =================================================================
LAMBda    Lambda prime.  For calculating <exp-beta(E(lambda')-E(lambda)>
          and  related  quantities.   No  meaning  if  TI  is specified.
          Level 0 warning issued.  (See *Note Details: (pdetail).).

ONE       Indicates  that lambda is exactly 1.  Overides input lambda in
          file.  This  is  used  only  in TI.  In the case on non-linear
          lambda  dependence  the  derivatives due to reactant terms are
          identically  zero.   This provides a solution to the lambda ->
          zero catastrophe.

ZERO      Indicates  that lambda is exactly 0.  Overides input lambda in
          file.   Commands  ONE  and ZERO are mutually exclusive and are
          used  only  in  the TI post processor. Level 0 warning issued.
          This  command  is  used only in TI.  In the case on non-linear
          lambda  dependence  the  derivatives  due to product terms are
          identically  zero.   This provides a solution to the lambda ->
          zero catastrophe.

UMBRella  A flag to indicate that umbrella sampling is used.  A check is
          made of a parameter line in each data file to see if  umbrella
          sampling was indeed used.

EAVG      calculate <Etot.> and uncertainty for this value of lambda
          ignored if TI.

SKIP      skip first nskip records.

NMAX      maximum number of points to read. skips skip number of values
          first.
========= =================================================================

The following parameters are only used when processing internal coordinate
perturbation data.

========= =================================================================
BEGI      specifies number of first dataset to use in accumulating
          averages.

STOP      specifies number of last dataset to use in accumulating
          averages.
========= =================================================================


Processing Chemical Perturbation Data
-------------------------------------
 
The PROCess  command is usually issued several times.  When using
the  perturbation  method  one  would  issue  it  at least once for every
lambda.  In all of our work so far, we have employed double-wide sampling
in that for each value of lambda whereupon dynamics are run, we "perturb"
both  up  and  down  from lambda to lambda prime (i.e. both less than and
greater  than  lambda,  definitely  see  :doc:`Detail <pdetail>`.). The
program  rewinds  the files after each PROCess  command.  For each lambda
->  lambda prime perturbation, a separate PROCess command is issued.  For
the   TI   method,   one   PROCess   command   per  lambda  is  used  and
<dE(lambda)/dlambda> is calculated.


 
Processing Internal Coordinate Perturbation Data
------------------------------------------------

The PROC command is used to specify the processing of perturbation 
data from a single window.  The PROC command is usually used several times 
and the results from the various windows are usually constructed into 
thermodynamic surfaces.

The user may specify that a subset of all the data read is to be 
used in the calculation of the averages and thermodynamics.  This option 
is useful for examining the convergence of the thermodynamic properties.  
The integers following the BEGI and STOP keywords are the numbers of the 
first and last datasets (not dynamics steps), respectively, to be used for 
processing the data for a given window.  By default, all of the data is 
used.  If the limits are set using the BEGI and STOP keywords, they are 
only used on the data processed by the particular PROC command which set 
the limits (e.g. the defaults are reinstated after each PROC command).

The integer following the BINS keyword is the number of datasets per 
batch, n, used in the calculation of the statistical uncertainties by the 
the method of batch averages (see :doc:`Detail <pdetail>`.).  This number 
must be specified as there is no default value.  We typically use BINS 
100.

.. _perturb_pend:

END
===

This  command  terminates the post-processing.  When this command
is  received  the  averages generated by the issuing the PROCess commands
are  combined  and  total  values  of  the  thermodynamic  properties are
computed  and  output.

For chemical perturbations,  if the  TI method is used,  a spline 
polynomial is fit to the averages and  integrated over  limits determined
by the minimum and maximum lambda's.  The lambda values do not have to be
processed  in  order  since the program will sort them.  It is the user's
responsibility to cover the whole range  for  lambda = 0 to 1 (if that is
the  intention).   Though  there  are  cases  where a range that does not
include  those  two  endpoints may be useful  (e.g.  mixing linear TI and
linear perturbation :doc:`Detail <pdetail>`.), discontinuous gaps in the
lambda curve do not make sense.  In the section :doc:`Detail <pdetail>`.
Input examples and explanations of the different methods are given.

For internal coordinate perturbations, if surface construction has 
been requested by including the SURF keyword on the TSM POST IC command 
line, the average internal coordinate and thermodynamic values (free 
energy, internal energy, and entropy differences) are sorted for construc-
tion of the surfaces.  The sorting is done according to increasing values 
of the average internal coordinate of the first perturbed internal coordi-
nate (e.g. the i.c. specified in the first MOVE command issued in the 
dynamics input file).  Then the thermodynamic surfaces are constructed by 
combining the differences in the thermodynamic properties.  Finally, the 
surfaces are printed to the output file of the post-processing run.  If 
more than one internal coordinate is perturbed, an identical set of sur-
faces is printed out as functions of each perturbed internal coordinate.  
The surfaces may be simply cut from the output file using an editor for 
subsequent plotting (e.g. using the PLT2 program).
