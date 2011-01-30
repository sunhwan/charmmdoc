.. py:module:: tps

========================
Transition Path Sampling 
========================

The Transition Path Sampling (TPS) methods introduced by Chandler and   
co-workers to sample rare events (see References) are implemented as 
extensions of the RXNCoor, DYNAmics, and USER commands in CHARMM.  The 
TPS keyword must be included in pref.dat for the code to be compiled.


.. _tps_syntax:

Syntax required to invoke TPS
-----------------------------

::

   RXNCoor [ standard RXNCoor keywords ] [ set-spec ] [ bas-spec ] [wri-spec]

   set-spec ::= SET [ NRXN     1 ] NRXN{ name }

                where NRXN{ name } is an NRXN-long list of the names of 
                the order parameters to calculate.

   bas-spec ::= BASIn NRXN{ name alo ahi blo bhi }

   wri-spec ::= TPUNit NRXN{ name unit }

   [Syntax DYNAmics]

   DYNAmics [ RTRJ ] [ standard DYNAmics keywords ] -
            PATH [ mod-spec ] [ tps-spec ] [ sht-spec ] [ trj-spec ] [ hsa-spec ]

   mod-spec ::= [ ISVFrequency 0 ]

   tps-spec ::= [ NTPAth      0 ] [ NSAVP        0 ] [ NPRAccept         0 ] -
                [ ITPRint     0 ] [ ITPUnit STDOUT ] [ ACCU         STDOUT ] -
                [ USER        0 ] [ PSHOot     1.0 ] [ IMXShift          1 ] 
                [ SDUNit      0 ] [ SDINit       0 ] 

   sht-spec ::= [ VFRAction 0.0 ] [ TFRAction  1.0 ] [ NTFRaction   NTPAth ] -
                [ IFSHoot    -1 ] [ IRST         0 ] [ PHALf             0 ] -
                [ ISLO        0 ] [ ISLN         0 ]  

   trj-spec ::= [ NUNIt       1 ] [ IFIRst      -1 ] [ VFIRst IFIRst+NUNIt ] -
                [ BEGIn       0 ] [ SKIP        -1 ] [ STOP              0 ] -
                [ PNSAve   NSTEP / NSAVC + 1]

   hsa-spec ::= [ HSAM        0 ] [ IHUN    STDOUT ] [ IHFR              0 ] -
                [ NHSV       -1 ] [ IHPR         1 ] [ NHST              0 ] 
             

.. _tps_description:

RXNCoor                              
-------

The RXNCoor command has been modified to facilitate the definition of the
basins that constrain the endpoints of the path.  

Multiple order parameters can be specified using the NRXN keyword
following the SET keyword (see umbrel.doc for a description of the
latter).  However, note that the tree structure associated with the
RXNCoor command is not dynamically allocated, so that MAXNOD in
rxncom.fcm must be increased with for larger numbers of order parameters.

The BASIn keyword is used to specify the boundaries of the basins (A
and B).  For each order parameter, one must give the name followed by
four real numbers: alo, ahi, blo, and bhi (the order is important).
Here, alo (blo) is the lower bound for basin A (B) and ahi (bhi) is
the upper bound for basin A (B).

The TPUNit keyword designates units to which values of the order
parameters are written.  Values are written for each saved trajectory.
A unit number less than 1 suppresses writeout.

In addition, a steered molecular dynamics (SMD) has been implemented 
through RXNCOR and integrated with transition path sampling (TPS) to allow 
generation of reactive trajectories de novo.  In this bias-annealing method, 
SMD allows transitions between basins (defined in the RXNCOR module) of the 
reaction coordinates where symmetric (harmonic) biasing potentials are 
advanced in a ratchet-like manner.  TPS shooting moves can then be made with 
progressively milder biasing potentials until unbiased reactive trajectories 
could be obtained.  This facility is activated by including the SMDDel keyword
(along with a non-zero KUMB value) in the UMBRella subcommand.  This facility
can be used with DYNAmics without TPS so long as the basins are defined.


DYNAmics
--------

The DYNAmics command has been modified to provide a looping structure to
calculate multiple trajectories (paths).  TPS is invoked by including the PATH
keyword.

The additional keywords that are specific to TPS are:

=========== ===================================================================
NTPAth      The number of paths to calculate.

NSAVP       The frequency of saving paths to the trajectory and velocity 
            files.

USER        Whether to user USERSB to calculated the order parameters
            used to determine the stable states. 

ITPRint     The frequency with which to write order parameter values
            at basin evaluations.

ITPUnit     The unit number on which to write order parameter values 
            at basin evaluations.

SDUNit      The unit number to which to write  random number generator 
            seeds for TPS with Langevin dynamics.

SDINit      The unit number from which to read random number generator 
            seeds for TPS with Langevin dynamics.

NPRAccept   The frequency of printing acceptance statistics.

ACCUnit     The unit number to which to write acceptance statistics
            as a function of saved structure.  The columns are the same
            as those printed by TPSACC to the CHARMM output file.

PSHOot      The fraction of moves that are shooting moves.

IMXShift    The maximum number of saved phase space points by which to 
            reptate the path in a shifting move.   In other words, a 
            shift can be up to IMXS*NSAVC molecular dynamics steps long.

VFRAction   The amount to perturb the velocities in shooting moves.  A 
            random vector is chosen from a Gaussian (Maxwell-Boltzmann) 
            distribution and then scaled by VFRAction.  The scaled vector 
            is added to the current velocity vector and the result is 
            scaled to conserve kinetic energy after correcting for SHAKE 
            if necessary.  This procedure has the effect of rotating the 
            3N-dimensional velocity vector without changing its magnitude.

TFRAction   The amount by which to scale the kinetic energy in each 
            shooting move if annealing is desired.

NTFRaction  The number of ACCEPTED shooting moves in which to scale the 
            kinetic energy by TFRAction.  

IFSHoot     The first point from which to shoot in units of NSAVC.  A 
            value of -1 indicates that IFSHot is chosen randomly if RTRJ 
            is specified and it is set to the middle of the path
            if shooting from a structure.

IRST        The saved phase space point to save to the restart file.  It 
            is best if IRST is chosen to be close to the transition state.
            Otherwise, numerical errors can prevent one from regenerating 
            a valid path using the saved phase space point.

PHALf       The probability of shooting half a trajectory.  A stochastic
            element should be included if PHALf is greater than
            zero.  For example, see discussion about Langevin dynamics.

ISLO        The lowest saved structure from which to shoot.

ISLN        The number of saved structures from which to shoot (i.e.,
            the last saved structure shot from is ISLO + ISLN - 1).

RTRJ        If this keyword is present (in place of STARt or RESTart), 
            an entire trajectory is read at the beginning of a restart.  
            The keywords BEGIn, SKIP, and STOP have their usual meaning.
            
NUNIt       The number of trajectory files to read if RTRJ is specified.

IFIRst      The first trajectory file to read if RTRJ is specified.

VFIRst      The first velocity   file to read if RTRJ is specified.

PNSAve      Allow TPS to read shorter trajectories into longer ones in
            order to lengthen the allowed transition time. PNSAve specifies
            the number of phase points in the short trajectories (default
            is the number of phase points to be saved in the new trajectories).

HSAMple     If this keyword is present, paths are accepted if they 
            start in basin A and ever go through basin B.  Also, the 
            probability that the system is in basin B as a function of 
            time is calculated [<h_B(t)>].

NHSTart     The path at which to start calculating <h_B(t)>.

IHUNit      The unit number to which to write <h_B(t)>.

IHFRequency The frequency with which to write <h_B(t)>.

NHSV        The frequency of evaluating whether the path is in 
            basin B (h_B[x(t)]).  A value of -1 sets NHSV to NSAVC.

IHPRint     h_B[x(t)] is printed every IHPRint*NHSV steps.
=========== ===================================================================
     

In addition, note that the meaning of the ISVFrequency keyword is changed
during TPS.  It refers to the number of PATHS, not the the number of 
molecular dynamics steps, between writes to the restart file.  

Note that if PHALf is greater than 0, shooting moves are carried out in 
which the path is only updated in one direction.  In this case a stochastic 
element should be included in the integration, such as Langevin dynamics. 
It is possible to apply Langevin integration to only the periphery of the 
simulation using the RBUF keyword.  When using Langevin dynamics with TPS, 
the random number seed used to generate the random forces is recorded for 
every saved structure.  This is necessary to regenerate the appropriate 
displacement vectors from the coordinates and velocities in 2-step dynamics 
during a shooting move. When writing and reading trajectories, these seeds 
can be written/read as designated by the keywords SDUNit and SDINit. If 
seeds are not read in with a trajectory, seeds are generated randomly for 
the initial shooting move, which results in a different displacement vector 
than in the original structure. If SHAKe is used in conjunction with Langevin 
dynamics, the new displacement vectors have some velocity components along 
the constrained bond.  These components are zeroed, and the overall kinetic 
energy will be reduced for that step.  Also note that, if SHAKe is used, 
there is a small error in the regeneration of  displacement vectors for 
atoms that have different values of FBETA and are connected by a shaken bond.


.. _tps_references:

REFERENCES
----------

All studies that employ TPS in CHARMM should reference:
     
*  Hagan, M. F., Dinner, A. R., Chandler, D. and Chakraborty, A. K. (2003) 
   Atomistic understanding of kinetic pathways for single base-pair binding 
   and unbinding in DNA.  Proc. Natl. Acad. Sci. USA 100, 13922-13927.

In addition, studies that employ SMD based on RXNCOR should reference:

*  Hu, J., Ma, A. and Dinner, A. R. (2006) Bias annealing:  A method for obtaining
   transition paths de novo.  J. Chem. Phys., submittted.

Additional references on TPS:

*  Dellago, C., Bolhuis, P., Csajka, F. and Chandler, D. (1998)
   Transition Path Sampling and the Calculation of Rate Constants.
   J. Chem. Phys . 108, 1964.

*  Dellago, C., Bolhuis, P. and Chandler, D. (1998) Efficient
   Transition Path Sampling:  Application to Lennard-Jones Cluster
   Rearrangements.  J. Chem. Phys. 108, 9236.

*  Dellago C., Bolhuis, P. G., Geissler, P. L. (2002) Transition path
   sampling. Adv. Chem. Phys. 123, 1.

*  Bolhuis, P. G., Chandler, D., Dellago, C. and Geissler, P. (2002)
   Transition Path Sampling: Throwing ropes over mountain passes, in
   the dark.  Ann. Rev. Phys. Chem. 59, 291.
