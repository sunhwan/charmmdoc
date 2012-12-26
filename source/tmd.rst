.. py:module:: TMD

===========================
Targeted Molecular Dynamics
===========================

The Targeted Molecular Dynamics (TMD) method introduces a holonomic
constraint that reduces the rmsd with a predefined target at each MD
step. Three flavors of the method are available: the original TMD
method (J. Schlitter, M. Engels, P. Kruger, E. Jacoby and A. Wollmer,
Mol. Sim. 10, 291 (1993)), the zeta-TMD method (Q. Cui, to be
published), and the restricted perturbation - TMD method (A. van der
Vaart and M. Karplus, J. Chem. Phys. 122, 114903 (2005)).
The methods are implemented for the LEAP integrator, and Berendsen's
thermostat must be used; CHARMM needs to be compiled with the "TMD"
keyword present in the pref.dat file.


.. _tmd_syntax:

Syntax for the Targeted Molecular Dynamics commands
---------------------------------------------------

::

  TMDInitialize { [ INRT integer ]  }   [ atom-selection ]  [ atom-selection ]
                { [ DINC real ]     }
                { [ FRMS real ]     }
                { [ ITMD integer ]  }
                { [ FTMD integer ]  }
                { [ ENER integer ]  }
                { [ MAXF real ]     }
                { [ MAXB real ]     }
                { [ SUMP integer ]  }
                { [ ZETA ]          }
                { [ CZETa real ]    }
                { [ ZTOL real ]     }
                { [ ZMIT integer ]  }

  INRT    1000     Number of step that one needs to get rid of artificial
                   rotational motion in TMD simulation.

  DINC    0.0      RMS increment in TMD simulation (ignored in RP-TMD
                   simulation).

  FRMS    1.0D-6   Stop dynamics when a rmsd of FRMS with the target is reached.

  ITMD    -1       Write analysis file to this unit (-1: don't write).**

  FTMD    -1       Frequency of writing analysis file.**

  ENER     0       Write potential energy difference due to the TMD constraint
                   to the analysis file if >0 (very expensive!).**

  MAXF    -1.0     If positive: perform a RP-TMD simulation with the value of
                   MAXF as the maximum perturbation.
                   If negative: perform an "original" TMD simulation.**

  MAXB    -1.0     If MAXF<0: ignored.
                   If MAXF>0 and MAXB>0: Use MAXF as the maximum perturbation
                   when C>0, use MAXB as the maximum perturbation when C<0,
                   where C = Sum ( |P_i| Cos(p_i F_i)).
                   If MAXF>0 and MAXB<0: Always use MAXF as the maximum
                   perturbation.**

  SUMP     0       If MAXF<0: ignored.
                   If 0: the maximum perturbation is the maximum atomic
                   perturbation (MAXF = max |p_i|).
                   Else: the maximum perturbation is a sum over all
                   atoms (MAXF = Sum |p_i|).**

  ZETA             Indicates Zeta form of TMD is being used (default: not active)

  CZETa   1.0      Zeta form expotential factor.

  ZTOL    1.0E-10  Tolerance used in calculating Zeta-TMD constraint.
                   Positions are solved for (Zeta-Zeta0) < ZTOL.

  ZMIT    1000     Number of iterations allowed in the Zeta-TMD constraint
                   subroutine.

  1st atom-selection   Apply TMD fit to the selected atoms only.

  2nd atom-selection   Apply TMD perturbation to the selected atoms only.


  ** For TMD/RP-TMD only.

.. _tmd_description:

Description of the Targeted Molecular Dynamics Commands
-------------------------------------------------------

To invoke TMD, the TMDInitialize command should be given before the
DYNAmics command. After TMDInitialize, the target structure should be
read in with the "READ COOR TARG" command. All TMD variables and the
target coordinates are cleared after the DYNAmics command; before
a restart you should invoke the TMDInitialize command (and the "READ
COOR TARG" command) once again.

Two atom selections are used with the TMDInitialize command.  The first
selection is used to define the atoms used in fitting both targets
to the current structure (done every INRT steps).  The second selection
is used to define the atoms which the TMD constraint is applied.
If only one selection is given in TMDI, this selection will be used for
both fitting and applying the constraint.  To run TMD in its original
form, one must use 'select all end' for both selections. One should
perform an overlay of the structures before the simulation.

Please note that the "standard" TMD and RP-TMD methods are parallellized;
the 'Zeta-TMD' has not been parallellized.

A) Original TMD method.

   For doing Targeted Molecular Dynamics (TMD), one needs to define
   a moving coordinate and a target coordinate. You slowly pull the
   moving structure towards the target structure by gradually decreasing
   the RMS distance between two. The 'pulling' speed is defined
   by the user.

   Commands:

   ::

     OPEN UNIT 88 WRITE CARD NAME tmd.dat
     TMDINITIALIZE ITMD 88 FTMD 10 FRMS 1.2 INRT 10 DINCRE 0.0004 -
     SELE ALL END SELE ALL END

     OPEN READ UNIT 2 CARD NAME target.crd
     READ COOR UNIT 2 CARD TARG
     CLOSE UNIT 2

     DYNA RESTART LEAP TCONST TCOUPL 0.5 TREFER 300.0 ...


B) Zeta-TMD method.

   To constrain dynamics between two target structures (between a
   starting and ending structure of a conformational transition,
   for example), the 'Zeta' form of the constraint function is used:

   ::

       Zeta(t) - Zeta0(istep) = 0  (contraint),

   where

   ::

       Zeta(t)  = -1/(1+EXP(-CZETA*RMSD1(t))) + 1/(1+EXP(-CZETA*RMSD2(t)))
       Zeta0(istep) = Zeta0(istep-1) - DINC
       RMSD1 = (mass weighted) root-mean-squared difference (RMSD)
               between the current structure and TARG (using atoms
               defined by second selection in TMDI)
       RMSD2 = RMSD between the current structure and TAR2

   The two target structures are read using READ COOR TARG and READ COOR TAR2,
   respectively.  The sign of DINC (incrementation of the contraint function)
   determines which structure the molecule is pulled towards, and which one is
   pushed away from, during dynamics run. For DINC > 0, TMD pulls the molecule
   towards TAR2.  For DINC < 0, TMD pulls the molecule towards TARG.

   The starting value of Zeta0 is based on the coordinates at the start of
   the dynamics run; istep is the current step number in the dynamics run.
   The ZETA keyword must be used in the DYNA command line for this.  If two
   targets are read in, but ZETA is not specified, then only the one TARG
   structure is used in the TMD algorithm and the Zeta form of the
   constraint is not used.

   The Zeta form is useful, since it is more effective at pulling molecules
   towards target structures than other relative constraint forms, such as
   ((RMSD1 - RMSD2) - rho) = 0, where the difference in RMSDs may be well
   defined, but the current structure may be far from both target structures.
   Also, transitions are not limited to paths which only allow for the RMSD
   to one target structure to decrease monotomically.

   ZTOLerance and ZMITerations are used in the minimization scheme for
   calculating the coodinates which satisfy the TMD constraint. They are
   similar to the cooresponding terms in the SHAKE algorithm.

   The constrained RMSD for one-target TMD is not allowed to go below zero.
   Similarily, the restriction |Zeta| <= -1/(1+EXP(-CZETA*RMSD0))+1/2 is
   is used, where RMSD0 is the RMS Difference between the two target
   structures.  Once these values are reached during dynamics, the contraint
   value for RMSD (or Zeta) is held at this limiting value.

   This subroutine outputs RMSD1, RMSD2, and the actual Zeta value (which
   is within +/- ZTOL of Zeta0(istep)), with PRNLEV >= 5.  For each dynamics
   step, this is likely to print out a few times, due to the iterative scheme
   used between this subroutine and the SHAKE subroutine.

   Commands:

   ::

     TMDINITIALIZE INRT 1 DINC -0.0003  -
     ZETA  CZETA 1.0  ZTOL 1.0E-8  ZMIT 1000 -
     SELECT ALL END SELECT ALL END

     OPEN READ UNIT 2 CARD NAME target.crd
     READ COOR UNIT 2 CARD TARG
     CLOSE UNIT 2

     OPEN READ UNIT 2 CARD NAME init.crd
     READ COOR UNIT 2 CARD TAR2
     CLOSE UNIT 2

     DYNA RESTART LEAP TCONST TCOUPL 0.5 TREFER 300.0 ...


C) Restricted perturbation - TMD method.

   In this method, the coordinate displacement (perturbation) is limited
   to a preset value; given this displacement, the rmsd with the target
   is minimized at each step. This procedure prevents the crossing of
   large energy barriers, and may increase the efficiency of the calculation.
   Either the total perturbation Sum |p_i| or the maximum atomic
   perturbation Max |p_i| can be restricted.
   The function C = Sum ( |p_i| Cos(p_i,F_i) ) is a good indicator of
   barrier crossings: when C is negative, a barrier has probably been
   crossed. To reduce barrier crossings, the perturbation can be decreased
   when C<0 (this will increase the simulation time).
   Note that in this method, the rmsd fluctuates along the trajectory and
   the length of the simulation may vary (simulations may get "stuck" when
   very small perturbations are used).
   See J. Chem. Phys. 122, 114903 (2005) for a more detailed discussion
   of the algorithm, and a comparison of the RP-TMD method with the
   standard TMD method.

   Commands:

   ::

     OPEN UNIT 88 WRITE CARD NAME tmd.dat
     TMDINITIALIZE ITMD 88 FTMD 10 FRMS 1.2 INRT 10 MAXF 0.001 -
     MAXB 0.0008 SUMP 1 -
     SELE ALL END SELE ALL END

     OPEN READ UNIT 2 CARD NAME target.crd
     READ COOR UNIT 2 CARD TARG
     CLOSE UNIT 2

     DYNA RESTART LEAP TCONST TCOUPL 0.5 TREFER 300.0 ...

   Examples: tmdtest32.inp (serial & parallel), and tmd_zeta.inp (serial).



