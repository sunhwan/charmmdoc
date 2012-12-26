.. py:module::mts

=================================
Multiple Time Scales Method (MTS)
=================================

In CHARMM, multiple time scales method (MTS) algorithm is similar
to code of the algorithm described in the paper by Tuckerman, Berne,
and Martyna [J.C.P., 97, 1990 (1992)]. Please refer to this paper for 
details of derivations of this MTS-RESPA method. In addtion, more details 
can be seen in J. Chem. Phys. 99, 8063 (1993) and J. Phys. Chem., 99, 5680 
(1995) by  M. Watanabe and M. Karplus. In this new release, MTS method can
be called under parallel platforms. All modules under MTS should work in
parallel. To run CHARMM in parallel, please refer to parallel.doc.

The MTS method can be combined with Langevin dynamics via the
LN algorithm, described by Barth and Schlick [J.Chem.Phys., 1998, in press].
This version includes the slow forces via extrapolation and is expected to
allow larger timesteps than reversible MTS-RESPA. See
general notes at the end of this documentation file.
LN algorithm was implemented in CHARMM by Eric Barth (8/97) and
Adrian Sandu (7/98).    

In this documentation we refer to the rRESPA code as MTS-RESPA
(performing Newtonian dynamics) and to the LN code as MTS-LN
(performing Langevin  dynamics).     


.. _mts_syntax:

Syntax for the Multiple Time Scaled Method (MTS)
================================================

In this Multiple Time scaled method in CHARMM, a reversible RESPA
(Reference System Propagator Algorithm) is modified. (See Tuckerman's
paper about a reversible RESPA.) Two types of reversible RESPA methods
can be used in CHARMM.
 
1) Single reversible RESPA (Two-time-scale propagator)
 
   ::
   
     MTSmethod  1  I
     mts-spec
     END

      I       :: multiple time scales (Integer) - See Description
   
      mts-spec:: selection of hard forces or fast time scaled force
          [BOND]      ! Forces from all Bond-stretching motions  
          [ANGL]      ! Forces from all angle-bending motions
          [DIHE]      ! Forces from all diheral motions and all improper 
                      !   torsional motion
          [ALL ]      ! Forces from all internal motions which are
                      !   defined in CHAMM force fields
          [MASS] 1 M  ! Nonbonded forces involving atoms whose masses are
                      !   less than M. M is a mass weight.
          [SLFG]      ! Nonbonded forces separation by long- and short-
                      !   range forces  (see more detail in the text)
          [CLEA]      ! Clear MTS module and assignments

 
2) Double reversible RESPA (Three-time-scale propagator)

   ::
   
     MTSmethod  I  J
     mts-spec
     END

      I and J  :: multiple time scales (Integers) - See Description

      mts-spec:: selections of fast and medium time scaled forces
      
          [BOND]  K    ! Forces from all Bond-stretching motions  
          [ANGL]  K    ! Forces from all angle-bending motions
          [DIHE]  K    ! Forces from all diheral motions and all improper
                       !   torsional motion
          [ALL ]  K    ! Forces from all internal motions which are
                       !   defined in CHAMM force fields
          [MASS] K M   ! Nonbonded forces involving atoms whose masses are
                       !   less than M. M is a mass weight.
          [SLFG]       ! Nonbonded forces separation by long- and short-
                       !   range forces  (see more detail in the text)     
          [CLEA]       ! Clear MTS module and assignments

          K is 1 or 2 - 1 - force considered as a short time scaled
                        2 - force considered as a medium time scaled

   In both single and double RESPA methods, all interaciton forces, which
   are not selected by MTS command, are  considered as a long time scaled 
   degree of freedom

   You can see the more details in J. Chem. Phys. 99, 8063 (1993) and
   J. Phys. Chem., 99, 5680 (1995) by  M. Watanabe and M. Karplus.

3) MASS (Atomic Mass force separation)

   If MASS is selected, separating mass should be given.

   ::
   
      MTS
      .
      .
      MASS  K   M
      END

   K is the force considered as a short or medium time scaled. This
   should be 1 or 2 ( 1 for short and 2 for medium )
   M is the atomic mass. Separate the force contributions by the atomic
   mass. Contributions from the atom less than M mass is considered as the
   faster time scaled motions than those from the atom more than M mass.
   
   See the more details about this in J. Phys. Chem. 99, 5680 (1995),
   by M. Watanabe and M. Karplus.

4) SLFG (Short-Long Range force separation)

   If SLFG is selected, Short range force cutoff distance (RSCUT), 
   switching function healing length (RHEA), and Buffer healing 
   length (BUFF) should be given by the following way:

   ::
   
      MTS
      .
      SLFG   RSCUT [number (6.0)]  RHEA [number (1.0)]  BUFF [number (1.0)]
      END

      [number] is the distance in the unit of Angstrom.
      (..) is the default values.

   See the more details about these lengths in J. Phys. Chem. 98, 6885, 
   (1994) by D.D.Humphreys, R.A.Friesner, and B.J. Berne.

.. _mts_desc:

Description of MTS Dynamics Commands
====================================

MTS method approach is effective for special system where a separation
between the fast and slow time components is natural. The nature of
CHARMM force field allow us to separate some time scales. But in gener
there will be coupling between those motions, so this leads the
limitation of time scales.

a. In Multiple time scale, I and J are the number of cycle that you
   want to calculate short time scaled and medium time scaled motions,
   respectively, before calculating long time scaled motion.

   ::
   
      Delta t = J * Dtau2 = I * J * Dtau1 
   
   where Dtau1 and Dtau1 are the integral time step for short and
   medium time scaled motions respectively and Delta t is the integration
   time step of long time scaled motions. Dtau1 is defined in DYNAmic
   module as TIME.

b. MTS-RESPA method uses the velocity Verlet algorithm.
   MTS-LN algorithm solves the simple Langevin equation and
   relies on position Verlet and on constant extrapolation.
   For MTS-RESPA (Newtonian dynamics) use "DYNA VVER ..."
   For MTS-LN  (Langevin dynamics) define FBETA and use "DYNA LNX ..."      
   
c. Energy and forces from Urey-Bradely term is incoporated with bond
   command.

d. MTS-RESPA method is interacted with Nose-Hoover method. 
   In order to call Nose-
   Hoover method, you have to use Nose-Hoover module (see nose.doc and
   testcase for MTS method.)

e. If you are using with SHAKE, i.e. treat water molecules by SHAKE and
   other molecules without SHAKE, you MUST specify SHAKE BEFORE calling
   MTS command.

f. If you are going to use IMAGE module, you also have to specify the
   module before calling MTS module.

g. MASS selection only works with ATOM nonbonded selection. (See
   :doc:`nonbond`)

h. all NONBOND and UPDATE options listed before calling MTS module 
   are highly recommended.

i. If SLFG selection is used, update frequency of nonbond list may be
   specifically assigned in order to achieve the better energy 
   conservation of the system instead of using the automatic update
   frequency (INBFRQ -1).

j. ENERGY module has to be called before calling MTS or after calling
   dynamics. Otherwise, ENERGY cannot be calculated.
   

.. _mts_note:

Energy routine and selections of MTS method
===========================================

Selections of MTS method depend on energy routines and nonbond
options. Followings are lists of possible selections. (In the lists,
SHAKE means whether SHAKE method can be applied to the only environment, 
such as water molecules.)

1) GROUP selection

   ::
   
      S - Single reversible RESPA method
      D - Double reverisble RESPA method

      Scalar include fast and slow routine. 
      (SLFG only works with fast scalar routine.)

      X - acceptable selection
      ----------------------------------------------
      Selection | Vector    |  Scalar   |  SHAKE
                |  S  |  D  |  S  |  D  |  S  |  D
      ----------------------------------------------
       BOND     |  X  |     |  X  |  X  |  X  |  X
       ANGL     |  X  |     |  X  |  X  |  X  |  X
       DIHE     |  X  |     |  X  |  X  |  X  |  X
       ALL      |  X  |     |  X  |  X  |  X  |  X
       MASS     |     |     |     |     |     |
       SLFG     |     |     |  X  |  X  |  X  |  X 
      ----------------------------------------------
       MASS selection dosen't work with GROUP option.


2) ATOM selection 

   ::
   
      S - Single reversible RESPA method
      D - Double reverisble RESPA method

      X - acceptable selection

      ----------------------------------------------
      Selection | Vector    |  Scalar   | SHAKE
                |  S  |  D  |  S  |  D  |  S  |  D
      ----------------------------------------------
       BOND     |  X  |  X  |  X  |  X  |  X  |  X
       ANGL     |  X  |  X  |  X  |  X  |  X  |  X
       DIHE     |  X  |  X  |  X  |  X  |  X  |  X
       ALL      |  X  |  X  |  X  |  X  |  X  |  X
       MASS     |  X  |  X  |  X  |  X  |     |
       SLFG     |     |     |  X  |  X  |  X  |  X
      ----------------------------------------------

.. _mts_exam:

Examples of using MTS-RESPA method
==================================

The followings are examples of MTS method. Also you can check two
testcases for MTS method.

::

   ! Multiple Time Scaled Method Start
   ! Part I - Single Reversible RESPA method

   MTS 6   ! Integration time steps of 0.5fs for short time scale motions
           !   and 3.0fs for long time scales are used.
   BOND    ! Forces from all bond-stretching and bond bending motions
   ANGL    !   are treated as hard forces or short time scale degree of 
           !   freedom
   END

   DYNA VVER STRT NSTEP 1000 TIME 0.0005 -
        NPRINT 100   IPRFRQ  1000   -
        INBFRQ 20   IHBFRQ 0 FIRSTT 200.0 -
        IUNREA -30  IUNWRI 31  IUNCRD -32  IUNVEL -33 -
        KUNIT -34  IUNO -41 NSAVC 5  NSAVV 5  NSNOS 10 ISVFRQ 1000 

   ! Part II
   ! Double Reversible RESPA method
   ! Mass-scaling selection

   MTS 5 2      ! Integration time steps - 
                !     0.5fs for short time scale motions
                !     2.5fs for medium time scale motions and
                !     5.0fs for long time scale motions
   MASS 2  3.0  ! Nonbonded forces acting on atoms whose mass is less
                !   than 3.0g are treated as medium time scale.
   BOND 1       ! Forces from all bond-stretching and angle-bending
   ANGL 1       !   motions are considered as short time scale motions
   DIHE 2       ! Forces from all dihedral and imporper torsion motions
                !   are considered as medium time scale.
   END          ! Rest of force contributions are considered as long time
                !   scaled motions

   DYNA VVER REST NSTEP 1000 TIME 0.0005 - 
        NPRINT 100   IPRFRQ  1000   - 
        INBFRQ 20   IHBFRQ 0 -   
        IUNREA 30  IUNWRI 31  IUNCRD -32  IUNVEL -33 - 
        KUNIT -34  IUNO -41 NSAVC 5  NSAVV 5  NSNOS 10 ISVFRQ 1000 


   !
   ! PartIII
   ! Short-long range selection of MTS
   !
   MTS 2  3    ! Integrated time step
               ! 0.5fs for the fast time scale
               ! 1.0fs for the medium time scale
               ! 3.0fs for the slow time scale
   BOND  1     ! Forces from all bond-streching and angle-bending
   ANGLE 1     !    motions are consideered as fast time scale motions
   DIHE  2     ! Forces from all dihedral and improper torsion motions
               !    are considered as medium time scale
   SLFG  RSCUT  6.0   RHEA 2.0  BUFF 1.0
               ! Short-long range forces selection
               ! Short range forces are cut-off at 6.0A
   END


   DYNA VVER STRT NSTEP 200 TIME 0.0005 -
        NPRINT 10   IPRFRQ  1000   -
        FIRSTT 300.0  IUNREA -30  IUNWRI -31  IUNCRD -1  IUNVEL -1 -
        KUNIT -1  IUNO -1 NSAVC 5  NSAVV 5  NSNOS 10 ISVFRQ 1000  -
        TSTRUC 300 

   STOP

.. _mts_line:

Examples of using MTS-LN method
===============================

::

   ! PartIV
   ! LN algorithm - MTS with Langevin dynamics
   ! 
   ! Langevin coupling parameter needs to be defined
   SCALAR FBETA SET 20 SELE ALL END 

   ! Shake, if needed
   SHAKE BONH PARAM TOL 1.0D-9 SELE RESNAME TIP3 END

   MTS 4  24   ! Integrated time step
               ! 0.5  fs for the fast timestep
               ! 2.0  fs for the medium timestep
               ! 48.0 fs for the slow timestep
   BOND  1     ! Forces from all bond-streching and angle-bending
   ANGLE 1     !    motions are consideered as fast time scale motions
   DIHE  1     ! Forces from all dihedral and improper torsion motions
               !    are considered as medium time scale
   SLFG  RSCUT  6.0  RHEA 2.0  BUFF 1.0
               ! Distance class definition
               ! Short range forces are cut-off at 6.0A
   END


   DYNA LNX STRT NSTEP 200 TIME 0.0005 -
        NPRINT 10   IPRFRQ  1000   -
        FIRSTT 300.0  IUNREA -30  IUNWRI -31  IUNCRD -1  IUNVEL -1 -
        KUNIT -1  IUNO -1 NSAVC 5  NSAVV 5  NSNOS 10 ISVFRQ 1000  -
        TSTRUC 300 

   STOP

.. _mts_ln:

Parameter settings for LN
=========================

The algorithm relies on existing CHARMM force splitting routines under 
the MTS command. The  LN slow forces are incorporated via extrapolation as
opposed to "impulses" as in the MTS-RESPA method. This alleviates severe
resonance problems and permits larger outer timesteps to be used for
additional speedup.

The LN algorithm is compatible with SHAKE and with the use of boundary 
conditions. Any other combinations of options have not been tested.

There are several parameters that are set for a LN simulation:

1. Langevin parameter gamma (FBETA in CHARMM notation), the damping constant:

   * Recommended value = 5 to 20 ps^(-1)

   Too small a value will render the simulation unstable. On the other
   hand, the larger gamma is, the greater the overdamping of the 
   low-frequency modes. The above recommendation reflects a balance 
   found by experimentation. Gamma can also be simulation-goal dependent.
     
2. Timestep Protocol for force splitting:

   ::
   
      Dt(fast)   = inner TIMESTEP  for updating the "fast" forces

        * Recommended value = 0.5 -- 1 fs (no shake)
                               1  -- 2 fs (with shake)

      Dt(medium) = K1*Dt(fast), update  frequency for "medium" forces

        * Recommended value =  1  -- 3 fs            

      Dt(long)   = K2*Dt(medium), update frequency for "slow" forces

        * Recommended value =  6  -- 200 fs 
     
   Larger computational savings can be realized with a larger Dt(long).
   However, the speedup is limited and reaches an asymptotic value
   since the evaluation of medium forces becomes increasingly costly. 

   The asymptotic maximum speedup can be reached for outer timesteps of 
   24 or 48 fs, for example, but the precise value depends on the timestep 
   protocol employed and the application system. This should be tested 
   carefully by the user for the problem at hand.

3. Definition of the force splitting classes:

   Recommended Protocol --: 

   ::
   
     * Fast forces = BOND 1, ANGL 1, DIHE 1

     * Medium forces = Nonbond cutoff 
            cutoff distance = 6 A - 8 A
            healing region  = 1 A - 3 A
            buffer region   = 1 A - 3 A
       SLFG RSCUT [cutoff distance] RHEA [healing region]  BUFF [buffer region]

     * Longrange forces = remaining terms
   
   Nonboned pairlists are currently updated in the LN code
   every outer timestep; it is possible (but more costly)
   to attempt the updating every medium timestep.
    
   The nonbonded pairlist is updated in the current 
   implementation every Dt(long).

4. The GROUP electrostatics option works much better 
   than ATOM electrostatics.
   Use of the latter is discouraged based on our test problems.

.. note::
   All the LN parameters above can be sensitive to the 
   specific protocol used for the dynamics simulations
   and are problem dependent (see discussion of results
   in the LN papers).

   For further guidance, feel free to contact Tamar
   Schlick at the email:    schlick@nyu.edu


