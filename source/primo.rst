.. py:module:: primo

======================
PRIMO Module of CHARMM
======================

Recently developed PRIMO model(proteins2010, 78:1266) consists of interaction
sites made up of one to three atoms. The recontruction from PRIMO to atomistic
phase is based on analytic relations based on the standard bonding geometries
observed in atomistic models.

Since PRIMO energy surface is smoother compared to the rugged one of atomistic
phase, coarse-grained conformations are tolerant of small local[bond/angle]
deformations which are forbidden in atomistic phase. This is a serious issue,
since the analytic reconstruction would not guarantee good energetic structure
in atomisitic phase. To alleviate this issue, a serious of virtual atoms are
constructed and used to restrict the sampling of coarse-grained particles.
Three types of virtual atoms(vs1,vs2,vs3) are supported as sketched below.


::

        vs1                   [4]
          .                  /
           .            q2  /
        b1  .           t2 /
             .            /
             [1]--------[3]
      t1    /             .
      q1   /               . b2
          /                 .
         /                   .
       [2]                   vs2-------[5].......vs3


Here [1],[2],[3],[4],[5] are the coarse-grained(CG) sites. Virtual site vs1 is
constructed using the bond distance vs1-[1], angle vs1-[1]-[2] and dihedral
vs1-[1]-[2]-[3]. Virtual site vs2 is constructed using bond vs2-[3], angle
vs2-[3]-[1] and dihedral vs2-[3]-[4]-vs1. Virtual site vs3 is estimated as
2*[5]-vs2.

Harmonic distance and angle potential involving virtual sites are supported.

.. _primo_syntax:

Syntax of the PRIMO Command
---------------------------

::

  PRIMO  RESN residue-selection [CSTEP]
         {  VS1 vs1-build-specs vs1-potential-specs
           [VS2 vs2-build-specs vs2-potential-specs]
           [VS3 vs3-build-specs vs3-potential-specs]
         }

  vs1-build-specs :: = { B1 real T1 real Q1 real ATM1 atom-selection
                         ATM2 atom-selection ATM3 atom-selection }

  vs2-build-specs ::= { vs1-build-specs ATM4 atom-selection B2 real T2 real
                        Q2 real}

  vs3-build-specs ::= { vs1-build-specs vs2-build-specs ATM5 atom-selection }

  vs1-potential-specs ::= { DIS1 DSEL1 atom-selection KDT1 real MND1 real
                           [ANG1 ASEL1 atom-selection ASEL2 atom-selection
                            KTH1 real MNT1 real] }

  vs2-potential-specs ::= {[DIS2 DSEL2 atom-selection KDT2 real MND2 real]
                           [ANG2 ASEL3 atom-selection ASEL4 atom-selection
                            KTH2 real MNT2 real] }

  vs2-potential-specs ::= {[DIS3 DSEL3 atom-selection KDT3 real MND3 real
                           [ANG3 ASEL5 atom-selection ASEL6 atom-selection
                            KTH3 real MNT3 real] }


.. _primo_function:

Function
--------

+-----+-------------------------------------------------------------------+
|RESN | constructs virtual atom for residue-selection. The residue        |
|     | selection should contain one residue(ALA,PHE..) and must be       |
|     | recognizable by RES(I) array of charmm module                     |
+-----+-------------------------------------------------------------------+
|CSTEP| evaluate potential during every CSTEP                             |
+-----+-------------------------------------------------------------------+
|VS1  | specifies virtual site vs1 construction                           |
+-----+-------------------------------------------------------------------+
|ATM1 | atom selections for building vs1, specified as CG particle name   |
|ATM2 | prefixed by '+', '-' or '' to indicate selection belongs to       |
|ATM3 | either next, previous or same residue as vs1 respectively.        |
+-----+-------------------------------------------------------------------+
|B1   | Equilibrium bond length: vs1 and ATM1                             |
+-----+-------------------------------------------------------------------+
|T1   | Equilibrium bond angle:  vs2, ATM1 and ATM2                       |
+-----+-------------------------------------------------------------------+
|Q1   | Equilibrium dihedral: vs3, ATM1, ATM3 and ATM4                    |
+-----+-------------------------------------------------------------------+
|VS2  | specifies virtual site vs2 construction                           |
+-----+-------------------------------------------------------------------+
|ATM4 | atom selection for building virtual site vs3 (in addition to      |
|     | ATM1,2,3..see above)                                              |
+-----+-------------------------------------------------------------------+
|B2   | Equilibrium bond length: vs2 and ATM3                             |
+-----+-------------------------------------------------------------------+
|T2   | Equilibrium bond angle:  vs2, ATM3 and ATM4                       |
+-----+-------------------------------------------------------------------+
|Q2   | Equilibrium dihedral: vs2, ATM3, ATM4 and virtual vs3             |
+-----+-------------------------------------------------------------------+
|VS3  | specifies virtual vs3 construction                                |
+-----+-------------------------------------------------------------------+
|ATM5 | atom selection for building CG-atom (in addition to ATM1,2,3,4)   |
+-----+-------------------------------------------------------------------+
|DIS1 | distance based harmonic potential between vs1 and CG particle     |
+-----+-------------------------------------------------------------------+
|DSEL1| atom selection for DIS1                                           |
+-----+-------------------------------------------------------------------+
|KDT1 | spring constant for DIS1                                          |
+-----+-------------------------------------------------------------------+
|MND1 | Equilibrium distance for DIS1                                     |
+-----+-------------------------------------------------------------------+
|ANG1 | angle based harmonic potential between vs1 and CG particles       |
+-----+-------------------------------------------------------------------+
|ASEL1| atom selections for ANG1 potential.                               |
|ASEL2|                                                                   |
+-----+-------------------------------------------------------------------+
|KTH1 | spring constant for ANG1                                          |
+-----+-------------------------------------------------------------------+
|MNT1 | Equilibrium angle for ANG1                                        |
+-----+-------------------------------------------------------------------+

The meaning DIS2, DSEL2, KDT2, MND2 and DIS3 DSEL3 KDT3 MD3 are similar to the
above definitions and are valid for vs2 and vs3 sites respecitively.

Similarly ANG2, ASEL3, ASEL4, KTH2, MNT2 are for vs2 angle potential and ANG3,
ASEL5, KTH3, MNT3 are for vs3 site respectively.


.. _primo_examples:

Examples
--------

* Example (1)

  Construct virtual C-atom for ala and restrain it with respect to PRIMO
  N-particle. This emulates N-CA-C angle term.

  ::

    primo resn ala vs1 b1 0.616 t1 41.50 q1 0.0 -
    atm1 CO atm2 +N  atm3 CA -
    dis1 dsel1 sele (resname ala .and. atom * * N) end kdt1 30.0 mnd1 2.472

* Example (2)

  Maintain the planarity of phe by applying an angle term between SC1,
  virtual-CG atom and SC2 particles. Apply the potential every 5th md step

  ::

    primo resn phe cstep 5 vs1 b1 0.616 t1 41.50 q1 0.0 -
               vs2 b2 1.552 t2 111.0 q2 -122.0 vs3 -
               atm1 CO atm2 +N atm3 CA atm4 N atm5 SC1 -
               ang3 asel5 sele (resname phe .and. atom * * SC1) end -
               asel6 sele (resname phe .and. atom * * SC2) end kth3 10.0 mnt3 150.0

* Example (3)

  The above task could also be accomplished by using a distance potential between
  virtual-CB and SC1/SC2 atoms of phe.

  ::

    primo resn phe cstep 5 vs1 b1 0.616 t1 41.50 q1 0.0 -
               vs2 b2 1.552 t2 111.0 q2 -122.0 atm1 CO atm2 +N atm3 CA atm4 N -
               dis2 dsel2 -
               sele (resname phe .and.  (atom * * SC2 .or. atom * * SC3) ) end -
               kdt2 30.0 mnd2 3.814

In all the cases the value of spring constant, must be carefully chosen as not
to adversely affect the CG conformational sampling.

Application of PRIMO module as residual PRIMO force field terms is given in
the CHARMM test case.

