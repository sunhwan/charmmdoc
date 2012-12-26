.. py:module:: pressure

============================================
Constant Pressure/Temperature (CPT) Dynamics
============================================

Two types of constant pressure/temperature dynamics are available
in CHARMM.  The weak coupling method for temperature and pressure
control described in the paper by Berendsen et al. (JCP 81(8) p3684
1984) was the first constant pressure and temperature algortihm
implemented in CHARMM.  Extended system constant pressure and
temperature algorithms have now been implemented based on the work
of Andersen (JCP 72(4) p2384 1980), Nose & Klein (Mol Physics 50(5)
p1055 1983), Hoover (Phys. Review A 31(3) p1695 1985).  Additionally,
a variant on the extended system method which treats the control
variables by means of a Langevin equation is available (Feller, Zhang,
Pastor & Brooks, JCP, 103, 4613 (1995)).

Shape matrix propagation and coordinates scaling for triclinic
unit cell is done according to D. Brown and J.H.R. Clarke in
Computer Physics Comm. 62 (1991) 360-369.

A constant surface tension algorithm is included which is useful
for studying interfacial systems where one wishes to allow the area
to change dynamically during the simulation.  The dynamical equations
and statistical ensemble are discussed in (Zhang, Feller, Brooks &
Pastor, JCP, 103, 10252 (1995)).

.. pressure_syntax:

Syntax
======

::

  [Syntax DYNAmics CPT]

  DYNAmics CPT ... cpt-spec

  cpt-spec::=  [ pressure-spec ] [ temperature-spec ] [ surface-tension-spec ]


  pressure-spec::= PCONST {[PINTernal]} {BEREndsen berensen-spec} ref-pressure-spec
                          { PEXTernal } { langevin-piston-spec  } [ IUPTEN int ]

  temperature-spec::= { TCONst [TCOUpling real] [TREFerence real] } ! Berendsen
                      {                                           }
                      { HOOVer    [TMASs real] [ REFT real]       } ! Hoover


  berensen-spec::= [COMPressibility real] [PCOUpling real]

  langevin-piston-spec::= piston-mass-spec [PGAMMA real] [TBATH real]

  surface-tension-spec::= [SURFace] [TENSion real]

  piston-mass-spec::= [PMASs real]
       [PMXX real] [PMYY real] [PMZZ real] [PMXY real] [PMXZ real] [PMYZ real]

  ref-pressure-spec::=  [PREFerence real] [PREFInitial real] [PREFFinal real]
     [PRXX real] [PRYY real] [PRZZ real] [PRXY real] [PRXZ real] [PRYZ real]
       [PIXX real] [PIYY real] [PIZZ real] [PIXY real] [PIXZ real] [PIYZ real]
         [PFXX real] [PFYY real] [PFZZ real] [PFXY real] [PFXZ real] [PFYZ real]
           [VOLUME real]


.. _pressure_description:

Description of CPT Dynamics Commands
====================================

Only a few changes are needed to a standard CHARMM dynamics input
file to run a CPT MD simulation. There are a few things to note :

a. The CPT algorithm is invoked by the CPT keyword.

b. It's not possible to use LANGEVIN dynamics with the constant
   pressure or temperature algorithm.

c. All the non-Langevin dynamics keywords have the same meaning as
   the in a standard dynamics input file. This includes the
   keywords STRT and REST.

e. The CPT specific keywords (apart from CPT itself) are :

   1. PCONstant    - do a constant pressure calculation.  Extended system
                     algorithm is the default, weak-coupling is available
                     with BEREndsen keyword.
   2. TCONstant    - do a constant temperature calculation with the weak-
                     coupling algorithm.  HOOVer constant temperature is
                     only available with PCONstant simulations.

f. The CPT module is only available for use with the leap-frog integrator

   3. To be used with Berendsen algorithm:

      ::

        COMPressibility <real> - the isothermal compressibility
                                 (atmospheres**-1).
        PCOUple         <real> - the pressure coupling constant
                                 (picoseconds).

      To be used with extended system algorithm:

      ::

        PMASs           <real> - the mass of the pressure piston (amu)
        PGAMma          <real> - Langevin piston collision frequency (1/ps)
        TBATh           <real> - Langevin piston bath temperature
        TENSion         <real> - reference surface tension (dyne/cm)
        IUPTEN          <int>  - unit number, P tensor at every step

        To be used with either constant pressure algorithm:
        PREFerence      <real> - the reference pressure (atmospheres).
                                 (for isotropic pressure)
        PREFInitial     <real> - initial reference pressure tensor (atmospheres).
                                 (for isotropic pressure)
        PREFFnitial     <real> - final reference pressure tensor (atmospheres).
                                 (for isotropic pressure)

        PRXX,PRYY,PRZZ  <real> - the reference pressure tensor (atmospheres).
        PRXY,PRXZ,PRYZ           (for anisotropic pressure)

        PIXX,PIYY,PIZZ  <real> - initial reference pressure tensor (atmospheres).
        PIXY,PIXZ,PIYZ           (for anisotropic pressure)

        PFXX,PFYY,PFZZ  <real> - final  reference pressure tensor (atmospheres).
        PFXY,PFXZ,PFYZ           (for anisotropic pressure)

        PREFI,PREFF,PIXX...,PFXX... - are used for linear pressure ramping


   4. To be used with Berendsen algorithm

      ::

        TCOUple         <real> - the temperature coupling constant
                                 (picoseconds).
        TREFerence      <real> - the berendsen reference temperature (K).

      To be used with extended system (HOOVer) algorithm

      ::

        TMASs           <real> - the mass of the thermal piston (kcal*mol^-1*ps^2).
        REFT            <real> - the hoover reference temperature (K).

   .. note:: for full descriptions of these parameters and the suggested
      values to use see the reference given above.

.. _pressure_notes:

Other Points
============

Suggested values for solvated systems:

::

      COMPressibility (beta) = 4.63e-5 /atm (for proteins)
      PCOUple                = 5.0 ps (or more)
      PREF                   = 1.0 atm (default)
      PMASs                  = 500 amus (default is infinity)
      TCOUPle                = 5.0 ps (or more)
      TREF                   = 298.0 K (default)
      REFT                   = 298.0 K (default)
      TBATh                  = 298.0 K (default)
      TMASs                  = 1000.0 kcal ps^2 (default is infinity)
      TENSion                = 0.0 (default, results in regular constant pressure)

Other Points
------------

1. Although the heating and equilibration commands are the same as
   for standard dynamics it is possible to use the CPT algorithm
   to perform both without velocity modification (c.f. Langevin
   dynamics).

2. The algorithm requires the use of the CHARMM CRYSTAL facility
   for constant pressure dynamics. If only Berendsen constant temperature
   is requested, then the crystal code need not be used.

3. For the Berendsen algorithm, when a reference pressure term
   ( PRXX,PRYY,PRZZ,PRXY,PRXZ,PRYZ is
   set to a very large negative number (less than -9999.0), then this
   component of the pressure will not be considered.  For orthorhombic
   simulations, the particular box length will not change (for example,
   if the command says: PRXX 1.0  PRYY -10000. PRZZ -10000.  then
   only the box length in the x direction will change during dynamics).
   If a cubic simulation is performed, then the pressure terms corresponding
   to the large negative values are not considered in the calculation of
   the instentaneous pressure.

4. For the extended system pressure algorithm, setting any component of the
   piston mass array to zero will cause the corresponding simulation cell
   length to remain constant.  For example when using the orthormobic
   cell, setting pmxx=pmyy=0 results in only the zdirection changing during
   the dynamics (and it changes according to the z component of the pressure
   tensor).  This is the standard method for interfacial NPAT systems.

5. A discussion of Hoover temperature control can be found in the documentation
   file nose.doc.  The temperature control implemented in the velocity verlet
   integrator is very similar to the one used in the leapfrog integrator.
   NOTE:  Hoover temp control only works in conjunction with constant pressure

6. The Berendsen pressure/temperature control scheme may not be appropriate
   for inhomegenous systems (protein in water, aqueous membrane, interfacial
   systems).  This is especially true if SHAKE constraints are used on one
   component.  A full discussion is given in the paper by Feller, Zhang, Pastor
   and Brooks (JCP, 9/15/95).

7. The extended system pressure algorithm can be run with temperature control
   (resulting in isothermal-isobaric ensemble) or without (resulting in
   isoenthalpic ensemble).  The Berendsen pressure method must be run with the
   constant temperature control (the ensemble for these methods is unknown).

8. If PINTernal is used (default), then the pressure is determined by the
   internal virial and the atoms in the box are instantaneously scaled in a
   homogeneous response to the altered box dimensions.  If PEXTernal is used,
   then the external virial (related to the force required to maintain the
   symmetry constraint) is used and the atom positions are not instantaneously
   scaled by box size changes.  The PEXTernal option is normally used for
   minimization, but may be used with molecular dynamics.  It is not recommended
   for large systems.

9. The pressure tensor (extended system) on every integration step is saved
   to the formatted file indicated by unit number specified by IUPTEN; this
   allows calculation of viscosity using the Green-Kubo relationship.  The
   viscosity is computed from the integral of the autocorrelation function
   of the off-diagonal elements of the P tensor, scaled by V/kT; see
   J. Phys. Chem. 100:17011-17020 (1996).  The column order of the data is:

   ::

     Time  PIXX  PIXY  PIXZ  PIYX  PIYY  PIYZ  PIZX  PIZY  PIZZ

.. _pressure_examples:

Examples
========

Examples of Constant Pressure Usage
-----------------------------------

1.  Basic constant pressure, appropriate for a cubic simulation cell, box
    length is coupled to the trace of the pressure tensor, using langevin on
    pressure piston degree of freedom.  This also works for a non cubic
    cell, but in that case each length moves independently to maintain a
    constant pressure tensor.  Constant volume is the limit pmass -> infinity
    (implementation in CHARMM: Set pmass = 0 for constant V).

    ::

      dynamics cpt leap restart time 0.001 nstep 10000 iseed 314159 -
              pconstant pmass 400.0 pref 1.0 pgamma 20.0 -
              tbath 300.0

2.  Constant normal pressure, constant area.  Appropriate for orthorombic
    cell where only the z direction is allowed to change.  The box
    length in z direction is coupled to the z component of the pressure tensor.

    ::

      dynamics cpt leap restart time 0.001 nstep 10000 iseed 314159 -
              pconstant pmzz 225.0 pmxx 0.0 pmyy 0.0 pref 1.0

3.  Constant pressure (stress) tensor.  Each box length moves independently
    to maintain the desired pressure tensor.

    ::

      dynamics cpt leap restart time 0.001 nstep 10000 iseed 314159 -
              pconstant pmass 225.0 przz 1.0 prxx 2.0 pryy 3.0

4.  Constant pressure, constant surface tension.  Z direction moves independently
    of x and y and is coupled to the bulk pressure (z component of pressure
    tensor).  X and y box lengths move to maintain constant surface tension.
    Note:  this is only appropriate for interfacial systems where the interface
    is perpendicular to the z axis.

    ::

      dynamics cpt leap restart time 0.001 nstep 10000 iseed 314159 -
              pconstant pmass 225.0 pref 1.0 surface tension 50.0

5.  Constant pressure and temperature (NPT)

    ::

      dynamics cpt leap restart time 0.001 nstep 10000 iseed 314159 -
              pconstant pmass 400.0 pref 1.0 pgamma 20.0 -
              tbath 300.0 tcons hoover reft 300. tmass 1000.

Examples of Constant Temperature Usage
--------------------------------------

1.  Basic constant temperature using the Hoover method

    ::

      dynamics cpt leap restart time 0.001 nstep 10000 iseed 314159 -
                 HOOVer   TMASs 1000.0  REFT 298.0

2.  Constant T with calculation of pressure data; this will also print
    the surface tension in the output log, useful for studying
    interfacial systems. The optional IUPTEN keyword will store the
    pressure tensor data for every timestep.

    ::

      open unit 29 card write name dyn.ptn
      dynamics cpt leap restart time 0.001 nstep 10000 iseed 314159 -
            pcons pmass 0.0 pint pref 1. iupten 29 -
            hoover tmass 1000.0 reft 293.0

.. _pressure_pressure:

The PRESsure command
====================

Process the pressure commands for the system. There are three
modes :

* Mode 1 : Initialise all pressure arrays.

  ::

     Syntax:

     PRESsure INITialise

* Mode 2 : Calculate and print the instantaneous pressures for
  a system.

  ::

     Syntax:

     PRESsure INSTantaneous TEMPerature <Real> VOLUme <Real> -
                            NDEGf <Integer> NOPRint

  The external isotropic pressure and tensor are calculated
  if a volume is present. The isotropic internal pressure is
  calculated if a volume is present and a temperature has been
  given. If no degrees of freedom are specified then a value
  of 3*NATOM is taken be default. The virials are always printed.
  NOPRint will suppress all printing.

  .. note: a previous call to energy is required so that the
     virials (and volume if not specified) are available
     in ENERGY.FCM. The command also accumulates the
     average and the square of the pressure variables.

* Mode 3 : Print the averages and fluctuations.

  ::

     PRESsure STATistics

