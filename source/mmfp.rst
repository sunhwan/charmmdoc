.. py:module::mmfp

=============================
Syntax of basic MMFP commands
=============================

.. index:: mmfp; syntax
.. _mmfp_syntax:


Syntax of basic MMFP commands
-----------------------------

::

   GEO reset

   GEO [MAXGEO integer] [shape_specification] [position_spec] [RCM] 
                 [potential_spec] [atom_selection] [ DISTANCE atom_selection]
                           [ ADISTANCE atom_selection atom_selection ] [PERP] 
                           [ ANGLE atom_selection atom_selection ]
                           [ DIHEDRAL 3 X atom_selection ]
                           [ IUMMFP unit]

   SSBP reset

   SSBP [atom_selection] [atom_selection] [ssbp_specification]

   BHEL [atom_selection] 

   SHEL [atom_selection] [shell_options-specification] 

   VMOD RESEt

   VMOD INIT MXMD integer UMDN integer [NATM integer] -
        KROT real KCGR real UOUT integer NSVQ integer [atom selection]

   { VMOD ADD IMDN integer KMDN real QN real }

   { VMOD CHANge NREStraint integer KMDN real QN real }


   shape_specification:==  { [SPHERE] } [XREF real] [YREF real] [ZREF real] -
                                        [TREF real]
                           { CYLINDER } [XDIR real] [YDIR real] [ZDIR real] 
                           { PLANAR   }  

   potential_spec:== { HARMonic } { INSIDE    } [FORCE real] -
                                                [DROFF real] [DTOFF real]
                     { QUARtic  } { OUTSIDE   } [P1 real] [P2 real]  
                     { EXPOnent } { SYMMETRIC }
                     { GAUSsian }
                     { SAWOod   }
                     

   ssbp_pecification:== KIRKWOOD NMULT [integer (15)]  [DIEC real] [RADI real] 
                          [DRDIE real] CAVITY HSR ANGU [EMP1 real] [EMP2 real]

   shell_options-specification:== DRSH [real (5.0)] RELA [real (0.0005)]
                                  FREF [real (0.9)]
                                  RSLV [real (0.0)] FOCO [real (3.0)]

   atom-selection:== (see *note select:(chmdoc/select.doc).)


.. _mmfp_details:

Details of basic MMFP commands
------------------------------

1) GEO RESET

   Cancels all restraints in GEO free all space allocated on the HEAP

2) GEO [MAXGEO  int]

   Allocate space on the HEAP to be used for all subsequent GEO potential terms.
   By default, MAXGEO is set to NATOM unless specified. The MMFP subroutine calls
   WRNDIE if there is not enough space allocated. 
   The keyword IUMMFP followed by a unit number will cause the position R,
   or the angle (degrees) for that constraint to be written to that file. 
   By default no writing is performed for each GEO restraint. 

3) RCM key word

   With the keyword RCM any restraints is not applied to each individual
   atoms of a selection but applied to the center of mass of the selected atoms.


4) [shape_specification]

   The shape of the potential is chosen from SPHERE, CYLINDER or PLANE key words,
   SPHERE is the default.  The shape specification gives the origin 
   (XREF, YREF, ZREF, or TREF for angles) and the orientation (XDIR, YDIR, ZDIR) 
   of a vector such that a sphere, plane or cylinder may be defined.  Using the 
   shape_specification the potential is calculated from the 
   general distance from a (x,y,z) reference point (SPHERE), distance from 
   an axis (CYLINDER) or distance from a plane (PLANE).  By default, all 
   values are zero and the origin of the boundary is at (0.0,0.0,0.0).
   If the shape of the boundary requires a unit vector (true for cylinder 
   and plane), and no values are given the subroutine will call WRNDIE.

5) [potential_spec]

   ::
   
      [potential_spec] [HARMonic]
                       [QUARtic]
                       [EXPOnential]
                       [GAUSsian]
                       [SAWOod]
                       
   The potential specification has a number of 
   parameters: 
   
   * [FORCE real] is the amplitude of the potential term
   * [P1 real] is a parameter used in the quartic, the gaussian,
     the exponetntial and the Saxon-Wood-type potential
   * [P2 real] parameter used in the Saxon-Wood-type potential 
   * [DROFF real] is an offset distance such that GEO(r) = 0 if r<droff
   * [DTOFF real] is an offset angle such that GEO(theta) = 0
     if theta<dtoff 
   * [INSIDE] the potential used only for r-droff<0
   * [OUTSIDE] the potential used is only for r-droff>0
   * [SYMMETRIC] the potential used is for \|r-droff\|
               
   They determine which kind of potential function will be used in combination
   with the geometrical shape.  The default is a harmonic potential.  A fourth
   order polynomial can be used with the key word QUARTIC, the potential has
   the form: ``GEO(r) = FORC*DELTA**2*(DELTA**2-P1)``, with DELTA=(R-DROFF).
   Using the parameters [FORCE 0.2 P1 2.25] the QUARTIC potential can be used
   to setup a spherical boundary potential with a well depth of -0.25 kcal/mol
   at r=DROFF+1 followed by a smoothly rising repulsion. Such potential is
   appropriate for a water sphere of radius DROFF+1.5 and is  very similar
   to that used in SBOUND, see :doc:`sbound`.

   The key word EXPO defines a exponential potential to mimic interfacial
   solvation effects:

   ::
   
           = HALF*FORC*EXP(-DELTA/P1),       for r > DROFF
           = FORC*(1 - HALF*EXP(+DELTA/P1),  for r < DROFF

   When defined in combination with PLANE shape_specification, this potential
   reproduces the "hydrophobic" potential used for  transmembrane polypeptide
   by O. Edholm.  and F. Jahnig, Biophys. Chem. 30, 279-292 (1988).
   The key word GAUSS defines a similar gaussian potential to mimic interfacial
   solvation effects. The parameter P1 gives the width of the interface.

   The keyword SAWO defines an exponential Saxon-Wood-type flat-bottom 
   potential of the form:

   ::
   
           = FORC/( 1 + Exp((P2-DELTA)/P1) ) - V(0)    for r > DROFF
           = FORC/( 1 + Exp((P2+DELTA)/P1) ) - V(0)    for r < DROFF

   where P1 is responsible for the steepness of the potential and P2
   determines the width (the distance between the two inflection points)
   of the restraint. V(0) is an offset correction to ensure a value of
   zero at the equilibrium point.
   
   This restraint should be helpful e.g., for binding free energy 
   difference calculations (it doesn't perturb the potential energy 
   landscape of the system within an adjustable range).

6) DISTANCE key word 

   With the keyword DISTANCE a restraint is setup
   between two sets of atoms or between their center of mass if the key
   word RCM is used.  A second atom selection must be specified.

7) ADISTANCE key word 

   With the keyword ADIS a restraint is setup
   between one atom set, and two other sets of atoms, such that the position
   of the first selection is constrained at some distance parallel to 
   the axis joining the centres of mass of the second and third atom selections.  
   A second and third atom selection must be specified. The keyword PERP
   will instead constrain the first atom selection at a distance 
   perpendicular to the axis vector. 

9) ANGLE keyword 

   With the keyword ANGLE a restraint is setup between 3
   sets of atoms or their center of masses if the keyword RCM is
   used. Three sets of atom selections must be made, note that the force
   constant is per radian**2 and NOT per degree**2 even though the TREF 
   (theta-reference, equivalent to DROFF of v29) 
   variable (angle constraint) is to be specified in degrees. 
   Specification of DTOFF variable can allow shifting of the potential 
   away from TREF, as is useful in the INSIde restraint. 

10) DIHEDRAL keyword

    With the keyword DIHEDRAL a restraint is setup
    between 4 sets of atoms or their center of masses if the keyword RCM
    is used. Four sets of atom selections must be made, note that the
    force constant is per radian**2 and NOT per degree**2 even though the
    TREF (equivalent to DROFF) variable (dihedral constraint) is to be
    specified in degrees.
    An offset of DTOFF may also be used for this restraint. 

11) SSBP key word

    Stands for Spherical Solvent Boundary Potential.  Current implementation of
    the method described in Beglov & Roux, J. Chem. Phys., 100:9050 (1994).
    The method follows from a rigorous reduction of the multi-dimensional
    configuration integral from N solvent molecules (10**23) to "n" solvent
    molecules (e.g., 1 to 1000).

    The SSBP potential corresponds to a constant temperature and constant
    pressure system.  The non-bonded interactions must be treated with EXTENDED
    electrostatics otherwise the system is unstable.  There are several
    contributions to the boundary potential of mean force:  HSR (hard sphere
    restriction) is a term setting the external pressure and surface tension;
    CAVITY ressembles to the standard stochastic boundary potential and
    corresponds to the van der Waals interactions; KIRKWOOD is the multipolar
    expansion for the reaction field due to a dielectric continuum surrounding a
    cavity containing a charge distribution;  ANGU is an angular correction that
    works for three sites water models and is used to restore the isotropic
    angular distribution near the edge of the sphere.  EMP1 and EMP2
    are two parameters for empirical gaussian potential (Deng, Y and Roux B.
    J. Phys. Chem. B, 108 (42), 16567--16576).
    The magnitude of the gaussian is controlled by EMP1,
    which has a default value of 1.1 kcal/mol.
    
    The width of the gaussian potential is controlled by EMP2, which
    has a default value of 0.008 angstrom^-2.  The empirical correction
    reduces the pressure in the simulation sphere, which is essential
    for correct free energy simulations.  The variable radius of
    the sphere is calculated on the fly and does not need to be specified. 
    The first atom selection flags the atoms for which the VDW and the ANGU 
    potentials are applied.  It also determines the radius of the boundary sphere.
    The second selection is optional.  If present it flags those atoms that 
    determine the radius of the boundary sphere.  By default, only the first
    flags everything; the second selection is there if one wants to remove
    some part of the system to determine the radius of the boundary sphere
    (such as a large part of a protein in an active site simulation).
    For bulk water sphere simulations, the first atom selection for should  
    be "select type OH2 end".  The second atoms selection is optional and 
    could be "select type OH2 end" or could be "select (.not. type H*) end".
    In NO CASE should the second selection includes the water hydrogens, since
    the results were NOT parametrized for this selection. 

12) BHEL  key word
 
    Stands for defining the boundary of the primary shell model as described in 
    Beglov & Roux, Biopolymers 35: 171-178 (1995).  This method is useful
    to provide one layer of solvent around a flexible polypeptide.
    The selection should be that of the protein or peptide heavy atoms only. 

13) SHEL  key word
 
    Stands for defining the solvent heavy atoms for the primary shell model.
    Other options allow to modify the effective force reference (analogous to
    the pressure (FREF).

14) VMOD  key word

    The VMOD facility  (David Perahia, Sylvain Frederic & Charles H. Robert
    2002-2008) is used to add one or more terms to the potential energy, each
    corresponding to a restraint to a given mrms projection on a normal
    mode or other 3N-dimensional vector. An appropriate reference is Floquet
    et al. (2006) FEBS Lett. 580, 5130-6. This facility is compatible with
    parallel operation.

    The command has several forms: initialization, adding a specification
    of a mode restraint, changing the restraint parameters, printing data
    concerning the current structrue, and resetting (to free the heap).

    VMOD INIT performs the initialization of the VMOD facility:
    
    * MXMD maximum number of mode restraints to add
    * KROT harmonic force constant for rotational restraint of the system
    * KCGR harmonic force constant for translational restraint of the system
    * UMDN the unit number of the open modes file. The keyword CARD can be used
      to specify a card-formatted modefile, the default is a binary file
    * NATM number of atoms in the modefile (defaults to number in current PSF)
    * UOUT unit number (formatted output) to write normal mode mrms coordinates
      and restraint energies at a given step
    * NSVQ frequency in terms of minimization or dynamics steps for writing
      detailed data to UOUT
      
    An optional atom selection permits restricting the restraint force to the
    desired subset of atoms present in the mode file.
    
    Note: The reference structure (e.g., structure for which the modes were
    calculated) must be in the main coords when invoking this command!

    VMOD ADD restraint statement(s) must (each) specify the following
    
    * IMDN is the mode number in the mode file
    * KMDN is the harmonic force constant (kcal/mol-A) for the mrms restraint
      QN is the desired target mrms value

    VMOD CHANge restraint statement(s) must (each) specify the following
    
    * NREStraint is the constraint (not mode) number (i.e., 1...MXMD)
    * KMDN is the new harmonic force constant (defaults to current value)
    * QN is the new desired target mrms value (defaults to current value)

    VMOD PRINt summarizes the mode projections and energies for the current structure

    VMOD RESET removes all existing VMOD restraints. It will give an error
    unless a VMOD INIT command has already been executed.

    In minimization or dynamics runs, the total restraint energy (Emode+Etrans+Erotat)
    is reported in the "MINI MMFP2>" or "DYNA MMFP2>" output, while more detailed
    data is written to the UOUT file at the desired frequency as specified in the
    VMOD INIT statement.


.. index:: mmfp; examples
.. _mmfp_examples:

Examples of MMFP GEO subcommnads
--------------------------------


1) To setup a harmonic spherical restraint on all oxygens around the origin
   (by default is harmonic potential and a sphere centered at the origin)

   ::
   
      MMFP
      GEO  force 100.0 select type O* end
      END

      The entirely equivalent detailed command would be
      MMFP
      GEO  sphere harm xref 0.0 yref 0.0 zref 0.0 force 100.0 select type O* end
      END

2) The spherical quartic potential is very similarly to SBOUND potential
   (Suitable for a sphere of radius of 13.0 angstroms centered at the origin)

   ::
   
      MMFP
      GEO  sphere quartic -
           force 0.2 droff 13.0 p1 2.25 select type OH2 end
      END

3) To impose a harmonic restraint on the center of mass of carbon alpha around
   (x,y,z) = (1.0,2.0,3.0)

   ::
   
      MMFP
      GEO  sphere  RCM -
           xref 1.0 yref 2.0 zref 3.0 -
           force 10.0 droff 0.0 select type CA end 
      END

4) To apply a harmonic cylindrical tube constraint of 8 angstroms radius, 
   the axis of the cylinder is directed along ydir 1.0 and passes through the 
   point: xref=4.0,yref=5.0,z=6.0)

   ::
   
      MMFP
      GEO  cylinder -
           xref 4.0 yref 5.0 zref 6.0 xdir 0.0 ydir 1.0 -
           force 100.0 droff 8.0 select type CA end
      END

5) To apply a planar harmonic constraint with normal in zdir 1.0

   ::
   
      MMFP
      GEO  plane -
           xref 7.0 yref 8.0 zref 9.0 zdir 1.0 -
           force 100.0 droff 0.0 select type N* end
      END

6) To fix the distance between the center of mass of two subset of atoms
   (e.g., two domains of a protein, two amino acids, etc...)

   ::
   
      MMFP
      GEO  sphere  RCM  distance -
           harmonic symmetric force 10.0 droff 5.0 -
           select bynu 1:10 end    select bynu 11:20 end
      END

7) To constrain the distance along an axis vector joining the center of mass 
   of two subset of atoms
   (e.g., and ion between two domains of a protein, two amino acids, etc...)

   ::
   
      MMFP
      GEO  ADIStance  sphere  RCM SELE RESName POT end -
           harmonic symmetric force 10.0 droff 5.0 -
           select bynu 1:10 end    select bynu 11:20 end
      END

8) To constrain the angle between the center of mass of 3 subset of atoms
   (e.g., 3 domains of a protein, 3 amino acids, etc...)

   ::
   
      MMFP
      GEO  sphere  RCM  angle -
           harmonic symmetric force 1000.0 tref 5.0 dtoff 0.0 -
           select bynu 1:10 end    select bynu 11:20 end   select bynu 21:30 end
      END

   Thus, the TREF variable specifies the reference angle value while the DTOFF
   variable specifies the offset to be used if necessary.
   The previous implementation (till c30 version) used DROFF to specify reference 
   angle/dihedral value with no provision for specifying flat-bottom harmonic 
   potential with an offset, the previous command is still valid but is not
   recommended.

9) To constrain the dihedral angle between the center of mass of 4 subset of
   atoms (e.g., 4 domains of a protein, 4 amino acids, etc...)

   ::
   
      MMFP
      GEO  sphere  RCM  dihedral -
           harmonic symmetric force 1000.0 tref 5.0 dtoff 0.0  -
           select bynu 1:10 end    select bynu 11:20 end -
           select bynu 21:30 end   select bynu 31:40 end
      END

10) In using VMOD to constrain the system to a given mrms value along a
    normal mode, the modes must have been calculated previously and the binary
    file opened for reading before entering the MMFP facility. Further, when
    the initialization command is given the current coordinates must be those
    of the minimum-energy configuration used for the mode calculation.

    To restrain the system to an mrms value of 0.1 along the first vibrational
    mode (mode 7), the following sequence is appropriate:

    ::
    
      MMFP
      VMOD INIT MXMD 1 KROT 1000 KCGR 1000 UMDN 10 UOUT 11 NSVQ 10
      VMOD IMDN 7 KMDN 100 QN 0.1
      END

    Force constants must of course be adapted to the problem at hand.

11) To reset all GEO potentials to zero and deallocate the HEAP space

    ::
    
      MMFP
      GEO reset
      END


.. index:: mmfp; substitutions
.. _mmfp_substitutions:
   
MMFP Substitution Parameters
----------------------------

There are several different variables that can be substituted in 
titles or CHARMM commands that are set by some of the MMFP commands 
(see :doc:`miscom`).  Here is a summary and description of each variable.

* :sub:`GEO`

   The total energy contribution of the GEO restraining potentials.

* :sub:`XCM`, :sub:`YCM`, :sub:`ZCM`

  The position of the center of mass of the last set of atom is returned.

* :sub:`XCM2`, :sub:`YCM2`, :sub:`ZCM2`

  The position of the center of mass of the second set of atoms is
  returned if the key word DISTANCE ADISTANCE or ANGLE or DIHEDRAL was issued.

* :sub:`XCM3`, :sub:`YCM3`, :sub:`ZCM3`

  The position of the center of mass of the third set of atoms is
  returned if the key word ADISTANCE, ANGLE or DIHEDRAL was issued.

* :sub:`XCM4`, :sub:`YCM4`, :sub:`ZCM4`

  The position of the center of mass of the fourth set of atoms is
  returned if the key word DIHEDRAL was issued.

* :sub:`RGEO`

  The distance/angle/dihedral used in the last potential calculation 
  is returned. Set if a MMFP constraint with the keyword DIST, ADIS or
  ANGLE or DIHEDRAL was used.

* :sub:`RADI`

  The instantaneous sphere radius for the SSBP method.

* :sub:`SSBPLRC`

  long-range free energy correction for SSBP.  Only set in PERT
  calculation with SSBP

* :sub:`SSBPLRCS`

  standard deviation of SSBP long-range correction.  Only set in PERT
  calculation with SSBP

Future developments
-------------------

1. The SSBP potential will be implemented for active site solvation (in
   which a large part of the protein lies outside the spherical region).  

2. A  primary shell model for the solvation of polypeptides will be
   implemented in the coming year.  For details, see Beglov & Roux, Biopol.
   (1995, in press).

   The method is used for providing a first shell of waters around a
   markedly non-spherical system.  The boundary potential is flexible and
   variable.  It adapts dynamically to the shape of the polypeptide during a
   dynamics.
