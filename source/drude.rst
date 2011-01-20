.. py:module:: drude

===========================
New Drude Oscillator Format
===========================

::

   by   Benoit Roux            (roux@uchicago.edu)
   and  Guillaume Lamoureux    (Guillaume.Lamoureux@umontreal.ca)
   and  Edward Harder          (eharder@uchicago.edu)
   and  Alex MacKerell Jr.     (alex@outerbanks.umaryland.edu)

As of July 2007 the old "Drude Oscillator Command" and "ANISOTROPY 
command" (see below) used to generate the Drude polarizable model will 
be replaced by a new format that encodes the Drude model from the RTF.  

In the present implementation the drude particles are
generated for all atoms for which polarizabilities (ie. via ALPHA) are
specified in the RTF file.  This allows for the drudes to be generated
automatically when a molecule is generated in CHARMM.  In addition,
code has been developed to allow for inclusion of atom-based Thole
scale factors, atom-based anisotropic polarizabilities and the
addition of lone pairs to selected atoms at the RTF level.  These
enhancements allow for all the information for the polarizable drude
force field to be included in the RTF and parameter files.

For each atom with a specified polarizability, a "Drude oscillator" 
is created by attaching to the atom an additional particle 
(using a fictitious chemical bond of length zero and of force constant 
'KDRUDE = k/2').  Each Drude particle is given a mass and a charge, 
taken from the mass and the charge of its atom (so that the total mass 
and charge are conserved for the "atom-Drude" pair).

As a whole, each "atom-Drude" pair has a charge 'Q', unchanged from
the partial charge the non-polarizable atom had prior to calling the
DRUDE command.  The "atom-Drude" pair forms a dipole 'q*d', where 'q'
is the charge on the Drude particle and 'd' is the displacement vector
going from the atom to its Drude particle.  Any external field 'E'
creates a net displacement 'd = q*E/k', and thus the "atom-Drude" pair
behaves as a point charge 'Q' with a polarizability 'alpha = q**2/k'.
The polarizabilities (in Angst**3) are read from WMAIN, and converted
into charges 'q', assuming a force constant 'k = 2*KDRUDE'.

See J. Chem. Phys. 119, 3025-3039 (2003) for more details.

In this implementation of the Drude oscillator, the force constant of the 
spring is a diagonal rank-2 tensor with components KDRUDE. This leads to 
an isotropic atomic polarizability, 'alpha = alpha_{11} = alpha_{22}
= alpha_{33} = q**2/k'. The ANISOTROPY command modifies the components of 
the Drude force constant tensor allowing for anisotropic atomic 
polarizabilites.

See J. Chem. Theory Comput. 2, 1587-1597 (2006) for more details.

The non-bonded lists are constructed such that, if the "real" atoms are 
in a 1-2, 1-3, or 1-4 relationship, their corresponding Drude particles
will also be in a 1-2, 1-3, or 1-4 relationship, respectively. 
 
The bonded interactions are modified to allow 1-2 and 1-3 screened
dipole-dipole interactions, as proposed by Thole [see Chem.Phys. 59,
341 (1981)].  If two atoms were not "seeing" each others in the
non-polarizable force field, their dipoles (and only their dipoles,
not their net partial charges) will "see" each others in the
polarizable force field.  The screening function strength is modulated 
by a parameter "a" that has been generalized to depend on the atom
type associated with the interacting pair of Drude oscillators 
(i.e. a = a_i + a_j).   These variable thole parameters "a_i" are 
encoded in the RTF through the THOLE command. 

.. _drude_description:

Description of new DRUDE model RTF
----------------------------------

1) Items in RTF file

   A) The AUTOGENERATE command now includes a flag for the drude
      particles.  This should be specified at the beginning of the RTF.
      
      ::

         AUTO ANGLES DIHE DRUDE

   B) MASS statements must be included for Drude particle types.
   
      ::

         MASS   104 DRUD      0.00000 DD ! drude particle

      The default drude type is DRUD.  However, alternate drude types may be
      specified.  These additional drude types must be i) specified with
      MASS statements, ii) the ATOMs to which they are applied must include
      a "TYPE drudetype" statement and iii) the NONBond parameters must be
      specified. For example, to include a special drude type on water the
      following would be needed in the RTF

      ::
      
         i)   MASS   153 DOH2      0.00000 DD ! water Drude
         ii)  ATOM   OH2  ODW      0.00000 ALPHA -0.97825258 TYPE DOH2
         iii) (NONBOND section)    D*       0.0   -0.0000    0.0000
                                              ! Wildcard for Drudes and dummy atom
         iv)  (via NBFIX) MAGD    DOH2    -0.01000  3.10000

      Term iv) is an NBFIX that allows for the interactions between
      magnesium and the water drude to be treated with a special LJ term.
      Indeed, the motivation for additional drude types are such special
      terms.

   C) ALPHA and THOLE parameters are set in the RTF after the charge.
      The presence of the "ALPHA" keyword flags that atom as being
      polarizable and a drude particle is assigned to it upon generation of
      the molecule.  An atom based Thole scale factor may also be specified
      via the "THOLE" keyword.  If the THOLE term is missing, a default
      value of 1.3 is assigned to that atom. This corresponds to a parameter 
      a = a_i + a_j = 2.6 which is the THOLE parameter in the old Drude 
      command syntax.

      ::
      
         ATOM C    C      0.630  ALPHA  -1.104  THOLE 1.073

   D) Specific drude types for a polarizable atom may be specified as
      follows using the "TYPE" keyword.  See 1A) above, sections i to iv for
      more details.

      ::
      
         ATOM C    C      0.630  ALPHA  -1.104  TYPE DC

   E) Lone pairs may now be set in the RTF using a syntax similar to the
      standard lonepair command (lonepair.doc).  All options for the
      generation of lonepairs (FIXed, CENTer, COLOcate etc.) as specified
      in the lonepair documentation may be used although the atom selection
      specification is simplied as shown in the following example.

      ::
      
         LONEPAIR relative LPA O  C  CL distance 0.3 angle 91.0  dihe 0.0
         LONEPAIR relative LPB O  C  N  distance 0.3 angle 91.0  dihe 0.0

   F) Anisotropic atomic polarizabilities.  Atom based anisotropic
      polarizabilities may be assigned to selected atoms via the ANISOTROPY
      syntax in the RTF.  the first atom selects the atom on which the
      polarizability is to be made anisotropic.  The second atom along with
      the first defines the "11" direction and the 3rd and 4th atom
      selections define the "22" direction, with the anisotropy defined by
      fractional polarizability components.  Accordingly, only the 11 and 22
      direction are needed since the trace of the A tensor is set to ONE. An
      example follows

      ::
      
         ANISOTROPY  O    C   CL     N  A11 0.697 A22 1.219

2) Items in the Parameter file

   A) The KDRUDE force constant for all the atoms is set by a wildcard in
      the bond parameters.  The term must be included for all drude types
      included in the model.  The wildcard terms can be overwritten by putting
      chemical type specific bond paramters following the wild card term.
      
      ::
      
         DRUD X     500.000     0.0000
         DRUD O     487.740     0.0000

   B) NONBOND parameters must be included for all drude types. In
      addition, NBFIX terms for the drudes may be included as needed (see 1A
      above).

3) Command Line Items

A) A "DRUDE" flag must included with the generate command to specify
   the inclusion of drudes on the molecule.  In addition, the drude mass
   is also set in this command using the DMASS keyword followed by the
   mass in AMU.  The default value for the drude mass is XXXX.

   generate NMA first none last none setup warn DRUDE DMASS 0.4


Drude Oscillator Command
------------------------

::

   by   Benoit Roux          (Benoit.Roux@med.cornell.edu)
   and  Guillaume Lamoureux  (Guillaume.Lamoureux@umontreal.ca)

The DRUDE command generates a polarizable system by modifying the
topology and parameters of an existing non-polarizable system.  For
each selected atom, it creates a "Drude oscillator" by attaching to
the atom an additional particle (using a fictitious chemical bond of
length zero and of force constant 'KDRUDE = k/2').  Each Drude
particle is given a mass and a charge, taken from the mass and the
charge of its atom (so that the total mass and charge are conserved
for the "atom-Drude" pair).

As a whole, each "atom-Drude" pair has a charge 'Q', unchanged from
the partial charge the non-polarizable atom had prior to calling the
DRUDE command.  The "atom-Drude" pair forms a dipole 'q*d', where 'q'
is the charge on the Drude particle and 'd' is the displacement vector
going from the atom to its Drude particle.  Any external field 'E'
creates a net displacement 'd = q*E/k', and thus the "atom-Drude" pair
behaves as a point charge 'Q' with a polarizability 'alpha = q**2/k'.
The polarizabilities (in Angst**3) are read from WMAIN, and converted
into charges 'q', assuming a force constant 'k = 2*KDRUDE'.

See J. Chem. Phys. 119, 3025-3039 (2003) for more details.

The bonded lists are modified so that, if the "real" atoms are in a
1-2, 1-3, or 1-4 relationship, their corresponding Drude particles
will also be in a 1-2, 1-3, or 1-4 relationship, respectively.  (This
is done by creating additional fictitious bonds of force constant zero
between the particles.)

For a single atom (charges in parenthesis):

::

           DRUDE           (q)
     A    -------->      A~DA
    (Q)              (Q-q)

For a diatomic molecule:

::

     A1     DRUDE        A1~DA1
     |     -------->     |
     A2                  A2~DA2

                         1-2 pairs: A1-A2, A1-DA2, DA1-A2, DA1-DA2

For a triatomic molecule:

::

     A1                   A1~DA1
       \       DRUDE        \
        A2    -------->      A2~DA2
       /                    /
     A3                   A3~DA3

                          1-2 pairs: A1-A2, A1-DA2, DA1-A2, DA1-DA2,
                                     A2-A3, A2-DA3, DA2-A3, DA2-DA3
                          1-3 pairs: A1-A3, A1-DA3, DA1-A3, DA1-DA3

The bonded interactions are modified to allow 1-2 and 1-3 screened
dipole-dipole interactions, as proposed by Thole [see Chem.Phys. 59,
341 (1981)].  If two atoms were not "seeing" each others in the
non-polarizable force field, their dipoles (and only their dipoles,
not their net partial charges) will "see" each others in the
polarizable force field.

In this implementation of the Drude oscillator, the force constant of the 
spring is a diagonal rank-2 tensor with components KDRUDE. This leads to 
an isotropic atomic polarizability, 'alpha = alpha_{11} = alpha_{22}
= alpha_{33} = q**2/k'. The ANISOTROPY command modifies the components of 
the Drude force constant tensor allowing for anisotropic atomic 
polarizabilites.


.. _drude_syntax:

Syntax of the DRUDE command
---------------------------

::

   DRUDe RESEt

   DRUDe [MASS real] [KDRUde real]  -
         [VTHOle logical] [THOLe real [SHAPe integer]]  -
         atom-selection


   atom-selection::= (see *note select:(chmdoc/select.doc).)


.. _drude_description:

Description of the DRUDE command
--------------------------------

   =======  =========  =================================================
   Keyword  Default    Purpose
   =======  =========  =================================================
   RESET               Desactivates the Drude particles.  After a reset,
                       the user should delete all Drude particles.

   MASS        0.0     Mass of the Drude oscillators (in amu).  The
                       default zero value causes the masses to be read
                       from the topology file.  Any nonzero value will
                       override the topology file.

   KDRUDE      0.0     Force constant of the atom-Drude bonds
                       (in kcal/mol/Angst**2).  The default zero value
                       causes the bond force constants to be read from
                       the parameter file.  Any nonzero value will
                       override the parameter file.

   VTHOLE              Uses variable thole parameters for dipole-dipole 
                       interactions. 
 
   THOLE       2.6     The screening factor for dipole-dipole interactions
                       between atoms excluded from the non-bonded
                       interactions.  To have no dipole-dipole
                       interactions between these bonded atoms, use
                       THOLE = 0.

   SHAPE       1       Specifies the shape of the dipole-dipole
                       screening function.
   =======  =========  =================================================

1) KDRUDE

   KDRUDE is the force constant (in kcal/mol/Angst**2) for the bond
   between each atom and its Drude particle the user wishes to use.  It
   is overriding any bond constant found in the parameter file.  For
   highly polar molecules like water, the recommended value for KDRUDE is
   500 kcal/mol/Angst**2.

   The atomic polarizabilities (in Angst**3) are read from the WMAIN
   array:
   
   ::
   
       alpha = abs(WMAIN)

   The charge on every Drude particle is computed using the following
   formula:
   
   ::
   
       q = sqrt( 2*KDRUDE * alpha / CCELEC ) * sign(WMAIN)

   The charges are given the signs of the WMAIN values.  As long as
   KDRUDE is large enough, the Drude particles will stay very close to
   their atoms, and the sign of 'q' is irrelevant.


2) THOLE, SHAPE

   (For 1-2 and 1-3 atoms only. All other interactions are regular
   Coulomb.)

   The THOLE parameter is a dimensionless factor that specifies the
   extent of the smearing of the charge 'q' on the Drude oscillators and
   of a contribution '-q' on the "real" atom.  The default value of THOLE
   is 2.6, that is, the 1-2 and 1-3 dipole-dipole interactions are turned
   on by default.  To turn the interactions off, set THOLE to zero.

   The default constant THOLE of 2.6 can be replaced by variable thole
   parameters using the VTHOLE keyword on the DRUDE command line.  
   The THOLE parameter between oscillators I and J is given by 
   THOLE = THOLEI + THOLEJ.  The THOLEI parameters for each atom are fit 
   along with the charges against ab initio data.  Values of THOLEI must 
   be given in the WCOMP array for all Drude oscillator containing atoms
   prior to the DRUDE command.  

   Because the dipole are explicitly made of two charges, the screened
   dipole-dipole interaction between two polarizable atoms (that is, two
   "atom-Drude" pairs) is actually the sum of the following four screened
   charge-charge interactions:
   
   ::
   
       ('q1' on Drude 1) - ('q2' on Drude 2)
       ('q1' on Drude 1) - ('-q2' on atom 2)
       ('-q1' on atom 1) - ('q2' on Drude 2)
       ('-q1' on atom 1) - ('-q2' on atom 2)

   The screened charge-charge interaction has the form:
   
   .. math::
      
      U(r_{12}) = \mathrm{CCELEC} * q_1 * q_2 * \frac{S(u_{12})}{r_{12}}
       
   where 'u12' is the normalized distance:
   
   .. math::
   
      u_{12} = r_{12} * \frac{ \mathrm{THOLE} }{ ( \alpha_1 * \alpha_2 )^\frac{1}{6} }

   'S' is a screening function defined by the SHAPE parameter:

      =====   =============================   ====================
      SHAPE   Screening function              Charge distributions
      =====   =============================   ====================
      1       S(u) = 1 - (1+u/2)*exp(-u)      Slater-Delta
      2       S(u) = erf(u)                   Gaussian
      =====   =============================   ====================

   The default value of SHAPE is 1, which is also the only shape
   currently implemented.  SHAPE=2 is reserved for Gaussian-Gaussian
   distributions.

   Two "atom-Drude" pairs have dipole-dipole interactions if the
   following conditions are met:
   
   1. The THOLE parameter is nonzero.
   2. In the non-polarizable force field, the two atoms where
      in the nonbonded exclusion list.

   To see if all the desired atoms have dipole-dipole interactions, use
   PRNLEV > 7.  Each call to the energy will print the atom numbers,
   polarizabilities and Drude particles's charges of each interacting
   pair.

   The energy from the dipole-dipole interactions is added to the ELEC
   energy term, and "SKIP ELEC" will skip the Thole interactions as well.


.. _drude_toppar:

Effect on the topology and parameters
-------------------------------------

1) New atoms

   The Drude particles are inserted immediately after their corresponding
   atoms.  For an atom type 'CA1', the DRUDE command will assign the atom
   type 'DCA1' for the Drude particle.  Since no regular atoms have names
   starting with a 'D', the Drude oscillators can be selected with
   "SELECT TYPE D* END".

2) Masses

   The masses for the selected atoms are modified so that the total mass
   of the atom-Drude pair corresponds to the atomic mass.  Try "SCALAR
   MASS SHOW" before and after calling the DRUDE command.

3) Charges

   The charges for the selected atoms are modified so that the total
   charge of the atom-Drude pair corresponds to the atomic partial
   charge.  Try "SCALAR CHARGE SHOW" before and after calling the DRUDE
   command.

4) Bonded interactions

   In addition to the atom-Drude bonds, zero force bonds are created to
   maintain between the atom-Drude pairs the same 1-2, 1-3, and 1-4
   relationships that were existing previously to the DRUDE call.  For
   two bonded atoms A1 and A2, with Drude particles DA1 and DA2

   ::
   
       DA1      DA2
        |        |
        A1 ----- A2

   zero force bonds are created between DA1 and DA2, between DA1 and A2,
   and between A1 and DA2, so that any particle of the A1-DA1 pair is 1-2
   bonded to any particle of the A2-DA2 pair.  Since the force constants
   of these fictitious bonds are zero, the computational overhead is
   minimal.

5) Non-bonded interactions

   Weither the Drude particles have Lennard-Jones parameters or not, the
   Lennard-Jones parameters of the selected atoms are kept unchanged.
   Since the polarizable force field is built from the same "ingredients"
   as the non-polarizable force field, all the NBONDS options can be used
   as before (notably the PMEWALD method).

.. _drude_warning:

To be aware of when using the DRUDE command

1) Call the DRUDE command after all the atoms are built

   Otherwise, some zero-force bonds between the Drude particles and
   neighboring atoms (as discussed in the previous section) may be
   missing.  And since bad contacts are difficult to resolve with a
   polarizable force field, it is probably safer to minimize/equilibrate
   the system using first a non-polarizable force field.

2) Call the DRUDE command before SHAKE and LONEPAIR

   The preferred call to SHAKE is:

   ::
   
      COOR COPY COMP
      SHAKE BONH PARAM TOLERANCE 10E-9 -
            NOFAST -
            SELECT .NOT. TYPE D* END -
            SELECT .NOT. TYPE D* END

2) Always delete all the Drude particles after a RESET

   Treat this sequence of commands a single command:

   ::
   
       DRUDE RESET
       DELETE ATOMS SELECT TYPE D* END

   The "DRUDE RESET" command puts back the mass and the charge of the
   Drude particles on the heavy atoms, and erases the distinction between
   a Drude particle and a regular atom.

3) Beware of atom names conflicts

   Since the atom names don't have more than four letters, atoms with
   different names may end up having Drude particles with the same
   name:

   ::
   
      C210 --> DC21
      C211 --> DC21

4) MASS and KDRUDE are overriding the toppar files

   Any nonzero value for MASS and KDRUDE specified by the user is
   overriding the values from the topology and parameter files.

.. _drude_examples:

Usage examples of the DRUDE command
-----------------------------------

1) Polarizable benzene
   (see test/c30test/drude_benzene.inp)

   After reading the standard, non-polarizable topology and parameter
   files, a standard benzene molecule is generated:

   ::
   
      READ SEQUENCE BENZ 1
      GENERATE BENZ SETUP FIRST NONE LAST NONE

   The polarizabilities on all carbon atoms are set to 1.5 Angst**3:

   ::
   
      SCALAR WMAIN SET +1.5 SELECT .NOT. TYPE H* END
      DRUDE SELECT .NOT. TYPE H* END

   The selection contains the atoms CG, CD1, CD2, CE1, CE2, and CZ.  The
   DRUDE command will look for atom types DCG, DCD1, DCD2, DCE1, DCE2,
   and DCZ.  If these types are unknown, the program will crash.  For
   this reason, it is necessary to append the atom types of the Drude
   particles when reading the topology:

   ::
   
      OPEN READ CARD UNIT 1 NAME @TOPPAR/top_all22_drude.inp
      READ RTF CARD APPEND UNIT 1

   Similarly, the program will crash if bond parameters are missing, and
   the additional bond parameters should be appended to the parameters:

   ::
   
      OPEN READ CARD UNIT 1 NAME @TOPPAR/par_all22_drude.inp
      READ PARAM CARD APPEND UNIT 1

   The structure is minimized:

   ::
   
      MINI SD STEP 0.001 NSTEP 1000 NPRINT 100

   Since benzene is a nonpolar molecule, the Drude particles are not
   significantly moved from their heavy atoms.  To find the induced
   atomic dipoles for a given structure, one should use ``CONS FIX SELECT
   .NOT. TYPE D* END`` before calling MINI.
   
   The molecular polarizability is obtained using the VIBRAN command with
   the DIPOLES keyword:
   
   ::
   
       VIBRAN
       DIAGONALIZE
       PRINT NORMAL VECTORS DIPOLES SELECT ALL END
       END
   
   The total polarizability is an anisotropic tensor similar to the
   experimental results for benzene [J.Chem.Phys. 95, 5873 (1991)].  The
   strong anisotropy comes from the 1-2 and 1-3 dipole-dipole
   interactions.  Desactivating these interactions by using THOLE=0, the
   polarizability tensor is almost isotropic.

2) SWM4-DP water
   (see test/c30test/swm4.inp)

   See J. Chem. Phys. 119, 5185-5197 (2003) and Chem. Phys. Lett. 418, 
   245-249 (2005) for a complete description of the model.
   After reading the topology and parameter files, the model is built as
   following:

   ::
     
      READ SEQUENCE SWM4 ...
      GENERATE WAT SETUP NOANGLE NODIHEDRAL
     
      READ COOR CARD ...
     
      SET ALPHAO = 1.042520
      SET DOM    = 0.238080
     
      SCALAR WMAIN SET @ALPHAO SELECT ( SEGID WAT .AND. TYPE OH2 ) END
      DRUDE SELECT ( SEGID WAT .AND. TYPE OH2 ) END
     
      COOR COPY COMP
      SHAKE BONH PARAM TOLERANCE 1.0E-9 -
            NOFAST -
            SELECT ( SEGID WAT .AND. .NOT. TYPE D* ) END -
            SELECT ( SEGID WAT .AND. .NOT. TYPE D* ) END
     
      LONEPAIR BISECTOR DIST @DOM ANGLE 0.0 DIHE 0.0 -
               SELECT ATOM WAT * OM END  SELECT ATOM WAT * OH2 END - 
               SELECT ATOM WAT * H1 END  SELECT ATOM WAT * H2  END
     

   The molecular dynamics for polarizable water is explained in
   :doc:`vv2`.

ANISOTROPY Command
------------------

::

   by   Edward Harder        (eharder@uchicago.edu)
   and  Benoit Roux          (roux@uchicago.edu)

In the above implementation of the Drude oscillator, the force constant 
of the spring is a diagonal rank-2 tensor with components KDRUDE. This 
leads to an isotropic atomic polarizability, 'alpha = alpha_{11} = 
alpha_{22} = alpha_{33} = q**2/k'. The ANISOTROPY command modifies the 
components of the Drude force constant tensor allowing for anisotropic 
atomic polarizabilites.

Syntax of the ANISOTROPY command
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

   ANISOTROPY [MASS real] [KDRUde real]  -
              [THOLe real [SHAPe integer]]  -
              atom-selection


   atom-selection::= (see *note select:(chmdoc/select.doc).)

Description of the ANISOTROPY command
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   =======  ========   ===============================================
   Keyword  Default    Purpose
   =======  ========   ===============================================
   K11,      KDRUDE    Components of the force constant tensor 
   K22,                (in kcal/mol/Angst**2).  The default value
   K33,                is the KDRUDE parameter assigned from the DRUDE      
                       command.

   VERBOSE             Prints atoms involved in selection and 
                       components of the force constant tensor.  
   =======  ========   ===============================================

1) K11, K22, K33

   K11, K22, K33 are the components of the force constant tensor
   (in kcal/mol/Angst**2).  The directions 1, 2 and 3 are defined in the 
   molecular reference frame using atom selections that follow the 
   assignment of these variables.  The first atom selection contains the 
   anisotropic Drude oscillator.  The vector connecting the first atom to 
   the second atom selection defines the 1 direction.  The vector between 
   the 3rd and 4th atom selections defines the 2 direction.  The 3 direction 
   is orthogonal to 1 and 2.  The charge on the Drude for anisotropic 
   oscillators is:

   ::
   
      q = sqrt( 2*K33 * alpha / CCELEC ) * sign(WMAIN)


To be aware of when using the ANISOTROPY command
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1) Call the ANISOTROPY command after the rest of the charge model 
   has been built (i.e. DRUDE and LONEPAIR)

   Otherwise, inconsitencies in the particle charges and intended 
   polarizabilites may result.

Usage examples of the ANISOTROPY command
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1) Polarizable NMA

   ::
   
      ANISOTROPY -
      K11 700 K22 400 K33 450 select segid NMA .and. type O show end -
                              select segid NMA .and. type C show end -
                              select segid NMA .and. type CL show end -
                              select segid NMA .and. type N show end  VERBOSE

   The atomic polarizability of the carbonyl oxygen in NMA is made 
   anisotropic.  The 1 direction is parallel to the carbonyl bond. 
   The 2 direction is parallel to the vector between the CL carbon
   and nitrogen.  The 3 direction is normal to the plane spanned by 
   1 and 2.

