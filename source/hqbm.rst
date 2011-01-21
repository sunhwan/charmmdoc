.. py:module:: hqbm

=========================
The HQBM Module of CHARMM
=========================

By Emanuele Paci, 1997/2000

HQBM is an external perturbation designed induce conformational
changes in macromolecules. The time dependent perturbation is designed
to introduce a very small perturbation to the short time dynamics of
the system and does not affect the conservation of the constants of
motion of the system (the conservation of the total energy or of the
suitable conserved quantity when an extended Lagrangian is used can
then be used as a check of the correctness of the forces).

The external perturbation needs:

- a reference (or target) structure
- a reaction coordinate which defines a "distance" from the 
  reference structure


.. _hqbm_syntax:

Syntax
------

- read the reference structure 

  ::
  
     OPEN UNIT 1 READ FORMATTED NAME coor0.crd
     READ COOR CARD COMP UNIT 1 
     CLOSE UNIT 1

- call the perturbation choosing a coupling constant [ALPHA], a
  reaction coordinate (see summary below), and a selection of atoms 
  which define the reaction coordinate. Several biases may be 
  in operation in any time: each must be set up by a separate
  HQBM command. The general form of the setup command is:

  ::
  
     HQBM [RC1 | RC2 | RC3...] ALPHA real [IUNJ integer] [XIMAX real] -
           [ANAL FIRSTU integer NUNIT integer] -
           coord-specific-options

       coord-specific-options are listed below for each coordinate

- energy NO LONGER NEEDS TO BE CALLED after HQBM !!
  this won't affect anything, just increase the step number by 1 each time.
  necessary in order to have multiple reaction coordinates & keep the
  output synchronous.

- reset all HQBM biases, i.e. EHQBM = 0.0 always
  
  ::
  
   HQBM RESET

- only change the coupling constant (ALPHA); useful for equilibration
  
  ::
  
   HQBM RCX UPALPHA real
  
  RCX is RC1, RC2, etc. - which ever reaction coordinate needs
  alpha updated

  'real' is a new value for the coupling constant.

- coord-specific-options:
  
  A description of each coordinate, and the options is given in the
  section Function. Also, this will surely be out of date rapidly,
  so the source is the best recourse.
  
  ::

    RC1: [AWAY] [SMD GAMMA real] [FIX] [NOEN] [READLIST integer] -
         [READREF integer] [IUNK integer] atom-selection

    RC2: [AWAY] [SMD GAMMA real] [FIX] [NOEN] [READLIST integer] -
         [IUNK integer] atom-selection
    
    RC3/PHI: [AWAY] [SMD GAMMA real] [FIX] [NOEN] [COMB] [AVEP AALPHA real] -
         [READLIST integer] [IUNR real] [BETA real] [EXCL real] -
         [RCUT real] [TOL real] [ZERO] [IUND integer] -
         IUNP real atom-selection
    
    RC4/HX: [AWAY] [SMD GAMMA real] [FIX] [NOEN] [IUNK integer] [IUND integer] -
         [EEF1] [NHCON] [SPLIT] [NONN [CUTON real] [CUTOF real]] [BETA real] -
         [BETC real] [BETH real] [EXCL real] [RCUT real] [HCUT real] [ZERO] -
         IUNP integer -
         atom-selection1 atom-selection2 atom-selection3 atom-selection4

    RC5: [NOEN] [TARGET real] [READLIST integer] atom-selection

    RC6/NOE: [AWAY] [SMD GAMMA real] [FIX] [ZERO] [SIXT | LINE] [NOEN] -
         [IUND integer] IUNN integer

    RC7/RDC: ... not done yet ...

    RC8/S2 ***: [IUND integer] [FIX] IUNS integer 

    RC9/J3: [IUND integer] [IUNK integer] [ZERO] [NOEN] J3UNIT integer

    RC10/PSI ***: [IUND integer] [IUNK integer] [FIX] [ZERO] [BETA real] 
          [RCUT real] [TOL real] IUNP integer

    *** These coordinates can ONLY be used in the replica/ensemble version.


.. _hqbm_function:

The following section describes the keywords of the HQBM command.

HQBM introduces a half quadratic perturbation on a given reaction
coordinate (see below)

Meaning of the HQBM parameters
------------------------------

General Parameters & Parameters common to many coordinates 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

(check syntax to see whether a given option is supported with the reaction
coordinate of interest)

* AWAY drive the system away from the reference coordinate.
  As an example, if the reaction coordinate measures the deviation from
  a reference conformation, the perturbation will increase it.

* ALPHA is the force constant of the half harmonic potential.

* RC1, RC2, RC3/PHI, RC4/HX, RC5, RC6/NOE, RC7/RDC, RC8/S2, RC9/J3, RC10/PSI
  will select other reaction coordinates (descriptions below)

* atom-selection

  some coordinates require an atom selection -
  only the selected atoms will be used to define the coordinate.
  See below for more specific definitions. 

* IUNJ

  write the output (istep rc(t) max(rc)) on unit IUNJ

* FIX
  
  make the target value of the reaction coordinate the initial value.

* ZERO

  make the target value of the reaction coordinate ZERO (same as FXRG).

* IUND integer

  a unit to dump calculated phi-values, protection factors to
  at regular intervals during the trajectory

* IUNK integer

  a unit to dump initial contact lists to.

* SMD

  use schulten style "steered molecular dynamics". This requires
  a speed to move the target reaction coordinate, given by the
  GAMMA option.

* NOEN

  when using the ensemble version of the code (see: ensemble.doc)
  this will force a particular reaction coordinate NOT to use
  the ensemble averaged form.

* BETA real

  the value of beta in the smooth function for counting native
  contacts 1.0/(1+exp(beta(r-rcut))).

* RCUT real

  see entry for BETA above.

* TOL real

  When counting native contacts in non-native structures, allow
  an extra TOL angstroms (i.e. rcut is increased by TOL).

* EXCL integer

  Do not count contacts between residues separated by fewer
  than EXCL.

Description of each coordinate and its specific parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*  RC1

   A reaction coordinate based on the mean square difference from the
   target coordinates. If the target coords are all set to zero
   (e.g. with SCALAR), the reaction coordinate is like a radius
   of gyration (it is in fact the square of the radius of gyration
   over the selected atoms assuming equal masses). If only two atoms 
   are selected, the reaction coordinate is the distance between them.
   
   ::
        
        [READLIST integer] read a list of atom index pairs specifying native
                                contacts, i.e. in the format:
                                i1 j1
                                i2 j2
                                ...
         [READREF integer] read a list of atom index pairs specifying native
                                contacts, AND distances between them, i.e.:
                                i1 j1 r1
                                i2 j2 r2
                                ...

*  RC2

   Works exactly like RC1, except that instead of ``rho = \sum_ij (r_{ij}-r_{ij}^{ref})^2``,
   ``rho = sum_{ij} exp(((r_{ij}-r_{ij}^{ref})/r_{ij}^{ref})^2)``.

*  RC3/PHI

   Drive system to satisfy experimental phi-values, defined as a residue-based
   fraction of native contacts.
   
   ::
   
      [COMB] : if specified, the native contact list will be constructed
             by making all possible combinations of the atom selection.
             Used for hydrophobic clustering in unfolded state (Julia Wirmer).
      [AVEP AALPHA real] : ONLY works with ensemble code. As an ensemble,
             the replicas are driven to satisfy the expt phi-values; the
             AVEP bias ensures that each replica will also satisfy the 
             average phi value, AALPHA being a separate coupling constant
             for this. Only one HQBM invokation is needed for both the
             standard phi and the average phi (by default average phi is off).
      [READLIST integer] : read native contacts from a file:
                             i1 j1
                             i2 j2 
                             ...
      [IUNR real] ????
      [IUNP real]: unit with phi-values:
                     res1 phi1
                     res2 phi2
                     ...
      atom-selection: the atoms to use for counting native contacts if
             not reading native contact list from a file.

*  RC4/HX

   Hydrogen exchange bias. System driven to satisfy experimental 
   protection factors. Protection factors defined as logP = Bc*Nc+Bh*Nh

   ::
   
         atom-selection1: defines heavy atom contacts
         atom-selection2: oxygen selection (for hbonds) 
         atom-selection3: nitrogen selection (only for EEF1 - otherwise ignored)
         atom-selection4: hydrogen selection (for hbonds)

         [EEF1] - this ONLY works in analysis mode. The EEF1 energy of nitrogen
                atom is used for the burial term (Nc). Uses third atom selection.
         [NHCON] - used HN_i --- heavy atom contacts for burial 
                default is heavy_atoms_i --- other heavy atoms
         [SPLIT] - when writing to IUND file, separate hydrogen bonding and 
                burial contributions to the protection factor.
         [NONN [CUTON real] [CUTOF real]] - Use all contacts, not just native
                ones, for burial. Requires a cutoff function for efficiency.
                cutof must be larger than cuton.
         [BETC real] = bc above
         [BETH real] = bh above
         [HCUT real] - cutoff for counting hydrogen bonds (default = 2.4 Angstrom
                        O-H distance)
         [IUNP integer] - unit with protection factors:
                        res1 logP1 type1
                        res2 logP2 type2
                        ....
                    The protection factor "type" is one of 0, -1, or 1:
                    0: protection factor must be satisfied exactly
                    -1: protection factor must be smaller than value given
                        (for residues exchanged in dead time)
                    1: protection factor must be larger than value given
                        (for global exchange data)

*  RC5

   Works like RC1, except drives system towards target value specified by TARGET
   and holds it there.

*  RC6/NOE

   Drives system towards experimental NOE values.
   
   ::
   
         [SIXT | LINE] - type of averaging. Default is <r^{-3}>^{-1/3}
                        SIXTh specifies <r^{-6}>^{-1/6}
                        LINEar is normal (linear) averaging

         IUNN integer - unit with noe's, format:
                N
                i1 j1 lbound1 ubound1
                i2 j2 lbound2 ubound2
                ...
                iN jN lboundN uboundN

*  RC7/RDC: not implemented

*  RC8/S2

   Order parameter bias. Drives an ensemble of configurations
   to satisfy experimental order parameters. Obviously, this
   ONLY works for the ENSEMBLE code (see ensemble.doc).

   ::
   
        IUNS integer - unit with order parameters, format:
                N
                i1 j1 S2_1
                i2 j2 S2_2
                ...
                iN jN S2_N

*  RC9/J3: Drive system to satisfy scalar coupling restraints

   ::
   
        J3UNIT integer  - unit with couplings, format:
                i1 j1 k1 l1 A1 B1 C1 D1 J1
                i2 j2 k2 l2 A2 B2 C2 D2 J2
                ...
                where i,j,k,l are the atom indices defining
                the dihedral, and A, B, C and D are the
                karplus parameters using the form of the equation:
                J(phi) = A*cos^2(phi+D) + B*cos(phi+D) + C
                Ref: Chou et al. JACS, 125, 8959-8966 (2003)

*  RC10/PSI

   Drive system to satisfy psi-values (sosnick papers)
   (not finished...)

The method is  described in
E. Paci and M. Karplus.  Forced unfolding of fibronectin type 3
modules: An analysis by biased molecular dynamics simulations.
J. Mol. Biol., 288: 441-459, 1999.

TESTCASES (in test/c32test)
---------------------------

*  hqbm_single_test.inp

   This is a test of the single copy versions of 
   RC1, RC2, RC3, RC4, RC6 & RC9
   It may be run in the test directory by invoking:
   ./test.com arch output bench 32
   which will run this + all the other c32 testcases

* hqbm_rc3_ens_test.inp: Ensemble test of RC3/PHI -- see below for how to run
* hqbm_rc4_ens_test.inp: Ensemble test of RC4/HX -- see below for how to run
* hqbm_rc8_ens.inp: Test of RC8 (only ensemble)  -- see below for how to run

To run ensemble tests, use the following command in the test directory:

::

   ./test.com E arch 
   
in this case the optional fourth command specifying target will be ignored.

