.. py:module:: ace

============================================
Analytical Continuum Solvent (ACS) Potential
============================================

Purpose: calculate solvation free energy and forces based on
a continuum description of the solvent, in particular the analytical
continuum electrostatics (ACE) potential.

Please report problems to Michael Schaefer at schaefer@piaf.u-strasbg.fr

.. warning::
   The module is still being developed and may change in the future.

.. warning::
   Note on ACE2: the version 2 of ACE as of Jan 2002 is not yet fully
   parameterized; it yields reasonably stably MD trajectories of native
   proteins when using param19 (united atom parameters), but is         
   unreliable with all-hydrogen parameters.                             


REFERENCES:
*  M. Schaefer & M. Karplus (1996) J. Phys. Chem. 100, 1578-1599.
*  M. Schaefer, C. Bartels & M. Karplus (1998) J. Mol. Biol. 284, 835-847.
*  N. Calimet, M. Schaefer & T. Simonson, (2001) Proteins 45, 144-158
*  M. Schaefer, C. Bartels, F. Leclerc& M. Karplus (2001),
   J. Comp. Chem. 22, 1857-1879.


.. _ace_syntax:

Syntax
------

::

   Syntax: The ACE specifications can be specified any time the nbond 
           specification parser is invoked, e.g., 
   	ENERgy [other-spec] [ace-spec]

   ace-spec::=
           [ ACE ] [ IEPS real ] [ SEPS real ] [ ALPHa real ]
           [ SIGMa real ] [ IDEAl | CURRent ] [FVSCaling real]
           [ ACE2 [ MXBSolv real ] [TBSOlv real ] [ TBSHydrogens real ]]
 

.. _ace_defaults:

Defaults
--------

The defaults for the ACE potential are:

===== ====
IEPS	1.0
SEPS	80.0
ALPHa	1.3
SIGMa	0.0
IDEAl true
FVSCa 1.0
===== ====

The additional defaults for the ACE2 potential are:

=====   ====
MXBSo   14.0
TBSOl   8.4
TBSHy   3.85
=====   ====

In the current implementation, ACE should be used with united atom parameters,
ALPHa set equal to 1.3, the standard PARAM19 parameter file param19.inp and
Voronoi volumes as given in acepar19.inp (toppar and test/data).


.. _ace_function:

Introduction
------------

The analytical continuum solvent (ACS) potential is introduced to
perform molecular dynamics/minimization calculations with a continuum
approximation of the solvent.

Two solvent contributions to the effective (free) energy of a solute
are included: the electrostatic solvation free energy, and the
non-polar (i.e., non-electrostatic) solvation free energy.
The first (electrostatic) contribution (G_el) is calculated using an
analytical approximation to the solution of the Poisson-equation
called ACE (from: analytical continuum electrostatics).
The non-polar solvation free energy (G_np) is approximated by a pairwise
potential which yields results that are very similar to the well-known
surface area approximations of the hydrophobic (solvation) energy
(e.g., Wesson and Eisenberg, Prot. Sci. 1 (1992), 227--235; see
the ASP potential in CHARMM).

Restriction
-----------

The ACE solvation potential has to be used together with no cutoff or with
atom based switching.

Compatibility
-------------
1. ACE can be used with BLOCK (but: the diagonal elements of the BLOCK
   matrix MUST NOT be zero).

2. ACE can be used with fixing atoms (CONS FIX); the resulting energy and
   forces are an approximation, because all the interaction-dielectric terms
   of the potential (eq (47) in Schaefer & Karplus, JPC 100 (1996), 1578)
   which involve two fixed atoms are neglected, despite the fact that they
   exist and that they are not invariant!

Meaning of the ACE parameters
-----------------------------

1.  IEPS 

    Dielectric constant for the space occupied by the atoms that are treated
    explicitly, e.g., the space occupied by the protein.

2.  SEPS

    Dielectric constant for the space occupied by the solvent that is treated
    as a continuum (i.e., the complement of the space occupied by the protein).

3.  ALPHa

    The volumes occupied by individual (protein) atoms are described by
    Gaussian density distributions. The factor ALPHa controls the width of these 
    Gaussians. The net volume of the individual atom Gaussian distributions is
    defined in the volume table in the parameter file acepar19.inp.
    The volumes in the acepar19.inp file are expected to work best
    for an ALPHa of 1.3.

4.  SIGMa

    The ACE solvation potential includes a hydrophobic contribution
    which is roughly proportional to the solvent accessible surface area.
    The factor SIGMa scales the hydrophobic contribution. For peptides
    with about 10-15 residues, a SIGMa factor of 3.0 results in hydrophobic
    contributions that are approximately equal to the solvent accessible 
    surface area multiplied by 8 cal/(mol*A*A).

4.  IDEAl | CURRent
    
    As of c29a2, the ACE potential considers the distances between atoms
    in the nonbonded exclusion list as invariant. This is consistent with
    the assumption that the forces involving these atoms are governed by
    the internal energy terms (bond, angle, and some 1-4 atom pairs in
    aromatic ring systems). Note that solvation forces still apply to
    pairs of these atoms, considered as a polar group.
    
    With the IDEAl option (default), ACE calculates the nonbonded exclusion
    list distances from ideal bond length and angles where possible; the
    distances for 1-4 atom pairs in the exclusion list are calculated
    from the current atom positions at the first ACE energy call.
    With the CURRent option, all the distances between atoms in
    the nonbonded exlusion list are calculated from the current
    coordinates of the atoms. These distances are considered invariant
    for all subsequent energy calls, during minimization and dynamics.
    Recalculation of the nb-exclusion list atom pair distance is
    enforced only when toggling IDEAl on/off, fixing/unfixing atoms,
    or a change of the psf (e.g., REPLica).
    
4.  FVSCal

    One major problem with ACE1 (and gneralized Born methods in general)
    is the overestimation of the desolvation by the pairwise de-screening
    function ESELFIK (see ace.src). One way to reduce the impact of this
    systematic error is to reduce the volume that is assigned to the atoms
    by a constant factor FVSCal < 1 as proposed in Calimet et al., Proteins
    45 (2001), 144-158. The default value for FVSCal is 1.0, though a value
    of 0.9 appears reasonable in conjunction with param19 and volumes
    in acepar19, using the ACE1 potential (work in progress). Note that
    the modified treatment of the self energy (de-screening) potential
    in ACE2 is aimed at fixing the overestimation problem of ESELFIK
    such that the re-scaling of volumes becomes obsolete (work in progress).
    
4.  ACE2

    The ACE2 keyword implies ACE (no need to specify both). It invokes
    a modified treatment of the Born solvation radii which are limited
    by un upper bound --- MXBSolv (see below). This takes account of the
    overestimation of the desolvation of charges by the pairwise de-screening
    potential in ACE1.
    
4.  MXBSolv

    The Born solvation radii of all atoms (charges) are limited
    by the upper bound parameter MXBSolv (default 14.0 Angstrom).
    
4.  TBSOl

    In the ACE2 potential, the conventional conversion of the atomic
    solvation to the Born solvation radii is applied until a Born radius
    of TBSOlv is obtained ("turning point"). After that, atomic solvation
    energies (i.e., the de-solvation) is converted in a way that prevents
    the Born solvation radii from exceeding the imposed maximum.
    Details will be given in an upcoming publication. 
    
4.  TBSHyd

    This parameter has the same meaning as TBSOl, but applies
    to hydrogens, which are most susceptible to an overestimation
    of the desolvation by neighboring atoms (volumes). The smaller
    the TBSOl and TBSHyd, the more the over-desceening is counter-
    acted (parametrization in progress).


.. _ace_examples:

Examples
--------

To set up simulations/minimizations with the ACE solvation potential,
read the standard CHARMM topology and parameter files and the corresponding
ACE parameter file using

::

   read ACEParameters card unit IUN

e.g., the file acepar19.inp with param19 parameters.
The following energy call is expected to be adequate for most cases,
including proteins:

::

   ENERgy ATOM ACE2 IEPS 1.0 SEPS 80.0 ALPHa 1.3 SIGMa 2.5 SWITch -
          VDIS VSWI CUTNB 13.0 CTONNB 8.0 CTOFNB 12.0

When you run molecular dynamics or minimization with ACE, you get
two more lines in the log file printout with energy terms, e.g.,

::

   DYNA DYN: Step         Time      TOTEner        TOTKe       ENERgy  TEMPerature
   DYNA PROP:             GRMS      HFCTote        HFCKe       EHFCor        VIRKe
   DYNA INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers
   DYNA EXTERN:        VDWaals         ELEC       HBONds          ASP         USER
   DYNA PRESS:            VIRE         VIRI       PRESSE       PRESSI       VOLUme
   DYNA ACE1:      HYDRophobic         SELF    SCREENing      COULomb 
   DYNA ACE2:        SOLVation  INTERaction 
    ----------       ---------    ---------    ---------    ---------    ---------
   DYNA>        0      0.00000  -3423.29671      0.00000  -3423.29671      0.00000
   DYNA PROP>          4.45310  -3423.12228      0.52327      0.17442   -532.70519
   DYNA INTERN>        6.58717     60.43092      0.00000     56.00750      7.32144
   DYNA EXTERN>     -380.26218  -3173.38156      0.00000      0.00000      0.00000
   DYNA PRESS>         0.00000    355.13679      0.00000      0.00000      0.00000
   DYNA   ACE1>      109.04469  -3829.20991   2750.59427  -2203.81062
   DYNA   ACE2>    -1078.61564    546.78365
    ----------       ---------    ---------    ---------    ---------    ---------

and the same during minimization (MINI...) or after
an energy calculation (ENER...).

The terms in lines with ACE1 and ACE2 are:

===========  ==================================================================
HYDRophobic  Hydrophobic potential, equivalent to a surface based
             solvation term proportional to the sigma input parameter;

SELF         Self contribution to electrostatic solvation free energy,
             Delta-E_self, first term of eq(8) (i.e., sum over all atomic
             solvation energies, Delta-E_self_i, eq(28));

SCREENing    Interaction contribution to electrostatic solvation free energy,
             i.e., screening of Coulomb interactions, eq(38) (sum over all
             atom pairs, including bonded and 1-3 atom pairs!);

COULomb      Coulomb energy with constant dielectric of EPSI (sum over
             all atom pairs for the first term in eq(36) -- excluding
             bonded and 1-3 atom pairs, and 1-4 atom pair contributions
             scaled with E14FAC);

SOLVation    Electrostatic (!) solvation free energy, sum of SELF and
             SCREENing;

INTERaction  Electrostatic interaction, sum of SCREENing and COULomb
             (eq(36), but taking account of the bonded, 1-3, and 1-4
             exclusion in the Coulomb term, see above).
===========  ==================================================================

The term "ELEC" in line "DYNA EXTERN>..." is the total electrostatic energy: 

===========  ==================================================================
ELEC:        Sum of SELF, SCREENing, COULomb.
===========  ==================================================================

Equation numbers refer to Schaefer & Karplus, J. Phys. Chem. 100 (1996), 1578.

See also: test cases c27test/ace1.inp and c29test/ace_v2.inp.

