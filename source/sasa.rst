.. py:module::sasa

=================================
The SASA implicit solvation model
=================================

Characteristics of the SASA model
---------------------------------

The SASA model is a fast implicit solvation model that is useful to
simulate structured peptides and miniprotein motifs [1]. The polar and
non-polar contributions of each atom to the free energy of solvation are
assumed to be proportional to their solvent accessible surface areas.
The SASA model uses only two surface-tension like solvation parameters
(constants of proportionality) and approximates the solvent accessible
surface area of each solute atom with a simple analytical function that
is easily derivable. The electrostatic screening between solute charges
is accounted for by using a distance dependent dielectric function and
by neutralizing the formal charges (Asp, Glu, Arg, Lys, and the termini)
as in the EEF1 model [2].

The SASA model has been successfully applied to peptides, removing the
major artifacts of in vacuo simulations and reproducing reversible
folding [3]. Benchmarks indicate that a simulation with SASA is only
about 50% slower than an in vacuo simulation.

Range and limitations
---------------------

The SASA model has been applied to structured peptides, see for instance
[3]. However, it should not be used for large proteins mainly for two
reasons. Firstly, it has been parameterized for small proteins [1] and
secondly, the dielectric function does not take different environments
into account, i.e., it does not distinguish whether or not the
interacting partial charges are buried or on the protein surface.

Theoretical aspects
-------------------

The potential energy of the system consisting of the solute and the
solvent can be decomposed in three parts: the intra-solute potential
energy U(X), the intra-solvent potential energy V(Y), and the
interaction potential energy of the solute and the solvent W(X,Y), where
X denotes the degrees of freedom of the solute and Y the degrees of
freedom of the solvent. Integrating out all solvent degrees of freedom
one obtains the potential of mean force W(X), also called the effective
energy. It can be written as the sum of the intra-solute potential
energy U(X) and the free energy of solvation, or mean solvation term,
DW(X) that describes all solvent induced effects. This is a rigorous
result from statistical mechanics. For more details consult, for
instance, the review [4].

Any implicit solvation model based on the solvent accessible surface
area approximates the major contribution to the free energy of solvation
using the first solvation layer, and the screening effect of the
solvent. Following this idea it is assumed that the mean solvation term
is a sum of atomic contributions proportional to the solvent accessible
surface area of each atom, plus an energy term that accounts for the
screening. The constants of proportionality are the surface-tension like
solvation parameters.

Technical aspects of the SASA model
-----------------------------------

There are three important aspects for any implicit solvation model that
is based on the solvent accessible surface area: (A) how to calculate
the atomic solvent accessible surface areas, (B) how many
surface-tension like solvation parameters to use, and (C) how to account
for the screening of solute charges. 

(A) Calculation of the solvent accessible surface area

    In the SASA model the solvent accessible surface area is approximated by
    a probabilistic approach. For details please refer to [5]. The base is a
    simple formula for the probability to hit the accessible surface of atom
    i if N atoms are present by choosing randomly a spot on the solvation
    shell of atom i. However, this formula is only true under the assumption
    that all the atoms are distributed randomly. This is of course not the
    case because of covalent geometry and the Pauli exclusion principle. Now
    instead of elaborating in sophisticated probability calculations the
    pragmatic approach of [6] is adapted where the original formula is
    parameterized by including two sets of parameters: a set of
    probabilistic parameters (called atom type parameters in [1,6]) and a
    set of connectivity parameters. The probabilistic parameters depend on
    the atom type and correct for systematic errors primarily due to
    hybridization. The connectivity parameters distinguish bound atoms from
    more distant ones. Also, the radii used to calculate the approximated
    solvent accessible surface were optimized for this purpose. (These radii
    are used only for the calculation of the approximated area and they are
    different from the radii used for the CHARMM van der Waals energy term.)
    The optimal values for the radii, probabilistic parameters, and connectivity
    parameters were taken from [6]. The atom types in [6] and the CHARMM
    atom types do not match exactly, so the most reasonable assignments were
    chosen, analogous to the choices in [7].
  
    The only internal degree of freedom considered by this approximation is
    the interatomic distance. Therefore the major defect is the absence of a
    better correlation between exact and approximated areas upon changes in
    internal coordinates like dihedral angles.

(B) Solvation parameters

    By default, the SASA model uses only two non-vanishing surface-tension
    like solvation parameters: one for hydrophobic and one for hydrophilic
    groups. The solvation of explicit hydrogen atoms is neglected. However,
    this can be changed by the user (see below). The solvation parameters
    were adjusted for the EEF1-modified CHARMM 19 polar hydrogen parameter
    set by Philippe Ferrara in a trial and error approach in 1999. The
    criterion was to minimize the root mean square deviation from the native
    state for six small proteins by performing molecular dynamics
    simulations of 1ns at 300K [1].

(C) Screening of solute charges

    For the screening of solute charges the SASA model uses a distant
    dependent screening function, eps(r)=2r, and neutralizes the charged
    groups of polar amino acids in exactly the same way as it is done in the
    EEF1 model of Lazaridis and Karplus [2].

Implementation in CHARMM
------------------------

To every CHARMM atom type the SASA model assigns a surface-tension like
solvation parameter (zero by default for explicit hydrogen atoms and
non-zero for hydrophobic and hydrophilic groups), a radius optimized for
the approximation of the solvent accessible surface area, and a
probabilistic parameter. Additionally, a connectivity parameter is
assigned to every pair depending on whether the two atoms of the pair
form a 1-2 or a more distant pair. These parameters are called the SASA
parameters [1].

When initializing SASA, i.e., by invoking the SASA command, all the SASA
parameters are printed to the CHARMM output file. A value of -999.000
means that this specific parameter hasn't yet been determined for SASA
and therefore the corresponding CHARMM atom type can't be used by
default for the SASA calculations - the user would have to assign a
meaningful value. (Currently, SASA does not support the following CHARMM
atom types by default: HA, HT, LP, CT, CM, NP, OH2, OM, OT, OS, and FE.)
All SASA parameters can be changed by the user. For the corresponding
syntax, see below.

To evaluate the approximative formula for the solvent accessible surface
area of an atom i, all neighbors of atom i within a certain cutoff are
required. This cutoff is calculated by 2*(2.365A+1.4A)=7.53A, where
2.365A is the largest van der Waals radius in the CHARMM parameter set
19, and 1.4A is the radius of the solvent probe sphere. Most, but not
all of these neighbors, are included in the nonbond pair list in CHARMM.
The missing pairs are stored in a new pair list, the SASA pair list.
More precisely, the SASA pair list contains all pairs that (1) are not
in the nonbond pair list, that (2) do not belong to the fixed exclusions
(as given in the topology file), and that (3) are within the above
mentioned cutoff. The SASA pair list assures correct operation of the
SASA model for any nonbond exclusion mode (-5 to 5). This means that for
any nonbond exclusion mode from 1 to 5, the SASA energy and its
derivatives are identical. The same is true for any nonbond exclusion
mode from -5 and 0. The differences in the SASA energy and its
derivatives between the nonbond exclusion modes regions of -5 to 0 and 1
to 5 stem from the fact that the fixed exclusions are treated
differently: For an exclusion mode from -5 to 0 they are included in the
nonbond pair list, opposed to an exclusion mode from 1 to 5 where they
are excluded from the nonbond pair list.

Caveat
------

Please note that the radii used for the calculation of the approximated
solvent accessible surface area are labeled 'van der Waals radii' in the
SASA output of the CHARMM output file. This is consistent with the
terminology used in [6]. However, these radii are different from and do
not replace in any way the CHARMM default van der Waals radii used to
calculate the van der Waals energy term.

Additional Input Files
----------------------

Two additional files are needed to use SASA (taken from EEF1):

(1) toph19_eef1.inp : This is a modification of toph19.inp where ionic
    side chains and termini are neutralized and contains
    an extra parameter type (CR).
(2) param19_eef1.inp: This is a modification of param19.inp which includes
    the extra parameter type (CR).

These files can be found in test/data/.

.. _sasa_syntax:

Syntax of the SASA command
--------------------------

There is only one SASA command:

::

   SASA atom-selection [S<number> <real>] [R<number> <real>] [P<number> <real>] [fcon <real>] [ncon <real>] [surf] [infx] [newp]


   atom-selection:== (see *note select:(chmdoc/select.doc).)

   <number>:== number corresponding to a CHARMM atom type from param19.inp
               or param19_eef1.inp according to the following list:

   number = 001 for H
   number = 002 for HC
   number = 003 for HA
   number = 004 for HT
   number = 005 for LP
   number = 006 for CT
   number = 007 for C
   number = 008 for CH1E
   number = 009 for CH2E
   number = 010 for CH3E
   number = 011 for CR1E
   number = 012 for CM
   number = 013 for C1ES
   number = 014 for N
   number = 015 for NR
   number = 016 for NP
   number = 017 for NH1
   number = 018 for NH2
   number = 019 for NH3
   number = 020 for NC2
   number = 021 for O
   number = 022 for OC
   number = 023 for OH1
   number = 024 for OH2
   number = 025 for OM
   number = 026 for OT
   number = 027 for OS
   number = 028 for S
   number = 029 for SH1E
   number = 030 for FE
   number = 031 for CR

The SASA command sets up the SASA model for a simulation. All the values
have to be given on one command line. Invoking the SASA command a second
time reinitializes all values either to the default or to the user
specified values.

------------------ --------------------------------------------------------
atom-selection     This determines the atoms to be used for the SASA
                   calculations. All atoms that are not included in
                   this selection are treated by SASA as if not
                   existent. Their solvation free energies are not
                   considered, i.e., are set to zero, and the decrease
                   in the solvent accessible surface areas of the
                   selected atoms due to the not selected ones is
                   neglected.

S<number> <real>   Changes the surface-tension like solvation parameter
                   of the CHARMM atom type corresponding to <number>
                   (see list above) from the default value to <real>.

R<number> <real>   Changes the radius used by SASA (for the calculation
                   of the approximated solvent accessible surface areas)
                   of the CHARMM atom type corresponding to <number>
                   (see list above) from the default value to <real>.

P<number> <real>   Changes the probabilistic parameter of the CHARMM
                   atom type corresponding to <number> (see list above)
                   from the default value to <real>.

fcon <real>        Changes the connectivity parameter for 1-2 pairs from
                   the default value to <real>.

ncon <real>        Changes the connectivity parameter for more distant
                   than 1-2 pairs from the default value to <real>.

surf               The approximated atomic solvent accessible surface
                   areas are stored in WMAIN.

infx               Includes the fixed exclusion pairs in the SASA pair
                   list. By default, the fixed exclusion atoms are not
                   considered in the SASA surface calculations for
                   historical reasons. This means that, for instance,
                   for the CG of the residue PHE (see the topology
                   file), CZ is not considered to calculate its
                   accessible surface. To include all neighbors use the
                   keyword infx, but note that SASA was not
                   parameterized with this option.

newp               Triggers the use of the new (Haberthuer 2002) radii,
                   probabilistic parameters, and connectivity
                   parameters (see below). This keyword must be used
                   with infx.
------------------ --------------------------------------------------------

The SASA standard setup [1,3] looks like this:

::

   nbond nbxmod 5 atom rdiel shift vatom vdistance vshift -
         cutnb 8.0 ctofnb 7.5 ctonnb 6.5 eps 2.0 e14fac 0.4 wmin 1.5

   sasa selection (.not. hydrogen) end


The nonbond options are the default nonbond options from the default
param19.inp file with the exception of the eps value that is set to 2
instead of 1 and rdiel is used instead of cdiel. The default solvation
parameters are -0.06 for hydrophilic groups (N, NR, NH1, NH2, NH3, NC2,
O, OC, and OH1) and 0.012 for hydrophobic groups (C, CH1E, CH2E, CH3E,
CR1E, S, SH1E, and CR) and 0.0 for explicit hydrogen atoms. Please note
that the param19_eef1.inp file has different default nonbond options
(especially the cutoffs are different) that were not used to
parameterize SASA and therefore should not be used or used only with
care in simulations with SASA since consistency is lost.

If you select all atoms for SASA instead of excluding the (explicit)
hydrogens, all atoms will be considered for calculating the solvent
accessible surface area of each atom. However, since the solvation
parameter for hydrogen atoms is zero by default, the solvation energy of
the hydrogen atoms is also zero by default. If you insist on including
the solvation energy due to the solute hydrogen atoms, you have to
assign a non-zero solvation parameter to the hydrogen atoms by yourself,
using 'S001 <real> S002 <real>' in the CHARMM input file. Be aware that
this is not the default.

Solvation Parameters for Proteins (not fully tested)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A second set of surface-tension like solvation parameters has been
derived: -0.144 for hydrophilic groups, 0.024 for hydrophobic groups and
0.0 for explicit hydrogen atoms with eps(r)=r. It is not the default. It
seems to work better for large proteins but no results have been
published up to date (June 2004).

More Than One Chain
^^^^^^^^^^^^^^^^^^^

If you run a simulation with several molecules (so that you have more
than one segment identifier), make sure that you invoke the SASA command
after the generation of the last segment since any molecule generated
after the last use of the SASA command is not included in the SASA
calculations. You don't have to invoke the SASA command after every
generation of a segment, it is sufficient to use the SASA command once
after the generation of the last segment.

Accessing the Solvation Energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The SASA solvation energy is stored in the variable 'SASL'. Use '?SASL'
in the CHARMM input file to access the value.

New Surface Parameters
^^^^^^^^^^^^^^^^^^^^^^

A new set of surface parameters (radii, probabilistic parameters, and
connectivity parameters) was derived in 2002 by Haberthuer. They give a
better correlation with exact analytical surfaces than the original
Hasel and Still surface parameters. (A set of 20 structures [10 native
and 10 unfolded conformations] was used for the calibration.) The new
surface parameters are not the default because the surface-tension like
solvation parameters were optimized by Ferrara with the original Hasel
and Still surface parameters. The following example illustrates how to
setup SASA with the new surface parameters.

::

   nbond nbxmod 5 atom rdiel shift vatom vdistance vshift -
         cutnb 8.0 ctofnb 7.5 ctonnb 6.5 eps 2.0 e14fac 0.4 wmin 1.5

   sasa infx newp

.. _sasa_references:

References
----------

[1] Ferrara, P.; Apostolakis, J.; Caflisch, A.; Evaluation of a Fast
    Implicit Solvent Model for Molecular Dynamics Simulations; Proteins
    2002; 46; 24-33.

[2] Lazaridis, T.; Karplus. M.; Effective Energy Function for Proteins
    in Solution; Proteins 1999; 288; 477-487.

[3] Ferrara, P.; Caflisch, A.; Folding Simulations of a Three-Stranded
    Antiparallel Beta-Sheet Peptide; Proc. Natl. Acad. Sci. USA; 2000;
    97; 10780.

[4] Roux, B.; Simonson, T.; Implicit Solvent Models; Biophysical
    Chemistry; 78; 1999; 1-20.

[5] Wodak, S. J.; Janin, J.; Analytical Approximation to the Solvent
    Accessible Surface Area of Proteins; Proc. Natl. Acad. Sci. USA;
    1980; 77; 1736.

[6] Hasel, W.; Hendrickson, T. F.; Clark, S. W.; A Rapid Approximation
    to the Solvent Accessible Surface Area of Atoms; Tetrahedron
    Computer Methodology 1988; Vol. 1; No. 2; 103-116.

[7] Fraternali, F.; van Gunsteren, W. F.; An Efficient Mean Solvation
    Force Model for Use in Molecular Dynamics Simulations of Proteins
    in Aqueous Solution; J. Mol. Biol. 1996; 256; 939.

.. _sasa_example:

Examples
--------

Check test/c29test/sasa.inp for a more complex example with
minimization, equilibration and dynamics. Here is a short version:

::

   * Example input file for the SASA implicit solvation model.
   *



   ! --- Begin generation procedure ---

   open read card name toph19_eef1.inp unit 30
   read rtf card unit 30
   close unit 30

   open read card name param19_eef1.inp unit 30
   read parameter card unit 30
   close unit 30

   open read card name filename.crd unit 30
   read sequence coor unit 30
   close unit 30

   generate main warn setup

   ! --- End generation procedure ---



   ! --- Begin reading coordinates ---

   open read card name filename.crd unit 30
   read coordinate card unit 30 
   close unit 30

   ! --- End reading coordinates ---



   ! --- Begin setting up SASA ---

   ! Use the SASA standard setup.

   nbond nbxmod 5 atom rdiel shift vatom vdistance vshift -
         cutnb 8.0 ctofnb 7.5 ctonnb 6.5 eps 2.0 e14fac 0.4 wmin 1.5

   sasa selection (.not. hydrogen) end

   ! --- End setting up SASA ---



   ! --- Begin minimization ---

   minimize sd   nstep 300 nprint 20 tolgrad 0.1
   minimize conj nstep 200 nprint 20 tolgrad 0.1

   ! --- End minimization ---



   stop
