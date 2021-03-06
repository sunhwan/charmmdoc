.. py:module:: mcma

=========================================
Monte Carlo Minimization/Annealing (MCMA)
=========================================

The MCMA commands modify the system coordinates, facilitiating 
conformational searches for the global energy minimum of a macromolecular
system.  Rigid-body translation and/or rotation of a subset of atoms is
supported for use in docking applications.  Also supported are rotations
about single bonds, biased to favor conformations observed in 
high-resolution crystal structures of proteins.  As of August 2004, this
biased rotation about single bonds has been tested more extensively than
have the docking moves, and this biasing is implemented only for phi, psi,
chi1, chi2, and pre-proline omega angles of amino acids.  Other dihedral
angles are changed without bias.

The MCMA commands only change atomic coordinates (internal and/or
Cartesian).  They were designed for use with a CHARMM script 
(e.g., mcma.inp included in the test cases) that evaluates the energy 
and accepts/rejects trial conformations.  Consequently, the MCMA commands
can be used with any energy function implemented in CHARMM.  These commands
assume that the "main" IC table contains only those dihedral angles to be
modified by MCMA moves (e.g., single bonds).  As in mcma.inp, the "saved"
IC table can be used to store all ICs, from which the entire structure 
can be built.

.. _mcma_syntax:
    
MCMA commands
-------------

The initialization of the MCMA routines must be performed once
before any MCMA moves can be performed.  The form of the command is:

::

   MCMA INIT  char1  char2  real1  real2  real3  real4  real5  int     (1)

The two character and first four real values are only relevant in
ligand-docking applications: char1 and char2 specify the CHARACTER*4
segid (SEGMove) and resid (RESMove) of the ligand, respectively.  The
first three real values contol changes made to the Euler angles (phi,
psi,theta) describing the ligand's orientation.  Changes to cos(theta)
are uniformly distributed between +- real1.  Changes to the ligand phi
angle (radians) are uniformly distributed between +- real2.  Changes to
the ligand psi angle (radians) are uniformly distributed between +- real3.
Changes to the x, y, and z coordinates of the ligand are each uniformly
distributed between +- real4.  In biased searches of peptide conformations,
real5 specifies a scale factor applied to the standard deviations of the 
Gaussians used to bias dihedral-angle changes.  setting real5 = 1.0 
(recommended) invokes the biasing of Abagyan and Totrov, JMB, 1994.  
Finally, the integer value specifies the seed used in random-number
generation.   

Upon invoking the MCMA INIT command, all entries in the "main" IC
table are classified as phi, psi, omega, or chi_n angles (n=1,2,3,4).  
NOTE: Identification is done based on atom type.  That is, it is assumed
that CA = alpha carbon, CB = beta carbon, etc.  The results of this
classification are output.

Conformational changes are made by the following MCMA commands:

::

   MCMA ROTA
   MCMA TRAN
   MCMA ROTR
   MCMA BIAS int1
   MCMA SIDE int1
   MCMA MAIN int1  int2  real
   MCMA BACK int1  int2  real
   MCMA ALL

The first three moves change the coordinates of a ligand (segid = SEGMove,
resid = RESMove) by a rigid-body ROTAtion, TRANslation, or both, respectively.
The magnitude of the coordinate changes is governed by the maximum values 
read in the MCMA INIT command (above).  

The remaining moves modify dihedral angles of one or more entries in the IC
table.  The BIAS move randomly selects int1 entries in the IC table.
For each entry, biased random changes are made to it and to its partner angle,
if any.  That is, phi-psi and chi1-chi2 angles are changed simultaneously
as in Abagyan and Totrov, 1994.  A value of 1 is recommended for int1,
as acceptance decreases when more than one angle or pair of angles is changed
simultaneously.

The SIDE move is a special-case BIAS move, restricted to select int1
side-chain torsions.

The MAIN move randomly selects a residue as the first of int1 consecutive
residues for which the (phi,psi) angles are changed to the same
values, biased according to the identity of the first residue.  By moving
int1 consecutive residues to the same point in the Ramachandran plane, MAIN
moves can facilitate the creation of alpha helices.

The BACK move randomly selects a residue as the first of int1 consecutive
residues for which the (phi,psi) angles are changed, each biased
independently according to the identity of the respective residue.

For BACK and MAIN moves, int2 (int2 .ge. 0) gives the total number of steps of
molecular dynamics used to heat and cool the system following the
dihedral-angle changes and prior to accepting or rejecting the move.
The maximum temperature to be reached (TMAX) is given by the real value
(e.g., 450.0).  The system is heated to about TMAX in 0.1*int2 steps and
cooled to about 200 K in 0.9*int2 steps.

NOTE:  The use of SHAKE is assumed in BACK and MAIN moves for which int2 > 0;
a 2-fs step is used during annealing.

NOTE: the BIAS, SIDE, MAIN, and BACK moves change the internal and Cartesian
coordinates. 

For side-chain moves, the Cartesian coordinates of the entire side chain are
rebuilt.  For main-chain moves, the given residue and all preceding or
following it are rebuilt, depending on which choice involves the fewest
residues.  In future implementations, it may be desirable (e.g., for loop
modeling) to allow more control over the atoms to be rebuilt following any
of the MCMA moves.

The ALL move performs a BIAS move for each IC entry (or pair of entries).
Unlike the BIAS, MAIN, BACK, and SIDE moves, the ALL move only changes the
IC values, allowing the user to specify the coordinates to be IC-BUILt after
the move.  This could be useful when modeling loops: the IC table would
involve only the loop residues and the coordinates to be "IC BUILt" would
be those of the loop atoms.


.. _mcma_assumptions:

Assumptions made in current implementation
------------------------------------------

(as of August, 2004)
 
The MCMA module was written with two primary applications in mind,
ligand docking and peptide-structure prediction.  The ligand-docking moves
have not been tested as much as the other moves.
The following assumptions/restrictions are imposed.

1. Identification of dihedral angles (in the IC table) is done based on atom
   type.  Atom types used in CHARMM19 and CHARMM22 are assumed: CA = alpha carbon,
   CB = beta carbon, etc.  Any IC entries not recognized as a protein phi, psi,
   omega, or chi angle are changed without bias.  The results of this
   identification are output.

2. When analyzing the IC table, it is assumed that psi and chi2 ICs
   immediately follow the corresponding phi and chi1 ICs, respectively.

3. The use of SHAKE to constrain the lengths of bonds to hydrogen atoms is
   assumed in BACK and MAIN moves for which the number of annealing steps (int2)
   is greater than zero.  A 2-fs integration time step is used during simulated
   annealing by molecular dynamics.

4. The BIAS, MAIN, BACK, and SIDE moves modify the internal and Cartesian
   coordinates.  For side-chain moves, the Cartesian coordinates of the entire
   side chain are rebuilt.  For main-chain moves, the given residue and all
   preceding or following it are rebuilt, depending on which involves the fewest
   residues.  The ALL move, which modifies all dihedrals in the IC table, was
   included so as to allow the user to specify the atoms to be rebuilt
   (see mcma.inp).

5. Arrays X0, Y0, and Z0 are used for temporary storage of coordinates.
   As of c32a1 (August 2004), these arrays have two conflicting uses:
   
   i) For docking: ligand coordinates relative to ligand center of mass, and
   ii) For peptide structure prediction: storage in MAIN and BACK moves.
   
   Therefore, the code must be modified to do docking and main-chain moves
   simultaneously.


.. _mcma_example:

Example
-------

See mcma.inp in the test cases for a peptide-folding run.


.. mcma_references:

References
----------

(1) P J Steinbach, Exploring Peptide Energy Landscapes: A Test of Force Fields
    and Implicit Solvent Models; Proteins, in press (2004)

(2) R Abagyan and M Totrov, J Mol Biol (1994)

(3) Li and H A Scheraga, ; Proc Natl Acad Sci USA  (1987)
