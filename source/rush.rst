.. py:module::rush

==================================================================
RUSH: A simple implicit-solvent force-field for protein simulation
==================================================================

- Olgun Guvench             (oguvench@post.harvard.edu)
- Charles L. Brooks III     (brooks@scripps.edu)

RUSH is a simple implicit-solvent force-field that adds terms
to the bonded portion (bond + angle + dihe + impr + urey) of the
all-atom CHARMM22 force field to account for volume-exclusion (_R_epulsion),
the hydrophobic effect (_U_nburied _S_urface), and intra-molecular
and protein-solvent hydrogen-bonding (_H_ydrogen-bonding) (hence
_R_ _U_ _S_ _H_).

See:

::

  Guvench and Brooks. "Folding the Trp-cage mini-protein to atomic
  resolution: Searching parameter space for a free-energy minimum".
  J. Chem. Phys. ??? (2006)

.. _rush_overview:

Usage overview: The steps involved in using the RUSH force field
----------------------------------------------------------------

The basic steps to use the RUSH force field are:

1) read in RUSH-modified versions of the CHARMM22 protein topology and 
   parameter files and the RUSH cmap parameter file (the appropriate
   non-bond list values are set in the parameter file)
2) build the appropriate PSF and read in / build the coordinates
3) use the RUSH keyword to initialize the subroutines and specify
   parameters other than the defaults
4) use "skipe all excl bond angl urey dihe impr rrep rpho rhbn rbdo rbac 
   raro cmap" to turn on the proper combination of energy terms
5) do energy, minimization (first derivative methods _only_), and
   dynamics calculations

To turn off RUSH and use the CHARMM22 force field:

1) use "RUSH off" to reset the subroutines and release any memory
   used by RUSH
2) use "SKIPE none" to reactivate Lennard-Jones, electrostatic, etc
   energy terms
3) delete all atoms
4) read in the original CHARMM22 force field topology and parameter 
   files to replace the RUSH topologies and parameters and non-bond
   parameters
5) proceed as usual

.. _rush_syntax:

Syntax of the RUSH command
--------------------------

Initialization:

::

    RUSH [PHOB <real>] [HBND <real>] [BDON <real>] [BACC <real>] -
         [KARO <real>] [CARO <real>] [DRMX <real>]

Clean-up:

::

    RUSH OFF

.. _rush_description:

Description
-----------

======== ======== ==========================================================
Keyword  Default  Purpose
======== ======== ==========================================================
PHOB     0.00072  multiplier for the hydrophobic surface-area term 
                  (kcal/(mol*angstrom**7/2)

HBND     -5.1     intra-molecular hydrogen bond strength (kcal/mol)

BDON      0.04    burial penalty for hydrogen-bond donors (kcal/mol)
                  
BACC      0.04    burial penalty for hydrogen-bond acceptors (kcal/mol)
                  
KARO      0.00    multiplier for force-shifted Coulomb attraction between
                  aromatic moieties (unitless)

CARO      0.00    cutoff for force-shifted Coulomb attraction between
                  aromatic moieties (angstrom)
                  
DRMX      0.50    minimum atomic motion for update of RUSH-associated
                  neighbor-lists (angstrom)

OFF       n/a     release heap used by RUSH, reset the RUSH-associated
                  global variables (including energies), and skip RUSH 
                  energy terms during subsequent energy/derivative
                  calculations
======== ======== ==========================================================

.. _rush_restrictions:

Restriction
-----------

The following will not work or give incorrect results:

- anything that requires 2nd derivatives
- periodic boundary conditions, images, etc.
- the free energy modules (pert, tsm, block)
- parallel CHARMM

.. _rush_notes:

Notes
-----

The hydrogen masses are set to 12.01 amu, a holdover from
early in the development of the force field when a single-sided
harmonic potential was used to model volume exclusion ( as opposed
to the current implementation, which employs the repulsive part
of the Weeks-Chandler-Andersen decomposition of the CHARMM22 
Lennard-Jones term) and required the increased hydrogen mass to
allow for energy conservation with a 2-fs timestep. This does not
affect the thermodynamic properties of the system owing to the lack
of atomic mass in the configurational part of the partition function.
The dynamic properties are affected, but this is a moot point since
the solvent is modeled implicitly. From a practical perspective, using
SHAKE is optional: a 2-fs timestep will achieve energy conservation
with or without SHAKE owing to the heavy hydrogens. From a philosophical
perspective, some may choose to use SHAKE since hydrogen has significant
quantum-mechanical character at room temperature. All development and
testing was done without SHAKE.

The surface-area based hydrophobic term is discontinuous in the first
derivatives. This occasionally causes problems with conjugate gradient
minimization; stick to steepest-descent. Energy conservation during
MD is not a problem despite the discontinuity. Using a 2-fs timestep,
the average total energy for a constant energy simulation of the trpzip2
designed beta hairpin after equilibration to 298 K is 631.7 kcal/mol
for the first 20-ps interval and 632.5 kcal/mol for the last 20-ps
interval of a 100-ps simulation.

PSF files from the standard CHARMM22 all-atom force field are NOT
compatible with RUSH. Use the rush topology file to build the PSF.

The aromatic energy term RARO was added subsequent to the publication of
the original model. Its purpose is to maintain the proper geometry of
interacting aromatic residues (Guvench and Brooks. "Tryptophan side chain
electrostatic interactions determine edge-to-face vs parallel-displaced
tryptophan side chain geometries in the designed beta-hairpin tripzip2".
J. Am. Chem. Soc. 127:4668-4674 (2005)) It is simply a force-switched
electrostatic term (Steinback and Brooks. "New spherical-cutoff methods
for long-range forces in macromolecular simulation". J. Comput. Chem. 15:
667-683 (1994)) that is applied to all atoms in the topology file with
non-zero charge, which currently comprise only the atoms of the
aromatic moieties of the Phe, Trp, and Tyr sidechains and the Tyr -OH
group. The default for KARO is 0.0, so that this term is not included
in the calculation.

.. _rush_examples:

Examples
--------

::

   * 1) read in the RUSH topology, parameter, and cmap files
   * 2) read in the protein psf and coordinates
   * 3) turn on RUSH
   * 4) do a minimization and energy calculation
   * 5) turn off RUSH and exit
   *

   open unit 1 read form name @TOPPAR/rush/top_rush_058.inp
   read rtf card unit 1
   close unit 1

   open unit 1 read form name @TOPPAR/rush/par_rush_058.inp
   read param card unit 1
   close unit 1

   open unit 1 read form name @TOPPAR/rush/par_rush_058_a.cmap
   read param append card unit 1
   close unit 1

   open unit 1 read form name ./@PDB.psf
   read psf card unit 1
   close unit 1
        
   open unit 1 read form name ./@PDB.pdb
   read coor pdb unit 1
   close unit 1

   rush

   skipe all excl -
    bond angl urey dihe impr rrep rpho rhbn rbdo rbac raro cmap

   minimize sd nstep 1000 tolgrd 0.1

   energy

   rush off

   stop
