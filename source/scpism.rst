.. py:module::scpism

===========================================================
Screened Coulomb Potentials Implicit Solvent Model (SCPISM)
===========================================================


The SCPISM is a continuum model for solvated macromolecules that treats
implicitly the effects of water. At the current stage of development the model
only incorporates a continuum description of electrostatic effects and,
optionally, a cavity term to account for hydrophobic interactions. The
electrostatic component of the model is based on a superposition of screened
Coulomb potentials (SCP) as derived from basic theory of polar liquids [2,5,6].

The model corrects hydrogen bond (HB) energies between all amino acid
pairs (sidechain-sidechain, sidechain-backbone, and backbone-backbone) due to
the competing effects of water molecules in their hydration shells. This HB
correction is treated independently of the electrostatics [3,4].

The current implementation is suitable for dynamics (plain or Langevin)
simulations, energy evaluation and minimization of peptides and proteins [1].
The models can be used with the all-atom representation (either PAR22 or CMAP)
or (since version c35) with PAR19. 
 
In the current implementation there is one parameter per atom type. 
The parameterization was done based on hydration energies of amino acid side
chain analogs and restrictions imposed to the space of parameters, as reported
[2] (the parameterization is not based on reproducing PB data. However, tests
have shown that the SCP model correlates well with results from PB
calculations [2,7]).


.. _scpism_syntax:

SCPISM commands
---------------

An effort has been made to minimize the number of input options available 
to the user (i.e., no parameters can be modified from the input file since
the physics of the system is already incorporated into the model and then
hardwired into the algorithm). To activate the model the following command
line is used:

::

   SCPIsm UISM [int]                                                 (1)

where UISM is the unit number for reading the SCP parameters. These
parameters are stored in scpism.inp and must be opened for reading
in the usual way before the model is requested (see examples below), i.e.,

::

   OPEN READ UNIT int CARD NAME "scpism.inp"

Once the model is requested, any options for energy, minimization 
and dynamics calculations are supported. Electrostatic interactions
are truncated using a shift function (any other specification of cutoff
of electrostatics is automatically disabled once the SCPISM is requested
(a shift function was shown to be the best option for the functional form
used in the SCPISM).

Any restraining option that is introduced as an additional term in the 
potential is supported. Use of the CONS FIX constraints is discouraged
because they are inconsistent with the way the self-energies are
calculated: the calculation of Born radii in the SCPISM requires the
positions of all the atoms around each atom to be available, in particular
the positions of the atoms to be fixed.  

It is possible (but not required) to deactivate the model once any 
of the required tasks (energy evaluation, minimization or dynamics) is 
completed. To exit the model use

::

   SCPIsm END

in this case CHARMM will return to the default (vacuum) electrostatics  
options (the model can be turned on again by using the command line 
(1) above). 

Note that, although the current implementation of the SCPISM describes
only electrostatic effects, to be consistent with earlier applications
of the model, an overall hydrophobic term proportional to the total
solvent accessible surface area (SASA), also based on a contact model 
approach, can be requested by using the subcommand HYDRophobic in the
command line (1) above (see examples). If this term is not activated, 
other model accounting for hydrophobic interactions should be used.

The CPU time of the serial version of the SCPISM is between 2 and 3 times
slower than vacuum (depending on platform and compilation optimization;
see [8] 

for a comparison with other continuum models in CHARMM). From version c34 
the models has been parallelized (MPI routines added to the scpism.src source
code). The best current performance is for 2 or 4 processors (this should be
improved in the future)

A major change has been implemented in the developmental version c35b1 to 
account more accurately for HB energies. These HB energies has been estimated
through systematic calculations of potentials of mean force (PMF), as reported
[4]. These PMF were obtained with explicit TIP3P water molecules. To use this
new implementation the SCPISM must be requested by the following command (see 
EXAMPLE 2 below):

::

   SCPIsm UISM (int) UIHB (int)                                            (2)

which is similar to the syntax in (1) but a new unit number UIHB allows to read
parameters for the HB corrections. Prior to c35 a simplified version of these
parameters was hardwired into the code, but (2) allows: i) more flexibility to
adjust specific HB interactions if other models other than TIP3 are used to
calculate PMF), and ii) HB interactions between amino acid pairs to be exactly
reproduced according to the PMF reported in [4]. (Note: To use this new HB
parameters in c35b1 the variable NEWOLD in scpism.src must be set to 'YN'
instead of the default 'NY'; the latter allows to use the OLD HB
implementation; the former allows to use the NEW implementation)

Two files must now be opened:

::

   OPEN READ UNIT int CARD NAME "scpism.inp"
   OPEN READ UNIT int CARD NAME "scphb.inp"

scpism.inp contains parameters for the electrostatics (same as previous
versions), and scphb.inp contains parameters for the HB corrections. In c35b1
parameters for both par22 and par19 are available, so these files are actually
scpism-par22.inp and scphb-par22.inp, and scpsim-par19.inp and scphb-par19.inp,
respectively. As in previous versions a hydrophobic term can be requested with
the keyword HYDR. (Note: For consistency between par19 and par22, parameters
in scpism-par22.inp are slighlty different than in scpism.inp of previous
versions; experimental hydration energies are however not affected)

.. _scpism_background:

Structure of the parameter input files
--------------------------------------

All electrostatic parameters required for the model are atom-type based
and are collected in a single file scpism.inp. Determination of these
parameters was described in [1,2,3].

The scpism.inp parameter file has the following format (note that not all of
the fields in this file are parameters for electrostatics):

::

   H   0.4906  0.6259  0.5000  0.7004  0.3700  0.0052 PH ! polar H
   HC  0.5300  9.6746  0.5000  0.7280  0.3700  0.0052 PH ! N-ter H
   HA  0.5300  0.5000  0.5000  0.7280  0.3700  0.0052    ! nonpolar H
   HT  0.4906  2.3800  0.5000  0.7004  0.3700  0.0052 PH ! TIPS3P WATER HYDROGEN
   HP  0.5274  0.5000  0.5000  0.7262  0.3700  0.0052    ! aromatic H
   HB  0.4858  0.5000  0.5000  0.6970  0.3700  0.0052    ! backbone H
   HR1 0.4651  0.5000  0.5000  0.6820  0.3700  0.0052    ! his he1, (+)his HG,HD2
   HR2 0.4580  0.5000  0.5000  0.6768  0.3700  0.0052    ! (+) his HE1

   Col. 1: Atom type defined in PAR22
   Col. 2: Alpha_i controls slope of D(r) around atom-type i
   Col. 3: for PH interaction with PA this value controls the Born radius
           of PH to modulate hydrogen bonding strength; the increase or descrease 
   	of the PH Born radius is defined by the product Col.3(PH)*Col.3(PA)
   	[2,3]
   	As a rule, the larger the value of this product, the weaker the HB
   	interaction.
   Col. 4: Extension of Born radius R_iw to obtain R_ip, i.e.,
           R_ip = R_iw + Col.4 [2,3]; only one value for atoms considered
   Col. 5: SQRT(alpha_i); this is needed for alpha_ij=alpha_i alpha_j [2,3]
   Col. 6: covalent radius for each atom (only 5 values considered, i.e.,
           for C,N,O,S,H)  
   Col. 7: gamma_i in hydrophobic energy = summatory over i of gamma_i SASA_i
           (note that all coefficients are equal in this first release) 
   Col. 8: denotes atoms involved in H-bonding (PH = polar proton; PA =
           proton acceptor) 

The scphb.inp parameter file has the following format:

::

   N   AD    R    A'D'

*N* is the index of the correspoding amino acid pair.

*AD* denotes the amino acid pair, with *A* the acceptor and *D* the donor
(one-letter nomenclature except j, which denotes the unprotonated form of His);
a letter 'b' means that the corresponding acceptor or donor is a backbone atom;
last entry 'bb' corresponds to backbone-backbone interactions.

*R* is the parameter that controls the particular HB strength 

*A'D'* is an amino acid pair not explicitly used in the PMF calculations but
assumed to have same HB energy as *AD* 


Theoretical background
----------------------

The SCPISM is based on a superposition of screened potentials. Screening
functions are derived from the Lorentz-Debye-Sack theory [2,5,6]. The model
uses screenings D(r) that modulate the potential phi(r) rather than dielectric
functions eps(r) that modulate the electric field E(r). The relation between
both functions is given by the definition of potential, E(r) = -Grad[phi(r)]

Both D(r) and eps(r) are signoidal functions of r. Once eps(r)
is known from theory or experiments. D(r) is obtained by integration [2,5]
Based on these results, the SCPISM uses atom type-dependent sigmoidal
functions. 

In the SCPISM the standard electrostatic component of the force field
is replaced by terms that describe both the electrostatic interaction energy, 
and the self energy. The screening functions are continuous functions of the 
position and describe a dielectric medium that permeates all of space. For the 
solvated protein, D(r) approaches bulk screening only far from the protein
(see [5,6,7] for discussion). Therefore, the SCPISM does not introduce a
priori an internal and an external dielectric constant and, then, there is no
boundary that separates the protein from the solvent, as used in true
macroscopic electrostatics. At the microscopic scale of biomolecules the
physically meaningful way to account for this in a continuum approximation is
through the position-dependent density of water [6]. The model is being
constructed incrementally to incorporate all the physical effects that are
removed when the explicit solvent is removed. 

Because the effective screening functions that characterize the overall
modulation of the electrostatics are obtained from properties of bulk solvent, 
the short-range (typically < 5 Ang.) interactions that characterize the HB
strengths must be corrected [3,4]. This is so because the forces between two
solutes in close proximity are controled not only by bulk electrostatic but
also by the structure of water around the solutes (see discussion in [6]).
To obtain a reasonable representation of HB strength, each and all HB
interactions between amino acids pairs are calibrated individually via the
self-energy terms (by adjusting the Born radii of the shared polar H [3]).
For the current implementation (suitable for dynamic simulations) this
approach was simplified with respect to the original development (for MC
simulations) [2,3], so the HB interactions depend only on the atom type,
while no directionality has been included. However, the geometry of HB was
shown to be important in MC for structure calculation to avoid artifacts
such as multiple HB for a single donor or acceptor atom. Therefore, MC
simulations are not recommended with the current implementation of the model.


.. _scpism_references:

References
----------

[1] S A Hassan, E L Mehler, D Zhang and H Weinstein, Molecular Dynamics
    Simulations of Peptides and Proteins with a Continuum Electrostatic Model
    based on Screened Coulomb Potentials; Proteins 51, 109 (2003)

[2] S A Hassan, F Guarnieri and E L Mehler, A General Treatment of Solvent 
    Effects based on Screened Coulomb Potentials; J Phys Chem B 104, 6478 
    (2000)

[3] S A Hassan, F Guarnieri and E L Mehler, Characterization of Hydrogen 
    Bonding in a Continuum Solvent Model; J Phys Chem B 104, 6490 (2000)

[4] S A Hassan, Intermolecular Potentials of Mean Force of Amino Acid Side
    Chain Interactions in Aqueous Medium; J Phys Chem B 108, 19501 (2004)

[5] S A Hassan and E L Mehler, From Quantum Chemistry and the Classical
    Theory of Polar Liquids to Continuum Approximations in Molecular Mechanics
    Calculations; Int. J. Quant. Chem. 102, 986 (2005)

[6] S A Hassan, Liquid Structure Forces and Electrostatic Modulation of
    Biomolecular Interactions in Solution; J. Phys. Chem. B 111, 227 (2007)

[7] S A Hassan and E L Mehler, A Critical Analysis of Continuum Electrostatics:
    The screened Coulomb potential-implicit solvenmt model and the study of
    the alanine dipeptide and discrimination of misfolded structures of
    proteins; Proteins 47, 45 (2002)
 
[8] B R Brooks et al., CHARMM: the biomolecular simulation program; 
    J. Comp. Chem. 30, 1545 (2009) 


.. _scpism_example:

SCPISM example input files
--------------------------

EXAMPLE 1 (versions c34 and previous)

::

   * energy, minimization and dynamics with the SCPISM
   *

   open read card unit 10 name top22.inp
   read rtf unit 10 card
   close unit 10 

   open read card unit 10 name par22.inp
   read para unit 10 card 
   close unit 10 

   open read unit 10 card name SCPISM.inp

   ! define system
   generate main setup

   open read unit 2 card name filename.crd
   read coor card unit 2 
   close unit 2

   ! set up all options for the run (nbond cutoff distance, shake, vdW, etc)
   ! these options can also be specified in the 'energy' command line, as usual

   ! request SCPISM 

   SCPI HYDR UISM 10

   ! this will activate electrostatics and a simple hydrophobic term;
   ! if other model is used to describe hydrophobic interactions then
   ! replace the above command line by: SCPI UISM 10

   ! now calculate energy, minimize structure and run dynamics

   ener 

   mini [options] 

   dyna [options]

   ! exit SCPISM

   SCPI END

   ! if needed, the model can be requested again at any point

   stop


EXAMPLE 2 (from developmental c35b1 onward; for c35b1 changed NEWOLD to 'YN')

::

   * energy, minimization and dynamics with the SCPISM
   *

   open read card unit 10 name top22.inp
   read rtf unit 10 card
   close unit 10

   open read card unit 10 name par22.inp
   read para unit 10 card
   close unit 10

   ! below open the scp and hb parameters files
   ! if par19 is used instead of par22, then open the corresponding files

   open read unit 10 card name scpism-par22.inp
   open read unit 11 card name scphb-par22.inp

   ! define system
   generate main setup

   open read unit 2 card name filename.crd
   read coor card unit 2
   close unit 2

   ! set up all options for the run (nbond cutoff distance, shake, vdW, etc)
   ! these options can also be specified in the 'energy' command line, as usual

   ! request SCPISM

   SCPI HYDR UISM 10 UIHB 11

   ! this will activate electrostatics and a simple hydrophobic term;
   ! if other model is used to describe hydrophobic interactions then
   ! replace the above command line by: SCPI UISM 10 UIHB 11

   ! now calculate energy, minimize structure and run dynamics

   ener

   mini [options]

   dyna [options]

   ! exit SCPISM

   SCPI END

   ! if needed, the model can be requested again at any point

   stop


