.. py:module::sccdftb

=====================================================================================
Combined Quantum Mechanical and Molecular Mechanics Method Based on SCCDFTB in CHARMM
=====================================================================================

by  Qiang Cui and Marcus Elstner
(cui@chem.wisc.edu, elstner@phys.upb.de)

The approximate Density Functional program SCCDFTB (Self-
consistent charge Density-Functional Tight-Binding) is interfaced with
CHARMM program in a QM/MM method.  

This method is described in 

``Phys. Rev. B  58 (1998) 7260``,
``Phys. Stat. Sol. B 217 (2000) 357``,
``J. Phys. : Condens. Matter.  14 (2002) 3015``.

The QM/MM interface in CHARMM has been described in 
``J. Phys. Chem. B 105 (2001) 569``

The GHO-SCC-DFTB/MM boundary treatment has been described
in ``J. Phys. Chem. A 108 (2004) 5454``.

The extension of the SCC-DFTB method to work with the Replica Path 
and the Nudged Elastic Band methods has been described in the following
paper and should be cited when applied... 

H. L. Woodcock, M. Hodoscek, and B. R. Brooks Exploring SCC-DFTB Paths
for Mapping QM/MM Reaction Mechanisms J. Phys. Chem. A; 2007; 111(26)
5720-5728.


.. _sccdftb_description:

Description
-----------

The SCCDFTB QM potential is initialized with the SCCDFTB command

::

   [SYNTAX SCCDFTB]

   SCCDFTB   [REMOve] [CHRG] (atom selection) [GLNK atom-selection]
             [TEMPerature] [SCFtolerance] 
             [CUTF] [EWAD EOPT KAPPA kappa KMAX kmax KSQMAX ksqmax]
             [UPDT 0]
             [MULL]
             [DISP]
             [HBON GAUS]
             [DHUB DHGA]
          

   REMOve:  Classical energies within QM atoms are removed.

   CHRG:    Net charge in the QM subsystem.

            The atoms in selection will be treated as QM atoms.

   GLNK atom-selection: contains a list of atoms that are selected as
                        GHO boundary atoms.

   TEMPerature:  Specifies the electronic temperature (Fermi distribution).
                 Can be used to accelerate or achieve SCF convergence 
                 (default =0.0).
              
   SCFtolerance: Convergence criteria for the SCF cycle. As default
                 a value of 1.d-7 is used.

   CUTF: a flag that turns on cut-off, will use the same scheme used for
         MM interactions (allows atom-based fshift and fswitch)

   EWAD: a flag that turns on ewald summation for QM-QM and QM-MM 
         electrostatic interactions
   EOPT: performs an internal optimization for kappa and kmax as well as
         real-space sum. NOT recommended - very inefficient. incompatible
         with CUTF

   KAPPA, KMAX, KSQMAX: parameters used for QM-QM and QM-MM ewald 
         contributions.

   UPDT: Whether to update the box during a MD run. Default is 0 (not update!)
         Do NOT forget to specify this to 1 when running CPT calculations

   MULL: Transfer the mulliken charges to the CG array such that they can be
         printed out by transferring the CG array to any vectors (e.g.,
         scalar wmain = charge; print coor)

   DISP: Dispersion interactions among QM atoms can be calculated using an 
         empirical formular (Elstner et al. J. Chem. Phys. 114, 5149, 2001).
         One needs to specify a set of parameters in the DISPERSION.INP file
         (see the above ref for details).

   HBON: The short-range behavior of XH gamma function is modified with a damping 
         function. This significantly enhances the hydrogen bonding interactions
         (Elstner, Cui, unpublished). The dampling exponent needs to be specified
         in the sccdftb.dat file. However, with the modified gamma, the repulsive
         potential has to be adjusted accordingly. This can be done in an
         empirical fashion by including a Gaussian (constrained to operate in
         a range by a switching function) in the relevant repulsive potential;
         the parameters for the Gaussian and the switching function need to be
         specified in the spl file and turned on using the GAUS keyword.

   DHBU: To improve proton affinities, which depend much on the charges in the
         protonated DHGA and deprotonated molecules, the SCC-DFTB is expanded
         to the 3rd order. Currently only on-site terms have been included,
         which were observed to have a major impact on the calculated PAs for
         many molecules. The relevant parameters are the derivative of the
         Hubbard parameters (related to chemical hardness), which need to be
         specified in the sccdftb.dat file. The DHGA keyword includes further
         flexibility in the behavior of the Hubbard parameters as a function of
         charge (Elstner, Cui, unpublished).

In the SCCDFTB program the atomtypes are represented by consecutive 
numbers. The definition of SCCDFTB atom numbers has to be accomplished 
before invoking the SCCDFTB command. The numbers are stored in WMAIN.
If the QM system e.g contains only O, N, C and H atoms,
the the numbering can be executed as follows:

::

   scalar WMAIN set 1.0 sele type O*  SHOW end
   scalar WMAIN set 2.0 sele type N*  SHOW end
   scalar WMAIN set 3.0 sele type C*  SHOW end
   scalar WMAIN set 4.0 sele type H*  SHOW end

Now, the O atoms are represented by 1.0, the N atoms by 2.0 etc. 

Link atom may be added between an QM and MM atoms with the
following command:

::

   ADDLinkatom  link-atom-name  QM-atom-spec  MM-atom-spec

         link-atom-name ::= a four character descriptor starting with QQ.

         atom-spec::= {residue-number atom-name}
                      { segid  resid atom-name }
                      { BYNUm  atom-number     }

When using link atoms to break a bond between QM and MM
regions bond and angle parameters have to be added to parameter file
or better use READ PARAm APPEnd command.

If define is used for selection of QM region put it after all
ADDLink commands so the numbers of atoms in the selections are not
changed. Link atoms are always selected as QM atoms.

.. _sccdftb_usage:

SCCDFTB input files
-------------------

SCCDFTB needs to read in the parameter files, which have 
a two-body character. Therefore, the interaction parmeters
for all pairs of atoms have to be read in.
These files are named like oo.spl, on.spl, oc.spl, no.spl etc.,
where oo.spl contains the two-center integrals for the O-O interaction,
on.spl the  two-center integrals for the O-N interaction etc.
DFTB needs these parameters for the O-N and N-O interaction,
similarily for all other pairwise interactions.
The file sccdftb.dat contains the paths to these parameters, as:

::

   'potential:atom-1-atom-1'
   'potential:atom-1-atom-2'
   'potential:atom-1-atom-3'
   ... \\
   'potential:atom-1-atom-N' 
   'potential:atom-2-atom-1' 
   ... \\
   'potential:atom-2-atom-N' 
   'potential:atom-N-atom-1' 
   ... \\
   'potential:atom-N-atom-N'

where atom-1 is the atom defined by 1.0, as described above,
atom-2 defined in WMAIN by 2.0 etc.

For the example of the system containing O N C and H, sccdftb.dat would
contain:

::

   'PATH/oo.spl'
   'PATH/on.spl'
   'PATH/oc.spl'
   'PATH/oh.spl'
   'PATH/no.spl'
   'PATH/nn.spl'
   'PATH/nc.spl'
   'PATH/nh.spl'
   'PATH/co.spl'
   'PATH/cn.spl'
   'PATH/cc.spl'
   'PATH/ch.spl'
   'PATH/ho.spl'
   'PATH/hn.spl'
   'PATH/hc.spl'
   'PATH/hh.spl'

where PATH specifies the path to the directory where the data files 
are located. Be careful, an error in the sequence or a wrong assingnment
 of parameters to atoms (coordinates) will make results meaningless.
Parameter files can be requested from Marcus Elstner 
(elstner@phys.upb.de).


SCCDFTB output files (currently disabled)
-----------------------------------------

::

   SPE.DAT : contains the Kohn-Sham energies with occupations numbers.
   CHR.DAT : contains the atomic  (Mulliken) charges of the atoms 
             (first row) and for the orbitals (s, px,py,pz,dxx.. ) in the 
             following columns.
   REST.DAT: contains dipolemoment (D), calculated from the Mulliken 
             charges (not a reliable estimate of Dipolemoment in general!)

.. _sccdftb_installation:

Installation of SCCDFTB
-----------------------

The source code of SCCDFTB ist distributed with CHARMM.
To compile the SCCDFTB method as the quantum part:

::
   
   ./install machine size T

   T invokes the SCCDFTB 
   
The parameter files have to be reqeusted and stored in a directory,
which can be reached by 'PATH' (see up).

Diagonalization routines
------------------------

As default, the library routine dsygv.f (LAPACK) is used for 
the diagonalization of the hamiltonian matrix.
This is called by ``chmdir/source/scctbint/scctbsrc/ewevge.f``.
A faster (about factor 2) solution is given by the dsygvd.f routine 
(but less stable), which is called by ewevge-dsygvd.f: 
copy ewevge-dsygvd.f to ewevge.f and recompile to invoke this option.
Contact Marcus Elstner for more details or questions. 


.. _sccdftb_fep:

Free energy perturbations with SCC-DFTB/MM
------------------------------------------

The code currently allows dual-topology based SCC-DFTB/MM free energy
perturbation calculations; since all scaling related to the QM component
of the free energy derivative is done inside SCC-DFTB, the FEP 
calculations do not have to use BLOCK.

As discussed in JPC, 107, 8643 (2003), a practical problem of using
FEP with QM/MM potentials is that the structure of the QM region 
undergoes significant distortions at end-points if one scales the entire
QM molecule; there is no such problem if one chooses to scale only
QM/MM interactions, but that requires calculation of new terms. The
general solution is to add harmonic constraints on the QM part -
either only at the end-points and then re-weight the calculated free 
energy derivatives - or, more elegantly, add harmonic constraints 
as "chaperones" throughout the "alchemy" simulation and compute 
corrections based on local configuration integrals. See W. Yang et al. 
J. Chem. Phys. 2004.

For the special case where the two end-states have very similar chemical
structures - such as in redox, metal-exchange and pKa applications, 
which we believe are scenarios where QM/MM treatment is useful, a simple
dual-topology-single-coordinate (DTSC) approach has been introduced. As
the name implies, one uses only one set of coordinates for the two 
states (e.g., reduced and oxidized states). Due to the fact that the 
free energy is path-independent, such an approach is formally exact. In 
practical applications, error might arise due to SHAKE - i.e. X-H 
distances are assumed to the same in the two states - which usually has 
negligible effects.

At each configuration (hence single-coordinate) along the trajectory, 
two electronic structure calculations are carried out (dual topology) 
and the free energy derivative with respect to the coupling parameter 
is evaluated and averaged on the fly. 

With minor modifications, the algorithm also works for pKa prediction
for a specific group in large molecules. For more details, refer to the
following publications:

-	M. Formaneck, G. Li, X. Zhang, Q. Cui, J. Thero. Comput. Chem.
	1, 53-68 (2002) 
-	G. Li, X. Zhang, Q. Cui, J. Phys. Chem. B (2003) 107, 8643
-	G. Li, Q. Cui, J. Phys. Chem. B, (2003) 107, 14521

NOTE BENE: It MUST be used with "FAST OFF" because only generic 
atom-atom codes have been modified so far (made default).

Due to the fact that ALL QM related components are handled within SCC
(including GSBP and eWald, see next section), FEP (such as pKa) 
calculations can be used with both eWald and GSBP - provided that it is 
the QM part that undergoes "alchemical" mutation. 
The code has NOT been extensively tested in which both QM and MM 
undergo changes.

Two examples are given to illustrate computational details; the
first one deals with redox potential calculations for FAD in Cholesterol
oxidase, and the second one concerns pKa calculations of ethanethiol
(Ch3CH2SH) in water. The test files are scc_fep_dtsc.inp and 
scc_pka1/2.inp. 

Example 1
^^^^^^^^^

For redox potential calculations, the following set-up is used,

::

   ......

     SCCDFTB LAMDa [REST] STOP [CUTF] OUTPut int -
     [REMOve] (atom selection 1) [CHRG] [TEMPerature] [SCFtolerance] -
     INIT @lam PASS int STEP int TIAV int -
                                 [CHRG] [TEMPerature] [SCFtolerance] -
              (atom selection 2) -
              (atom selection 3) 

   LAMD:  invoke the TI method to perform free energy calculations
   REST:  Restart option for accumulating statistics concerning <DU/DL>
          i.e., necessary values will be read in from dynamics restart
          file
   STOP:  employ the dual-topology-single-coordinate approach
   CUTF:  invoke cutoff for QM/MM electrostatic interactions
   OUTP:  unit number for storing the free energy derivative <DU/DL>
   INIT:  the current lamda value
   PASS:  numbers of MD steps to be skipped when accumulating <DU/DL>
   STEP:  the frequency of collecting statistics for <DU/DL> 
   TIAV:  the frequency of computing the average of DU/DL.

   atom selection 1: Reactant+Product to set up MM list for QM atoms
   atom selection 2: Reactant state
   atom selection 3: Product  state

For pKa calculations, two free energy simulations are in 
principle required; in the first step, the protonated state is mutated
into the ionized state as the acidic proton is mutated into a dummy atom
in the second step, the dummy atom is transferred into the gas phase. 
Test calculations indicate that the contribution from the 2nd step is 
likely to be small. 

Example 2a
^^^^^^^^^^

In the first step, BLOCK is used together with SCCDFTB

::

   ......

   BLOCK 3
   SCCDFTB STOP PKAC ISTP 1
   CALL 2 SELE qm1  END
   CALL 3 SELE qmh  END
   CALL 1 SELE .not. (qm1 .or. qmh) END
   COEF 1  1  1.0
   COEF 1  2  1.0
   COEF 1  3  1.0
   COEF 2  2  0.0
   COEF 2  3  @lam
   COEF 3  3  0.0
   END

   SCCDFTB PKAC ISTP 1 HYGN int [CUTF] OUTPut int -
   [REMOve] (atom selection 1) [CHRG] [TEMPerature] [SCFtolerance] -
   INIT @lam PASS int STEP int TIAV int -
                               [CHRG] [TEMPerature] [SCFtolerance] -
   atom selection 2 -
   atom selection 3
   ......

In the BLOCK section :

::

   the SCCD keyword is used to set up coefficent matrix for calculating 
   bonded contribution involving the dummy atom to <DU/DL>. 
   ISTP 1: the first step in pKa calculations
   qm1 is the ionized state (e.g., CH3CH2S-); 
   qmh is the acidic proton.

In the SCCDFTB section:

::

   PKAC  : invoke pKa calculation
   ISTP 1: the first step in pKa calculations
   HYGN  : atomic index (number) of the acidic proton in the psf 
   atom selection 1:   protonated state (CHRG:   protonated state)
   atom selection 2:   protonated state (CHRG: deprotonated state)
   atom selection 3: deprotonated state

Example 2b
^^^^^^^^^^

In the second step for pKa calculations, the dummy atom is transferred 
into vacuum,

::

   ......
   calc 1mlam 1.0-@lam

   BLOCK 3
   SCCDFTB STOP PKAC ISTP 2
   CALL 2 SELE qm1  END
   CALL 3 SELE qmh  END
   CALL 1 SELE .not. (qm1 .or. qmh) END
   COEF 1  1  1.0
   COEF 1  2  1.0
   COEF 1  3  @1mlam
   COEF 2  2  0.0
   COEF 2  3  0.0 bond 1.0 angl 1.0 dihe 1.0
   COEF 3  3  0.0
   END

   SCCDFTB PKAC ISTP 2 HYGN int [CUTF] OUTPut int -
   [REMOve] [CHRG] (atom selection 1) [TEMPerature] [SCFtolerance] -
   INIT @lam PASS int STEP int TIAV int 
   ......

Note that with BLOCK, the coefficient matrix is different in the 
second step: we are only scaling the non-bond (vdW) interaction between
the environment and the dummy atom (1 and 3). The bonded terms between
the QM and the dummy atom (2 and 3) is kept (coefficient as 1.0) and
will be taken out with local configuration integrals.

In SCC-DFTB, atom selection 1: deprotonated state


.. _sccdftb_electrostatics:

Since 2004, electrostatics in SCC-DFTB/MM simulations can be 
treated in several ways for both spherical and periodic conditions:

i) As for other QM packages, the default is no cut-off for QM/MM 
   electrostatic interactions. This is NOT recommended when cut-off 
   is used for MM; the imbalance will cause over-polarization of the media
   (e.g.,see discussion in classical simulations by Woods, J. Chem. Phys. 
   103, 6177, 1995). A useful option is to use extended electrostatics for
   MM.

ii) Cut-off is introduced for QM/MM electrostatics, similar to MM
    interactions; i.e., the same scaling factors are the same as those for 
    MM interactions. Simply add "CUTF" to the SCC-DFTB command line.
    Currently only supports energy/force-shifts based on atoms

iii) For spherical boundary conditions, the
     GSBP approach can now be used with SCC-DFTB. The current implementation
     takes GSBP contributions into the SCF iteration, although for a large
     inner region, this may not be necessary. Further tests are being carried
     out. The code will be extended to other boundary conditions and QM 
     methods in the future. If one uses sorting (i.e., truncate size of
     basis in GSBP), make sure a SCC-DFTB/MM 
     energy calculation with MULL (save Mulliken charge) is carried out 
     before issuing GSBP, since the Mulliken charges are used to estimate
     contributions from various basis functions to the QM related terms. 
     See test cases for examples.

iv) For PBC simulations, one can use either cut-off or eWald sum for
    SCC-DFTB/MM interactions. No PME has been implemented for the QM/MM 
    interactions although it may not be too unreasonable to use PME for the
    expensive MM part and ewald for QM/MM interactions. The current QM/MM
    implementation allows in principle all cell shapes.
    For eWald, one can either let the code optimize the exponent to get the
    best balance between real space sum and the reciprocal space sum (EOPT)
    or one can specify a set of parameters (Kappa, KMAX, KSQMAX).
    The real space sum is done till convergence is met with EOPT or 
    without CUTF (so more expensive); EOPT is incompatible with CUTF.
    With cutoff (CUTF), the real sum is limited to atoms within the cutoff- 
    which is recommended (much more efficient). In any case, one should
    carefully test kappa, KMAX to ensure the convergence of energy and, 
    more importantly, force from SCC-DFTB/MM calculations. 

A sample command line would be:

::

   ......
   SCCDFTB remove CHRG 2 SELE resn @m END TEMP 0.00 SCFT 0.00000001 EWAD -
   CUTF Kappa 0.45  KMAX 6 KSQMAX 100
   ......

The eWald code is not as efficient as one might hope for at this stage.
Typically QM/MM-eWald is about 5-8 times slower than a QM/MM calculation
without eWald. 

An important point for PBC simulations is that all image must be used 
with "UPDAte IMAL ". This is because symmetry operations have not been
considered in the SCC-DFTB/MM code - which obviously needs to be fixed
in the future.

.. _sccdftb_status:

Status
------

The current implementation has analytical first derivative and thus
allows energy minimizations, reaction path search (e.g., travel) and
molecular dynamics simulations; SCC-DFTB/MM also works with Monte Carlo.
Replica can also be used, which makes it possible to use replica path
and related approaches (such the nudged elastic band) for determining
reaction path with the SCC-DFTB/MM potential; along the same line,
path integral simulations can be carried out as well, although only for
equilibrium properties at this stage.

Several aspects of the code will be improved in the near future, 
and new functionalities will be added:

1. Interface with centroid path-integral simulations and Tsallis
   statistics.
2. More flexible interface with BLOCK for general free energy 
   simulations.
3. Better methods for open-shell systems; constrained density functional
   theories.
4. Time-dependent treatment for electronically excited states; non-adiabtic MD.
5. Integration with polarizable force field models (Drude).

.. _sccdftb_glnk:

Description of the GLNK Command
-------------------------------

::

   [GLNK atom-selection]

   atom-selection: contains a list of atoms that are boundary atoms.

   Restrictions: see the correponding entry for GLNK in quantum.doc

   Description:  see the correponding entry for GLNK in quantum.doc

Limitations: The present implementation allows up to 5 QM-boundary
atoms. To improve the geometry for the QM/MM boundary bond, an
empirical correction (Ecor) term is added. Currently, Ecor parameters
are only available for cases where the QM/MM partition cuts a C-C,
a C-O, or a C-S bond.  For other cases, no empirical corrections
will be included. Unrestricted GHO-SCC-DFTB for open-shell system
is not implemented.

Reference: Reference made to the following paper, which contains
a more thorough description and discussion of test cases, is appreciated.

Jingzhi Pu, Jiali Gao, and Donald G. Truhlar,
J. Phys. Chem. A 108, 5454-5463 (1998). "Combining Self-Consistent-Charge
Density-Functional Tight-Binding (SCC-DFTB) with Molecular Mechanics by
the Generalized Hybrid Orbital (GHO) Method."
