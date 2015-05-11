.. py:module:: mrmd

=============================================================
Multi-Surface Adiabatic Reactive Molecular Dynamics (MS-ARMD)
=============================================================

::

        Multi-Surface Adiabatic Reactive Molecular Dynamics (MS-ARMD)

  by 
  Tibor Nagy           (tibornagy@chem.elte.hu)  [coder]
  Juvenal Yosa Reyes   (juvenal.yosa@unibas.ch) 
  Markus Meuwly.       (m.meuwly@unibas.ch)

  References:
  [1] Multi-Surface Adiabatic Reactive Molecular Dynamics 
      Tibor Nagy*, Juvenal Yosa Reyes, Markus Meuwly*
      submitted to JCTC (Oct, 2013).
  [2] State-selected ion-molecule reactions with Coulomb-crystallized 
      molecular ions in traps
      Xin Tong, Tibor Nagy, Juvenal Yosa Reyes, Matthias Germann, 
      Markus Meuwly*, Stefan Willitsch*
      Chemical Physics Letters, Volume 547, 21 September 2012, Pages 1�8
       
  The Multi-Surface Adiabatic Reactive Molecular Dynamics (MS-ARMD) method allows 
  the construction of global reactive potential energy surface from standard 
  force fields and running molecular dynamics or Monte Carlo simulations on it. 
  The effective surface is always the lowest energy surface, except for the region
  where several surfaces have the same low energy, where it switches smoothly 
  between them by changing their weights. The algorithm is based on an energy 
  difference-based switching method, which conserves total energy during dynamics.
  The code can be run also with a single state allowing the easy usage of some 
  advanced parametrization (Morse potential, MIE potential) also for non-reactive 
  simulations. The CROSS module of CHARMM is based on the ARMD method, which uses
  an alternative, time-dependent switching function. The ARMD method has been used
  successfully for the simulation of reactions in large molecules (see cross.doc 
  for references). The comparison of MS-ARMD and ARMD methods and their advantages
  features are discussed in detail in reference [1]. Looking at the example 
  provided as a test case with CHARMM (mrmd_h2so4.inp with mrmd_h2so4.par MRMD
  parameter file, see Example Section) can help in understanding this 
  documentation greatly.

  * Menu:

  * Syntax::                
                            Syntax of the MRMD command
  * Description::           
                            Description on the MRMD command
  * Parameter file::         
                            Description of multiple states and their connection 
  * States::
                            Declaration of states: leveling and their connection
  * Nonbonded::
                            (Re)parametrization of nonbonded interactions for 
                            each state
  * Bonded::
                            (Re)parametrization of bonded interactions for each 
                            state
  * Output::                   
                            Detailed description of output generated during normal 
                            termination
  * Example::
                            Water elimination from vibrationally highly 
                            excited sulfuric acid molecule                         
  * Troubleshooting::
                            Hint for troubleshooting
  * Code::
                            Essential code related notes

  
  File: MRMD, Node: Syntax,Previous: Top,Up: Top, Next: Description

            Syntax for the Reactive Molecular Dynamics commands

                      |---------------optional arguments---------------| 
  MRMD   UPAR integer [ UCRG integer ] [ PRDY integer ] [ PRCA integer ]
  MRMD   RSET

  UPAR    Int    Unit containing the mrmd parametrization for all surfaces.
                 This unit must be opened (FORMatted) for reading before MRMD
                 is called.
  UCRG     -1    Unit to which a geometry will be written (FORMatted write) in PDB 
                 format at each crossing point (default=-1 => no writing). 
  PRCA     -1    Period in module calls at which MRMD prints out the energy  
                 and the weight of each surface and the energy of the effective
                 surface when dynamics is not active (default=-1 => no writing).
  PRDY     -1    Period in time-steps at which MRMD prints out the energies  
                 and the weights of each surface and the energy of the effective
                 surface during dynamics. (default=-1 => no writing)
  RSET           Deactivate currently active MRMD parametrization. Should not be
                 used together with other MRMD keywords.
                 
  
  File: MRMD, Node: Description, Previous: Syntax, Up: Top, Next: Parameter file

          Description of the Reactive Molecular Dynamics command

  -----------------------------
  MRMD command in CHARMM input:
  -----------------------------
  MS-ARMD force field is invoked by the MRMD command. On the first call the 
  content of MRMD parameter file is fully processed. This contains information on 
  one or more force fields by providing the changes from the PSF/parameter file 
  by means of modifying/removing/adding FF potential terms. The PSF and parameter 
  arrays in CHARMM are never modified! MRMD will affect only the total energy 
  and the forces on individual atoms. The parameter file also contains information
  on the switching function, that is how the surfaces are connected.

  ----------------------------------------------------------------------------
  Position of the MRMD command during input relative to other CHARMM commands:
  ----------------------------------------------------------------------------
  Before MRMD is called, it is important that the bonded parameter list
  is up-to-date, therefore it is recommended to execute the UPDAte
  command before calling MRMD. It requires that the PSF is already
  built.

  Before executing the DYNA command it is important that crossing
  geometry (UCRG) file is opened.

  All commands (MINI, VIBR, DYNA) which are supposed to use the force
  field parametrization provided by MRMD should be preceded by
  calling MRMD.

  A typical input sequence for an MRMD simulation is as follows:
  "
  ...read or generate PSF...

  UPDATE
  OPEN UNIT 9 READ FORMATTED NAME mrmd.par
  OPEN UNIT 14 WRITe FORMatted NAME crossings.pdb
  MRMD UPAR 9 UCRG 14 PRCA 10 PRDY 1000

  MINI ...

  OPEN UNIT 12 WRITe FORMatted NAME new.res
  OPEN UNIT 13 WRITe UNFORMatted NAME  mrmd_traj.dcd

  DYNA ....
  "
  ----------------------------------------------
  Compatibility with other major CHARMM commands
  ----------------------------------------------
  - PSF modifying routines -
  Calling MRMD assumes that the PSF will not be modified later. If it is
  modified then MRMD has to be called again either with the same or a
  new parameter file (UPAR).

  - VIBRan -
  MRMD does not provide (analytical) second derivatives, therefore all
  routines which requires second derivatives of the energy can work only
  if they have an option for numerical second derivatives (i.e. VIBR
  with FINITE).

  - IMAGe -
  The nonbonded interaction between atoms with modified non-bonded
  parameters and image atoms are not updated. If the atoms with modified
  nonbonded parameters remain far from the edges of the center box then
  this will not cause a problem at all. It is planned to provide this
  functionality in the future.

  - NBONded -
  compatible with NBONded for the following settings:
      ELEC   : electrostatic interaction
      CDIElec: constant dielectric
      ATOM   : atom-based electrostatic interaction
      SHIFt  : shifting cut-off method for electrostatic interaction 
      VDW    : van der Waals interaction by Lennard-Jones potential
      VSWItch: switching cut-off method for van der Waals interaction 
      VATOm  : atom-based van der Waals interaction
      NBXMOD : +1,+2,+3,+4,+5 are all supported

  not compatible with NBONded for the following settings:
      GROUp   : group-based electrostatic interaction
      SWITch  : switching cut-off method for electrostatic interaction
      FSWItch : force-switching cut-off method for electrostatic interaction
      VGROup  : group-based van der Waals interaction
      VSHIft  : shifting cut-off method for van der Waals interaction
      FVSWitch: force-switching cut-off method for van der Waals interaction      

  - SKIP -
  ARMD is compatible with SKIPe commands
  SKIPe BOND : skips bond energy correction defined in BOND HARM
               skips bond energy correction defined in BOND MORS
  SKIPe ANGL : skips energy correction defined in angle part of ANGL HARM
  SKIPe UREY : skips energy correction defined in Urey-Bradley part of ANGL HARM
  SKIPe DIHE : skips energy correction defined in DIHE FOUR
  SKIPe IMPR : skips energy correction defined in IMPR HARM
  SKIPe ELEC : skips energy correction defined in point charge part of NBON ATOM
  SKIPe VDW  : skips energy correction defined in Lennard-Jones part of NBON ATOM
               skips energy correction defined in NBON GVDW

  
  File: MRMD, Node: Parameter file, Previous: Description, Up: Top, Next: States

  The user has to provide an additional external parameter file beyond the usual 
  CHARMM parameter file to describe one or more potential energy surfaces, which 
  together define a global surface. This (re)parametrization is based on atom 
  indices instead of atom types, therefore allows specific reparametrization of 
  any part of the system.

  -------------------------------------------------------
  Determining effective reactive potential energy surface
  -------------------------------------------------------
  First the energy correction to the CHARMM energy for each surface is calculated 
  based on the (re)parametrization given in this file. Then, the energy shift is 
  added to each surface. Then a smooth connection of the low-lying energy surfaces
  is calculated using the energy difference based switching functions. Finally, 
  the Gaussian*polynomial functions of pairwise energy differences are added 
  optionally to adjust the shape of the crossing regions between surfaces. 
  The MRMD module in CHARMM corrects the total energy, and does not modify the
  individual energy terms (BOND, VDW, etc). MRMD energy correction will be added 
  only to the total energy.

  The parameter file should always have all the following blocks even
  if the parameters are not changed/used:
  SURF:      defines the number of surfaces and their level shifts
  SWCH:      defines switching function and switching related parameters
  SHAP GAPO: defines shaping functions to crossing regions between surfaces 
  NBON ATOM: (re)parametrization of charges and Lennard-Jones potentials of atoms
  NBON GVDW: parametrization of pair-wise generalized Lennard-Jones potentials 
  BOND HARM: (re)parametrization of harmonic bond potentials 
  BOND MORS: parametrization of Morse potentials
  ANGL HARM: (re)parametrization of harmonic angle potentials 
  DIHE FOUR: (re)parametrization of Fourier dihedral angle potentials  
  IMPR HARM: (re)parametrization of harmonic improper dihedral angle potentials 

  
  File: MRMD, Node: States, Previous: Parameter file, Up: Top, Next: Nonbonded

  Description of the individual blocks:
  ================================================================================
  1. Block SURF
  ================================================================================
  ROLE: 
  Defines the number of surfaces (nsurf) and the energy shifts (VshiftI). At each 
  call to MRMD, the energy correction of all surfaces are determined, and based 
  on the energy difference based switching formula, corresponding weights are 
  determined and the effective correction is calculated.

  --------------------------Structure of block SURF-------------------------------
  SURF nsurf
  1 Vshift1
  2 Vshift2
  ...
  I VshiftI
  ...
  nsurf Vshiftnsurf
  --------------------------------------------------------------------------------
  SHIFTING THE SURFACES (kcal/mol):
  dVI(r)=dV0I(r)   +VshiftI
  VI(r) =VCHARMM(r)+dVI(r)

  dVI0: 
  Energy correction of surface I to the CHARMM energy based on the
  (re)parametrization of force field terms given in the MRMD parameter.

  dVI: 
  Total energy correction of surface I to CHARMM energy. Their weighted sum plus 
  the shaping functions determine the effective correction to the CHARMM energy.

  VCHARMM:
  CHARMM energy calculated without MRMD module

  VI:
  Energy of MRMD surface I. Their weighted sum plus the shaping functions 
  determine the effective global surface. 
  --------------------------------------------------------------------------------
  nsurf ({integer,>0})
  Number of surfaces. At least one is required. If only one surface is given then 
  no surfaces crossing possible.

  VshiftI ({real,kcal/mol}): 
  Energy shift for the I-th surface. It is defined in separate lines for each 
  surface in increasing order of the surface index.  Energy shifting is used to 
  adjust the level of each surface to a global reference scale (eg. ab initio 
  energies), which is necessary to reproduce reaction energetics and predict 
  dividing surfaces between them.

  ================================================================================
  2. Block SWCH
  ================================================================================
  ROLE: 
  Defines the switching function and switching related parameters.
  MS-ARMD switching concept:
  1. Normalized energy difference from the lowest-energy surface is determined.
  2. Raw weights(wi0) are defined using simple mathematical switching functions.
  3. Raw weights are renormalized to give mixing weights (wi). 
  4. The weighted linear combination of surfaces plus the shaping functions
     with the weights give the MS-ARMD surface.

  --------------------------Structure of block SURF-------------------------------
  SWCH
  switching_function deltaV
  --------------------------------------------------------------------------------
  switching_function (string): mathematical switching function 
  deltaV({real,>0,kcal/mol}) :  switching parameter
  Note, input for SWCH has to be provided even if there is only one surface!

  MS-ARMD switching concept:
  1. Normalized energy difference from the lowest-energy surface is determined:
  Vi                                     ! energy of the surface i
  Vmin(r)    =  min(Vi(r),i=1...n)       ! energy of the lowest-energy surface
  deltai(r)  = (Vi(r)-Vmin(r))/deltaV    ! normalized energy difference of surface 
                                         ! i from the lowest-energy surface

  2. Raw weights wi0 are defined using switching functions (f).
  Currently available switching functions:
  A, JOHNSON7: 7th order switching function by B.R.Johnson 
  B, EXPDECAY: exponential decay switching function (recommended)
  w0i(r)=f(deltai(r))

  3. Raw weights are renormalized and the linear combination of surface
  energies are formed using the renormalized weights.
  wi(r)     = w0i(r)/(w01(r)+w02(r)+...) ! normalized weights
  Veff0(r)  = sum(wi(r)*Vi(r),i=1...n)   ! effective, smooth surface
  Veff0(r)  = sum(wi(r)*Vi(r),i=1...n)=VCHARMM+sum(wi(r)*dVi(r),i=1...n)
  Veff(r)   = Veff0(r)+shaping functions
  --------------------------------------------------------------------------------
  7th order switching function by B.R. Johnson 
  B.R. Johnson: J. Chem. Phys. 83, 1204 (1985)
  --------------------------------------------------------------------------------
  JOHNSON7 delta
  --------------------------------------------------------------------------------
  w0i(r) = f(deltai(r))=fJ(1-deltai(r))
  where:
        { 0                              if x<0  
  fJ(x)={ -x**4*(((20*x-70)*x+84)*x-35)  if 0<x<1
        { 1                              if x>1  
  deltaV ({real,>0,kcal/mol}): 
  if the energy of a surface is less than the energy of the lowest-lying surface 
  +delta, then it will have a non-zero weight, and start contributing to the 
  effective surface. When two or more lowest lying surfaces are crossing, their
  weights change between 0-1.

  NOTE: while the Johnson's switching function is 3 times continuously
  differentiable, the final MS-ARMD surface is not differentiable if the 
  three lowest surfaces are within the deltaV energy of each other.

  --------------------------------------------------------------------------------
  EXPDECAY switching function
  (See JCTC 2013 MS-ARMD paper)
  --------------------------------------------------------------------------------
  EXPDECAY delta
  --------------------------------------------------------------------------------
  w0i(r) = exp(-deltai(r))
  --------------------------------------------------------------------------------
  deltaV ({real,>0,kcal/mol}): 
  The raw weight of surfaces drops exponentially with increasing energy.
  This exponential decay has a characteristic energy of deltaV, that is
  an increase in energy by deltaV causes a drop by a factor of
  e~=2.7 in the weight. EXPDECAY switching function is recommended over JOHNSON7 
  switching as it always provides a smooth (analytical, infinite times 
  differentiable) effective surface.

  ================================================================================
  3. Block SHAP GAPO
  ================================================================================
  ROLE: 
  Defines the parameters of Gaussian*Polynomial (GAPO) shaping functions to adjust 
  the crossing region for pairs of surfaces.

  --------------------------Structure of block SURF-------------------------------
  SHAP GAPO ngapo
  surf11 surf12 n1 x10 x11 a10 a11 a12 ... a1n1
  surf21 surf22 n2 x20 x21 a20 a21 a22 ... a2n2
  ...
  surfI1 surfI2 nI xI0 xI1 aI0 aI1 aI2 ... aI2nI
  ...
  ngapo...
  --------------------------------------------------------------------------------
  Energy correction:
  Veff(r)=Veff0(r)+sum_I[ {wsurfI1+wsurfI2}*VGAPOI(VsurfI2(r)-VsurfI1(r)) ]
  --------------------------------------------------------------------------------
  GAUSSIAN*POLYNOMIAL FUNCTION (kcal/mol):
  xI=VsurfI2(r)-VsurfI1(r)               
  VGAPOI(xI) = exp(-(xI-xI0)^2/(2 xI1^2))*
              *[aI0+aI1*(xI-xI0)+aI2*(xI-xI0)^2+...+aInI*(xI-xI0)^nI]
  where r(angstroem) is the coordinates of all atom.
  --------------------------------------------------------------------------------
  ngapo ({integer,>=0}): 
  Total number of GAPO functions used for adjusting the crossing regions.

  surfI1, surfI2 ({integer,nsurf>=,>=1}):
  The indices of the two surfaces of the I-th GAPO function.

  nI ({integer,>=0})
  Polynomial order of the I-th GAPO function.

  xI0 ({real,1/(kcal/mol)})
  Shift of the I-th GAPO function.

  xI1 ({real,>0,(kcal/mol)})
  Standard deviation parameter of the Gaussian of the I-th GAPO function.

  aI0,aI1,...,aInI ({real,(kcal/mol)^(1-j) where j=0,1,2,3})
  0-th, 1-st, ..nI-th order polynomial coefficients of the I-th GAPO function.

  
  File: MRMD, Node: Nonbonded, Previous: States, Up: Top, Next: Bonded

  ================================================================================
  4. Block NBON ATOM
  ================================================================================
  ROLE: 
  Reparametrizes LJ interaction and point-charge electrostatic interaction on
  each surface.

  --------------------------Structure of block ATOM-------------------------------
  NBON ATOM natom
  atom1   q11 eps11 rmh11 xeps11 xrmh11    q12 eps12 rmh12 xeps12 xrmh12   ...
  atom2   q21 eps21 rmh21 xeps21 xrmh21    q22 eps22 rmh22 xeps22 xrmh22   ...
  ...
  atomI   qI1 epsI1 rmhI1 xepsI1 xrmhI1    qI2 epsI2 rmhI2 xepsI2 xrmhI2   ...
  ...
  atomnatom ...
  --------------------------------------------------------------------------------
  ELECTROSTATIC INTERACTION: COULOMB POTENTIAL WITH SHIFTING (kcal/mol):
  V(r)=fshift(r)/4pi/EPS0/EPS*q*q/r
  where r is the distance of two atoms.

  S1(r)=fshift(r)={  0                    if r >roff }
                  { (1-(r/roff)^2)^2      if r<=roff } 

  V(r)=V(r)*E14FAC for atoms in 1-4 position if NBXMOD=1,2,3,5
  where roff=CTOFNB
  NBXMOD,CTOFNB,EPS,E14FAC are described in nbond.doc
  --------------------------------------------------------------------------------
  VAN DER WAALS INTERACTION: LENNARD-JONES POTENTIAL WITH SWITCHING (kcal/mol): 

  eps=sqrt(eps1*eps2)

  rmin=rmh1+rmh2

  V(r;eps,rmin)=fswitch(r)*eps*[(rmin/r)^12-2*(rmin/r)^6]
  where r(angstroem) is the distance of atomI1 and atomI2.

             {1                                                     if r<ron     }
  fswitch(r)={(roff^2-r^2)^2*(roff^2+2r^2-3ron^2)/(roff^2-ron^2)^3  if ron<r<roff}
             {0                                                     if roff<r    }

             where:roff=CTOFNB,ron=CTONNB
  CTONNB,CTOFNB are described in nbond.doc
  --------------------------------------------------------------------------------
  natom ({integer,>=0})
  Number of atoms with reparametrized nonbonded interaction. 

  atomI ({integer,>0}): 
  PSF index of atom I to be reparametrized, no duplicate definition is
  allowed for the same atom.

  qIJ ({real, elementary charge unit}): 
  Charge of atom I on surface J

  epsIJ ({real,>=0,kcal/mol}): 
  Well depth of Lennard-Jones potential for atom I on surface J.

  rmhIJ ({real,>0,angstroem}): 
  Half of the separation corresponding to the minimum of the Lennard-Jones 
  potential well ("Rmin/2") for atom I on surface J (usual CHARMM convention).

  xepsIJ ({x} or {X} or {real,>=0,kcal/mol}): 
  epsIJ value for "special van der Waals interaction" (CHARMM jargon), which acts 
  between atoms in 1-4 position and described by Lennard-Jones potential, which is
  active in the case of NBXMOD=+5. Opposite to CHARMM convention, this is always 
  given as non-negative value. If letter x or X is given, then the value copied 
  from epsIJ.

  xrmhIJ ({x} or {X} or {real,>0,angstroem}): 
  rmhIJ value for "special van der Waals interaction" (CHARMM jargon), which acts 
  between atoms in 1-4 position and described by Lennard-Jones potential, which is
  active in the case of NBXMOD=+5. If letter x or X is given, then the value 
  copied from rmhIJ.

  ================================================================================
  5. Block NBON GVDW
  ================================================================================
  ROLE: 
  Adds general-exponent Lennard-Jones potential (equivalent to the Mie-potential)
  and removes 6-12 Lennard-Jones potential for pair of atoms on each surface.
  General-exponent Lennard-Jones potential can improve the description of the van 
  der Waals interaction between atoms in the reaction region and it also better 
  captures the energetics for the reactant and product complexes and around the 
  transition state.
  --------------------------------------------------------------------------------
  VAN DER WAALS INTERACTION: GENERAL-EXPONENT LENNARD-JONES POTENTIAL (kcal/mol)

  V(r;eps,rmin,rep,att)=att/(rep-att)*eps*[(rmin/r)^rep-rep/att*(rmin/r)^att]

  --------------------------Structure of block GVDW-------------------------------
  NBON GVDW ngvdw
  atom11 atom12  eps11 rmin11 att11 rep11   eps12 rmin12 att12 rep12    ....
  atom21 atom22  eps21 rmin21 att21 rep21   eps22 rmin22 att22 rep22    ....
  ...                                                                
  atomI1 atomI2  epsI1 rminI1 attI1 repI1   epsI2 rminI2 attI2 repI2    ....
  ...
  atomngvdw1 atomngvdw2 ...
  --------------------------------------------------------------------------------
  ngvdw ({integer,>=0})
  Number of atom pairs with reparametrized van der Waals interaction.

  atomI1,atomI2 ({integer,>0}): 
  PSF indices of atom pair I. 
  No duplicate definition is allowed for the same pair of atoms (ij=ji).

  epsIJ ({x} or {X} or {real,>=0,kcal/mol}): 
  Well depth of general-exponent Lennard-Jones potential for atom pair I on 
  surface J. If letter x or X is given for epsIJ, then following values of rminIJ,
  attIJ and repIJ are not used, and the 6-12 Lennard-Jones potential will be used 
  for the corresponding atom pair on surface J.

  rminIJ ({real,>0,angstroem}): 
  Distance of atoms at the minimum of general-exponent Lennard-Jones potential I
  on surface J. If letter x or X is given for epsIJ, then following values of 
  rminIJ, attIJ and repIJ are not used, and the 6-12 Lennard-Jones potential from 
  the CHARMM parameter file will be used for the corresponding atom pair on 
  surface J.

  attIJ ({real,>0,repIJ>}): 
  Attractive exponent of general-exponent Lennard-Jones potential I on surface J. 
  If letter x or X is given for epsIJ, then following values of rminIJ, attIJ and 
  repIJ are not used, and the 6-12 Lennard-Jones potential will be used for the 
  corresponding atom pair on surface J.

  repIJ ({real,>0,>attIJ}): 
  Repulsive exponent of general-exponent Lennard-Jones potential I on surface J. 
  If letter x or X is given for epsIJ, then values of rminIJ, attIJ and repIJ are 
  not used, and the 6-12 Lennard-Jones potential will be used for the 
  corresponding atom pair on surface J.

  
  File: MRMD, Node: Bonded, Previous: Nonbonded, Up: Top, Next: Output 

  ================================================================================
  6. Block BOND HARM
  ================================================================================
  ROLE: 
  Reparametrizes, adds or removes harmonic bond potentials on surfaces.

  --------------------------Structure of block HARM-------------------------------
  BOND HARM nharm
  atom11 atom12    fch11 req11    fch12 req12    ...
  atom21 atom22    fch21 req21    fch22 req22    ...
  ...
  atomI1 atomI2    fchI1 reqI1    fchI2 reqI2    ...
  ...
  atomnharm1 atomnharm2 ...
  --------------------------------------------------------------------------------
  HARMONIC BOND POTENTIAL (kcal/mol):
  V(r;fch,req)=fch*(r-req)^2
  where r(angstroem) is the length of bond atomI1-atomI2.
  --------------------------------------------------------------------------------

  nharm({integer,>=0}): 
  Number of those atom pairs between which the bond is (re)defined.

  atomI1, atomI2 ({integer,>0}): 
  PSF index of the two atoms forming the I-th bond.  No duplicate definition is 
  allowed for the same pair of atoms (ij=ji). The definition for the same pair 
  of atoms is not allowed in the MORS block. 

  fchIJ ({x} or {X} or {real,>=0, kcal/mol/angstroem^2}): 
  Half force constant of bond potential I on surface J (usual CHARMM convention)
  If letter x or X is given then the bond is considered as broken in this surface 
  (even if it is in the PSF file) If letter x or X is given for fchIJ, then the 
  following value of reqIJ is not used.

  reqIJ ({>0,real,angstroem}): 
  Equilibrium bond length of bond I on surface J. If letter x or X is given for 
  fchIJ, then the following value of reqIJ is not used.

  ================================================================================
  7. Block BOND MORS
  ================================================================================
  ROLE: 
  Adds Morse bond potential, removes harmonic bond potential on surfaces.
  If the same bond exists in PSF, then it is automatically replaced with this.

  --------------------------Structure of block MORS-------------------------------
  BOND MORS nmors
  atom11 atom12   de11 req11 beta11   de12 req12 beta12    ...
  atom21 atom22   de21 req21 beta21   de22 req22 beta22    ...
  ...
  atomI1 atomI2   deI1 reqI1 betaI1   deI2 reqI2 betaI2    ...
  ...
  atomnmors1 atomnmors2 ...
  --------------------------------------------------------------------------------
  MORSE-POTENTIAL (kcal/mol):
  V(r)=de*[1-exp(-beta*(r-req))]^2
  where r(angstroem) is the length of bond atomI1-atomI2.
  --------------------------------------------------------------------------------
  nmors({integer,>=0}): 
  Number of redefined Morse bonds.

  atomI1, atomI2 ({integer,>0}): 
  PSF indices of atom pair I. 
  No duplicate definition is allowed for the same pair of atoms (ij=ji). 
  The definition for the same pair cannot be used in the HARM block. 

  deIJ ({x} or {X} or {real,>=0, kcal/mol}): 
  Well depth of Morse bond I on surface J.
  If letter x or X is given the bond is broken in this surface (even if it is a 
  PSF harmonic bond). If letter x or X is given for deIJ, then following values 
  of reqIJ and betaIJ are not used.

  reqIJ ({>0,real,angstroem}):
  Equilibrium bond length of Morse bond I on surface J. If letter x or X is given 
  for deIJ, then following values of reqIJ and betaIJ are not used.

  betaIJ ({>0,real, 1/angstroem}):
  Beta parameter of Morse bond I on surface J.
  If letter x or X is given for deIJ, then folowing values of reqIJ and betaIJ are
  not used.

  ================================================================================
  8. Block ANGL HARM
  ================================================================================
  ROLE: 
  Reparameterizes, adds or removes harmonic angle potentials and Urey-Bradley 
  potentials on surfaces.

  --------------------------Structure of block ANGL-------------------------------
  ANGL HARM nangl
  atom11 atom12 atom13 fch11 phieq11 ufch11 ureq11 fch12 phieq12 ufch12 ureq12 ... 
  atom21 atom22 atom23 fch21 phieq21 ufch21 ureq21 fch22 phieq22 ufch22 ureq22 ... 
  ...
  atomI1 atomI2 atomI3 fchI1 phieqI1 ufchI1 ureqI1 fchI2 phieqI2 ufchI2 ureqI2 ...
  ...
  atomnangl1 atomnangl2 ...
  --------------------------------------------------------------------------------
  HARMONIC ANGLE POTENTIAL (kcal/mol): 
  V(phi;fch.phieq)=fch*(phi-phieq)^2
  where phi(radian) is the angle formed by atomI1-atomI2-atomI3.
  --------------------------------------------------------------------------------
  UREY-BRADLEY POTENTIAL (kcal/mol):
  V(r;ufch,ureq)=ufch*(r-ureq)^2
  where r(angstroem) is the distance between atomI1 and atomI3.
  --------------------------------------------------------------------------------
  nangl({integer,>=0}): 
  Number of (re)defined harmonic angle and Urey-Bradley potentials.

  atomI1, atomI2, atomI3  ({integer,>0}): 
  PSF indices of the three atoms forming the I-th angle (atomI1-atomI2-atomI3). 
  No duplicate definition is allowed for the same three atoms in the give order
  (ijk=kji).

  fchIJ ({x} or {X} or {real,>=0,kcal/mol/radian^2}): 
  Half force constant of angle potential I on surface J (usual CHARMM convention).
  If letter x or X is given the angle potential is not present on this surface due
  to missing bond. If letter x or X is given for fchIJ, then following values of 
  phieqIJ, ufchIJ, ureqIJ are not used.

  phieqIJ ({>0,real,degree}): 
  After reading it in, it is immediately converted to radian. Equilibrium angle of
  angle potential I on surface J. If letter x or X is given for fchIJ, then 
  following values of phieqIJ, ufchIJ, ureqIJ are not used.

  ufchIJ ({x} or {X} or {real,>=0,kcal/mol/angstroem^2}): 
  Half force constant of Urey-Bradley potential I on surface J 
  (usual CHARMM convention). 
  If letter x or X is given the Urey-Bradley potential is not present on this 
  surface and the following values of phieqIJ, ufchIJ, ureqIJ are not used. 

  ureqIJ ({real,>0,degree}): 
  Equilibrium distance of Urey-Bradley potential I on surface J. If letter x or X 
  is given for fchIJ, then following values of phieqIJ, ufchIJ, ureqIJ are not 
  used.

  ================================================================================
  9. Block DIHE FOUR
  ================================================================================
  ROLE: 
  Reparameterizes, adds or removes proper dihedral potentials on surfaces.

  --------------------------Structure of block DIHE-------------------------------
  DIHE FOUR ndihe
  atom11 atom12 atom13 atom14 per1    amp11 phi011   amp12 phi012    ...
  atom21 atom22 atom23 atom24 per2    amp21 phi021   amp22 phi022    ...
  ...
  atomI1 atomI2 atomI3 atomI4 perI    ampI1 phi0I1   ampI2 phi0I2    ...
  ...
  atomndihe1 atomndihe2 atomndihe3 atomndihe4 ...
  --------------------------------------------------------------------------------
  PROPER DIHEDRAL POTENTIAL (kcal/mol): 
  per=1,2,..,6 => V(phi)=amp*[1+cos(per*phi-phi0)]    
  per=0        => V(phi)=amp*(phi-phi0)^2
  where phi(radian) is the dihedral angle formed by atomI1-atomI2-atomI3-atomI4.
  Expected bond connectivity of atoms: atomI1-atomI2-atomI3-atomI4
  --------------------------------------------------------------------------------
  ndihe({integer,>=0}): 
  Number of (re)defined proper dihedral potentials.

  atomI1, atomI2, atomI3, atomI4  ({integer,>0}): 
  PSF indices of the four atoms forming the I-th dihedral angle. No duplicate 
  definition is allowed for the same 4 atoms in the given order (ijkl=lkji) 
  with the same periodicity.

  perI ({integer:1,2,3,4,5,6}): 
  Periodicity of cosine dihedral potential I. For the same group of atoms multiple
  declarations with different periodicity are allowed (Fourier series:1,2,3,4,5,6)
  No duplicate definition is allowed for the same 4 atoms in the given order
  (ijkl=lkji) with the same periodicity.
                              
  ampIJ ({x} or {X} or {real,kcal/mol},{real,>0,kcal/mol}): 
  Amplitude of dihedral potential I on surface J. If letter x or X is given for 
  ampIJ, then the dihedral potential does not exist on the corresponding surface 
  and the following value of phi0IJ is not used.

  phi0IJ ({x} or {X} or {real,degree}): 
  Zero phase for dihedral potential. After reading it in, it is immediately 
  converted to radian within the code. 
  The location of extrema are:
  if ampI>0:  
   cos(n*phi_max-phi0) maximal <=> n*phi_max-phi0=2*k*pi       k is integer
   =>maxima: phi_max=(2*k*pi+phi0)/n                           k is integer
   cos(n*phi_min-phi0) minimal <=> n*phi_min-phi0=(2*k+1)*pi   k is integer
   =>mixima: phi_min=((2*k+1)*pi+phi0)/n                       k is integer   
   
  if ampI<0: 
   -cos(n*phi_min-phi0) minimal <=> n*phi_min-phi0=2*k*pi      k is integer
   =>minima: phi_min=(2*k*pi+phi0)/n                           k is integer
   -cos(n*phi_max-phi0) maximal <=> n*phi_max-phi0=(2*k+1)*pi  k is integer
   =>maxima: phi_max=((2*k+1)*pi+phi0)/n                       k is integer   

  ================================================================================
  10. Block IMPR HARM
  ================================================================================
  ROLE: 
  Reparameterizes, adds or removes improper dihedral potentials on surfaces.

  --------------------------Structure of block IMPR-------------------------------
  IMPR HARM nimpr
  atom11 atom12 atom13 atom14 per1    fch11 phi011   fch12 phi012    ...
  atom21 atom22 atom23 atom24 per2    fch21 phi021   fch22 phi022    ...
  ...
  atomI1 atomI2 atomI3 atomI4 perI    fchI1 phi0I1   fchI2 phi0I2    ...
  ...
  atomnimpr1 atomnimpr2 atomnimpr3 atomnimpr4 ...
  --------------------------------------------------------------------------------
  IMPROPER DIHEDRAL POTENTIAL (kcal/mol): 
  perI=0        V(phi)=fch*(phi-phi0)^2
  where phi(radian) is the dihedral angle formed by atomI1-atomI2-atomI3-atomI4.
  Expected bond connectivity of atoms: 
  atomI1 - atomI4
    |   \
  atomI2 atomI3

  IMPORTANT:
  CHARMM topology should use the same topology definition for improper dihedrals
  in the reaction center otherwise dihedrals will not be successfully removed
  on MRMD surface if a bond gets broken. Improper dihedrals to be removed or
  reparameterized have to have the same or reverse atom order.
  ------------------------------------------------------------------------------
  nimpr({integer,>=0}): 
  Number of (re)defined improper dihedral potentials.

  atomI1, atomI2, atomI3, atomI4  ({integer,>0}): 
  PSF indices of the four atoms forming the I-th improper dihedral angle. 
  No duplicate definition is allowed for the same group of atoms in the 
  given order (ijkl=lkji)

  perI ({integer:0}): 
  Periodicity of dihedral potential I. Only perI=0 is allowed and it implies 
  quadratic angle potential.

  fchIJ ({x} or {X} or {real,kcal/mol/radian^2}): 
  Half force constant of quadratic dihedral potential I on surface J (usual CHARMM
  convention). If letter x or X is given for fchIJ, then the improper potential
  does not exist on the corresponding surface and the following value of phi0IJ is
  not used.

  phi0IJ ({x} or {X} or {real,degree}):
  Equilibrium angle of quadratic dihedral potential I on surface J. After reading 
  it in, it is immediately converted to radian within the code.  
                  
  
  File: MRMD, Node: Output, Previous: Bonded, Up: Top, Next: Example

  Output to standard output:

  ----------------------------------------------
  1. At the call of MRMD command in CHARMM input
  ----------------------------------------------
  Subroutine MRMD_INIT will be executed, which will print:
  Printed lines start with:
  'MRMD_INIT>...'

  if PRNLEV>=5 the following will be printed:
  - the arguments with which the MRMD command was called
  - progress of reading various blocks of MRMD parameter file 
  - after reading the parameter file, the number of records in each block

  if PRNLEV>=6 the following will be printed:
  - all the parameters that are read from the MRMD parameter file.
  - changes in 1-2, 1-3, 1-4, 1-(more than 4) neighbourship lists for each surface

  -----------------------------------------------------
  2. At energy call if certain conditions are fulfilled
  -----------------------------------------------------
  Subroutine EMRMD will be executed, which calls several other subroutines.
  Printed lines start with:
  'EMRMD>...'
  'MRMD_ENERG>...'
  'MRMD_SURF>...'

  -------------------------------
  2A.Printing of surface crossing:
  -------------------------------
  - PRNLEV>=4: time of crossing and the corresponding two surfaces

  ---------------------------------------------
  2B.Printing energetics at various detail level
  ---------------------------------------------
  - During dynamics after every integration step defined after PRDY argument.
  - When dynamics is not active after every PRCA calls.
  - Initial call of EMRMD routine. 
  - When surface crossing is detected.

  PRNLEV>=4 
  MRMD energy correction.

  PRNLEV>=5
  MRMD energy correction for each surface.

  PRNLEV>=6
  Summed energy correction of each type 
  (NBON ELEC,NBON VDW,NBON GVDW,BOND HARM,BOND MORS,
  ANGL HARM, ANGL UREY,DIHE FOUR,IMPR HARM) for each surface.

  PRNLEV>=7
  Each individual MRMD energy term, except for nonbonded interactions.
  Each GAPO energies.

  PRNLEV>=8
  Each individual nonbonded energy term is also printed.

  -------------------------------------------
  2C. Printing forces at various detail level
  -------------------------------------------
  PRNLEV>=9
  All nonzero MRMD force corrections for each atom are printed.

  PRNLEV>=10
  All forces of all potential energy surfaces are printed.

  
  File: MRMD, Node: Example, Previous: Output, Up: Top, Next: Troubleshooting

  The mrmd_h2so4.inp test case with the mrmd_h2so4.par parameter file are 
  distributed with the package. This example simulation describes the water 
  elimination from a vibrationally highly excited sulfuric acid molecule.  
  Two surfaces are defined: one for sulfuric acid molecule and one for water 
  and sulfur-trioxid molecules. The initial conditions are prepared so that 
  water elimination takes place within a few timesteps. All the force field 
  parameters of the H2SO4 molecule are redefined in the parameter file 
  (mrmd_h2so4.par) in order to demonstrate most of the features of the 
  MRMD module.

  
  File: MRMD, Node: Troubleshooting, Previous: Errors, Up: Top, Next: Code

  Before running production calculations, it is highly recommended to check 
  whether all the expected individual energy terms are cancelled out or/and 
  added for each surfaces. Increasing the PRNLEV parameter up to 8, a detailed 
  listings is done, showing also the values of all used parameters and geometric 
  variables (ie. bond length) used for the calculation of the removed and added 
  terms. If total energy seems to be wrong, this listing can help in identifing 
  the energy terms which are responsible for it.

  
  File: MRMD, Node: Code, Previous: Troubleshooting, Up: Top

  ----------------------------------
  Interfaces, subroutines, functions
  ----------------------------------
  memory (de)(re)allocation interface
  MRMD_ALLOC using subroutines:
  MRMD_ALLOC_REAL1,MRMD_ALLOC_REAL2,MRMD_ALLOC_INTG1,MRMD_ALLOC_INTG2,
  MRMD_ALLOC_INTG3,MRMD_ALLOC_LOGI1,MRMD_ALLOC_LOGI2,MRMD_ALLOC_CH16_1

  reads and processes parameter file and sets up all global variables
  MRMD_INIT

  energies, forces for effective and individual surfaces
  EMRMD,MRMD_ENERG,MRMD_SURF

  switching related
  MRMD_SWITCH_FUNCT,MRMD_SWITCH_WEIGHTS

  energy terms:
  MRMD_GAPO,MRMD_ELEC,MRMD_VDW,MRMD_GVDW,MRMD_HARM,MRMD_MORS,MRMD_ANGL,MRMD_DIHE

  write out pdb file
  MRMD_PDB

  converting strings to uppercase:
  function MRMD_UPPERCASE

  complete deallocation
  MRMD_DEALLOC

  ----------------------------------
  Global variables in other routines
  ----------------------------------
  variable for the preprocessor:
  RMD: logical variable to include MRMD module or not into the code to be compiled

  Fortran variables, constants:
  logical MRMD_ACTIVE  = whether MRMD module is active (reactive surface loaded) 
  integer MRMD         = 103
  string  CETERM(MRMD) = 'ERMD' 
                       => use ?ERMD to request ERMD energy correction
  logical QETERM(MRMD) = whether MRMD modul is compiled
                       => use ?MRMD to request whether MRMD module is compiled
                       ?MRMD .eq. 1 => yes
                       ?MRMD .ne. 1 => no
  energy  ETERM(MRMD)  = energy correction to PSF

  --------------------------------------
  Other modules refering to MRMD module:
  --------------------------------------
  charmm/miscom.src : calls MRMD_INIT
  energy/energy.src : calls EMRMD
  energy/energym.src: defines MRMD=103
  energy/eutil.src  : defines CETERM(MRMD)='ERMD'                        
  misc/genetic.src  : calls EMRMD
  pert/epert.src    : calls EMRMD 
  pert/icpert.src   : calls EMRMD
