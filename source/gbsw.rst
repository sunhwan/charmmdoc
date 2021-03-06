.. py:module::gbsw

=============================================================================================================
Generalized Born with a simple SWitching (GBSW) (Electrostatic + Nonpolar) Solvation Energy and Forces Module
=============================================================================================================


.. note::

   **Questions and comments regarding GBSW should be directed to**
   
   - Wonpil Im (wonpil@ku.edu)
   - Charles L. Brooks, III (brooks@scripps.edu)

   **References for GBSW:**
   
   1. W. Im, M.S. Lee, and C.L. Brooks III
      "Generalized Born Model with a Simple Smoothing Function."
      J. Comput. Chem. 24:1691-1702 (2003). 

   2. W. Im, M. Feig, and C.L. Brooks III
      "An Implicit Membrane Generalized Born Theory for the Study of 
      Structure, Stability, and Interactions of Membrane Proteins."
      Biophys. J. 85:2900-2918 (2003).

   3. W. Im, J. Chen, and C.L. Brooks III
      "Application of a rotationally invariant procedure to
      a Generalized Born Model"
      in preparation (2005).


.. _gbsw_description:

The GBSW module provides the (electrostatic + nonpolar) solvation
energy and forces.  A Generalized Born method is used for the
electrostatic part and the solvent-exposed surface ares for the
nonpolar part with a phenomenological surface tension coefficient.
Based on volume integration schemes used in the GBMV module [M.S. Lee,
F.R. Salabury, Jr., and C.L. Brooks III, J. Chem. Phys., 116, 10606
(2002)], we have recast the calculation of the self-electrostatic
solvation energy to utilize  a simple smoothing function at the
dielectric boundary. The GBSW model is formulated in  this manner to
provide consistency with the Poisson-Boltzmann (PB) theory previously
developed to yield numerically-stable electrostatic solvation forces
based on finite-difference methods [W. Im, D. Beglov, and B. Roux,
Comp. Phys. Comm., 111, 59 (1998)].  However, it is also possible to
mimic the PB results with the molecular surface by reparametrizing two
adjustable parameters, a_0 to modulate the Coulomb field term and a_1
to include a correction term beyond Coulomb field.

The GBSW module takes the influence of biological membranes into
account. Consistent with continuum Poisson-Boltzmann (PB)
electrostatics, the membrane is approximated as an
solvent-inaccessible infinite planar low-dielectric slab. The membrane
GB model closely reproduces the PB electrostatic solvation energy
profile across the membrane.
    
The GBSW module works with the IMAGE facility. The GBSW calculations
are about 4 times slower than the corresponding vacuum
calculations. Using the simple smoothing function makes the present GB
model roughly 2-3 times faster than the GBMV module.

The backbone phi/psi cross-term (CMAP) and the atomic input radii were 
recently re-optimized specifically for GBSW implicit solvent, to balance the 
solvation and intramolecular interactions and to reproduce experimental 
conformational equilibria of a range of peptides.  See Examples for detailed
description of the optimal settings (Ref: Chen, Im and Brooks, JACS, 2006).


.. _gbsw_syntax:

::

   GBSW [SW real] [AA0 real] [AA1 real] [MOLSURF] [GBENER] [ROTINV] -
        [NANG integer] [NRAD integer] [RMAX real] [DGP real] [RBUFFER real] - 
        [EPSP real] [EPSW real] [CONC real] [TEMP real] [SGAMMA real] - 
        [TMEMB real] [MSW real] -
        [IGBFRQ integer]

   GBSW RESET  ! reset GBSW

======== ========= ==========================================================
SW       [0.3]     half of smoothing length in Ang.
                   (default value is changed to 0.2 when MOLSURF is issued.)
AA0      [aa0(sw)] coefficient for the Coulomb Field Approximation term
AA1      [aa1(sw)] coefficient for the correction term
                   (optimized default values for aa0(sw) and aa1(sw)
                   are given below) 
MOLSURF  [FALSE]   approximation to PB with molecular surface
GBENER   [FALSE]   calculate and print the solvation energy
                   (No cutoff is used for GB electrostatic solvation energy.)
ROTINV   [FALSE]   rotationally invariant numerical quadrature procedure
NANG     [38]      number of angular integration points
NRAD     [0]       number of radial integration points
                   (default value means the use of optimized 24 radial
                   integration points)
RMAX     [20.0]    maximum distance for radial integration in Ang.
DGP      [1.5]     grid spacing for lookup table in Ang.
RBUFFER  [0.0]     buffer length for lookup table in Ang.
EPSP     [1.0]     dielectric constant of both protein and reference state 
EPSW     [80.0]    solvent dielectric constant
CONC     [0.0]     salt concentration in M
TEMP     [300.0]   temperature in K (only necessary with CONC)
SGAMMA   [0.0]     nonpolar surface tension coefficients in kcal/(molxA^2)
TMEMB    [0.0]     thickness of low-dielectric membrane slab centered
                   at Z=0 (in Ang.)
MSW      [sw]      half of membrane switching length in Ang.
IGBFRQ   [1]       updating frequency of effective Born radii
======== ========= ==========================================================


.. _gbsw_function:

General discussion regarding the GBSW module
--------------------------------------------

1. Volume Integration

   The GBSW module uses the numerical quadrature method for the
   volume integration.  The integration points and weights for the radial
   component are generated by the Gaussian-Legendre quadrature those for
   the angular component by the Lebedev quadrature.  The default values
   for the integration (NANG, NRAD, and RMAX) should be appropriate
   for most cases.  However, one can specify NANG, NRAD, and RMAX
   independently. Note that NANG should be one of 26, 38, or 50.

   A grid-based lookup table is used to increase the efficiency of the
   integration.  Keywords DGP and RBUFFER are related with the lookup
   table. The current default value should be optimal for most case.
   However, one can check the efficiency and optimize those by performing
   short MD runs.


2. Choice of SW
   In prinicple, one can choose any SW. However, it should be noted that
   GBSW calculations take more time as SW increases. As default, SW=0.3
   is recommended for the smooth boundary and SW=0.2 for the molecular
   surface.  The optimized coefficients a_0 and a_1 are shown below for  
   each SW. Those coefficients were obtained by minimizing the error
   between GB and PB self-electrostatic solvation energies.

   * default A_0 and A_1 for the smoothed dielectric boundary
   
     === ======= ======= ======= ======= ======= ======= ======= ======= ======= =======
     SW     0.1     0.2     0.3     0.4     0.5     0.6     0.7     0.8     0.9     1.0
     A_0 -0.0811 -0.1481 -0.1801 -0.1680 -0.1542 -0.1731 -0.2279 -0.3064 -0.3943 -0.4820
     A_1  1.6000  1.7292  1.8174  1.8560  1.8864  1.9453  2.0359  2.1472  2.2645  2.3801
     === ======= ======= ======= ======= ======= ======= ======= ======= ======= =======

   * default A_0 and A_1 for the molecular surface

     ==== ======= ======= ========
     SW      0.1     0.2     0.3  
     A_0   1.2642  1.2045  1.1177
     A_1   0.0593  0.1866  0.3406
     ==== ======= ======= ========


3. Physical Parameters
   
   It should be noted that a_0 and a_1 were optimized with EPSP=1.0 and
   EPSW=80.0. Therefore, one should be careful when other values for
   EPSP and EPSW are used.  In other words, the electrostatic solvation
   contribution may not be optimal.  The influence of salt is taken into
   account based on the formalism of [J. Srinivasan, M.W. Trevathan,
   P. Beroza, and D.A. Case, Theor. Chem. Acc., 101, 426-434 (1999)].

   The nonpolar solvation contribution is considered only when non-zero
   SGAMMA is issued.  Note that the dimension is kcal/(molxA^2), and 0.01
   to 0.04 might be suitable for SGAMMA.

4. Low-dielectric slab for membrane
   
   The influence of membrane hydrophobic core as the low dielectric
   medium is approximately captured in the GBSW module (see reference 2
   for details).  Note that the membrane switching function is applied in
   the following region;

   * Z > 0 :  Tmemb/2.0 - MSW to  Tmemb/2.0 + MSW
   * Z < 0 : -Tmemb/2.0 + MSW to -Tmemb/2.0 - MSW


.. index:: gbsw; examples
.. _gbsw_examples:

Usage Examples
--------------

The examples below illustrate some of the uses of the GBSW module. 
(See also c30test/gbsw.inp)

There are two requirements for running GBSW;

1. "SWITCH" should be chosen in NBOND specifications.
2. WMAIN should be filled with a proper set of radii.  It is
   recommended to use the optimized PB radii
   (~charmm/test/data/radius.str) for the GBSW module.

.. note::

   1. A self-consistent GBSW force field is optimal for peptide and protein
      simulations (as of 2005). For this, WMAIN should be filled with a new set of
      radii ("radius_gbsw.str"). In addition, a special CMAP should be used for
      optimal treatment of peptide backbone ("par_all22_prot_gbsw.inp"). Both files
      locate in ~charmm/toppar/gbsw/.  (see Example 5 for illustration).

   2. For optimal performance in folding simulations, the following GBSW
      command-line options should be used with the self-consistent GBSW force field:
  
      ::
      
         GBSW sgamma 0.005 nang 50


Example 1
^^^^^^^^^

::

  !To perform a single-point energy calculation with infinite cutoffs: 

  prnlev 0
  stream radius.str
  prnlev 5
  scalar wmain statistics select .not. type H* end
  define check select (.not type H* ) .and. ( prop wmain .eq. 0.0 ) show end
  if ?nsel ne 0  stop       !some heavy atom have a zero radius

  GBSW sw 0.3 sgamma 0.03 dgp 1.5 GBenergy


Example 2
^^^^^^^^^

::

  !To perform a minimization or dynamics with cutoffs

  prnlev 0
  stream radius.str
  prnlev 5
  scalar wmain statistics select .not. type H* end
  define check select (.not type H* ) .and. ( prop wmain .eq. 0.0 ) show end
  if ?nsel ne 0  stop       !some heavy atom have a zero radius

  GBSW sw 0.3 sgamma 0.03 dgp 1.5 GBenergy

  NBOND atom switch cdie vdw vswitch -
        ctonnb 16 ctofnb 16 cutnb 20
  ENERGY

  (minimization or dynamics) 

Example 3
^^^^^^^^^

::

  !To perform a minimization or dynamics with images

  (image definition)

  NBOND atom switch cdie vdw vswitch -
        ctonnb 20 ctofnb 20 cutnb 24 cutim 24   ! should be before GBSW

  prnlev 0
  stream radius.str
  prnlev 5
  scalar wmain statistics select .not. type H* end
  define check select (.not type H* ) .and. ( prop wmain .eq. 0.0 ) show end
  if ?nsel ne 0  stop       !some heavy atom have a zero radius

  GBSW sw 0.3 sgamma 0.03 dgp 1.5 GBenergy

  ENERGY

  (minimization or dynamics) 


Example 4 
^^^^^^^^^

::

  !To perform a minimization or dynamics with membrane

  prnlev 0
  stream radius.str
  prnlev 5
  scalar wmain statistics select .not. type H* end
  define check select (.not type H* ) .and. ( prop wmain .eq. 0.0 ) show end
  if ?nsel ne 0  stop       !some heavy atom have a zero radius

  GBSW sw 0.3 sgamma 0.03 dgp 1.5 tmemb 35.0 msw 2.5

  NBOND atom switch cdie vdw vswitch -
        ctonnb 16 ctofnb 16 cutnb 20
  ENERGY

  (minimization or dynamics) 

Example 5
^^^^^^^^^

::
  
  !To setup for using the self-consistent GBSW force field

  ! read in the CMAP topology file (standard)
  open read card unit 10 name @toppar/top_all22_prot_cmap.inp
  read rtf  card unit 10
  close unit 10

  !read in the parameter file that contains GBSW specific CMAP
  open read card unit 10 name @topar/par_all22_prot_gbsw.inp
  read para card unit 10
  close unit 10
 
  ...

  ! read in the new input radii
  prnlev 0
  stream @toppar/radius_gbsw.str
  prnlev 5

  ! verify that all heavy atoms have non-zero radii
  scalar wmain statistics select .not. type H* end
  define check select (.not type H* ) .and. ( prop wmain .eq. 0.0 ) show end
  if ?nsel ne 0  goto diehard      !some heavy atom have a zero radius

  ! ativate GBSW energy term
  gbsw sgamma 0.005 nang 50

  nbond atom switch cdie vdw vswitch -
        ctonnb 16 ctofnb 16 cutnb 20

  energy

  ...

  (minimization and/or dynamics) 

