.. py:module:: gbim


Generalized Born Solvation Energy Module with Implicit Membrane
---------------------------------------------------------------

GBIM is a modification of the  GENBORN module  that includes
Implicit Membrane in the calculations of the  electrostatic contribution
to solvation energy. The non-polar region of the membrane is approximated
as a planar dielectric slab having the same dielectric constant as inside
the molecule.  It permits the calculation of the Generalized Born solvation
energy and forces following the formulation of the Qui & Still pairwise
GB approach in  linearized version of B. Dominy and C.L. Brooks, III
(see genborn.doc).

The Generalized Born model with Implicit Membrane is described in
Spassov et al., 2002 (see below).

In the  GBIM  module the polarization energy is computed following the
equation:

.. math::

   G_{pol} = -C_{el} \left( \frac{1}{\epsilon_m} - \frac{1}{\epsilon_{slv}} \right) \left[ \frac{1}{2} \sum_{i=1}^{N} \sum_{j=1}^{N} \frac{q_i q_j}{ ( r_{ij}^2 + \alpha_i \alpha_j e^{-D_ij} )^\frac{1}{2} } \right]
   
:math:`\epsilon_m` is the dielectric constant of the reference medium and :math:`\epsilon_{slv}` is 
the dielectric constant of the solvent. 

If the membrane is present, the effective Born radii are calculated as:

.. math::

   \alpha_i = - \left( \frac{1}{\epsilon_m} - \frac{1}{\epsilon_{slv}} \right) \frac{ C_{el} }{ 2G_{pol,i} }
        
where

.. math::

   G_{pol,i} = \left( \frac{1}{\epsilon_m} - \frac{1}{\epsilon_{slv}} \right) \left[ GAM( Z(i), R(i), L, \lambda, \gamma ) 
                      + \sum_\mathrm{bonds} \left( P2 \cdot \frac{V(j)}{r_{ij}^4} \right)
                      + \sum_\mathrm{angles} \left( P3 \cdot \frac{V(j)}{r_{ij}^4} \right)
                      + \sum_\mathrm{non-bonded} \left( P4 \cdot V(j) \cdot \frac{C_{ij}}{r_{ij}^4} \right) \right]
  
The self term  (1/eps_m -1/eps_slv)* GAM(i)  approximates the polarization
energy of a single ion in the presence only of membrane.  Z(i) is the 
distance of the atom from the membrane midplane and L is the membrane
thickness.

::

            slv     cntr    slv
  GAM(i) = g    + (g     - g    ) / {1 + exp[ Gamma * (Z(i) + R(i) - L/2)] }
            i       i       i

where

::

     slv
    g    = -Cel/(2*Lambda*R(i)) + P1*Cel/(2*R(i)^2 
     i

and

::

     cntr     Cel ln(2)
    g    = - ----------
     i           L
   

The rest of the variables are the same as in :doc:`genborn`.
 
The gradient of polarization energy is also computed so GBIM can be used
in energy minimization and dynamics.

The combined use of GBIM and ASPENRMB (:doc:`aspenrmb`) can be used for
calculations of solvation energy in the frames of GBSA/IM
(Generalized Born - Surface Area model with Implicit membrane) approach.


REFERENCES:

V.Z. Spassov, L. Yan and S. Szalma. Introducing an Implicit Membrane in
Generalized Born / Solvent Accessibility Continuum Solvent Models.
J. Phys. Chem. B, 106,8726-8738 (2002)         

.. _gbim_syntax:

Syntax
------

Syntax of the Generalized Born model with Implicit Membrane commands

::

   GBIM { P1 <real> P2 <real> P3 <real> P4 <real> P5 <real> 
           [LAMBda <real>] [EPSILON <real>] [EPSMOL <real>]
           [TMEMB <real>] [ZMDIR (or XMDIR or YMDIR) ] [CENTER  <real>] 
           [CUTAlpha <real>] [WEIGht] [ANALysis] }
           { CLEAr                                                     }


.. _gbim_function:

Parameters of the Generalized Born Model with Implicit Membrane
---------------------------------------------------------------

========= ===========================================================================
P1-P6     The parameters P1, P2, ..., P5 specify the particular parameters
          controlling the calculation of the effective Born radius for
          a particular configuration of the biomolecule.
       
          ::

             Alpha(i) = [  GAM( Z(i),R(i), L, Lambda, Gamma )

                      + Sum {P2*V(j)/rij^4} + Sum {P3*V(j)/rij^4}
                       bonds                 angles

                      + Sum {P4*V(j)*Cij/rij^4}]^(-1)*(-CCELEC)/2
                     non-bonded

             with Cij = 1 when (rij^2)/(R(i)+R(j))^2 > 1/P5

             and Cij = 1/4(1-cos[P5*PI*(rij^2)/(R(i)+R(j))^2])^2 otherwise.


          .. note::
             R(i), V(i) correspond to the vdW radius and volume respectively,
             CCELEC is the conversion from e^2/A to kcal/mol, rij is the separation
             between atom i and atom j.

Lambda    This is the scaling parameter for the vdW radius in the :math:`g^{slv}`
          term
          of GAM function. It has the  same meaning, as in genborn.doc.

          .. note::
             The parameters P1-P5 and Lambda correspond to parameters for a 
             particular CHARMM parameter/topology set.

          .. warning::
             The parameters P1-P5 and Lambda are required input

EPSILON   This is the value of the dielectric constant for the solvent medium.
          The default value is 80.0

EPSMOL    This is the value of the dielectric constant for the
          reference medium.  The default value is 1.0.

TMEMB     Membrane thickness  

ZMDIR     Membrane normal is along Z axis (or XMDIR or YMDIR)

CENTER    Position of membrane midplane  ( Z coordinate, if ZMDIR)  
     
Gamma     Empiric parameter regulating the slope of GAM function
          A good accuracy  for the charmm19 force field is achieved
          with Gamma =  0.55 [A^(-1)].

CUTAlpha  This is a maximum value for the effective Alpha for any atom
          during the calculation for a particular conformation of the
          biomolecule.  It is necessary because in some instances the
          expression above for Alpha(i) can take on negative values
          of numerical problems with the expression for very buried atoms
          in large globular biopolymers.  The default for this value is
          10^6.

WEIGht    This is a flag to specify that you want the vdW radii for the
          atoms to be taken from the wmain array instead of the parameter
          files (from Rmin values).  The default is to use the parameter
          values.  These values are used for the R(i) and V(i) noted above.

ANALysis  This flag turns on an analysis key that puts the atomic contributions
          to the Generalized Born solvation energy into an atom array (GBATom)
          for use through the scalar commands.

CLEAr     Clear all arrays and logical flags used in Generalized Born 
          calculation.
========= ===========================================================================

.. _gbim_examples:

Examples
--------

The examples below illustrate some of the uses of the generalized Born
model with charmm19.  See also c31test/gbsaim.inp

Example (1) 
^^^^^^^^^^^

Calculates the generalized Born solvation energy using atomic radii
from the wmain array (the example illustrates the useage but simply
uses the same radii as would be employed w/o the "Weight" option). The
membrane is present as a 30. Angstrom  dielectric slab. The membrane
normal is along Z and membrane midplane has a coordinate Z = 0. 
A value of 2.0 is used for the  molecular & membrane dielectric constant
and  80. for the water solvent.

::

   !  Test use of radii from wmain array
   scalar wmain = radii
   !  Now turn on the Generalized Born energy term using the param19 parameters
   Gbim P1 0.415 P2 0.239 P3 1.756 P4 10.51 P5 1.1 Lambda 0.730 -
        Epsilon 80.0 Epsmol 2. -
        Tmemb 30. Zmdir  Center 0.0  Gamma 0.55  Weight


   ! Now calculate energy w/ GB
   energy cutnb 20 ctofnb 16 ctonnb 14

   GBIM Clear


Example (2)
^^^^^^^^^^^

Use of the ANALysis key to access atomic solvation energies.

::

   Gbim P1 0.415 P2 0.239 P3 1.756 P4 10.51 P5 1.1 Lambda 0.730 -
        Epsilon 80.0 Epsmol 2. -
        Tmemb 30. Zmdir  Center 0.0  Gamma 0.55   Analysis
   energy cutnb 990

   !  What are the current Generalized Born Alpha, SigX, SigY, SigZ and T_GB 
   !  and atomiuc solvation contribution (GBATom) values?
   skipe all excl GbEnr
   energy cutnb @cutnb
   scalar GBAlpha show 
   scalar SigX show 
   scalar SigY show 
   scalar SigZ show 
   Scalar T_GB show 
   Scalar GBAtom show  ! One can now use these individual contributions
   GBIM Clear


Example (3)
^^^^^^^^^^^

Do a minimization (could be dynamics too)

::

   !  Minimize for 1000 steps using SD w/ all energy terms.
   skipe none

   Gbim P1 0.415 P2 0.239 P3 1.756 P4 10.51 P5 1.1 Lambda 0.730 -
        Epsilon 80.0 Epsmol 2. -
        Tmemb 30. Zmdir  Center 0.0  Gamma 0.55   Analysis


   mini sd nstep 1000 cutnb 20 ctofnb 18 ctonnb 18 -
   elec cdiel Eps 2 switch

Note, that Eps must be equal to Epsmol for consistent results!
