.. py:module:: facts

=======================================================
FACTS: Fast Analytical Continuum Treatment of Solvation
=======================================================

.. note::

   **Questions and comments regarding FACTS should be directed to**

   Amedeo Caflisch (caflisch@bioc.uzh.ch)

   **Reference for FACTS:**

   1. Haberthuer and Caflisch, J. Comput. Chem., 29(5): 701-715, 2008
      DOI: 10.1002/jcc.20832

.. _facts_description:

Description
-----------

FACTS is an efficient generalized Born implicit solvent model [1].Because
of its speed it is particularly useful for MD simulations. It is based on
the fully analytical evaluation of the volume and spatial symmetry of the
solvent that is displaced from around a solute atom by its neighboring 
atoms. The two measures of solvent displacement are combined in empirical
equations to approximate the atomic (or self) electrostatic solvation
energy and the solvent accessible surface area. The former directly yields
the effective Born radius, which is used in the generalized Born formula
to calculate the solvent-screened electrostatic interaction energy.
The solvent accessible surface area is used to approximate the non-polar
contribution to solvation. FACTS is only four times slower than using the
vacuum energy in molecular dynamics simulations of peptides and proteins.

FACTS can be used with either force field "param19" or "param22/27", but
parameters were derived for protein atoms only. Parameters for unknown 
atom types can be extrapolated from known ones with the "TAVW" option
(see below).

Periodic boundary condition can be set up with the CHARMM IMAGE facility
and the FACTS "TAIM" option.

The FACTS algorithm has been parallelized in 2007 by F.Marchand
(fmarchand@bioc.uzh.ch).

.. index:: facts; syntax
.. _facts_syntax:

Syntax of the FACTS Command
---------------------------

::

   FACTs {TCPS 19} {TEPS 1.0} [TKPS real] [GAMM real] 
         {TCPS 22} {TEPS 2.0} [TAIM] [TAVW] [TPSL] 

.. _facts_function:

Description of the FACTS keywords and options
---------------------------------------------

=======  =======  =============================================================
Keyword  Default  Purpose
=======  =======  =============================================================
TCPS        19    Force-field specific, either 19 or 22. TCPS 22 also works
                  with param27, although only protein atoms are parametrized.

TEPS       1.0    Solute dielectric constant. One can specify only TEPS = 1.0
                  or TEPS = 2.0, since only these values were used as interior
                  dielectric constant for the finite-difference Poisson
                  calculation of atomic solvation energies used to
                  parametrize FACTS.

TKPS       4.0    Kappa, the factor in the denominator of the exponential
                  function in the Generalized Born formula. It usually
                  ranges between 4.0 (Still's original equation) and 8.0.

GAMM      0.015   Nonpolar surface tension coefficients in kcal/(molxA^2).

TAIM      false   Flag to include contribution from image atoms. The FACTS
                  module works with the CHARMM IMAGE facility. This flag
                  should be included in the FACTS command line if images (PBC)
                  are used.

TAVW      false   FACTS distinguishes atoms solely based on their VdW radius
                  and only atoms found in proteins are currently parametrized.
                  If a system contains atoms whose radii are not among those
                  already parametrized in FACTS, then the "TAVW" option will
                  extrapolate new parameters from already parametrized atoms.

TPSL      false   Flag to print detailed information about the atomic
                  self-electrostatic contribution to the energy.
=======  =======  =============================================================


.. _facts_examples:

Examples
--------

PARAM 19
^^^^^^^^

REMARK: It was observed that (for param19) FACTS performs better with
        the dielectric constant set to 2.0 rather than to 1.0.

::

   ! ------------------------------------------------------------------
   !
   !  large, globular systems (e.g., protein in folded state)
   !  epsilon=2.0 and gamma=0.025

   set diele 2.0

   nbond nbxmod 5 atom cdiel eps @diele shift vatom vdistance vshift -
         cutnb 9.0 ctofnb 7.5 ctonnb 6.5 e14fac 0.4 wmin 1.5

   scalar wmain = radius
   scalar wmain set 1.0 selection (type h*) end

   facts tcps 19 teps @diele gamm 0.025

   ! ------------------------------------------------------------------
   !
   !  Reversible folding simulations of structured peptides
   !  epsilon=2.0 and gamma=0.015

   set diele 2.0

   nbond nbxmod 5 atom cdiel eps @diele shift vatom vdistance vshift -
         cutnb 9.0 ctofnb 7.5 ctonnb 6.5 e14fac 0.4 wmin 1.5

   scalar wmain = radius
   scalar wmain set 1.0 selection (type h*) end

   facts tcps 19 teps @diele gamm 0.015

   ! ------------------------------------------------------------------
   !
   !  Unstructured peptides and peptide aggregation 
   !  epsilon=2.0 and gamma=0.0075

   set diele 2.0

   nbond nbxmod 5 atom cdiel eps @diele shift vatom vdistance vshift -
         cutnb 9.0 ctofnb 7.5 ctonnb 6.5 e14fac 0.4 wmin 1.5

   scalar wmain = radius
   scalar wmain set 1.0 selection (type h*) end

   facts tcps 19 teps @diele gamm 0.0075

   ! ------------------------------------------------------------------
   !
   !  Periodic Boundary Conditions
   !  (images should be previously set up with the IMAGE facility)
   !  Print detailed informations about the atomic self-electrostatic 
   !  contributions to the energy
   !  epsilon=1.0 and gamma=0.015

   set diele 2.0

   nbond nbxmod 5 atom cdiel eps @diele shift vatom vdistance vshift -
         cutnb 9.0 ctofnb 7.5 ctonnb 6.5 e14fac 0.4 wmin 1.5

   scalar wmain = radius
   scalar wmain set 1.0 selection (type h*) end

   facts tcps 19 teps @diele gamm 0.015 taim tpsl


PARAM 22
^^^^^^^^

::

   !  epsilon=1.0 and gamma=0.015

   set diele 1.0

   nbond nbxmod 5 atom cdiel eps @diele shift vatom vdistance vswitch -
         cutnb 14.0 ctofnb 12.0 ctonnb 10.0 e14fac 1.0 wmin 1.5

   scalar wmain = radius

   facts tcps 22 teps @diele gamm 0.015


See also: test cases 

* ~/charmm/test/c35test/facts_p19.inp
* ~/charmm/test/c35test/facts_p22.inp
