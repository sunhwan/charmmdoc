.. py:module::rdfsol

============================
Radial Correlation Functions
============================


The RDFSOL command computes radially resolved correlation
functions, such as radial distribution functions or orientational
correlation functions. The function of interest is computed either between
pairs of atoms from two atom selections, or between a pair consisting of
one atom selection and a reference point, which can be either a fixed point
in space or the center of mass of a set of atoms. As the name rdfSOL
suggests, the routine allows for the special treatment of solvent molecules
(TIP3 water is supported by special routines, others need the use of the
PROTotype facility (see see :ref:`Selection <proto_selection>`). 


.. _rdfsol_syntax:

Syntax for the RDFSOL command
-----------------------------

::

   RDFSOL  [ RDF int ]  [ DDIP int ]  [ QDIP int ]  [ HD int ] -

                     [ setA-spec ]  [ setB-spec ]  [ around-spec ] -

                     [ SAME ]  [ RMAX real ]  [ NBIN int ] -

                     [ VOLUme real ]  [ PRECise ]  [ BRUTe ]  [ MINI ]  -

                     [ traj-spec ]  [ SPFAc int]



         setA-spec:: [ SITE ]  [ atom-selection ]
                     [ XREF real ]  [ YREF real ]  [ ZREF real ]
                     [ WATEr ]
                     [ PROTo int ]


         setB-spec:: [ SITE ]  [ atom-selection ]
                     [ WATEr ]               
                     [ PROTo int ]


         around-spec:: AROUnd  [ RAROund real ]  [ LOCAl ] -

                            [ atom-selection ]
                            [ XREF real ]  [ YREF real ]  [ ZREF real ]


         traj-spec:: [ FIRStu int ] [ NUNIt int ] [ BEGIn int ] -

                     [ STOP int ] [ SKIP int ]


         atom-selection::= see *note Selection:(chmdoc/select.doc)

.. _rdfsol_general:

General overview and options
----------------------------

RDFSOL calculates radially resolved pair-distribution or angular
correlation functions between two sets of atoms (setA and setB).
The type of function computed is selected by keywords:

===== ================================================================
RDF   Radial Distribution Function
      If one of the two sets is WATEr then distribution functions
      for the oxygen and hydrogens are computed.
      If both sets are WATEr then the hydrogen-hydrogen distribution
      function will also be evaluated.
      If the first set is not water and the keyword SITE is present,
      the center of mass of the set is taken as single center. Else
      the average over all points in setA is taken.
     
DDIP  Dipole-Dipole correlation function. If one or both sets are
      not WATEr the center of mass and dipole moment of this set
      is used (no matter whether the keyword SITE is present or
      not). In this case setA must not be a fixed point in space
      since the dipole moment is not defined in this case.
     
QDIP  Charge-Dipole correlation function. As with RDF setA can
      either be a SITE or the average of all points.
===== ================================================================

The integer after each function to be calculated gives the unit number
the respective function is to be written to.


.. _rdfsol_sets:

SiteA/B specifications
----------------------

SetA and SetB are two sets for which the chosen function is evaluated
for all pairs A-B. Both can be WATEr, in which all TIP3 residues present
will be included. In this case the oxygen positions will be used as the
centers of the molecules. In both cases the center of mass (and set
dipole moment if needed) can be used if the keyword SITE is present.
SetA can be a fixed point in space: (XREF/YREF/ZREF) (if SITE is present
but no atom selection (0/0/0) will be used as default).
Finally for each of the two sets a previously defined prototype set
(see :ref:`Prototypes <proto_prototypes>`) can be used. In this case
the center of geometry (or mass with keyword MASS) and dipole of each
individual set member will be used in the requested functions.
For both sets WATEr is the default.

.. _rdfsol_limit:

Limiting Sets
-------------

If only a subset which is localized around a certain point should be
used in each frame this can be achieved by the AROUnd keyword. If it is
present setA will be re-selected in each frame. If the keyword LOCAl is
present setB will also be re-selected. RAROund <real> is the radius
around the selected center within which an atom must lie to be available
for evaluation in this frame. The center itself can either be a fixed
point in space (XREF/YREF/ZREF) or an atom selection of which the center
of mass will be used.

.. _rdfsol_options:

Other options
-------------

=============== =========================================================
SAME            if this keyword is present, only setA is used for both
                sets, thus calculating auto-functions (this algorithm
                should be faster than the general one if setA and setB
                use the same selection)

RMAX <real>     the maximum distance up to which a pair A-B is evaluated
                (default: 7.5A)

NBIN <int>      the number of bins used to sample (each bin is RMAX/NBIN
                wide)
                (default: 150)

VOLUme <real>   the volume of the total system, necessary for the
                normalization. If not specified by the user and crystal
                is in use, the resulting cell volume will be used.
                Finally, if crystal is not used and no volume
                is specified, and if both sets are localized (see AROUnd), 
                the volume of the limiting sphere will be used.

PRECise         if RDFs are calculated and one or both sets contain
                WATEr, some pairs including water hydrogens will be
                missed since only oxygen distances are evaluated.
                If PRECise is present, these pairs are also included
                which results in a slightly diminished efficiency of the
                cubing algorithm

BRUTe           use a simple double loop algorithm rather than a cubing
                algorithm

MINI            use 'real' minimum image conventions. Currently only one
                function can be calculated if MINI is used. Its major use
                is the computation of the distance dependent Kirkwood
                G-factor (with DDIP, second column). Here, one needs to go
                'into the corners' (i.e. sqrt(3)/2 * L for a cubic box)
                without counting pairs twice.
                (caution: needs lots of memory)

SPFAc           if images are present, the number of total
                atoms/pairs/cubes may change from frame to frame. 
                So an estimate of the needed space needs to be made
                before reading the trajectory so SPFAc times the actual
                values is allocated
                (default: 3)
=============== =========================================================

.. _rdfsol_traj:

Trajectory specifications
-------------------------

These are the usual specs. The trajectory is read NUNIt units starting
with FIRSTu reading from frame BEGIn to STOP where SKIP frames are
skipped between reading.

.. _rdfsol_caveats:

Caveats and Comments
--------------------

- When computing dipole-dipole correlations for a set which is not WATEr,
  only its center of mass and dipole moment for primary atoms will be
  evaluated. So if a part of a large molecule which is re-centered
  bysegment (e.g. a protein) and "sticks out" of the primary box, some
  pairs may not be sampled.

- Normalization of RDFs differs slightly from that used in COOR ANAL

- no excluded volume correction

- point-point (e.g. two SITEs or DDIP for two non-WATEr sets...) not yet
  implemented

.. _rdfsol_examples:

Examples
--------

(See also test/c30test/rdfsol.inp test/c30test/rdfsol2.inp testcases)

::

   RDFSOL RDF 10 SETA WATER SAME RMAX 7.5 NBIN 150 PRECISE -
          FIRSTUNIT 11 NUNIT 1

This will calculate g_OO, g_OH and g_HH for all waters in the simulated
system up to 7.5 A into 150 bins. One trajectory file will be read from
unit 11 and the result output to unit 10.

----------------------------------------------------------------------

::

   RDFSOL RDF 10 SETA WATER SAME RMAX 7.5 NBIN 150 PRECISE -
          AROUND RAROUND 7.5 LOCAL SELECT ATOM PROT 1 NH END    -
          FIRSTUNIT 11 NUNIT 1

The same as above but only waters around the NH of residue 1 of segment
PROT will be considered.

----------------------------------------------------------------------

::

   RDFSOL RDF 10 QDIP 11 DDIP 12 SETA WATER SAME RMAX 7.5 NBIN 150 PRECISE -
          FIRSTUNIT 13 NUNIT 1

Same sets as in the first example but here all three functions are
calculated at once (i.e. the trajectory is only read once).

