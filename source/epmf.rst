.. py:module:: epmf

=====================
EPMF Module of CHARMM
=====================

EPMF module implements an empirical energy term based on
an angle and distance (HBPMF) or a distance (1DPMF).  1D PMF is defined
as function of distance between two atoms(DONR and ACCP), whereas HBON PMF
is function of distance (DONR-ACCP) and cosine of angle (DONR-HYD-ACCP).

::

                      (H)                                  (H)
         |             |                     |              |
         |       {B}  d|    {A}              |         BLEN |
         |             |    a                |              |
   ----(accp).......(donr)-----(atom2)     (accp)........(donr)
         |             |                     |      F1   /     \
         |             | b  {D}              |      F2  /       \
         |             |                     |         /         \
                     (atom1)                        (atom1)    (atom2)

                HBON/DEFA                           HBON/GEOM

A sketch of Hydrogen atom construction and interactions is shown above.
Few comments on the notation.

(1) 1D PMF is function of distance between (donr) and (accp) atoms.

(2) Empirical hydrogen bonding provided by the DEFA option is applicable
    only when atom2, atom1 and donr atoms are planar.Then position of the hydrogen
    atom bonded to donr can estimated using the relation

    ::

                        _                                _
                 -|d|  |   sin(B)          sin(A)         |
      vec(d) =  -----  |   -----  vec(a) + ------- vec(b) |
                sin(D) |_   |a|             |b|          _|

   where A,B,D are the angles defined by (H-donr-atom2), (H-donr-atom1)
   and (atom1-donr-atom2) respectively.

   A more flexible way of estimating hydrogen atom is also provided by the GEOM
   option which uses the bonding information: hydrogen atom-donr bond, hydrogen
   atom-donr-atm1 angle and hydrogen atom-donr-atm1-atm2 dihedral.


.. _epmf_syntax:

Syntax of the EPMF Command
--------------------------

::

  EPMF  [UPFR] int [CLEA]
        {DIST   dist-option-spec dist-PMF-spec selection-spec}
        {HBON   ATM1 atom-selection F1 real
                ATM2 atom-selection F2 real
                BLEN real hbon-option-spec  hbon-PMF-spec selection-spec}

  dist-option-spec ::= {  PMF0 unit  [PMF1 unit] [PMF2 unit] [PMF3 unit]
                         [PMFN unit] [PMFM unit] [PMFP unit]
                         [SCL0] real [SCL1] real [SCL2] real [SCL3] real
                         [SCLN] real [SCLP] real [SCLM] real
                       }

  hbon-option-spec ::= {  DEFA defa-specs
                          GEOM geom-specs
                         [PMF1 unit] [PMF2 unit] [PMF3 unit] PMFN unit
                         [SCL1] real [SCL2] real [SCL3] real [SCLN] real
                         [EMIN] real
                       }

  selection-spec :: = { DONO atom-selection ACCP atom-selection}


.. _epmf_function:

Purpose of the various EPMF variables
-------------------------------------

+---------+-------------------------------------------------------------------+
|Variable | Explanation                                                       |
+---------+-------------------------------------------------------------------+
|CLEA     | clears the memory allocated for EPMF module data structures.      |
+---------+-------------------------------------------------------------------+
|DIST     | evaluates interaction energy based on 1D PMF for the              |
|         | given selections The donor and the acceptor atoms are             |
|         | specified by keywords DONOr and ACCEptor respectively.            |
+---------+-------------------------------------------------------------------+
|HBON     | evaluates interaction energy based on 2D PMF for the              |
|         | given selections. This option in addition to DONOr atom,          |
|         | needs two other atoms (ATM1, ATM2)for estimating position of      |
|         | H-atom bonded to the Donor atom. Further three additional         |
|         | real arguments are needed. For DEFA option these factors are      |
|         |                                                                   |
|         | ::                                                                |
|         |                                                                   |
|         |   F1 =  ( -|d|/sin(D) * sin(B)/|a| )  (see the sketch             |
|         |   F2 =  ( -|d|/sin(D) * sin(A)/|b| )   in top section)            |
|         |   BLEN = |d|                                                      |
|         |                                                                   |
|         | Note that the hydrogen atom construction for this case remains    |
|         | accurate only when angles A, B and D are around 120 degrees.      |
|         |                                                                   |
|         | To alleviate this issue, an additional robust method GEOM, which  |
|         | relay on bonding features is provided. For GEOM, option F1 becomes|
|         | the equilibrium angle Hydrogen-Donor-ATM1 and F2 becomes the      |
|         | equilibrium dihedral Hydrogen-Donor-ATM1-ATM2 . The BLEN option   |
|         | remains the same as in previous case.                             |
+---------+-------------------------------------------------------------------+
|PMF0     |                                                                   |
|PMF1     | PMF0..N  specify the PMF data file used for                       |
|PMF2     | the (i, i+n(n=1..N)) residue interactions.                        |
|PMF3     | See PMF-file section for the format description of                |
|PMFN     | these files. PMF0, PMFP(i+1 residue) and PMFM(i-1 residue)        |
|PMFP     | are valid only for DIST1                                          |
|PMFM     |                                                                   |
+---------+-------------------------------------------------------------------+
|UPFR     | specifies the update frequency for the (i,i+n)                    |
|         | donor-acceptor list. Default value is 25                          |
+---------+-------------------------------------------------------------------+
|SCL0     |                                                                   |
|SCL1     |                                                                   |
|SCL2     | Scaling of PMF0..N energies. Default scaling factors are 1.0      |
|SCL4     | SCL0, SCLM and SCLP are valid only for DIST PMF                   |
|SCLN     |                                                                   |
|SCLM     |                                                                   |
|SCLP     |                                                                   |
+---------+-------------------------------------------------------------------+
|EMIN     | specifies the minimum allowed energy of EPMF interactions defined |
|         | per residue. If the EPMF energy for particular residue becomes    |
|         | less than the desired value, the subsequent EPMF calculations     |
|         | for that residue is reset to ZERO. By default EMIN is set to      |
|         | -10 kcal/mol                                                      |
+---------+-------------------------------------------------------------------+
|ATM1     | These selections are used only for HBON PMF. The selections       |
|ATM2     | are specified by atom name, prefixed by '+', '-' or ''.           |
|         | If DONO selection belongs to jth residue, than '+' indicates      |
|         | ATM1/ATM2 belongs to j+1th residue, '-' ATM1/ATM2 belongs to      |
|         | j-1th residue. If prefix is '', ATM1/ATM2 belong to jth           |
|         | residue                                                           |
+---------+-------------------------------------------------------------------+

::

  atom-selection:== (see *note select:(chmdoc/select.doc).)

.. _epmf_pmf_file

PMF-File
--------

The PMF-file must contain either two(1D PMF/DIST) or three(HBOND) columns.
The last column always correspond to energy and other column(s) are the grid
points of PMF co-ordinate(s). In addition to data, the PMF-file must contain
a header section describing type of PMF, maximum and minimum values of grid
points and number of X/Y entries in the file.

For 1D PMF the header section of PMF-file must be

::

  <XMAX>     <XBINS>   <XMIN>    <XPTS>
  4.50   0.1    0.00
  0.00   0.00
  0.1    0.05
  ...
  ...
  ...
  4.49   0.95
  4.50   1.00

For HBOND PMF the header section is

::

  <XMAX>     <XBINS>   <XMIN>    <YPTS>
  <YMAX>     <YBINS>   <YMIN>    <XPTS>
  4.50   0.1    3.00    21
  1.0    0.05   0.00    16
  3.00   0.00   0.00
  3.00   0.05  -0.15
  ...
  ...
  ...
  4.50   0.95   0.15
  4.50   1.00   0.00

.. _epmf_example:

Examples
--------

* Example (1)

  ::

    EPMF DIST PMF3 19 PMFN 23 DONO select type (OE* .or. OD*) end ACCP select
    ND .or. NE* end

  A simple distance based PMF to mimick salt bridges for some polar and charged
  residues. The PMF data for i+/-3 interactions are read from unit 13 and i+/-n
  from unit 23

* Example (2)

  ::

    EPMF HBON defa atm1 CA f1 -.6736 atm2 -C f2 -.7627 blen .997 -
         PMF1 13 PMFN 17  DONO select type N end ACCP select type O end

  This example calculates putative hydrogen bond interactions between backbone O
  and N atoms. Position of hydrogen is estimated using DEFA option of HBOND
  potential. Additional atoms CA (ith residue) , -C (i-1th residue) and
  corresponding factors F1,F2,BLEN  are needed for estimation of position of
  hydrogen atom bonded to N of (ith residue). Also needed are the PMF data files
  given by unit 13 (for i+/-3) and unit 17 ( for i+/-n, where n>3).


* Example (3)

  ::

    EPMF HBON geom atm1 CA atm2 -C  f1 116.0 f2 180.0  blen 0.997 -
         PMF1 13 PMFN 17  DONO select type N end ACCP select type O end

  Same as Example(2), but using GEOM option of HBON for estimating Hydrogen
  atom position.

Application of DIST and HBON PMFs in context of PRIMO force field is given as
test case of charmm c36a451 version.

