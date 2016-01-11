.. py:module:: gopair

=======================================
The Go pairwise Energy Module of CHARMM
=======================================

::

                         By Charles L. Brooks III, 2014


The Go Pair facility in CHARMM was added to complement and extend the
functionality of the Karinoclas and Brooks Go model implementation in
CHARMM (ETEN) to more easily accomodate the situation where individual
protein/nucleic chains are treated as KB-specific Go models and the
inter-chain interactions are modeled as either KB-like Go interactions
or some generaic form of coarse-grained interaction model, e.g., if
one wants to consider the interactions of two proteins of known
structure, for which one can build a KB Go model, and have the
inter-chain/protein interactions occur with a general Mizawa-Jernigan
pairwise Ca-based interaction. In this case, the GoPair facility
accounts for the non-generic intra-protein/chain interactions via the
pairwise modified MJ KB Go interactions and the generic repulsive
terms within the chain are treated via the normal non-bonded routines
as are the inter-chain interactions.

.. note::

   In this implementation the intra-chain KB-based Go interactions
   are not subject to periodic boundary conditions and all pairwise
   interactions are considered, i.e., no cutoffs. Whereas, the
   non-Go/non-specific interactions are treated via periodic boundary
   conditions with cutoffs.

The ETEN functionality is available and can be turned on or off
independently from that govering the non-specific pairwise
interactions. Thus one can do "mixed" models where intra-protein
interactions are treated with the KB form and inter-chain are treated
with standard LJ, or visa-versa (remember, the replusive intra-chain
interactions are treated via the normal ETEN ON/OFF model here.

How it works
------------

The Go Pair functionality reads a list of pairwise interaction
parameters, i.e., emin and rmin ananlgous to the NBFIX pairwise
parameters, and sets up the list structure to calculate the energy and
forces on the specified pairs using the specified functional form
(ETEN ON/OFF). These pairwise interactions are put on the non-bonded
exclusion list so they are not double counted. In doing so, as noted
above, the set of pairwise interactions specified by the GOPAIR
commands ARE NOT SUBJECT TO CUTOFFS OR PEREIODIC BOUNDARY CONDITIONS,
and hence are only appropriate for single chains or multi-chain
systems that are confied by some means other than periodic boundary
conditions.

References:
(1) Karanicolas & Brooks, Protein Science, 11, 2351 (2002).

.. _gopair_syntax:

Syntax
------

::

   [INPUT GoPair command]

     GOPAIR [READ UNIT integer] [ON/OFF] [CLEAR] [ETEN ON/OFF]


.. _gopair_function:

Function
--------

::

   READ  : Read the Go-model pairwise interactions from unit <integer>, 
           default is to read from input stream

   ON/OFF : Turn the Go Pair functionality on or off (note reading the
            parameters automatically turns it on and clearing the
            data structures automativally turns it off.

   CLEAr : Turn off GoPair functionality and delete data structures

   ETEN  : Turns the KB ETEN functional form on <ON> or off <OFF>. 
           Default is ETEN OFF.


.. _gopair_example:

Example usage
-------------

Consider a system comprised of two chains (A & B) for which a standard
KB Go model has been constructed for each chain A and chain B (Note:
in the near future the Go Server will produce an additional file for
the GoPair intrachain interactions separately to facilitate usage of
the GoPair modeling facility, let's assume it's called
model_kb-gopair.param.) and the "generic" intra-chain repulsive
interactions, the generic Mizawa-Jernigan or otherwise determined
inter-chain non-bonded parameters (specified as NBFIXes for the
N_A*(N_B-1)/2 pairs of interacting sites), the "standard" NBFIX KB
parameters for intra-chain A-A and B-B, bonded, angle and torsion
parameters are contained in the parameter file model_kb-go.param. The
topology is described in model_kb-go.top and the initial Ca-only
structures are in model_a.pdb and model_b.pdb.

::

   Files assumed: model_kb-go.top - contains the required topology information
                  model_kb-go.param - contains the associated parameters
                  model_a.pdv, model_b.pdb - contains coordinates for chains
                  A and B.
                  model_kb-gopair.param - contains the intra-chain non-bonded
                  KB Go parameters.


   ***************
   model_kb-go.top:
   * Topology for Go model of 1arr
   *
      20   1
   MASS 1   A1       131.000000
   MASS 2   A2       128.000000
   MASS 3   A3       57.000000
            .
            .
            .
   MASS 103 B50     157.000000
   MASS 104 B51     113.000000
   MASS 105 B52     57.000000
   MASS 106 B53     71.000000

   DECL +CA
   AUTOGENERATE ANGLE DIHEDRAL

   RESI A1         0.0
   GROU
   Atom  CA  A1       0.0
   Bond CA +CA

   RESI A2         0.0
   GROU
   Atom  CA  A2       0.0
   Bond CA +CA
            .
            .
            .
   RESI B52       0.0
   GROU
   Atom  CA  B52     0.0
   Bond CA +CA

   RESI B53       0.0
   GROU
   Atom  CA  B53     0.0
   Bond CA +CA

   END

   *****************
   model_kb-go.param:
   * Parameters for Go model of 1arr
   *

   BOND
   A1      A2        378.000000  3.841480
   A2      A3        378.000000  3.832207
            .
            .
            .
   B50    B51      378.000000  3.817184
   B51    B52      378.000000  3.740986
   B52    B53      378.000000  3.807921

   ANGLE
   A1      A2      A3         75.600000 87.002943
   A2      A3      A4         75.600000 92.769692
            .
            .
            .
   B49    B50    B51       75.600000 81.795665
   B50    B51    B52       75.600000 100.811114
   B51    B52    B53       75.600000 98.404340

   DIHEDRAL
   A1   A2   A3   A4       0.070661 1  148.427948
   A1   A2   A3   A4       0.642645 2  247.750476
   A1   A2   A3   A4       0.131763 3  98.732133
   A1   A2   A3   A4       0.076565 4  20.955060
   A2   A3   A4   A5       0.155810 1  253.803724
   A2   A3   A4   A5       0.433367 2  21.748974
   A2   A3   A4   A5       0.116055 3  221.349291
   A2   A3   A4   A5       0.169406 4  13.133496
            .
            .
            .
   B50 B51 B52 B53     0.083619 1  189.775808
   B50 B51 B52 B53     0.746189 2  228.220978
   B50 B51 B52 B53     0.184920 3  106.554957
   B50 B51 B52 B53     0.071306 4  353.723408

   NONBONDED NBXMOD 3 ATOM CDIEL SWITCH VATOM VDISTANCE VSWITCH -
     CUTNB 399.0 CTOFNB 398.5 CTONNB 395.5 EPS 1.0 WMIN 1.5

   A1     0.0  -0.000546  2.474456  !The following are generic repulsive
   A2     0.0  -0.000012  4.738648  ! interactions beween atoms.
   A3     0.0  -0.000224  6.214139
            .
            .
            .
   B50     0.0  -0.000155  3.489662
   B51     0.0  -0.000654  2.905577
   B52     0.0  -0.000224  3.193040
   B53     0.0  -0.000272  4.391179

   NBFIX
   A1      A4         -2.36758  3.843863  ! These are the intra-chain KB Go
   A1      A6         -0.924458  4.937081 ! parameters for chain A
   A2      A5         -0.391404  6.159995
   A14     A18        -2.522816  6.612486
            .
            .
            .
   B1     B4        -2.36758  3.843863  ! These are the intra-chain KB Go
   B1     B6        -0.924458  4.937081 ! parameters for chain B
   B2     B5        -0.391404  6.159995
   B14     B18        -2.522816  6.612486
   B14     B19        -2.288776  6.526276
            .
            .
            .
   A1    B1       -0.257896    8.490313  ! These are the generic MJ non-bonded
   A1    B2       -0.117140    8.439674  ! parameters for all N_A*(N_B-1)/2
   A1    B3       -0.160122    4.500000  ! pairs of interactions.
   A1    B4       -0.257896    8.490313
            .
            .
            .
   A53    B47       -0.061876    6.957982
   A53    B48       -0.071323    6.762163
   A53    B49       -0.109110    4.500000
   A53    B50       -0.086438    7.992829
   A53    B51       -0.216330    6.743982
   A53    B52       -0.109110    4.500000
   A53    B53       -0.128476    5.530169

   END

   **********************
   model_kb-gopair.param:

   198 bynu  ! Note this format specifies atom pair selection by number (bynu)
   1 4 -2.36758 3.843863    ! This specfies atom 1 and atom 4 emin and rmin.
   1 6 -0.924458 4.937081   ! These are the same as the intrachain NBFIXes
   2 5 -0.391404 6.159995   ! in the model_kb-go.param file.
            .
            .
            .
   102 105 -2.36758 6.726333
   102 106 -1.775686 6.685909
   103 106 -0.68216 5.186839

   Note: there are 53 residues (Ca atoms in each chain, thus 106 total atoms.

   There is an alternative form for specifying the pairwise interactions in
   this file, where specific atom selection syntax is used.

   **********************************
   Alternative model_kb-gopair.param:

   198   ! Note that there is no bynu here so full selection syntax below.
   sele bynu 1 end sele bynu 4 end -2.36758 3.843863
   sele bynu 1 end sele bynu 6 end -0.924458 4.937081
   sele bynu 2 end sele bynu 5 end -0.391404 6.159995
            .
            .
            .
   sele bynu 102 end sele bynu 105 end -2.36758 6.726333
   sele bynu 102 end sele bynu 106 end -1.775686 6.685909
   sele bynu 103 end sele bynu 106 end -0.68216 5.186839


   *********************
   Example input stream:

   The following CHARMM command script provides an example of 
   setting up energy calculations using the Go Pair facility.

   * Test input script for model system 1arr
   * dimer represented by a mixed KB-Go/MJ interaction
   * model.
   *

   ! Read the general RTF and parameter files
   read rtf card name model_kb-go.top
   read param card name model_kb-go.param

   ! Set-up the PSF by reading sequence, generating and 
   ! reading coordinates for each chain.
   read sequ pdb name go_a.pdb
   generate proa autogenerate angle dihedral
   read coor pdb name go_a.pdb

   read sequ pdb name go_b.pdb
   generate prob autogenerate angle dihedral
   read coor pdb name go_b.pdb

   ! Turn on periodic boundary conditions using
   ! images for a cubic volume
   set boxsize = 120

   read image card
   * IMAGE FILE FOR CUBIC TRANSFORMATION
   * BOX SIZE IS @boxsize X @boxsize X @boxsize ANGSTROMS
   *
   SCALE  @boxsize @boxsize @boxsize
   IMAGE  X
   TRANS    1.0       0.0       0.0
   IMAGE  A
   TRANS   -1.0       0.0       0.0
   IMAGE  XY
   TRANS    1.0       1.0       0.0
            .
            .
            .
   IMAGE  XBC
   TRANS    1.0      -1.0      -1.0
   IMAGE  BC
   TRANS    0.0      -1.0      -1.0
   IMAGE  ABC
   TRANS   -1.0      -1.0      -1.0
   END

   IMAGE BYSEGID XCEN 0.0 YCEN 0.0 ZCEN 0.0

   ! We will use the ETEN model in this section
   eten on

   energy cutnb 25 ctonnb 25 ctofnb 25 cutim 25
   set e_kbgo = ?ener   ! Set energy variable for later comparison

   ! Now set-up GoPair model
   open unit 1 read form name model_kb-gopair.param
   gopair read unit 1 eten on

   update  ! doing update here ensures that exclusion list is built

   energy

   set e_gopair = ?ener

   calc diff = abs ( @e_kbgo - @e_gopair )  ! value of diff should be zero

   ! Now turn gopair off and calculate energy, should match @e_kbgo
   gopair off

   update  ! rebuild corrected exclusion list

   energy

   calc diff = abs ( @e_kbgo - ?ener ) ! Again, this should be 0.

   ! Turn GoPair back on and check
   gopair on

   update

   energy

   calc diff = abs ( @e_kbgo - ?ener ) ! Again, this should be 0.

   ! Now turn off gopair and switch to ETEN OFF

   gopair off

   eten off

   update

   energy

   set e_kbgo = ?ener

   gopair on eten off

   update

   energy

   set e_gopair = ?ener

   calc diff = abs ( @e_kbgo - @e_gopair )  ! value of diff should be zero

   stop

   From here one may move on to do dynamics or any of the other molecular 
   mechanics manipulations with the GoPair model.
