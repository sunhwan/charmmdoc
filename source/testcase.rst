.. py:module:: testcase

================
CHARMM Testcases
================

The CHARMM test cases are designed to test features of CHARMM and some of
error handling.  Though the test cases are not designed as a tutorial and
some used options and input parameters are not recommended, the test cases
are a valuable learning tool in setting up input files for CHARMM.


.. _testcase_overview:

Notes about the Testcase Suite
------------------------------

Testcases are reformed.  All testcases before version 22 are collected
in ~/cnnXm/test/c20test and new tests are written while we develop
CHARMM.  The new testcases are gathered in ~cnnXm/test/c22test,
c23test, ...  Note the following.

(1) In the new testcase suite, we use formatted I/O for
    topology/parameter files in order for all testcases to run
    independently each other.
(2) We make testcases self-contained whenever possible.  If external
    data files are required to run the test, they are in ~/cnnXm/test/data.
(3) CHARMM command parameter 0 and 9 are reserved to point
    directories for data and scratch files.  These are defined in
    the stream file ~/cnnXm/test/datadir.def.  The default
    definitions are
    
    ::
    
        set  0  data/
        set  9  scratch/
        
(4) We do not include benchmark suite in the distribution package
    any more.  Benchmarks produced on some Harvard local machines are
    available upon request.
(5) A developer must run all the testcases and confirm the benchmarks
    before he send in his developmental version.  If he develops a
    new function, he has to provide proper testcases.


.. _testcase_instruction:

How to Run Testcases
--------------------

In order to run testcases, you have to locate ~/cnnXm/test/datadir.def
and review the content.  It directs the directories for data and
scratch files needed for testcases.  Since the stream file is
stream'ed at the beginning, you may specify other job characteristics
in the file, e.g., bomb level (by the BOMLEV command), FASTer level, etc.

A C-shell script is provided in ~/cnnXm/test/test.com, which expects
at least one argument.

::

   test.com host_machine_type [ output_dir benchmark_dir target_executable ]

(1) host_machine_type can be convex, cray, ibmrs, sgi, stardent, sun, etc.
    and must be specified.
(2) output_dir is an optional argument which indicates the directory
    where you collect output files.  The default is output.
(3) benchmark_dir is an optional and points the directory that
    contains benchmark files.  If you specify 0 for this argument,
    then no comparison is done.
(4) target_executable is the executable you are testing.  If not specified
    the default executable, ~/cnnXm/exec/{host_machine_type}/charmm, is
    used.

For example, if your host is a Convex and you want run under the
default setting, your command is

::

   test.com convex

When test.com is finished, you may find the report file that
records differences between the test output and the target benchmark.
Normally the report file is named as {output_dir}.rpt.  The review
may indicate small numeric discrepancies due to different machine
precision and compiler options used.

For non-vector machines such as SGI, SUN, etc., ignore the warning
message ``"***** No VECTOR code compiled."``


.. _testcase_c20test:

Testcases from CHARMM20
-----------------------

~charmm/test/c20test contains 48 testcases inherited from CHARMM20
test set (the VAX version CHARMM20 and previous).  All external file
I/O are done in CARD format and binary I/O testing is performed by
stdtest in c22test.

::

   (1) brbtest.inp
       File : toph8.rtf/param3.prm
       Model: AMN-CBX-AMN-CBX
       Test : (1) IC PARAM, IC EDIT, IC SEED, IC BUILD
              (2) COOR COPY, COOR ORIENT MASS
              (3) HBOND
              (4) MINI CONJ NSTEP 100
              (5) MINI NRAPH NSTEP 40
       NOTE : I/O in card format
              BRBTEST.CRD inserted.

   (2) constest.inp
       File : toph9.rtf/param6.prm
       Model: TRP
       Test : Various CONStraint commands
       NOTE : BOMBLEV set to -2 to read param6.inp 
              Modified to test proper dihedral constraints by Arnaud Blondel
              (15-Feb-94)

   (3) control.inp
       File : None
       Model: None
       Test : CMDPAR and flow control

   (4) djstest.inp
       File : toph10.rtf/param7.prm, bpti.crd
       Model: PTI (58aa) and 4 waters
       Test : ABNER minimization and various electrostatic options
       Note : Was reduced to hepta peptide with one water model.
              Recover the original model since it talks 1.44 CPU min.
              on Convex C220.  The coordinates are read from bpti.crd.
              djstest.crd is no longer required.

   (5) dyntest1.inp
       File : toph8.rtf/param3.prm
       Model: TRP
       Test : DYNAmics commands with various CONStraints

   (6) dyntest2.inp
       File : toph8.rtf/param3.prm
       Model: TRP
       Test : DYNAmics LANGevin commands with various CONStraints

   (6a)dyntest3.inp (from R. Stote)
       File : toph8.rtf/param3.prm
       Model: TRP
       Test : TESTS A NUMBER OF DYNAMICS RESTART CALCULATIONS

   (7) enbtest.inp
       File : None
       Model: CO-NH
       Test : Various nonbonded interaction options using VIBRAN commands
              on selected energy terms 
              (a) Hydrogen bond energy, (b) van der Waals energy
              (c) Atom electrostatics,  (d) Group electrostatics
              (e) Extended electrostatics
       Note : FASTer OFF in order to use SKIPE commands

   (8) enertest.inp
       File : toph10.rtf/param8.prm, bpti.crd
       Model: PTI and 4 waters
       Test : ENERGY/GETE commands for selected energy terms (SKIPE)
              under various FAST/nonbond options.
              Forces are also examined.

   (9) genertest.inp
       File : toph8.rtf/param3.prm
       Model: GLY-PRO, PRO-GLY, PRO-PRO
       Test : some of the generation and patching routines

   (10) h2otst.inp
       File : toph8.rtf/param3.prm
       Model: two waters
       Test : HBONd and NBONd commands then MINI ABNR/NRAPH
              Runs a water dimer to convergence and a true minimum.

   (11) hbondtest.inp
       File : toph9.rtf/param6.prm
       Model: RNase beta sheet part 1 (PRO VAL ASN THR PHE VAL HSC)
              RNase beta sheet part 2 (SER ILE THR ASP CYS ARG GLU)
              and 5 waters
       Test : HBUIld and HBOND commands
       Note : ANAL command testing commented out

   (12) hbuildst2.inp
       File : toph9.rtf/param6.prm
       Model: RNase beta sheet part 1 (PRO VAL ASN THR PHE VAL HSC)
              RNase beta sheet part 2 (SER ILE THR ASP CYS ARG GLU)
              and 5 ST2 waters
       Test : HBUIld and HBOND commands

   (13) ictest.inp
       File : toph8.rtf/param4.prm
       Model: AMN-CBX
       Test : internal coordinate and coordinate manipulation commands

   (14) imbetash.inp
       File : toph9.rtf/param6.prm
       Model: ALA-ALA betasheet
       Test : build a beta sheet by IMAGe/IMPAtch commands
              ABNR 25 step minimization

   (15) imh2otest.inp
       File : toph8.rtf/param3.prm, cubic.img, wat125.crd
       Model: Box of 125 OH2 waters
       Test : 10 step dynamics in cubic periodic boundary condition
       Note : imh2otest.img renamed to cubic.img
              imh2otest.crd renamed to wat125.crd

   (16) imst2test.inp
       File : toph9.rtf/param6.prm
       Model: Box of 125 ST2 waters
       Test : 10 step dynamics in cubic periodic boundary condition
       Note : imst2test.img renamed to cubic.img
              imst2test.crd renamed to st2125.crd
              dynamics performed.

   (17) imtest.inp
       File : toph9.rtf/param6.prm
       Model: ALA 9-mer betasheet
       Test : Build C2 rotated image and 100 step ABNR minimization
       Note : imtest.img and imtest.crd are inserted

   (18) langtest1.inp
       File : None
       Model: 4 extended atom butane
       Test : 2500 step Langevin dynamics with FBETA 6.657235

   (19) langtest2.inp
       File : None
       Model: 4 extended atom butane
       Test : 100000 step Langevin dynamics with FBETA 100.0
    
   (20) lsqptest.inp
       File : toprna10r.rtf/pardna10.prm
       Model: GUA-CYT
       Test : COOR LSQP (least-squares-plane) commands

   (21) maatest.inp
       File : None
       Model: N-methyl alanyl acetamide
       Test : Dihedral constraint ABNR minimization to mimic Uray-Bradley
              terms.  Then, release the constraint and further ABNR-NRAPH
              minimize.  Perform VIBRAN.
       Note : maatest.crd inserted.

   (22) nbondtest.inp
       File : toph19.rtf/param19.prm, bpti.crd
       Model: PTI and 4 waters
       Test : FASTer ON/OFF nonbond interaction energy

   (23) noetest.inp
       File : toph19.rtf/param19.prm, bpti.crd
       Model: PTI and 4 waters
       Test : NOE distance restraints
       Note : bpti.crd is used instead of noetest.crd (identical)

   (24) partest.inp
       File : toprna10r.rtf/pardna10.prm, partest.prm
       Model: GUA-CYT dimer
       Test : parameter file I/O
       Note : partest.par renamed to partest.prm

   (25) patchtest.inp
       File : toph9.rtf/param6.prm
       Model: seven HSC residue segment
       Test : patch

   (26) powelltes.inp
       File : toph9.rtf/param6.prm, bpti.crd
       Model: PTI and 4 waters
       Test : MINI POWEll with SHAKE and HARMonic constraints
       Note : powelltes.man and powelltes.sol are replaced by bpti.crd
              TIP3 waters are used instead of ST2 waters.
              (ST2 minimization is found in rigidst2.inp)

   (27) psftest.inp
       File : toph9.rtf/param5.prm, rtftest.psf
       Model: PRO-PRO-PRO and ALA-ALA-ALA
       Test : PSF I/O

   (28) quasi.inp
       File : toph9.rtf/param5.prm
       Model: Ethanol
       Test : VIBRAN, quasiharmonic dynamics

   (29) rgyrtest.inp
       File : toprna10r.rtf/pardna10.prm
       Model: ADE-CYT-GUA-URI
       Test : COOR RGYR

   (30) rigidst2.inp
       File : toph9.rtf/param6.prm
       Model: two ST2 waters
       Test : HBUILD, TEST FIRST, Constraint/unconstraint minimization
              and dynamics

   (31) rtftest.inp
       File : param5.prm, rtftest.rtf, toprna10r.rtf
       Model: 21 amino acid sequence
       Test : RTF I/O

   (31a) rtf2.inp [Ryszard Czerminski, 30-Apr-92]
       File : all *rtf files from data directory
       Test : removing old residues when reading rtf append
       Note : does not work properly (yet)

   (32) sbdtest1.inp
       File : toph10.rtf/param7.prm, sbdtest1.pot
       Model: TIP3 water
       Test : SBOUNDARY and energy due to the boundary potential

   (33) sbdtest2.inp
       File : toph10.rtf/param7.prm, sbdtest2.pot
       Model: ST2 water
       Test : SBOUNDARY and energy due to the boundary potential

   (34) sbpgentst.inp
       File : sbpgentst.sbt
       Model: None
       Test : SBOUNDARY POTENTIAL

   (35) simp.inp
       File : toph9.rtf/param5.prm
       Model: AMN-CBX, AMN-CBX
       Test : COOR ORIENT, Q commands
       Note : simp.crd is inserted

   (36) st2test.inp
       File : toph9.rtf/param6.prm, st2125.crd
       Model: 125 ST2 waters
       Test : 15 step Verlet dynamics
       Note : st2test.crd renamed to st2125.crd

   (37) surftst.inp
       File : toph9.rtf/param6.prm, bpti.crd
       Model: PTI and 4 waters
       Test : COOR SURFace and COOR CONTact commands
              Checks the accessible surface calculation.
       Note : use full BPTI structure instead of a shortened one.
              surftst.chr is no longer needed.

   (38) test.inp
       File : toph8.rtf/param4.prm, bpti.crd
       Model: PTI and 4 waters
       Test : HBUILD, 10 step Verlet dynamics, INTERaction, MINI CONJ
       Note : use full BPTI structure instead of a shortened one.
              test.crd is no longer needed.

   (39) testcons.inp
       File : top9.rtf/param6.prm, lysozyme.crd
       Model: Lysozyme (129aa)
       Test : Harmonic atom constraints
       Note : testcons.src renamed to lysozyme.crd

   (40) testsel2.inp
       File : toph8.rtf/param4.prm
       Model: part of PTI sequence (60 atoms)
       Test : atom and tag selection, define command
       Note : test.crd is inserted

   (41) tipstest.inp
       File : tip125.crd, cubic.img
       Model: box of 125 TIP3 waters
       Test : energy with and without minimum image convention
       Note : topwat.inp and parwat.inp inserted
              tipstest1.crd replaced by tip125.crd
              tipstest1.img replaced by cubic.img
              tipstest1.inp renamed to tipstest.inp

   (42) trnphi.inp
       File : toph8.rtf/param3.prm
       Model: TRP
       Test : 500 step dynamics, MONITOR, MERGE, CORREL commands
       Note : trnphi.crd is inserted

   (43) vibpafl.inp
       File : None
       Model: All hydrogen methane
       Test : VIBRAN DIAG, READ, EDIT and PAFL commands

   (44) vibran.inp
       File : toph9.rtf/param6.prm
       Model: TRP
       Test : VIBRAN commands

   (45) vibrtst.inp
       File : toph8.rtf/param3.prm
       Model: AMN-CBX, AMN-CBX
       Test : VIBRAN commands
       Note : vibrtst.crd inserted

   (46) vibwat.inp
       File : None
       Model: a water
       Test : VIBRAN commands

   (47) voltest.inp
       File : toph10.rtf/param8.prm, bpti.crd
       Model: BPTI
       Test : SCALAR and COOR VOLUME commands
       Note : voltest.crd is replaced by bpti.crd

   (48) xray.inp
       File : toph10.rtf/param8.prm, bpti.crd
       Model: BPTI
       Test : WRITE XRAY command
       Note : enertest.crd is replaced by bpti.crd


.. _testcase_c22test:

New Testcases in CHARMM Version 22
----------------------------------

The following tests are written during the CHARMM22 development
period.  Note that most tests are self-contained and only lengthy
data files are left out in the data directory.  We also need to use
toph19.rtf/param19.prm.  stdtest test most CHARMM commands supported
in the version.

General Tests
^^^^^^^^^^^^^

::

   (0) stdtest.inp [Ryszard Czerminski, 20-Dec-91]
       File : toph19.rtf/param19.prm
       Model: ALA-TRP
       Test : most CHARMM commands

   (1) block1.inp [Bruce Tidor]
       File : topnah1r.rtf/parnah1r.prm, gal11.crd
       Model: GAL11 (dimer of ADE-ADE-GUA-THY-GUA-THY-GUA-ADE-CYT-ADE-THY)
       Test : patch and  block commands
       Note : used to be self-contained.  RTF, PARAM and CRD are separated
              out to topnah1r.rtf, parnah1r.prm and gal11.crd respectively.

   (2) block2.inp [Ryszard Czerminski, 11-Dec-91]
       File : None
       Model: Methanol (Me-OH to HO-Me mutation)
       Test : BLOCK FREE energy calculation
       Note : command line parameters 1-8 in use

   (3) cortst.inp
       File : None
       Model: ACE GLY GLY GLY GLY GLY GLY GLY GLY GLY GLY CBX
       Test : 100 step dynamics and CORRelation commands
       Note : modified from the VAX version cortst.inp.
              self-contained.  NOT working on some machines.

   (4) covaritst.inp [Charlie L. Brooks III, 09-Dec-91]
       File : None
       Model: deca-alanine
       Test : 50 step dynamics and cross correlation calculation
              COOR COVAriance command

   (5) ewions.inp [Roland Stote and Stephen Fleischman, 04-Dec-91]
       File : None
       Model: 108 Na(+)Cl(-) in a cubic box
       Test : Ewald summation energy under various FAST options.

   (6) exsgtst.inp [Ryszard Czerminski, 11-Dec-91] ! TO BE REMOVED !!!
       File : toph19.rtf/param19.prm, bpti.crd
       Model: BPRI
       Test : UPDATE/ENERGY EXSG subcommand

   (6a)exsg.inp [Ryszard Czerminski, 26-Mar-92]
       File : None
       Model: four hydrogen atoms
       Test : UPDATE/ENERGY EXSG subcommand

   (7) fshake1.inp [Stephen Fleischman, 04-Dec-91]
       File : cubic.img, tip125.crd
       Model: Glycerol in 125 water periodic box
       Test : SHAKE FAST against normal SHAKE

   (8) fshake2.inp [Stephen Fleischman, 04-Dec-91]
       File : cubic.img, tip125.crd
       Model: 125 water periodic box
       Test : SHAKE FAST against normal SHAKE

   (9) icfix.inp [Charlie L. Brooks III, 09-Dec-91]
       File : None
       Model: Three methane molecules
       Test : TSM ic constraint commands

   (10) icpert.inp [Charlie L. Brooks III, 09-Dec-91]
       File : tip125.crd, cubic.img
       Model: ACE-ALA-CBX in 125 TIP3 water periodic boundary box
       Test : the internal coordinate constraint and TSM commands

   (11) mewtest.inp [Charlie L. Brooks III, 09-Dec-91]
       File : mewtest.crd, cubic.img
       Model: Methane in 245 TIP3 water rectangular periodic boundary box
       Test : non-linear lambda scaling for methane -> nothing perturbation

   (12) slowgr.inp [Charlie L. Brooks III, 09-Dec-91]
       File : None
       Model: Ethanol -> Propane
       Test : TSM slow growth free energy simulation example

   (13) solanal.inp [Charlie L. Brooks III, 09-Dec-91]
       File : tip216.crd, cubic.img
       Model: 216 water molecules in a periodic box
       Test : solvent analysis on water

   (14) window.inp [Charlie L. Brooks III, 09-Dec-91]
       File : None
       Model: Ethanol -> Propane
       Test : TSM window/TI free energy simulation example

   (15) path.inp [Ryszard Czerminski, 11-Dec-91]
       File : None
       Model: Alanine Dipeptide
       Test : PATH between minima
       Note : NOT working, YW 17-Dec-91

   (16) pert.inp [Ryszard Czerminski, 11-Dec-91]
       File : None
       Model: Methanol (Me-OH to HO-Me mutation)
       Test : Free energy perturbation calculation by PERT command

   (17) travel.inp [Stefan Fischer, 20-Jun-91]
       File : chair.crd, boat.crd
       Model: Cyclohexane
       Test : TRAVEL commands

   (18) umbrella.inp [Jeyapandian Kottalam & Youngdo Won, 10-Dec-91]
       File : None
       Model: Cyclohexane
       Test : RXNCOR commands

   (19) xtlala1.inp [Martin J. Field, 22-Nov-90]
       File : None
       Model: Alanine crystal
       Test : COOR CONVERT and CRYSTAL commands.  Crystal optimization.
       Note : xtl_ala[1-4].inp are merged into this testcase

   (20) xtlala2.inp [Martin J. Field, 22-Nov-90]
       File : None
       Model: Alanine P1 crystal
       Test : CRYSTAL commands.  Crystal vibration and phonon analysis.
       Note : xtl_ala[5-6].inp are merged into this testcase

   (21) xtlala3.inp [Martin J. Field, 24-Jan-91]
       File : None
       Model: Alanine P1 crystal
       Test : two 100 step CPT dynamics
       Note : xtl_ala7.inp renamed to xtlala3.inp


Energy Tests
^^^^^^^^^^^^

::

   (1) cuttest1.inp [Stephen Fleischman, 04-Dec-91]
       File : toph19.rtf/param19.prm, bpti.crd
       Model: BPTI with four crystal waters
       Test : energy and force under various ATOM electrostatic cutoff
              options and VATOM VSWITCH with truncated cutoff (i.e., 
              CTONNB = CTOFNB).
           
   (2) cuttest2.inp [Stephen Fleischman, 04-Dec-91]
       File : toph19.rtf/param19.prm, bpti.crd
       Model: BPTI with four crystal waters
       Test : energy and force under various ATOM electrostatic cutoff
              options and VATOM VSHIFT.
           
   (3) cuttest3.inp [Stephen Fleischman, 04-Dec-91]
       File : toph19.rtf/param19.prm, bpti.crd
       Model: BPTI with four crystal waters
       Test : energy and force under various ATOM electrostatic cutoff
              options with van der Waals interactions skipped.
           
   (4) cuttest4.inp [Stephen Fleischman, 04-Dec-91]
       File : toph19.rtf/param19.prm, bpti.crd
       Model: BPTI with four crystal waters
       Test : energy and force under various ATOM electrostatic cutoff
              options and VATOM VSWITCH.
           
   (5) cuttest5.inp [Stephen Fleischman, 04-Dec-91]
       File : toph19.rtf/param19.prm, bpti.crd
       Model: BPTI with four crystal waters
       Test : energy and force under various VATOM van der Waals options
              with  electrostatic interactions skipped.
           
   (6) ew14test.inp [Stephen Fleischman, 04-Dec-91]
       File : tip125.crd, cubic.img, ew14test.str
       Model: Glycerol in 125 TIP3 water periodic box
       Test : Ewald energy calculation with exclusions

   (7) ewh2oderiv.inp [Stephen Fleischman, 04-Dec-91]
       File : tip125.crd, cubic.img, ewh2oderiv.str
       Model: 125 TIP3 water periodic box
       Test : Ewald derivative and energy with various van der Waals 
              cutoff options

   (8) ewh2oexcl.inp [Stephen Fleischman, 04-Dec-91]
       File : tip125.crd, cubic.img
       Model: 125 TIP3 water periodic box
       Test : Ewald nonbond exclusions


Dynamics Test
^^^^^^^^^^^^^

::

   (1) ewtipdyn.inp [Stephen Fleischman, 04-Dec-91]
       File : tip125.crd, cubic.img
       Model: 125 TIP3 water periodic box
       Test : DYNA VERLET 50 step with EWALD and SHAKE (FAST)
              FASTer VECTOR, SCALAR and OFF with either VSWITCH or VSHIFT

QUANTUM
^^^^^^^

::

   (1) quantum1.inp [Jeff Evansec, 28-May-92]
       File : none
       Model: monohydrated acetone
       Test : quantum mechanics / molecular mechanics


.. _testcase_c23test:

New Testcases in CHARMM Version 23
----------------------------------

The following tests are written during the CHARMM23 development period.
New features and major modifications are tested with the testcases.

::

   (1) clustst.inp [Mary E. Karpen, 09-Jan-93]
       File : clustst.hex
       Model: YPGDV peptide
       Test : CLUSter

   (2) cmdpar.inp [Leo Caves, 18-Jan-1994]
       File : None
       Model: None
       Test : command line parameters: assignment,substitution and manipulation.
       Note : error handling tested. one temporary file created to test import
              of parameter from external file.

   (3) mmfptest.inp [Benoit Roux, 31-Jan-94]
       File : toph19.inp, param19.inp
       Model: tripeptide ASP-ALA-ARG
       Test : miscelaneous boundary and restraints
       Note : command line parameters 1 in use

   (4) mtsm1.inp [Masa Watanabe, 18-Aug-1993]
       File : toph19.rtf and param19.prm
       Model: Met-enkephalin
       Test : Multiple time-step Method (MTS)

   (5) mtsm2.inp [Masa Watanabe, 18-Aug-1993]
       File : toph19.rtf and param19.prm
       Model: Met-enkephalin
       Test : MTS method with Nose-Hoover heat bath

   (6) nmrtest1.inp [Benoit Roux, 31-Jan-94]
       File : toph19.inp, param19.inp
       Model: tripeptide ASP-ALA-ARG
       Test : generate trajectory and calculate NMR properties
       Note : command line parameters 1 in use

   (7) nose1.inp [Masa Watanabe, 18-Aug-1993]
       File : toph8.rtf and param3.prm
       Model: TIP3P water in a box
       Test : Single Nose-Hoover Dynamic Method

   (8) nose2.inp [Masa Watanabe, 18-Aug-1993]
       File : None
       Model: A ethane molecule in TIP3P water
       Test : Nose-Hoover Method with multiple heat bath

   (9) replica.inp [Leo Caves, 18-Aug-1993]
       File : toph19.inp, param19.inp
       Model: alanine dipeptide
       Test : replication of PSF; energy,forces and nonbonded exclusions.
       Note : nonbonded exclusions cannot be tested for all list generation 
              routines on a given machine (eg. CONVEX specific FNBL). 

   (10) rism.inp [Georgios Archontis, 18-Aug-1993]
       File : None
       Model: pure solvent: tip3p water 
              solute 1: extended-carbon with weak charge
              solute 2: diatomic
       Test:  solvent-solvent calculation
              solute-solvent calculation
              solute-solute calculation
              calculation of chemical potential of solvation for
              the two solutes and  decomposition to 
              energy and entropy of solvation

   (11) zmat.inp [Benoit Roux, 31-Jan-94]
       File : toph19.inp, param19.inp
       Model: TIP3P water dimer
       Test : construct the optimized configuration for water dimer
       Note : difficult to do with the IC table


.. _testcase_c24test:

New Testcases in CHARMM Version 24
----------------------------------

The following tests are written during the CHARMM24 development period.
New features and major modifications are tested with the testcases.

::

   (1) autogen.inp [Rick Venable, 04-Aug-1995]
       File : none
       Model: alanine tetrapeptide
       Test : automatic regeneration of angles and dihedrals
              designed for use after multiple PATCh statements
              compare PSF before and after; should be identical
       Note : patches may now be written w/o ANGL and DIHE terms
              the PSF should not contain any water molecules


   (2) bcdtest.inp [Wonpil Im, 02-Aug-95]
       File : None
       Model: Beta-cyclodextrin with 8 crystal waters
       Test : Crystal build
       Note : introduced to check the unit cell rotation bugfix

   (3) block3.inp [Stefan Boresch, 01-Aug-95]
       File : datadir.def; data/tip125.crd; data/cubic.img;
              scratch files (trajectories+restart files) produced.
       Model: ethane/methanol hybrid in small TIP3 water box
       Test : Tests BLOCK in combination with IMAGE module
       Note : Runs a few steps of dynamics and tests the supported
              post-processing options

   (4) calc.inp [Benoit Roux, 15-Feb-94]

   (5) dihtest1.inp [Arnaud Blondel, 15-Feb-94]
       File : top_all22_na.inp / par_all22_na.inp
       Model: ADE : CYT
       Test : Various dihedral energy routines.
       NOTE : Introduced to check the correspondance between old
              and new dihedral energy routines.
              Designed to test various parts of the code.

   (6) dihtest2.inp [Arnaud Blondel, 15-Feb-94]
       File : None
       Model: Extended atom butane
       Test : None planar equilibrium dihedral energy terms.
              TEST SECOnd command (second derivatives).
       NOTE : New dihedral energy routine only.

   (7) dimb1.inp [Herman van Vlijmen, 15-Feb-95]
       File : None
       Model: Deca-alanine
       Test : DIMB and DIMB DWIN normal mode calculations

   (8) dimb2.inp [Herman van Vlijmen, 15-Feb-95]
       File : None
       Model: Deca-alanine
       Test : Reduced basis diagonalization with compressed Hessian

   (9) dyn4Dtest.inp [Carol B. Post, 15-Feb-95]

   (10) mmfptest2.inp [Benoit Roux, 15-Feb-94]

   (11) mmfptest3.inp [Benoit Roux, 15-Feb-94]

   (12) nptdyn.inp [Scott Feller, 01-Aug-1995]
       File : tip125.crd
       Model: cubic water box
       Test : extended pressure system; CPT dynamics
             isotropic system with Langevin piston
       Note : minimal test, 40 steps

   (13) pbound1.inp [Charles L. Brooks, III and William A. Shirley, 15-Aug-95]
       File : toph19.rtf, param19.prm, cubic.img, tip125.crd
       Model: 125 TIP3P water molecules
       Test : The explicit (simple) periodic boundaries (BOUND command)
       Note : This test case checks the energy and dynamics of a simple
              implementation of explicit periodic boundaries against the
              IMAGE facility values.  As of c24b1, works only with the
              FAST Scalar routine.

   (14) pert2.inp [Stefan Boresch, 01-Aug-95]
       File : none; on none Unix platforms /dev/null may have to be
              changed appropriately.
       Model: one bond-length of triatomic symmetric molecule is
              changed from 1 to 2 A (gas phase).  Various protocols
              are used
       Test : Tests PERT in combination with SHAKE; constraint correction
       Note : The free energy differences for this system can be
              calculated analytically; these results are listed at
              the end of the input file together with the results
              I obtained on one of our HP's.  The differences between
              simulation and analytical result should give some feeling
              as to what to expect on other platforms.

   (15) tntest1.inp [P. Derreumaux 28-Jan-94]
       File : None
       Model: Alanine dipeptide
       Test : Truncated Newton Minimizer TNPACK


.. _testcase_c25test:

New Testcases in CHARMM Version 25
----------------------------------

The following tests are written during the CHARMM25 development period.
New features and major modifications are tested with the testcases.

::

   (1) allxtl.inp [Bernard R. Brooks, 15-JUL-97, c25b1]
       Test : Crystal symmetry test case.  Test ALL of the crystal types.

   (2) anal.inp [Bernard R. Brooks, 15-JUL-97, c25b1]
       Test : long/short energy print-out with the anal command

   (3) cortst25.inp [Bernard R. Brooks, 15-JUL-97, c25b1]
       Model: alfa decaglycine

   (4) cwat.inp [Paul Lyne, 01-Sep-1995]
       File : libfil.dat, modpot.dat and cwat.str
       Model: two water molecules
       Test : QM(CADPAC)/MM(CHARMM)
              First water is calculated with TIP3P and the second water is
              QM-STO3G.
       Note : The Hamiltonian and other CADPAC control commands are found in
              cwat.str


   (5) ewald_atom.inp [Bernard R. Brooks, 15-JUL-97, c25b1]

   (6) ewald_grp.inp [Bernard R. Brooks, 15-JUL-97, c25b1]

   (7) ewald_pert.inp [Bernard R. Brooks, 15-JUL-97, c25b1]

   (8) fastest.inp [Bernard R. Brooks, 15-JUL-97, c25b1]
       Test : fast options with various non-bond schemes
       Model: BPTI coordinates with all but 4 water removed

   (9) gmstst.inp [Milan Hodoscek]
       File : none
       Model: Alalnine 
       Test : QM(GAMESS)/MM(CHARMM)
              CTERM is QM and the rest is MM.  Link atom is between CA and C
       Note : Runs ~ 3 min on HP-735

   (10) hba1.inp [Lennart Nilsson]
       File : top_all22_prot.inp and par_all22_prot.inp
       Model: three water molecules
       Test : hydrogen bond analysis facility

   (11) helix.inp [Bernard R. Brooks, 15-JUL-97, c25b1]
       Test : helix analysis code
              makes a duplicate of the perfect octamer of CA atoms and
              performs coordinates transformations to align it as a helix
              in a parallel fashion.

   (12) hrbestfit.inp [Bernard R. Brooks, 15-JUL-97, c25b1]
       Test : Harmonic restraints best-fit test case
       Model: acetamide
       Files: toph19.rtf, param19.prm

   (13) mtsm3.inp [Masa Watanabe, 09-Feb-1996]
       Model: Met-enkephalin
       Test : multiple time scaled method

   (14) pull.inp [Lennart Nilsson]
       Test : application of external forces to the system

   (15) quiet.inp [Bernard R. Brooks, 15-JUL-97, c25b1]
       Test : a number of dynamics calculations

   (16) resdtest.inp [Bernard R. Brooks, 15-JUL-97, c25b1]
       Test : restrained distances

   (17) rpath1.inp [Bernard R. Brooks, 15-JUL-97, c25b1]
       Test : Replica Path method

   (18) td149.inp [Bernard R. Brooks, 15-JUL-97 c25b1]
       Model: truncated dodecahedron box of 149 water molecules

   (19) vibwat25.inp [Bernard R. Brooks, 15-JUL-97, c25b1]
       Test : water normal modes

   (20) xtloct1.inp
       Model: octane crystal
       Test : cryst building code using X-cryst fractional coordinates
              as starting structure
       Ref. : H. Mathisen and N. Norman, Acta Chemica Scandinavica (1961) 15, 1747


.. _testcase_c26test:

New Testcases in CHARMM Version 26
----------------------------------

The following tests are written during the CHARMM26 development period.
New features and major modifications are tested with the testcases.

::

   (1) block4.inp [Thomas Simonson, 24-JUL-97, c26a1]
       Test : BLOCK enhancement that allows different coefficients for different
              energy terms in free energy simulations
       Model: hybrid Asn-Asp in vacuum


   (2) cftigas.inp [Krzysztof Kuczera, 22-Mar-1997, c26a1]
       Test : one-dimensional conformational thermodynamic integration
       Model: butane in vacuum
       Files: top_all22_prot.inp and par_all22_prot.inp

   (3) cftmgas.inp [Krzysztof Kuczera, 22-Mar-1997, c26a1]
       Test : CFTM protocol:
              (a) Perform MD with holonomic constraints, save gradient file
              (b) Perform elementary analysis
       Model: ALA10 with all phi and psi constrained, vacuum
       Files: top_all22_prot.inp and par_all22_prot.inp

   (4) conmin.inp [Krzysztof Kuczera, 22-Mar-1997, c26a1]
       Test : energy optimization with holonomic constraints
       Model: ALA10 with all phi and psi fixed
       Files: top_all22_prot.inp and par_all22_prot.inp

   (5) luptst.inp [Krzysztof Kuczera, 22-Mar-1997, c26a1]
       Test : Generate butane and run LUP protocol to create path between
              trans and gauche minima
       Files: top_all22_prot.inp and par_all22_prot.inp

   (6) mbtest18.inp [Robert Nagle, 10-JUL-97, c26a1]
       Test : MBOND/CHARMM TEST #18 
              Test MBOND dynamics body-based production runs with different
              conditions and different substructuring.
              Begins from an body based equilibration.

   (7) mbtest19.inp [Robert Nagle, 10-JUL-97, c26a1]
       Test : MBOND/CHARMM TEST #19
              Test MBOND dynamics body-based production runs with thermostat on.

   (8) mbtest24.inp [Robert Nagle, 10-JUL-97, c26a1]
       Test : MBOND/CHARMM TEST #24
              mode generation/storing

   (9) mbtest25.inp [Robert Nagle, 10-JUL-97, c26a1]
       Test : MBOND/CHARMM TEST #25 MTS

   (10) pathint.inp [Benoit Roux, 10-JUL-97, c26a1]
       Test : Classical trajectories for acetylacetone close to
              the transition state
       Files: top_all22_prot.inp and par_all22_prot.inp

   (11) pbeqtest1.inp [Benoit Roux, 10-JUL-97, c26a1]
       Test : (1) the Poisson Boltzmann Equation solver for alanine dipeptide
              (2) the Atomic Born Radii
              (3) the Solvation Forces
       Files: top_all22_prot.inp, par_all22_prot.inp and radius.str

   (12) pbeqtest2.inp [Benoit Roux, 10-JUL-97, c26a1]
       Test : (1) the Poisson Boltzmann Equation solver with membrane
              (2) the Solvation Forces with membrane
       Files: top_all22_prot.inp, par_all22_prot.inp and radius.str

   (13) pm6test1.inp [Benoit Roux]


   (14) whamtest.inp [Benoit Roux]
       Test : Perturbation calculation with WHAM post-processing


.. _testcase_nbondtest:

NonBOND Testcases in CHARMM Version 25
--------------------------------------

The following tests are written during the CHARMM25 development period.
New features and major modifications of nonbonding energy routines are tested
with these testcases.

::

   (1) coul_test.inp [Bernard R. Brooks, 15-JUL-97, c25b1]
       Test : calculate coulomb interaction of 2 protons

   (2) form.inp [Bernard R. Brooks, 14-Jul-1997, c25b1]
       Test : calculate group-group energy (nonbonded) of 2 formamides

   (3) form_ewald.inp [Bernard R. Brooks, 14-Jul-1997]

   (4) form_ewald_m.inp [Bernard R. Brooks, 14-Jul-1997]

   (5) form_fsw.inp [Bernard R. Brooks, 14-Jul-1997]

   (6) form_mm_m.inp [Bernard R. Brooks, 14-Jul-1997]

   (7) form_simp_m.inp [Bernard R. Brooks, 14-Jul-1997]

   (8) ion_fsw.inp [Bernard R. Brooks, 14-Jul-1997]

   (9) vdw_test.inp [Bernard R. Brooks, 14-Jul-1997]


.. _testcase_mmfftest:

MMFF Testcases in CHARMM Version 25
-----------------------------------

The following tests CHARMM/MMFF features.  Common external files are
in the test/data directory:

::

        mmffang.par
        mmffbond.par
        mmffchg.par
        mmffdef.par
        mmffoop.par
        mmffstbn.par
        mmffsup.par
        mmffsymb.par
        mmfftor.par
        mmffvdw.par
        mmff_setup.str

::

   (1) mmff.inp [Ryszard Czerminski, 11-May-1993]
       Test : MMFF parameter reader, energy & derivatives

   (2) mmff_amino.inp [Ryszard Czerminski, 11-May-1993]
       Test : new MMFF ring perception code

   (3) mmff_append.inp [Ryszard Czerminski, 11-May-1993]
       Test : 'read merck ... append' command

   (4) mmff_c60.inp [Tom Halgren, 11-May-1993]
       Test : PARAMETERS of C60 releated molecular moieties

   (5) mmff_clpert.inp [Ryszard Czerminski, 11-May-1993]
       Model: chloromethane
       Test : PERT command (slow growth method) to calculate free energy
              perturbation for migrating -Cl atom in chloromethane
              (CH3-Cl -> Cl-CH3)

   (6) mmff_cutoff.inp [Jay L. Banks, 02-Dec-1993]
       Test : MMFF cutoff schemes on small molecule 

   (7) mmff_gener.inp [Ryszard Czerminski, 11-May-1993]
       Test : MMFF parameter reader, energy & derivatives

   (8) mmff_h2o.inp [Ryszard Czerminski, 11-May-1993]
       Test : MMFF parameter reader

   (9) mmff_icpert.inp [Jay L. Banks, 13-Apr-94]
       Files: icala.mrk
       Test : internal coordinate TSM with MMFF force field

   (10) mmff_pep27.inp [Thomas, 11-Feb-1995]
       File : top_all22_prot_mmff.inp

   (11) mmff_pert.inp [Ryszard Czerminski, 08-Sep-1994]
       Files: sg.15K.punit lambda schedule file
       Test : PERT command (slow growth method) with MMFF to calculate free
              energy perturbation for migrating -OH group in methanol
              (CH3-OH -> OH-CH3)

   (12) mmff_ring.inp [Ryszard Czerminski, 11-May-1993]
       Test : new MMFF ring perception code

   (13) mmff_rtf.inp [Ryszard Czerminski, 11-May-1993]
       Test : new RTF keywords (SINGLE, DOUBLE & TRIPLE) with ala2


   (14) mmff_solanal.inp [Jay L. Banks, 13-Oct-1993]
       Test : MMFF water model using solvent analysis

   (15) mmff_vib.inp [Ryszard Czerminski, 30-Sep-1993]
       Test : VIBRAN facility working with MMFF


.. _testcase_graftest:

New Graphics Testcases in CHARMM Version 24
-------------------------------------------

As of c24b1, CHARMM graphics is substantially enhanced and a testcase
suite is also developed (test/cgrftest).  Not all testcases in
cgrftest may be valid on a given system; grftest.com will select
and run the appropriate testcases based on the graphics keyword found
in pref.dat.

::

             [ XDISPLAY GLDISPLAY NODISPLAY NOGRAPHICS APOLLO ]

   grfapo.inp     [Rick Venable, 05-Aug-1995]
   grfgldsp.inp
   grfnodsp.inp
   grfnowin.inp
   grfxwin.inp
       Files: toph19.rtf, param19.prm, bpti.crd
       Model: BPTI crystal structure
       Test : exercise main graphics features
              separate tests for display window and derived files
       Note : for grfxwin, the env var DISPLAY must be set
              PostScript, other output files are created in cgrftest
              see following table for testcase/keyword matchup

       --------------------------------------------------------------
         Testcase                    Keyword in pref.dat
       --------------------------------------------------------------
                          XDISPLAY   GLDISPLAY  NODISPLAY    APOLLO
         grfapo.inp                                             +
         grfgldsp.inp                    +
         grfnodsp.inp                    +          +           +
         grfnowin.inp        +
         grfxwin.inp         +
       --------------------------------------------------------------

