.. py::module: mmff

====================================
Merck Molecular Force Field (MMFF94)
====================================

.. _mmff_usage:

Usage
=====

In order to use MMFF in CHARMM, the user has to issue the following
commands:

::

   1. use mmff force field
   2. <read mmff parameter files>
   3. (a) read rtf name <MMFF-capable rtf file>, or
      (b) read merck name <file_name>
      (c) read mol2 name <file_name>
      (d) read db mol_name name <file_name>
   4. read sequence  ! if input is via the rtf route (step 3 (a))
   5. generate  ! note that there may be multiple segments in one .mrk file
   6. patch     ! if input is via rtf/sequence route, apply appropriate patches
                ! to force a new mmff_setup; either include the keyword "mmff" 
                ! on the final patch or follow the final patch by the command:
                ! "use mmff atom types"
   7. read coord, or ic build  ! if input is via the read rtf/sequence route.  

Steps 1 & 2 can be done by streaming the file "mmff_setup.STR."  An example
of this file is shown below.  Documentation on the contents and usage of the
MMFF parameter files may be found in mmff_params.doc.

Step 3a requires a MMFF-capable rtf file.  This means a file in which 
BOND records have been replaced by analogous DOUBLE or TRIPLE records for
cases in which the chemical structure (or any valid Kekule representation)
has a double or triple bond.  Mass records in a MMFF-capable rtf file 
must also be augmented to add the atomic symbol for each CHARMM atom type
after the atomic mass entry.  Note that MMFF-capable rtf files are *back 
compatible*.  That is, such rtf files can equally well be used for 
calculations that utilize the CHARMM force field.  Thus, it is *not* necessary 
to maintain two versions of the rtf files.

Format of .mrk file optionally read in step 3b
----------------------------------------------

Merck-format files consist of one or more consecutive molecular-data entries.

.. note::
   when embedded in a CHARMM input script, a mrk file must be followed by
   a card reading "END" in columns 1:3.

::

   FILE_FORMAT:

   	All entries in a Merck-format (.mrk) file have the format:
 
                Line #     # of Lines      Use
 
                  1            1            Header_1
                  2            1            Header_2
                  3            1            Number of atoms and bonds
                  4            n            Data on n atoms
                 4+n           k            Bonding data on the 5*k bonds 
                                            in the structure  (Each line 
                                            contains data on five bonds)
 
   Header_1:

   	The format for the first header line is:

   		(A70,A10)

   	Each field contains the following information:

                   column  Description of use
 
                    1-70    User defined title
                   71-72    Present Year (YY)
                   73-75    Present Date (DDD)
                   76-79    Time of Day (HHMM), e.g., "1709" for 5:09 pm
                      80    Must be a "1" for the file to be valid
 
   Header_2:

   	The second header line has the following format:

   		(A4,A8,X,A1,X,A65)
	
   	Each of the fields has the following information:

                   column  Description of use
 
                    1- 4    The string "MOL "
                    5-12    User name
                      14    Source of file : (e.g., E for MOLEDIT, C 
   			 for Cambridge, D for Distance Geometry etc.)
                   16-80    Column used by other programs such as the 
   			 Cambridge Programs and OPTIMOL
 
   Number_of_atoms_and_bonds:

   	The format for this record is:

   		(I5,X,I5)

   	Each of the fields has the following information:

                   column  Description of use
 
                    1-5     NATOM
                    7-11    NBND
 
   Data_on_atom_n:

   	The format for the atom records is:

           (3(F10.4,1X),I5,1X,I2,1X,I1,1X,I5,1X,3A4,F8.4,6X,A4)

   	Each of the fields has the following information:
 
           Columns     Field               Description
 
             1-10       X                  X coordinate of the atom
            12-21       Y                  Y coordinate of the atom
            23-32       Z                  Z coordinate of the atom
 
            34-38       Atomic Number      (I5) field containing the type
                                            of atom. (i.e. -- 6 for Carbon;
                                            8 for Oxygen; etc...) A value 
                                            of 0 indicates a lone pair.
 
            40-41        Atom Subtype       (I2) field: on output, contains the
                                            MMFF atom type; is not read on input
 
               43        Charge Code        Formal charge code of the atom.
 
            45-49       Sequence Number     (I5) field containing the unique
                                            number by which every atom in 
                                            the structure can be identified.
                                            Note: in the CHARMM implementation,
                                            these quantities are not actually 
                                            read.  However, the atoms are
                                            expected to be numbered consecutively
                                            from 1 to NATOM and to correspond to 
                                            the numbers used in the bond_data 
                                            records defined below.
 
            51-54        Atom Name          Left justified (A4) field. 
                                            Should be unique inside a
                                            given residue. (Examples -- "C24 ",
                                            "NH  ", etc...).
 
            55-58        Residue Name       Right justified (A4) field.
                                            (Examples -- " 123", "123A",
                                            etc...).
 
            59-62        Residue Type       Left justified (A4) field.
                                            (Examples -- "TRP ", "LYS ",
                                            etc...). 
 
            63-70        Partial Charge     (F8.4) field containing the partial
                                            charge of an atom in proton units.
                                            Note: this entry is written on output,
                                            but is not read on input.

            77-80        Segment ID         Left justified (A4) field containing 
                                            a one to four character segment ID
                                            identifier.
 
   Note: if any of the A4 fields specified above are blank, the file reader will
   construct a default name.
 
   Charge_code:
 
           The valid charge codes are:
 
                   Code            Charge Code
 
                     0              Neutral
                     1               +1
                     2               -1
                     3              Radical
                     4               +2
                     5               -2
                     6               +3
                     7               -3
                     8               +4
                     9               -4

   Bond_data:

   	The block of data at the end of the .mrk file contains the bonding 
   	information.  Each line of bond data can contain a maximum of five 
   	bond definitions.  The format for the bond data is:
		
   		5(I5,X,I5,X,I2,2X)
		
           For each bond definition,
 
   		Field       Description
 
                   IFROM       (I5) Sequence number of the starting 
                               atom of the bond 
 
   		ITO	    (I5) Sequence number of the terminating
                               atom of the bond 
 
   		ITYPE       (I2) Order of the bond. (i.e. 1 for a single 
   			    bond, 2 for a double bond, etc.)
                               Bond orders are always integral


end of mrk format specification
-------------------------------

As noted, the .mrk file reader in CHARMM can read concatenated .mrk files. 
It should also be possible to 'read merck ... append'.
These two input routes should be equivalent as far as final the data 
structure is concerned.

.. note::

   1. no binary parameter files are supported for MMFF.
   2. MMFF is an all hydrogen force field -- i.e., extended atoms
      are not supported

Format of .mol2 file optionally read in step 3c
-----------------------------------------------

SYBYL MOL2-format files provides a complete representation of a molecule for
use with software from Tripos Inc. (including SYBYL). Details of the format
can be found in documentation from Tripos Inc.
Note: when embedded in a CHARMM input script, a mol2 file must be followed by
a card reading "END" in columns 1:3.

::

   FILE_FORMAT:

   The exact content of MOL2 files generated by SYBYL may vary based on
   different processing of the molecules. However, it should at least contain
   the following records:

      @<TRIPOS>MOLECULE
      @<TRIPOS>ATOM
      @<TRIPOS>BOND
      @<TRIPOS>SUBSTRUCTURE

   These four sections provide different information about the molecule
   and are necessary to reconstruct the molecule.

   @<TRIPOS>MOLECULE section

   Format:

       mol_name
       num_atoms num_bonds num_subst num_feat num_sets
       mol_type
       charge_type

       mol_name:
       This entry indicates the name of the molecule and has a string format.

       num_atoms:
       This indicates the number of atoms in the molecule. Integer format.

       num_bonds:
       This indicates the number of bonds in the molecule. Integer format.

       num_subst:
       This indicates the number of substructures in the molecule. Integer format.

       num_feat:
       This indicates the number of features in the molecule. Integer format.

       num_sets:
       This indicates the number of sets in the molecule. Integer format.

       mol_type:
       This indicates the molecule type.

       charge_type:
       This indicates the type of charges associated with the molecule.

   @<TRIPOS>ATOM section

       The format of this section contains the following information

       (atom_id atom_name x y z atom_type subst_id subst_name charge)

       and has the following format:

       (I8,A4,4X,3(F10.4),1X,A4,3X,I4,1X,A4,6X,F8.4)

       Each of the fields has the following information:

            column  Field       Description of use

             1- 8   atom_id     the ID number of the atom at the time the mol2
                                file was created
             9-12   atom_name   the name of the atom
            17-26   x           the x coordinate of the atom
            27-36   y           the y coordinate of the atom
            37-46   z           the z coordinate of the atom
            48-51   atom_type   the SYBYL atom type for the atom
            55-58   subst_id    the ID number of the substructure containing
                                the atom
            60-63   subst_name  the name of the substructure containing the atom
            70-77   charge      the charge associated with the atom

   @<TRIPOS>BOND section

       The format of this section contains the following information

       (bond_id origin_atom_id target_atom_id bond_type)

       and has the following format:

       (1X,3I5,1X,2A)

       Each of the fields has the following information:

            column  Field         Description of use

             2- 6  bond_id        the ID number of the bond at the time the mol2
                                  file was created
             7-11  origin_atom_id the ID number of the atom at one end of the bond
            12-16  target_atom_id the ID number of the atom at the other end
                                  of the bond
            18-19  bond_type      the SYBYL bond type

   @<TRIPOS>SUBSTRUCTURE section

       The data line contains the substructure ID, name, root atom of the
       substructure, substructure type, dictionary type, chain type, subtype,
       number of inter substructure bonds, SYBYL status bits, and user defined
       comment. Information contained in this section is not read nor used by
       the MMFF module. The format is open for this section.


Format of .mol2 file optionally read in step 3d
-----------------------------------------------

SYBYL MOL2 database files have a format identical to that described in
step 3c. If the database is read in as an external file, there is no need
to put "END" at the end of every mol2 molecule.

end of mol2 format specification
--------------------------------

(1) Each atom in the MOL2 file should have a unique atom name in order
    for the MMFF bond types to be assigned properly.
(2) For external database reading capability, the maximum length of
    a molecule name in the MOL2 database file is currently set to be
    a string of 20 UPPERCASE characters. A molecule name is read in
    the line of mol_name in @<TRIPOS>MOLECULE section.
(3) Due to the fact that bonds are not explicitly typed in the MOL2
    format, a conversion of MOL2 non-integer bond type (e.g. ar and am)
    into MMFF recognizable type was made.  The type of an amide bond
    is always set to be 2. For aromatic bonds within an aromatic ring,
    they are assigned to be alternating single and double bonds.
    The algorithm first separates aromatic bonds (and the associated
    atoms) from any integer-type bond. It arbitrarily sets the first
    aromatic bond to be a single bond and then starts a loop of
    aromatic bond assignment. During the course of assignment,
    the surrounding connectivity information of an atom with aromatic
    bond type is taken into account.  However, problems may still occur
    during this step. The authors welcome reports of any problematic
    molecules.

Examples of MMFF usage in CHARMM are given in mmff*.inp files in the test 
directory.

.. _mmff_quanta:

Here is the current "prescription" for to use MMFF in CHARMm from
QUANTA.

(1) In the CHARMm menu, select "MMFF" within the "CHARMm MODE" menu item.

(2) Proceed as you normally would; until an alternative MODE is selected,  
    all requests for CHARMm energy services will use the MMFF force field.

.. note::

   QUANTA communicates with CHARMm by writing a .mrk (Merck-format) file
   named .charmm_mmff.  Because MMFF does not recognize special "aromatic" or
   "resonant" bond orders (e.g., 7), a translation to a 'Kekule' structure is
   made as the .mrk file is being written.  On some ocassions, the routines in
   QUANTA that make this translation (as of February 1996) do so incorrectly.  
   It is therefore safest - and sometimes *necessary* - for the QUANTA user to
   first employ the Molecular Editor to change the structure to Kekule format,
   to examine it visually, and to repair incorrect bonding if needed.

   Quanta also sends the requisite "stream mmff_setup.STR" and "read Merck"
   commands to CHARMm.  A typical mmff_setup.STR file is shown below:

   ::

      mmff_setup.STR
      --------------
      * setup of MMFF in CHARMM
      *
      use mmff force field

      open read form unit 1 name "$CHM_DATA/MMFFSUP.PAR"
      read parameter card mmff SUPP unit 1
      read parameter card mmff PROP name "$CHM_DATA/MMFFPROP.PAR"
      read parameter card mmff SYMB name "$CHM_DATA/MMFFSYMB.PAR"
      read parameter card mmff DEFI name "$CHM_DATA/MMFFDEF.PAR"
      read parameter card mmff BNDK name "$CHM_DATA/MMFFBNDK.PAR"
      read parameter card mmff HDEF name "$CHM_DATA/MMFFHDEF.PAR"
      read parameter card mmff AROM name "$CHM_DATA/MMFFAROM.PAR"
      read parameter card mmff VDW  name "$CHM_DATA/MMFFVDW.PAR"
      read parameter card mmff BOND name "$CHM_DATA/MMFFBOND.PAR"
      read parameter card mmff CHRG name "$CHM_DATA/MMFFCHG.PAR"
      read parameter card mmff PBCI name "$CHM_DATA/MMFFPBCI.PAR"
      read parameter card mmff ANGL name "$CHM_DATA/MMFFANG.PAR"
      read parameter card mmff STBN name "$CHM_DATA/MMFFSTBN.PAR"
      read parameter card mmff DFSB name "$CHM_DATA/MMFFDFSB.PAR"
      read parameter card mmff OOPL name "$CHM_DATA/MMFFOOP.PAR"
      read parameter card mmff TORS name "$CHM_DATA/MMFFTOR.PAR"
      close unit 1

      return

.. _mmff_status:

Status of MMFF implementation into CHARMM (February 1996)
=============================================================

This implementation of MMFF in CHARMM is principally due to Ryszard
Czerminski (MSI) and Jay Banks (first of MSI, later a consultant
to Merck and to NIH), working in conjunction with Tom Halgren (Merck).

Features currently supported in CHARMM/MMFF


(1) energy, first & second derivatives
(2) minimization
(3) dynamics
(4) most ATOM based cutoff options (force switch is not implemented for
    vdW interactions; for vdW force shift, a generalized version is used 
    with beta=4 -- see Steinbach and Brooks, J. Comput. Chem., 15, 667-683 
    (1994)). 
(5) fast routines, implelented using the "PARVEC" paradigm
(6) the multiple time step algorithm (should work, if it does not use custom 
    calls for energy services)
(7) PERT, BLOCK, and TSM free energy methods, but only for a limited range 
    of problems.  The current MMFF setup code requires that the input 
    structure be a valid chemical species (e.g., no more than four bonds 
    to carbon), and therefore does not allow for dummy atoms.  However, 
    it should be possible to use TSM for internal-coordinate perturbations 
    and BLOCK for perturbations in which the blocks are not interbonded 
    (examples are given in the mmff*pert*.inp scripts that may be found in 
    the test directory).  For PERT, it is also possible to use rtf/sequence 
    input and to add dummy atom(s) after the "generate" command has done a 
    MMFF  setup on the original data structure.  This would be accomplished 
    by applying one or more patches and then, without repeating the MMFF 
    setup (e.g., without again giving the generate command), using scalar 
    commands to set the MMFF atom types and partial charges.  See the 
    mmff_pert.inp script that may be found in the test directory (if it is 
    up to date). In this case, parameters for the dummy atom(s) are read 
    from the MMFFSUP.PAR supplementary-parameters file. An example of such 
    a file is shown below:
    
    ::

      -------------------------- MMFFSUP.PAR ------------------------------------
          1    1    0    0    1    0    0    2
      !  NV,  NS, MUA,  NQ,  NB,  NO, NSB,  NT
      !
      ! NV    - supplementary VDW parameters
      ! NS    - supplementary BOND strech parameters
      ! MUA   - not used
      ! NQ    - supplementary CHARge parameters
      ! NB    - supplementary ANGL bending parameters
      ! NO    - supplementary OOPL parameters
      ! NSB   - not used
      ! NT    - supplementary TORSional parameters
      !
      VDW
         0.25      0.2       12.       0.8        0.5
         99     0.100     0.100     0.100     0.000 - DUMMY
      BOND
      0   5   99     1.000     0.500   parameters for dummy atoms
      ANGLE
      0   1    5   99     0.100   120.000   parameters for dummy atoms
      TORSION
      0  99    5    1    5   0.000   0.000   0.100   parameters for dummy atoms
      0  99    5    1    6   0.000   0.000   0.100   parameters for dummy atoms
      ----------------------------------------------------------------------------


Major features NOT currently implemented in CHARMM/MMFF:

(1) bonds between primary atoms and image atoms.
(2) Some cutoff options.  In particular,
    group-based cutoffs are not supported.
(3) Fast multipoles.


Other known limitations:

(1) correlation analysis tools have not been implemented for MMFF specific 
    energy terms -- e.g. it is not possible to calculate the correlation
    function for an out-of-plane bending angle, etc ...     
(2) .mrk files do not have group information -- i.e. residues = groups
(3) only all-atom models (no extended atoms)

There are probably other problems/limitations/bugs. Your comments about 
limitations of the current MMFF implementation in CHARMM (and bugs) will be 
very valuable.

Similarly, comments about deficiencies (as well as of particular strengths!) 
of the current MMFF parametrization would be very valuable for Tom Halgren, 
the author of MMFF.

Please direct comments to:

::

   Ryszard Czerminski, MSI
   e-mail: ryszard@msi.com
   phone:  (617)229-8875 x 217

   Tom Halgren, Merck Research Laboratories.
   e-mail: halgren@merck.com
   phone: (908) 594-7735

KNOWN BUGS:


.. mmff_thoery:

Theory
======

::

                      The Merck Molecular Force Field (MMFF94)

          A Broadly Parameterized, Computationally Derived Force Field 
                     for Organic and Bio-organic Systems 

                              Thomas A. Halgren
 
                 Merck Research Laboratories, Rahway, New Jersey 07065

                                February, 1996


1. Introducing The Merck Molecular Force Field.

   The Merck Molecular Force Field (MMFF) represents a systematic attempt
   to combine the best features of such well-regarded force fields as MM3,
   OPLS, AMBER, and CHARMM into a *single* force field that is equally
   adept in small-molecule and macromolecular applications.  In particular,
   MMFF strives for MM3-like accuracy for small molecules in a force field
   that can be used with confidence in condensed-phase simulations.  

   References to five papers introducing MMFF94 are given elsewhere within
   this documentation.

2. The Basis and Motivation for the Formulation of MMFF.

   Ideally, a single molecular mechanics/dynamics force field would reproduce 
   all of the following, and other, molecular properties accurately both in 
   gas-phase and in condensed-phase simulations: 

   * molecular geometries
   * conformational and stereoisomeric energies
   * torsional barriers and torsion-profile energies
   * intermolecular-interaction energies
   * intermolecular-interaction geometries
   * vibrational frequencies
   * heats of formation

   Because of their relatively simple construction, however, current force 
   fields necessarily make a variety of compromises.  MMFF94 focusses on 
   accurately reproducing conformational and intermolecular-interaction energies.
   It also regards molecular geometries, torsional barriers, and intermolecular-
   interaction geometries as being relatively important.  Vibrational 
   frequencies have been parameterized against a combination of theoretical
   and experimental data, but are regarded as being less important.  Heats of 
   formation are not normally needed to understand such qunatities as differences 
   in ligand-enzyme binding energies, and are not addressed in MMFF. 

   To be widely applicable, MMFF could not be parameterized against experimental 
   data because far too little data of high quality are available, especially for 
   conformational and intermolecular-interaction energies.  Instead, MMFF has 
   been derived almost solely from computational data, though experimental data
   have been used liberally in its validation.

   Many of the processes we wish to model at Merck occur in condensed phases.  
   Like many other well-known force fields, MMFF therefore employs effective pair 
   potentials that reflect in an averaged sense the enhancement of the charge 
   distribution in a high-dielectric medium due to molecular polarizability; a 
   better, but still future, approach would of course be to include polarizability 
   explicitly.  


3. Discussion

   The principal distinguishing feature of MMFF is that it is primarily
   computationally derived.  This approach is made possible because of recent
   increases in computing power; it is made necessary because pertinent
   experimental data are lacking for many of the chemical structures a force field
   suitable for general use in chemical and pharmaceutical applications must be
   prepared to handle.  MMFF's parameterization utilizes a large amount of
   high-quality computational data -- ca. 500 molecular structures optimized at
   the HF/6-31G* level, 475 structures optimized at the MP2/6-31G* level, 380
   structures evaluated at the composite "MP4SDQ/TZP" level using MP2/6-31G*-
   optimized geometries, and 1450 structures evaluated in single-point 
   calculations at the MP2/TZP level. This core has been significantly expanded 
   by using data from approximately 2800  Cambridge Structural Database 
   structures in conjunction with additional computational data and with a 
   series of carefully calibrated empirical rules and default-parameter 
   assignment procedures.  This expanded parametrization embraces nearly all 
   stable organic compounds in a systematic, objective, and consistent way, 
   making "missing parameters" virtually a thing of the past.

   The computationally derived "core" MMFF parameters cover a broad range
   of functional groups.  Among "monofunctional" chemical families, MMFF has been
   parameterized for alkanes, alkenes, alcohols, phenols, ethers, aldehydes,
   ketones, ketals, acetals, hemiketals, hemiacetals, amines, amides, peptides,
   ureas, imides, carboxylic acids, esters, carboxylate anions, ammonium cations,
   thiols, mercaptans, disulfides, halides (chlorides and fluorides), imines,
   iminium cations, amine N-oxides, hydroxylamines, hydroxamic acids, amidines,
   guanidines, amidinium cations, guanidinium cations, imadazolium cations,
   aromatic hydrocarbons, and heteroaromatic compounds.  The structural coverage
   is quite broad for many of these chemical families, but still is somewhat 
   limited for others. 

   Many of the bifunctional compounds included in the parameterization are
   unsaturated analogs of families listed above, i.e.: conjugated alkenes and
   aromatic hydrocarbons (e.g., styrenes); alpha,beta-unsaturated variants of
   amides, imines, aldehydes, ketones, carboxylic acids, esters, and carboxylate
   anions; vinylic ethers, alcohols, amines and esters; and allylic aldehydes,
   ketones, amines and alcohols.  Other bifunctional compounds include:
   beta-ketoacids; beta-hydroxyesters; dicarboxylic acids; 1,2-diols, 1,2-diamines
   and 1,2-dithiols; and nonconjugated dienes.  A limited selection of alkanes,
   amines, ketones, halides and ethers containing 4- or 5-membered rings has also
   been included. Compounds containing SO2 and phosphate groups have been 
   parameterized as a part of the extension of MMFF's parameterization mentioned 
   above. 

   Another important advantage of MMFF is that nearly all of its parameters
   have been determined in a mutually consistent fashion from the full set
   of computational data.  In most other force fields, parameters are
   determined for one functional group at a time, and then frozen before
   moving on to the next functional group.  This approach fails to allow
   for correlations that can make one subset of the parameters
   inappropriate for fitting data on subsequent functional groups.  MMFF's
   derivation, in contrast, simultaneously employed all data (e.g., on
   conformational energies) in determining the associated parameters (e.g.,
   torsion).  Furthermore, the parameter derivation procedures were
   iterated between three and four times, in order to allow each class of
   parameters (e.g., bond and angle reference values, quadratic force 
   constants, charges, torsion parameters) to be determined in a mutually 
   consistent fashion in the context of successively refined values for 
   parameters belonging to other classes.

   The reliance almost solely on computational data, the quality and
   quantity of the supporting *ab initio* calculations, and the methodology
   used in deriving mutually consistent values for most classes of
   parameters, together with novel elements of its functional form, combine
   to make MMFF's derivation unusual and possibly unique.  They also
   combine to produce a force field that by contemporary standards
   performs very well.  MMFF reproduces the computational data used in its
   parameterization with rms deviations of 0.006 angs for bond lengths,
   1.16 deg for bond angles, 5 deg for most torsion angles, 0.31 kcal/mol
   for conformational energies, and 0.50 kcal/mol for comparisons of
   relative energies along torsion profiles.  Crucially important
   intermolecular-interaction energies and geometries closely adhere to
   benchmarks established using *ab initio* calculations on small-molecule dimers. 
   Molecular charge distributions are also described reasonably well: rms 
   deviations are 0.39 D for HF/6-31G* molecular dipole moments and 5.5 deg 
   for dipole directions.  

   In addition, MMFF predicts experimental bond lengths, bond angles, and 
   vibrational frequencies essentially as accurately as does MM3, and 
   reproduces conformational energies and rotational barriers to 
   0.4 kcal/mol rms, about as well as can be expected given the disparate
   nature and uncertain accuracy of the experimental results. These results 
   are encouraging, because they demonstrate that fitting MMFF to high-quality 
   theoretical data has simultaneously conferred the ability to fit experiment.  
   In contrast to experimentally derived force fields, MMFF's great strength is 
   that it can be expected to perform equally well for the wide range of systems 
   for which it has been parameterized but for which no experimental data are 
   available. 

   I expect a computational approach like the one employed for MMFF to be
   indispensable in future efforts to derive still more accurate force
   fields which, for example, may explicitly incorporate polarizability and
   represent the electrostatic potential more accurately than is possible
   using only atom-centered charges.  Fortunately, further improvements in
   computer technology can be expected to make it increasingly feasible
   both to utilize the more complex force fields which result and to employ
   even more rigorous computational models to generate the data needed to
   parameterize them.  I doubt that any other approach will be capable of
   producing a physically superior force field which not only performs
   accurately in condensed-phase simulations but is parameterized sufficiently 
   broadly to support the full range of significant pharmaceutical, organic and 
   biochemical applications.


.. _mmff_funcform:

MMFF Functional Form
====================

The MMFF energy expression can be written as

.. math::

   E_{MMFF} = \sum E^{bond}_{ij} & + \sum E^{angle}_{ijk} + \sum E^{bend}_{ijk} + \sum E^{oop}_{ijk;l} \\
                                 & + \sum E^{torsion}_{ijkl} + \sum E^{vdW}_{ij} + \sum E^{Q}_{ij} \quad \mbox{(1)}

where the constituent terms, each expressed in kcal/mol, are defined as 
shown below.

1. Bond Stretching.  MMFF employs the quartic function:  

   .. math::
   
      E^{bond}_{ij} = 0.5 * 143.9325 * k^{bond}_{IJ} * \Delta r_{ij} *(1 + cs * \Delta r_{ij}^2 + \frac{7}{12} cs^22 * \Delta r_{ij}^2) \quad \mbox{(2)}
      
   where :math:`k^{bond}_{ij}` is the force constant in md/angs, :math:`\Delta r_{ij} = r_{ij} - r^{ref}_{ij}` is the
   difference in Angstroms between actual and reference bond lengths, and :math:`cs = -2 \AA^{-1}` is the "cubic stretch" constant.  This function corresponds to an
   expansion through fourth order of a Morse function with an "alpha" of :math:`2 \AA^{-1}`.
   Results published in a recent high-level *ab initio* study [1]_ show
   this value for alpha to be a representative one.  Note: throughout this
   Account, the indices i, j, k, ... represent atoms and I, J, K, ... denote the
   corresponding numerical MMFF atom types. 

2. Angle Bending.  MMFF normally uses the cubic expansion:

   .. math::
   
      E^{angle}_{ijk} = 0.5 * 0.043844 * k^{angle}_{IJK} * \Delta \theta_{ijk}^2 * (1 + cb * \Delta \theta_{ijk}) \quad \mbox{(3)}
 
   where :math:`k^{angle}_{IJK}` is the force constant in md-ang/rad**2, :math:`\Delta \theta_{ijk} = \theta_{ijk} - \theta_{IJK}^{ref}` is
   the difference between actual and reference bond angles in degrees, and
   :math:`cb = -0.007 deg^{-1}` is the "cubic-bend" constant.  

   For linear or near-linear bond angles, MMFF instead employs the well-behaved 
   form used in DREIDING [2]_ and UFF [3]_:

   .. math::
   
      E_{ijk}^{angle} = 143.9325 * k^{angle}_{IJK} * (1 + \cos(\theta_{ijk})) \quad \mbox{(4)}

3. Stretch-Bend Interactions.  MMFF employs the form:

   .. math::
      
      E_{ijk}^{bend} = 2.51210 * ( k^{bend}_{IJK} * \Delta r_{ij} + k^{bend}_{KJI} * \Delta r_{kj}) * \Delta \theta_{ijk} \quad \mbox{(5)}
   
   where :math:`k^{bend}_{IJK}` and :math:`k^{bend}_{KJI}` are force constants in md/rad which couple the i-j and 
   k-j  stretches to the i-j-k bend, and :math:`\Delta r_{ij}`, :math:`\Delta r_{jk}` and :math:`\Delta \theta_{ijk}` are as defined
   above. Stretch-bend interactions are omitted for linear bond angles. 

4. Out-of-Plane Bending at Tricoordinate Centers.  MMFF uses the form:

   .. math::
   
      E^{oop}_{ijk;l} = 0.5 * 0.043844 * k^{oop}_{IJK;L} * X_{ijk;l}^2 \quad \mbox{(6)}
      
   where :math:`k^{oop}_{IJK;L}` is the force constant in md-angs/rad**2 and :math:`X_{ijk;l}` is the
   Wilson angle [4]_ in degrees between the bond j-l and the plane i-j-k. Because
   it uses eq 3 for the "in-plane" angles, MMFF is able to properly describe the
   nonplanar centers found, e.g., in enamines, sulfonamides, and even amides. 

5. Torsion Interactions.  MMFF uses the three-fold representation employed 
   in MM2 and MM3, where W is the i-j-k-l dihedral angle:

   .. math::
   
      E_{ijkl}^{torsion} = 0.5 * (V_1 (1 + \cos W) + V_2 (1 - \cos {2W}) + V_3 (1 + \cos {3W})) \quad \mbox{(7)}
      
6. Van der Waals Interactions.  MMFF employs the recently developed
   "Buffered 14-7" form (eq 8) together with an expression which relates the
   minimum-energy separation R*II to the atomic polarizability aI (eq 9), a
   specially formulated combination rule (eqs 10, 11), and a Slater-Kirkwood
   expression for the well depth epsIJ (eq 12) [5]_: 

   ::
   
      Evdwij  =  epsIJ*{1.07R*IJ/(Rij+0.07R*IJ)}**7 *  
                       {1.12 R*IJ**7/(Rij**7 + 0.12R*IJ**7) - 2}                (8)

      R*II = AI * aI**(0.25)                                                    (9)

      R*IJ =  0.5 * (R*II + R*JJ) * (1 + 0.2 (1 - exp(-12*gIJ**2)))            (10) 

      gIJ = (R*II - R*JJ)/(R*II + R*JJ)                                        (11)

      eIJ =  181.16*GI*GJ*aIaJ/[(aI/NI)**0.5 + (aJ/NJ)**0.5]*R*IJ**(-6)        (12)

   Most vdW well depths and radii conform to simple systematic trends 
   adduced from high-quality experimental data on vdW interactions of rare-
   gas atoms and of small molecules with one another [5]_

7. Electrostatic Interactions.  MMFF uses the buffered Coulombic 
   form

   .. math::
   
      E^{Q}_{ij} = 332.0716* \frac{ q_i q_j }{D*(R_{ij} + d)} \quad \mbox{(13)}

   where :math:`q_i` and :math:`q_j` are partial atomic charges, :math:`R_{ij}` is the internuclear
   separation in angs, d = 0.05 angs is the "electrostatic buffering" constant, and D is
   the "dielectric constant" (normally taken as D = 1, though use of a distance-
   dependent dielectric constant is also supported).  Partial atomic charges :math:`q_i`
   are constructed from initial full or fractional formal atomic charges (usually
   zero, but, e.g., -0.5 for carboxylate oxygens) by adding contributions from
   bond charge increments wKI which describe the polarity of the bonds to atom i
   >from attached atoms k.  Specifically, MMFF computes :math:`q_i` as 

   .. math::
      
      q_i = q_{0i} + \sum wKI \quad \mbox{(14)}

   where wIK= - wKI. 1,4-interactions are scaled by a factor of 0.75.  Distance 
   buffering (d > 0) prevents infinite attractive electrostatic energies from 
   overwhelming the bounded repulsive vdW interaction given by eq 8 as 
   oppositely charged atomic centers approach.  

   Unlike MM2 and MM3, MMFF employs a unit dielectric constant, and
   thereby allows straightforward application to condensed-phase simulations
   employing explicit solvent molecules.  Like AMBER [6]_, CHARMM [7]_, OPLS [8]_ and
   other force fields used in molecular dynamics simulations, MMFF describes
   hydrogen bonding interactions as being essentially electrostatic in nature,
   whereas MM2 (1987 parameters and later) and MM3 in some cases attribute a
   significant portion of the stabilization energy to an attractive vdW term which
   would not be attenuated upon immersion in a high-dielectric medium.  This
   difference, too, may serve to make MMFF more readily applicable to
   condensed-phase simulations. 


References:

.. [1] Orozco, M.; Luque, F. J. J. Comput. Chem. 1993, 881-894.

.. [2] Mayo, S. L.; Olafson, B. D.; Goddard III, W. A. J. Phys. Chem. 1990, 94,
       8897. 

.. [3] Rappe, A. K.; Casewit, C. J.; Colwell, K. S.; Goddard III, W. A; Skiff, W.
       M. J. Am. Chem. Soc. 1992, 114, 10024-10035, and references therein. 

.. [4] Wilson, E. B., Jr; Decius, J. C.; Cross, P. C., Molecular Vibrations;
       Dover: New York, 1955, Chapter 4. 

.. [5] Halgren, T. A. J. Am. Chem. Soc. 1992, 114, 7827-7843.

.. [6] Weiner, S. J.; Kollman, P. A.; Nguyen, D. T.; Case, D. A. J. Comput. Chem. 
       1986, 7, 230-252;  Weiner, S. J.; Kollman, P. A.; Nguyen, D. T.; Case, D. A.; 
       Singh, U. C.; Ghio, C.; Alagona, G.; Profeta, S.; Weiner, P. J. Am. Chem. Soc. 
       1984, 106, 765-784. 

.. [7] Brooks, B. R.; Bruccoleri, R. E.; Olafson, B. D.; States, D. J.;
       Swaminathan, S.; Karplus, M. J. Comput. Chem. 1983, 4, 187-217. 

.. [8] Jorgensen, W. L.; Tirado-Rives, J. J. Am. Chem. Soc. 1988, 110, 1657-
       1666, and references therein.

.. mmff_refs:

References
==========

The following five papers introduce the MMFF94 force field:

[1] "Merck Molecular Force Field. I. Basis, Form, Scope, Parameterization, and
    Performance of MMFF94," Thomas A. Halgren, J. Comput. Chem., 17, 490-519 
    (1996).

[2] "Merck Molecular Force Field. II. MMFF94 van der Waals and Electrostatic
    Parameters for Intermolecular Interactions," Thomas A. Halgren, J. Comput. 
    Chem., 17, 520-552 (1996)

[3] "Merck Molecular Force Field. III. Molecular Geometries and Vibrational 
    Frequencies for MMFF94," Thomas A. Halgren, J. Comput. Chem., 17, 553-586 
    (1996).

[4] "Merck Molecular Force Field. IV. Conformational Energies and Geometries 
    for MMFF94," Thomas A. Halgren and Robert B. Nachbar, J. Comput. Chem., 17, 
    587-615 (1996).

[5] "Merck Molecular Force Field. V. Extension of MMFF94 Using Experimental
    Data, Additonal Computational Data, and Empirical Rules," Thomas A. Halgren, 
    J. Comput. Chem., 17, 616-641 (1996).



