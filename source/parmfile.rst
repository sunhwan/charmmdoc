.. py:module::parmfile

==========================================
CHARMM Emprical Energy Function Parameters
==========================================

::

   CHARMM Element doc/parmfile.doc 1.2
   
   File: Parmfile, Node: Top, Up: (chmdoc/commands.doc), Previous: (chmdoc/usage.doc)Standard Files, Next: Overview

                     CHARMM Emprical Energy Function Parameters

           This section describes parameters in the CHARMM empirical
   energy function.

   * Menu:

   * Overview::      Overview of CHARMM parameter file
   * Multiple::      Rules for the use of multiple dihedrals in CHARMM22
   * Conversion::    Rules for conversion of old nucleic acid rtf and
                     param to CHARMM22 format 
   * PARMDATA::      Description of Parameter Files available for general use.

   
   File: Parmfile, Node: Overview, Up: Top, Previous: Top, Next: Multiple

                       Overview of CHARMM parameter files

                    By Alexander D. MacKerell Jr., July 1997
                          Updated; January 2010

           This section of the documenation contains a brief description
   of the contents of a parameter file.  The CHARMM parameter file
   contains the information necessary to calculate energies etc. when
   combined with the information from a PSF file for a structure.
   Information on the keywords found in the parameter file is in IO.DOC.

   (A)   * CHARMM example parameter file
         *
   (B)   BOND
         H   O   500.0  1.00
   (C)   ANGLe (THETa)
         H   O   H  100.0  104.51  20.0  1.70
   (D)   DIHEdral (PHI)
         HT  CT  CT  HT    10.0   3    180.0
         X   CT  CT  X     10.0   3    180.0
   (E)   IMPH
         O   C   CT  N      5.0   1      0.0
         X   C   CT  X      5.0   1      0.0
         X   X   CT  N      5.0   1      0.0
         O   X   X   N      5.0   1      0.0
   (F)   CMAP 
         CT1  N   CT1  C  N  CT1  C  CT1 24
         cmap data: 2D array of energy correction values
   (G)   NBONDed nonbond-spec
         H    0.00  -0.046  0.2245   0.00  -0.023  0.2245
         O    0.00  -0.120  1.8000   0.00  -0.060  1.8000
   (H)   NBFIX
         H    O   -0.30  1.50   -0.15  1.50
   (I)   HBONDs  hbond terms (IO.DOC)
         H    O   -0.00  1.00
   (J)   END

           The parameter file starts with a title (A) which contains
   information on the origins and applicability of that file.

           Section (B) BONDs, contains information on all bond force
   constants and equilibrium geometries.  In this as well as the
   remainder of the parameter file the bonds etc. are specified by the
   atom type associated with each IUPAC atom in the topology file.

           Section (C) ANGLes or THETas, are specified by 3 atom types
   followed by the force constant and equilibrium geometry.  If a
   Urey-Bradley term is desired between the 1 and 3 atom types of the
   angle a second U-B force constant and equilibrium geometry are
   included.

           Section (D) DIHEdrals (PHI), contains the 4 atom types
   specifing a dihedral followed by the force constant, the multiplicity
   of the dihedral and the minimium geometry of the dihedral.  With
   dihedrals wildcards, X, as shown may be included for the terminal
   atoms.  Also, multiple dihedrals of different multiplicities may be
   specified for a single dihedral as outlined below. As of version 25
   any value of the phase (besides 0 and 180) is possible.  However, the
   use of phases other than 0 and 180 is discouraged, as the resulting
   potential is asymmetric, and hence unphysical in most (if not all)
   cases.  See io.doc for more information.

           Improper dihedrals (E) IMPH, used for out of plane motions are
   specified in the same fashion as dihedrals.  The use of wildcards, X,
   is also allowed in a number of variations.  Multiple improper
   dihedrals are not supported.  Ordinarily, improper dihedrals are given
   a multiplicity of 0, which imposes a harmonic restoring potential
   instead of a cosine function.  In this case, the central atom must be
   either the first or the last atom in the definition, else the
   resulting potential will be asymmetric and exhibit a discontinuity at
   equilibrium!  By convention, the first atom of an improper dihedral
   (type A-X-X-B or A-B-C-D) should be the central atom.  See io.doc for
   more information.

           Cross-term energy correction map, CMAP (F), is a 2D grid
   correction where the 2 dimensions correspond to two dihedral angles,
   where the first four atom types correspond to the first dihedral and
   the next four atom types to the second dihedral angle.  A CMAP
   specification is required in the topology file (see rtop.doc) and
   properly generated for this energy term to be applied.  The integer
   following the 8 atom types represents the dimensions of the 2D grid.
   In the present example, 24 indicates that for the 360x360 deg surface,
   grid data in 15 deg increments is being used. The data in the CMAP
   grid is presented starting from -180 down to 0 and then increased to
   the final value (165 for the present example).

           Parameters for (G) NONBonded VDW parameters may be specified
   in two ways.  Initially the Tanford-Kirkwood Formula was used where
   the atom polarizabilities, Number of effective electrons, and (minimum
   radius)/2 were required.  In this formulation the first term following
   the atom type is the atom polarizability, the second term is the
   number of effective electrons and must be positive in order to specify
   the Tanford-Kirkwood Formula and the third term is the (minimum
   radius)/2.  If the second term is negative, then the first number is
   ignored, the second term is the well-depth (epsilon) and the third
   term is the (minimum radius)/2.  Both formulations use the
   Lennard-Jones 6-12 formula to determine the VDW interactions, in the
   first method the Tanford-Kirkwood Formula is used to calculate the
   well-depth (epsilon) and in the second method it is used directly.
   With both formulations a second set of 3 numbers may be specified to
   indicate the VDW parameters to be used for the calculation of 1-4
   nonbonded interactions. Wildcards (*, %, etc. see MISCOM.DOC) may be
   used with the NONBond as well as the NBFIX and HBOND sections of the
   parameter file.

           The NBFIX section (H) allows VDW interactions between specific
   atom pairs to be modified.  This is done by specifing the 2 atom types
   followed by the well depth and the minimum radius (not (minimum
   radius)/2 as in NBOND).  A second well depth and minimum radius may be
   specified to determine the 1-4 interactions.

           The final section (I) contains the hydogen bond well depths
   and minimum radii for various atom pairs.  In current versions of the
   CHARMM parameter sets (PARAM19, CHARMM22 protein and nucleic acid
   parameters) hydrogen bonding is included implicitly in the
   electrostatic and VDW interactions.  Thus, the HBOND well depth is set
   to -0.00 and in most calculations IHBFRQ should be set to 0 to avoid
   updating the hydrogen bond lists. This facility is still supported to
   allow calculations using the Lennart Nilsson nucleic acid parameters,
   AMBER parameters and for analysis of hydrogen bond geometries.  It
   should be noted that both the NBOND and HBOND keywords are followed by
   a number of keywords dictating truncation schemes, 1-4 interaction
   treatments and dielectric constants, amoung others.  These
   specifications are of the upmost importance for relabile calculations
   and deviations from the default values supplied with the parameter
   files should be done with the utmost caution. As of version 27, a
   special parameter files, par_hbond.inp, was added to the toppar
   directory; it contains information to calculated hydrogen bonds. Note
   that the resulting energetics are meaningless and are only included to
   aid in the analysis of hydrogen bonds.

   
   File: Parmfile, Node: Multiple, Up: Top, Next: Conversion, Previous: Overview

             Rules for the use of multiple dihedrals in CHARMM24


    1) The association of 1 or more dihedrals with different 
       multiplicities to a specfic dihedral type (as specified 
       by atom types) is specified by the presence of 2 or 
       more dihedral parameters in the parameter file.  When 
       multiple dihedrals are read in the parameter file and if the
       warning level (wrnlev) is 6 or more CHARMM
       will list those dihedrals in the output file (Note: the following
       type message "PARRDR> Multiple terms for dihedral type: INDEX  427
       CODE31141959     CT3 -OS  -CD  -OB" indicates that the multiple
       dihedral has been successfully read). 

    2) If dihedral angles are AUTOGENERATED, then the RTF should
       not specify them again. Additional dihedrals in the RTF will 
       be ignored and warnings given.

    3) Without AUTOGENERATE, each dihedral should appear only once 
       in the RTF.  Multiple listings of a dihedral will be ignored
       and warnings given.

    4) The order of the multiple dihedral entries associated with
       a specific dihedral is important; they must be placed
       sequentially in the parameter file.  If they are not sequential
       errors will be given.  This is new in C24B1 and later versions.
       For example:
       P    ON2  P2   ON2      0.03    2     0.0  
       P    ON2  P2   ON2      0.03    3     0.0  
       will place both a 2-fold and a 3-fold term on the P-ON2-P2-ON2
       dihedral.

    5) Wildcards may be used in the parameter file to specify multiple
       dihedrals(ie. X  C1  C2  X), however, all the dihedrals in the
       parameter file associated with that dihedral type must be
       wildcards.  Use of wildcards with multiple dihedrals is NOT
       recommeded.

    6) Specific dihedral entries always override wildcard entries.
       For example:
       X  C2  C3  X   100.0 1 180.0
       C1 C2  C3  C4  100.0 2 180.0
       X  C2  C3  X   100.0 3 180.0
       will assign the 2-fold term to C1-C2-C3-C4 while 1-fold and
       3-fold terms would be assigned to C5-C2-C3-C6 and any other
       dihedral centered about the C2-C3 bond.  This assignment of
       the multiple terms to a number of dihedrals is why the use
       wildcards for the specification of multiple dihedrals in NOT
       recommeded.  The preferred method is as follows:
       X  C2  C3  X   100.0 2 180.0
       C5 C2  C3  C6  100.0 1 180.0
       C5 C2  C3  C6  100.0 3 180.0
       will assign the 1-fold and 3-fold terms to C5-C2-C3-C6 and the
       2-fold term to C1-C2-C3-C4 and any other dihedral centered about 
       the C2-C3 bond.  This limits the potential for multiple dihedrals
       being mistakenly assigned to a dihedral centered on the C2-C3 bond.
       Thus, it is advised that when creating a multiple dihedral all
       4 atom types be explicitly stated and, if necessary, new atom
       types be created to avoid conflicts.

    7) This design is such that previous CHARMM topology and parameter
       files for proteins are compatible with CHARMM24.  However, due
       to complexities in the multiple dihedral setup for the nucleic
       acid sugars (ribose and deoxyribose) the nucleic acid topology 
       and parameter files are NOT compatible with CHARMM22.  In order
       to make them compatible the following alterations must be 
       performed.  Alternatively, the altered files may be obained 
       from Alexander D. MacKerell Jr.

   
   File: Parmfile, Node: Conversion, Up: Top, Previous: Multiple, Next: PARMDATA

    Rules for conversion of old nucleic acid rtf and param to CHARMM22 format

    The following conversion rules apply to CHARMM22.  Compatability with
    C24B1 and later versions will be insured if the multiple dihedrals in
    the converted parameters are sequential, as disscused above.

           ALL-HYDROGEN

           Protocol for conversion of all-hydrogen nucleic acid topology and
   parameter files (topnah*.inp and parnah*.inp) from a CHARMM21 or
   previous format to a format compatible with CHARMM22.  This change is
   due to a new methodology for the treatment of multiple dihedrals in
   CHARMM22.

           In Topology File (TOPNAH1.INP, TOPNAH1E.INP, TOPNAH1R.INP)

       1)  Create a new atom type, OSS

       2)  Convert the atom type of all O4' atoms to OSS

           In Parameter File (PARNAH1.INP)

       1)  Copy all OS parameters (bonds, angles, dihedrals etc.)
           and in the copy change OS to OSS.  Be sure that the
           original OS parameter remains.  Some OS to OSS copies
           can be avoided (such as OS  P terms), however, one
           must be careful that all the necessary OSS parameters
           relating to O4' are present.  Creating extra OSS 
           parameters which are unused is not a problem. One
           exception occurs with the dihedral OS CH CH OS, where
           only one of the terminal OS atom should be converted
           to OSS.

       2)  In the DIHEDRAL (PHI) parameters under the heading
           "WILMA OLSON SUGAR MODEL" the following steps must be
           performed once all the OSS dihedral parameters are
           created.

       A)  In all the explicit OS terms which don't include
           wildcards (X) or P atom types and have both 2 and
           3-fold periodicities (2nd of 3 numbers following the
           dihedral) the 2nd 3-fold term must be commented out 
           with a !.

       B)  Of the new explicit OSS terms the following
           3-fold terms must be commented out with a !.

           OSS  CH   CH   OS       1.4000    3    0.0000 
           OH   CH   CH   OSS      1.4000    3    0.0000 


           Lastly, when generating the structure be sure only the AUTOGENERATE
   ANGLE term is used. (i.e. do NOT use AUTOGENERATE DIHEDRAL).
      
           At this point the topology and parameter files should be
   compatible with CHARMM22 (but not CHARMM21 or a previous version of
   CHARMM).  A test should be performed on a (deoxy)ribose containing
   containing compound.  In this test the energies should be calculated
   1) using CHARMM21 or a previous version using the original, unmodified
   topology and parameter files and 2) with CHARMM22 using the modified
   OSS containing topology and parameter files. These energies should be
   equivalent.


           EXTENDED (UNITED) ATOM

           Protocal for the conversion of extended (united,explicit) atom
   nucleic acid topology and parameter files from CHARMM21 or previous
   format to a format compatible with CHARMM22.  This change is due to a
   new methodology for the treatment of multiple dihedrals in CHARMM22.

           In Topology File (TOPRNA10 or TOPRNA10R)

       1)  Create 2 new atom types, OSS and OST

       2)  Convert the atom type of all O4' atoms to OSS except
           in the the patch PRES DEOX where it must be changed
           to atom type OST.  This conversion to OST must also be 
           performed in any residue, such as RESI DRIB, in which 
           deoxyribose is used explicitly.

       3)  In the patch PRES DEOX add the line: 

           ATOM O4'   OST   -0.30  ! (check the charge)

           before the GROUP statement and comment out the terms

     !DELETE DIHE O4'  C4'  C3'  O3' ! WE NEED THIS AS A MULTIPLE TERM IN DEOXY
     !DIHE O4'  C4'  C3'  O3' ! threefold
     !DIHE O4'  C4'  C3'  O3' ! twofold

           such that no alterations in the dihedral setup are made.

           In Parameter File (PARDNA10.INP)

       1)  Copy all OS parameters twice (bonds, angles, dihedrals etc.);
           in the first copy change OS to OSS and in the second change OS
           to OST.  Be sure that the original OS parameter remains.  Some 
           OS to OSS(OST) copies can be avoided (such as terms in which OS 
           is adjacent to P), however, one must be careful that all the 
           necessary OSS(OST) parameters relating to O4' are present.  
           Creating extra OSS(OST) parameters which are unused is not a 
           problem. One exception occurs with the dihedral OS CH CH OS, 
           where only one of the terminal OS atom should be converted
           to OSS(OST).

       2)  In the DIHEDRAL (PHI) parameters under the heading
           "WILMA OLSON SUGAR MODEL" the following steps must be
           performed once all the OSS(OST) dihedral parameters are
           created.

       A)  In all the explicit OS terms which don't include
           wildcards (X) or P atom types and have both 2 and
           3-fold periodicities (2nd of 3 numbers following the
           dihedral) the 2nd term must be commented out 
           with a ! (mostly 3-fold terms and 1 or 2 2-fold term).

       B)  Of the new explicit OSS terms the following 
           3-fold terms must be commented out with a !.

           OSS  CH   CH   OS       1.4000    3    0.0000 
           OH   CH   CH   OSS      1.4000    3    0.0000 
      
       C)  Maintain all of the OST dihedral terms.

           An example of the additions/alterations to pardna10.inp are
   listed below.

     BOND
     HO   OSS    450.0000    0.9600
     HO   OST    450.0000    0.9600
     OSS  CH     292.0000    1.4300 
     OSS  C2     292.0000    1.4300 
     OST  CH     292.0000    1.4300 
     OST  C2     292.0000    1.4300 
     C3   OSS    292.0000    1.38   
     C3   OST    292.0000    1.38   
     C    OSS    292.0000    1.43  
     C    OST    292.0000    1.43  

     THETA
     OSS  C2   C3     150.5000  111.0000  
     OSS  C2   CH      70.0000  112.0000
     OSS  C2   C2      82.0000  112.0000  
     OST  C2   C3     150.5000  111.0000  
     OST  C2   CH      70.0000  112.0000
     OST  C2   C2      82.0000  112.0000  
     C2   CH   OSS     46.5000  111.0000
     C2   CH   OST     46.5000  111.0000
     C3   CH   OSS     46.5000  111.0000
     C3   CH   OST     46.5000  111.0000
     CH   CH   OSS     46.5000  111.0000
     CH   CH   OST     46.5000  111.0000
     OSS  CH   NS      46.5000  111.0000
     OSS  CH   NH2E    46.5000  111.0000
     OST  CH   NS      46.5000  111.0000
     OST  CH   NH2E    46.5000  111.0000
     C2   OSS  C2      82.0000  111.5000  
     CH   OSS  CH      46.5000  111.5000
     HO   OSS  CH      46.5000  107.3000
     HO   OSS  C2      46.5000  107.3000
     C2   OST  C2      82.0000  111.5000  
     CH   OST  CH      46.5000  111.5000
     HO   OST  CH      46.5000  107.3000
     HO   OST  C2      46.5000  107.3000
     CH   OSS  C3      46.5     107.3
     CH   OST  C3      46.5     107.3
     C    OSS  C3      46.5     120.5 
     C    OST  C3      46.5     120.5 
     O    C    OSS     70.0     120.0
     O    C    OST     70.0     120.0
     CH   C    OSS     70.0     125.3 
     NA   C    OSS     70.0     120.0 
     CH   C    OST     70.0     125.3 
     NA   C    OST     70.0     120.0 
     OSS  CH   CS      46.5     111.0
     OST  CH   CS      46.5     111.0

     PHI
     X    CH   OSS  X        0.9000    3    0.0000  
     X    CH   OST  X        0.9000    3    0.0000  
     X    C2   OSS  X        0.5000    3    0.0000
     X    C2   OST  X        0.5000    3    0.0000
     ! OSS SUGAR TERMS
     OSS  CH   CH   OS       0.5000    2    0.0000 
     !OSS  CH   CH   OS       1.4000    3    0.0000 Should be commented out
     OH   CH   CH   OSS      0.5000    2    0.0000 
     !OH   CH   CH   OSS      1.4000    3    0.0000 Should be commented out
     OSS  CH   CH   CH       0.5000    2    0.0000 
     OSS  CH   CH   CH       1.4000    3    0.0000
     OSS  CH   C2   CH       1.0000    2    0.0000 
     OSS  CH   C2   CH       1.4000    3    0.0000
     OSS  CH   CH   C2       1.4000    3    0.0000 
     OSS  CH   CH   C2       0.5000    2    0.0000
     OSS  C2   C2   C2       1.4       3    0.0    
     OSS  C2   C2   C2       0.5       2    0.0    
     ! OST SUGAR TERMS
     OST  CH   CH   OS       0.5000    2    0.0000 
     OST  CH   CH   OS       1.4000    3    0.0000 
     OH   CH   CH   OST      0.5000    2    0.0000 
     OH   CH   CH   OST      1.4000    3    0.0000 
     OST  CH   CH   CH       0.5000    2    0.0000 
     OST  CH   CH   CH       1.4000    3    0.0000
     OST  CH   C2   CH       1.0000    2    0.0000 
     OST  CH   C2   CH       1.4000    3    0.0000
     OST  CH   CH   C2       1.4000    3    0.0000 
     OST  CH   CH   C2       0.5000    2    0.0000
     OST  C2   C2   C2       1.4       3    0.0    
     OST  C2   C2   C2       0.5       2    0.0    
     ! additional terms for tRNA
     OSS  CH   CS   CF       1.5       3       0.0
     OST  CH   CS   CF       1.5       3       0.0
     C2   CH   C    OSS      1.5       3       0.0
     C2   CH   C    OST      1.5       3       0.0
     X    C    OSS  X        1.8       2       180.00 
     X    C    OST  X        1.8       2       180.00 
     ! THE FOLLOWING TERMS UNDER THE HEADER
     ! "WILMA OLSON SUGAR MODEL":
     ! SHOULD BE COMMENTED OUT
     !OS   CH   CH   OS       1.4000    3    0.0000 
     !OS   CH   CH   CH       1.4000    3    0.0000
     !OH   CH   CH   OS       1.4000    3    0.0000 
     !OS   CH   C2   CH       1.4000    3    0.0000
     !OS   CH   CH   C2       0.5000    2    0.0000
     !OS   C2   C2   C2       0.5       2    0.0     

     IMPHI
     OSS  X    X    CH      31.5000 0  35.2600
     OST  X    X    CH      31.5000 0  35.2600
     CH   OSS  C2   NS      31.5000 0  35.2600
     CH   OSS  CH   NS      31.5000 0  35.2600
     CH   OSS  C2   NH2E    31.5000 0  35.2600
     CH   OSS  CH   NH2E    31.5000 0  35.2600
     CH   OST  C2   NS      31.5000 0  35.2600
     CH   OST  CH   NS      31.5000 0  35.2600
     CH   OST  C2   NH2E    31.5000 0  35.2600
     CH   OST  CH   NH2E    31.5000 0  35.2600

     NBONDED
     OSS      0.64      7.0       1.6
     OST      0.64      7.0       1.6

           Lastly, when generating the structure be sure only the
   AUTOGENERATE ANGLE term is used. (i.e. do NOT use AUTOGENERATE
   DIHEDRAL).

           At this point the topology and parameter files should be
   compatible with CHARMM22 (but not CHARMM21 or a previous version of
   CHARMM).  A test should be performed on a (deoxy)ribose containing
   containing compound.  In this test the energies should be calculated
   1) using CHARMM21 or a previous version using the original, unmodified
   topology and parameter files and 2) with CHARMM22 using the modified
   OSS containing topology and parameter files. These energies should be
   equivalent.

   
   File: Parmfile, Node: PARMDATA, Up: Top, Previous: Conversion, Next: Top

   Description of topology and parameter files in version c31 and
   subsequent versions.

   In version 31 a major restructuring of the topology and parameter
   files was undertaken. This was performed to create a more modular
   approach to the files, thereby avoiding the problem of files becoming
   increasingly large.  The restructuring was done such that the primary
   topology and parameter files for biomolecules can be used as they were
   previously.  However, the topology and parameter information for the
   majority of model compounds used in the parameter development,
   additional molecules, including coenzymes, and patches parametrized to
   be compatible with the CHARMM force fields have been moved to toppar
   stream files.  These toppar stream files have to be streamed in a
   CHARMM input script typically following reading of the parent topology
   and parameter files.  They include both the topology and parameter
   information for the selected molecules, using the "read rtf card
   append" and "read param card append" commands to append the additional
   information to the topology and parameter lists.  As of version c31 it
   is still necessary to maintain all the MASS atom lists in the parent
   topology and parameter files.

   Files in Version C35/C36

   Topology files
    top_all22_prot.inp       all hydrogen RTF for proteins, CHARMM22 with CMAP
    top_all35_carb.rtf       all hydrogen RTF for sugars
    top_all32_lipid.rtf      all hydrogen RTF for lipids with alkane dihedral update
    top_all36_lipid.rtf      all hydrogen RTF for lipids (C36 and on)
    top_all27_na.rtf         all hydrogen RTF for nucleic acids
    top_all32_na_lipid.rtf   all hydrogen RTF for nucleic acids and lipids
    top_all36_na_lipid.rtf   all hydrogen RTF for nucleic acids and lipids (C36 and on)
    top_all27_prot_na.rtf    all hydrogen RTF for proteins and nucleic acids
    top_all32_prot_lipid.rtf all hydrogen RTF for proteins and lipids
    top_all36_prot_lipid.rtf all hydrogen RTF for proteins and lipids (C36 and on)
    top_all35_ethers.rtf     all hydrogen RTF for ethers
    top_all30_cheq_prot.inp  all hydrogen RTF for protein charge equilibration polarizable model
    toph19.inp               extended atom RTF for proteins
    toprna10r_22.inp         extended atom RTF for nucleic acids

   Parameter files
    par_all22_prot.inp       all hydrogen parameters for proteins with CMAP
    par_all35_carb.prm       all hydrogen parameters for sugars
    par_all32_lipid.prm      all hydrogen parameters for lipids with alkane dihedral update
    par_all36_lipid.prm      all hydrogen parameters for lipids (C36 and on)
    par_all27_na.prm         all hydrogen parameters for nucleic acids
    par_all32_na_lipid.prm   all hydrogen parameters for nucleic acids and lipids
    par_all36_na_lipid.prm   all hydrogen parameters for nucleic acids and lipids (C36 and on)
    par_all27_prot_na.prm    all hydrogen parameters for proteins and nucleic acids
    par_all32_prot_lipid.prm all hydrogen parameters for proteins and lipids
    par_all36_prot_lipid.prm all hydrogen parameters for proteins and lipids (C36 and on)
    par_all30_cheq_prot.inp  all hydrogen parameters for protein charge equilibration polarizable model
    par_all35_ethers.prm     all hydrogen parameters for ethers
    param19.inp              extended atom parameters for proteins
    pardna10_22.inp          extended atom parameters for nucleic acids

   (C) Toppar stream files (see stream subdirectory) listed under the parent 
   topology and parameter files required for the individual files.

   Parent files: can be used with prot, na and lipid files

    toppar_dum_nobel_gases.str: dummy atom, helium and neon
    toppar_hbond.str: stream file to estimate hydrogen bond interactions

   Parent files: top_all22_prot.inp, par_all22_prot.inp
             (or top_all22_prot_cmap.inp, par_all22_prot_cmap.inp)

    toppar_all22_prot_model.str: model compounds used in protein parameter development
                                 as well as additional compounds
    toppar_all22_prot_aldehydes.str: small molecule aldehydes
    toppar_all22_prot_aliphatic_c27.str: extends all22 protein force field to include
                                 all27 alkane parameters        
    toppar_all22_prot_fluoro_alkanes.str: optimized fluoroalkanes, requires
                                 toppar_all22_prot_aliphatic_c27.str
    toppar_all22_prot_heme.str: heme, O2, CO, CO2 and related patches
    toppar_all22_prot_pyridines.str: various substituted pyridines
    toppar_all22_prot_retinol.str: retinol, model compounds, Schiff's bases

   Parent files: top_all27_na.rtf, par_all27_na.prm

    toppar_all27_na_model.str: model compounds used in na parameter development including
                                 individual bases etc.
    toppar_all27_na_base_modifications.str: various chemical modifications of bases
    toppar_all27_na_carbocyclic.str: constrained bicyclic sugars
    toppar_all27_na_nad_ppi.str: NAD, NADH, ADP, ATP and others.  Useful in combination with
                                 protein force field via the top_all27_prot_na.rtf
                                 and par_all27_prot_na.prm


   Parent files: top_all32_lipid.rtf, par_all32_lipid.prm 

    toppar_all27_lipid_model.str: model compounds used in lipid parameter development, including
                                  alkenes
    toppar_all27_lipid_cholesterol.str: cholesterol and related model compounds

   Parent files: top_all36_lipid.rtf, par_all36_lipid.prm 

    toppar_all36_lipid_model.str: model compounds used in lipid parameter development, including
                                  alkenes
    toppar_all36_lipid_cholesterol.str: cholesterol and related model compounds

   Parent files: top_all27_prot_na.rtf, par_all27_prot_na.prm

    toppar_prot_na_all.str: all compounds that require both protein and nucleic acid toppar information
                                 includes phosphorylated tyrosine, serine and threonine and some
                                 coenzymes (SAH)
    toppar_all27_na_bkb_modifications.str: various chemical modificaiton of the na backbone
                                 including abasic variants and phosphoramidate

   Parent files not needed

    toppar_water_ions.str: contains TIP3P water model and ions.  All of these
                                 are also included in the prot, na and lipid topology
                                 and parameter files.
    toppar_amines.str: highly optimized neutral aliphatic amines.


   (D) Parameters for the polarizable force field based on a classical
   Drude oscillator.  As this force field is currently under development
   such that the parameters have been placed in the "drude" subdirectory.
   Parameters for water, alkanes, ethers, aromatics, alcohols and amides
   are available as of January 2007 and this list will be expanding.  See
   the 00readme for details and the approriate references.

   (E) Parameters for selected silicate and aluminosilicate surfaces have
   been developed.  These parameters are designed to be compatible with
   the CHARMM22 and 27 force fields allowing for biological
   molecule-silicate surface interactions.  As use of these parameters
   requires creation of the surface, which entails creation of the
   necessary patches, the parameters are included in the "silicates"
   subdirectory.  This directory also includes examples and code to
   create the extended surfaces. See the 00readme file for more details.

   ref: Lopes, P.E.M., Murashov, V. Tazi, M. Demchuk, E. MacKerell,
   A. D., Jr. "Development of an Empirical Force Field for
   Silica. Application to the Quartz-Water Interface," Journal of
   Physical Chemistry B, 110: 2782-2792, 2006.

   (F) Additional topology and parameter files from various sources are
   included in subdirectories of the toppar directory.  A description
   follows:

   non_charmm: Contains toppar files for AMBER, Bristol-Myers Squibb
   (BMS) and OPLS force fields along with a stream file for the SPC and
   SPC/E water models.  These files have been tested to the extent that
   they may be considered reliable representations of the original force
   fields, though potentially not exact representations.  These files are
   NOT maintained and, thus, use at your own risk.  See the 00readme
   files and note that AMBER requires a special version of CHARMM as
   described in the 00readme file.

   tamdfff: An internal coordinate force field (ICFF) that was built
   based on the CHARMM 22 protein force field.  Specifically, it provides
   a backbone covalent geometry suitable for torsion angle molecular
   dynamics (TAMD) and the necessary CMAP cross-term corrections to
   suppress distortions of the potential energy surface due to rigid
   covalent geometry.  Additional details can be found in tamd.doc.

   Ref. J. Chen, W. Im and C. L. Brooks III, J. Comp. Chem. 2005, 26,
   1565-1578.

   rush: A simple implicit-solvent protein force-field that adds terms to
   the bonded portion (bond + angle + dihe + impr + urey) of the all-atom
   CHARMM22 force field to account for volume-exclusion (_R_epulsion),
   the hydrophobic effect (_U_nburied _S_urface), and intra-molecular and
   protein-solvent hydrogen-bonding (_H_ydrogen-bonding) (hence _R_ _U_
   _S_ _H_).  Usage instructions are in doc/rush.doc

   gbsw: Optimized protein backbone parameters (par_all22_prot_gbsw.inp)
   and atomic input radii (radius_gbsw.str) for a balanced GBSW implicit
   solvent force field.  The backbone phi/psi cross-term (CMAP) and the
   atomic input radii have been re-optimized specifically to balance the
   solvation and intramolecular interactions and to capture experimental
   conformational equilibria of both helical peptides and
   beta-hairpins. Additional information can be found in gbsw.doc.

   Ref. J. Chen, W. Im and C. L. Brooks III, J. Am. Chem. Soc. 128,
   3728-36 (2006).

   Description of topology and parameter files prior to version c31.
   These files can be accessed via the toppar_history subdirectory of the
   toppar directory.

   (A) Topology files
        top_all22_prot.inp       all hydrogen RTF for proteins
        top_all22_model.inp      all hydrogen RTF for protein model cmpds
        top_all22_sugar.rtf      all hydrogen RTF for sugars
        top_all27_na.rtf         all hydrogen RTF for nucleic acids
        top_all27_lipid.rtf      all hydrogen RTF for lipids
        top_all27_na_lipid.rtf   all hydrogen RTF for nucleic acids and lipids
        top_all27_prot_na.rtf    all hydrogen RTF for proteins and nucleic acids
        top_all27_prot_lipid.rtf all hydrogen RTF for proteins and lipids
        toph19.inp               extended atom RTF for proteins
        toprna10r_22.inp         extended atom RTF for nucleic acids

   (B) Parameter files
        par_all22_prot.inp       all hydrogen parameters for proteins
        par_all22_sugar.prm      all hydrogen parameters for sugars
        par_all27_na.prm         all hydrogen parameters for nucleic acids
        par_all27_lipid.prm      all hydrogen parameters for lipids
        par_all27_na_lipid.prm   all hydrogen parameters for nucleic acids and lipids
        par_all27_prot_na.prm    all hydrogen parameters for proteins and nucleic acids
        par_all27_prot_lipid.prm all hydrogen parameters for proteins and lipids
        param19.inp              extended atom parameters for proteins
        pardna10_22.inp          extended atom parameters for nucleic acids
        par_hbond.inp            hydrogen bond parameters for analysis only

   Phased out: The CHARMM22 all-atom nucleic acid and lipid topology and
   parameter files are no longer included in the toppar directory due to
   their becoming obsolete.  Note that they are included in the
   toppar_history directories.

   The CHARMM all-hydrogen topology and parameter sets may be
   considered to be stable, however, minor bug fixes may be performed as
   required.  Additions may also occur leading to an expanding set of
   parameters which are compatible across proteins, nucleic acids,
   lipids, and, ultimately, carbohydrates.  The carbohydrate(sugar)
   parameter work is still in progress by John Brady and coworkers; the
   number of sugar types should expand in the future.  See the file
   toppar_all.history for a listing changes in the files over time.
   top_all22_model.inp includes the majority of model compounds used in
   the protein parameterization and is to be used in conjunction with
   par_all22_prot.inp.  

   Three sets of combined topology and parameter files are included
   for use with 1) proteins and nucleic acids, 2) protein and lipids
   and 3) nucleic acids and lipids.  In all cases the CHARMM22 protein
   parameters and the CHARMM27 nucleic acid or lipid parameters are
   used.  The designation all27 for these files is based on the
   use of the most recent nucleic acid or lipid parameters.  Test
   calculations using these combined files have yielded good results.

   Added as of July 1997 was the parameter file par_hbond.inp, which has
   been renamed to stream/toppar_hbond.str.  This file is included for
   the analysis of hydrogen bonds; it includes information to calculate
   h-bond energies, but these are basically meaningless.  The hydrogen
   bonds should NOT be used for energy, minimization and dynamics
   calculations with the CHARMM all-hydrogen topology and parameter sets.

   Ions in the all22 files are from two sources.  Mg and Ca are from
   Prodhom and Karplus and were optimized specifically for the all22
   parameters.  The remaining cations are from Benoit Roux (see his
   thesis).  They were optimized to be consistent with Param19, however,
   MD studies in a number of groups have shown them to work well.  Note
   the presence of a variety of NBFIXES for the ions.  These were
   initially optimized based on the proteins and later transferred to the
   lipids and nucleic acids based on analogy (by ADM Jr.).  Ions in the 
   all27 files have been optimized based on free energies of solvation
   by Roux and coworkers.  As of August 1999 there were no NBFIXes
   used with these ions.

   The extended atom parameters for proteins are the same as those
   included with CHARMM20 which are based on Wally Reiher's thesis.  They
   have been included in the supplement material of a recent publication
   (see suggested citations below).  For the extended atom nucleic acid
   parameters those of Nilsson and Karplus, J. Comp.  Chem.  7:591-616,
   1986 are used which were also included in the CHARMM20 release and are
   the only set to include explicit hydrogen bonding terms.  Some
   alterations of the extended atom nucleic acid topology and parameter
   files have been made in order to maintain compatibility with the
   multiple dihedral scheme in CHARMM22.
  
   Please send all remarks and suggestions to the CHARMM web page at
   www.charmm.org, Parameter Set Discussion Forum

   ADM Jr., July 2008
   www.charmm.org or
   http://www.pharmacy.umaryland.edu/faculty/amackere/

   References

   EXTENDED ATOM NUCLEIC ACID PARAMETERS

   Nilsson, L. and Karplus,M. Empirical Energy Functions for Energy
   Minimizations and Dynamics of Nucleic Acids J. Comp.  Chem.
   7:591-616, 1986

   PARAM19 PROTEIN PARAMETERS

   Reiher, III., W.E. Theoretical Studies of Hydrogen Bonding, Ph.D.
   Thesis, Department of Chemistry, Harvard University, Cambridge, MA,
   USA, 1985

   and

   Neria, E., Fischer, S., and Karplus, M.  Simulation of Activation Free
   Energies in Molecular Systems, Journal of Chemical Physics, 1996, 105:
   1902-21. 

   CHARMM22 PROTEIN PARAMETERS

   MacKerell, J., A.D.; Bashford, D.; Bellott, M.; Dunbrack Jr., R. L.;
   Evanseck, J.; Field, M. J.; Fischer, S.; Gao, J.; Guo, H.; Ha, S.;
   Joseph, D.; Kuchnir, L.; Kuczera, K.; Lau, F. T. K.; Mattos, C.;
   Michnick, S.; Ngo, T.; Nguyen, D. T.; Prodhom, B.; Reiher, I., W. E.;
   Roux, B.; Schlenkrich, M.; Smith, J.; Stote, R.; Straub, J.; Watanabe,
   M.; Wiorkiewicz-Kuczera, J.; Yin, D.; Karplus, M.  All-hydrogen
   Empirical Potential for Molecular Modeling and Dynamics Studies of
   Proteins using the CHARMM22 Force Field.  Journal of Physical
   Chemistry B, 1998, 102, 3586-3616.

   CMAP 2D correction surface

   MacKerell, A.D., Jr,. Feig, M., Brooks, C.L., III, Extending the
   treatment of backbone energetics in protein force fields: limitations
   of gas-phase quantum mechanics in reproducing protein conformational
   distributions in molecular dynamics simulations, Journal of
   Computational Chemistry, 25: 1400-1415, 2004.

   FOR PHOSPHOTYROSINE

   Feng, M.-H., Philippopoulos, M., MacKerell, Jr., A.D., and Lim, C.
   Structural Characterization of the Phosphotyrosine Binding Region of a
   High-Affinity SH2 Domain-Phosphopeptide Complex by Molecular Dynamics
   Simulation and Chemical Shift Calculations, Journal of the American
   Chemical Society, 1996, 118: 11265-11277

   CHARMM27 NUCLEIC ACID PARAMETERS

   Foloppe, N. and MacKerell, Jr., A.D. "All-Atom Empirical Force Field for
   Nucleic Acids: 2) Parameter Optimization Based on Small Molecule and
   Condensed Phase Macromolecular Target Data.   2000, 21: 86-104.

   MacKerell, Jr., A.D. and Banavali, N. "All-Atom Empirical Force Field for
   Nucleic Acids: 2) Application to Molecular Dynamics Simulations of DNA
   and RNA in Solution.  2000, 21: 105-120.

   CHARMM27/32 LIPID PARAMETERS

   Feller, S. and MacKerell, Jr., A.D. An Improved Empirical Potential
   Energy Function for  Molecular Simulations of Phospholipids, Journal
   of Physical Chemistry B, 2000, 104: 7510-7515.

   Feller, S. and MacKerell, Jr., A.D. An Improved Empirical Potential
   Energy Function for  Molecular Simulations of Phospholipids, Journal
   of Physical Chemistry B, 2000, 104: 7510-7515.

   Klauda, J.B., Brooks, B.R., MacKerell, A.D., Jr., Richard M. Venable,
   R.M. and Pastor, R.W., An Ab Initio Study on the Torsional Surface of
   Alkanes and its Effect on Molecular Simulations of Alkanes and a DPPC
   Bilayer, Journal of Physical Chemistry B, 109; 5300- 5311, 2005

   CHARMM36 LIPID PARAMETERS

   !above references and personal communication from the following
   Jeffery Klauda
   Joseph O'Connor
   Marcus Hadle
   Richard Venable
   J. Alfredo Freites
   Douglas Tobias
   Carlos Mondragon-Ramirez
   Igor Vorobyov
   Alexander D. MacKerell, Jr
   Richard W. Pastor

   POLYUNSATURATED LIPIDS

   Feller, S.E., Gawrisch, K. and MacKerell, Jr., A.D. "Polyunsaturated
   Fatty Acids in Lipid Bilayers: Intrinsic and Environmental
   Contributions to their Unique Physical Properties,: Journal of the
   American Chemical Society, 2002, 124:318-326

   NAD+, NADH and PPI

   Pavelites, J.J., Bash, P.A., Gao, J. and MacKerell, Jr., A.D. A
   Molecular Mechanics Force Field for NAD+, NADH, and the Pyrophosphate
   Groups of Nucleotides, Journal of Computational Chemistry, 1997, 18:
   221-239.

   CHARMM22 NUCLEIC ACID PARAMETERS

   MacKerell Jr., A.D., Wiorkiewicz-Kuczera, J. and Karplus, M. An
   all-atom empirical energy function for the simulation of nucleic
   acids, Journal of the American Chemical Society, 1995,
   117:11946-11975. 

   CHARMM30 PROTEIN FLUCTUATING CHARGE POLARIZABLE MODEL

   Patel, S., MacKerell, A.D., Jr., Brooks, C.L., III, CHARMM fluctuating
   charge force field for proteins: II Protein/solvent properties from
   molecular dynamics simulations using a nonadditive electrostatic
   model, Journal of Computational Chemistry, 25: 1504-1514, 2004

   CHARMM35 ether parameters

   Vorobyov, I., Anisimov, V.M., Greene, S., Venable, R.M., Moser, A.,
   Pastor, R.W., and MacKerell, A.D., Jr. "Additive and Classical Drude
   Polarizable Force Fields for Linear and Cyclic Ethers," Journal of
   Chemical Theory and Computing, 3: 1120-1133, 2007

   ! O-C-C-O torsion modified

   Hwankyu Lee, Richard M Venable, Alexander D MacKerell Jr., Richard W Pastor
   Molecular dynamics studies of polyethylene oxide and polyethylene glycol:
   Hydrodynamic radius and shape anisotropy
   Biophysical J., 95: 1590-1599, 2008

   CHARMM35 Carbohydrate Additive Parameters

   ! pyranose monosaccharides
   Guvench, O., Greene, S.N., Kamath, G., Brady, J.W., Venable, R.M.,
   Pastor, R.W., MacKerell, Jr., A.D. “Additive empirical force field for
   hexopyranose monosaccharides” Journal of Computational Chemistry, 29:
   2543-2564, 2008. PMID: 18470966

   ! linear sugars, sugar alcohols, and inositol
   Hatcher, E., Guvench, O., and MacKerell, Jr., A.D. “CHARMM Additive
   All-Atom Force Field for Acyclic Polyalcohols, Acyclic Carbohydrates
   and Inositol,” Journal of Chemical Theory and Computation, 5:
   1315-1327, 2009, DOI: 10.1021/ct9000608.

   ! hexopyranose glycosidic linkages
   Guvench, O., Hatcher, E. R., Venable, R. M., Pastor, R. W., MacKerell,
   A. D. Jr. “Additive Empirical CHARMM Force Field for glycosyl linked
   hexopyranoses,” Journal of Chemical Theory and Computation, 5,
   2353–2370, 2009, DOI: 10.1021/ct900242e

   ! furanose monosaccharides
   Hatcher, E. R.; Guvench, O.; MacKerell, Hatcher, E., Guvench, O., and
   MacKerell, Jr., A.D. “CHARMM Additive All-Atom Force Field for
   Aldopentofuranose Carbohydrates and Fructofuranose.” Journal of
   Physical Chemistry B. 113:12466-76, 2009, PMID: 19694450

   !CHARMM GENERAL FORCE FIELD, CGenFF
   Vanommeslaeghe, K., Hatcher, E.,Acharya, C., Kundu, S., Zhong, S.,
   Shim, J., Darian, E., Guvench, O., Lopes, P., Vorobyov, I. and
   Mackerell Jr., A.D., J. Comput. Chem., DOI: 10.1002/jcc.21367


