.. py:module:: zerom

========
Z Module
========

The Z Module is a general facility for carrying out conformational
searches based on rigid-geometry mapping or so-called "zero-order" 
minimization.  It also includes 1st-order minimization methods (Steepest
Descent and Conjugant Gradient), but the fundamental structure of the method
is nonetheless grid-based.

Reference: 

RJ Petrella, M Karplus. "A Versatile Deterministic Method for Conformational
Searches:  Application to CheY" (to be published).

.. _zerom_method:

The Z Method
------------

The Z method depends on a partitioning of the conformational space of
a system into subsystems or "subspaces", where a subspace is a subset of all
the degrees of freedom in the system.  Once defined, the subspaces can be 
searched independently or in combinations.  Specifically, the facility searches
all N!/(N-n)!n! combinations of N subspaces taken n at a time, where N and n 
are any positive integers, and n<=N. Each grid or combination searches a
set of n subspaces having C(1),C(2),...C(n) conformations across a total of 
C(1)*C(2)*...C(n) conformations. For example, if there are 4 subspaces, each 
containing 20 conformers, taking 2 subspaces at a time would produce 6 grids
(4!/2!2!), each with 400 grid-points or energy calculations.  The overall 
search results in a subspace, the "product subspace", that is the union of all 
the starting or "reactant" subspaces. Hence, when n > 1, the product subspace 
is always larger than each of the reactant subspaces.  For each grid, the 
subspaces not being searched are held fixed at their initial conformations 
(internal coordinate values occurring in the IC table prior to any of the 
grid searches having been done).  The conformers resulting from the search 
may be written to a .dcd ("trajectory") file or to an output conformer file. 
The latter is formatted in such a way as to allow it to be used as an input 
conformer file in subsequent searches.

In the current implementation, the degrees of freedom are internal
coordinates, although they could, in principle, include Cartesian coordinates
as well. The degrees of freedom--proper or improper dihedral angles, bond
angles, or bond lengths--must be defined.  Then, conformers for each subspace
are read from one or more input conformer files, which contain lists of the
degrees of freedom, their possible values, and corresponding subspace numbers.
Hence, the input conformer files define the subspaces as well as specifying
their possible conformations.

Before beginning the search, the N subspaces to be searched must be
"loaded."  This loading feature allows for searching various parts of the
system without having to redefine substructures or change the numbering of
subspaces--i.e. all the substructures/subspaces may be defined once and then
subsets can be selected for inclusion in the searches by use of the loading
feature.  

A subspace has associated with it an atomic "substructure", which is 
comprised of the atoms that are to be included in the energy calculations.
For each overall search, the set of atoms included in the calculation is the
union of the sets of atoms in all of the loaded subspaces.  Hence, for a
given subspace combination, or grid, within the larger search, atoms
corresponding to subspaces not being searched are also included in the
calculation. Note that an atom can be (and usually is) a part of multiple
substructures, but in general, a degree of freedom should not be a part
of multiple subspaces.

Subspaces can be searched in multiple regions.  For the sake of 
simplicity, this is done essentially by replicating subspaces, rather than by 
the use of a separate "region" layer in the hierarchical formalism of the 
method. This means that rather than assigning a subspace two search regions, 
the subspace is simply replicated and the two resulting subspaces are then 
assigned to different categories or types.  A search can be done over all 
N!/M!m!(N-M-m)! combinations of the N subspaces, taken M + m at a time, where 
M is the number of subspaces being searched in the first region and m is the 
number of subspaces being searched in the 2nd (or "minor") region.  This is 
discussed separately below under Supplementary.

In an attempt to optimize its utility, the Z module is designed to be 
able to be entered and exited easily and often over the course of a single
charmm script. This allows sections of Z module commands in CHARMM scripts 
to be conveniently interspersed with other CHARMM commands and functions.
Note, however, that some standard miscellaneous CHARMM commands, such as 
the DEFIne, SET, and CALC commands, and the OPENing of units for writing 
output, are supported within the Z module, itself.

Distance and degree-of-freedom constraints may be defined, and then
selected ones may be used (loaded) in a given search. For each new 
conformation, the loaded distance constraints (if they exist) are checked 
against the structure. If the conformer satisfies all of the distance 
constraints, then the loaded degree-of-freedom constraints are checked.  If 
the degree of freedom constraints are all satisfied, then the energy 
calculation on that conformer is done.

Z module searches result in changes in the main coordinate set as
well as in the main IC table.

Other Notes:

- At the start of the Z module run, the ZMEMory command must be given 
  so as to allocate the memory for the run.
- The Z module uses the CHARMM internal coordinate tables.  Hence, in 
  order for the module to function, 1) an IC table must be present, having
  either been generated (GENErate ...  SETup command) or read in
  (READ IC command); and 2) the IC FILL command must be given.
  This command fills the IC table with internal coordinate values derived
  from the Cartesian coordinates in the main coordinate set (see also
  intcor.doc).
- The Z module requires the use of the BYCC list-builder (BYCC keyword in
  the NBONDs or UPDAte command).
- Note that the Z module modifies the active atom lists, so that the
  NBACtive and BACTive commands (which specify the atoms to be considered
  in the non-bonded and bonded energies, respectively) should be invoked
  after exiting the module and before calculating CHARMM energies or forces
  with non-Z-module commands. The same is true for the fixed atom list:
  the CONS FIX command should be invoked if the CFIX keyword has been used
  in the Z module searches (i.e. in the ZSEArch command).
- For the Z module to function properly, The CHARMM code must be compiled
  with the ACTBOND and ZEROM keywords in the precompiler keyword list file 
  (build/xxx/pref.dat).  The NOBYCC and GENETIC precompiler keywords cannot be
  used with the Z module (ZEROM).  The GENETIC and ACTBOND keywords are
  incompatible.
- The Z module currently supports CMAP, but not for SD or CONJ minimizations.
  It does not support the implicit solvation models, (ACE, GB, etc.) except
  for EEF1.  
  It is not currently parallel.


.. _zerom_syntax:

Syntax of the ZMOD command
--------------------------

::

   ZMOD        enter the Z module
   END         exit the Z module


   Subcommands Keywords

   ZMEMory     NDOFs  {int} NSUBspaces {int} NCONFormers {int} NVALues {int}  
               NATOm  {int}  CSIZe  {int} NDIStances NDAToms NDFConstraints
 
   ZDEFine     DOFR  {int} internal coordinate specification
        
                                       {ATOM SELECTION X 2}  REVErse
                                     /
                  SUBStructure {int}             or
                                     \
                                        ALIAs {int}

               DIST {int} {ATOM SELECTION X 2}  GTHAn {REAL}  LTHAn {REAL}
 
               DFCO {int} DOF {int}  GTHAn {REAL}  LTHAn {REAL} 

   ZREAd       CONF RUNI {int} WUNI {int} MISSing REDUndant SSDIsorder

   ZLOAd  NULL SUBS {int1 int2 int3 ... }  CLEAR
               DIST {int1 int2 int3 ... }  CLEAR 
               DFCO {int1 int2 int3 ... }  CLEAR

   ZSPECifications 

   ZSEArch     TAKE {int} WRUN {int} ECUT {real} VCUT {real} 
               WCOM MINA SPECifics CFIX MSD NOOR WDUN NDST
               SD {int} CONJ {int} 

   ZMIN        LOAD


.. _zerom_description:

ZMOD Subcommand Description
---------------------------

*  ZMOD
*  END

   ::
  
      ZMOD 
      END

   Enter (ZMOD) and exit (END) the Z module.


*  ZMEMory 

   Allocates memory for the entire module. It and several keywords must be
   specified before any other ZMOD subcommands are given. NDOFs is the number of 
   degrees of freedom to be defined in the system; NSUBspaces is the number of
   subspaces/substructures to be defined; NCONFormers is the total number of
   conformers in all subspaces to be read in via conformer files; NVALues is the 
   cumulative length of all conformer files (total number of degree of freedom
   values, one value per line) to be read in; NATOM is the total number of atoms 
   in all substructures (since atoms can occur in multiple substructures, this 
   number may exceed the total number of atoms in the system); CSIZe is the 
   maximum number of degrees of freedom that occur in any one conformer 
   (optional). NDIST is the number of distance constraints (required only
   if distance constraints are being used). NDAToms is the total number of atoms 
   used in all of the distance constraints (required only for distance 
   constraints). NDFConstraints is the number of degree-of-freedom constraints 
   (required only for use of DOF constraints  the degrees of freedom (DOFR), 
   the substructures (SUBS), distance constraints (DIST), and internal coordinate
   constraints (DFCO) which must all be numbered serially from 1.  The degrees
   of freedom (DOFR) syntax is similar to that for editing the IC table
   (see IC EDIT), and involves either bond distances (BOND), bond angles (ANGL),
   or improper or proper dihedral angles (DIHE).  The IC table should have been
   filled "IC FILL" prior to the ZDEFine DOFR command.  Importantly, when a degree
   of freedom (DOFR keyword) is thus defined, its initial value from the IC table
   is stored and used as a default value when the subspace containing that degree
   of freedom is included among the N subspaces in a search but is not one of
   the n subspaces being searched on the current grid. 
   Although the degree-of-freedom definitions can be changed easily with the 
   ZDEFine DOFR command, for convenience, the design of the Z module is 
   intended to allow the user to define all the necessary degrees of freedom 
   initially, and then to perform all searches on the basis of that initial set of
   definitions.


*  ZDEFine

   The substructure definition contains two atom selections: the first
   defines the atoms that make up the substructure, and the second defines 
   the moving atoms--i.e. the atoms whose positions are to be initialized and 
   rebuilt during the searches. By default, the Cartesian coordinates are built 
   from the internal coordinates in the order of the atom numbers in the 
   substructure, lowest to highest.  The REVErse keyword causes the substructure 
   to be built up in the opposite sense, from higher atom number to lower.  This 
   reversal is applied only to those substructures to which the keyword is 
   assigned (and their aliases).  The "REVErse" keyword is often necessary, for 
   example, in building up loop structures from N-terminal and C-terminal halves.

   (The alternative to explicitly defining a substructure is to define it
   as an alias of another one by use of the keyword "ALIAs {int}."
   This is useful in combination with the ZCATegorize MINOR and ZSEARCh MINOR
   commands to search the same subspace in two different regions.
   See Searching Subspaces in Two Regions, below.)

   The distance constraint definitions also contain two atom selections
   each. Each constraint is applied to the distance between the centers of
   geometry of the two selections.  The distance constraint number
   (serially ordered from 1), two atom selections, and upper and lower distance
   bounds must be given for each distance constraint.  

   The internal coordinate or DOF constraint (DFCO) definitions require the 
   specification of the DOF constraint number as well as the integer corresponding
   to the actual degree of freedom as defined in the ZDEFine DOF commands.  The 
   upper and lower bounds for the constraints (in degrees) must also be specified.
   The constraints are currently implemented for dihedral, improper, and bond 
   angles (not bond lengths).  The bounds must be in the interval [-360,360]. This
   feature is intended for applying constraints to degrees of freedom that are not
   explicitly being searched.  A good example occurs in loop searches, in which 
   both halves of the loop are being searched together for low-energy 
   combinations.  The dihedral between the two searched halves of the loop is not 
   explicitly being searched (it doesn't exist in either half of the loop), but it
   may need to be constrained in a simultaneous search of both halves.

   ZEDEfine must be invoked once for each DOFR, SUBS, DIST definition.

   Example 1:

   ::
   
      ZDEFINE DOFR 11 DIHE 44 CA 44 CB 44 CG 44 OD1

   This defines the 11th  degree of freedom as a dihedral involving atoms CA, 
   CB, CG, and OD1 of residue 44.

   Example 2:

   ::
   
      ZDEFine SUBStructure 2 sele resi 20 .around. 15  end -
         SELE RESI 20 END

   This defines the second substructure as containing all atoms within 15 
   angstroms of residue 20 and specifies that residue 20 will be rebuilt 
   in the searches.

   Example 3:

   ::
   
      ZDEFine DFCO 1 DOF 11 GTHAn -100 LTHAn 100 

   This defines the first internal coordinate constraint as being
   applied to degree of freedom 11 (the dihedral angle as defined above),
   with bounds of [-100,100] degrees.


*  ZREAd 

   Reads in input conformer files.  CONF reads in a conformer file from unit 
   RUNI and writes the "compressed"  conformer file to unit WUNI. The conformer 
   file is compressed in the sense that there is no redundant information when the
   file is read from start to finish--i.e. only the degrees of freedom that change
   from one conformer to the next are written.  This form of the file, which is 
   the form stored internally and used by the module, decreases space and memory 
   requirements and increases speed.  ZREAd also sorts each conformer by degree of
   freedom internally. ZREAd expects the first conformer of each subspace to be 
   "complete"--i.e. all its degrees of freedom must be present.  Subsequent 
   conformers in the subspace may be present in compressed form.
   
   The format of the conformer file is as follows (fortran format 
   I14,I14,I14,F14.7):
   
   ::

             1             1             7    18.4857895  
             1             1             3   180.7483736  
             1             1             1    50.1098635  
             1             2             7    24.2836786    
             1             2             3   180.7483736    
             1             2             1    50.1098635    
             1             3             7    30.09877     
             1             3             3   180.168350     
             1             3             1    60.87635
             2             4             8    64.27840  
             2             4             9    80.287635  
             2             4            11    17.981386  
             2             5             8    56.976233    
             2             5             9    80.287635
             2             5            11     1.8923759    

   The first column gives the subspace number, the second column gives the 
   conformer number, the third column gives the degree of freedom, and the fourth 
   column gives the numerical value taken by the degree of freedom.  This example 
   file has not been compressed.  The compressed file would appear as follows:

   ::
   
             1             1             1    50.1098635
             1             1             3   180.7483736
             1             1             7    18.4857895
             1             2             7    24.2836786
             1             3             1    60.8763500
             1             3             3   180.1683500
             1             3             7    30.0987700
             2             4             8    64.2784000
             2             4             9    80.2876350
             2             4            11    17.9813860
             2             5             8    56.9762330
             2             5            11     1.8923759

   Note that parts of conformers 2 and 5 have been omitted in the compressed file 
   because they contain redundant information. 

   By default, ZREAd expects the input file to have perfectly ordered (and 
   consecutive) conformer and subspace numbers, beginning with 1. This behavior 
   can be overridden with the MISSing and SSDisorder keywords, respectively.  
   MISSing allows for "gaps" in conformer numbers (they still need to be ordered),
   and SSDIsorder allows for disorder in the subspace numbers.  This means that, 
   for example, subspace 2 can occur before subspace 1 in the file.  However,
   each subspace must be made up of consecutive input lines--i.e. subspaces cannot
   be "split" in a given conformer file.  ZREAd also expects degrees of freedom 
   not to occur in multiple subspaces.  This can be overridden with the keyword 
   REDUndant, which is necessary when subspaces overlap or when the same subspace
   is defined twice (either with or without aliasing).  If two consecutive 
   conformers in a conformer file are exactly the same, ZREAd will complain, but 
   it does not check the entire list for conformer redundancies.  ZREAd expects 
   the first subspace in a file to be subspace 1; this again can be overridden 
   with the SSDIsorder keyword. 

   If multiple conformer files are read in, each requires a separate ZREAd 
   CONFormer command, and the subspaces of the files will be renumbered 
   internally, so that they are consecutive.  For example, if the last subspace of
   the first file is 10, the first subspace of the second file will be called 
   subspace 11, regardless of the number appearing in the subspace column.  In 
   general, if multiple conformer files are being read, the subspaces should be 
   well-ordered in all files--i.e. no "gaps" or inverted orders within each given 
   file.  The SSDIsorder keyword may not help in the case of multiple conformer 
   files with disordered subspace numbers.


*  ZLOAd 

   "Loads" subspaces or constraints to be used in searches.  This means the
   specified subspaces or constraints are selected out of the sets of all possible
   defined subspaces or constraints, and will be used in the searches.

   If the "SUBS" keyword is used, the selected subspace/substructure will be 
   used in the searches.  A subspace/substructure must be loaded in order to be 
   searched.  (Substructure/subspace aliases should also be loaded if they 
   are to be searched.)  Note that if more than one substructure/subspace
   is loaded, the "product" subspace/substructure--i.e. the one appearing
   in the output conformer file--will generally be larger than each of the
   "reactant" subspaces/substructures, since the product is the union of all
   the reactant (loaded) subspaces/substructures.

   The "NULL" keyword allows for the loading of empty subspaces--i.e. ones 
   containing no conformers or ones for which no conformer file has been read in. 
   This is usefull in cases where, for example, a single subspace corresponds to 
   (the union of) 2 previously defined substructures.  When it is used, the
   "NULL" keyword must occur as the first keyword in the ZLOAd command.

   If the "DIST" keyword is used, the selected distance constraints will be
   applied in the searches. 

   The "DFCOnstraints" keyword is analogous to the "DIST" keyword, except it
   is used for loading the internal coordinate (DOF) constraints.

   The CLEAr keyword causes the previously loaded subspaces (or distance 
   constraints) to be unloaded prior to the parsing of any other keyword in the 
   ZLOAd command.


*  ZSPEcifications

   This command causes the currently stored (user defined) search 
   specifications to be written to standard output.


*  ZSEArch

   Carries out N!/(N-n)!n! grid searches of the N loaded subspaces taken n at
   a time, where n <= N.   It eliminates conformers that are above specified 
   thresholds in energy, or outside of distance constraints, and can write 
   results to a product conformer file or a .dcd file (or both).  The output 
   conformer file is formatted exactly as the input conformer file, except that, 
   in addition to the subspace, conformer, degree-of-freedom, and d.o.f.
   numerical value information, the energies are also written in an additional
   column to the right (in 1X,F19.7 fortran format).  For convenience, the energy
   of each conformer is written next to each of its component degrees of freedom.
   If requested, the MSD's (mean square deviations from a target structure) are 
   written to a sixth column.  If writing to a .dcd file is specified, all atoms
   in the system are included.

   TAKE {int} specifies the number of subspaces to be searched in 
   combination--i.e. the dimensionality of the grid ('n' in the expression above).

   WCOM partially compresses the information in the output conformer file. 
   For each grid, only the first outputted conformer is written in its entirety; 
   the remaining conformers include only degrees of freedom that are being varied 
   in that grid.  Since the ZREAd command requires that the first conformer be 
   "complete," and since the conformers are only accessed serially from the 
   start of the file, care should be taken when modifying a compressed conformer 
   output file prior to its use as an input file for a subsequent search.  In 
   general, such modification is not recommended.

   ECUT (real} and VCUT {real} are the fixed and variable energy cutoffs, 
   below which conformers will not be written to the output conformer file. If 
   both are specified, both cutoffs are used. ECUT is an absolute energy cutoff 
   value. VCUT is a tolerance above the running energy minimum and hence varies 
   over the course of the calculation.  For example, ECUT 2000 VCUT 30 would 
   result in the exclusion of all conformers with energies > 2000 kcal/mol OR
   >30 kcal/mol above the running minimum (but not necessarily the absolute
   minimum).

   MINA assigns the structure (main coordinate set) to the global minimum
   after the search is completed.

   WRUN the unit number to which the output conformer file is to be written.

   MWRU is the write unit number to which the minimum-energy conformer file
   is to be written.

   TAG {int} is the numerical (integer) name given to the product subspace 
   and written in the first column of the output conformer file.

   CFIX keyword causes the non-bonded energy contributions of fixed-fixed 
   atom pairs not to be calculated (bonded energy contributions remain).
   Fixed atoms are defined as ones that are not specified as moving for any
   of the substructures included in a search.  I.e. the moving atoms in a search
   are taken as the union of the second atom selections of all the ZDEFINE
   SUBStructure commands corresponding to the loaded substructures.  The fixed
   atoms are all the atoms in the system that are not taken as moving.
   Since a significant fraction of the atoms included in the searches may not
   be moving, the use of CFIX may speed up the calculations without affecting
   the relative energies.  Note that CFIX will also fix the atoms during
   minimization.

   MSD calculates the MSDs (mean squared deviations, or squared RMSDs) of the
   conformers from the structure in the comparison coordinates.  The MSD's are 
   calculated for all atoms in all the loaded substructures (union of all first 
   atom selections in the ZDEFine SUBStructure command).  Least-squares-fit 
   reorientation is done prior to rmsd calculation unless the NOORientation 
   keyword is also specified.  The MSD keyword currently is not compatible with 
   WCOMpression.

   SPECifics keyword causes the user-defined specifications that will
   be used in the search to be written to standard output.

   WDUN unit for writing out structures to .dcd ("trajectory") files.

   NDSTeps total number of frames to be written to .dcd [default = 1000]

   NOEN causes the energies not to be calculated.  This is useful,
   for example, when using only distance or internal coordinate constraints,
   or when reading in a conformer file only for the purposes of writing out
   the corresponding .dcd file.

   SD {int} or CONJ {int} cause a minimization of the coordinates to be 
   carried out at each gridpoint in the search.  The coordinates of all
   atoms involved in the search (union of all loaded substructures) are minimized,
   unless constraints are set prior to the minimization (e.g. CONS FIX can be
   used as usual, outside of the Z module). Note that after
   each minimization, the unminimized coordinates are restored before the next 
   point on the grid.  This is necessary to preserve the internal coordinates that
   are not being explicitly searched but are nontheless affected by minimization. 
   The non-bonded list is always updated at the start of the minimizations, 
   provided the update frequency (INBF) is > 0.  The energies written to the 
   conformer file and the structures written to the .dcd file correspond to the 
   minimized coordinates.  The assigned minimum-energy structure (which is 
   obtained with the MINA keyword in the ZSEArch command) corresponds to the 
   unminimized coordinates. 
    
*  ZMIN  

   assigns the structure to the minimum indicated.  If LOAD keyword
   specified, it assigns the structure to the minimum over all searches done
   since the last ZLOAD command.

.. _zerom_examples:

ZMOD Command Usage Examples
---------------------------

::

   ZMOD  !enter module
   ZMEMOry NDOF 100 NSUBSPACE 20 NATOM 25000 NCONF 2000 NVAL 10000 !allocate memory

   ZDEFINE DOFR  1 DIHE 10 C 11 N 11 CA 11 C
   ZDEFINE DOFR  2 DIHE 11 N 11 CA 11 C 12 N
   ZDEFINE DOFR  3 DIHE 11 C 12 N 12 CA 12 C
   ZDEFINE DOFR  4 DIHE 12 N 12 CA 12 C 13 N
   ZDEFINE DOFR  5 DIHE 29 C 30 N 30 CA 30 C
   ZDEFINE DOFR  6 DIHE 30 N 30 CA 30 C 31 N
   ZDEFINE DOFR  7 DIHE 30 C 31 N 31 CA 31 C
   ZDEFINE DOFR  8 DIHE 31 N 31 CA 31 C 32 N
   ZDEFINE DOFR  9 DIHE 35 C 36 N 36 CA 36 C
   ZDEFINE DOFR  10 DIHE 36 N 36 CA 36 C 37 N

   ZDEFine SUBStructure 1 SELE ALL END SELE IRES 1:13 END
   ZDEFine SUBStructure 2 SELE ALL END SELE IRES 1:32 END
   ZDEFine SUBStructure 3 SELE ALL END SELE IRES 1:38 END
   ZDEFine SUBStructure 4 SELE ALL END SELE IRES 1:50 END

   open unit 12 write card name echodata
   open unit 10 read card name short.conf  !conformer file

   ZREAD REDU CONF RUNI 10 WUNI 12  !we are allowing redundancy of dofs in the 
     !conformer file

   END ! exit the Z module

   !other charmm commands here

   ENERGY
   MINIMIZE SD
   .
   .
   .

   ZMOD !reenter Z module

   ZLOAD SUBS CLEAR  2 3 4   !load 3 out of 4 defined subspaces 

   open unit 15 write card name ener.file
   !specify searches that take 2 subspaces at a time
   !write to unit 15, assign the minimum to the main coordinate
   !set after the search, set product subspace tag to one (for output),
   !save (and write) only the conformers within 500 kcal/mol of the running
   !minimum energy, write search specifications to standard output

   ZSEArch TAKE 2 WRUN 15 MINAssign TAG 1 VCUT 500 SPEC  
   close unit 15

   END !exit Z module

   !other charmm input:
   COOR RMS ORIENT SELE ALL END

   OPEN UNIT 10 WRITE CARD NAME crd/final.crd
   WRITE COOR CARD UNIT 10
   CLOSE UNIT 10

.. _zerom_supplementary:

More Complex Searches: Searching Subspaces in Two Regions
---------------------------------------------------------

In some cases, the user may want to be able to search the same 
conformational subspace in two different regions in a series of grids.  Take, 
for example, the case of 3 rigid helices for which the optimal packing 
arrangement is sought.  The most efficient search algorithm may be one in 
which each helix is allowed to sample a large region of conformational space, 
while the other two helices make local adjustments.  Hence, each helix would 
have two defined regions: one for the large changes and one for the local 
adjustments--i.e. a "major" and a "minor" region.  As helix A was searched in 
its major region, helices B and C would be searched in their minor regions, and
likewise for the other combinations. 

Rather than introduce a "region" formalism into the structure of the 
Z module--i.e. having each subspace be composed of multiple regions--which 
would be complicated for both users and developers, the decision was made to 
allow subspaces to have aliases.  The alias of a subspace has exactly the same
degrees of freedom and atomic substructure, but it is called separately and can
be assigned different conformations from its primary subspace counterpart
through the conformer file(s).  The CATEgorize MINOr command is necessary to
categorize the aliased subspace as "minor", or secondary, so that it is treated
as such in the ZSEARch command.  The "MINOR {int}" keyword must be used in the
ZSEArch command to indicate how many minor subspaces or regions are to be 
searched at a time.

ZDEFine SUBSPACE {int} ALIAS {int} is used as an alternative to explicitly
defining an atomic substructure. The integer after ALIAS specifies the number
of the primary subspace to which the alias corresponds. Aliasing is necessary
when searching two regions of the same subspace (i.e. a subspace and its alias)
in a set of searches, because it tells the ZSEArch command that the two loaded
subspaces correspond. This allows the Z module to exclude search grids in
which the same subspace effectivley appears more than once. Without aliasing,
subspaces are combined in searches without regard to whether they are the same
or not, so that results may be redundant or incomplete. The atom selections
and REVErsal specifications are not necessary (and not allowed) when ALIAsing
a subspace.

Syntax for additional commands in multiple-region searches
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

   ZCATegorize MINOr CLEAR NONE {int1 int2 int3 ... }

   ZSEArch     MINO {int} RESS


Description of additional commands
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*  ZCATegorize 

   Categorizes a subspace, effectively creating a different search region for
   the space.  Currently, only one alternative region can be searched for a given
   subspace, which is assigned the category (keyword) "MINOr."  The keyword 
   "minor" indicates that the subspace is to be included in the set of "MINOr" 
   or secondary regions that are to be searched by the ZSEARch command (it implies
   nothing about size).  A set of integers must be specified which refer to the 
   numbers of the subspaces being categorized as minor.  The CLEAR or NONE 
   keywords clear any previous categorizations of subspaces as minor.


*  ZSEARCH MINOR {int} RESS

   MINOr {int} specifies the number of subspaces, out of TAKE searched
   subspaces, that should be searched in their "minor" or alternative regions
   [default=0].,

   The RESS ("REpeat SubSpace") keyword allows multiple versions of the same
   subspace, as defined by the ZDEFine ALIAS command, to be searched at the same
   time (i.e. on the same search grid).  This is not generally recommended and the
   user employs this keyword at his own risk.

Example 
^^^^^^^

An example of a search using multiple subspace regions, modified from 
an example above

::

   ZMOD  !enter module
   ZMEMOry NDOF 100 NSUBSPACE 20 NATOM 25000 NCONF 2000 NVAL 10000 !allocate memory

   ZDEFINE DOFR  1 DIHE 10 C 11 N 11 CA 11 C
   ZDEFINE DOFR  2 DIHE 11 N 11 CA 11 C 12 N
   ZDEFINE DOFR  3 DIHE 11 C 12 N 12 CA 12 C
   ZDEFINE DOFR  4 DIHE 12 N 12 CA 12 C 13 N
   ZDEFINE DOFR  5 DIHE 29 C 30 N 30 CA 30 C
   ZDEFINE DOFR  6 DIHE 30 N 30 CA 30 C 31 N
   ZDEFINE DOFR  7 DIHE 30 C 31 N 31 CA 31 C
   ZDEFINE DOFR  8 DIHE 31 N 31 CA 31 C 32 N
   ZDEFINE DOFR  9 DIHE 35 C 36 N 36 CA 36 C
   ZDEFINE DOFR  10 DIHE 36 N 36 CA 36 C 37 N


   ZDEFine SUBStructure 1 SELE ALL END SELE IRES 1:13 END
   ZDEFine SUBStructure 2 SELE ALL END SELE IRES 1:32 END
   ZDEFine SUBStructure 3 SELE ALL END SELE IRES 1:38 END
   ZDEFine SUBStructure 4 SELE ALL END SELE IRES 1:50 END
   ZDEFine SUBStructure 5 ALIAs 1               ! 5 is same as 1 !***************
   ZDEFine SUBStructure 6 ALIAs 2               ! 6 is same as 2 !***************

   ZCAT CLEAR MINOR  5 6  ! categorizing 5 and 6 as minor subspaces (regions) !***

   open unit 12 write card name echodata
   open unit 10 read card name short.2.conf  !conformer file, now contains 
     !conformers for subspaces 5 and 6 also

   ZREAD REDU CONF RUNI 10 WUNI 12  ! allowing redundancy of dofs in conformer file

   END ! exit the Z module

   !other charmm input

   ENERGY
   MINIMIZE SD
   .
   .
   .

   ZMOD !reenter Z module

   ZLOAD SUBS CLEAR  2 3 4 5 6   !load 5 out of 6 defined subspaces !***********

   open unit 15 write card name ener.file
   !specify searches that take 3 subspaces at a time, one of which is searched
   !in its minor region,  write to unit 15, assign the minimum to the main 
   !coordinate set after the search, set product subspace tag to one (for output),
   !save (and write) only the conformers within 500 kcal/mol of the running
   !minimum energy, write search specifications to standard output

   !include minor keyword *****************
   ZSEArch TAKE 3 MINO 1 WRUN 15 MINAssign TAG 1 VCUT 500 SPEC  !***************
   close unit 15

   END !exit Z module

   !other charmm input:
   COOR RMS ORIENT SELE ALL END

   OPEN UNIT 10 WRITE CARD NAME crd/final.crd
   WRITE COOR CARD UNIT 10
   CLOSE UNIT 10
