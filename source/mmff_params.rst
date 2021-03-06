.. py::module: mmff_params

==============================================
The MMFF94 Setup Procedure And Parameter Files
==============================================

.. _mmff_params_mmffsymb:

1. MMFFSYMB.PAR

   Starting from the input atomic species, connectivity, 
   and formal bond orders (for aromatic systems, for example, a Kekule 
   structure having alternating single and double bonds must be supplied), the 
   MMFF structural perception code automatically "sets up" the calculation by 
   perceiving and classifying rings, detecting aromaticity, and creating 
   appropriate lists of bond, angle and torsional interactions.  The atom typing 
   procedures (currently overseen by subroutines XTYPE, HTYPE and RGTYPE) 
   then assign a 4-character symbolic atom type to each atom.  Finally, the 
   entries in MMFFSYMB.PAR are used to translate the symbolic atom types into 
   the numeric atom types actually employed in assigning force-field 
   parameters to the force-field interaction terms.  Note that this is a many-
   to-one translation, in that more than one symbolic atom type usually 
   corresponds to a given numeric atom type.  The descriptive text strings to 
   the right of the columns bearing the symbolic and numeric atom types are 
   given for informational purposes only.

.. _mmff_params_mmffarom:

2. MMFFAROM.PAR

   The assignment of symbolic atom types is actually a 
   two-stage process.  For the carbons in benzene, for example, the first stage, 
   which uses only the local connectivity, would assign a symbolic atom type of 
   CE, meaning unsaturated or olefinic carbon.  The second stage (currently 
   controlled by subroutine RGTYPE) perceives that these carbons lie in an 
   aromatic six-membered ring and uses the elements of MMFFAROM.PAR to 
   assign the correct aromatic atom type.  Some lines from the version of this 
   parameter file as of 2/96 are listed below:

   ::
   
      *
      *     Copyright (c) Merck and Co., Inc., 1994, 1995, 1996
      *                    All Rights Reserved
      *
      * MMFF AROMATIC SYMBOLIC TYPES VS ORIGINALLY ASSIGNED TYPES 
      * OLD	AROM	AT    RING	       IM      N5
      * TYPE	TYPE	NUM   SIZE     L5      CAT    ANION
         C*	 CB	6	6	0	0	0     
         N*	 NPYD	7	6	0	0	0     
         NCN+  NPD+	7	6	0	0	0     
         N+=C  NPD+	7	6	0	0	0     
         N=+N	 NPD+	7	6	0	0	0     
         N2OX	 NPOX	7	6	0	0	0     
         *
         C*	 C5A	6	5	2	0	0
         C*	 C5B	6	5	3	0	0
         C*	 C5	6	5	4	0	0
         N*	 NPYL	7	5	1	0	0
         N*	 N5A	7	5	2	0	0
         N*	 N5B	7	5	3	0	0
         N*	 N5	7	5	4	0	0
         CNN+	 CIM+	6	5	2	1	0
         CNN+	 CIM+	6	5	3	1	0
         CNN+	 CIM+	6	5	4	1	0


   As in all the MMFF parameter files, lines beginning with a "*" are comment 
   lines.  The column labeled "L5" refers, in the case of 5-ring systems, to the 
   position of the atom in question relative to the unique pi-lone-pair 
   containing heteroatom (which itself occupies position "1"); a "4" is an 
   artificial entry which is assigned when no such unique heteroatom exists, as 
   for example occurs in imidazolium cations and in tetrazole anions.  An entry 
   of "1" in the "IM CAT" or "N5 ANION" column must also be matched for such 
   ionic species to convert the "OLD" (first stage) to "AROM" (second stage; 
   aromatic) symbolic atom type.  Note: in matching the "OLD" symbolic atom 
   types, an "exact" match is first attempted.  If this match fails, a wild-
   carded match, using for example "C*" is then employed.

.. _mmff_params_mmffhdef:

3. MMFFHDEF.PAR

   The symbolic atom type for a hydrogen atom is 
   assigned by reference to that for its parent non-hydrogen atom, using the 
   entries in this parameter file.  Note that this approach assumes that a given 
   hydrogen is bonded to a single parent and thus disallows bridging hydrogens.  
   Note also that the code treats molecular hydrogen as a special case, in as 
   much as there is then no "parent" atom.  

.. _mmff_params_mmffdef:

4. MMFFDEF.PAR

   The procedures outlined above assign of both a symbolic 
   and a numeric atom type for each atom.  For the purpose of matching force 
   field interaction terms (such as the stretching interaction between atoms i 
   and j), the code first uses the so-generated primary atom type in seeking a 
   match with an entry in the relevant (e.g., MMFFBOND.PAR) file.  In this 
   process, the parameter files, which are kept in "canonical order" based on 
   indices derived from the numerical atom types, are processed using a rapid 
   binary search algorithm.  If present, the fully qualified parameter 
   corresponding to the precise set of atom types (supplemented, in ambiguous 
   situations, by defined bond, angle, stretch-bend, or torsion interaction types 
   as described below) is used.  For vdW, bond stretching, stretch-bend 
   interaction, and bond charge increment parameters, no equivalences are 
   recognized.  For angle bending, out-of-plane bending, and torsion 
   interactions, however, whenever the fully qualified parameter is not found, 
   the code executes a staged "step down" procedure in which increasingly 
   generic values are sought.  This protocol is governed by the entries in 
   MMFFDEF.PAR, where the primary ("Level 1") atom types define the fully 
   qualified parameters.  Entries from Levels 2 - 5 are employed as needed in 
   subsequent searches.  Those at Level 5, always "0", serve as wild cards.  Such 
   wild card values are used only for "wing atoms" in an angle bending or 
   torsional interaction or for non-central atoms in an out-of-plane interaction.
   Level 4 generally corresponds to the atomic species, Level 3 to atomic 
   species plus hybridization.  Currently, the first two levels employ identical 
   numerical atom types; additional sets of unique values for Level 1 may be 
   specified later if we find that certain atom types need to be defined more 
   specifically.  The protocol used in the step-down procedure depends on the 
   interaction type (angle, torsion, out-of-plane)  Thus, for bending of the
   i-j-k angle, a five-stage process based in the level combinations 1-1-1,
   2-2-2, 3-2-3, 4-2-4 and 5-2-5 is used.  For i-j-k-l torsion interactions,
   a five-stage process based on level combinations 1-1-1-1, 2-2-2-2, 3-2-2-5,
   5-2-2-3, and 5-2-2-5 is used, where stages 3 and 4 correspond to "half default"
   or "half wild card" entries.  For out-of-plane bending ijk;l, where j is the 
   central atom, the five-stage protocol 1-1-1;1, 2-2-2;2, 3-2-3;3, 4-2-4;4, 
   5-2-5;5 is used.  The final stage provides wild-card defaults for all except 
   the central atom.  Finally, if no parameter is found, a default rule may be 
   invoked (see below).  Note that the parameter lines beginning with a "*" are 
   comment lines and are not used, and that the symbolic atom types in the first 
   column and the abbreviated text-string definitions in the right-most portion 
   of the parameter line are for informational purposes only.

.. _mmff_params_mmffprop:

5. MMFFPROP.PAR

   This file defines atomic properties associated with the MMFF
   numeric atom types.  As described below, these properties are used in setting
   up the MMFF energy expression and in determining values for missing parameters
   from defined empirical rules.  The first two columns list the MMFF numeric atom
   type and the atomic number for that atom type. Column 3 ("crd) specifies the
   requisite number of bonded neighbors, and column 4, labeled "val", gives the
   total number of bonds made to that atom type.  Where an entry shows a two-digit
   number, such as "34" for atom types 54 and 55, the valency may be either of the
   two co-joined values -- e.g., either 3 or 4. Entries of "1" in column 5
   ("pilp") define those atom types which have a "pi lone pair" capable of
   partaking in resonance interactions with, say, an adjacent multiple bond. 
   Column 6, headed "mltb", specifies cases in which double ("2") or triple ("3")
   bonds are expected to be made to an atom having the listed atom type.  This
   information is used, for example, in the empirical rule for reference bond
   lengths discussed below in connection with the MMFFBOND.PAR file.  Entries of
   "1" in this column designate cases in which intermediate hybdrization and
   valence based bond shortenings are employed in that empirical rule, and roughly
   correspond to cases in which strong delocalization to a neighboring atom can be
   expected.  Column 7 then specifies, through an entry of "1", the set of
   "aromatic" atom types, and column 8 ("lin") in the same manner delineates those
   (necessarily dicoordinate) atom types which involve linear bond angles, i.e.,
   bond angles near 180 deg. Entries of "1" in the final column, headed, sbmb for
   "single-bond--multiple-bond", designate the set of atom types which can be
   involved in multiple bonds but alternatively can also be involved in
   delocalized, formally single, bonds with a sp2 (or sp) hybridized neighbor. 
   The latter instance, for example, would occur for the bond between the two
   central carbons of butadiene, each of numeric atom type 2.  These entries
   control the setting of the "bond type index" discussed below.  Specifically, an
   nonstandard bond type index BTIJ of "1" is assigned whenever a single bond 
   (formal bond order 1) is found a) between non-aromatic atoms i and j of types I
   and J for which "sbmb" entries of "1" appear in the MMFFPROP.PAR file or b)
   between aromatic atoms belonging to different aromatic rings (as in the case of
   the central C-C bond in biphenyl). 

.. _mmff_params_mmffbond:

6. MMFFBOND.PAR

   The form for bond stretching used in MMFF is:

   ::
   
      EBij  =  0.5*143.9325*kbIJ*Drij**2*(1 + cs*Drij + 7/12 cs**2*Drij**2),  (1)

      where: 	kbIJ =	bond stretching force constant in md/angstroms for bonded 
      		atoms i  and j of types I and J .
      	Drij =	rij - roIJ, the difference in angstroms between actual and
      		reference bond lengths between bonded atoms i and 
      		j of types I and J.
      	cs = 	-2 angs**(-1), the "cubic stretch" constant.

   Note: throughout this description, the indices i, j, k, ... represent atoms;
   I, J, K, ... denote the corresponding numerical MMFF atom types (or, 
   occasionally, the atomic species).

   A few entries taken from the MMFFBOND.PAR file are shown below:

   ::
   
       		At. Types	
            BTIJ	I	J	kbIJ    roIJ	Origin/Comment
      	-------------------------------------------------------
      	0	1	1	4.258	1.508	C94
      	0	1	2	4.539	1.482	C94
      	0	1	3	4.190	1.492	C94
      	0	1	5	4.766	1.093	C94
      	0	1	6	5.047	1.418	C94
      	0	1      10	4.664	1.436	C94
      	0	2	2	9.505	1.333	C94
      	1	2	2	4.814	1.430	#C94
      	1	2	3	4.565	1.468	C94

         C94  = CORE MMFF94 parameter - obtained from ab initio data
         X94  = EXTD MMFF94 parameter - fit to additional ab initio data
         E94  = r0 from fit to X-ray data, kb from empirical rule
         #C94 = r0 lies between C94 and E94 values, kb from empirical rule
         #X94 = r0 lies between X94 and E94 values, kb from empirical rule
         #E94  = r0 and k both from empirical rules

   The first entry, labeled BTIJ, is the "bond type" index discussed above in 
   connection with the MMFFPROP.PAR file.  As described there, this index 
   normally takes the value 0 but occasionally takes the value 1.  Thus, the 
   parameter whose three numerical indices are "0 2 2" describes a normal 
   double bond between olefinic carbons (reference bond length, 1.333 angs), 
   whereas the next-listed parameter with numerical indices "1 2 2" 
   describes a delocalized but formally single bond between two such 
   carbons, e.g., as for the central carbon-carbon bond in butadiene (reference 
   bond length, 1.430 angs).  The last column, labeled "Origin/Comment", 
   currently (as of 2/96) distinguishes core MMFF parameters derived in the 
   extended MMFF parameterization from ab initio data (labeled "C94") 
   from parameters derived from crystallographic data (labeled "E94").  
   The latter parameters take their reference bond lengths from fits to 
   experimentally determined structures, and then use the empirical rule 
   discussed below in connection with the MMFFBNDK.PAR file to assign the 
   force constant.  The "#C94" designation for the "1 2 2" parameter indicates
   that the reference bond length for this parameter reflects an intentional 
   compromise between the ab-initio derived value (e.g., for butadiene, which
   gives 1.447 angs) and (shorter) values more typically found in crystallographic
   structures (which tend to involve "push-pull" delocalization that increases
   the pi-bond order across the central, formally single C-C bond). 

  The canonical ordering reflected in the table is defined through an 
  integer index, CXB, which is computed as

  ::
  
    CXB = MC*(I*MA +J) + BTIJ,

   where BTIJ is the bond-type index, MC is an integer equal to at least the 
   maximum permissible bond-type index +1, MA is an integer equal to the 
   maximum numerical atom-type index + 1, and atom type I is less than or 
   equal to J.  Thus, the first atom type, I, changes least rapidly and the bond-
   index BTIJ changes most rapidly, such that parameters which reference the 
   same atom types but different bond types appear consecutively in the file.  
   As in the case of other force-field interaction terms discussed below, this
   canonical index is used in the binary-search algorithm used to detect the
   presence or absence of a specific parameter (i.e., a specific BTIJ, I, J
   combination).

.. _mmff_params_mmffbndk:

7. MMFFBNDK.PAR

   The entries in this file govern the empirical rule 
   used for bond-stretching force constants when a parameterized value derived 
   from fits to ab initio second derivatives is unavailable.  MMFF uses an 
   inverse sixth-power dependence1 to relate the desired force constant to the 
   force constant values tabulated in MMFFBNDK.PAR for a specific bond length, 
   usually that expected for a single bond between atoms i and j of the 
   associated atomic numbers.2  Specifically, MMFF takes 

   ::
   
      kbIJ = kbref (roref /roIJ)**6						(2)

   where kbref and roref are the force constants and associated reference bond 
   lengths given in MMFFBNDK.PAR and roIJ is the reference bond length 
   retrieved from the MMFFBOND.PAR parameter file or calculated from an
   empirical-rule procedure contained in the code (see below).  This rule was used
   to generate the force constants tabulated in MMFFBOND.PAR for the extended-
   parameterization entries labeled "E94". 

   The format of MMFFBNDK.PAR is as shown below:

   ::
   
             atomic number
      	i	j   	roref      kbref  	Source
      	----------------------------------------------
      	1	6	1.084	    5.15	C94
      	1	7	1.001	    7.35	C94
      	1	8	0.947	    9.10	C94
      	1	9	0.92	   10.6		E94
      	1      14       1.48	    2.3		E94
      	1      15	1.415	    2.95	C94
      	1      16	1.326	    4.30	C94

         C94 =  Fitted to ab-initio derived core MMFF94 force constants
         E94 =  Based on Herschberg/Laurie parameterization of Badger's rule

   Note that kbref and roref depend only on the atomic species for i and j rather 
   than on their associated MMFF numeric atom types.  The reference force 
   constants kbref labeled C94 were determined by requiring that the force
   constants kbIJ given by eq 2 closely reproduce the core, computationally
   derived force constants tabulated in MMFFBOND.PAR when the associated
   computationally derived reference bond lengths roIJ are used in the empirical
   rule -- i.e., with a rms deviation of 0.75 md/angs, as against a rms value of
   5.81 md/angs. (Badger's rule also performed reasonably well, but gave a
   somewhat higher rms deviation of 1.06 md/angs).  The force constants kbref for
   the entries labeled E94 could not be calibrated against known reference values
   in a comparable manner.  Instead, these force constants were calculated from
   the tabulated roref values using Herschbach and Laurie's parameterization of
   Badger's rule [3]_. Indeed, the empirical rule employed in MMFF would employ
   this implementation of Badger's rule directly, with roIJ given as shown below,
   should a case arise for which the atomic parameters required for eq 2 were
   unavailable. 

   Reference values not supplied in the parameter tables are assigned using 
   a modified form of the Schomaker-Stevenson rule [4]_ similar to that described
   by Blom and Haaland [5]_.  Specifically, MMFF takes 

   ::
   
      roIJ = roI + roJ - c |xI - xJ|**n - d					(3)

   where roI and roJ represent bond-order adjusted covalent radii for atom types I
   and J, xI and xJ are Pauling-scale electronegativities defined by Allred and
   Rochow and depend only on the atomic species [6]_, c is a proportionality
   constant, n is the power dependence of the electronegtive term, and d is a
   "shrinkage" factor which recognizes that the reference bond lengths in a force
   field calculation tend to be somewhat shorter than the optimized bond lengths
   for which the comparison with observed bond lengths is to be made. The
   quantities in eq (3) are provided explicitly in the program code. 

.. _mmff_params_mmffang:

8. MMFFANG.PAR

   Most angle-bending interactions are described using a 
   cubic expansion:

   ::
   
      EAijk  =  0.5 * 0.043844 * kaIJK * DTijk**2 * (1 + cb*DTijk),             (4)

      where: 	kaIJK =   angle bending force constant in md-angs/rad**2 for the
      		  angle between atoms i, j and k of atom types I, J and K.
 
   	DTijk =   Tijk - ToIJK, the difference between actual and 
   		  reference i-j-k bond angles in degrees.

   	cb = 	  -0.007 deg**(-1), the "cubic-bend" constant.


   For linear or near-linear bond angles such as those which occur in alkynes,
   nitriles, isonitriles, azides, and diazo compounds, we employ the form used 
   in DREIDING[7]_ and UFF[8]_:

   ::
   
      EAijk  =  143.9325 * kaIJK * (1 + cos(Tijk))                              (5)

   where kaIJK and Tijk are as defined above.  Angle terms i-j-k are included 
   in all cases in which atoms i and k are bonded to a common central atom j.  

   The format of MMFFANG.PAR is depicted below:

   ::
   
            ---- Atom Types ----
      ATIJK  I	J	K      kaIJK     ToIJK	   Origin/comment
      -------------------------------------------------------------------
       0	0	1	0	0.000	108.100	   0:*-1-* MMFF94 DEF
       0	1	1	1	0.851	109.608	   C94
       0	1	1	2	0.736	109.445    C94
       0	2	1	5	0.632	110.292	   C94
       0	2	1	6	1.074	108.699	   C94
       0	2	1      10	1.160	107.963	   E94
       0	3	1	3	0.974	111.746	   E94
       0	3	1	5	0.650	108.385	   C94
       1	2	2	2	0.747	121.550	   C94
       2	2	2	2	0.796	126.284	   E94
 
         C94  = CORE MMFF94 parameter - obtained from ab initio data
         X94  = EXTD MMFF94 parameter - fit to additional ab initio data
         E94  = theta0 from fit to X-ray data, ka from empirical rule
         #E94 = theta0 and ka both from empirical rules

   The seven columns define the angle-type index ATIJK, the atom types for the 
   wing, central, and wing atoms in the angle, the force constant in 
   md-angs/rad**2, the reference bond angle in degrees, and the origin of the
   parameter, where "C94" and "E94" have meanings analogous to those discussed 
   above. In particular, the latter combine reference angles obtained by fitting 
   to crystallographic geometries with empirical-rule force constants
   calculated as noted below. 

   The angle-bending parameters employ angle-type indices ATIJK ranging 
   between 0 and 8.  Their meanings are as defined below:

   ::
   
      ATIJK		Structural significance
       ---------------------------------------------------------------------------
         0		The angle i-j-k is a "normal" bond angle
         1 		Either bond i-j or bond j-k has a bond type of 1
         2		Both i-j and j-k have bond types of 1; the sum is 2.
         3		The angle occurs in a three-membered ring
         4		The angle occurs in a four-membered ring
         5		Is in a three-membered ring and the sum of the bond types is 1
         6		Is in a three-membered ring and the sum of the bond types is 2
         7		Is in a four-membered ring and the sum of the bond types is 1
         8		Is in a four-membered ring and the sum of the bond types is 2

   The canonical ordering of the angle parameters is such that the atom-type 
   index J changes least rapidly and the angle-type index ATIJK changes most 
   rapidly.  Atom type I is always less than or equal to K and changes less 
   rapidly than K.  All angle interactions having a common central atom type 
   thus appear together in the file.  The canonical-order index, CXA, is computed 
   as:

   ::
   
     CXA = MC*(J*MA**2 + I*MA + K) + ATIJK

   where MA is again the maximum permissible atom type +1, and MC is at least 
   one greater than the maximum permissible angle-type index.  

   The procedure used to assign force constants kaIJK and the reference 
   angles ToIJK for a particular i-j-k interaction is analogous to that described 
   above for bond stretching.  Note that the MMFF angle bending parameters also 
   include default parameters which have zero values for atom types I and K 
   (and a zero value for the listed force constant); the first parameter listed 
   above, for example, is of this type.  These parameter entries are used to 
   assign reference bond angles for interactions I-J-K of angle-type M = ATIJK 
   when the parameter file contains neither the fully-qualified M:I-J-K  
   parameter nor any related parameter obtained by successively equivalencing 
   atom types I, J and K to simpler atom types in the manner described above in 
   connection with the MMFFDEF.PAR file.  In such a case, the "M:0-J-0" entries 
   are used to assign the reference angle, and procedure based on a previously 
   published an empirical rule (but employing slightly different parameters) is 
   then employed to calculate the force constant [9]_.  (If no match on M and J is 
   found either, information relating to hydridization and ring size is used to 
   assign the reference value.)  These default reference values were generated 
   by averaging over I and K for all M:I-J-K reference bond angles determined 
   either from the MP2/6-31G* structures used in core MMFF or from the 
   crystallographic structures used in the extended CSD parameterization.  


.. _mmff_params_mmffstbn:

9. MMFFSTBN.PAR

   The form used for stretch-bend interactions is:

   ::
   
      EBAijk   = 2.51210 (kbaIJK*Drij + kbaKJI*Drkj)*DTijk,	 		(6)

      where: 	kbaIJK =   force constant in md/rad for i-j stretch coupled to 
      		   i-j-k bend.

   	kbaKJI =   force constant for k-j stretch coupled to i-j-k bend.

   and Dr and DT are defined as above.  Currently, stretch-bend interactions are 
   omitted when the i-j-k interaction corresponds to a linear bond angle.  

   The format of the MMFFSTBN.PAR file is illustrated below:

   ::
   
    	        --- Atom Types ---
        SBTIJK      I	J	K    kbaIJK    kbaKJI	Origin/Comment
     ------------------------------------------------------------------
   	0	1	1	1     0.206	0.206	C94
   	0	1	1	2     0.136	0.197	C94
   	0	1	2	1     0.250	0.250	C94
   	0	1	2	2     0.203	0.207	C94
   	2	1	2	2     0.222	0.269	C94
   	2	1	2	3     0.244	0.292	C94
   	0	1	2	5     0.215	0.128	C94
   	1	2	2	2     0.250	0.219	C94
   	2	2	2	3     0.155	0.112	C94
   	0	2	2	5     0.207	0.157	C94

      C94 - CORE MMFF94 parameter, from fits to HF/6-31G* 2nd D's
      X94 - EXTD MMFF94 parameter, also from fits to HF/6-31G* 2nd D's

   All the parameters listed in the file represent values derived by fitting to 
   HF/6-31G* energy derivatives.  Note that the format is similar to that 
   employed for angle bending.  Two force constants are listed, however. The 
   first, kbIJK, couples i-j-k bending with the stretching of the i-j bond, 
   whereas the second, kbaKJI, couples i-j-k bending to k-j stretching.  The 
   requisite reference bond lengths and bond angle are taken from the respective 
   bond and angle parameter sets, with the stretch-bend type SBTIJK listed in 
   the first column serving with the atom types to establish the proper
   connection.

   The stretch-bend types are defined in terms of the constituent bond types BTIJ 
   and BTJK and the angle type ATIJK as shown below: 

   ::
   
           Stretch-Bend     Angle	      ---- Bond Type ----
               Type	      Type	      I-J	      J-K
       -----------------------------------------------------------
      	  0		0		0		0
      	  1		1		1		0
      	  2		2		0		1
      	  3		2		1		1
      	  4		4		0		0
      	  5		3		0		0
      	  6		5		1		0
      	  7		5		0		1
      	  8		6		1		1
      	  9		7		1		0
        	 10		7		0		1
      	 11		7		1		1

   Canonical order is established through an index, CXBA, which is computed 
   from the stretch-bend type and the atom types in much the same maner as was 
   described above for the index CXA.  The ordering of the tabulated parameters 
   therefore follows the same prescription.

   The matching of stretch-bend interactions with tabulated parameters
   proceeds much as described above, i.e. by comparing the CXBA index for the
   M:I-J-K combination of interest (where M = SBTIJK) to the indices stored with
   the entries in the parameter file.  In this case, however, default parameters
   are used whenever the fully qualified M:I-J-K parameter is not found in the
   initial search.  (Thus, no atom-type equivalences are employed.)  The default
   parameters are taken from the MMFFSBDF.PAR file, as described below.


.. _mmff_params_mmffdfsb:

10. MMFFDFSB.PAR

   The format of this file is illustrated below:

   ::
   
                    Row in 
             -- Periodic Table --    Force constant     
      	IR	JR	KR     kbaIJK   kbaKJI 
        ------------------------------------------------
      	0	1	0	0.15	 0.15
      	0	1	1	0.10	 0.30
      	0	1	2	0.05	 0.35
      	0	1	3	0.05	 0.35
      	0	1	4	0.05	 0.35
      	0	2	0	0.00	 0.00
      	0	2	1	0.00	 0.15
      	0	2	2	0.00	 0.15
      	0	2	3	0.00	 0.15
      	0	2	4	0.00	 0.15
      	1	1	1	0.30	 0.30
      	1	1	2	0.30	 0.50

   Note that the tabulated force constants depend only on the rows in the 
   periodic table to which the constituent atoms i, j and k belong.  Wherever 
   possible, they represent averages of the core-MMFF stretch bend parameters 
   for each given row combination; for new row combinations, they reflect 
   trends in the core-MMFF stretch-bend force constants.


.. _mmff_params_mmffoop:

11. MMFFOOP.PAR.  For out-of-plane bending, MMFF uses the quadratic form

   ::
   
      EOOPijk;l  =  0.5 * 0.043844 * koopIJK;L * Xijk;l**2,                     (5) 

      where: 	koopIJK;L =  out-of-plane bending force constant in md-angs/rad.

      	Xijk;l    =  angle in degrees between the bond j-l and the 
      		     plane i-j-k, where j is the central atom.

   The angle Xijk;l is computed for all tricoordinate centers using the Wilson 
   definition [10]_.  Note that three such out-of-plane angles arise at each 
   trigonal center j, as any one of the attached atoms i, k and l can serve as 
   the "out of plane" atom.  In MMFF, all three such angles are assigned the 
   same koop force constant.  

   The format of the MMFFOOP.PAR file is depicted below:
 
   ::
   
           ------  Atom Types ------
           I       J       K       L     koopIJK;L	   Origin/Comment
        -------------------------------------------------------------------
   	0	2	0	0	0.020	   *-2-*-* C94 DEF
   	1	2	1	2	0.030	   C94
   	1	2	2	2	0.027	   C94
   	1	2	2	3	0.026	   C94
   	1	2	2	5	0.013	   C94
   	2	2	2	5	0.013	   C94
   	2	2	3	5	0.012	   C94
   	2	2	5	5	0.006	   C94
   	2	2	5	6	0.027	   C94

      C94  - CORE MMFF94 parameter, from fits to HF/6-31G* 2nd D's
      #C94 - Value adjusted from CORE MMFF94 value
      E94  - Value assigmed by analogy

   The canonical index, CXO, is computed as:

   ::
   
   	CXO = J*MA**3 + I*MA**2 + K*MA + L, 

   where MA is as defined previously, L is greater than or equal to K, and K is 
   greater than or equal to I.  This ordering groups together in the file all out-
   of-plane entries which share a common central atom type J.  The entries 
   for which atom types I, K and L are zero serve as default values; they 
   match out-of-plane interactions for which the fully qualified I-J-K-L 
   parameter is not found.  Wherever possible, the default force constants 
   approximate the average of the force constants for the central atom type J
   derived in the fit to the HF/6-31G* second derivatives in the core- MMFF
   parameterization.  The entries labeled "C94" were derived from the ab initio
   data , whereas entries labeled "E94" were assigned by analogy to values
   found for similar atom types for which computationally derived values were
   available. 
 

.. _mmff_params_mmfftor:

12. MMFFTOR.PAR

   MMFF employs the three-fold representation also used in
   MM2 [11]_ and MM3 [12]_: 

   ::
   
      ETijkl = 0.5 * (V1*(1 + cosW) + V2*(1 - cos2W) + V3*(1 + cos3W))          (8)

   where W is the i-j-k-l dihedral angle in degrees.  The constants V1, V2 and 
   V3 depend on the atom types I, J, K and L for atoms i, j, k and l, where i-j, 
   j-k and k-l are bonded pairs and i is not equal to l.  

   The format of the parameter file is illustrated below:
   
   ::
   
              --- Atom Types ---
      TTIJKL I     J     K     L      V1       V2        V3      Origin/Comment
      ---------------------------------------------------------------------------
       0     0     1     1     0     0.000     0.000	  0.300    C94 0:*-1-1-* Def
       5     0     1     1	 0     0.200    -0.800	  1.500    C94 5:*-1-1-* Def
       0     1     1     1	 1     0.103     0.681    0.332    C94
       5     1     1     1     1     0.144    -0.547    1.126    C94
       0     1     1     1	 2    -0.295     0.438	  0.584    C94
       0     1     1     1	 3     0.066    -0.156	  0.143    C94
       0     1     1     1     5     0.639    -0.630    0.264    C94
       0     1     1     1     6    -0.688     1.757    0.477    C94
       5     1     1     1     6     0.000     0.000    0.054    C94

         C94  - CORE MMFF94 parameter - from fits to conformational energies
         X94  - EXTD MMFF94 parameter - also from fits to conformational E's
         E94  - EXTD MMFF94 parameter - from empirical rule
         #E94 - Adjusted from empirical rule value

   The first column gives the value of the torsion type index, TTIJKL.  This 
   index normally takes the value "0", but is "1" when the j-k bond has a bond 
   type index BTJK of 1, is "2" when BTJK is "0" but BTIJ and/or BTKL is "1", is 
   "4" when i, j, k, and l are all members of the same 4-membered ring, and 
   is "5" when the four atoms are members of a 5-membered ring which is 
   not aromatic and contains no unsaturation.  

   The torsion parameters are ordered using the canonical index 

   ::
   
      CXT = MC*(J*MA**3 + K*MA**2 + I*MA + L) + TTIJKL

   In this case, as of 2/96, MC (the maximum permissible torsion-type index plus
   1) and MA (the maximum permitted atom type plus 1) are set to 6 and 136,
   respectively, to insure that CXT will fit within a 32-bit integer word.  Thus,
   J changes least rapidly and K next least rapidly.  For I, L and TTIJKL, the
   effect on the ordering can be seen in the portion of the MMFFTOR.PAR file
   displayed above.

   Requests for one-, two- and three-fold torsion parameters which can not 
   be resolved using the primary numeric MMFF atom types and the associated 
   torsion-type index TTIJKL are first handled via the step-down procedure
   outlined above in connection with the MMFFDEF.PAR file, and then, if necessary,
   by invoking an empirical rule distinct from but patterened after those used by
   Mayo, Olafson and Goddard in DREIDING [7]_ and by Rappe and co-workers in their
   successor force field, UFF [8]_. The request may terminate with fully wildcarded
   M:0-J-K-0 parameters (where for simpilicity M denotes the torsion-type index
   TTIJKL) or, in a few cases, may access "half-wildcard" M:0-J-K-L or M:I-J-K-0
   parameters. For those J-K bond types covered in the core MMFF parameterization,
   the wildcarded parameters M:0-J-K-0, M:0-J-K-L, or M:I-J-K-0 were generated by
   averaging all explicitly determined M:I-J-K-L parameters over I and L,
   respectively.  For M:I-J-K-L parameters not covered in core MMFF but required
   for the extended CSD parameterization, wildcarded parameters of the type
   M:0-J-K-0 were generated from the empirical rule just cited.  Thus, the
   additional torsion interactions covered in the extended MMFF parameterization
   utilize default parameters which depend only on the torsion type and the atom
   types of the central two atoms. 


.. _mmff_params_mmffvdw:

13. MMFFVDW.PAR

   For van der Waals interactions, MMFF employs the 
   recently developed "Buffered 14-7" form (eq 9) together with an 
   expression which relates the minimum-energy separation R*II to the 
   atomic polarizability aI (eq 10), a specially formulated combination rule 
   (eqs 11, 12), and a Slater-Kirkwood expression for the well depth eIJ (eq 
   13) [13]_: 

   ::
   
      Evdwij  =  eIJ*{1.12R*IJ/(Rij+0.12R*IJ)}**7 *  
                       {1.07 R*IJ**7/(Rij**7 + 0.07R*IJ**7) - 2}                (9)

      R*II = AI * aI**(PEXP)                                                    (10)

      R*IJ =  0.5 * (R*II + R*JJ) * (1 + AFACT(1 - exp(-BFACT*gIJ**2)))         (11) 

      gIJ = (R*II - R*JJ)/(R*II + R*JJ)                                         (12)

      eIJ =  181.16*GI*GJ*aIaJ/[(aI/NI)**0.5 + (aJ/NJ)**0.5]*R*IJ**(-6)         (13)

   The first non-comment line in the parameter file contains five floating 
   point numbers which define the variables PEXP, AFACT, BFACT, DARAD, and 
   DAEPS, respectively.  PEXP (currently 0.25), AFACT (currently 0.2) and BFACT 
   (currently 12.0) are used in the equations shown above;  DARAD and DAEPS 
   are used as explained below.  The format of the remainder of the parameter 
   file is illustrated below:  

   ::
   
         At. Type I     aI     NI       AI     GI     D_AI   Symb.   Origin/Comment
       ----------------------------------------------------------------------------
      	1	1.050	2.490	3.890	1.282    -     CR      E94
      	2	1.350	2.490	3.890	1.282    -     C=C     E94
      	3	1.100	2.490	3.890	1.282    -     C=O     E94
      	5	0.250	0.800	4.200	1.209    -     HC      C94
      	6	0.700	3.150	3.890	1.282    A     OR      C94
      	7	0.650	3.150	3.890	1.282    A     O=C     C94
             10	1.000	2.820	3.890	1.282    A     NC=O    C94
             21	0.150	0.800	4.200	1.209    D     HOR     C94
             24	0.150	0.800	4.200	1.209    D     HOCO    C94 

         E94  - From empirical rule (JACS 1992, 114, 7827) 
         C94  - Adjusted in fit to HF/6-31G* dimer energies and geometries
         X94  - Chosen in the extension of the parameterization for MMFF94 

   Column 1 contains the numerical atom type, and column 2 contains the 
   associated atomic polarizability aI in units of angs**3.   Note that the
   effective-electron numbers N and the scale factors A and G depend only on the
   atomic species for the interacting atoms i and j, not on the full MMFF numeric
   atom type. Column 6, labeled "D_AI" specifies whether the atom type in question
   is considered to be a hydrogen bond donor ("D"), acceptor ("A"), or neither
   ("-"). When either i or j is a hydrogen-bond donor, MMFF uses the simple
   arithmetic mean (0.5 (R*II + R*JJ)) instead of eq 11 to obtain R*IJ.  If the
   i-j interaction is a donor-acceptor interaction, MMFF also scales R*IJ as given
   by eq 11 by DARAD (currently 0.8) and eIJ as given by eq 13 by DAEPS (currently
   0.5).  For informational purposes only, column 7 lists a representative MFFF
   symbolic atom type. 

.. _mmff_params_mmffchg:

14. MMFFCHG.PAR

   MMFF uses the buffered 
   Coulombic form

   ::
   
      EQij  =  332.0716*qi*qj/(D*(Rij + d)),                                     (14)

   where qi and qj are partial atomic charges, Rij is the internuclear separation
   in angs, d = 0.05 angs is the "electrostatic buffering" constant, and D is the
   "dielectric constant" (normally taken as D = 1, though use of a distance-
   dependent dielectric constant is also supported).  In MMFF, 1,4- interactions
   are scaled by a factor of 0.75.  The distance buffering (d > 0) prevents
   infinite attractive electrostatic energies from overwhelming the bounded
   repulsive vdW interaction given by eq 9 as oppositely charged atomic centers
   approach. 

   The partial atomic charges qi used in eq 14 are constructed from 
   initial full or fractional formal atomic charges q0I (usually zero, but, e.g., 
   +1/3 for guanidinium nitrogens) by adding contributions from bond charge 
   increments wKI which describe the polarity of the bonds to atom i from 
   attached atoms k.  Thus, wKI is the contribution to the total charge on 
   atom i of atom type I accumulated from, and at the expense of, its bonded 
   neighbor k of atom type K.  Specifically, MMFF computes qi as

   ::
   
      qi = (1 - nI*uI)*q0I + uI*Sum(K)q0K + Sum(K)wKI			(15)

   where wKI = - wIK and where the sums on the right hand side run over the nI =
   crd(I) atoms k of MMFF atom type K directly attached to atom i.  (Crd(I) comes
   from MMFFPROP.PAR.) In this equation, q0I and q0K are the formal charges
   assigned in the atom typing procedure (usually, by subroutine XTYPE), and the
   sum of the first two terms gives the "effective" fractional formal atomic
   charge residing on atom i.  This approach allows a formal atomic charge
   initially affixed by the atom-typing procedure (e.g., q0I) to be shared in a
   prescribed manner with the neighbors bonded to the atom in question.  For
   example, for the series PO4(-3), HPO4(-2), H2PO4-, H3P04, it allows allows the
   partial charges on the terminal oxygens (each represented by the same numerical
   atom type, "32") to vary in a way which properly reflects the partial charges
   obtained from fits to the 6-31G* electrostatic potential.  In particular, the
   difference between the resultant charges qi calculated for the single terminal
   oxygen in H3PO4 and for the four equivalent terminal oxygens in PO4(-3) comes
   to -0.375, half (because u32 = -0.5) the difference of -0.75 in the q0K charges
   (i.e., 0.00 and -0.75, respectively) and reasonably in accord with the
   difference of -0.42 found by fitting the electrostatic potential. 

   The bond charge wKI increments  are taken from MMFFCHG.PAR, and the 
   disbursement-fraction factors uI come from MMFFPBCI.PAR.  The format of 
   the MMFFCHG.PAR file is as shown below:

   ::
   
          BTIK     I       K	 wIK	   Origin/Comment
         --------------------------------------------
   	0	1	1	0.0000	   #C94
   	0	1	2      -0.1382	   C94
   	0	1	3      -0.0610	   #C94
   	0	1	5	0.0000	   #C94
   	0	1	6      -0.2800	   #C94 
   	0	1      10      -0.3001	   C94
   	0	2	2	0.0000	   #C94
   	1	2	2	0.0000	   #C94
   	1	2	3      -0.0144	   C94

       C94  - CORE MMFF94 parameter, obtained from fit to dipole moments
       #C94 - CORE MMFF94 parameter, either fixed or adjusted to fit 
              HF/6-31G* interaction energies
       X94  - EXTD MMFF94 parameter, obtained from fit to dipole moments
       #X94 - EXTD MMFF94 parameter, fixed in the fit to dipole moments
       E94  - derived from partial bond charge increments (empirical rule)

   The BTIK are the previously discussed bond types.  Note that the "0 1 2" 
   line, for example, indicates that a (carbon) atom of type "2" in a "1 - 2" C-
   C bond gains a charge of -0.1382 from the polarity of the bond, while the 
   atom of type "1" gains an equal and opposite charge of +0.1382.  The 
   parameter lines identified by the label "#C94" list bond charge 
   increments which were held fixed in the fit to molecular dipole moments used
   in the derivation of the "C94" values.

.. _mmff_params_mmffpbci:

15. MMFFPBCI.PAR

   The format of the MMFFPBCI.PAR file as illustrated below:

   ::
   
          Atom			
         type I       pI        uI     Origin/Comment
         ---------------------------------------------
   	1	 0.000	    0.00    E94	
   	2    	-0.135	    0.00    E94		
   	3     	-0.095	    0.00    E94		
          30	-0.166	    0.00    E94		
          31	 0.161	    0.00    E94		
          32 	-0.732	    0.50    E94		

   The uI parameters, which usually are zero, are used in the equation given 
   above for qi.  The pI are used to estimate missing bond-charge increments as

   ::
   
      wKI = pI - pK 								(16)

   Because of the additive decomposition they give of the bond charge increments
   wKI, the pI are called "partial bond charge increments".  They have been shown
   to reproduce the ab-initio parameterized wKI with an rms error of about 20% and
   were used to generate the bond charge increments labeled "E94" in the
   MMFFCHG.PAR file. 

.. _mmff_params_mmffsup:

16. MMFFSUP.PAR

   This file allows additional MMFF parameters of various 
   types to be added, and also allows existing parameters to be redefined (i.e.,
   overwritten, replaced).  The supplementary parameters are read from a file
   assigned to unit "INP".  The supplementary-parameters file consists of a header
   card, followed by sets of cards specifying supplementary parameters of various
   types.  

   Note: the current (February 1996) CHARMM implementation is incomplete.  In 
   this implementation, the subtype sections *in some cases* can appear in any
   order.  In particular, the VDW, BOND, ANGLE, and TORSION sections are 
   initiated by a header record of that name.   These headers, which *must* be 
   present in these four cases, allow CHARMM to position the supplementary-
   parameters file at the first "real" data record for that block. CHARMM then 
   reads just the requested number of data records (e.g., NB  records of "BOND" 
   type).  For readability, it is therefore permissible to close the block 
   with a record like "END_VDW."   Alternatively, or in addition, any number of 
   comment records of arbitrary format can be added after the "real" data.  

   Unfortunately, the current implementation treats blocks of bond-charge-
   increment and out-of-plane bending parameters differently.  Here, *no* header
   record is allowed.  Rather, as in other implementations of the MMFF parameter
   reading facility, the code assumes that the supplementary-parameters file will
   be read sequentially and will be positioned properly when supplementary
   parameters of any given subtype, including these two, are to be read.  

   Fortunately, there is a work-around.  First, order the reading of the main
   MMFF parameter files as is shown in the mmff_setup.STR example in the mmff.doc
   file.  Second, specify the supplementary parameters in the order shown below.
   Third, when the associated NV, NB, NA, and NT count variables are nonzero, 
   include VDW, BOND, ANGLE, and TORSION headers, but *do not* supply an END
   card, or any other data, following the parameter data cards.  In this way, if
   NQ, for example, is nonzero but NV or NB is not, the supplementary 
   parameters file will have been left positioned at the first bond-charge-
   increment parameter line at the time the reading of this data is called for.

   To determine when the implementation for reading supplementary parameters has
   been made complete, go to the mmff source directory and issue the command:
   "grep -i find_loc \*.src".  You will see the file-positioning find_loc commands 
   bearing the arguments "VDW", "BOND", "TORSION" AND "ANGLE".  If you also see 
   others  bearing such arguments as "CHARGE" and "OUTPLANE", you will know that 
   these blocks also require header records.  At this point, it will also be
   possible to add "trailer" records and to position the subblocks in arbitrary
   order.  


   a. First read header card for the supplementary-parameters file:
   
      ::

          READ(INP,'(8I5)') NV, NB, NQ, NP, NA, NO, NBA, NT

          NV -- NT are numbers of supplementary parameters:

             Quantity and definition                     Block header
             ------------------------------------        ---------------------------
             NV -- van der waals parameters              VDW 
             NB -- bond stretch paramaters               BOND 
             NQ -- bond charge increment parameters      
             NP -- partial bond charge increment  
                   parameters (not yet implemented)
             NA -- angle bending parameters              ANGLE 
             NO -- out-of-plane bending parameters
             NBA - stretch-bend interaction (not yet 
                   implemented)
             NT -- torsion parameters                    TORSION 

      The block identifiers need to be specified in addition to the numerical
      data defined in the fllowing sections.


   b. If NV is nonzero, first read new values for various vdw control parameters:

      ::
      
          READ(INP,'(5F10.5)') PEXP, AFACT, BFACT, DARAD, DAEPS

          NOTE: each nonzero value replaces the value supplied in the first line of 
          the MMFFVDW.PAR file (zero values are ignored).  See the description of
          the MMFFVDW.PAR file for the meaning of these parameters.

          Then read NV supplementary vdW parameters:

          DO I=1,NV
            READ(INP,'(I5,4F10.3,1X,A1)') I, aI, NI, AI, GI, D_AI
          ENDDO

      These six quantities provide additional or replacement values for the
      corresponding items on the MMFFVDW.PAR parameter lines.

      ::
      
         I     -  MMFF numeric atom type
         aI    -  polarizability for an atom of type I
         NI    -  Slater-Kirkwood effective electron number
         AI    -  "A" value
         GI    -  "G" value
         D_AI  -  Donor ("D"), acceptor ("A"), or "neither" ("-") character string
         
   c. If nonzero, read "NB" supplementary stretching parameters: 

      ::
      
         DO I=1,NB
            READ(INP,'(I1,I4,I5,2F8.3,A)') BTIJ, I, J, kbIJ, roIJ, COMMENT
         ENDDO 

         BTIJ    - bond-type index for this bond between atoms of types I and J
         I, J    - MMFF atom types for atoms i and j
         kbIJ    - stretching force constant in md/angstrom
         roIJ    - reference bond length in angstroms
         COMMENT - optional comment field

         
   d. If nonzero, read "NQ" supplementary bond charge increments:

      ::
      
          DO K=1,NQ
             READ(INP,'(I1,I4,I5,2F8.3,A)') BTIK, I, K, wIK
          ENDDO
      
          BTIK    - bond-type index for this bond between atoms of types I and K.
          I, K    - MMFF atom types for atoms i and k.
          wIK     - bond-charge increment in electron units, added to the atom of
                    higher atom type and subtracted from the atom of lower atom type


   e. If nonzero, read "NP" partial bond charge increments (not yet implemented
      as of February 1996):

      ::
      
         DO K=1,NP
            READ(INP,'(I1,I4,I5,F10.5)') I, pI, uI
         ENDDO
      
         I      - MMFF atom type 
         pI     - partial-bond-increment change in electrons, used
                  to estimate "missing" bond charge increments
         uI     - formal charge distribution factor
         
   f. If nonzero, read "NA" supplementary angle bending parameters: 

      ::
      
          DO I=1,NA
             READ(INP,'(I1,I4,2I5,4F8.3,A)') ATIJK, I, J, K, kaIJK, ToIJK, COMMENT
          ENDDO 
      
          ATIJK      -   angle-type index
          I, J, K    -   MMFF atom types for atoms i, j, k
          kaIJK      -   bending force constant in md-angs/rad**2
          ToIJK      -   reference bond angle in degrees
          COMMENT    -   optional comment field

         
   g. If nonzero, read "NO" supplementary out-of-plane parameters: 

      ::
      
          DO I=1,NO
              READ(INP,'(4I5,F10.3,A)') I, J, K, L, koopIJK;L, COMMENT
          ENDDO

          I, J, K, L  -  MMFF atom types for atoms i, j, k and l, where j is the
                         central atom
          koopIJK;L   -  out-of-plane bending force constant in md-angs/rad**2
          COMMENT     -  optional comment field
      
         
   h. If nonzero, read "NBA" supplementary stretch-bending parameters (not yet 
      implemented as of February 1996): 

      ::
      
          DO I=1,NBA
             READ(INP,'(I1,I4,2I5,4F8.3,A)') SBTIJK, I, J, K, kbaIJK, kbaKJI, COMMENT
          ENDDO
      
          SBTIJK       -   stretch-bend-type index
          I, J, K      -   MMFF atom types for atoms i, j, k
          kbaIJK       -   stretch-bending force constant in md/rad
                           coupling i-j stretching with i-j-k bending
          kbaKJI       -   stretch-bending force constant in md/rad
                           coupling k-j stretching with i-j-k bending
          COMMENT      -   optional comment field


   i. If nonzero, read "NT" supplementary torsion angle parameters: 

      ::
      
          DO I=1,NT
             READ(INP,'(I1,I4,3I5,3F8.3,A)') TTIJKL, I, J, K, L, V1, V2, V3, COMMENT
          ENDDO

          TTIJKL         - torsion-type index
          I, J, K, L     - MMFF atom types for atoms i, j, k, l
          V1, V2, V3     - one-fold, two-fold, three-fold barrier constants in 
                           kcal/mol
          COMMENT        - optional comment field


References
----------

.. [1] Roy, R. S. Proc. Phys. Soc., Ser. 2, 1968, 1, 445-448.
.. [2] Peterson, M. R.; Csizmadia, I. G. J. Mol. Struct. 1985, 123, 399-412; cf. 
       Table 1. 
.. [3] Herschbach, D. R.; Laurie, V. W. J. Chem. Phys. 1961, 35, 458-463.  
.. [4] Schomaker, V.; Stevenson, D. P. J. Am. Chem. Soc. 1941, 63, 37-40.
.. [5] Blom, R.; Haaland, A. J. Mol. Struct. 1985, 128, 21-27.
.. [6] Allred, A. L.; Rochow, E. G. J. Inorg. Nucl. Chem. 1958, 5, 264.
.. [7] Mayo, S. L.; Olafson, B. D.; Goddard III, W. A. J. Phys. Chem. 1990, 94, 8897. 
.. [8] Rappe, A. K.; Casewit, C. J.; Colwell, K. S.; Goddard III, W. A; Skiff, W.
       M. J. Am. Chem. Soc. 1992, 114, 10024-10035, and references therein. 
.. [9] Halgren, T. A. J. Am. Chem. Soc. 1990, 112, 4710-4723.
.. [10] Wilson, E. B., Jr; Decius, J. C.; Cross, P. C., Molecular Vibrations;
        Dover: New York, 1955, Chapter 4.
.. [11] (a)Allinger, N. L. J. Am. Chem. Soc. 1977, 89, 8127.  (b) Bukert, U.; 
        Allinger, N. L. Molecular Mechanics; American Chemical Society: 
        Washington, DC, 1982.  (c) Allinger, N.L.; Yuh, Y. QCPE 1980, 12, 395.
.. [12] Allinger, N. L.; Yuh, Y. H.; Lii, J.-H. J. Am Chem. Soc., 1989, 111, 8551.
.. [13] Halgren, T. A. J. Am. Chem. Soc. 1992, 114, 7827-7843.
 