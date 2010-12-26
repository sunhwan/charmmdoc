.. py:module:: io

####################################
Input-Output Commands
####################################

The commands described here are used for reading and writing
data structures used in the main part of CHARMM. Some of data structures
used in the analysis facility may also be read and written.

.. index:: io; read
.. _io_read:

Reads Data from External Sources
--------------------------------

This command reads data into the data structures from external
sources. The external sources can be either card image files or binary
files. The fortran unit number from which the information is read, is
specified with the unit-spec.

The precise format of all these files is described only in the
source code as that serves as the only definitive, accurate, and up to
date description of these formats. The description of the data
structures provides pointers to the subroutines which should be
consulted, see :ref:`Data Structure <usage_data>`.

.. index:: read; syntax
.. _io_read_syntax:

Syntax of READ Command
^^^^^^^^^^^^^^^^^^^^^^

::

   READ { RTF {   CARD  [APPEnd] [PRINt] }                } [ UNIT integer
                                                               | NAME filename ]
        {     { [ FILE ]                 }                }

        { PARAmeter {  CARD  [PRINt]  [APPEnd]  [FLEX] }  }
        {           { [FILE] [NBON ]            [MMFF] }  }

        { IC  { [ CARD ] } [ APPEnd ] [ SAVEd ]           }
        {     {   FILE   }                                }

        { SEQUence  { [ CARD ]                   }        }
        {           {   COOR  [ RESId ]          }        }
        {           {   PDB   [ RESId ] [CHAIn <CHAR>] }  }
        {           {   TIPS  integer            }        } ! TIP3P water model
        {           {   ST2   integer            }        } ! ST2 water model
        {           {   WATEr integer            }        } ! OH2 residue model
        {           {   DUM   integer            }        } ! Dummy atoms
        {           {   resname integer          }        } ! Any RESI in the RTF

        { HBONd  { [ FILE ] }                             }
        {        {   CARD   }                             }

        { PSF    { [ FILE ] }  [APPEend]                  }
        {        {   CARD   }                             }

        { CONStraint { [ CARD ] }                         }

        { NBONd [ FILE ]                                  }

        { TABLe [ FILE ]                                  }

        { TRAJectory [ COMP ]                             }

        { IMAGes [ CARD ] [ INIT ]                        }

        { XRAY                                            }

        { UNIVersal-coordinate-format                     }

        { COORdinate coor-spec [ COMP ]                   }

        { SEGId segment  { PDB     }  [BUILd [SETUp]]     }
        {                { CARD    }                      }
        {                { FREE    }                      }
        { NAMD FILE "filename" **                         }
    ** no other options are available

    coor-spec ::= { FILE      [IFILE int]                 } coor-option
                  { CONTinue                              }
                  { CARD      [OFFS  int] [ RESI ]        }
                  { PDB    [OFFS  int] [MODEL int] [OFFI] }
                  { UNIVersal [OFFS  int] [ RESI ]        }
                  { IGNOre                                }
                  { DYNR      CURR|DELT|VEL               }
                  { TARG                                  }
                  { TAR2                                  }

    coor-option ::=  [APPEnd] [INITial] [FREEfield] atom-selection

   Syntactic ordering: The second field must be specified as shown.
   The file to be read can be specified either through the UNIT number (the
   same number as in a preceding OPEN statement) or through the NAME keyword.


Reading the sequence and coordinates of segments simultaneously
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This command provides convenient way to transform a system in PDB file
format into new CHARMM segments with given coordinates.  When read in segments
from a PDB file, one can specify :chm:`BUILd` to generate all atom connectivities and
atom types. If there are missing atoms in the PDB file, one can specify :chm:`SETUp`
to generate an internal coordinate table of the segments to be used to 
generate the coordinates of those missing atoms.  Each chain in the PDB file
will form a new segment named as the given :chm:`SEGId` followed by its segment 
number. These generated segments are well qualified CHARMM segments and
can be used for atom based simulation. This is a very convenient way to 
generate simulation systems from PDB files. However, It requires that all
residue and atom names in the input file are consistent with that in the
CHARMM RTF file.

For example: 

::

          open read unit 10 card name 1b5s.pdb
          read segid b5s PDB build setup unit 10
          
This command can be used to create a new segment from either a
PDB file (PDB), a CHARMM coordinate file (:chm:`CARD`), or a free format coordinate
file (:chm:`FREE`). If :chm:`BUILd`  option is not specified, the generated 
segment contains only atoms listed in the input PDB file but no atomic 
connectivities are generated.  Such a segment can be used to generate a map 
object needed in the :doc:`EMAP module <emap>`. With this command, a map 
object can be quickly converted from a PDB structure.


.. index:: read; sequence

Specifying a sequence of residues for a segment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The specification of :chm:`SEQUence` causes the program to accept a
sequence of residue names to be used to generate the next segment in the
molecule. Unless the :chm:`WATEr`, :chm:`TIPS`, or :chm:`ST2` option is used, the sequence is
specified as follows:

::

        title
        number of residues
        repeat(residue names)

The form of the title is defined in the syntactic glossary,

.. note::

   :ref:`Syntactic Glossary <usage_syn>` The number of
   residues is specified on the line following the title in free field
   format. If the number of residues you specify is less than zero,
   CHARMM will read residues until it encounters a blank line or end of
   file. If the number is greater than zero, it will also stop once it
   has read at least as many residues as you've specified. If the number
   you specify is zero, you will get a warning message as one common
   error is to forget the number entirely. In this case, the first
   residue name will be consumed as the number and converted to zero.
   
The residue names are specified as separate words, each no
longer than 4 characters, on as many lines as are required for all the
residues. This sequence may be placed immediately following the READ
command if the unit number is the stream or may be placed in a separate file.

When reading is complete, CHARMM will list all the residues it
has read, and tell you which residues it thinks can be titrated.

The :chm:`WATEr` option allows a sequence of water molecules to be
specified. This will give the old 3-center water model (not recommended).
The integer which follows the keyword gives the number of waters. The TIP3P
water model may be specified with the :chm:`TIPS` option. Likewise, the :chm:`ST2`
option allows ST2 waters to be specified. Obviously, no sequence on
separate lines need be given. The topology file must contain the residue
named (OH2,TIP3,ST2); otherwise, the :chm:`GENErate` command invoked subsequently
will fail.

The :chm:`COOR` option will read the sequence from a CHARMM format card
coordinate file. The residue numbers are ignored except that when a change
occurs, a new residue is added. If the :chm:`RESId` keyword is also present,
then the resid's are obtained from the resid field of the coordinate file.
For the PDB option resids are always read from the resSeq (resid) field.
This is useful when one wants to specify residue names (rather than use
the number representation). No other information is read from the coordinate
file during this process. To read the sequence for a specific chain in a PDB
file the :chm:`CHAIn <char>` option can be used; `<char>` is the one letter PDB chain
id in position 22 of ATOM/HETATM records.


.. index:: read; coordinates

Reading coordinates
^^^^^^^^^^^^^^^^^^^

The reading of coordinates is done with the :chm:`READ COOR` command,
and there are several options (which may change over in future versions).
Coordinates may be read into the main set or the comparison coordinate set
using the COMP keyword.

There are three possible file formats that can be used to read
in coordinates. They are coordinate binary files, dynamics coordinate
trajectories, and coordinate card images. In addition, NAMD program
binary restart coordinates(and velocities) files can be read (only
into main set). Protein Data Bank (PDB) formatted files can also be
read. PDB files do however require some editing first. All the HEADER
and other junk before the actual coordinate section has to be removed
and optionally replaced by a standard CHARMM title. There should be no
line with NATOM (= number of atoms) preceding the actual coordinates.
CHARMM does no translation whatsoever of residue or atom names, so you
would either have to rename some entries in the PSF or in the
coordinate file in case there are differences. The MODEL option reads
the specified MODEL number from an NMR style multiple coordinate set
PDB file.

For all formats, a subset of the atoms in the PSF may be selected
using the standard atom selection syntax. For binary files, This is a
risky maneuver, and warning messages are given when this is attempted.
Only coordinates of selected atoms may be modified. When reading binary
files, or using the IGNOre keyword, coordinate values are mapped into
the selected atoms sequentially (NO checking is done!).
Selection of atoms does not work with NAMD binary files (example:
`read namd file "myfile.coor.rst"`
)

The reading of the first two file formats is specified with the
:chm:`FILE` option. The program reads the file header to tell which format it
is dealing with. The coordinate binary files have a file header of
:chm`COOR` and contain only one set of coordinates. These are created with a
:chm:`WRIT COOR FILE` command. The dynamics coordinate trajectories have a file
header of :chm:`CORD` and have multiple coordinate sets. These files are
created by the dynamics function of the program. To specify which
coordinate set in the trajectory to be read, the :chm:`IFILE` option is
provided. One specifies the coordinates position within the file. The
default value for this option will cause the first coordinate set to be
read. If the :chm:`IFILE` value is negative, then the next file (other than
the first one) will be read. This will only work if a set has already been
read from the file with a positive :chm:`IFILE` value.

For binary files, the :chm:`APPEnd` command will 'deselect' all atoms
up to the highest one with a known position. This is done in addition
to the normal atom selection. This is useful for structures with several
distinct segments where it is desirable to keep separate coordinate
modules.

The :chm:`CARD` file format is the standard means in CHARMM for
providing a human readable and writable coordinate file. The format is
as follows:

* Normal format for less than 100000 atoms and PSF IDs with less than
  five characters
   
  ::
   
         title
         NATOM (I5)
         ATOMNO RESNO   RES  TYPE  X     Y     Z   SEGID RESID Weighting
           I5    I5  1X A4 1X A4 F10.5 F10.5 F10.5 1X A4 1X A4 F10.5

* Expanded format for more than 100000 atoms (up to 10**10) and with
  up to 8 character PSF IDs. (versions c31a1 and later)
  
  ::
  
         title
         NATOM (I10)
         ATOMNO RESNO   RES  TYPE  X     Y     Z   SEGID RESID Weighting
           I10   I10 2X A8 2X A8       3F20.10     2X A8 2X A8 F20.10

The title is a title for the coordinates, see :ref:`Syntactic Glossary <usage_syn>`,
for details. Next comes the number of coordinates. If this number is zero or too large,
the entire file will be read. Finally, there is one line for each coordinate.

``ATOMNO`` gives the number of the atom in the file. It is ignored
on reading. `RESNO` gives the residue number of the atom. It must be
specified relative to the first residue in the PSF. The :chm:`OFFSet` option
should be specified if one wishes to read coordinates into other positions.
The :chm:`APPEnd` option adds an additional offset which points to the
the residue just beyond the highest one with known positions. This option
also 'deselects' all atoms below this residue (inclusive).
For example, if one is reading in coordinates for the second segment of a
two chain protein using two card files, and the :chm:`APPEnd` option is used,
``RESNO`` must start at 1 in both files for the file reading to work
correctly.

It should also be remembered that for card images, residues are
identified by ``RESNO``. This number can be modified by using the
:chm:`OFFSet feature`, which allows coordinates to be read from a different PSF.
Both positive and negative values are allowed. The :chm:`RESId` option will
cause the residue number field to be ignored and map atoms from ``SEGID``
and ``RESID`` labels instead.

``RES`` gives the residue type of the atom. ``RES`` is checked against
the residue type in the PSF for consistency. ``TYPE`` gives the IUPAC name
of the atom. The coordinates of an atom within a residue need not be
specified in any particular order. A search is made within each residue
in the PSF for an atom whose IUPAC name is given in the coordinate file.

The :chm:`RESId` option overrides the residue number and fills coordinates
based on the SEGID and RESID identifiers in the coordinate file.
This is the recommended method where different PSF's are used.

The :chm:`IGNORE` option allows one to read in a card coordinate file
while bypassing the normal tests of the residue name, number, and atom
name. When :chm:`IGNORE` is specified in place of card, the identifying
information is ignored completely. Starting from the first selected
atom, the coordinates are copied sequentially from the file.

The PDB option works very much like the CARD option, but expects the
actual file format to be according to Protein Data Bank standards:

::

  text IATOM  TYPE  RES  IRES      X  Y  Z    W
   A6   I5  2X A4   A4    I5  4X     3F8.3 6X F6.2

The :chm:`OFFI` option enforces the official PDB format. The `SEGID` (chain id)
has to be one character in length on read  and it is truncated
to one character on write.

Normally, the coordinates are not reinitialized before new values
are read, but if this is desired, the :chm:`INITialize` keyword, will cause the
coordinate values for all selected atoms to be initialized. Note that only
atoms that have been selected, will be initialized (9999.0). The :chm:`COOR INIT`
command provides a more general way to initialize coordinates.

The :chm:`READ COOR DYNR` variant reads a full coordinate set from a dynamics
restart file. It **REQUIRES** a matching PSF and allows no selections or
other manipulations. A restart file (usually) contains three sets of
atom data, and you chose which one to read in with keywords: 

   ======   ================================================================
   CURR     the current coordinates
   DELT     the displacement to be taken from the current coordinates
   VEL      the current velocities (in AKMA units)
   ======   ================================================================

.. note::

   The restart file written after a crash may be slightly different,
   at present (c28a2) it contains the previous coordinates instead of velocities.

The :chm:`READ COOR TARG` and :chm:`READ COOR TAR2` commands read in the coordinates of the 
target for Targeted Molecular Dynamics (see :doc:`TMD <tmd>`)


.. index:: read; universal

Reading coordinates from nonstandard formats
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The reading of coordinates is done with the :chm:`READ COOR` command,
and there are several options.  One such option is the :chm:`READ COOR
UNIVersal` command which will read using a previously specified format.
The Universal format is specified by the :chm:`READ UNIVersal` command.  This
reads the specification from the input stream or from a specified
file.

::

   READ UNIVersal

The following commands clear the translation table and sets up
default specifications for the file format.

   =======  ================================================================
   CHARMM   setup standard CHARMM format (default)
   PDB      setup brookhaven format
   AMBER    setup standard AMBER  format 
   UNKNown  setup null format (everything must be specified)
   =======  ================================================================

The following commands specify the field locations of various items
When reading free-of-field, the starting values are sorted to determine
the ordering of parsing.

::

   SEGID start length
   RESID start length
   TYPE  start length
   RESN  start length
   IRES  start length
   ISEQ  start length
   X     start length
   Y     start length
   Z     start length
   W     start length

The following commands specify how input lines should be considered.

::

   PICK  start length  string    ! choose only line that match one or more of these 
   EXCL  start length  string    ! exclude any line that contains one of these 
   TITL  start length  string    ! add any line containing one of these to the title 

The following commands specify character translation upon reading the file.

::

   TRANslate { SEGID external-segid internal-segid                         }
            { RESID external-resid internal-resid match-segid             }
            { RESN  external-resn  internal-resn  match-segid             }
            { TYPE  external-type  internal-type  match-resn  match-segid }

   END   ! terminate reading universal file format


.. index:: read; parameter

The Format of Parameter Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

   READ { PARAmeter {  CARD  [PRINt]  [APPEnd]  [FLEX] }  }
        {           { [FILE] [NBON ]            [MMFF] }  }

The :chm:`CARD`/:chm:`FILE` keywords specify a card (readable) or binary file format.

The :chm:`PRINT` and :chm:`NBON` options determine the extend of printing while
reading parameters.  The :chm:`NBON` will list the ``NATVDW*(NATVDW+1)/2`` vdw table.

The :chm:`APPEnd` keyword will add the new parameters to the existing parameter
set. APPEnd does not work with binary files, MMFF, CFF, SPAS.  Also,
only parameters of the same type (e.g. both :chm:`FLEXible`) may be appended.

The :chm:`MMFF` keyword invokes the Merck Force Field parameter reader (see :doc:`MMFF <mmff>`)

The :chm:`FLEX` keyword specifies the new flexible parameter format.  This is
the same as the standard CHARMM parameter format, but;

1. allows general wildcarding for all terms
2. allows parameter substitution for missing parameters
3. does not require a previously read RTF (no global MASSES list required)
4. allows the definition of parameter equivalence groups.


Parameters can be read from cards or binary modules by the routine
PARRDR.  After the title, card file data is divided into sections beginning
with a keyword line and followed by data lines read free field:

::

        ATOM                             (Flexible parameters only)
         MASS   code   type   mass       (Flexible parameters only)

        EQUIvalence                      (Flexible parameters only)
         group  atom [ repeat(atom) ]    (Flexible parameters only)

        BOND
         atom atom force_constant distance

        ANGLe or THETA
         atom atom atom force_constant theta_min UB_force_constant UB_rmin

        DIHE or PHI
         atom atom atom atom force_constant periodicity phase

        IMPRoper or IMPHI
         atom atom atom atom force_constant periodicity phase

        CMAP
         atom atom atom atom atom atom atom atom resolution
         <...cmap data...>
         
        NBONd or NONB  [nonbond-defaults]
         atom* polarizability  e  vdW_radius -
              [1-4 polarizability  e  vdW_radius]

        NBFIX
         atom_i* atom_j*  emin rmin [ emin14 [ rmin14 ]]


        HBOND [AEXP ia] [REXP ir] [AHEX ih] [AAEX iaa] [hbond-defaults]
         donor-heavy-atom* acceptor-heavy-atom* well_depth distance


       ( SPAS only paramter types )

            FLUC
             atom chi_value zeta_value prin_integer chma_value
         
            KAPPa
             atom atom atom atom atom atom force_constant 
         
            LCH2
             atom atom atom atom atom force_constant 
         
            14TG
             atom atom atom atom trans_const gauche_const
         

        PRINt [ON ]
              [OFF]


     where '*' allows wildcard specifications:
      *  matches any string of characters (including none),
      %  matches any single character,
      #  matches any string of digits (including none),
      +  matches any single digit.

   ---------------------------------------------------------------------------
   nonbond-defaults::= [NBXMod int] [CUTNB real] [CTOFNB real] [CTONNB real]
                             [WMIN real] [E14Fac real] [EPS real]

     [ATOM ] [CDIElectric] [SHIFt  ] [VATOm ] [VSWItch ] [BYGRoup] [GEOMetric ]
     [GROUp] [RDIElectric] [SWITch ] [VGROup] [VSHIft  ] [BYCUbe ] [ARIThmetic]
                           [FSWITch]          [VFSWitch]
                           [FSHIft ]

   hbond-defaults::= [ ACCEptor   ] [ HBEXclude   ] [ BEST ]
                     [ NOACceptor ] [ HBNOexclude ] [ ALL  ]

           [CUTHB real] [CTOFHB real] [CTONHB real]
               [CUTHA real] [CTOFHA real] [CTONHA real]
                   [REXP int(def12)] [AEXP int(def10)]
                       [HAEX int(def4)] [AAEX int(def2)]
   ---------------------------------------------------------------------------

Sections end with the occurrence of the next keyword line, or a line with
the word END, the latter terminating parameter reading.

Errors in the input file will result in warning messages but not
termination of the run.

No wildcard usage is allowed for bonds and angles. For dihedrals,
two types are allowed; A - B - C - D (all four atoms specified) and
X - A - B - X (only middle two atoms specified). Double dihedral
specifications may be specified for the four atom type by listing a
given set twice. When specifying this type in the topology file, specify
a dihedral twice (with nothing intervening) and both forms will be used.

There are five choices for wildcard usage for improper dihedrals;

1. A - B - C - D  (all four atoms, double specification allowed)
2. A - X - X - B
3. X - A - B - C
4. X - A - B - X
5. X - X - A - B

When classifying an improper dihedral, the first acceptable match (from
the above order) is chosen. The match may be made in either direction
( A - B - C - D = D - C - B - A).

The periodicity value for dihedrals and improper dihedral terms
must be an integer. If it is positive, then a cosine functional form is used.
Only positive values of 1,2,3,4,5 and 6 are allowed for the vector, parallel
vector and cray routines. Slow and scalar routines can use any positive
integer and thus dihedral constrains can be of any periodicity.

Reference angle 0.0 and 180.0 degree correspond to minimum in staggered 
and eclipsed respectively. Any reference angle is allowed. The value
180 should be preferred over -180 since it is parsed faster and more
accurately. When the periodicity is given as zero, for OTHER THAN THE
FIRST dihedral in a multiple dihedral set, then a the amplitude is a
constant added to the energy. This is needed to effect the
Ryckaert-Bellemans potential for hydrocarbons (see below). 

The normal dihedral energy equation is:

::

      E = K * ( 1.0 + cos( periodicity * phi - phase ) )

When the periodicity is given as zero, then a harmonic restoring potential
in (phi - phi_min) is used. The phase value gives phi_min for this option.
This functional form is identical to that reported in the CHARMM paper,
except that either functional form (referred to as proper and improper)
may be used for dihedrals and improper dihedrals. The distinction between
these terms is that separate lookup tables are kept and the default atom
choices are still different. For dihedrals, the selection is usually based
on the middle two atoms, and for improper dihedrals, the selection is based on
the outer two atoms. For either terms, all 4 atoms may be required.

The HBOND line can be used to specify exponents for the hbond function,
with ia and ir being the attractive and repulsive radial terms and
ih and iaa the cosine exponents on the angular terms at the h and a
respectively. Defaults 4, 6, 4, and 2 respectively.

For atom types with no NBOND parameters given, no van der Waals
interactions will be calculated.  You will be warned, but be careful.

The nbond parameters for 1-4 interactions can be specified by placing the
extra set of parameters after the first.  By default the same parameters
will be used for 1-4 and all other interactions.

NON-BOND parameter combination rules depend on how the parameters are listed.
If the second number is negative, it is used as Emin, and

::

        Emin(ij)=-sqrt(Emin(i)*Emin(j)).
        
If the second number is positive, it is used as Neff, and the Slater Kirkwood
formula is used to compute Emin(ij).

The PARRDR card field ,NBFIX, allows individual atom type
van der Waals pair interactions to be specified. Subsequent lines must have;

::

        atom_i atom_j  emin rmin [ emin14 [ rmin14 ]]

If emin is positive, a severe warning is issued. The wildcard "X" may
be given. In the case where both atoms are wildcards, the entire
nbond parameter set will be modified.

If emin14 and rmin14 are not specified, then the value of emin
and rmin will be used. NOTE: The previous value will not be used.

NBFIXes are processed in order. For that reason, wildcard
usage should come first. In the case of duplicate specifications,
there is no check, and the last specification will be used.

The maximum number of NBFIX entries is currently set at 150.
The space for this is allocated in PARMIO.

PARAMETER I/O ADDENDUM:
^^^^^^^^^^^^^^^^^^^^^^^

In order to calculate the Ryckaert-Bellemans torsional potential for butane
and other extended atom hydrocarbons, the following terms should be included
in the parameter file:

::

   V = gamma[1.116 - 1.462cos(phi) - 1.578 cos**2(phi) + 0.368 cos**3(phi)
                   + 3.156 cos**4(phi) - 3.788 cos**5(phi)]
       and gamma = 1.987 kcal/mol

J. P. Ryckaert and A. Bellemans, Chem. Phys. Lett. 30, 123 (1975).
J. P. Ryckaert and A. Bellemans, Disc. Farad. Soc. 66,  95 (1978).

::

   PHI
   ! Ryckaert Bellemans has trans = 0.0
   ! since cos is an even function cos(-phi)=cos(phi), invert the 
   ! sign of the coefficients with odd power of cos(phi)
   CH3E CH2E CH2E CH3E   0.470467 5   0.0
   CH3E CH2E CH2E CH3E   0.783947 4   0.0
   CH3E CH2E CH2E CH3E   2.53516  3   0.0
   CH3E CH2E CH2E CH3E   1.56789  2   0.0
   CH3E CH2E CH2E CH3E   2.34787  1   0.0
   CH3E CH2E CH2E CH3E  -4.70368  0   0.0

The potential should be used with SHAKE bonds and angles or bonds only
as required.  The zero periodicity (constant) term should NOT be the
first in the set, otherwise it will be treated as an improper torsion.


The Format of a Residue Topology File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here is a description of what is currently (24-May-1982) in
residue topology files (as they are stored in ascii files). You may use
this format if you specify the CARD option in the READ command. The
format of binary files depends on the current implementation of the RTF
data structure (see RTF.FCM).

The purpose of residue topology files is to store the
information for generating a representation of macromolecule from its
sequence. These files are read by RTFRDR a subroutine in RTFIO which
should be be consulted for formats and the final word on what is
actually done with these files.

The residue topology files are named RTOP... .  There are two
forms, binary module (.MOD) and card format (usually .INP).  The card
format files are used only for creating binary modules and therefore
are structured as input files for CHARMM, beginning with a run title
and the command READ RTF CARD, followed by the actual topology file.

The first section of the topology files is a title section in
the usual format of up to ten lines delimited by a line containing only
a * in column 1.

The remaining information is read in free field format as
commands to define the RTF. The ordering of the commands is important
in that some information is needed to define others (i.e. the atoms
of a residue must be defined before the bonds between them).
The recommended structure of this file is:

::

        Initial setup:
                MASS specification for each atom type
                DECLarations of out of segment definitions
                DEFAults for patching on the fist and last residues
                AUTOgenerate anlges or dihedrals
        For each residue:
                RESIdue name and total charge specification
                        (or PRESidue if this is a patch)
                ATOM definitions within this residue
                GROUping dividers between atom definitions
                BOND specification
                ANGLe specifications
                DIHEdral angle specifications
                IMPRoper dihedral angle specifications
                CMAP dihedral angle specifications, resolution
                DONOr specifications
                ACCEptor specifications
                IC information
                PATChing residues to use if defaults are not desired
        Closing:
                END statement
        Display control:
                PRINT option

The format above is not rigid. In particular, The 'out of
residue declarations' may be augmented and redefined at any point.
These declarations are checked against all 'out of segment' atom
references. This is done to avoid potential problems where atom names
are misspelled. The number following the declaration is ignored, and
is for the users own reference (or debugging).

The syntax of all subcommands are as follows:

::

     MASS atom-type-code atom-type-name mass

     DECLare out-of-residue-name
      This adds names to be considered for possible connections
      to the previous or next residues.  This is done as a spelling
      check. Any atoms names not contained with in the residue nor
      on this list of declarations will be flagged as an error.
      Use the symbol "-" as an atom name prefix to denote the previous
      residue, use "+" for the subsequent residue. Use "#" as a prefix
      for the (n+2) residue.

     DEFAults [ FIRSt { name } ] [ LAST { name } ]
                      { NONE }          { NONE }

     AUTOgenerate [ ANGLes ] [ DIHEdrals ]
                  [NOANgles] [NODIhedrals]

     { RESIdue  } name [total-charge]
     { PRESidue }
     Residues labled PRES may only be used for patching. Residues
     defined with RESI may not be used as a patch.

     ATOM iupac atom-type-name charge repeat(exclusion-names)

     GROUp
      This keyword divides the structure into specific electrostatic
      groups.  These are used with explicit group-group electrostatic
      options and are used to make the atom-atom list generation
      more efficient.  If a RESIdue does not start with a GROUp command,
      then any ATOMs defined will belong to the last group of the
      previous residue.  Also, the maximum number of atoms allowed in
      any group is currently set at 1000 (MAXING in dimens.fcm).
      As a general guide, and electrostatic group should be roughly neutral
      or have unit charge.  A group should generally be a rigid group of
      atoms, and should not have heavy (non-hydrogen) atoms in a 1-5
      arrangement.  Hydrogens should always be in the same group as its
      bonded partner.  A group should NEVER include two or more groups
      of atoms that are not covalently linked.

     BOND repeat(iupac iupac)

     { ANGLe } repeat(iupac iupac iupac)
     { THETa }

     { DIHEdral } repeat(iupac iupac iupac iupac)
     { PHI      }

     { IMPRoper } repeat(iupac iupac iupac iupac)
     {  IMPHi   }

     { CMAP } repeat(iupac iupac iupac iupac iupac iupac iupac iupac)

     DONOr [ hydrogen ] [ heavy-atom ] [ antecedent-1 antecedent-2 ]
           [ BLNK     ] [ hydrogen   ]

     The antecedents are not required unless hydrogen position
     generation is desired.

     ACCEptor iupac [iupac [iupac] ]
     The first antecedents is required if and angle dependence about
     the acceptor atom is desired. The second antecedent is unused.

     {  IC   }
     { BILD  } name name name name bond angle phi angle bond
     { BUILd }
     BLNK may be used to indicate a missing atom name.

     DELEte   { ATOM             }  iupac  [COMBine iupac]
              { BOND             }  (iupac iupac)
              { THETa | ANGLe    }  (iupac iupac iupac)
              { DIHEdral | PHI   }  (iupac iupac iupac iupac)
              { IMPHi | IMPRoper }  (iupac iupac iupac iupac)
     Deletions are allowed only in patch residues (PRES); the optional
     COMBine keyword for ATOM deletions allows passing part of the IC
     data for the deleted atom to the "combine" atom, i.e. stereochemistry
     of atoms bonded to the deleted atom.  In order to use the COMBine
     option, both atoms must be present in the PSF and it must be invoked
     from the PATCh command (not the GENErate command).

     PATChing [ FIRSt { name } ] [ LAST { name } ]
                      { NONE }          { NONE }

     PRINt { ON  }
           { OFF }


     The PRINt command may be used to control the display of lines as
     they are read by the RTF reader. The initial setting for printing is
     controlled by the READ command itself. If PRINT is specified, then
     printing will initially be enabled; otherwise, the commands will not
     be echoed. PRINT ON turns on echoing of RTF specifications; PRINT OFF
     turns them off. This command is useful for debugging an addition to a
     previously tested topology file.

A small sample RTF card file follows:

::

   *  title for documentation example
   *
      18    1
   MASS     1 H      1.00800
   MASS    11 C     12.01100
   MASS    12 CH1E  13.01900
   MASS    13 CH2E  14.02700
   MASS    14 CH3E  15.03500
   MASS    31 N     14.00670
   MASS    38 NH1   14.00670
   MASS    51 O     15.99940
   MASS    56 OH2   15.99940

   DECL -C
   DECL -O
   DECL +N
   DECL +H
   DECL +CA

   DEFA FIRS NTER LAST CTER

   RESI ALA     0.00000
   GROU
   ATOM N    NH1    -0.35
   ATOM H    H       0.25
   ATOM CA   CH1E    0.10
   GROU
   ATOM CB   CH3E    0.00
   GROU
   ATOM C    C       0.45
   ATOM O    O      -0.45
   BOND N    CA        CA   C         C    +N        C    O         N    H
   BOND CA   CB
   THET -C   N    CA             N    CA   C              CA   C    +N
   THET CA   C    O              O    C    +N             -C   N    H
   THET H    N    CA             N    CA   CB             C    CA   CB
   DIHE -C   N    CA   C         N    CA   C    +N        CA   C    +N   +CA
   IMPH N    -C   CA   H         C    CA   +N   O         CA   N    C    CB
   CMAP -C  N  CA  C   N  CA  C  +N
   DONO H    N    -C   CA
   ACCE O C
   BILD -C   CA   *N   H      0.0000    0.00  180.00    0.00   0.0000
   BILD -C   N    CA   C      0.0000    0.00  180.00    0.00   0.0000
   BILD N    CA   C    +N     0.0000    0.00  180.00    0.00   0.0000
   BILD +N   CA   *C   O      0.0000    0.00  180.00    0.00   0.0000
   BILD CA   C    +N   +CA    0.0000    0.00  180.00    0.00   0.0000
   BILD N    C    *CA  CB     0.0000    0.00  120.00    0.00   0.0000


   RESI OH2     0.00000
   GROUP
   ATOM OH2  OH2    -0.40000     H1   H2
   ATOM H1   H       0.20000     H2
   ATOM H2   H       0.20000
   BOND OH2  H1        OH2  H2
   THET H1   OH2  H2
   DONO H1   OH2 -O -O
   DONO H2   OH2 -O -O
   ACCE OH2
   PATC FIRS NONE LAST NONE

   END

.. note::

   The use of improper dihedrals for the PSF is unrelated
   to the use of improper dihedrals for the internal coordinate tables.

   :: 
   
                         L
      PSF usage:         |
                         |
                         I
                        / \
                       /   \
                 -----J---- K------ 


      IC table usage: 

                    I      L
                     \    /
                      \  /
                       *K
                        |
                        |
                        J

Note that for PSF usage the first atom is the central atom,
and the last atom is the atom to be restrained relative to
the axis defined by the middle pair of atoms.  For the IC table
usage, the central atom is in the third position, but the 
axis is again defined by the middle pair of atoms. 


Reading data other than the sequence or coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The parameter files (PARA) and internal coordinate files (IC)
and hydrogen bond (HBONd) data files can be read as card images or binary
files. Specifying CARD signifies card image input; specifying FILE
signifies binary file input. Please note that topology file must be read
in before the parameters can be read.

Protein structure files (PSF) files and non bonded lists (NBONd)
can only be read as binary files. The constraints (CONStraint) which
includes dihedral restraints may only be read as formatted file (card).

There are two types of IC card files (residue number vs. resid's).
The residue number option is the default, and atom assignments are based
on residue number. This is the low precision form. The resid option
is the high precision form and atom assignments are based on SIGID's and
RESID's. This is also useful where different homologies are used.

The Image file (IMAGes) containing transformation information can
only be read in card image format (see :doc:`images`).
The INIT keyword will remove all existing image data. Without the
INIT keyword, any existing image items (such as bonds) would be kept.
This allows one to modify the crystal geometry without the necessity
of regenerating all image items.

The TABLe file contains the nonbond energy lookup information.
Once read in, The effects cannot be reversed. The nonbond energy
evaluation is now under control of the table routines.


.. index:: io; write
.. _io_write:

Writes Data Structures to External Files
----------------------------------------

::

         WRITe { { PSF        } [FILE]              }  UNIT unit-number |
                                                       NAME filename
               {                [CARD]    [XPLOr]   }
               { { RTF        }                     }
               { { PARAmeter  }                     }
               { { NBONd      }*                     }
               { { TABLe      }                     }
               {                                    }
               { { COORdinate coor-spec } [CARD]    }
               {                          [PDB [MODEL int [FIRSt|LAST]] [OFFI]}
               {                          [DUMB]    }
               { { IC  [RESId]  [SAVEd] } [FILE]    }
               { { HBONd [ANAL]         }           }
               {                                    }
               { { IMAGes imag-spec} [CARD]         }
               { { ENERgy             }             }
               { { CONStraint [PSF 0] }             }
               { { TITLe              }             }

          title

               { NAMD FILE "filename" **            }
      ** no other options are available

        coor-spec:== [COMP]  [OFFS int] [IMAGes]  atom-selection

        imag-spec::= [ TRANsformations ] [ FORCes ] [ PSF ]

      *: The NBOND list can only be WRITten in binary (FILE) form. Use PRINt to get
         formatted output.


Function
^^^^^^^^

The primary purpose of this command to save some of CHARMM's
data structures. The coordinate and internal coordinate data structures
can be written in formatted form so that they be edited independent
of CHARMM using a text editor. The option, FILE, specifies that a file
is to be written in unformatted form (binary).  The option, CARD,
specifies that a file is to written in formatted form.  For the
coordinate and internal coordinate file, CARD is the default.  The
coordinate option PDB gives a file in Protein Data Bank format, with
just the ATOM records; the MODEL N option writes a PDB file in the NMR-style
multiple coordinate set format (note that for this to work the file has to 
be specified as UNIT <int>, not as NAME <string>): 

============================= ============================================================
MODEL 0 (or no MODEL keyword) just write standard PDB file
MODEL 1                       writes beginning of multicoordinate file (title, MODEL 1,
                              coor, TER, ENDMDL)
MODEL N (N>1)                 appends just coordinates for MODEL N  (MODEL N, coor,
                              TER, ENDMDL)
MODEL N (N<0)                 appends last coordinate set, and END (MODEL \|N\|, coor,
                              TER, ENDMDL, END)
============================= ============================================================                              

Keyword FIRSt forces writing of title even if N.NE.1, LAST forces
writing of END line. 

The XPLOr option of WRITe PSF produces an XPLOR style PSF file (atom
names are used instead of atom numbers) 

The selection of "PSF 0" in the WRITe CONS only works with PERT and
writes data for the lambda=0 PSF. 

A set of title lines must follow the WRIT command. This title
will be written at the start of the file and serves to document the
file. For your protection, one should always make good use of this
title, as it may be the only documentation for the file.

The UNIT keyword specifies what Fortran unit the output should be
written to. It cannot be omitted unless the filename is provided with the NAME
keyword.


.. index:: io; print
.. _io_print:

Writes information to output file (unit 6)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

         PRINt { PSF         [XPLOr]  }
               { RTF                  }
               { CONStraint [PSF 0]   }
               { PARAmeter   [USED]   }
               { RESIdue              }
               { COORdinate coor-spec }
               { IC        [ SAVEd ]  }
               { HBONd       [ ANAL ] }
               { NBOND                }
               { IMAGes imag-spec     }
               { TITLe                }
               { ENERgy               }

        coor-spec::= [COMP] [OFFS int] [IMAGes] atom-selection

        imag-spec::= [ TRANsformations ] [ FORCes ] [ PSF ]


        Syntactic ordering: All commands must be typed in the order shown.

Function
^^^^^^^^

This command is used to list information contained in data
structures used by the program. The information must already have been
created through use of a READ, GENE, HBON, etc., command. The printable
output is sent to unit 6.

The XPLOr option of PRINt PSF produces an XPLOR type PSF
listing.  Atom names are printed instead of atom numbers.

The selection of "PSF 0" in the PRINt CONS only works with PERT
and prints data for the lambda=0 PSF.

For printing parameters, the USED option causes the print of only
the parameters that were used in the most recent energy evaluation. This
option is PSF dependent.

For hydrogen bonds, ANAL gives a geometrical and energy analysis
of the hydrogen bonds. Representing the hydrogen bond as
A2-A1-X-H....Y-, the distances X-Y, H-Y, the angle (180 - <X-H-Y ), the
dihedral angle A2-A1-X-H and the hydrogen bond energy contribution are
listed. A more versatile hbond analysis facility is provided by 
COOR HBOND (see :doc:`corman`).


.. index:: io; titles
.. _io_titles:

Specifying and manipulating titles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Titles are optional. All title lines MUST begin with a "*".
If no title is specified, the title will be untouched. This is useful when
a series of titles are needed. Titles are terminated with a line containing
only a "*" in the first colunm. There may be up to 32 lines contained
in any title.

The titles are read using RDCMND, thus parameter substitutions are
allowed.

A command TITLe has been added to CHARMM which can be used to specify
a title to be used by subsequent write commands.

For interactive use, A title is always required (no backspace can
be done) when RDTITL is called.

The date,time, and user is added at the end of the title when
a title is written to a file. If a date and time is already present,
it will be superseded. For the print option, the date and time
information is left as it was.

A second title array TITLEB has been added to CTITLA.FCM
TITLEA is to be used for writing, and TITLEB must be used for reading
from data files. In this way, the main title is never destroyed by reading
a data file. For any write command, TITLEA can be modified by specifying
a title. Any further writes will use that title, unless a new title is
specified.

As it is now, title lines should not end in "-" and any characters
beyond a "!" will not be included in the title.

Titles may begin with a "#" as well as "*". The pound sign
is converted to a "*" upon reading. When the first title line begins
with "#", the old title is not destroyed. All entered title lines
supersede any previous title lines. Obviously, if more title lines are entered
than were previously present, then there will be no difference in the two
methods. This option was added for cases where a series of identical
titles, except for a different first line, was needed.

The COPY keyword of the TITLe command will copy the current TITLB
(the reading title) to TITLA (the writing title) before reading the
subsequent title. If there is no subsequent title, then just a copy is done.

Normally, when titles are written to card files, the first column
"*"s are retained. With the WRITe TITLe command, several changes are made.
First, the first colunm of "*"s is suppressed. Second, no date and time
and username is added. Third, the file is not closed. This command is
primarily used for creating files for plotting. It is often used in
conjunction with looping and energy terms. Here is an example of possible
applications;

::

   OPEN WRITE CARD UNIT 23 NAME ENERGY.DAT        ! Open the file for plot data
   WRITE TITLE UNIT 23
   *  this file contains .....
   *  more message data  .....
   *
   SET 1 -180.0                   ! Set the initial dihedral angle value
   LABEL LOOP                     ! Here is the loop return point
   CONS DIHE ....... MIN @1       ! Introduce the desired dihedral constraint
   MINIMIZE .....                 ! Minimize
   CONS CLDH                      ! Remove the dihedral constraint
   SET 2 @1                       ! copy parameter one to parameter two
   TRIM 2 FROM 1 TO 10            ! Pad parameter two with blanks for formatting
                                  ! It will now be 10 characters long
   WRITE TITLE UNIT 23
   * DIHEDRAL = @2 ENERGY = ?ENER ! write this only this line to unit 23
   *
   INCREMENT 1 BY 15.0            ! Add 15 to parameter one
   IF 1 LT 180.1 GOTO LOOP
   CLOSE UNIT 23
   STOP
