.. py:module:: intcor

=============================================
The Internal Coordinate Manipulation Commands
=============================================

The commands in this section can be used to construct cartesian
coordinates from internal coordinate values. The internal coordinate
data structure can also be used for analysis purposes.
There are flexible editing commands for manipulating the data structure.
When these commands are used in conjunction with the Coordinate
Manipulation commands (see :doc:`corman`)  and the
I/O commands (see :doc:`io`), a rather complete model
building facility exists.

.. index:: intcor; syntax
.. _intcor_syntax:

Syntax of Internal Coordinates commands
---------------------------------------

::

        IC  { PARAmeters [ALL]                            }
            { FILL  [COMP] [APPEnd] [PREServe] [SAVEd]    }
            { GENErate   [THREe]    atom-selection        }
            { DIFFerences [COMP] [APPEnd] [SCALe real]    }
            { DERIvatives [COMP] [APPEnd] [DELTa real]    }
            { DYNAmics  dynamics-spec                     }
            { EDIT                                        }
            { BUILd   [COMP]  [SAVEd]                     }
            { SEED atom atom atom  [COMP]                 }
            { PURGe     [SAVEd]                           }
            { ADD       [SAVEd]                           }
            { SUBTract  [SAVEd]                           }
            { SCALe scale-spec [SAVEd]                    }
            { RANDom  [ISEEd int] [SAVEd]                 }
            { GAUSsian [ZMIX] UNIT int  atom atom atom    }
            { PUCKer 5x(atom) ANGLe real AMPL real        }
            {                                             }
            { { DELete } { BYNUM int [int]  } [SAVEd]     }
            { { KEEP   } { ic-selection     }             }
            {                                             }
            { SAVE     [PREServe ]                        }
            { RESTore  [OVERwrite]                        }
            {                                             }
            {  READ  [FILE] [APPEnd] UNIT int [SAVEd]     }
            {                                             }
            {  WRITe [FILE] [RESId]  UNIT int [SAVEd]     }
            {                                             }
            {  PRINt   [SAVEd]                            }

   atom::= { residue-number atom-name }
           { segid  resid atom-name   }
           { BYNUm  atom-number       }
           { next-from-atom-selection }

   dynamics-spec::= { [AVERages]   }    [FIRStunit int] [NUNIts int]
                    { FLUCtuations }      [BEGIn int] [STOP int] [NSKIp int]

   ic-selection::= {                                          } atom-selection
                   { [FIRSt] [SECOnd] [THIRd] [FOURth] [IMPR] }
                   {                                   [DIHE] }

   scale-spec ::= [ BOND real ] [ ANGLe real ] [ DIHEdral real ]
   atom-selection ::= see *note select:(chmdoc/select.doc).
   next-from-atom-selection ::= all atoms are selected from a single atom
                                selection in sequential order.


The syntax for the EDIT subcommands are:

::

   { DISTance   atom   atom   real                  [ADD]   }
   { ANGLe      atom   atom   atom   real           [ADD]   }
   { DIHEdral { atom   atom  [*]atom  atom }  real  [ADD]   }
   {          { ICNUmber integer           }                }
   { END                                                    }



.. _intcor_function:

Purpose of the various Internal Coordinate commands
---------------------------------------------------

These commands are used to setup, modify and process the
internal coordinates of the molecule. This operation is very useful
in setting up atom coordinates whenever they are not known. This
occurs when a protein structure is built from scratch or when an
existing structure is modified. The modification can be simply a
conformational change, or a change in the residue sequence through
replacement, insertion, or deletion. Many of these modifications can
be processed within the program as it currently stands. Other more
difficult modifications can be facilitated by editing the internal
coordinate card file by using external programs.

This facility is also useful as an analysis tool. Several support
program use the output from IC tables for conformational analysis
(phi-psi maps, ring pucker, pseudo-rotational angles, solvent structure,...).


Command ordering
^^^^^^^^^^^^^^^^

The Internal Coordinate commands (except EDIT and READ) can only be
used if internal coordinates exist (i.e. if the IC common is filled).
This can only be filled by reading an IC file, or by using the SETUp
keyword in the GENErate or PATCh commands. The information used to setup
is obtained from the residue topology file used in the generation process.

Subcommand interpretation
^^^^^^^^^^^^^^^^^^^^^^^^^

* PARAmeter [ALL] - Fill table with parameter values

  Fill the internal coordinates using standard values from
  the parameter file, unless otherwise specified in the residue topology
  file (see RTF:(IO)Rtf File Formats.). A value of zero for any bond or
  angle (not dihedral) indicates that this value should be obtained
  from the parameters. If the ALL keyword is specified, then all angle
  and bond values will be filled from the parameter set regardless of the
  existing values. Setting bond and angles values to zero with the IC edit
  command makes it possible to selectively use this command.

* FILL [COMP] [APPEnd] [PREServe] - Convert from cartesian to internal coordinates

  Fill the internal coordinate values wherever possible from the
  known atomic coordinates.  IC's for atoms that are not placed are zeroed
  unless the PREServe keyword is specified, in which case the entries are not
  modified.  If the COMP keyword is used, then The alternate coordinate set
  will be used to fill the IC data structure.  The APPEnd option will add
  the current values to the existing values of the table.

  For example, one way to see how the current coordinates match a reference
  ic table:
  
  ::
  
     ic scale bond -1.0 angle -1.0 dihe -1.0
     ic fill append
     print ic


* GENErate [THREe] atom-selection - Generate an IC table from connectivity data

  The IC GENErate command will generate additional IC table entries
  based on the bond list (connectivity data).  This command will not modify the
  existing IC table data entries.  It will attempt to add one IC entry for each
  selected atom, but for isolated molecules, it will generate NSEL-3, because the
  first three entries will be incomplete (and thus ignored).  If the "THREe" 
  keyword is specified, then it will also include any 3 atom (incomplete)
  IC entries from the search procedure.

  The main purpose for this command is to be able to easily produce an
  IC table for internal coordinate analysis.  The simple algorithm employed
  here does not do a very good job of guessing reasonable IC dihedral values
  for use in coordinate building such as;
  
  ::
  
      IC GENErate ...
      IC PARAM
      IC SEED ...
      IC BUILD
      
  For this usage, all dihedrals are set to trans-planar (180.0 degrees) and all
  improper types are set to 180, +/-120, or +/-90 degrees, depending on the number
  of bonds on the central atom.  Some editing of the table (see IC EDIT) may
  be essential before constructing coordinates.  Where the current algorithm
  has known problems:
  
  1. Ring structures - trans-planar is not a good starting guess.  Some may
     need to be edited to be cis-planar (0 degrees) or gauche (+/-60 degrees).

  2. Linear bonds - It may be necessary to add explicit IC entires for square
     planar configurations in order to avoid the linear bond problem with exact
     square planar configurations (i.e. the NA-FE-NC angle in the heme group)
     
  3. The algorithm always uses the most massive branch as the "mainchain".
     If the last residue (c-terminal) of a polypeptide chain is LYS, then
     there will be no -C-N-CA-C (phi) torsion angle, since the sidechain
     is more massive than the carboxyl group (instead you'll have -C-N-CA-CB).
     This can be "fixed" by temporarily setting the masses of the terminal
     atoms to large real values.  Likewise, an atom in a branch to be avoided
     as the "mainchain" can be temporarily set to a large negative value.

  4. Tetrahedral chiral centers - Currently the improper values are set based
     on atom order for non-mainchain atoms.  The opposite chirality can be
     achieved with IC EDIT, or an atom reordering.

.. note::

   The current algorithm (roughly 300 lines of code) is rather simple.  The
   problems listed above could be solved with a more complex method.  Whether this
   is ever done (and when) will depend on user demand - BRB - March 2, 1998


* DIFFerence [COMP] [APPEnd] - Fill table with the difference of two structures

  The DIFF command will cause the IC table entries to be filled
  with differences of internal coordinate values. Normally the
  values are filled (MAIN-COMP), but this is reversed if the COMP
  keyword is used. The APPEnd keyword will cause the differences
  to be added to the existing IC table values.

* DERIvative [COMP] [APPEnd] - Fill table with internal derivatives

  The DERIvative command will fill the IC table entries with the
  analytical internal derivatives associated with a particular vector
  (velocity, forces, or normal mode are typical examples). Normally, it is
  assumed that the vector is stored in the main coordinate set and the
  coordinates are stored in the comparison set. If the COMP keyword is
  specified, then their roles are reversed. The APPE keyword will cause
  the new values to be added to the existing table values.

* DYNAmics - Fill table with dynamic averages or fluctuations.

  The IC DYNAmics command generates averages or fluctuations for
  the IC table from a dynamics trajectory. The syntax is;

* IC DYNAmics  

  :: 
  
     IC DYNAmics  { [AVERages]   }     [FIRStunit int] [NUNIts int]
                  { FLUCtuations }        [BEGIn int] [STOP int] [NSKIp int]

  Either the averages, or the fluctuations about the current table values
  can be computed. The sequence;

  ::
  
        IC FILL
        IC DYNAmics AVERage ...
        PRINT IC
        IC DYNAMics FLUCtuations ...
        PRINT IC

  will print out the averages and fluctuations about the averages. For
  dihedrals, whether computing fluctuations or averages, a reference value
  is subtracted before summing (i.e. values are always within 180 degrees
  of the reference value), thus explaining the need for the IC FILL command
  preceding the first IC DYNAmics command.

* EDIT - Add to or modify the IC table elements

  Edit the internal coordinate file. This command causes the
  input stream to transfer to the IC edit mode. The edit mode
  commands are:
  
  ::
  
                DIST atom atom real               [ADD]
                ANGLE atom atom atom real         [ADD]
                DIHE atom atom [*]atom atom real  [ADD]
                END

         atom::= {residue-number atom-name}
                 { segid  resid atom-name }
                 { BYNUm  atom-number     }

  These commands will specify a particular internal coordinate value.
  All occurrences of the specified item will be modified.
  If the specified atoms have no corresponding IC table entry,
  then a new IC entry will be added for these specified atoms.
  For the ANGLe option when a new IC entry is added, the corresponding
  1-2 and 2-3 distances will be filled from other existing values
  (or left as zeros). For the DIHEdral option, an optional '*' on the third
  atom denotes that this is the central atom of an improper dihedral type.
  When adding a new IC entry for dihedrals, the associated bond and angle
  terms are filled from existing table values is possible, otherwise,
  they are added with zeros.
  
  The ADD option will add the specified value to the current
  corresponding value in the IC table.  An error will be issued if the
  IC table entry does not already exist.  For example, the command
  
  ::
  
      DIHE 15 N 15 CA 15 C 16 N  10.0  ADD
      
  will increase the psi angle of residue 15 by 10.0 degrees.

  The END command is used to exit from the edit IC mode.

* BUILd [COMP] [SAVEd] - Convert from internal to cartesian coordinates

  This command determines the cartesian coordinates for all
  unspecified atoms from the data in the IC file (wherever possible).
  The user is responsible to make sure that the designation for all atoms
  is unique. In the case that the system is over specified, An atom is
  placed on the first opportunity (no checking is done for currently placed
  atoms). If it is desired to modify the position of atoms with known
  coordinates, the coordinates for those atoms must be reinitialized using
  the COOR INIT command. If an IC element contains a zero bond length or
  angle (not dihedral), then it will not be used to place the terminal atom.
  This option is useful in cases where the system is over-specified and
  building is not desired for some IC's. For example;
  
  ::
  
      IC:  2 O4' 2 C2' 2 C1' 2 H1' 0.0 0.0 120.0 109.5 1.0
      
  can be used to place H1' but will not place atom O4'. Again, if the COMP
  keyword is used, then the alternate coordinate set will be used and modified.
  If the "SAVEd" keyword is used, then it will use the IC table generated
  by the most recent "IC SAVE" command in lieu of the normal IC table.

* SEED atom atom atom [COMP] - Place first three atoms for building reference

  When the cartesian coordinates are not specified for any atoms,
  the BUILd command cannot be used to generate positions since all positions
  are determined relative to known positions. The SEED command specifies the
  positions of the three atoms . It puts the first at the origin, the second
  on the x-axis, and the third in the xy-plane. The three atoms must have
  entries in the IC file corresponding to: dist 1-2, angle 1-2-3, dist 2-3.
  The COMP keyword causes the alternate coordinate set to be modified.

* DELEte

  :: 
  
      DELEte { BYNUM int [int] } - Delete selected elements from the table
             { ic-selection    }                                          

      ic-selection::= {                                          } atom-selection
                      { [FIRSt] [SECOnd] [THIRd] [FOURth] [IMPR] }
                      {                                   [DIHE] }

  This commands deletes a specified set of IC's from the data file.
  The delete can be by number (using the BYNUM keyword and a range), or by
  atoms selection. Any IC that contains a selected atom will be removed.
  By default, the atom can match in any position. However, a specific match
  may be requested by specifying one or more of (FIRSt,SECOnd,THIRd,FOURth).
  Specifying all of them is equivalent to the default.  Also, the keywords
  DIHE or IMPR may be used to select to delete only those which represent
  normal dihedrals (DIHE), or those of type improper (IMPR).
  
  ::
  
        IC DELE DIHE FIRST SELE TYPE CA END - This command will delete any IC
                                              element that is a dihedral type AND
                                              has a CA in the first position.
  
* KEEP

  ::
  
    KEEP { BYNUM int [int] } - Delete all non-selected elements from the table
         { ic-selection    }

    ic-selection::= {                                          } atom-selection
                    { [FIRSt] [SECOnd] [THIRd] [FOURth] [IMPR] }
                    {                                   [DIHE] }

  The keep command is almost the logical opposite of the DELEte
  command.  Its options are identical, except that the selected set of IC's
  is kept, and all of the remaining ones are deleted. As in the IC DELEte
  command, a positional match may be selected.  Also, the keywords DIHE or
  IMPR may be used to ADDITIONALLY keep all entries which represent normal
  dihedrals (DIHE), or those of type improper (IMPR).
  
  ::
  
        IC KEEP DIHE FIRST SELE TYPE CA END - This command will retain any IC
                                              element that is a dihedral type OR
                                              any improper type of IC element
                                              with a CA in the first position.

* PURGe - Clean up the IC table

  The PURGe command will cause all IC's that contain undefined atoms
  to be deleted. This is not automatic because sometimes it is desirable
  to keep partial IC table entries (where less than 4 atoms are defined).

* SCALe [BOND real] [ANGLe real] [DIHE real] - Scale table elements by a factor.

  The SCALe command will multiply all elements of a table by a
  constant factor. This is primarily used when the table contains IC
  differences or derivatives, and new structures are to be generated based
  on these values.  For example, the following sequence will generate
  a structure that is midway between two structures in internal
  coordinate space (Note: this is different from COOR AVERage):
  
  ::
  
        IC FILL COMP
        IC FILL APPEND
        IC SCALE BOND 0.5 ANGLE 0.5 DIHE 0.5
        COOR INIT SELE ALL
        IC SEED ...
        IC BUILD
        COOR ORIE RMS MASS

* RANDom [ISEEd int]- Randomize all dihedral values

  The RANDom command will randomize all dihedral values in the
  table.  It will use and modify the specified ISEED value.  To randomize
  a subset of dihedral values, the following procedure may be optimal:
  
  ::

      [read in RTF, PARAM, and sequence]
      GENErate MAIN   SETUp ! generate segment with IC table
      IC PARAM   ! fill zeroes in IC table with optimal parameter values
      IC SAVE    ! save the entire IC table
      IC DELEte IMPRoper ! get rid of improper IC terms
      IC DELEte .... ! get rid of dihedrals that will not be randomized
      IC KEEP   .... !  "
      IC PRINT    ! print to check if correct dihedrals are randomized
      IC RANDOM ISEED 12345678 ! Randomize all remaining IC dihedral values
      IC RESTore PREServe ! add all IC entries that are not randomized
      IC SEED ... ! start the build process
      IC BUILD    ! complete the build process
      COOR ORIENT MASS ! align the center of mass with the origin

* SAVE 

  :: 
  
    SAVE [PREServe ] - Save the current IC table
         [OVERwrite]

  The SAVE command will copy the main IC table to a second
  IC table for later retrieval.  If the PREServe keyword is specified,
  then any IC elements already in the second table will be unmodified
  and only IC entries from the main set not contained in the second
  set will be appended to the second set.
  
  If the OVERwrite keyword is specified, then all entries from the main
  set will be copied to the second set, however, any IC elements
  already in the second table and not contained in the main table will
  be unmodified.
  
  If neither PRESeve nor OVERwrite is specified, then a simple copy
  of the main set to the second set is performed (current second set
  data is lost).

* RESTore

  :: 
  
     RESTore [PREServe ] - Restore a previously saved IC table
             [OVERwrite]

  The RESTore command will copy the saved IC table (See IC SAVE
  command) to the current IC table.  If the PREServe keyword is specified,
  then any IC elements already in the main table will be unmodified
  and only IC entries from the second set not contained in the main
  set will be appended to the main set.

  If the OVERwrite keyword is specified, then all entries from the second
  set will be copied to the main set, however, any IC elements
  already in the main table and not contained in the second table will
  be unmodified.
  
  If neither PRESeve nor OVERwrite is specified, then a simple copy
  of the second set to the main set is performed (current main set
  data is lost).

* GAUSsian [ZMIX] UNIT int atom atom atom - Make a GAUSSIAN86 input file from 
  CHARMM coordinates

  The GAUSsian command will make a GAUSSIAN98 coordinate in Z-matrix 
  form for use with the popular ab initio program. The MAIN coordinates will 
  be used unless the COMP keyword is specified.  The first three atoms must 
  be specified (in the IC SEED format) and an output unit number must be 
  specified for a write access file. If internal coordinates are undefined
  the entire molecule will be printed in Z-matrix form. In case internal 
  coordinates are defined for a portion of molecule then this will 
  result in Z-matrix printed for these atoms only. The remaining atoms will
  not be shown.
  
  Mixed Cartesian-internal coordinate printing is activated by ZMIX keyword. 
  Atoms having internal coordinate records will be printed in Z-matrix form 
  while the other atoms will be printed in their Cartesian representation.
  The mixed Cartesian-internal coordinate representation is designed for
  performing partial geometry optimization of mutual orientation of two 
  interacting fragments (molecules) taken at their fixed geometry. Atoms of 
  the first fragment are defined by Cartesian coordinates. Atoms of the second 
  fragment are defined by CHARMM internal coordinates and these atoms will be 
  printed in Z-matrix form relatively to the first fragment.

* PUCKer - Set ring pucker values

  Set the ring pucker magnitude and phase IC table to a specified value.
  To force the conformation of the ring to these values, use "CONS IC.." and
  then minimize.
  
  To set the pucker for a DNA sugar, use the following command;

  ::
  
      IC PUCKer <I> c4' <I> o4' <I> c1' <I> c2' <I> c3' ANGLe <real> AMPLitude <real>

  Where <I> is the residue number (may be a loop variable) or SEGID/RESID
  combination (see :ref:`Syntax <intcor_syntax>`).  For this command to
  function without warning and/or errors, all 5 canonical torsion angles
  corresponding to the selected atoms must be part of the current IC table.

  .. note::
  
     The order of the atoms is significant!  


.. _intcor_structure:

Internal Coordinate concepts
----------------------------

Given the positions of any three atoms, the position of a fourth
atom can be defined in relative terms (internal coordinates) with three
values: a distance, an angle, and a dihedral specification. Where many
atoms are connected in a long sequence (as in proteins) it is easiest
to consider four atoms in a chain. If the positions of one end of the chain
is known, it is possible to find the positions of all of the remaining atoms
with a series of internal coordinate values. But in the more general case,
where some central portion of a molecule is known it is necessary to be able
work in both directions. This lead to the present form of the internal
coordinate data structure (five values for four atoms) where if either
endpoint is unknown and the other three atoms are determined, the position
of the end atom can be found.

The improper type of internal coordinate data structure was created
for branching structures (as opposed to simple chains). Since there are
roughly five values in the data structure for every atom it is clear that
the positions are over-specified. Keep this in mind when externally editing
IC files. The program will use the first acceptable value when building
a structure and ignore any redundancies. The EDIT commands will always
modify all occurrences of each edited parameter.

Normal IC table entry:

::

                I
                 \
                  \
                   J----K
                         \
                          \
                           L
        values (Rij),(Tijk),(Pijkl),(Tjkl),(Rkl)

Improper type of IC table entry

::

                I        L
                 \     /
                  \   /
                   *K
                   |
                   |
                   J
        values (Rik),(Tikj),(Pijkl),T(jkl),(Rkl)


Internal Coordinate file structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The internal coordinate file can be stored in either card or binary
form. for most purposes the card form will be used (since it can be edited).
There are two types of elements in the internal coordinate file, those that
correspond to normal dihedral angles and those that correspond to improper
dihedrals. They can be distinguished by the presence of a '*' just before
the IUPAC name of the third (K) atom (its presence denotes an improper
dihedral type).

For each element there are four atoms (referred to as I,J,K,L) and
five values. Elements of the IC file are symmetric with respect to
inverting the order of the atoms except that for improper types only atoms
I and L can be interchanged (also the sign of phi must be changed since
phi(IJKL)=-phi(LJKI)  ).
