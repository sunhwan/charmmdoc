.. py:module:: struct

============================================
Generation and Manipulation of the Structure
============================================

The commands described in this node are used to construct and
manipulate the PSF, the central data structure in CHARMM (see 
PSF.FCM).  The PSF holds lists giving every bond, bond angle, torsion
angle, and improper torsion angle as well as information needed to
generate the hydrogen bonds and the non-bonded list. It is essential
for the calculation of the energy of the system. A separate data
structure deals with symmetric images of the atoms.  See :doc:`images`.

There is an order with which commands to generate and manipulate
the PSF must be given.  First, segments in the PSF must be generated one
at a time.  Prior to generating any segments, one must first have read a
residue topology file, see :ref:`Read <io_read>`.  To
generate one segment, one must first read in a sequence using the READ
command, see :ref:`Sequence <io_sequence>`. Then, the :chm:`GENERATE`
command must be given.

Once a segment is generated, it may be manipulated. This can
be done in a very general way using the patch command. The patch
command allows, for instance, the addition of disulfide bridges,
changing the protonation state of a titratible residue or to make a
histidine heme crosslink.

The PSF can be saved with the "WRITE PSF" command.  A PSF may be
read with the "READ PSF" command.  The "READ PSF" command has an "APPEnd"
option that allows the merging of individual PSF files.  In addition, the
"DELETE" command allows the deletetion of atoms and all references to the
deleted atoms.

.. _struct_generate:

The Generate Command - Construct a Segment of the PSF
-----------------------------------------------------

::

   GENErate [segid] { generate-spec        } [SETUp]
                    {  DUPLicate segid     }

   generate-spec::= [FIRSt pres] [LAST pres] [WARN] [ ANGLe   ] [ DIHEdrals ]
                                                    [ NOANgle ] [ NODIhedral]


This command uses the sequence of residues specified in the last
:chm:`READ SEQUuence` command and the information stored in the residue
topology file to add the next segment to the PSF. Each segment contains a
list of all the bonds, angles, dihedral angles, and improper torsions
needed to calculate the energy. It also assigns charges to all the
atoms, sets up the nonbonded exclusions list, and specifies hydrogen
bond donors and acceptors. Any internal coordinate which references
atoms outside the range of the segment is deleted. This prevents any
unexpected bonding of segments.

The FIRSt and LAST specifications define what patch-residues
should be used for the terminating residues. If no specification is given,
then the default patching as specified in the topology file will be used.

The WARN keyword, will list all elements that were deleted due
to nonexistant atoms (usually references to the terminating residues).

The SETUp option will cause any internal coordinate table entries
(IC) from the topology file to be appended to the main IC table.

The ANGLe (NOANgle) and DIHEdral (NODIhedral) options override the
autogeneration option specified in the topology files. This may be done
to suppress unwanted additional terms, or to add terms for specific
residues.

.. note::

   The solvent residues (TIP3, ST2, WAT) must be generated
   with the NOANgle and/or NODIhedral qualifier. This is only necessary for 
   the files which use the AUTOgenerate ANGLes and/or DIHEdrals as a
   default. This also means that a protein residue sequence and
   water molecules may not be combined in the same generate command.
   Also, there is a special "READ SEQUence residue_type integer" command where
   integer is the number of resudies of residue_type (often water molecules).
   This avoids the need to list the number of residues followed by the
   specification of each TIP3 residue name individually as is done with a
   protein.

For the DUPLicate segment option, the generate command MUST NOT
be preceded by a READ SEQUence command. This option will create a new
segment which is identical (except for the segid) to an existing segment.
This option is mainly intended for the use in setting up small crystals
for viewing and other analysis.


.. _struct_nbx:

Nonbonded Exclusion List (NBX)
------------------------------

Some pairs of atoms are excluded from the nbond exclusion lists
because their interactions are described by other terms in the hamiltonian.
By default directly bonded atoms and the 1-3 atoms of an angle are excluded
from the nonbond calculation.  In addition the diagonal interactions of
the six membered rings in tyrosine and phenylalanine were excluded from
the nonbond calculation through CHARMM version 15 with RTOPH6. Hydrogen
bonds, and dihedral 1-4 interactions are not excluded (note that other
workers may differ from us on one or both of these points).

The list of nonbonded exclusion is generated in two steps.  First
a preliminary list is made at generation by GENIC using any information
that may be present in the topology file (for example, diagonal
interactions in rings).  The second step is an automatic compilation of
all the bond and angle interactions, followed by a sorting of the list,
performed in MAKINB.  The list is stored in the linked list pair IBLO14/INB14,
where IBLO14(i) points to the last exclusion in INB14 to atom i.  If the list
is modified after MAKINB, then either MAKINB should be called again to
resort the list, or care must be taken to see that the INB14 list is ascending
with all INB14 entries having higher atom numbers than i and that all atoms
have at least one INB entry.

MAKINB is called by default after any operation which changes
internal coordinates such as generate, patch, or edit.

The exclusion list can be specified in three ways. First, interactions
that are to be excluded can be placed in the topology file by listing
the excluded atoms after the charge.  Second,
NBXM mode can be specified as a qualifier to any of the commands which
change internal coordinates.  Third, the default NBXM value can be specified
in the parameter file.  The NBXM values and actions are (in the
following "include" refers to what is being kept (included) in the
exclusion list):

        ======== ===========================================================
        0        use the existing list (do nothing)
        1 or -1  include nothing extra
        2 or -2  include only 1-2 (bond) interactions
        3 or -3  also include 1-3 (angle) interactions
        4 or -4  also include 1-4 interactions automatically.
        5 or -5  include up to 1-3 interactions as exclusions and process
                 1-4 interactions using the 1-4 van der Waal parameters and
                 reduced elecrostatics (E14FAC).
        ======== ===========================================================

Negative values suppress the use of the information present in
the topology file.  Positive values add to the information that was in
the topology file. 


.. _struct_patch:

Patch command to modify PSF
---------------------------

Syntax (command level)

::

        PATCh <pres-name> segid1 resid1 [, segid2 resid2 [,...
                                         [, segid9 resid9]...]]
                                          [SORT]
                                           [SETUp]
                                            [WARN]

Syntax (corresponding patch residue in RTF)

::

        PRES <pres-name>

        [GROUp]
        [ATOM  <I><atomname>  <parameter type>   <charge> ]
        [DELEte ATOM <I><atomname>]

        [ [DELEte] BOND <I1> <I2> ]
        [ [DELEte] ANGLe <I1> <I2> <I3> ]
        [ [DELEte] DIHEdral <I1> <I2> <I3> <I4> ]
        [ [DELEte] IMPRoper <I1> <I2> <I3> <I4> ]
        [ [DELEte] DONOr  [<I1>] <I2> [[<I3> [<I4>]] ]
        [ [DELEte] ACCEptor  <I1> [ <I2> [ <I3> ]] ]

        [ IC  <I1> <I2> [*]<I3> <I4>   real real real real real ]
        [ DELEte IC <I1> <I2> [*]<I3> <I4> ]

     where I1, I2, I3, I4 refer to <I><atomname>.


Rules governing the patch procedure:

1) If an atom is being added via a PATCH at least one or more atoms
   already existing in the residue to which the patch is being added
   must be included in the PRES with an ATOM statement.  Unless 
   this(these) atoms are deleted using the DELEte ATOM command
   internal terms associated with this atom which are already present 
   in the residue should NOT be included in the PRES.

2) if no <I> is specified before <atomname> the patch procedure assumes
   that the atom should be in residue (segid1 resid1).

3) a '-', '+', '#' as a first letter in <atomname> tries to locate or add
   the atom <atomname> in the previous, next, next of the next, residue
   of residue (segid<I> resid<I>), respectively.

4) GROUP brackets in a patch residue have highest priority.

5) If no GROUP is specified, the group numbers of referenced, already
   existing atoms remain unchanged. Added atoms are placed in the last group 
   of the referenced residue.

6) A GROUP statement in a patch residue CAN enclose atoms in different
   referenced residues. However, if there is a conflict between
   sequential residue AND group boundaries new residues MIGHT be created
   with resid's and segid's referring to the referenced residues.
   These cases are indicated by a message from MAPIC that a negative number
   of residues were created. The user has to check the PSF explicitly
   to decide whether the modifications done by PATCH are appropriate.

7) Along with the PSF the coordinates, comparision coordinates, harmonic
   constraints, fixed atom list, internal coordinates (IC) are
   mapped correctly.

8) THERE IS NO MAP OF NBONDS, HBONDS, SHAKE, DYNAMICS ETC.
   THE ATOMNUMBERS ARE CHANGED.

9) Any bond, angle, etc referring to deleted atoms is itself deleted.
   The bond, angle, etc lists are compressed.

10) Even if the AUTOgenerate ANGLe and/or DIHEdral option has been
    invoked new angles and/or dihedrals have to be included in
    the PRES when that particular patch is being called after
    the GENErate statement.  The angles and/or dihedrals will
    be generated automatically for any patch which is called
    in the GENErate statement following the FIRSt or LAST 
    statements. NOTE: If angles and dihedrals are present in
    a PRES which is called in a GENErate statement in which
    AUTOgenerate ANGLes and/or DIHEdrals is being used those
    angles and/or dihedrals will be invoked twice in the PSF
    and, thus, be included twice when the energy is calculated.
 
    The AUTOgenerate command (next) can be used to circumvent the above
    problems, and removes the need for specifying angles and dihedrals
    as part of a PRES definition.


.. _struct_autogen:

Completely autogenerate all angles and/or dihedrals
---------------------------------------------------

::

   AUTOgen   {  ANGLes     [ DIHEdrals ]  }
             {  DIHEdrals  [ ANGLes    ]  }

Sets the angle and/or dihedral counts to zero in the PSF, and rebuilds the
indicated list(s) of energy terms.  Intended to simplify the development of
patches, since only bonding terms need to be specified in PRES definitions
which are followed by this command.  Note that at least one keyword is
required, but both may be specified, in either order.

.. warning::

   may be a problem if the PSF contains any water molecules.


.. _struct_delete:

Delete atoms or energy terms in the structure
---------------------------------------------

::

   DELEte  {   ATOMs        atom-selection                 } [SORT]
           {                                               }
           { { BONDs              } double-atom-selection  }
           { { ANGLes             }                        }
           { { DIHEdrals          }                        }
           { { IMPRoper-dihedrals }                        }
           { { CONNectivity       }                        }

The DELEte ATOM option deletes selected atoms and all references
to them in PSF.

.. note::

   THIS WILL CHANGE THE ATOM NUMBERING.

.. note::

   If PERT is currently in use, this command only affects the active
   (lambda=1) PSF.  The reference PSF (lambda=0) is only modified by the PERT
   command.

For the internal energy terms, any entry that has an atom selected
in both atom selections will be deleted. Note, if an atom is selected in
both atom selections, all connections to this atom will be deleted,
except for bonds. For a bond to be deleted, one of its atoms must
appear in each of the atom selections. The CONN (connectivity) option
will delete all bond, angles, dihedrals, and improper dihedrals.
This option avoids the necessity of running the DELEte command four times
when one wishes to break some connectivity.

The SORT option performs an optional sorting of the PSF after the
deleted atoms have been mapped out.

.. _struct_rename:

RENAme - rename portions of the current PSF
-------------------------------------------

RENAme is invoked only from the main command parser and it
includes the working PSF. Its syntax is;

::

        RENAme  { SEGId }  new-name  atom-selection
                { RESId }
                { RESN  }
                { ATOM  }

Any atoms selected will have the corresponding ID modified.
There is a check for duplicate SEGIDs, RESIDs, and atom names, but it
wont stop you if BOMLEV is negative. Renaming ST2 will not change their
status (except in the setup for SHAKE, which will be fixed soon).


.. _struct_join:

Joining Two Adjacent Segments
-----------------------------

For some operations, it is convenient to be able to join
two adjacent segments together. This process has no effect on the
energy terms, but just reorganizes naming and grouping of atoms
into segments. This is especially useful with IMAGES so that
all images in the PSF are identified only as a single segment.

::

   JOIN  first_segment  [second_segment]  [ RENUmber ]
 
The second segment must follow the first sequentially in the
PSF.  There is no checking for duplicate residue identifiers. The
RENUmber option sets the resid for each residue of the composite
segment to the relative index in that segment (just as it would have
during a generate command).  If only a single segment is specified
with the RENUmber option, then the resid's of this segment will be
numbered sequentially.
