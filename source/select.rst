.. py:module:: select

##############
Atom Selection
##############

Atom selection is used for many commands within CHARMM.
Its existence is one of the main factors in the versatility of CHARMM.

.. index:: select; syntax
.. _select_syntax:

Recursive Atom Selection Syntax
-------------------------------

::

   .... SELEction   <factor>   [SHOW]  END ....

Listed in priority order (low to high)
(operators not separated by a blank line are processed sequentially)

::

   <factor>:==  <factor> .OR. <factor>

                <factor> .AND. <factor>

                <factor> .AROUND. <real>
                <factor> .SUBSET. <int*>
                <factor> .SUBSET. <int1> : <int2>

                .NOT. <factor>
                .BONDED. <factor>
                .BYRES. <factor>
                .BYGROUP. <factor>

                (  <factor>  )

                <token>

                <keyname>

     <token>::= SEGId <segid>*
                SEGId <segid1> : <segid2>
                ISEG  <segnum1> : <segnum2>
                RESId <resid>*
                RESId <resid1> : <resid2>
                IRES  <resnum1> : <resnum2>
                RESName <resname>*
                RESName <resn1> : <resn2>
                IGROup  <grpnum1> : <grpnum2>
                TYPE <type>*
                TYPE <type1> : <type2>
                CHEMical <chem>*
                CHEMical <chem1> : <chem2>
                ATOM <segid>* <resid>* <type>*
                PROPerty [ABS]<prop><.LT.|.GT.|.EQ.|.NE.|.GE.|.LE.|.AE.><real>
                POINt <x-coor><y-coor><z-coor> [CUT <rmax>] [PERIodic]
                BYNUnumber <int>*
                BYNU <int1> : <int2>
                INITial
                LONE
                HYDRogen
                USER
                PREVious
                RECAll <integer>
                ALL
                NONE

.. note::

     where '*' allows wildcard specifications:
     
     * \*  matches any string of characters (including none),
     * %  matches any single character,
     * #  matches any string of digits (including none),
     * \+  matches any single digit.


.. _select_double:

Double atom selections
----------------------

Some commands allow (or require) a double atom selection.

::

    command ... SELE first-selection-spec END  [ SELE second-selection-spec END ]

If no atom selection is passed, both selection arrays have all
atoms selected. If only one atom selection is specified, both selection arrays
will contain that atom selection. If an error is found, both atoms selections will have no atoms
selected, and the appropriate error processing will occur.

.. _select_function:

Description of Atom Selection Features
--------------------------------------

.. index:: select; wildcard

Wildcard specification
^^^^^^^^^^^^^^^^^^^^^^

* \*  matches any string of characters (including none),
* %  matches any single character,
* #  matches any string of digits (including none),
* \+  matches any single digit.

.. index:: select; range

Range specification
^^^^^^^^^^^^^^^^^^^

Ranges are indicated by ':' and are defined by the lexigraphical
of :chm:`SEGID`, :chm:`RESName`, :chm:`TYPE`, a combination of numerical
and lexigraphical order for :chm:`RESId` 's (see routine :chm:`SPLITI`) and
by numerical order for :chm:`BYNUmber` specifications.

Keyname option
^^^^^^^^^^^^^^

The user may specify keynames with the :chm:`DEFIne` command (see :doc:`MISCOM`).
Each keyname corresponds to a particular atom selection. Keynames
are processed before tokens, so if there is a naming conflict, the
keyname will prevail. Keynames may not be abreviated.  Whenever
the PSF is modified and the number of atoms changes ALL keynames
are removed.

User specified selection
^^^^^^^^^^^^^^^^^^^^^^^^

:chm:`USER` represents a user selection which should be defined in the
user subroutine USRSEL.

.. index:: select; around

Around atom specification
^^^^^^^^^^^^^^^^^^^^^^^^^

The operation ``<factor> .AROund. <real>`` finds all atoms within a
distance ``<real>`` around the atoms specified in ``<factor>``. For this
operation it must be ``(QCOOR = .TRUE.)`` and all coordinates should
be known. Otherwise a warning message is printed.

.. index:: select; point

Around a spatial point option
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The token ``POINt <xcoor><ycoor><zcoor> CUT <rmax>`` selects all
atoms within a sphere around point ``(x,y,z)`` with radius ``rmax``.
The default value for ``rmax`` is 8.0.  If the keyword :chm:`PERIodic`
is present AND simple periodic boundary conditions are in effect
through the use of the :chm:`MIPB` command, the selection reflects
the appropriate periodic boundaries. (See :ref:`MIPB <images_mipb>`)

.. index:: select; byres

Selection by residues as a whole
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The function ``<factor> .BYRes.`` includes all atoms in a residue
which contains at least one atom selected in ``<factor>``.

Selection based on atom properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The token :chm:`INITial` selects all atoms with known coordinates.
The token :chm:`LONE` selects all lone pairs (based on a mass selection).
The token :chm:`HYDRogen` selects all hydrogens (based on a mass selection).
These selections use the CHARMM intrinsic routines INITIAL, HYDROG
and LONE.

Use of previous selection
^^^^^^^^^^^^^^^^^^^^^^^^^

The token :chm:`PREVious` will start from the current contents of the
atom selection array. This feature only works for commands
where atom selection storage is permanent, and is usually local
to a specific command type.  For example, the :chm:`PREVious` token within
graphics will be the last atom selection previously used in
graphics, even if an atom selection was requested later in a non
graphics application.  At some point, this command will be made to
be consistent, but for now, use it if it works.  If not, the :chm:`DEFINE`
command is better.

.. index:: select; bonded

Selecting based on connectivity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Atoms can be selected based on connectivity.  For example
if one wants to spin a methyl group by 30 degrees, the following
command sequence may be used:

::

    DEFINE CARBON SELE selected-carbon-atom END
    DEFINE TOMOVE SELE TYPE H* .AND .BONDED. CARBON END
    IF ?NSEL .NE. 3 GOTO ERROR
    DEFINE SECOND SELE .NOT. TYPE H* .AND. .BONDED. CARBON END
    IF ?NSEL .NE. 1 GOTO ERROR
    COOR AXIS SELE CARBON END SELE SECOND END
    COOR ROTATE AXIS SELE TOMOVE END PHI 30.0

.. index:: select; subset

Selecting subsets of atoms by index
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Atoms may be selected by index on a subselection.  For example,
if one wants to analyze a set of interesting atoms, the sequence
might be:

::

    DEFINE INTERESTING SELE interesting-atoms END
    SET N ?NSEL
    SET I 0
    LABEL LOOP
    INCREMENT I BY 1
    COOR DIST SELE INTERESTING .SUBSET. @I END SELE ALL END ...
    IF @I .LT. @N GOTO LOOP

.. index:: select; property

Selection based on atom properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The token ``PROPerty <prop>`` selects all atoms which have
the specified property relative to a selected value.
The allowed properties are: ``1`` and ``<keyname>``

* ``1``: The active weighting array (WMAIN or WCOMP)
  (This old construct works only when properties are actually present)
          
* keyname: An array keyname (from the SCALar command syntax)
  The currently allowed keynames include:

  ::
      
    X        Y        Z        WMAIn    XCOMp    YCOMp    ZCOMp    WCOMp   
    DX       DY       DZ       ECONt    EPCOnt   MASS     CHARge   CONStrai
    XREF     YREF     ZREF     FBETa    MOVE     TYPE     IGNOre   ASPValue
    VDWSurfa ALPHa    EFFEct   RADIus   RSCAle   FDIM     FDCOns   FDEQ
    SCA1     SCA2     SCA3     SCA4     SCA5     SCA6     SCA7     SCA8
    SCA9     ZERO     ONE

.. index::
   pair: subst; select

Atom selection substitution parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

During each invocation of the atom-selection, the following substitution
parameters are specified:

  ============== ==============================================================
  :sub:`NSEL`    Number of selected atoms from the most recent atom selection
  :sub:`SELATOM` Atom number of first selected atom
  :sub:`SELCHEM` Chemical type of first selected atom
  :sub:`SELIRES` Residue number of first selected atom
  :sub:`SELISEG` Segment number of first selected atom
  :sub:`SELRESI` Resid of first selected atom
  :sub:`SELRESN` Residue type of first selected atom
  :sub:`SELSEGI` Segid of first selected atom
  :sub:`SELTYPE` Atom name of first selected atom
  ============== ==============================================================

These may be used in any subsequent CHARMM command (NOT in the current command).
These definitions remain valid up to and including the next CHARMM command
that includes an atom selection.  For commands with a double atom selection
the variables are defined by the final atom selection.

.. index::
   pair: examples; select

Some examples of atom selections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

  sele atom * * CA end    !will include all C alphas in the list.

  sele .not. type H* end  !will include all atoms that are not hydrogens.

  sele atom MAIN 1 * .around. 5.0 end
                          !will include all atoms that are within a sphere
                          of radius 5.0 around any atom of the residue
                          MAIN 1.

  sele bynu 1 : 100 end   !will include atoms number 1 to 100.

  sele resid 1 : 10 .and. segid. MAIN -
       .and. .not. ( type H .or. type N .or. type O ) end
                          !will include all the atoms of reside 1 to 10
                          in the segment MAIN except atoms H, N, and O.

  sele bynu 1 .or. bynu 3 .or. bynu 5 .or. bynu 7 .or. bynu 8 -
       bynu 11 .or. bynu 13 .or. bynu 15 .or. segid SOLV end
                         !will include atoms number 1, 3, 5, 7, 8, 11, 13,
                          and 15, and the SOLV segment.

  ! to select side chain atoms in the polygen all atom parameter set
  ! where mb is the myoglobin protein 
  sele segid mb .and. .not. ( type n .or. type ca .or. type c .or. - 
         type o .or. type oct* .or. type ha* .or. type hn .or. type ht* ) end 

  sele prop abs charge .gt. 0.5 end
                  ! select all atoms with charge > 0.5 or charge < -0.5

  sele prop radius .gt. 3.0 end  ! select all large atoms (radius>3.0).

