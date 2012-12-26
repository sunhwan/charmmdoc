.. py:module:: corman

====================================
The Coordinate Manipulation Commands
====================================

The commands in this section are primarily used for moving
some or all of the atoms. There is a wide range of commands and options.
All of the commands may be used on either the main coordinate set, or
the comparison set. Some commands require both sets of coordinates.

.. index:: crystal; syntax
.. _corman_syntax:

Syntax of Coordinate Manipulation commands
------------------------------------------

::

   COORdinates { INITialize                       } [COMP] [DIMS] [atom-selection]
               { COPY                             }   [WEIGhting_array]
               { SWAP                             }   [IMAGes] [SECOnd]
               { AVERage  [ FACT real ]           }
               { SCALe    [ FACT real ]           }
               { MASS_weighting                   }
               { ADD                              }
               { SET  vector-spec                 }
               { TRANslate vector-spec            }
               { ROTAte vector-spec {PHI real}    }
               {                    {MATRix}      }
               { TWISt  vector-spec   RATE real   }
               { ORIEnt [MASS] [RMS] [NOROtation] }
               { RMS    [MASS]                    }
               { TMSCore                          }
               { UFSR                             }
               { DIFFerence                       }
               { FORCe  [MASS]                    }
               { SHAKe  [MASS]                    }
               { DRAW      draw-spec              }
               { DISTance  distance-spec [DIFF]   }
               { DIPOle [OXYZ] [MASS]             }
               { MINDist   distance-spec          }
               { MAXDist   distance-spec          }
               { READ  io-specification           }
               { WRITe io-specification           }
               { PRINt io-specification           }
               { RGYR [MASS] [FACT <real>]        }
               { OPERate image_name               }
               { STATistics [MASS]                }
               { VOLUme    {SPACe integer}        }
               {                                  }
               { DUPLicate { 2X(atom-selection) } }
               {           { PREVious           } }

   COORdinates HISTogram { X } [IUNIt int]  HMIN real HMAX real HNUM integer  -
                         { Y }  [HSAVe] [HPRInt] [HNORm real] [HDENsity real] -
                         { Z }   [COMP] [WEIGhting_array] atom_selection
                         { R }

   COORdinates { HBONd   }  [CUTHB <real>] [CUTHA <real>] [IUNIt <int>]  -
               { CONTact }  [BRIDge <resnam>] [VERBose] [TCUT real] -
                              2X(atom-selection) traj-spec -
                            [IRHI <int> [DRH <real> ][RHMAx <real>] ] -
                            [ITHI <int> [DTH <real> ][THMAx <real>] ] -
                            [PBC [CUBIC|TO|RHDO BOXL|XSIZE <real> -
                                   [YSIZE <real> [ZSIZE <real>] ] ]]

   COORdinates SECStructure [first-selection [second-selection]] -
                            [QUIEt | VERBose] [CUTH real] [CUTA real]

   COORdinates  DYNAmics  [COMParison]  [PAX]  [atom-selection] [NOPRint] -
                          traj-spec   [ORIENT [MASS] [atom-selection] ]

   COORdinates  PAXAnalysis [COMParison]  [atom-selection] [NOPRint]  [SAVE] -
                traj-spec

   COORdinates  SEARch { search-spec               } disposition-spec
                       { INVErt                    }
                       { KEEP xvalue yvalue zvalue }
                       { EXTEnd  RBUFf real        }

    search-spec ::   [atom-selection] [COMP] [IMAGe] [operation-spec]
                       [XMIN real] [XMAX real] [XGRId integer]
                         [YMIN real] [YMAX real] [YGRId integer]
                           [ZMIN real] [ZMAX real] [ZGRId integer]

    operation-spec ::=  {              }  { [VACUum] }  { [RESEt] }
                        { [RCUT  real] }  {  FILLed  }  {  AND    }
                        { [RBUFf real] }  {  HOLES   }  {  OR     }
                                                        {  XOR    }
                                                        {  ADD    }

    disposition-spec::= { [NOPRint]      } [NOSAve] [CREAte segid CHEM type]
                        {PRINt [UNIT int]} [ SAVE ]


   COORdinates   SURFace  [atom-selection] [WEIGhting] {  CONTact-area   }
                            [ACCUracy real]           { ACCEssible-area }
                               [RPRObe real]


   COORdinates   CONVert-from/to-unit-cell [ from | to ] -
                 [atom-selection] [COMP] [IMAGe] -
                 a  b  c   alpha   beta  gamma

                 [ from | to ] ::= [ FRACtional | SYMMetric | ALIGned ]


   COORdinates   AXIS  atom-selection [atom-selection] [MASS] [COMP] [IMAGEs]

   COORdinates   LSQP  [ NORM  ] [VERBose] [MASS] [COMP] [IMAGEs] [WEIGh] -
                       [ MAJOr ]
                       [ MINOr ]
                                 atom-selection

   COORdinates COVAriance traj-spec 2x(atom_selection) [UNIT_for_output int] -
                          [RESIdue_average_nsets integer] [MATRix] -
                          [ENTRopy [TEMP <real>] [DIAG] [RESI] [SCHL] ]

   COORDinates DMAT -
          [RESIdue_averaging] [NOE_weighting] [SINGle_coordinate_file] -
          [CUTOff <real>] [UNIT_for_output <int>] [TRAJectory] [CUTOff <real>] -
          [PROJect UPRJ <int>] [PROBability UPRB <int>] [TOLE <real>] MKPRoj -
          traj-spec 2x(atom_selection) [ [RELAtive] RMSF [DUNIt <int>]] [MATRix]

   COORdinates PUCKer [SEGId segid] RESId resid1 [TO resid2] [AS | CP]

   COORdinates HELIx atom-selection [atom-selection]

   COORdinate ANALysis {WATer} [RLP <int>] <atom-selection>  -
     {XREF <real> YREF <real> ZREF <real>} -  ! setup arbitrary analysis point
     {CROSs|SITE [MULTI] <atom-selection>} -  ! setup solute analysis site or
                                              ! cross terms for arbitrary solvent
     traj-spec -                              ! reading trajectories
     NCORs <int> RSPIn <real> RSPOut <real> - ! MSD/IVAC set-up
     RSPHere <real> DR <real>  MGN <int> -    ! g(r) setup
     RDSP <real> -                            ! cutoff for DENS,KIRK and DBF
     DENS <real> -                            ! userspecified bulk density
                                              ! (atoms/A**3)
                                              ! for normalization of g(r)
     {IMSD <unit>|IVAC <unit>} IDENs <unit> - ! output for  MSD, VAC and DENsity
     {IGDISt <unit> [IHH <unit>] [IOH <unit>]|ISDISt <unit>} - ! g(r) requests
       {BYGRoup|BYREsidue|BYSEgment}          ! discard distances WITHIN
                                              ! specified unit for g(r)
     IMRD                       ! Magnetic Relaxation Dispersion analysis
         RRES  cutoff radius for calculation of residence time. if 0 use shell
               beteween RSPIN, RSPOUT

     IKIRkg <unit> -                   ! Kirkwood g-factor (dipole correlations)
     RKIRk               ! distance dependent Kirkwood factor for water
                         ! iff a SITE MULTI selection containing
                         ! at least two atoms is
             given, then a unit-vector pointing from the first to
             the second site atoms  will be used in the
             scalar product with a unit vector along the water dipoles
     NKIRk   number of points in r-dimension for IKIR and RKIR
             from r=0 to r=RDSP

     XBOX <real> YBOX <real> ZBOX <real> - !PBC info for analysis
     IFDBF <unit> IFDT <unit>  RCUT <real> ZP0 <real> NZP <int> - ! DBF analysis
     IHIST <unit> IPDB <unit> [XMIN <real> XMAX <real> DX <real>] - !3D histogram
                              [YMIN <real> YMAX <real> DY <real>] -
                              [ZMIN <real> ZMAX <real> DZ <real>] -
                              [WEIGht] [CHARge] [DIPOle] -
                              [THREshold <real>] [NORM <real>] -
      IDIP <unit> [MIND <real>] [MAXD <real>] [NUMD <int>] -
                                                   ! dipole distribution
      EXVC <atom-selection> MCP <int> MCSH <int> - ! EXcludedVolumeCorrection
      RPRObe <real> ISEEd [WEIG] -

      RCOR <integer> -                 ! Rotational Correlation Time Analysis
      ROUT <unit>  TLOW <real>  TUP <real>  MAXT <integer> -

      IHYDn <integer>  RHYD <real>     ! Hydration numner

   COORdinates INERtia [atom-selection] -
                    [ENTRopy [TEMPerature <real>] [SIGMa <real>] ] -
                    [STANdard <SOLUtion|GAS>]

   COORdinates CONFormational { <resname> } [ PRINT ] [ READ io-speficication ] -
                    [atom-selection] [COMP]

   COORdinates PATH { NREP <int> } {NAME <character*>} [<PDB|FILE|UNFO|CARD|FORM>]

   atom-selection:== (see *note select:(chmdoc/select.doc).)

   distance-spec::=
          {  WEIGhting vector-spec               atom-selection            }
          {                                                                }
          { [UNIT int] [CUT real] [ENERGy [CLOSe]] 2X(atom-selection) -    }

                   { [Nonbonds] } { [NO14exclusions] } { [NOEXclusions] }  -
                   { NONOnbonds } {    14EXclusions  } {    EXCLusions  }

                [TRIAngle]   [ HISTogram HMIN real HMAX real HNUM integer  -
                                [HSAVe] [HPRInt] [HNORm real] [HDENsity real] ]


   vector-spec::= {  [XDIR real] [YDIR real] [ZDIR real]  } [DISTance real]
                    [XCEN real] [YCEN real] [ZCEN real]       [FACTor real]
                 {  AXIS                                 }

   draw-spec::= [DFACt real] [NOMO]  UNIT integer

   io-specification:== (see *note io:(chmdoc/io.doc).)

   traj-spec::= [FIRSt int] [NUNIts int] [NSKIp int] [BEGIn int] [STOP int]

.. corman_simple

Descriptions of the simple coordinate manipulation commands
-----------------------------------------------------------

All of these commands allow either the main coordinate set (default),
or the comparison set (:chm:`COMP` keyword) to be modified. The other coordinate
set is only changed by the :chm:`SWAP` command and the :chm:`ORIEnt RMS` command when
the specified atoms are not centered about the origin.

The DIMS coordinate set (DIMS keyword) is used with the
DIMS command (:doc:`dims`) and it is mainly used with COPY
to load the target structure:  'COOR COPY DIMS'. The DIMS set also works with
ORIENT, PRINT, and STAT, but not with any other operations. Copy the
DIMS set to the comparison set ('COOR COPY DIMS COMP') if other
operations on the target structure are required.

Each of these commands may also operate on a subset of the full
atom space. The selection specification should be at the end of the command.
The default atom selection includes all atoms.

If the :chm:`IMAGes` keyword is specified, then the operation will be
performed on the image atoms as well (if images are present).

The :chm:`SECOnd` keyword specifies that the second comparison set be used.
This keyword can be used with any command that uses a comparison set (e.g.
:chm:`COPY COOR COMP SECOnd` to copy coordinates to the second comparison set;
:chm:`COPY COOR SECOnd` to copy the coordinates from the second to the main set).
Use of this command requires compilation with the COMP2 precompiler keyword.


The INITialize command
^^^^^^^^^^^^^^^^^^^^^^

The :chm:`INITialize` command returns the coordinate values of the
specified atoms to their start up values (9999.0). The main use of
this command is in connection with the :chm:`IC BUILD` command, which may
only find coordinates for atoms with the initial value.


The COPY command
^^^^^^^^^^^^^^^^

The :chm:`COPY` command will copy the coordinate values into the
specified set FROM the other coordinate set.


The SWAP command
^^^^^^^^^^^^^^^^

The :chm:`SWAP` command will cause the coordinate values of the
specified atoms to be swapped with the comparison set.


the AVERage command
^^^^^^^^^^^^^^^^^^^

The :chm:`AVERage` command will generate a new coordinate set at a
point along the displacement vector between the present coordinate set
and the other set. The :chm:`FACTor` value determines the relative step along
this vector. Its default value is 0.5 (a true average). A :chm:`FACTor` value of
1.0 is equivalent to the copy command. Negative or greater than unit
positive values are also allowed.


The SCALe command
^^^^^^^^^^^^^^^^^

The :chm:`SCALe` command will cause the coordinate values for all
selected values to be scaled by a required scale factor. This option
is designed to work with coordinate displacement vectors. A scale
factor of zero will set the selected coordinate values to zero.
This option may also be useful in plotting.


The MASS_weighting command
^^^^^^^^^^^^^^^^^^^^^^^^^^

The :chm:`MASS_weighting` command will cause all selected coordinates
to be scaled by the mass of each atom. If the :chm:`WEIGht` option is specified,
the weighting array will be scaled.


The ADD command
^^^^^^^^^^^^^^^

The add command will add the main and the comparison
coordinate values and store the results in the selected coordinate set.
As with other commands, only selected atoms will be modified. If
an atom in either set is undefined, then the sum will also be undefined.
This option is designed for use in cases where one or both coordinate
sets contain coordinate displacement vectors.


The SET command
^^^^^^^^^^^^^^^

The :chm:`SET` command will set all coordinate values of selected
atoms to a specified value determined by the vector specified. This is
a simple manner in which to zero a coordinate set with the command;

::

        COOR SET XDIR 1.0 DIST 0.0

Note, the :chm:`XDIR` keyword value was included so that the vector has a nonzero
norm (required for all vector specifications).


The TRANslate command
^^^^^^^^^^^^^^^^^^^^^

The :chm:`TRANslate` command will cause the coordinate values of
the specified atoms to be translated. The translation step may be
specified by either X, Y, and Z displacements, or by a distance along
the specified vector. When no distance is specified, The :chm:`XDIR`, :chm:`YDIR`, and
:chm:`ZDIR` values will be the step vector. If the :chm:`AXIS` keyword is used, then
the translation will be along the axis defined by the previous :ref:`COOR AXIS <corman_coor_axis>`
command. For this option, a distance may be specified, but if it isn't,
then the translation distance will be the :ref:`COOR AXIS <corman_coor_axis>` vector length


The ROTAte command
^^^^^^^^^^^^^^^^^^

The :chm:`ROTAte` command will cause the specified atoms to be rotated
about the specified axis vector through the specified center. The vector
need not be normalized, but it must have a non-zero length. If the :chm:`AXIS`
keyword is used, then the axis and center information from the last
:chm:`COORdinates AXIS` command will be used. The :chm:`PHI` value gives the amount
of rotation about this axis in degrees.
Only the atoms specified will be rotated. If the :chm:`MATRix` keyword is used
the rotation will be made using an explicit rotation matrix, input in
free format on the three following lines (3 real numbers /line):

::

    U(1,1) U(1,2) U(1,3)
    U(2,1) U(2,2) U(2,3)
    U(3,1) U(3,2) U(3,3)

.. note::

   This command uses a LEFT HAND sense, not the usual right hand rule...
   It was a mistake, but this is kept for historical reasons (numerous scripts).
   The left hand sense is consistent with dihedral angles (i.e. if you define a
   vector along bond A-B (from A to B) and then rotate B (and its bonds) by a
   positive angle (in the left hand sense), then the dihedral angles will
   increase.  Other rotation angles in CHARMM (should) use the regular
   right hand rule (except for the :chm:`COOR TWISt` command).


The TWISt command
^^^^^^^^^^^^^^^^^

The :chm:`TWISt` command will cause the specified atoms to be rotated
about the specified axis vector through the specified center. The vector
need not be normalized, but it must have a non-zero length. If the :chm:`AXIS`
keyword is used, then the axis and center information from the last
:chm:`COORdinates AXIS` command will be used. The amount of rotation will depend
on the projected distance of the atom on the axis multiplied by the :chm:`RATE`
value (in degrees).

This command was designed to generate helical structures that are more or
less twisted than an initial helical structure.  This is an easy way to
homogeneously perturb a helix.  I can be also used to induce a twist in
planar structures.

.. note::

  this command uses a left handed sense, not the usual right hand rule...
  (see :chm:`ROTAte` above).


The ORIEnt command
^^^^^^^^^^^^^^^^^^

The :chm:`ORIEnt` command will modify the coordinate values of ALL of
the atoms. The select set of atoms is first centered about the origin,
and then rotated to either align with the axis, or the other coordinate set.
The :chm:`RMS` keyword will use the other coordinate set as a rotation reference.
The :chm:`MASS` keyword cause a mass weighting to be done. This will
align the specified atoms along their moments of inertia. When the :chm:`RMS`
keyword is not used, then the structure is rotated so that its principle
geometric axis coincides with the X-axis and the next largest coincides
with the Y-axis. This command is primarily used for preparing a
structure for graphics and viewing. It can also be used for finding
RMS differences, and in conjunction with the vibrational analysis.

The :chm:`NOROtation` keyword will suppress rotations. In this case,
only one coordinate set will be modified.


The RMS command
^^^^^^^^^^^^^^^

The :chm:`RMS` command will compute the RMS or mass weighted :chm:`RMS`
coordinate differences between the selected set of atoms just as they
lie. This differences from the :chm:`COOR ORIENT RMS` command in that no coordinate
modifications are made and no translation is done.


The DIFF command
^^^^^^^^^^^^^^^^

The :chm:`DIFF` command will compute the differences between the main
and comparison set (or the reverse) and store this difference in the
modified coordinate set. Undefined or unselected atoms result in a zero.
If the :chm:`WEIGht` keyword is invoked, then the WCOMP array is subtracted from
WMAIN and the coordinates are untouched.


The FORCe command
^^^^^^^^^^^^^^^^^

The :chm:`FORCe` command will copy the current forces (DX,DY,DZ)
of the selected atoms to the specified coordinate set. Atoms not selected
are given a value of zero. If the :chm:`MASS` keyword is specified, then the
forces will be divided by the mass. This would correspond to an
acceleration in dynamics.


The SHAKe command
^^^^^^^^^^^^^^^^^

This command will :chm:`SHAKE` the selected coordinate set with respect
to the other (as a reference). A mass weighting may be used. Any atoms
that are not selected are considered to be fixed (infinite mass).
In order to use this command, the :chm:`SHAKe` command must first be invoked
which sets up the shake constraints.

Lone pairs (:doc:`lonepair`) with undefined coordinates can be built
by :chm:`COOR SHAKE`.


The DIPOle command
^^^^^^^^^^^^^^^^^^

Calculates the dipole moment of selected atoms. If total charge
is not zero, the dipole moment is somewhat ill-defined and coordinate system
dependent; in this case the center of geometry of the selected atoms is used
as origin for the coordinate system in which the dipole moment is calculated.
This can be altered by the :chm:`MASS` keyword. If it is present the center of mass
will be used as origin of the relative coordinate system.

For the purpose of compatibility with Gaussian program this feature can be
disabled by adding :chm:`OXYZ` keyword, which forces calculation of dipole moment
relatively to the origin of Cartesian coordinate system.

Prints out dipole moment cartesian components and magnitude (in Debyes) and
the total charge. CHARMM variables :sub:`CHARGE`, :sub:`XDIP`, :sub:`YDIP`, :sub:`ZDIP` and :sub:`RDIP` (charge, x,y,z and magnitude of dipole) are set.

The UFSR command
^^^^^^^^^^^^^^^^

Compare two structures (working set versus comparison set)
with the Ultra Fast Shape Recognition algorithm by Ballester and
Richards (:ref:`Ballester 2007 <dims_references>`). This
algorithm is intended to differentiate two structures based on atomic
distributions. Notice that in this approach the score is normalized
and a value of 1 means two identical structures. The current
implementation is identical to the one proposed in their paper.

.. corman_function

Descriptions of the remaining corman commands
---------------------------------------------

See the descriptions of the simple commands for some background
information on these commands.

The DISTance command
^^^^^^^^^^^^^^^^^^^^

The :chm:`COOR DIST` command will either find distances between atoms
or the distances of atoms from a fixed point in space (:chm:`WEIGh` option).
This command can find distances within a single coordinate set, or
find distances between atoms in two coordinate sets (:chm:`DIFF` option).

The :chm:`DISTance` command can find all atom distances between two
atom selections. A unit number may be specified (default=6) and a
cutoff distance may be included as well (default=8999.0). If no selection
is specified, all atoms will be included! The delimiter :chm:`END`
must separate the two sets of atom selections. The van der Waal energy
may be requested with the :chm:`ENERgy` keyword, and if this option is used,
the list of pairs with a positive van der Waal energy may be selected
with the :chm:`CLOSe` keyword (i.e. only close contacts will be listed).
The :chm:`NEAR` option will list the nearest atom in the second atom selection
to the atoms in the first selection.

The :chm:`COOR DISTance` command doesn't gives distances between
excluded atoms unless the :chm:`EXCLusions` keyword is specified. This make
it much easier to search for bad contacts. Likewise, 1-4 interactions and
other interactions may be requested or omitted.

The command;

::

         COOR DISTance ENERgy CLOSe CUT 5.0 SELE ALL END SELE ALL END -
                14EXclusions NONBonds

will list all atom pairs that have a positive van der Waal energy.

The command;

::

         COOR DISTance ENERGY CUT 5.0 NONONbonds NOEXclusions 14EXCLusions -
                SELE ALL END SELE ALL END

will list all 1-4 interactions and energies (and nothing else).

The command;

::

         COOR DISTance ENERgy CUT 4.5 SELE RESID 23 END SELE ALL END

will list all contacts less than 4.5A that residue 23 has with the rest of
the system without considering 1-4 interactions or excluded pairs.

The 1-4 vdw terms, E14FAC, and EPS values other than 1.0 are recognized.

The :chm:`WEIGht` option puts the distance of all selected atoms from some
specified point. If no point is specified, then the origin is used. This
is most useful in computing magnitudes of forces or coordinate differences.
For example, the sequence;

::

        ENERGY ...
        COOR FORCE COMP  ! copy forces to the comparison coordinates
        COOR DIST WEIGH COMP  ! put magnitudes in the weighting array.
        PRINT COOR COMP SELE PROP WCOMP .GT. 5.0 END
           ! print atoms with large forces.
           ! Note that all operations were done on the comparison set.

The :chm:`DIFF` keyword causes the selection to work on different coordinate
sets, where the first selection corresponds to the set specified (:chm:`MAIN` or
:chm:`COMP`), and the second atom selection uses the other coordinate set.

The :chm:`HISTogram` option allows a histogram of distances to be produced.
With the histogram, the :chm:`HMIN` and :chm:`HMAX` (the range of the histogram in angstroms)
and the :chm:`HNUM` (the number of bins) must be specified.  The :chm:`HSAVe` keyword causes
the histogram values to be saved for subsequent :chm:`COOR DIST` commands.  In a loop,
this allows g(r) to be calculated from a dynamics trajectory.  The :chm:`HPRInt`
option will cause the final histogram values to be printed.  The :chm:`HNORm` value
will be used to normalize the histogram before printing (divide by :chm:`HNORm`).

A density value, :chm:`HDENS`, is also required, which is the number of selected
objects divided by the volume per object.  Also note: In order to get
this to work with with the crystal facility, the first atom selection
(in the loop) should only include primary atoms, and the second atom
selection should include both primary and image atoms.
The histogram will be scaled by the reciprocal of the distance squared

The histogram will also be scaled by the reciprocal of the distance squared
(to get normalized g(r) plots).  Three columns of numbers are output;
(1) the bin midpoint distance, (2) the normalized g(r), and (3) the total
number of pairs within the bin divided by the :chm:`HNORM` value.
A :chm:`PRNLEV` less than 5 will suppress the listing of distance pairs.
Example of use to get a distance distribution plot:

::

      update imgfrq 20 cutim 20.0
      traj ....
      prnlev 4
      set 1 1
      label loop
      traj read
      update inbf 0 IMALL cutim 10.5
      coor dist image sele segid main .and. type OH2 end sele type OH2 end -
             cut 10.5   HIST HMIN 0.0 HMAX 10.0 HNUM 50 HSAVE
      incr 1 by 1
      if 1 .lt. 1000.5 goto loop

      calc dens = 216.0/30.0  !  #waters/(volume/water)
      coor dist sele none end sele none end -
            cut 10.5  HIST HMIN 0.0 HMAX 10.0 HNUM 50 HNORM 1000.0 -
            HPRINT  HDENS @dens


The RGYR command
^^^^^^^^^^^^^^^^

The :chm:`RGYR` command can compute the Radius of GYRation, center-of-mass
and total mass of the specified atoms. By default the :chm:`RGYR`, uses a unit
weighting factor providing the rms distance from the center of geometry.
The current keywords are:

      ===========   ================================================================
      :chm:`MASS`   use mass weighting (otherwise use unit weight per selected atom)
      :chm:`WEIG`   use a weight array (WMAIN or WCOMP) for the weighting
      :chm:`FACT`   constant to be subtracted from each weight
      ===========   ================================================================

The weight arrays can be filled, by using :chm:`COOR` or :chm:`SCALAR` commands,
before invoking the :chm:`RGYR` routine. In this way almost any :chm:`RGYR` can be computed.


The LSQP command
^^^^^^^^^^^^^^^^

The :chm:`LSQP` command computes the least-squares-plane through the
selected atoms. Weighting can be done by the atom masses [:chm:`MASS`], by
the weighting array [:chm:`WEIG`], or not at all (default). Output is the
equation for the plane, the sum-of-squared distances (weighted) from
the plane (SSQ), and the center-of-mass of the selected atoms.

The keyword :chm:`VERBose` causes some additional output, most useful of
which is the distance from the plane for each atom.

The options; :chm:`NORM`, :chm:`MAJOr`, and :chm:`MINOr` select which vector is
stored as the :chm:`AXIS` (see :ref:`COOR AXIS <corman_coor_axis>` command for more details).  The default
is to not set the :chm:`AXIS` variables.


The OPERate command
^^^^^^^^^^^^^^^^^^^

The :chm:`OPERate` command processes the selected coordinates through
the image transformation specified by name. This command may only be
used if an image file has been read.  The image_name is one of the
image transformation names (:ref:`WRITE IMAGE TRANS <images_write>`).  This is also the SEGID
of the image atoms created by the image update procedure.


The MINDistance command
^^^^^^^^^^^^^^^^^^^^^^^

The :chm:`MINDistance` command computes the minimum distance between
selected coordinates. Usually this command is executed with a double
selection.  Note that the default distance-spec excludes bonded atoms and 1-4
interactions.  If only one selection is given, then it will give the minimum
distance of the selected coordinates between the MAIN and COMP set.


The MAXDistance command
^^^^^^^^^^^^^^^^^^^^^^^

The :chm:`MAXDistance` command computes the maximum distance between
selected coordinates. This command is executed with a double selection.


The STATistics command
^^^^^^^^^^^^^^^^^^^^^^

The :chm:`STATistics` command will print some simple statistics
regarding the selected atoms. The values :sub:`XMIN`, :sub:`YMAX`, :sub:`XAVE`,
:sub:`YMIN`, :sub:`YMAX`, :sub:`YAVE`, :sub:`ZMIN`, :sub:`ZMAX`, :sub:`ZAVE`,
:sub:`WMIN`, :sub:`WMAX`, :sub:`WAVE` are set when this command is executed. These
variable values may then be used un subsequent commands with the "?" symbol.
For example, the command sequence may be used to shift a structure so that
a single atom is in the X-Y plane (e.g. shift in the z-direction);

::

   COOR STATistics SELE desired-atom END
   COOR TRANS  ZDIR ?ZAVE  FACT -1.0

The :chm:`MASS` option will place the average values at the center of mass.


.. _corman_coor_axis:

The AXIS command
^^^^^^^^^^^^^^^^

The :chm:`AXIS` command generates a vector and saves it for subsequent use
for either command parsing, or for use as input in the :chm:`COOR SET`, :chm:`COOR ROTAte`,
:chm:`COOR TRANslate`, or :chm:`COOR DISTance WEIGhting` commands by using the :chm:`AXIS` keyword.
There are two modes for the :chm:`AXIS` command. With a single atom selection, the
stored vector is the defined from the origin to the center of geometry/mass
of all selected atoms. With two atom selections, the vector spans from the
center of the first set of selected atoms to the center of the second.
The :chm:`MASS` keyword invokes the usage of the center of mass.
The :chm:`AXIS` command sets the variables :sub:`XAXIs`, :sub:`YAXIs`, :sub:`ZAXIs`, :sub:`RAXIs`, :sub:`XCEN`, :sub:`YCEN`,
and :sub:`ZCEN`, which may be accessed with the "?" symbol. These values define
the actual vector, the length of the vector, and the center of the vector
(midpoint). For example, to use the distance between two atoms as a
criterion to terminating a run, the following command sequence could be used;

::

   SET 1  10.0
   COOR AXIS SELE first-atom END SELE second-atom END
   IF  1 GT ?RAXIs   STOP

For another example, to rotate the chi-1 torsion of a
specified residue BY 30 degrees, the command sequence would be appropriate;

::

   DEFINE BACK SELE TYPE O .OR. TYPE N .OR. TYPE H .OR. TYPE CA .OR. TYPE C END
   COOR AXIS SELE ATOM MAIN 23 CA END  SELE MAIN 23 CB END
   COOR ROTATE AXIS PHI 30.0  SELE RESID 23 .AND. .NOT. BACK END


The DUPLicate command
^^^^^^^^^^^^^^^^^^^^^

The :chm:`DUPLicate` command copies coordinates between atoms within
a structure.  The coordinates are copied FROM the first selection TO the
second selection. If the selections overlap, watch out!. The matching is
done by number within the selected coordinate sets. If the two selection
have a different number of atoms, a warning will be issued, and the smaller
number will be used. For example, if one needs to compute the relative
orientation between two alpha helices, the following input might be used;

::

   COOR COPY COMP
   COOR DUPL COMP SELE backbone of first END SELE backbone of second END
   COOR ORIE RMS MASS COMP SELE backbone of second END

This will give the RMS shift between these helices as well as the
coordinate transformation required to map one into the other.

The :chm:`PREVious` option may be used with a single atom selection.
This assigns the coordinate position of selected atoms to the value
of the previous atom (by number). This has been used with the command;

::

        COOR DUPLicate PREVious SELE TYPE H* END

to assign hydrogen atom positions to that of the associated heavy atom.

The :chm:`COMP` keyword causes only the comparison coordinates to be used and
modified.  Otherwise, the entire operation involves only the main coordinates.


The DYNAmics command
^^^^^^^^^^^^^^^^^^^^

The :chm:`COOR DYNAmics` command will read a (set of) dynamics trajectory
files and compute the average coordinates (stored in the selected
coordinate set) and the isotropic RMS fluctuations (stored in the weighting
array). The first unit number (:chm:`FIRSt`)(default 51), number of units (:chm:`NUNIts`)
(default 1), frequency of accepted coordinate sets (:chm:`NSKIp`)(default 1),
starting set (:chm:`BEGIn`)(default first set), last set (:chm:`STOP`)(default last set),
may be specified. Option values are not remembered with subsequent
:chm:`COOR DYNA` commands.  The :chm:`NOPRint` suppresses much of the output.
If the keyword :chm:`ORIENT` is present, all coordinate frames will be
RMS re-oriented with respect to the :chm:`COMParison` set (must be defined);
if the word :chm:`MASS` is also there the coordinates will be mass-weighted for
re-orientation; if a second atom selection is provided, only those selected
atoms will be used.

The :chm:`PAX` command causes the principal axis of the motion of each atom
to be computed and save.  The print out gives the direction and magnitude
of the fluctuation as well as the anisotropies.  The PAX data is saved for
a subsequent :chm:`COOR PAXAnal` command if further analysis is desired.


The PAXAnal command
^^^^^^^^^^^^^^^^^^^

The :chm:`COOR PAXAnal` command computes additional data regarding the
principal axis data (computed by the most recent :chm:`COOR DYNA PAX` command).
The trajectory must be reopened and reread, or a different trajectory
may be substituted.  This command prints data for each selected atom and
averages over the selected atoms.  The printout includes the skew and
kurtosis, anisotropies, as well as all of the low moments of the motion.
The :chm:`SAVE` option causes the PAX data structure (from the :chm:`COOR DYNA PAX` command)
to be saved (for subsequent :chm:`COOR PAXA` commands).


The SEARch command
^^^^^^^^^^^^^^^^^^

::

   COORdinates  SEARch { search-spec               } disposition-spec
                       { INVErt                    }
                       { KEEP xvalue yvalue zvalue }
                       { EXTEnd  RBUFf real        }

     search-spec ::   [atom-selection] [COMP] [IMAGe] [operation-spec]
                        [XMIN real] [XMAX real] [XGRId integer]
                          [YMIN real] [YMAX real] [YGRId integer]
                            [ZMIN real] [ZMAX real] [ZGRId integer]

     operation-spec ::=  {              }  { [VACUum] }  { [RESEt] }
                         { [RCUT  real] }  {  FILLed  }  {  AND    }
                         { [RBUFf real] }  {  HOLES   }  {  OR     }
                                                         {  XOR    }
                                                         {  ADD    }

     disposition-spec::= { [NOPRint]      } [NOSAve] [CREAte segid CHEM type]
                         {PRINt [UNIT int]} [ SAVE ]

The :chm:`SEARch` command generates and/or manipulates a grid of small volume
elements.

The :chm:`SEARch` command will search through a set of grid points
for vacuum space points (i.e. points outside the van der Waal radius of
any atom). In the default mode (:chm:`NOPRint`), only the relative volume of filled
and vacuum points are printed concerning the selected atoms.
The grid specifiers must be input (min, max, and grid) for each dimension.
(grid implies number of grid points. Hence

::

        XMIN -10.0 XMAX 10.0 XGRID 41

implies a half Angstrom sampling along the x direction)

The :chm:`FILLed` option will cause non-vacuum points to be listed or plotted.
The :chm:`PRINt` option will cause all found grid points to be listed on the
output unit specified (default 6).

For this command, the atom sizes (radii) are taken from the weighting
array.  To get van der Waal radii into the weighting array, the command;

::

        SCALar WMAIn = RADIus

may be used. If a hole big enough to stuff a water into is to be found,
then the command sequence;

::

        SCALar WMAIn = RADIus
        SCALAR WMAIN ADD 1.6
        SCALAR WMAIN MULT 0.85

would be probably the best to use.

If the :chm:`RCUT` or :chm:`RBUFf` value is set to a nonzero value, then the accessible
volume command is enabled.  When :chm:`RCUT` is set, this is the maximum radius.
When :chm:`RBUFf` is set, then the maximum radius is the weighting array plus the
:chm:`RBUFf` value.  The weighting array is returned with the fraction of free volume
in the shell from the atom radius to the maximum radius.

If the :chm:`HOLEs` keyword is set, only the grid points not connected to the
first point (point in the negative corner of the box) are considered.
In this way, the volume of just the holes can be analyzed and saved.

The :chm:`ADD` option for the :chm:`COOR SEARCH` command has been added to allow
the calculation of partial occupancy factors.  This allow holes in proteins
to be analyzed for flexibility and variability.

It is possible to use multiple :chm:`COOR SEARch` commands and to use boolean
operations to combine the results.  For example, the script sequence;

::

   COORdinates   SEARch  IMAGe -
         XMIN -10.0 XMAX 10.0 XGRId 20 -
         YMIN -10.0 YMAX 10.0 YGRId 20 -
         ZMIN -10.0 ZMAX 10.0 ZGRId 20 -
         NOPRINT VACUUM  SAVE
   ....
   SCALAR WMAIN ...
   ....
   COORdinates   SEARch  IMAGe -
         XMIN -10.0 XMAX 10.0 XGRId 20 -
         YMIN -10.0 YMAX 10.0 YGRId 20 -
         ZMIN -10.0 ZMAX 10.0 ZGRId 20 -
         AND PRINT UNIT 22  RBUFF 2.0 FILLED  NOSAVE

Note, the results of these two commands are computed and the
intersection (AND) is printed.  The first command needs a ":chm:`SAVE`" in order
for the results to be saved.  Also, the grids (if specified) must exactly match
(same number of grid points in all dimensions) for this operation to work.
The :chm:`COOR SEARch` command allocates space, if needed, and frees the space when
the :chm:`NOSAve` option is used.  Thus, if four :chm:`COOR SEARch` commands are needed for a
single computation, the first must have the :chm:`SAVE` option.  The only way
to free the space allocated by the :chm:`COOR SEARch SAVE` command is to run another
:chm:`COOR SEARch` command with the :chm:`NOSAve` option.

If the :chm:`CREAte` option is used then the specified grid points will be
added to the PSF as dummy atoms.  The chemical type of the dummy atom must
be specified and it must be present in the current RTF.  This option can be
used for graphics or for other hole analysis (shape,...).  This option
will add one segment to the PSF, one residue and atoms and groups equal to
the number of selected grid points.


The VOLUme command
^^^^^^^^^^^^^^^^^^

The :chm:`VOLUme` command will compute the volume of a selected set of
atoms.  Its operation is the same as that of the SEARch command, except
that only the volume is printed and the degree of exposure for each atom
is returned in the weighting array.  The SCALAR storage arrays must be filled
before using this command.  The first storage array [1] must contain
the radii of each atom (RMIN) and the second storage array must contain the
outer probe distance (RMAX) for each atom.  The free volume within the RMIN
to RMAX range and not within RMIN of any other atom will be returned in the
weighting array as a ratio of the maximum possible value.  For example a
completely exposed atom will return a value of 1.0 and an atom in the interior
of a protein would return a value of 0.0.  The :chm:`HOLEs` keyword feature
causes holes within the selected atoms to be filled before computing
the total volume and the accessible volume.

:chm:`SPACE` is a maximum number of cubic pixels
i.e. :chm:`SPACE` = :math:`x_{points} \times y_{points} \times z_{points}`
Larger :chm:`SPACE` value results in more accurate calculation but it takes more
memory an computer time. Number of points in x,y and z directions are
determined according to the formula:

::

    factor = ( SPACE / (a*b*c) ) ** (1/3)
    x_points = factor*a
    y_points = factor*b
    z_points = factor*c

where a, b and c are dimensions of the smallest rectangular box
enclosing the molecule.


The SURFace command
^^^^^^^^^^^^^^^^^^^

The :chm:`COOR SURFace` command computes the Lee and Richards surface for
selected atoms and stores the result in the appropriate weighting
array. If the :chm:`WEIGhting` keyword is used, the radii are obtained from
the weighting array (and then written over), otherwise the radii are
obtained from the parameter file values. The radius of the probe may
be specified (default 1.6) and the accuracy may be specified (default 0.05).
Either :chm:`ACCEssible` surface (default) or :chm:`CONTact` surface may be specified.
Contact surface is equivalent to Accessible surface if a zero probe
radius is used.  If the accuracy is not specified (or set to zero), then
the analytic result is provided.  If a nonzero accuracy is provided,
then the original Lee and Richard's (points on a sphere) algorithm
is used.


The HELIX command
^^^^^^^^^^^^^^^^^

The :chm:`COOR HELIx` command will analyze a single helix, or the relative
orientation of two helices.  The use this command, one or two atom
selections should be provided selecting ONLY the atoms which will be
used to define the helix.  The order of these atoms is important.
With a single atom selection, this command calculates the normalized
axis (A) and the perpendicular vector (R0) from the origin to A of
the cylinder most closely approximating a helix on which the selected
atoms best fit (Algorithm by J. Aqvist Computers & Chemistry
Vol. 10, pp97-99, (1986)).

With a double atom selection, this command also computes helix
axis and helix-helix structure analysis (Algorithm by Chotia, Levitt, and
Richardson JMB 145, P215-250 (1981)).


The CONVert command
^^^^^^^^^^^^^^^^^^^

The :chm:`COOR CONVert` command will cause the coordinates of all
defined and selected atoms to be transformed from the unit cell to
cartesian coordinates or back from cartesian to fractional coordinates.

Two orientations in cartesian coordinates are supported :

 ================ ==============================================
 :chm:`ALIGned`   in which b-vector is along y-axis and a-vector
                  in xy-plane (this is old charmm standard)
 :chm:`SYMMetric` in which shape matrix constructed from unit
                  cell vectors is symmetric
 ================ ==============================================

Two keywords in any order :chm:`[FRAC|ALIG|SYMM]` are required after :chm:`CONVert`.
Unit cell parameters (a,b,c,alpha,beta,gamma) follow in the same line.

The angle values are specified in degrees. See the routine CONCOR for
details concerning the transformation.

As an example, the following manipulations should have no net affect on the
coordinates,

::

      COOR COPY COMP
      COOR CONVERT SYMMETRIC  FRACTIONAL 5.6 12.2 5.4 80.0 95. 100.
      COOR CONVERT FRACTIONAL SYMMETRIC  5.6 12.2 5.4 80.0 95. 100.
      COOR CONVERT SYMMETRIC  ALIGNED    5.6 12.2 5.4 80.0 95. 100.
      COOR CONVERT ALIGNED    FRACTIONAL 5.6 12.2 5.4 80.0 95. 100.
      COOR CONVERT FRACTIONAL ALIGNED    5.6 12.2 5.4 80.0 95. 100.
      COOR CONVERT ALIGNED    SYMMETRIC  5.6 12.2 5.4 80.0 95. 100.
      COOR DIFF
      COOR STAT

When working with a triclinic system, the user should be aware of the form
of the coordinates.  Most of the data from crystallography is in fractional
(coordinates between zero and one) or in the aligned frame.

.. note::
   All of the internal use in CHARMM for energy calls, minimization,
   or dynamics ASSUMES that the coordinates are in the symmetric frame.


The COVAriance command
^^^^^^^^^^^^^^^^^^^^^^

The covariance command under coordinate manipulations
computes covariances of the spatial atom displacements of
a dynamics trajectory for selected pairs of atoms.

.. math::

   \mu_{JK} &= E( (R_J - E(R_J)) (R_K - E(R_K)) ) \\
            &= E( R_J R_K ) - E( R_J ) E( R_K )

and the normalized covariance matrix is given by

.. math::

   C_{JK} = \mu_{JK} / \sqrt{ \mu_{JJ} \mu_{KK} }

The command syntax and variables are as in the :chm:`coor dynamics` command.
The exceptions are the keywords:

   ====================== ==============================================================
   :chm:`SET1`            specifies the selection for the "J" groups in covariance
   :chm:`SET2`            specifies the selection for the "K" groups in covariance
   :chm:`UNIT_for_output` specifies unit for output of covarience matrix (ascii)
   :chm:`RESIdue_average` is a logical for computing the average over
                          residues in SET2 specification.  When followed by
   :chm:`NSETS`           equal to 2 the average is over both SET1 and SET2
                          giving a NRES1 x NRES2 covariance matrix.
   :chm:`MATRix`          gives output of just the covariance values in a matrix format
   :chm:`ENTRopy`         config. entropy [kcal/mol/K] using approximation S'' of
                          Andricioaei&Karplus (J. Chem. Phys 115,6289 (2001)) or
   :chm:`SCHL`            J. Schlitter's variation S'
                          (Chem. Phys. Lett. 215, 617 (1993)) on Karplus&Kushick.
                          See also Schafer et al  J. Chem. Phys. 113, 7809 (2000).
                          This approximation is an upper limit to the true entropy.
                          Sets CHARMM variable ENTROPY
                          It is recommended to remove translational(rotational) motion
                          before extracting the entropy (merge orient..[norot].);
                          for flexible molecules removal of rotation may be tricky.
                          NB! The covariance matrix used for this calculation is
                          not normalized and is 3N by 3N
   :chm:`TEMP`            temperature used in entropy calculation (default 298.15)
   :chm:`DIAG`            use only diagonal elements of covariance matrix,
                          mainly for testing purposes
   :chm:`RESI`            evaluate entropy using covariance for each residue only
   ====================== ==============================================================

Example:

::

   !Get configurational entropy at T=300K and save the unnormalized covariance
   !matrix, using all atoms in the PSF
   coor cova firstu 51 nunit 1 entropy matrix unit 61 temp 300.0
   ! Same without saving or printing the matrix and with output for each residue
   coor cova firstu 51 nunit 1 entropy unit -1 temp 300.0 resi


The DMAT command
^^^^^^^^^^^^^^^^

This command is accessed with the command :chm:`COOR DMAT` and provides some
general tools for the calculation, manipulation and storage/extraction of
distance matrix based properties.  This routine has some overlap with the
new distance command introduced by Bernie Brooks but also provides significant
complementarity in extending the range of properties computed.
The entire syntax is:

::

    COORdinates DMAT -
        RESIdue_average NOE_weighting -
        SINGle -
        FIRSt_unit <int> NUNIt <int> BEGIn <int> SKIP <int> -
        STOP <int> 2x<atom selection (SET1, SET2)> -
        UNIT_for_output <int>  TRAJectory CUTOff <real> -
        PROJect UPRJ <int> [MKPRoj] PROBability UPRB <int> TOLE <real> -
        [ [RELAtive] RMSF] [DUNIt <int>] [MATRix]

The command structure is like that of most other coordinate manipulation
commands other sub-parser keywords are:

    ======== ===================================================================
    UNIT     the distance matrix will be written to the unit
             number specified as an ASCII file unless the TRAJ
             keyword is specified, in which case a binary "trajectory" of
             the distance matrix will be written.
    RESIdue  this keyword specifies to compute the distance matrix
             for a center of geometry weighted average of residues
    NOE      this keyword denotes that the averaging over distances
             in the distance matrix should be inverse sixth power
             weighted.
    TRAJ     write a dynamic trajectory file of the distance matrix
    SINGle   process only a single coordinate file
    CUTOff   print only those values of the distance matrix which are
             smaller than cutoff value
    PROJect  project out a subset of contacts for printing
    UPRJ     read projection matrix from unit UPRJ
    MKPRoj   A projection matrix will be printed. Its elements are 1 if
             the distance is < CUTOff, 0 otherwise. To be used with subsequent
             PROJ UPRJ unit command. (If a standard DMAT is used as projection
             matrix the CUTOff in the PROJ command has to be squared)
    PROB     compute the contact probability based on differences
             from reference contact map read from UPRB and with
             an upperbound tolerance of TOLE
    RMSF     Computes the root mean square fluctuation in the distance
             matrix from the trajectory. Disables the printing of
             the binary file.
    RELAtive Divides the RMSF value by the distance
    DUNIt    Write distances to file open on the specified unit. This
             allows calculation of distance and (relative) fluctuation
             matrices in one pass.
    MATRix   Output is in the form of a rectangular matrix with just the
             z-values (distances or fluctuations)
    ======== ===================================================================

.. note::
   The binary file produced is analogous to the binary trajectory files and
   contain the following information:

   ::

                  WRITE(UNIT) HDRD,ICNTRL
                  CALL WRTITL(TITLEA,NTITLA,UNIT,-1)
                  WRITE(UNIT) NSET1,NSET2
                  WRITE(UNIT) (IND1(I1),I1=1,NSET1)
                  WRITE(UNIT) (IND2(I2),I2=1,NSET2)

   and then nframes of

   ::

                  WRITE(UNIT) ((CO(I1,I2),I1=1,NRES1),I2=1,NRES2)

   Where ICNTRL is a 20 element integer array with the following data:

   ::

                  ENDDO
                  ICNTRL(1) = (STOP - BEGIN)/SKIP
                  ICNTRL(2) = BEGIN
                  ICNTRL(3) = SKIP
                  ICNTRL(4) = STOP - BEGIN
                  ICNTRL(5) = NSAV
                  ICNTRL(8) = NDEGF
                  ICNTRL(9) = NATOM - NFREAT
                  CALL ASS4(ICNTRL(10),SKIP*DELTA)
                  IF(LNOE) THEN
                     ICNTRL(11) = 1
                  ELSE
                     ICNTRL(11) = 0
                  ENDIF
                  IF(LRESI) THEN
                     ICNTRL(12) = 1
                  ELSE
                     ICNTRL(12) = 0
                  ENDIF

   and NSET1[2] are the number of atoms comprising the two selections and
   IND1[2](NSET1[2]).  The distance matrix CO(NRES1,NRES2) is a 2-D array of
   size either NSET1 x NSET2 or NRES(NSET1) x NRES(NSET2) depending on
   whether the residue flag was used in processing the commands

Examples of usage:

1.  Compute the distance matrix for a single coordinate file (resident
in the main coordinate set) and print this matrix to a file linked to
fortran unit 1.

::

   open unit 1 write form name total.dmat
   COOR DMAT SINGLE UNIT 1 SELE ALL END SELE ALL END

2.  Compute the side chain-side chain center of geometry distance map
from a single coordinate file and print the distance matrix to unit 1
zeroing all elements of the matrix with distances greater than 6.5
angstroms

::

   define bb select ( type ca .or. type n .or. type c .or. typ o ) end
   define side select ( (.not. bb) .and. (.not. hydrogen) ) end

   open unit 1 write form name side.dmat

   coor dmat residue_average single unit 1 cutoff 6.5 select side end -
        select side end

3.  Compute the average hydrogen atom-hydrogen atom distance map from
a trajectory file on unit 10 and print the average distance matrix to
unit 1.  Use NOE inverse-sixth power weighting in the averaging and
"filter-out" all distances in the final map with values greater than
6.0 angstroms.

::

   open unit 10 read unform name trajectory.crd
   open unit 1 write form name noe.dmat

   coor dmat unit 1 cutoff 6.0 noe_weighting select hydrogen end -
        select hydrogen end -
        first_unit 10 nunit 1 begin 100 skip 100 stop 10000

4.  Compute the center-of-gemoetry distance matrix for side chains and
write this as a binary "trajectory" file to unit 1.  Read the
trajectory from unit 10.

::

   open unit 10 read unform name trajectory.crd
   open unit 1 write unform name side.dm-trj

   define bb select ( type ca .or. type n .or. type c .or. typ o ) end
   define side select ( (.not. bb) .and. (.not. hydrogen) ) end

   coor dmat residue_average unit 1 traj select side end select side end -
        first_unit 10 nunit 1 begin 100 skip 100 stop 10000

5.  Compute the center-of-geometry contact map probability based on a
precomputed distance matrix (e.g. from a PDB structure) based on a 6.5 A
cutoff. (This example is for the interdomain (helix-helix) contacts in
GCN4.  The two helices are segids zipa and zipb.)

::

   ! First contacts
   open unit 1 read unform name "traj/crdp/2zta/2zta_d1-60p.crd"
                          ! trajectory file to use to compute probability from
   open unit 2 write form name "distance_matrix/2zta_d1-60p.dmatp"
                          ! file to write contact probability matrix to
   open unit 3 read form name "distance_matrix/2zta_full.dmat
                          ! reference contact map

   coordinates dmat residue unit 2 -
           first 1 nunit 1 begin 100 skip 100 stop 600000 -
   	select side .and. ( segid zipa ) end -
           select side .and. ( segid zipb ) end -
           probability uprb 3 tole 0.3 cutoff 6.5

   close unit 1
   close unit 2
   close unit 3

6.  The following example shows the use of the dmat command to count the
number of contacts (native and non-native) throughout the course of a
trajectory using the distance matrix projection operator and the fact
that the number of contacts are accessible through the ?ncontact variable.

::

   label dotraj

   !  Now we loop over the trajectory and compute time dependent properties
   open unit 1 read unform name "traj/crdp/2zta/2zta_d1-60p.crd"
   open unit 10 write form name "distance_matrix/2zta_d1-60p.traj"
   write title unit 10
   *# Properties for Contacts
   *# trajectory 2zta_d1-60p.
   *# time(ps)   C(native)    C(total)
   *

   traj iread 1 nread 1 begin 500 skip 500 stop 600000
   set time 1.0
   set frame 1
   label loop

   trajectory read

   !  First get the contact information
   open unit 3 read form name "distance_matrix/2zta_full.dmatp"
                        ! reference distance matrix to use for projection
   open unit 2 write form name "distance_matrix/temp.dmat"
                        ! junk distance matrix
   coor dmat single residue unit 2 cutoff 6.5 -
        select ( side .and. segid zipa ) end  -
        select ( side .and. segid zipb ) end  -
        proj uprj 3

   set cnat ?ncontact

   open unit 2 write form name "distance_matrix/temp.dmat"
   coor dmat single residue unit 2 cutoff 6.5 -
        select ( side .and. segid zipa ) end  -
        select ( side .and. segid zipb ) end

   set ctot ?ncontact

   !  Write information to file
   write title unit 10
   * @time   @cnat    @ctot

   incr time by 1.0
   incr frame by 1
   if frame lt 1200 goto loop


The ANALysis command
^^^^^^^^^^^^^^^^^^^^

Analysis module for computing solvent averaged properties
It is accessed from the coordinate manipulation
part (CORMAN) of CHARMM and is used with the following syntax.  This
piece of documentation is still under development.  CLBIII 1/1/1990

.. note::
   Keyword syntax changed after c25a2!!
   Unit numbers for output to file have to be specified, and
   the trajectory is now specified in the usual way with BEGIN,SKIP,STOP
   LNI 11/11/96

Keywords:

  ================= ===============================================================================
  SOLVent           specifies analysis is to be of pure solvent, which means xref, yref
                    and zref, or site keywords are inappropriate, i.e., analysis all configurations
                    of solvent using all solvent molecules. OBSOLETE)
  WATEr             specifies the solvent is water (acutally any three-site molecule),
                    and forces all distinct g(r)'s to be computed, i.e., g_oo, g_oh and g_hh.
                    The first atom selection specifies the solvent atoms/molecules to be analyzed.
  SPECies           specifies the solvent species.  If SOLVent is active then all
                    solvent molecules to be analyzed should be specified here, e.g., all of them
                    present in the simulations.  This keyword is followed by the standard selection
                    syntax and is terminated with the FINIsh_solvent_specification keyword.
                    OBSOLETE)
  SITE              Specifies the collection of atoms around which you would like to compute
                    solvent properties, e.g., if you would like to analyze the solvent distribution
                    and velocity correlation function around the center of geometry of a trp
                    residue this keyword would be followed by the selection syntax which selects
                    that residue.
  XREF, YREF, ZREF  specifies that solvent analysis around a specific spatial
                    position, (xref, yref, zref) is to be carried out.  This is the same as the
                    site keyword, as far as the analysis of solvent configurations it invokes,
                    however, this site is static whereas the SITE keyword permits selection of a
                    dynamically evolving site. The above dimensions ar taken from trajectory stored
                    information for crystal runs (w/ charmm22 or later)
  CROSs             allows the selection of two subset of atoms for g(r) analysis
                    (a&b: 'a' are the atoms specified by the first selection and 'b' are the atoms
                    specified by the second selection).  The g(r) for a-vs-b and b-vs-b are
                    calculated and returned in units IOH and IHH respectively.
                    g(r) for a-vs-a will be returned in unit IGDIst.

                    Note that CROSs does not exclude form the analysis the couple of atoms
                    belonging to the same segid since it is design for the analysis of
                    independent subset of solvent molecules.
  ================= ===============================================================================

.. note::
   The keyword CROSs cannot be selected with the following options:
   WATer, SITE, IKIRkg, ISDIst, IFDBf.
   IVAC, IMSD, IFMIn were not tested with CROSs.
   IVAC cannot be combined with any analysis requiring coordinates
   IGDIST and ISDIST are mutually exclusive flags

     ========= ==================================================================
     NCORs     number of steps to compute vac or msd
     RSPIn     inner radius for vac,msd, analysis around REF (or SITE)
     RSPOu     outer radius for vac,msd, analysis around REF (or SITE)
     RDSP      radius of dynamics sphere, used for densities, kirkwood and dbf
     DENS      density (atoms/A**3) to use in normalization of g(r) if the value
               as calculated from the density within RDSP is not satisfactory
     DR        grid spacing for analysis of rdf's
     RSPHere   radius around REF to use for rdf analysis
     MGN       number of points in g(r) curve
     RCUT      radius of interaction sphere in dbf calculation
     ZP0       initial reference site - dynamics sphere origin separation
     NZP       number of separations to compute dbf
     TYP       for DBF calc 1=oxygen, 1=hydrogen

     IHIS      unit for output of 3Dhistogram data (in "DN6" format) or
     IPDB      unit for output of "atoms" where density exceeds THREshold
     ========= ==================================================================

   with options:

     ============ =============================================================
     WEIG         use WMAIN to weight points       !! Not tested
     DIPO         accumulate dipole vector density !! NOT working yet (June 98)
     CHARge       accumulate charge density        !! Not tested
                  default is to just accumulated number density of sel. atoms
     NORM value   densities are divided by this value (and by number of frames)
                  (default 1)
     XMIN,XMAX,DX
     YMIN,YMAX,DY grid dimension&spacing (default +/- 20A,0.5A spacing)
     ZMIN,ZMAX,DZ
     THREshold    value for density to output atoms in PDB file format
     ============ =============================================================

The atoms indicated by the solvent selection are analyzed. If dipole
data is to be analyzed the selection should contain 1 atom/group - the
groups define what atoms are to be used for the dipole calculation.
This could be automated; also need minimum image combined with orienting
function.

  ==== =======================================================================
  IDIP specifies a unit to which a simple dipole distribution will be plotted.
       This facility is intended for use with polarisable modelling of bulk
       solvent, and requires the FLUCQ compilation keyword for activation.
       (If IDIP is not specified, then no distribution is plotted.)

       ========= ====== ====================================================
       MINDipole  real  The minimum dipole (in Debye) to plot (default 0)
       MAXDipole  real  The maximum dipole to plot (default 4.0 Debye)
       NUMDipole  int   The number of sampling points to use (default 100)
       ========= ====== ====================================================
  EXVC EXcludedVolumeCorrection for use with ISDIST - the soulte-solvent g(r)
       is corrected for the volume excluded around the solute (ie the SITE)
       by the atoms in the selection following EXCV. This correction is
       computed using a Monte Carlo procedure with parameters:

       ========= ====== ====================================================
        MCP       int   Total number of points to use in the Monte Carlo
                        (default 1000)
        MCSHells  int   Total number of equal volume shells to spread
                        the MCP in (10)
        RPRObe    real  Probe radius (1.5A); a point is considered as excluded
                        if it is within RPRObe+VDWR(i) of any atom i in
                        the EXVC set
        ISEEd     int   Seed for random number generator (3141593)
        WEIG            Use WMAIN instead of the vdW radii
       ========= ====== ====================================================
  ==== =======================================================================

The following has been found to give good results even when looking
at g(r) for water hydrogens around a site:

::

   scalar wmain = radius
   scalar wmain mult 0.85
   coor anal ...... EXVC select segid pept end -
         MCPoints 20000 MCSHells 20 WEIG RPRObe 0.0

The key is to make sure that the a non-zero accessible volume is obtained
at the shortest distances where g(r) starts being non-zero.
The data file produced with EXCV contains two extra columns; column 4 contains
the uncorrected g(r) and column 5 contains the accessible volume fraction.

EXAMPLES: (See also the test/c27test/solanal2.inp testcase)
The following examples use a trajectory of a short peptide in a periodic
water box

::

   ! MeanSquareDisplacement of all watermolecules to estimate diffusion coeff
   open unit 21 read unform name @9pept500.cor
   open unit 31 write form name @9pept500.msd
   coor anal select type oh2 end  -     ! what atoms to look at
         firstu 21 nunit 1 skip 10 -    ! trajectory specification
         imsd 31 -                      ! flag to do the MSD analysis
         rspin 0.0 rspout 999.9 -       ! we are interested in ALL waters
         ncors 20 -                     ! compute MSD to NCORS*SKIP (0.04ps)steps
         xbox @6 ybox @7 zbox @8        ! and we did use PBC

   ! g(r) for the waters; the program defaults are used to calculate the density
   ! using selected atoms within 10A (RDSP keyword) of the reference point (0,0,0)
   ! (REF keyword)
   open unit 21 read unform name @9pept500.cor
   open unit 31 write form name @9pept500.goo
   open unit 32 write form name @9pept500.goh
   open unit 33 write form name @9pept500.ghh
   ! specify WATEr to get all three g(r) functions computed
   coor anal water select type OH2 end -
         firstu 21 nunit 1 skip 10 -    ! trajectory specification
         igdist 31 ioh 32 ihh 33 -      ! flag to do the solvent-solvent g(r)
         mgn 100 dr 0.1 -               ! comp. g(r) at MGN points separated by DR
         rsph 999.9  -                  ! use ALL waters for rdf calculation
         xbox @6 ybox @7 zbox @8        ! and we did use PBC

   ! g(r) backbone amide hydrogen -  water oxygens
   ! if a single solute atom is looked at the MULTi keyword is not necessary
   ! when several solute atoms are specified as the site, their average position
   ! will be used as the reference position if MULTi is not present
   open unit 21 read unform name @9pept500.cor
   open unit 31 write form name @9pept500.gonh
   coor anal select type oh2 end  -     ! Water oxygens
         site select type H end multi - ! and the amide hydrogens
         firstu 21 nunit 1 skip 10 -    ! trajectory specification
         isdist 31  -                   ! do the g(r) (here solute-solvent)
         mgn 100 dr 0.1 -               ! comp. g(r) at MGN points separated by DR
         rsph 999.9  -                  ! we use ALL waters for the calculation
         xbox @6 ybox @7 zbox @8        ! and we did use PBC

   ! g(r) for GLY3 NH - the water oxygens - with excluded volume correction
   open unit 21 read unform name @9pept500.cor
   open unit 31 write form name @9pept500.gn3ox1
   coor anal  select type OH2 end -
         site multi select atom pept 3 H end -
         EXVC select segid pept end -
         MCPoints 2000 MCSHells 20 RPRObe 1.7 -
         firstu 21 nunit 1 skip 50 -    ! trajectory specification
         isdist 31 -                    ! flag to do the solvent-solvent g(r)
         mgn 100 dr 0.1 -               ! comp. g(r) at MGN points separated by DR
         rsph 999.9  -                  ! we use ALL waters for the calculation
         xbox @6 ybox @7 zbox @8        ! and we did use PBC


Subcommand RCOR (Rotational Correlation Time of Water)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Calculation of rotational correlation times corresponding to the three
rotational motions of a water molecule has been added to the solvent
analysis code. The three rotational motions refer to motion around the
dipole axis (twist), around an axis perpendicular to the molecular
plane (rock) and around an axis parallel to the H-H vector (wag) (Ref 1).
The correlation time is calculated by fitting the exponentional decay part
of the corresponding time correlation function C(t) to an
exponentional function of the form C(t) = A exp(-t/tau) where tau is
the correlation time. The direct correlation functions were calculated
via FFT method using the CORFUNC subroutine in the CORREL.SRC. The
calculation can be invoked by assigning a non-zero integeer value to
the keyword RCOR.


Keywords for rotational correlational time calculation are:

  ==== ========= =============================================================
  RCOR <integer> if RCOR > 0, invokes rotational correlational time analysis
  ROUT <unit>    write the three correlation functions of selected waters
                 into a fortran unit
  TLOW <real>    lower limit of time for fitting, default is 1.0ps
  TUP  <real>    upper limit of time for fitting, default is 4.0ps (Ref 2)
  MAXT <integer> maximum number of time steps, default is 512
  P1             compute P1 dipole correlation instead of wag/twist/rock
                 (< u(t)u(t+tau)>, where u is unit vector along water dipole
                 output is to unit specified by ROUT
  P2             compute P2 dipole correlation instead of wag/twist/rock
                 (<P2( u(t)u(t+tau) )>, where u is unit vector along
                 water dipole; P2(x)=(3x**2-1)/2
                 output is to unit specified by ROUT
  ==== ========= =============================================================

For P1 and P2 the analysis may be performed in a shell defined by RSPIn
and RSPOut, and the minimum image  xbox,ybox,zbox is also accounted for

REFERENCE:

1. Johannesson, H. and Halle, B. J. Am. Chem. Soc. 1998, 120, 6859-6870
2. Wallqvist, A. and Berne, B. J. J. Phys. Chem. 1993, 97, 13841-13851


EXAMPLE: see test/c27test/solanal2.inp

::

   ! Rotational Correlation Time of Water
   open unit 21 read unform name @9pept500.cor
   open unit 31 write form name @9pept500.rcor
   coor anal sele .byres. (type oh2  -  ! select all three atoms of water
     .and. (resn asp .and. type od1) -
     .around. 3.5) show end    -
     firstu 21 nunit 1 skip 10 -
     rcor 1                    -    ! rot corr time calculation
     timl 1.0 timu 3.0         -    ! lower and upper time limits for linear fit
     rout 31                   -    ! corr coef to unit 31
     xbox @6 ybox @7 zbox @8        ! and we did use PBC


Subcommand IHYD: Hydration Number Calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is to calculate hydration number or, in general, the number of solvent
molecules within a specified distance of a multi atom or single atom site:

* number of solvent molecules (residues) withn RHYD of the solute
* number of solvent atoms within RHYD of the solute
* number of solvent atoms within RHYD of solute atoms (ie, if three water
  molecules are all within RHYD of a 7-atom solute this will be 63)

Sets CHARMM variables NHYDRR, NHYDAR and NHYDAA to the averages for these
three numbers.
If IHYDN>0 these numbers are written to unit IHYD every timestep.
At the end averages over the trajectory are printed in the output file.

Hydration number calculation is invoked by specifying a non-zero cutoff RHYD.
NB! You need keyword MULTi if the solute (the SITE) has more than one atom.

Keywords for hydration number calculation are:

  ==== ========= ===================================================
  IHYD <integer> if IHYDN > 0, output to unit IHYDN each timestep
  RHYD <real>    calculate hydration number at this distance from
                 each atom in the site
  ==== ========= ===================================================

Example:

::

   ! Calculate hydration no
   coor anal sele resn tip3 .and. type oh2 end -
         site select resn asp .and. type od1 show end multi -
         firstu 21 nunit 1 skip 5 -
         rhyd 3.0 -                    ! calculate hyd no at 3.0A
         xbox @6 ybox @7 zbox @8


The DRAW command
^^^^^^^^^^^^^^^^

The DRAW command (called directly from CORMAN, not to be
confused with the DRAW command found under the ANALysis command)
is useful for displaying molecules. The output is a command
file that can be read by various displaying and plotting programs.
This command file can be edited for different types of displaying.
In addition to atom positions and bonds, velocity and forces may
also be displayed. The current keywords are:

   ====== ====================================================
   NOMO   No molecule option (only velocities or derivatives)
   DFACt  Derivative factor                (default 0.0)
   DASH   Spacing of dashed line used for Hbonds (default .01)
   FRAMe  Specifies that a frame tag will be written first
          (default - dont specify frame)
   RETUrn Specifies which stream the plotting program will
          return to after plotting this section (default none)
   ====== ====================================================

An atom selection is also looked for. Any atom not selected will
not be considered. The default is to include all atoms.


The HBONd / The CONTact command
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The HBONd command analyses a trajectory, or the current coordinates,
for hydrogen bonding patterns.

The form COOR CONTact ... ignores the hydrogen bond donor/acceptor
definitions in the psf and looks for all contacts which satisfy the
distance cutoff criterion between all atoms in the two selections; possibly
bridged by a residue as defined by the BRIDge keyword. This is useful for
hydrophobic contact analysis, or for salt bridges. No angle cutoff can
be used with this form of the command.
Output and other options are as for the COOR HBONd variant.

The form COOR HBONd makes use of the DONOR/ACCEPTOR definitions in the psf.
For each acceptor/donor in the first selection the average number and average
lifetime (for trajectories only) of hydrogen bonds to any atom in the second
selection is calculated. A hydrogen bond is assumed to exist when two
candidate atoms are closer than the value specified by CUT (default 2.4A,
(reasonable criterion, DeLoof et al (1992) JACS 114,4028), and if a value
for CUTAngle is given the angle formed by D-H..A is greater than this CUTAngle
(in degrees, 180 is a linear H-bond); the default is to allow all angles.
The current implementation assumes that hbonding hydrogens are present in
the PSF and uses ACCEptor and DONOr information from the PSF to determine
what pairs are possible. If output is wanted to a separate file the IUNIt
option can be used. If the BRIDge option is used the routine calculates average
number and lifetime of bridges formed between all pairs of atoms in the
two selections; a bridge is counted when a residue of the type specified with
the BRIDge <resnam>  hydrogen bonds (using same criteria as for direct
hbonding) to at least one atom in each selection. The typical
use of this would be to find water bridges. Here again, results are presented
for each atom in the first selection.

If FIRSTunit is not specified the current (MAIN) coordinates are analyzed.

Periodic boundary conditions are taken into account using the hardwired
minimum image code (:ref:`images_mipb`) if keyword PBC is
given. Supported geometries are:

===================== ========= ========================== ======================
Geometry              Keyword   Required information       Auxiliary information
===================== ========= ========================== ======================
"Orthogonal"          CUBIC     BOXL (or XSIZE)            YSIZE, ZSIZE if
                                                           different from XSIZE
Truncated octahedron  TO        BOXL (crystal A parameter)
Rhombic dodecahedron  RHDO      BOXL (crystal A parameter)
===================== ========= ========================== ======================

If crystal information is present in the trajectory it will be used to
set the actual box dimensions (overriding the value(s) specified on the
COOR command line). The minimum image code is turned off when the command
exits, which means that a previous BOUND command will no longer be in effect.

Keyword :chm:`VERBose` provides a more detailed output:

For trajectory analysis the duration and endtime (ps) of each H-bond,
or bridge, together with a specification of the atoms involved is output;
potentially very large amounts of data! Only hbonds/bridges with a lifetime
longer than the value specified by keyword TCUT (default 0.0 ps) are included
here and in the summary.

.. note::

   TCUT (and NSKIP) may influence the results, since hbonds with
   a duration < TCUT are not counted, and for the lifetime analysis a quick
   fluctuation in hbond distance may with one choice of NSKIP result in the
   hbond being perceived as broken at that instant, whereas with a longer NSKIP
   the event would not have been noticed, resulting in a longer lifetime
   being reported.

For single coordinate set analysis the VERBose keyword results in a more
detailed listing giving all atoms involved, and also the geometry for
direct hbonds.

For each donor/acceptor in the first selection the trajectory analysis outputs
the AVERAGE NO. of hydrogens bonds this atom has had during the trajectory
(aveno=sum over frames(number of hbonds formed by this atom)/(number of frames)
the average lifetime is defined as
avelife=
sum over hbonding events(duration of hbond between two atoms)/(number of
different hbonds formed by these atoms)
(ie, hbonds that have been broken for at least one frame between events)
Note that the lifetime can be influenced by end-effects (ie hbonds
still active at end of trajctory are counted as being terminated then!)

Output can be directed to a separate file specified by IUNIT int.

The following charmm substitution parameters are set in the module:

  ============= ================================================================
  :sub:`NHBOND` total number of hydrogen bonds for selected atoms (timeaveraged)
  :sub:`AVNOHB` average number of hydrogen bonds over selected atoms (timeaver.)
  :sub:`AVHBLF` average lifetime of hydrogen bonds
  ============= ================================================================

Note that these averages are over the selected atoms, which may include
a number of atoms with no hbonds > TCUT!

Distance and lifetime histograms can be computed for all (putative) hydrogen
bonds encountered in the analysis; ie, the distance histogram will in general
contain non-zero data also for bins > CUT. For bridges the lifetimes are those
of the bridging events, but the distances are computed from all individual
hydrogen bonds.

The three columns in the output are:

::

   distance (or time)   counts     counts/NSTEP

where NSTEP is the number of frames that have been analyzed from the
trajectory.

   ======= ========== ==========================================================
   Keyword  default    meaning
   ======= ========== ==========================================================
   IRHI      -1        unit to which distance histogram will be written
   DRH      0.05       bin size for distance histogram (A)
   RHMAx    10.0       distance in maximum bin (collects all distances >= RHMAx)
   ITHI      -1        unit to which lifetime histogram will be written
   DTH      5.0        bin size for lifetime histogram (ps)
   THMAx    1000.0     time in maximum bin (collects all times >= THMAx)
   ======= ========== ==========================================================


The HISTogram command
^^^^^^^^^^^^^^^^^^^^^

This command computes a histogram along the X,Y,Z or Radial directions
for the selected atoms.
The histogram can either be a simple count of the number of atoms
contained in each bin (specified by the HNUM=number of bins between
HMIN,HMAX keywords), or if the WEIGhting keyword is present the WMAIN
array is summed for the atoms in each bin.
HSAVe specifies that the histogram should be saved and incremented at
the next invocation of COOR HIST. HPRInt specifies that the resulting
histogram should be printed. For X,Y,Z histograms the output is
the accumulated density/HNORM (default=1.0) in each bin. If HDENS>0.0
(default=0.0) there is also a third column for R histograms containing
the accumulated density/(volume of shell containing this bin)/DENS.

The COMParison keyword results in XCOMP,YCOMP,ZCOMP,WCOMP being used.

The variable ?NCONFIG is set to the number of configurations (frames)
that have been accumulated so far.

The results may be output to a file specified by IUNIt int.

EXAMPLE:
To average the charge density in spherical shells from a trajectory
could be done in the following way:

::

   scalar wmain=charge

   traj iread ....

   set i 1
   label loop
   traj read
   !if you are reading velocities, you may want to convert to A/ps
   ! (and then you wouldn't use the weighting option like this)
   ! scalar x divi ?TIMFAC
   ! scalar y divi ?TIMFAC
   ! scalar z divi ?TIMFAC
   coor hist R hnum 50 hmin 0.0 hmax 10.0 hsave weig
   incre i by 1
   if i .lt. 100 goto loop

   ! you could also normalize for number of selected atoms
   ! set scale ?NSEL
   ! mult scale by ?NCONFIG
   ! then use @scale instead of ?NCONFIG below
   bomblevel -1 ! to get by the zero atom selected warning below
   coor hist R hnum 50 hmin 0.0 hmax 10.0 select none end hprint -
    hnorm ?NCONFIG [ hdens 0.03 (some reasonable bulk density/A**3) ]


The PUCKer command
^^^^^^^^^^^^^^^^^^

::

   COORdinates PUCKer [SEGId segid] RESId resid1 [TO resid2] [AS | CP]

The sugar pucker phase and amplitude, as defined by
Altona&Sundaralingam (default, keyword AS)  or (CP) Cremer&Pople (JACS 1975),
are calculated for the (deoxy)ribose of the specified residue(s);
the first segment is the default. A range of residues from resid1 TO resid2
can be analyzed.


The INERtia command
^^^^^^^^^^^^^^^^^^^

::

   COORdinates INERtia [atom-selection]

Principal moments of inertia I_xx, I_yy, I_zz are calculated and
the eigenvectors of the inertia tensor are printed. Normally atom selection
should not be used and the command

example:

::

   COOR INER

is sufficient, since all ithe atoms are selected by default. The units for
principal moments of inertia are

:math:`amu \cdot A^2`,  where amu - atomic mass unit (Carbon is 12), and A stands
for Angstrom.


The INERtia ENTRopy command
^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

   COORdinates INERtia [atom-selection] ENTRopy
                    [TEMPerature <real>] [SIGMa <real>] -
                    [STANdard <SOLUtion|GAS>]

Entropy calculation is an extension to the INERtia command.
In addition to calculation of principal moments of inertia the rotational
and translational entropy components will be evaluated. Calculation of
these two entropy terms is very fast. See :doc:`vibran.doc <vibran>` to see how to
calculate the vibrational entropy term.

Default value for TEMPerature is 298.15 K. Default SIGMa value is 1.0.
SIGMa is symmetry number which is 1 for non-symmetric molecule and some
low symmetry groups. For symmetric molecules one should enter a correct
value for sigma (see, for example, C.J.Cramer, "Essentials of Comp.Chem.",
2002,p.327).

Translational component of entropy depends on the defition of standard state.
There are two definitions: solution (1M) and ideal gas. The default is solution.
They differ by a constant of 6.35236 kcal/mol, with higher entropy in gas state.
See details inTidor and Karplus, J Mol Biol (1994) vol. 238 (3) pp. 405-14

example:

::

  COOR INER ENTRopy
  COOR INER ENTRopy TEMPerature 298.15 SIGMa 1
  COOR INER ENTRopy TEMPerature 298.15 SIGMa 1 STANdard SOLUtion
  COOR INER ENTRopy TEMPerature 298.15 SIGMa 1 STANdard GAS

  VIBRan
  DIAGonalize ENTRopy TEMP 298.15 SIGM 1
  DIAGonalize ENTRopy TEMP 298.15 SIGM 1 STANdard SOLUtion
  DIAGonalize ENTRopy TEMP 298.15 SIGM 1 STANdard GAS
  END

testcase in c32test/entropy.inp

The units for entropy are :math:`cal/(mol \cdot K)`. Rotational, translational, vibrational, and
total entropies can be accessed in CHARMM input file as ?SROT, ?STRA ?SVIB, and ?SSUM
substitution parameters.


The SECondaryStructure command (SECS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Computes secondary structure of residues in first-selection in the context of
the second-selection; eg, a beta-strand in the first-selection will be
rcognized as such if it forms appropriate hydrogen bonds to residues in the
second-selection. If no second-selection is given it is the same as the first
(which defaults to all). A residue is included if any atom in it is selected,
and amino acids are recognized by the presence of atoms named N,C and CA. The
amide hydrogen can be named either H or HN. Only operates on main coordinates.

Currently using Kabsch&Sander (Biopolymers 22, 1983, 2577) definition of
alpha-helix and beta-strand.

Sets CHARMM variables ?NALPHA and ?NBETA to number of residues in alpha/beta
structures, and ?ALPHA and ?BETA are set to fraction of residues with that type
of structure. The fraction is computed from number of peptide residues in the
first selection. On return Calphas have WMAIN-array set to 0, 1 (alpha), 2
(beta)

The default H-bond criterion is CUTH=2.6, slightly longer than the default
2.4A used in coor hbond (from DeLoof et al JACS 1992); this is to be slightly
more generous in defining secondary structures. CUTA can be used to define an
angle cutoff for the N-H..O angle (default is not to use this criterion).

Keywords QUIEt/VERBose control the amount of output


The CONFormational command
^^^^^^^^^^^^^^^^^^^^^^^^^^

::

  COORdinate CONFormational { <resname> } [ PRINT ] [ READ io-speficication ] -
                   [atom-selection] [COMP]

Current methods for generating transition paths between macromolecules e.g.,
the TMD and TREK modules, rely on the Cartesian coordinates of a subset of
atoms in a protein.  Although several residue types possess symmetry (e.g.
planar symmetry of a PHE ring), so that the conformation of such a residue is
invariant with respect to a rotation around the symmetry axis, rendering
certain groups of atoms effectively indistinguishable, topology files must
distinguish between these atoms (e.g. PHE CD1 vs. PHE CD2). Given two different
coordinate sets for a macromolecule, any two-set path generation method that
makes use of the Cartesian coordinates of atoms that belong to residues with
symmetry decides arbitrarily the correspondence between the `indistinguishable'
atoms. For example, performing TMD using coordinates of the ring atoms of a
PHE, will force the position of atom CD1 in the initial set to move to the
position of atom CD1 in the target set, although the movement from CD1 to CD2
is also possible. In such transitions, it is likely that there exist a path
with a high energy barrier (e.g. flipping of a PHE ring in a tightly-packed
protein interior) that can be avoided by making use of symmetry. The current
method, CONFormational consistency, is an algorithm for renaming certain atoms
to minimize rotation and flipping of the involved residues during path
generation.

The algorithm is heuristic and is as follows.  (Two coordinate sets are assumed
present, in the main and comparison sets).  For each residue in the optional
atom selection, the following procedure is performed.  The residue is
partitioned into three (non-disjoint) sets of atoms: swap atoms, orientation
atoms and test atoms.  Swap atoms are organized into pairs, which will be
swapped during the check. The residues in the two conformations are RMSD-
aligned based on the orientation atoms only. RMSD is computed between the test
atom positions in the two coordinate sets. The configuration of the swap atoms
that gives the lesser test-atom-RMSD value is accepted.  Positions of any
hydrogen atoms that are bonded to swap atoms are initialized, and can be
regenerated with HBUIld.

The three sets in the residue partitioning are defined by default for the
following residues (i.e. by default, {<resname>} can contain any number of
these)

::

  ARG ASP GLU HIS HSC HSD HSE HSP LEU PHE TYR VAL

Users can override pre-existing defaults for these residues, and declare new
residues in an optional input file.  In the following, the default residue
partitioning is shown for ARGinine (only the relevant atoms are shown):

::

                         HH11
                         |
            -- CD        NH1-HH12
                 \      //(+)
                  NE--CZ
                        \
                         NH2-HH22
                         |
                         HH21


  swap atoms:               NH1 NH2
  orientation atoms:        CZ NH1 NH2
  test atoms:               CD

Note that the HH* hydrogens will have undefined positions after the check is
complete, and can be redefined using HBUIld.  Also note that more than one
partitioning scheme may lead to the same results.

A custom residue partitioning file can be specified, following the READ
option.

For the twelve residue types supported by default, the equivalent partitioning
file is:

::

  12
  ARG 1 CD 1 CZ 1 NH1 NH2 0
  ASP 1 CA 2 CB CG 1 OD1 OD2 0
  GLU 1 CB 2 CG CD 1 OE1 OE2 0
  HIS 1 CA 2 CB CG 2 ND1 CD2 NE2 CE1 0
  HSC 1 CA 2 CB CG 2 ND1 CD2 NE2 CE1 0
  HSD 1 CA 2 CB CG 2 ND1 CD2 NE2 CE1 0
  HSE 1 CA 2 CB CG 2 ND1 CD2 NE2 CE1 0
  HSP 1 CA 2 CB CG 2 ND1 CD2 NE2 CE1 0
  LEU 1 CB 1 CG 1 CD1 CD2 0
  PHE 1 CA 3 CB CG CZ 2 CD1 CD2 CE2 CE1 0
  TYR 1 CA 3 CB CG CZ 2 CD1 CD2 CE2 CE1 0
  VAL 1 CA 1 CB 1 CG1 CG2 0

The first line specifies the number of lines to be read (number of residues)
Each subsequent line is organized as follows:

::

  <residue name> <# test atoms> <list of test atoms> -
                 <# orientation atoms that are not swapped> <list ...> -
                 <# PAIRS of orientation atoms that are swapped> <list...> -
                 <# swap atoms that are not part of the orientation set> <list...>

Note that the default residue partitioning file includes residues which do not
have any symmetry. These are histidine residues : HIS, HSD, HSE, HSP, and HSC.
In these cases the atoms ND1 and CD2 are assumed to be indistinguishable.

The optional PRINT command will print checking information for each tested
residue By default, the main comparison set is modified.  Specifying COMP will
cause the comparison set to be modified (note that this may lead to undefined
hydrogen atoms in the comparison set).

Finally, an atom selection may be specified. In this case, only the residues
for which at least one atom is selected will be tested.

Examples:

1)

   ::

     coor conf his arg phe tyr hsd glu asp print select all end

   will check the specified residues and, if needed, make modifications to
   the main set. Results for each residue will be printed. Default partitioning
   is used.

2)

   ::

     coor conf arg print select all end read
     * residue partitioning file
     *
     2
     ARG 1 CD 1 CZ 1 NH1 NH2 0
     ASP 1 CA 2 CB CG 1 OD1 OD2 0

   will check all arginines using the custom partitioning specified below the
   command line

   Testcase: c35test/confcons.inp

The PATH command
^^^^^^^^^^^^^^^^

::

  COORdinate PATH { NREP <int> } {NAME <character*>} [<PDB|FILE|UNFO|CARD|FORM>]

This command will create an interpolated path connecting two structures stored
in the main and comparison sets.  Currently, only linear interpolation in
Cartesian atom coordinates is implemented.

NREP specifies the number of replicas desired (this includes the two endpoints,
and must be at least three)

NAME specifies the base name of the file to which the interpolated coordinates
will be written.  An extension will be appended to the base name, which
consists of a number in the range [0.. NREP-1] followed by '.<ext>', in which
ext depends on the format specification as follows:

---------------- ---
format spec      ext
---------------- ---
PDB              PDB
FILE/UNFO/CARD   COR
---------------- ---

Example:

::

  coor path nrep 32 name output/conv card
  ! will create a linearly interpolated path of 32 replicas named
  ! output/conv0.cor, ..., output/conv31.cor
  ! in card format

Testcase: c35test/confcons.inp


.. _corman_substitution:

Coordinate Manipulation Values
------------------------------

There are several different variables that can be used in titles or
CHARMM commands that are set by some of the coordinate manipulation commands.
Here is a summary and description of each variable. See also :doc:`subst.doc <subst>` (which
may be more up-to-date).


* 'XAXI','YAXI','ZAXI','RAXI','XCEN','YCEN','ZCEN'

   A rotation axis vector and its length and the center of rotation.
   This data is set by the COOR AXIS, COOR LSQP, COOR ORIE, and COOR ORIE RMS
   commands.  These values may be used by any of the commands that uses the
   vector-spec with the AXIS keyword.

* 'XMIN','YMIN','ZMIN','WMIN','XMAX','YMAX','ZMAX','WMAX','XAVE','YAVE','ZAVE','WAVE'

   Statistics set by the COOR STAT command.

* 'THET'

   Angle of rotation set by the COOR ORIEnt command.

* 'XMOV','YMOV','ZMOV'

   Displacement of centers set by the COOR ORIEnt command.

* 'RMS'

   Resulting RMS value set by the COOR RMS, COOR ORIEnt, or COOR RGYR
   commands.

The TMSCore command
-------------------

Computes the TM-score between the selected sets of atoms.  The TM-score
(see Zhang, Y. and Skolnick, J. Proteins, 2004 57:702-710) is a scoring
function that quantifies the similarity between two structures, returning a
number between 0 and 1.  We assume that the sequences of the two structures
are identical.  The TM-score is computed as:

::

  TM-score = Max [ 1/N  sum_{i=1}^N  1/(1 + (di/d0)**2) ]

where di is the distance between the two structures of atom i, d0 is a
constant reference length that depends only on the number of residues in
the protein, N is the number of atoms selected, and the Max is computed
over many different alignment attempts of the two molecules (see Zhang and
Skolnick for more details).  The aim of the multiple alignments is to emphasize
the matching parts of the molecule.

After the command is executed, the TMScore, the TMScore with a cutoff of 10 A,
and the d0 value used to compute the TMScore are assigned to the variables
?tmscore, ?tm10 and ?tmd0, respectively.

Ex/

::

  coor tmsc sele type CA end
