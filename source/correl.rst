.. py:module::correl

=====================
Correlation Functions
=====================

The CORREL commands may be used to obtain a set of time series
for a given property from a trajectory. Once obtained, the time series
may be manipulated as required, saved or plotted, or to generate
correlation functions  ( :math:`C(\tau) = \langle A(t) \cdot A(t+\tau) \rangle` ). The correlation
functions may be manipulated, saved, plotted, and transformed to find
spectral density (Fourier transform of :math:`C(\tau)`), etc and determine the
correlation times.

Alternately, a covariance matrix may be computed for a collection
of time series. This option will compute the full matrix for use
in entropy calculations or for other applications.

Reorienting a coordinate trajectory is possible using the
COMPARE command. For details see :ref:`Merge <dynamc_reorient>`.

.. _correl_syntax:

Syntax for the CORREL command and subcommands
---------------------------------------------

::

   CORREL  [ MAXTimesteps int ]  [ MAXSeries int ]  [ MAXAtoms ] [ COVAriance]  -
                default 512         default 2        default 100

               [ nonbond-spec ] [ hbond-spec ] [ image-spec ] [NOUPdate]
               [  INBFrq 0    ] [  IHBFrq 0  ] [  IMGFrq 0  ]

   hbond-spec        *note Hbonds:(chmdoc/hbonds.doc).
   nonbond-spec      *note Nbonds:(chmdoc/nbonds.doc).
   image-spec        *note Images:(chmdoc/images.doc)Update.


   Subcommands:

   miscellaneous-commands

   COOR coordinate-manipulation-command

                   {   DUPLicate  time-series-name                         }
                   {                                                       }
   ENTEr  name     { [ BONDs   repeat(2x(atom-spec))        ] [ GEOMetry ] } c
                   { [ ANGLe   repeat(3x(atom-spec))        ] [ ENERgy   ] } c
                   { [ DIHEd   repeat(4x(atom-spec)) [NOTR] ]              } c
                   { [ IMPRo   repeat(4x(atom-spec))        ]              } c
                   {                                                       }
                   { [ ATOM ] [ X   ]  repeat(atom-spec) [ MASS ]          } e
                   { [ FLUC ] [ Y   ]                                      } c
                   {          [ Z   ]                                      }
                   {          [ R   ]                                      }
                   {          [ XYZ ]                                      }
                   {                                                       }
                   {   VECT   [ X   ]  repeat(2x(atom-spec))               } e
                   {          [ Y   ]                                      }
                   {          [ Z   ]                                      }
                   {          [ R   ]                                      }
                   {          [ XYZ ]                                      }
                   {                                                       }
                   {  ATOM DOTProduct repeat(2x(atom-spec)) [NORMal] [MASS]} e
                   {  FLUC DOTProduct repeat(2x(atom-spec)) [NORMal] [MASS]} e
                   {  VECT DOTProduct repeat(4x(atom-spec)) [NORMal]       } e
                   {                                                       }
                   {  ATOM CROSsproduct rep.(2x(atom-spec)) [NORMal] [MASS]} e
                   {  FLUC CROSsproduct rep.(2x(atom-spec)) [NORMal] [MASS]} e
                   {  VECT CROSsproduct rep.(4x(atom-spec)) [NORMal]       } e
                   {                                                       }
                   {  HBONd [4x(atom-spec)]*++       [ ENERgy  ]           } c
                   {                                 [ DISTance]           }
                   {                                 [ HANGle  ]           }
                   {                                 [ AANGle  ]           }
                   {                                                       }
                   {  DISTance repeat(2x(atom-spec))                       } c
                   {                                                       }
                   { [ GYRAtion ] [ CUT real ] [ MASS ]       s**          } c
                   { [ DENSity  ]                             s**          } c
                   {                                                       }
                   {   RMS  [ MASS ] [ ORIEnt ]               s**          } c
                   {   MODE mode-number                       s**          } c**
                   {   TEMPerature  [ NDEGF int ]             s**          } v
                   {   ENERGY                                              } c
                   {   HELIx [ SELE atom-selection END ]                   } c
                   {   INERtia  [ SELE atom-selection END ][ TRACe | ALL]  } c
                   {   PUCK  RESI <resid> [SEGI <segid>]                   } c
                   {                                                       }
                   {   PUCK ATOM [ Q ]      repeat(6x(atom-spec))          }
                   {             [ THETa ]                                 }
                   {             [ PHI ]                                   }
                   {                                                       }
                   {   USER user-value [ repeat(atom-spec) ]  s**          } e
                   {   SURF [ RPRPobe real] [ WEIG ] [ACCE|CONT] [RESI]    } c
                   {        [ SELE atom-selection END ]                    }
                   {                                                       }
                   {   CELL cell-spec                                      }
                   {   TIME [ AKMA ]                                       }
                   {   ZERO                                                }
                   {                                                       }
                   {   DIPO [ SELE atom-selection END ] [ OXYZ ] [ MASS ]  }
                   {   VECM repeat(2x(atom-spec))                          }
                   {                                                       }
                   {   SDIP [ SHELL int ] [ OXYZ ] [ MASS ]                } c***
                   {        [ BULK ]                                       }
                   {   SATM [ SHELL int ] atom-spec                        } c***
                   {        [ BULK ]                                       }
                   {   OMSC ETA friction-coefficient                       } c
                   {                                                       }
                   {   PRDI proto_num [ MASS ]                             } c****
                   {                                                       }
                   {   PRVE proto_num [ MASS ]                             } v*
                   {                                                       }

                       ( code: c-coordinates, v-velocities, e-either )
         c**  MODE time series is allowed only if CORREL is invoked from VIBRAN.

         s**  these utilize the first atom selection in the next TRAJ command.

         c*** needs a CHARMM executable with SHELL functionality
              see *note Shell:(chmdoc/shell.doc)

         c****,v*  needs a CHARMM executable with PROTo functionality
                   see *note Shell:(chmdoc/proto.doc)

         *++  Hydrogen bond atom order is one of:
                                  Donor,Hydrogen,Acceptor,Acceptor-antecedent
                                  Donor,Hydrogen,Acceptor
                                  Donor,Acceptor

         cell-spec::= one of { A B C ALPHa BETA GAMMa ALL SHAPe }

         atom-spec::= {residue-number atom-name}
                      { segid  resid atom-name }
                      { BYNUm  atom-number     }
                      { SELE atom-selection END} ***

         atom-selection::= see *note Selection:(chmdoc/select.doc)
         *** Note: If an atom-selection is used for atom-spec's, then
             all atom-spec's must be contained within one atom-selection

         *** WARNING: For angles and dihedrals, if SELE is used to
             specify atoms, then the order that the atoms are
             used to determine the angle value is the order that
             the atoms are in the psf/coord array. Recommend that
             BYNUm is used to specify the correct order of atoms.


   TRAJectory [ FIRStu int ] [ NUNIt int ] [ BEGIn int ] [ STOP int ]
                    [ SKIP int ] [ VELOcity ]  [first-atom-selection]
                        [ ORIEnt  [MASS]  second-atom-selection  ]


           { ALL                     } [P2] [UNIT int]
   SHOW    { time-series-name        }
           { CORRelation-function    }        (defines ?P2, ?AVER, ?FLUC)

           { ALL                     }
   EDIT    { time-series-name        }  edit-spec
           { CORRelation-function    }

           edit-spec::=  [INDEx int] [VECCod int] [CLASs int] [SECOnd int]
                               [TOTAl int] [SKIP int] [DELTa real]
                                   [VALUe real] [NAME new-name] [OFFSet real]

   READ  { time-series-name  } unit-spec edit-spec { [FILE]              }
         { CORRelation-funct }                     { CARD                }
                                                   { DUMB  [COLUmn int]  }

           { ALL                     }             { [FILE]              }
   WRITe   { time-series-name        }  unit-spec  { CARD                }
           { CORRelation-function    }             { PLOT                }
                                                   { DUMB [ TIME ]       }


   MANTIME time-series-name
               { DAVErage            } ! Q(t) = Q(t) - <Q(t)>,
                                         <Q(t)> implies time average
               { NORMal              } ! Q(t) = Q(t) / |Q(t)|
               { SQUAre              } ! Q(t) = Q(t) ** 2
               { COS                 } ! Q(t) = COS(Q(t))  (in degrees)
               { ACOS                } ! Q(t) = ACOS(Q(t)) (in degrees)
               { COS2                } ! Q(t) = 3*COS(Q(t))**2 - 1 (in degrees)
               { AVERage integer     } ! Q(t) = < Q(ti) >(ti=t-NUTIL+1,t)
               { SQRT                } ! Q(t) = SQRT(Q(t))
               { FLUCt name2         } ! print zero time fluctuations
               { DINItial            } ! Q(t) = Q(t) - Q(1)
               { DELN integer        } ! Q(t) = Q(t) - <Q(ti)>(ti=t-NUTIL+1,t)
               { OSC                 } ! print oscillations
               { COPY name2 [FIRSt int] [LAST int]        }
                                       ! Q(t) = Q2(t1), t1=FIRST,..,LAST
               { ADD  name2          } ! Q(t) = Q(t) + Q2(t)
               { RATIo name2         } ! Q(t) = Q(t) / Q2(t)
               { DOTProdcut name2    } ! Q(T) x-comp=Q(T).Q2(T)
                                         Q2(T)x-comp=angle Q(T) vs Q2(T) degrees
               { CROSproduct name2   } ! Q(T) = Q(T)xQ2(T)
               { KMULt name2         } ! Q(t) = Q(t) * Q2(t)
               { PROB integer        } ! Q(t) = PROB(Q(t))
               { HIST min max nbins  } ! Q(ibin) = Fraction of Q(t) values in ibin
               { POLY integer        } ! fit time series to polynomial (0-10)
                     [REPLace] [WEIGh name]
               { CONTinuous [real]   } ! make a (dihedral) time series continuous
                                         Q(t) = Q(t)+ n(t)*2*real, n(t)=integer
                                             (default real is 180.0)
               { MAP real1 real2     } ! Q(t)is mapped to [real1,real2]
                                       ! (typically [0,360] -
                                       ! default or both have to be specified!
               { LOG                 } ! Q(t) = LOG(Q(t))
               { EXP                 } ! Q(t) = EXP(Q(t))
               { IPOWer integer      } ! Q(t) = Q(t) ** integer
               { MULT   real         } ! Q(t) = real * Q(t)
               { DIVIde real         } ! Q(t) = Q(t) / real
               { SHIFt  real         } ! Q(t) = Q(t) + real
               { DMIN                } ! Q(t) = Q(t) - QMIN
               { ABS                 } ! Q(t) = ABS(Q(t))
               { DIVFirst            } ! Q(t) = Q(t) / Q(1)
               { DIVMaximum          } ! Q(t) = Q(t) / ABS(Q(MAX))
               { INTEgrate           } ! Q(t) = Integral(0 to t) (Q(t)dt)
               { MOVIng integer      } ! Q(t) = < Q(ti) >(ti=t-integer+1,t) (t)
               { TEST  real          } ! Q(t) = COS(2*PI*t*real/TTOT)
               { ZERO                } ! Q(t) = 0.0
               { DERIvative          } ! Q(t) = (Q(t+dt)-Q(t))/dt
               { SPHErical           } ! Q(t) = Q(t) 3-component vector series
                                       !      converted to spherical coord:
                                       !     (x,y,z)-> (r,phi,theta)

   CORFUN 2x(time-series-name)
             { [ PRODuct ]  [ FFT   ] [ LTC  ] [ P1 ] [ NONOrm ] } [ XNORm real ] [ TOTAl int ]
             {              [ DIREct] [ NLTC ] [ P2 ]            }
             {                                                   }
             {  DIFFerence                                       }

   SPECtrum  [FOLD] [RAMP] [SWITch] [SIZE integer]

   CLUSter time-series-name RADIus <real> [ MAXCluster <int> ] -
                            [ MAXIteration <int> ] [ MAXError <real> ] -
                            [ NFEAture <int> ] [ UNICluster <int> ] -
                            [ UNIMember <int> ] [ UNIInitial <int>] -
                            [ CSTEP <int> ] [ BEGIn <int> ] -
                            [ STOP <int> ] [ ANGLE ]

   END        ! return to main command parser



.. _correl_general:

General discussion regarding time series and correlation functions
------------------------------------------------------------------

The CORREL command invokes the CORREL subcommand parser.
The keyword values MAXTimesteps, MAXSeries, and MAXAtoms may be
specified for space allocation greater than the default options.
If there in insufficient virtual address memory for the space request,
it may be possible to achieve the desired results by removing the
nonbond lists before running the CORREL command.

The MAXTimesteps value is the largest number of steps any
time series will contain. The MAXSeries keyword is the largest number
of timeseries that will be contained at any time within CORREL.
A vector time series will counts as 3 time series in allocating space.
The MAXAtoms keyword allocates space for the atoms that are specified
in the ENTER commands (also duplicating a time series requires more space
for atoms). For bonds, angles, dihedrals, and improper dihedral
specifications, one extra value is needed for each entry to hold the
CODES value (so each bond uses 3 atom entries, 4 for angles...).

If the COVAriance keyword is given, no time series will be
computed, but instead, a complete equal time covariance matrix will
be computed. For this option, only one TRAJectory command is allowed.
The covariance matrix is then obtained by writing the time series, where
the elements are covariant with other time series.

The ENTER defines a time series. Many time series may be specified.
A time series is defined by the following items;

================  ==========================================================
Name              Each time series must have a unique (4 character) name.
Class code        The type of time series (BOND, USER, ATOM,...)
Number of steps   The number of time steps currently valid
Velocity code     Was the time series read from velocities?
Skip value        What multiple of delta do the time steps represent?
Delta             Integration time step
Offset            Time of first element
Secondary code    Depends on Class code (Geometry/Energy)(X/Y/Z...)
Vector code       1=simple time series, 3=vector, 0=Y or Z part of vector
Value             Utility series value, depends on Class code
Mass weighting    Are the elements to be mass weighted (only for ATOM)
Average           Time series average
Fluctuation       Time series fluctuation about the average
Atom pointer      Pointer into first specified atom in atom list
Atom count        Number atom entries given in the ENTER command
Time series       Series values from (1,NTOT)
================  ==========================================================

The TRAJectory command processes all of the time series which
have a NTOT (number of steps) count of zero. For this process,
the main coordinates are used for reading the trajectory. If fluctuations
are requested, the comparison coordinates MUST be filled with the
reference (or average) coordinates before invoking the TRAJectory
command. Allowing multiple TRAJectory commands separated by enter
commands make it possible to compute correlation function between
positions and velocities, or even for different trajectories.

The EDIT command allows the user to directly modify the time
series specifications.

The MANTIME command allows the user to manipulate the time
series values (and sometimes some of the specifications).

The SHOW command will display the specification data for all
of the time series.


.. _correl_enter:

Specifying time series
----------------------

The ENTER command defines a new time series. Each time series
specified by different enter commands must have a unique name (up to
4 characters). With this command, a time series may be defined and
then must be later filled with a TRAJectory command (or a MANTIME COPY,
or a READ time-series command). Alternatively, a time series may be retrieved
from an existing file, or duplicated from another time series that
currently exists.

The time series names "ALL" and "CORR" may not be used, and
are reserved for selecting all of the time series or the correlation
function respectivly.

The ENTER options are;

*  DUPLicate  time-series-name
   This causes an exact copy of an existing time series to be
   created (except with a different name). This may be useful where
   several different type of manipulations are required on a single
   time series.

*  READ  unit-number [CARD] [edit-spec]
   This causes a time series to be created and all data then
   read in from an existing time series file. All time series (up to the
   maximum allowed) will be read with this command.

*  Internal Coordinates
   ::

      [ BONDS   repeat(2x(atom-spec))        ] [ GEOMETRY ]
      [ ANGLE   repeat(3x(atom-spec))        ] [ ENERGY   ]
      [ DIHEd   repeat(4x(atom-spec)) [NOTR] ]
      [ IMPRo   repeat(4x(atom-spec))        ]

   These specifications cause a particular internal coordinate
   (or an average of several) to define the time series. It is not necessary
   that the specified atoms have a corresponding PSF entry, but if ENERGY is
   requested, the specified atoms must be able to produce a valid parameter
   code. The default is GEOMETRY. With geometry, any 4 atoms may be specified.
   A velocity trajectory should not be used to fill these types of time series.
   The NOTR option for dihedral prevents the analysis of dihedral transitions.

*  atom positions or velocities
   ::

      [ ATOM ] [ X   ]  repeat(atom-spec) [ MASS ]
      [ FLUC ] [ Y   ]
         [ Z   ]
         [ R   ]
         [ XYZ ]

   These ENTER commands define a time series, Q(t), based on atom
   positions or velocities. The ATOM option uses the (X,Y,Z,R,or XYZ) values
   directly.  The FLUCtuation option subtracts off the reference values
   (contained in the comparison coordinates). For example, if the average
   structure is desired as the reference value, then the command:

   ::

      COOR DYNA COMP trajectory-spec

   would be required before invoking the TRAJECTORY command.
   If more than one atom is specified, then Q(t) values are
   averaged over atoms.  If MASS is specified, then mass weighting is used in
   this averaging of Q(t) values.  The properties X,Y,Z, and R cause a scalar
   time series to be created with the requested property. The XYZ option causes
   a vector time series to be created.

   * ATOM:  Q(t) = X(t)
   * FLUC:  Q(t) = X(t) - Xref

*  Vector
   ::

      VECT   [ X   ]  repeat(2x(atom-spec))
       [ Y   ]
       [ Z   ]
       [ R   ]
       [ XYZ ]

   The VECTor command is similar to the ATOM and FLUCtuation
   commands listed above, except the values are given by the difference
   in position or velocity of 2 atoms. If more than one pair of atoms
   is specified, then the values for each vector are averaged.

   Q(t) = X1(t) - X2(t)

*  Vector product

   ::

      ATOM DOTProduct  repeat(2x(atom-spec))
      FLUC DOTProduct  repeat(2x(atom-spec))
      VECT DOTProduct  repeat(4x(atom-spec))

      ATOM CROSsproduct  repeat(2x(atom-spec))
      FLUC CROSsproduct  repeat(2x(atom-spec))
      VECT CROSsproduct  repeat(4x(atom-spec))

   These ENTER commands produce a scalar time series for
   velocities or positions with the following definitions;

   ::

      ATOM DOTP:  Q(t) =  ( r1(t) | r2(t) )
      FLUC DOTP:  Q(t) =  ( (r1(t)-r1(ref)) | (r2(t)-r2(ref)) )
      VECT DOTP:  Q(t) =  ( (r1(t)-r2(t)) | (r3(t)-r4(4)) )

   If more than one set of atoms is specified, then the vector values
   are averaged.  The dot product is then computed from the
   averaged vectors.  NOTE: the vectors are averaged, NOT the resultant
   dot products or cross products.   For the FLUC option, the reference
   coordinates must be in the comparison coordinate set.

*  Gyration

   ::

      [ GYRAtion ] [ CUT real ]
      [ DENSity  ]

   These commands define a scalar time series for a coordinate
   trajectory. The density calculation is based about the origin on all
   atoms within the CUT value; the radius of gyration is for all atoms
   within distance CUT of the geometric center of the molecule, and no
   mass weighting is applied.

*  MODE mode-number
   This option generates a scalar time series which is obtained
   by projecting the velocities onto the specified normal mode, or to
   project the coordinate displacement from the reference structure. The
   result is given by;

   * velocity:  Q(t) = < root(mass)*v(t) | q >
   * position:  Q(t) = < root(mass(i))*(r(t)-r(ref)) | q >

*  TEMPerature
   The time series is the temperature at each point.
   If NDEFG is specified as a positive value, then this is used instead of
   the NDEGF values from the trajectory file.  If a negative NDEGF value
   is specified, then NDEGF will be set to 3 times the number of selected
   atoms in the trajectory associated trajectory command.

*  HELIx atom-selection
   The x,y, and z components of the normalized vector defining the
   axis of a cylindrical surface best fitting the selected atoms.
   So you end up with a three-dimensional vector series.
   Intended for say alpha helices where the selection would be something
   like: ``SELE ATOM * * CA .AND. RESID 23:36 END``, to give the axis of
   an alpha helix running from residue 23 to residue 36.

*  INERtia atom-selection [ TRACe | ALL ]
   The x,y, and z components of the normalized vector defining
   the principal axis obtained from diagonalizing the moment of inertia
   tensor for the selected atoms at each time point.  The eigenvector
   corresponding to the smallest eigenvalue is returned, and 180 deg flips
   of the axis are explicitly prohibited (nonphysical).

   The optional TRACe keyword returns the sorted eigenvalues as a
   three column time series, instead of the principal axis vector.
   The optional ALL keyword (ALL and TRACe are mutually exclusive)
   returns all three principal axes as a vector with 9 components (x1,y1,z1,...)
   sorted with the main axis first.

   .. note::

      There may be problems, in particular for flexible systems, with
      exchange of the two minor axes; the code tries to correct for this
      (messages about this are printed at PRNLEV 7), but it may not always be
      right...

*  CELL  cell-spec
   If the cell-spec is one of the 6 unit cell parameters A, B, C,
   ALPHA, BETA, or GAMMA, then a single time series corresponding to that
   component is return.  The keyword ALL returns a 6 element time series,
   with the columns in the order given above.  The SHAPE keyword returns
   the shape matrix for the unit cell at each time point, in lower diagonal
   form.  The shape matrix has the angles as cosines, while ALPHA, BETA, and
   GAMMA are in degrees.

*  RMS  [ORIE]
   The RMS deviation from the COMPARISON coordinate set is
   computed for the atoms in the first selection on the TRAJ command,
   with a superposition to obtain a best fit to the same atoms in the
   COMParison coordinate set if ORIEnt is specified.
   If the TRAJ command also contains an ORIENT second_selection, this second
   selection will first have been used for a superposition onto the COMP
   coordinates.

*  PUCK RESI <resid> [SEGI <segid>]
   The sugar pucker phase and amplitude are calculated for
   the (deoxy)ribose of the specified residue; the first segment is
   the default. This gives a two-dimensional vector, with component 1
   being the phase (degrees) and component 2 the pucker amplitude
   (Angstroms), as defined by Cremer&Pople (JACS 1975).

*  ::

     PUCK ATOM [ Q ]      repeat(6x(atom-spec))
               [ THETa ]
               [ PHI ]

   Reports the Cremer & Pople puckering coordinates Q, THETa, and PHI for
   a six member ring of atoms. If Q, THETa, or PHI are not defined, all three
   coordinates are reported.

*  USER user-value [ repeat(atom-spec) ]
   The USRTIM routine is called for each coordinate or velocity
   set. The user value and atom list is also passed along. See the
   description in (USERSB.SRC)USRTIM for more details.

   Q(t) = Whatever you want!

*  SURF [RPRObe real] [WEIG] [ACCE|CONT] [RESI] [SELE atom-selection END]
   Computes the solvent accessible surface area vs time for the selected
   atoms in the context of the FIRST selection given to the TRAJ command. Uses the
   analytical method (see :doc:`SURF <corman>`).

   ========== ========= ===========================================================
   Keyword    Default   Meaning
   ========== ========= ===========================================================
   RPRObe     1.6       probe radius
   WEIG       .FALSE.   use WMAIN instead of LJ radii from parameter file
   ACCE|CONT  ACCE      accessible or contact surface
   RESI       .FALSE.   give ASA per residue in the selected set (creates a vector
                        time series with one component for each residue)
   ========== ========= ===========================================================

   Example:
   ::

      * Compute individual ASAs for 8 Trp residues in protein context given by all
      * residues with at least one atom within 8A of the Trp rings
      *
      ! r1 .. r8 are previously defined as 8 different Trp rings
      define trps sele r1 .or. r2 .or. r3 .or. r4 .or. r5 .or. r6 .or. r7 .or. r8 end
      define environment sele .byres. (segid cht .and. ( trps .around. 8.0 ) ) end
      long ! allows all ASA values at each time point to be written on one line
      correl maxseries 10 maxtime 50000 maxatom 200
      enter asa surf rprobe 1.4 sele trps end resi
      traj firstu 51 nunit 1 begin 100000 skip 500 sele environment end stop 125000
      write asa dumb time unit 21
      *hi
      *
      end


*  TIME [ AKMA ]
   The time is returned in picoseconds unless AKMA is specified.

   * Q(t) = t

*  ZERO

   A zero time series is specified ( Q(t)=0 ).
   This option is useful for cases where time series will be read with
   the DUMB option. For these cases, the EDIT command may also be needed
   to get desired results.

*  DIPO [ SELE atom-selection END ]
   Computes the dipole moment of all atoms specified in the atom
   selection. The OXYZ and MASS keywords have the same meaning as defined
   in COOR DIPO. See :doc:`corman` for further details.

*  VECM [ SELE atom-selection END ]
   Generates a series like VECT XYZ, but IMAGE aware (which need to
   be set up appropriately). If CUTIM is chosen appropriately (e.g., L/2
   for a cubic box), the vector in the time series will always represent the
   minimum image pair of the two atoms.

*  Dipole moment of a water/solvent shell
   ::

      SDIP [ SHELL int ]
      [ BULK ]

   Computes the dipole moment of a water/solvent shell. Returns
   X/Y/Z and the number of atoms in the shell.
   See :doc:`shell` for further details.
   The OXYZ and MASS keywords have the same meaning as defined in COOR DIPO.
   See :doc:`corman`. for further details.

*  Shell

   ::

      SATM [ SHELL int ] atom-spec
      [ BULK ]

   The series contains zero or one depending on whether the atom is
   in the specified shell (or the bulk). See :doc:`shell`.
   for further details.

*  PRDI int [ MASS ]
   This tree-dimensional time series contains the sum of all
   single dipole moments for each set in a given prototype set (see
   :doc:`proto`). This differs from the overall dipole moment
   for all sets only if the single sets carry a net charge. In this case
   the dipole moment of each set is calculated relative to a given
   reference point. If the MASS keyword is present, this point of
   reference is the center of mass of a given set, while in its absence
   the center of geometry is used. (Note: Almost equivalent functionality
   can be obtained with the DIPO series.)

*  PRVE int [ MASS ]
   Is similar to PRDI but calculates the sum of the center of
   geometry (or center of mass with keyword MASS) velocities of a given
   prototype set.

*  OMSC ETA real

   The series computes the cumulative Onsager-Machlup score
   (:ref:`Onsager-Machlup score <dims_onsager_Machlup_score>`).
   ETA is the friction coefficient of the dynamics (in 1/ps). As
   a first guess one may use the value used for the Langevin dynamics
   ('FBETA').

   OMSC can only be used as the single time series in a CORREL
   command. In particular, it is incompatible with RMS because they both
   use the same reference array for different things (RMS stores the
   comparison structure, OMSC the previous frame to compute velocities
   X(t) - X(t-1).)

   Output:
   The standard output (at PRNLEV 3 or higher) consists of lines

   ::

     OMSCORE> step-score normalized-cumulative-score

   The OM score for the first step is calculated and used to normalize
   all following scores. The numbers can become rather large and using
   the normalized score avoids using LONG in the output. Otherwise the
   output format overflows and only ******** would be printed.

   With the CORREL WRITE command, the normalized-cumulative-score for N-1
   steps is written to an CORREL output file. The first step contains the
   normalization factor s(t=0). You may have to postprocess the file
   (using for instance awk) after having written the output file
   omscore.dat with CORREL's WRITE name DUMB TIME ...:

   ::

     awk 'NR == 1 {s0 = $2}; {t=$1; s=s0*$2; print t," ",s}' \
         omscore.dat > omscore_nn.dat

   Note that it only makes sense to compare OM-scores for trajectories of
   the same system and of the same length.


.. _correl_trajectory:

Specification of the Trajectory Files
-------------------------------------

The TRAJectory command reads a number of trajectory files whose
Fortran unit numbers are specified sequentially. The first unit is given
by the FIRSTU keyword and must be specified. NUNIT gives the number of
units to be scanned, and defaults to 1.

BEGIN, STOP, and SKIP are used to specify which steps in the
trajectory are actually used. BEGIN specifies the first step number to
be used. STOP specifies the last. SKIP is used to select steps
periodically as follows: only those steps whose step number is evenly
divisible by STEP are selected. The default value for BEGIN is the first
step in the trajectory; for STOP, it is the last step in the trajectory;
and for SKIP, the default is 1.

The first atom selection in the TRAJectory command is meaningful
only for those time series that require an atom selection.  These are
time series defined by the following ENTER commands: GYRAtion, DENSity,
RMS, MODE, TEMPerature, and optionally USER.

General reorienting of a coordinate trajectory is possible using the
MERGE command. For details see :ref:`Merge <dynamc_reorient>`.
It is also possible to perform a simple rms best fit of each frame with the
reference coordinates (comparison set) using the ORIEnt option.  For this
option a second atom selection MUST be provided and a MASS keyword is an
option that allows for a mass weighting of the best fit. This superposition is
performed before any other manipulation on each frame to be analyzed.

If VELOcity is specified, a velocity trajectory will be looked
for. Otherwise, a coordinate trajectory is expected.

Any time series that has a zero count (NTOT=0) will be
filled by this command. The time series count will then be filled
with the total number of steps processed for each of these series.
Any time series with a nonzero count (NTOT>0) will not be affected
by this command. The count may be set to zero for a time series with
the EDIT command.

Upon conclusion, the average and fluctuation as well as some
other data is presented on each of the processed time series.

If any of the time series to be filled require a reference
coordinate set, then the comparison coordinates MUST be filled with the
reference (or average) coordinates before invoking the TRAJectory
command. Upon completion, the main coordinates contain the last coordinate
set read from the trajectory, and the comparison coordinates are unaffected.


.. index:: correl; edit
.. _correl_edit:

Editing a time series
---------------------

The EDIT command allows the time series specifications
to be modified directly.

.. warning::

   This command gives the user direct access to most time
   series specification. There is NO checking to see if what is being done
   makes sense. As such, this command is versatile and dangerous.

The EDIT command must be followed by a valid time series name.
All subsequent keywords will be based on that time series.
The series name "ALL" will cause the edit spec to operate on all
the time series. The name "CORR" will cause the edit to occur on the
correlation function.

The following may be specified for a time series;

=============== ============================================================
INDEx integer   May be specified to modify X,Y, or Z (1,2,3 resp)
                of a vector time series. Otherwise, all are modified.
                The index number is in fact an offset from the specified
                time series, where a value of 1 represents the selected
                time series. A value of 5 will cause the edit operation
                to modify the fourth time series from the specified.

CLASs integer   May be used to specify a class code (consult source).

TOTAl integer   The total number of valid steps may be altered, but
                none of the values are modified. By setting this
                value to zero, the time series is then ready again
                for the next TRAJectory command.

SKIP integer    May be specified to reset the SKIP value. This may be
                useful after reading an external time series.

DELTa real      May be specified to modify the basic time step. The
                actual time step for a series is (SKIP*DELTA).

OFFSet real     The time of the first element in the time series.

VECCod integer  User may specify a vector code. This may be useful
                in merging 3 separate time series into a vector
                time series (or the reverse). In fact any number of
                time series may be grouped together with this option.
                For example, if a table with 5 time series is desired,
                setting VECCOD to 5 for the first one and the writing
                this time series will output all 5.


VALUe real      This allows the user to modify the series utility
                value. Its function depends on the Class code.
                This value is currently used for (USER, GYRAtion,
                DENSity, MODE, and TIME)

SECOndary int   The secondary class code may be modified (consult source).
=============== ============================================================


.. index:: correl; mantime
.. _correl_mantime:

Manipulating the Time Series
----------------------------

The MANTIME command allows the user to manipulate selected
time series, Q(t), and performs the operation requested by the option
and leaves the resultant time series as the active time series.
This helps in performing various permutations of manipulations to increase
the options without increasing the number of ENTER commands.

The keyword ordering must be followed exactly.

=================== ===================================================================
DAVErage            subtracts the average of the time series from all elements.

NORMal              normalizes the vectorial time series.
                    (i.e. creates the unit vector by dividing all elements for
                    a given value of t by :math:`r(t) = \sqrt{x^2 + y^2 + z^2}` ).

SQUAre              squares all the elements

COS                 obtains the cosine of all elements.

ACOS                obtains the arc-cosine of all elements.

COS2                calculates 3*cos**2 - 1 for all elements.

AVERage integer     calculates the average for every <integer> consecutive
                    points and increases the time interval by a factor of
                    <integer>. Note: NTOT is divided by <integer>.

SQRT                obtains square root for all elements.
                    Negative elements are set to -SQRT(-q(t)).

FLUCt name          The Q(t) remains unchanged.
                    A second (b) time series must be specified.
                    The zero time fluctuations are computed and printed
                    out.  The following variables are computed:

                    * A = :math:`\langle Q_a(t) \cdot Q_b(t) \rangle`
                    * B = :math:`\sqrt {\langle Q_a(t)^2 \rangle}`
                    * C = :math:`\sqrt {\langle Q_b(t)^2 \rangle}`
                    * D = :math:`A/(B*C)`

DINItial            subtracts the value of the first element from all elements.
                    Q(t) = Q(t) - Q(1)

DELN integer        Q(I) = Q(I) - <Q(I)> I FROM 1 TO N, FROM N+1 TO N+N ETC.
                    (untested).

OSC                 counts the number of oscillations in Q(t) / unit time step.
                    The Q(t) remains unchanged.


COPY name           This copies the second time series to the first. NTOT
                    of the first is set to that of the second. If FIRSt or LAST is
                    specified, a subset (I=FIRST,,,LAST, with a total of
                    FIRST-LAST+1 points) of the second series is copied.
                    Defaults for FIRSt and LAST are 1 and NTOT of the second
                    series.

ADD name            Q(t) = Q(t) + Q2(t); add the second time series to the first

RATIo name          Q(t) = Q(t) / Q2(t)

CROSsprod name      :math:`Q(t) = Q(t) \times Q2(t)`; the 3D crossproduct of the two
                    3D vectors formed by the selected and named time series

DOTProd name        ::

                      Q(T) = x-comp of Q(T)= Q(T) . Q2(T)
                             x-comp of Q2(T) angle in degrees between the two vectors
                      NOTE! Modifies Q2 as well as Q
                      to get just the x-comp you may then edit the selected series:
                      EDIT series VECCOD 1

KMULt name          Q(t) = Q(t) * Q2(t)

PROB integer        give the probability to find a specific value of the
                    time series. <integer> subdivisions of the time series
                    are considered so that there are integer+1 values.

HIST min max nbins  Q(ibin) = Fraction of Q(t) values within ibin
                    This command replaces a time series with a
                    histogram of the time series divided into "nbins" with
                    a range from "min" to "max".  The histogram values sum to 1.

POLY integer        fit time series to polynomial. The order should
[REPLace]           be in the range of 0 to 10.
[WEIGh name]
                    * Order 0 will provide just the average,
                    * Order 1 will fit the time series to a straight line.
                    * Order 2 will fit to a quadratic function.

                    The REPLace option will replace the time series with
                    fitted one.  The WEIGht option will wait all data
                    by the values in a second time series.

CONTinuous real     Q(t) = Q(t) + n(t) , where n(t) is an integer such that
                    the ABS(Q(t)-Q(t-1))<=real

                    The default value is 180.0, which is appropriate for
                    making a dihedral time series continuous.  A different
                    positive value may be selected (such as a box size...).

LOG                 Q(t) = LOG(Q(t))

EXP                 Q(t) = EXP(Q(t))

IPOWer integer      Q(t) = Q(t) ** integer

MULT real           Q(t) = Q(t) * <real>

DIVI real           Q(t) = Q(t) / <real>

SHIFt real          Q(t) = Q(t) + <real>

DMIN                Q(t) = Q(t) - QMIN, QMIN being the minimum of the time series.

ABS                 Q(t) = ABS(Q(t))

DIVFirst            Q(t) = Q(t) / Q(1)

DIVMax              Q(t) = Q(t)/ ABS(Q(t) with max norm)

INTEgrate           Q(t) = Integral(0-t) [ Q(t) dt ]

MOVIng integer      Q(t) = Q(t) = < Q(ti) >(ti=t-integer+1,t) (t)
                    At each time, computed the moving average of the last
                    <integer> points.  It <integer> is zero or negative, the
                    moving average is taken over all the preceding points.

TEST real           Q(t) = COS ( 2 * PI * <real> / NTOT )

ZERO                Q(t) = 0
                    This option zeroes the specified time series.

DERIvative          Q(t) = (Q(t+dt)-Q(t))/dt, the last point is set to the one
                    before last
=================== ===================================================================


.. _correl_corfun:

Calculating a Correlation Function
----------------------------------

CORFUN: This option takes the specified time series and calculates the
desired correlation function from it.  The resultant correlation function
is saved in a time series named "CORR" which may then be used in subsequent
CORREL manipulation or write commands.  If multiple CORFUN commands are
requested, then the "CORR" time series is overwritten.
Command line substitution parameter CFNORM is set to the value that would be
used as the multiplicative normalization factor of the correlation function.

In the following, Qa and Qb refer to the time series that were
extracted using the CORREL command.

=============== ======================================================================
PRODuct         This option (default) generates a correlation function that is the
                product of the time series elements.
                C(tau) = < Q1(t)*Q2(t+tau) >

DIFFerence      The difference option is an alternative of the product option
                and it generates a function that is useful in calculating
                diffusion constants (slope at long tau).

                :math:`C(\tau) = \langle (Q_1(t) - Q_2(t+\tau))^2 \rangle`

FFT             This option is to calculate the correlation function using the FFT
                method.  There are certain limitations on the prime factors
                in the total number of points.

DIRECT          This option is to calculate the correlation function using the
                direct multiplication method.

P1              This option gives the direct correlation function, <Qa(0).Qb(t)>.
                If Qa and Qb are unit vectors, then this is also the first
                order Legendre Polynomial

P2              This is to obtain the correlation function of second order Legendre
                Polynomial, (3 <[Qa(0).Qb(t)]**2> - 1)/2.  For all applications
                that I can think of, Qa and Qb will be unit vectors. For P2, LTC = 0
                and NORM = 1

NLTC            no long tail correction.

LTC             long tail correction (subtracts <Qa>**2 if autocorrelation,
                <Qa>*<Qb> if cross correlation.  There is no LTC for P2
                so NLTC and LTC give same result.)
                This feature is to be used with care.  If the Qa and Qb are
                fluctuations from the mean (i.e. FLCT or MANTIME DELTA), then
                this can serve as a correction for roundoff error.  Otherwise,
                they are not centered about the mean, this correction causes
                the C.F. to be a less accurate calculation of fluctuations from
                the mean, i.e.

                .. math::

                   \langle Q_a(0) \cdot Q_b(t) \rangle - \mathrm{LTC}
                        =& \langle Q_a(0) \cdot Q_b(t) \rangle -
                           \langle Q_a \rangle \langle Q_b \rangle \\
                        =& \langle \Delta Q_a(0) \cdot \Delta Q_b(t) \rangle

NONORM          Correlations are not normalized. This is useful for adding
                correlations computed in different trajectories.
                (P2 is not normalized)

                The correlation functions are normalized unless NONORM is specified.

XNORm           Use this value if not zero as normalization factor (multiplies all
                values in correlation function). Overrides NONORM setting.

TOTAL integer   The TOTAL value determines the number of points to keep in
                the correlation function. The number of points may not be
                grater than the number of points in the time series. A reasonable
                value is about 1/4 to 1/3 the length of the time series.
                Correlation function values near the end have little weight.
                The default value is the nearest power of two less than half of
                the time series length.
=============== ======================================================================


The defaults are FFT, P1, NLTC.

.. note::

   The correlation time which is given by the program is calculated
   by an exponential fit to the first NTOT/8 points or up to the
   first crossing of the time axis.  This value should be considered
   a (poor) estimate, it is meaningful only for correlation functions
   which decay exponentially to zero with no oscillations.

For P1,
C(t) = (c(t) - ltc)/N
ltc and Normalization factors, N, are:

*  LTC, autocorrelation:

   .. math::

      \mathrm{LTC} &= \begin{cases}
         \langle Q_a \rangle ^2 &\; \text {for P1} \\
         0 &\; \text{for P2}
         \end{cases} \\
      N &= C(0) - \mathrm{LTC} \\
        &= \langle Q_a^2 \rangle - \mathrm{LTC}

*  LTC, cross-correlations:

   .. math::

      \mathrm{LTC} &= \langle Q_a \rangle \langle Q_b \rangle \\
      N &= \sqrt{ (\langle Q_a^2 \rangle - \langle Q_a \rangle ^2) (\langle Q_b^2 \rangle - \langle Q_b \rangle ^2) }


*  NLTC, autocorrelation:

   .. math::

      \mathrm{LTC} &= 0 \\
      N &= C(0)

*  NLTC, cross-correlations:

   .. math::

      \mathrm{LTC} &= 0 \\
      N &= \sqrt{ \langle Q_a^2 \rangle \langle Q_b \rangle ^2}


.. _correl_spectrum:

Generating a Spectrum from Correlation Functions
------------------------------------------------

There is a command, SPECtral-density, which may be used to generate
a spectrum from a correlation function. The syntax is;

::

   SPECtrum [SIZE integer] [FOLD] [RAMP] [SWITch]


.. _correl_cluster:

Clustering Time Series Data
---------------------------

This command clusters time series data obtained within the CORREL
facility.  The time series must first be defined using CORREL's ENTEr
command and the data read in via TRAJ or READ.  The CLUSter command
clusters these data into groups with similar time series values, with
each cluster being defined by a "cluster center".  The cluster centers are
output to UNICluster, and a list of time points and assigned clusters is
given in the cluster membership file (UNIMember).

For example, if you want to find similar conformations of a peptide
using dihedral angles, you would first define the set of dihedral angles to
be considered, say angle(1) -> angle(M), as M time series.  If the time series
were each N time steps long, then you would be clustering N "patterns", with
each pattern M "features" long.

Consecutive time series are clustered.  If the first time series
is, for example, "ts1" then the "veccod" of this time series can be
changed to the number of time series to be clustered:

::

   CORREL ...
       ENTE ts1 ...
       ENTE ts2 ...
       ...
       ENTE tsM ...
       EDIT ts1 veccod M
       TRAJ ... (or READ ...)
       CLUSTER ts1 ...
   END

Alternatively, NFEAture M can be specified in the CLUSter command line.
Note that vector time series count as three features.


The Clustering Algorithm
^^^^^^^^^^^^^^^^^^^^^^^^

ART-2' is a step-wise optimal clustering algorithm based on a
self-organizing neural net (Carpenter & Grossberg, 1987; Pao, 1989;
Karpen et al., 1993).  The algorithm optimizes cluster assignment subject
to a constraint on cluster radius, such that no member of a cluster is more
than a specified distance from the cluster center.  This optimization is
carried out as an iterative minimization procedure that minimizes the
Euclidean distance between the cluster center and its members.

A self-organizing net is created with each output node representing
a cluster.  The number of pattern features is equal to the number of input
nodes.  The weights of the connections between the input layer (layer i)
and the output layer (layer j) are denoted by b(j,i).  For each cluster j,
b(j,i), i = 1, nfeature, is the cluster center.  To create the net (which is
synonomous to learning the classification scheme or cluster centers) the
following algorithm is implemented:

1. To initialize the network, assign b(1,i) equal to the first
   pattern tq(1,i) for i = 1, nfeature.

2. For each pattern number k, calculate the Euclidean distance (rms)
   between the pattern tq(k,i) and all cluster centers b(j,i), where
   j is the cluster index.

   rms(j,k) = sqrt[sum [(b(j,i)-tq(k,i))**2] for i = 1, nfeature]

3. Find cluster j such that rms(j,k) < rms(i,k) for all i<>j.  If
   rms(j,k) <= Threshold, then update b(j,i):

   b(j,i) = ((m-1)*b(j,i) + tq(k,i))/m,

   where m is the number of prior updates of b(j,i).  Note that
   b(j,i) is the average of feature i for all patterns currently
   assigned to cluster j.

4. If rms > Threshold for all prior cluster centers (j=1,numclusters),
   then create a new cluster center by increasing the number of
   output nodes by one, and assign the weights b(numclusters,i) of
   this node the value of the pattern tq(k,i).

5. Repeat 2.-4. until all patterns have been input.

6. Compare the new set of cluster centers with the last set.  If
   the difference between them is less than MAXError, then halt
   clustering.

7. If the difference between the sets of cluster centers is greater
   than MAXError, then use the new set of cluster centers as the
   starting cluster centers, and repeat steps 2.-6.  Else, clustering
   is complete.

Note that the cluster centers currently being calculated in step 3
are only used for the comparison in step 2 during the first
iteration with no initial cluster centers.  Otherwise, the centers
calculated in the previous iteration (or read from UNIInit) are
used in the comparison in step 2.  Hence, in the initial "learning"
phase, cluster centers are recalculated as each new member is added.
In subsequent "refining" phases, cluster centers are not updated
until all conformations are read in and assigned.

References:

1) Carpenter, G. A., & Grossberg, S. 1987. ART 2: Self-organization of stable
   category recognition codes for analog input patterns. Appl. Optics  26:4919-
   4930.

2) Pao, Y.-H. 1989. Adaptive Pattern Recognition and Neural Networks, Addison-
   Wesley, New York.

3) Karpen, M. E., Tobias, D. T., & Brooks III, C. L. 1993. Statistical
   clustering techniques for analysis of long molecular dynamics trajectories.
   I: Analysis of 2.2 ns trajectories of YPGDV. Biochemistry  32:412-420.


CLUSter Parameters
^^^^^^^^^^^^^^^^^^

::

   CLUSter time-series-name RADIus <real> [ MAXCluster <int> ] -
                            [ MAXIteration <int> ] [ MAXError <real> ] -
                            [ NFEAture <int> ] [ UNICluster <int> ] -
                            [ UNIMember <int> ] [ UNIInitial <int>] -
                            [ CSTEP <int> ] [ BEGIn <int> ] -
                            [ STOP <int> ] [ ANGLE ]


1.  time-series-name: The name of the first time series (as defined by
    the ENTE command) to be clustered.

2.  RADIus: Maximum radius of cluster.  The rms cutoff or threshold for
    assignment to a cluster.

3.  MAXCluster:  Maximum number of clusters (default = 50).

4.  MAXIteration:  The maximum number of iterations allowed.  If the
    clustering has not converged by this number of iterations, all
    clusters are output (default = 20).

5.  MAXError:  If the rms difference between the position of the cluster
    centers for the last two iterations is less than maxerror, the system
    is considered converged and the clustering is halted (default = 0.001).

6.  NFEAture:  This variable gives the number of features in the input
    pattern, that is, the number of time series to be clustered at a time.
    The default is the veccod parameter associated with 'time-series-name'.
    NFEATure time series are clustered, starting with 'time-series-name'
    and continuing with the next nfeature-1 series specified in subsequent
    'ENTE' commands (default = veccod of time-series-name).

7.  UNICluster:  The unit number of the output cluster file.  If UNIC = -1
    (the default), the cluster parameters are output to the standard output.

8.  UNIMember:  The unit number of the output membership file.  This file
    lists each time point and the cluster(s) associated with the specified
    time series at that time point.  If UNIM = -1 (the default), the
    membership list is not output.

9.  UNIInit:  The unit number of the file with the initial cluster centers.
    If UNII = -1 (the default), no initial cluster centers are specified.

10. CSTEp:  This variable gives the spacing between time series in the
    input vector.  For each timepoint k, the set of patterns clustered is
    tq(k,1) -> tq(k,nfeature), tq(k,1 + cstep) -> tq(k,nfeature + cstep),
    ...,tq(k,nserie - nfeature + 1) -> tq(k,nserie) (default = nfeature).

11. BEGIn:  Indicates frame in time series where clustering begins
    (default = 1).

12. STOP:  Indicates the frame in the time series where clustering ends
    (default = minimum length (TOTAl in SHOW) of time series).

13. ANGLe:  A logical flag which when true specifies angle data is to be
    clustered, taking angle periodicity into account (default = .FALSE.).


Caveats
^^^^^^^

The clustering algorithm is initial-guess dependent, i.e., it is
dependent on the input order of the patterns.  The order of presentation
in CLUSter is simply the consecutive frames of the time series.  To check
for stable clustering, cluster centers can be calculated from time series
with the time frames randomized.  This is not currently implemented in
CHARMM, so the user will have to write a set of time series to a file
and then randomize row position outside of CHARMM.

It is relatively straight forward to compare features derived from
similar measures (i.e., time series with the same "class codes", for
example all DIHE/GEOM).  In some applications it may be desired to "mix"
units in the pattern, for example, cluster a set of time series derived
from both atomic positions and energies.  How best to compare "apples &
oranges" is a problem from measurement theory, and is application-specific.
Normalizing the variables such that they have unit variance is one
possibility, and this can be done by 1) determining the standard deviation
of the time series (FLUC given by the SHOW command), and 2) using this
value in the MANTim DIVI command.  Since only differences between features
are used in the clustering algorithm, shifting the time series to zero
mean is not necessary.

Duda & Hart have a good discussion of the issues involved in
clustering and normalization:

Duda, R. O., & Hart, P. E., Pattern Classification and Scene Analysis,
Wiley, New York, pp. (1973).


Cluster Output
^^^^^^^^^^^^^^

The following data are output to UNIC for each cluster:

* Cluster Index - The clusters are numbered starting with 1.
* No. of Members - Number of patterns assigned to the cluster.
* Cumulative No. of Members - The total number of patterns within the
  cluster radius.  This can be higher than the No. of Members due
  to patterns being within the maximum radius of more than one cluster.
* Standard Deviation of Patterns within Cluster -
  For cluster j with the number of features = Nfeature, this is
  sqrt(sum((tq(k,i) - b(j,i))**2)/Nfeature*N(j)) where the sum is
  over i = 1, Nfeature and over all k such that tq(k) is a member
  of j.  N(j) = the number of members in cluster j.  Note that
  b(j,i) = <tq(k,i)> (averaged over k in cluster j).
* Maximum Distance - the longest distance between the cluster center and
  an assigned pattern, normalized by sqrt(Nfeature).
* Cluster Centers - (b(j,i), i = 1, Nfeature)

The following data are output to UNIM:

* Cluster index of the assigned cluster
* Time series time step
* Time series index of first time series in pattern
* Distance of pattern from cluster center, normalized by sqrt(Nfeature)


.. _correl_io:
 
Input/Output of time series and correlation functions
-----------------------------------------------------

1) The SHOW command

   ::

              { ALL                     }
      SHOW    { time-series-name        }
              { CORRelation-function    }

   The SHOW command displays to print unit various data regarding
   the specified time series. This command is automatically run after the
   ENTER and EDIT commands as a verification of the last action.


2) The READ command

   ::

      READ  { time-series-name  } unit-spec edit-spec { [FILE]              }
            { CORRelation-funct }                     { CARD                }
                                                      { DUMB  [COLUmn int]  }

   The READ command allows a time series or correlation function to
   be directly read. The file formats for time series and correlation
   functions is identical. There are three basic methods by which time
   series may be read: FILE (default), CARD, and DUMB. The FILE and CARD
   options expect a file of specific type generated by the corresponding
   WRITE command. The DUMB option will read a free field card file with
   NO title or other header. The COLUmn option (default 1) may be specified
   to start reading the time series from any specified column. The DUMB
   option will usually include some edit specifications to properly set
   the time steps (etc.).

3) The WRITe command

   ::

              { ALL                     }              { [FILE]        }
      WRITe   { time-series-name        }  unit-spec   { CARD          }
              { CORRelation-function    }              { PLOT          }
                                                       { DUMB [ TIME ] }

   The WRITe command will write out time series or a correlation function.
   All of the write options expect a title to follow this command.
   There are several file formats; FILE (default), CARD, PLOT, and DUMB.
   The FILE and CARD options will write out all data regarding the specified
   time series with the expectation for later retrieval by CHARMM or another
   program. The PLOT option will create a BINARY file for plotting by PLT2.
   The first line of the title is used as the plot title, but this may be
   reset in PLT2.

   The DUMB options will simply write out the values with no title
   or header to a card file, one value to a line. If the TIME option is
   specified, the time value will precede the time series values (as needed
   for an X-Y plot). If the time series is a vector type, then all component
   values will be given on each line. Unless LONG (see miscom.doc) is in effect
   the output is limited to 8 values/line.  DUMB option is useful for making plot
   files, or for feeding the data to other programs.

   With the EDIT command, a user may merge 3 separate sequential
   time series into a vector time series (or the reverse). In fact any number
   of time series may be grouped together with this option.  For example,
   if a table with 5 time series is desired, setting VECCOD to 5 for the
   first one and the writing this time series will output all 5.


.. index:: correl; examples
.. _correl_examples:

Examples
--------

These examples are meant to be a partial guide in setting up
input files for CORREL. The test cases may be examined for a wider
set of applications.

Example (1)
^^^^^^^^^^^

::

   CORREL MAXSERIES 1 MAXTIMESTEPS 500 MAXATOMS 5
   ENTER AAAA  TORSION MAIN 28 N MAIN 28 CA MAIN 28 C MAIN 29 N    GEOMETRY
   TRAJECTORY FIRSTU 51 NUNIT 5 BEGIN 26000 STOP 31000 SKIP 10
   MANTIME AAAA DAVER
   WRITE AAAA UNIT 20 DUMB TIME
   * title
   *
   WRITE AAAA CARD UNIT 10
   * title for card
   * file containing the time series
   *
   CORFUN AAAA AAAA   FFT NLTC P0
   WRITE CORREL  UNIT 21 DUMB TIME
   * title
   *
   WRITE CORREL FILE UNIT 11
   * title for binary correlation function
   *

* Extracts the time series, PHI(t), for phi dihedral of residue 28.
* Makes the time series the fluctuation from the mean, delta PHI(t).
* Makes a plot file of delta PHI(t) vs. time.
* Makes binary file of delta PHI(t).
* Calculates C(t) = <delta PHI(0) . delta PHI(t)> / <PHI**2> by FFT
  with no long tail correction.
* Makes a plot file of C(t) vs. t.
* Makes a binary file of C(t).

Example (2)
^^^^^^^^^^^

::

   CORREL MAXSERIES 2 MAXTIMESTEPS 500 MAXATOMS 10
   ENTER PHI  TORSION MAIN 27 C  MAIN 28 N  MAIN 28 CA  MAIN 28 C  GEOMETRY
   ENTER PSI  TORSION MAIN 28 N  MAIN 28 CA MAIN 28 C   MAIN 29 N  GEOMETRY
   TRAJECTORY FIRSTU 51 NUNIT 5 BEGIN 26000 STOP 31000 SKIP 10
   MANTIME PHI DAVER
   MANTIME PSI DAVER
   CORFUN  PHI PSI  FFT NLTC P0 NONORM
   WRITE CORREL FILE UNIT 11
   * title for cross correlation binary file
   *
   WRITE CORREL PLOT UNIT 12
   * plot title
   *

* Extracts the time series PHI(t), for phi dihedral, and PSI(t), for
  the psi dihedral, of residue 28.
* Makes the time series the fluctuation from the mean.
* Calculates C(t) = <delta PHI(0) . delta PSI(t)> by FFT with no
  long tail correction.
* Makes a binary file of C(t).
* Makes a binary PLT2 file for plotting

Example (3)
^^^^^^^^^^^

Fluorescence Depolarization, for example

::

   CORREL MAXSERIES 6 MAXTIMESTEPS 500 MAXATOMS 8
   ENTER V1  VECTOR XYZ  MAIN 28 NE1 MAIN 28 CZ3 MAIN 28  NE1 MAIN 28 CE3
   ENTER V2  VECTOR XYZ  MAIN 28  CD1 MAIN 28 CH2 MAIN 28 CD1 MAIN 28 CZ3
   TRAJECTORY FIRSTU 51 NUNIT 5 BEGIN 26000 STOP 31000 SKIP 10
   MANTIME V1 NORMAL
   MANTIME V2 NORMAL
   SHOW ALL
   CORFUN  V1 V2   FFT P2
   WRITE CORREL PLOT UNIT 21
   * title for plot
   *

* Extracts the time series, consisting of the average of the vectors
  NE1 - CZ3 and NE1 - CE3 == V1(t) and of the average of CD1 - CH2 and
  CD1 - CZ3 == V2(t).
* Makes V1(t) and V2(t) unit vectors.
* Displays data regarding both time series
* Calculates P2(t) = (3< (V1(0)*V2(t))**2 > - 1) / 2
* Makes a binary plot file for PLT2
