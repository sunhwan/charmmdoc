.. py:module::vibran

==================
Vibration Analysis
==================

The vibrational analysis section of CHARMM has been designed
to be a general purpose normal mode generation and analysis facility.
Also included is an extensive set of vector analysis and comparison
features and entropy calculation.

Support programs such as the iterative diagonalization program,
or the restartable large matrix diagonalization program are compatable
with this facility. Also included are routine to generate trajectories
which can be used for examining modes on the picture system.

In order to process commands with the vibrational analysis routines,
The energy terms must all be defined, and the structure must be determined
(see :ref:`energy_needs`). At present, SHAKE and images 
(see :doc:`images`) are not supported. 
Systems with atoms fixed by the CONStraint FIX command can be treated 
by using the REDU FIX option.

Within the vibrational analysis command mode, all miscellaneous
(MISCOM), coordinate manipulation (CORMAN), and internal coordinate (IC)
commands are allowed.

Keywords used to define Hydrogen bonds and nonbonded interactions
may be included in the command that invokes VIBRAN.


.. _vibran_syntax:

Syntax for vibrational analysis command and subcommands
-------------------------------------------------------

Main command:

::

   VIBRan  [hbond-spec] [nbond-spec] [nmode-spec] [IHBFRQ 0] [INBFRQ 0]

   nmode-spec ::= NMODe integer
           (this specification allocates space on the heap and its    )
           (default is NATOM*3. For large systems, if normal modes are)
           (to be used, it should be set to the largest number needed.)
           (Its default value is set to 1 if NATOM is greater than 50)

   hbond-spec        !  see *note hbonds:(chmdoc/hbonds.doc).
   nbond-spec        !  see *note nbonds:(chmdoc/nbonds.doc).

Subcommands:

::

   miscellaneous-command-spec         ! see *note miscom:(chmdoc/miscom.doc).

   IC        ic-subcommand-spec       ! see *note intcor:(c22doc.intcor.doc).

   COOR      coor-subcommand-spec     ! see *note corman:(chmdoc/corman.doc).

   READ    { NORMal-modes  [FILE] unit-spec  [APPEnd] [mode-spec]   }

   WRITe   { NORMal-modes  [mode-spec] [CARD]                          } unit-spec
           { SECOnd-derivatives  [RAISe] [MASS]
                   [CARD [FINIt [STEP val] [TOL val] [atom-selection] ]] }
           { TRAJectory  mode-spec magnitude-spec [PHAS real] [SHAKe]
                   [SEQUential-files] [NCYC integer] [STEP real]
                   [SUPErpose] [RANDom] [ISEEd integer] }

   PRINt   NORMal-modes  [mode-spec] [magnitude-spec] [INTDer [FINIt]]
                   [VECTors] [DOTProducts] [DIPOles] [STATistics]
                                   [atom-selection]

   PROJect      mode-definition [mode-spec] [magnitude-spec]

   DIAGonalize  [NFREquencies integer] [NADD integer] [RAISe] [FINIte [STEP real] 
                [DSCF]] [ENTRopy [TEMPerature real] [SIGMa real]]

   QUASIharmonics [NUNIts integer]  FIRStunit integer [NSKIp integer]
                    atom-selection [FCUToff real]
                    [NFREquencies integer] [NADD integer] TEMP real [THERmo [RESI]

   REDUce         IUNBas int [IUNTrans int] [BIG] [NFREq integer] [RAISe]
                  [FIX] [CMPAct]

   RBQUasi        IUNBas int [IUNTrans int] [NFREq integer] TEMP real
                    [NUNIts integer]  FIRStunit integer [NSKIp integer]

   DIMB           IUNMode int [IUNRead int] [PARDim int] [NFREq int]
                  [NADD int] [RAISe] [CUTF1 real] [ITERation int]
                  [SAVFreq int] [TOLErance real] [STRENgth real] 
                  [DWINdow] [BIG] [SCILib]

   EXPLore        [mode-spec] [magnitude-spec]  unit-spec
                    [GRID int] [COMP] [SHAKe]
                      [ [ADJUst] [BOLTzmann temp] [GAUSs factor] ]

   FLUCtuations   { ATOM [atom-selection] } [mode-spec] [magnitude-spec]
                  { IC                     }       [QUANtum] [VERBose]
                  { USER [atom-selection] }

   PAFLuctuations { ATOM         } [atom-selection] [mode-spec] [magnitude-spec]
                  { GROUP [MASS] }   [QUANtum] [VERBose] [COORdinate] [SAVE]
                  { USER         }   [CONTinue]

   RAYLeigh       [mode-spec] [SAVE]

   PED            [mode-spec] [magnitude-spec] [TOL real]

   THERmodynamic  [mode-spec] [TEMPerature real] [STEP real] [FCUToff real]

   EDIT    { INCL  mode-definition  [ORTHog]  [TO mode]    } [atom-selection]
           { REMOve [mode-spec] mode-definition [NONOrm]   }
           { DELEte [mode-spec]                            }
           { ORTHogonalize [PURGe] [mode-spec] [TOL real]  }
           { SHAKe mode-spec                               }
           { ZERO mode-spec                                }
           { MOVE MODE n [TO m] [SCALe real]  [NONOrm]     }
           { ADD  MODE n [TO m] [SCALe real]  [NONOrm]     }
           { MULT MODE n SCALe real                        }
           { SET  MODE n SCALe real        [NONOrm]        }
           { COPY MODE n TO m                              }

   BASIS   { IC     { FIRSt BOND      }    }  [NOORthonorm]
           {        { FIRSt ANGLe     }    }
           {        { DIHEdral        }    }
           {        { SECOnd ANGLe    }    }
           {        { SECOnd BOND     }    }
           {                               }
           { TR atom-selection [BALAnce]   }

   FILL    {  DIFF } [mode-spec] [magnitude-spec] [APPE]
           {  COMP }

   CORREL  [ MAXTimesteps int ]  [ MAXSeries int ]  [ MAXAtoms ] [ COVAriance]


   MASS    ! turn on  mass weighting flag (default)
   NOMAS   ! turn off mass weighting flag

   SUBSystem   [atom-selection]

   END

   unit-spec ::= UNIT unit-number

   mode-spec ::= MODE integer [ THRU integer ]

   mode-definition ::=
           { { TRAN } { X }                } [NONOrm] [NOTR]
           { { ROTA } { Y }                }
           {          { Z }                }
           {                               }
           { SPHEre   { X  } [IX int]      }
           {           { Y  }  [IY int]    }
           {           { Z  }   [IZ int]   }
           {           { R  }    [IR int]  }
           {           { TX }              }
           {           { TY }              }
           {           { TZ }              }
           {                               }
           { COMP                          }
           { DIFF                          }
           { FORC                          }
           { USER integer                  }
           { BOND atom atom                }
           { ANGL atom atom atom           }
           { DIHE atom atom atom atom      }

   atom::= {residue-number atom-name}
           { segid  resid atom-name }
           { BYNUm  atom-number     }

   magnitude-spec ::= { TEMP real TFRE real }  [NONOrm]
                      { KCAL real TFRE real }
                      { RMS  real           }
                      { MRMS real           }
                      { FACT real           }

   miscellaneous-command-spec  *note misc:(chmdoc/miscom.doc).
   ic-subcommand-spec          *note ic:(chmdoc/intcor.doc).
   coor-subcommand-spec        *note coor:(chmdoc/corman.doc).


.. _vibran_normal_modes:

Normal Modes
------------

There are two ways that normal modes are stored internally
in CHARMM. The most common usage is as one double precision mass
weighted array. A series of such arrays usually span an orthonormal
basis (as would be the case upon diagonalization). The second method
is to represent a normal mode as three non mass weighted coordinate
displacement arrays, stored in double precision. The program automatically
converts between them as necessary. Whenever interconversion is to be
done, a "magnitude specification" may be given. This specification
requires a step type and step length. The valid step types are;

*        FACT (simple factor)
*        TEMP (put mode at this temperature)
*        KCAL (put in the specified Kcals)
*        MRMS (step along until this mass weighted RMS is obtained)
*        RMS  (step along until this RMS is reached).

When specifying TEMP or KCAL, a terminal frequency (cm-1) may be specified
(default TFRE is 5.0) which prevents excessive stepping along very low
frequency or translation-rotation modes.

The procedure used in going from double precion to coordinate
displacement arrays is;

1) Normalize double precision vector (unless NONOrm keyword is used)
2) Mass weight by dividing by root(mass)
3) Multiply by appropriate scale factor from step type and length

To convert from single precision into double precision, the
procedure is;

1) Save inner product (as initial step length)
2) Mass weight by multiplying by root(mass)
3) Normalize vector (unless NONO keyword is used)
4) Compute scale factor using initial step length and step type

Whenever a magnitude-specification is called for, some
interconversion will take place. The conversion from double precision
to coordinate displacements takes place in the subcommands;

::

                PRINt NORMal-modes
                WRITe TRAJectory
                PROJect
                FILL
                EXPLore
                FLUCtuations
                
The interconversion from coordinated displacements to double precision
takes place in the subcommands;

::

                PROJect
                EDIT NORMal-modes

The Normal Mode data structure of double precision arrays is local
to the Vibrational analysis section of CHARMM, and the storage space for
these arrays is released when exiting to CHARMM via the END command.


.. _vibran_io:

I/O For Normal modes
--------------------

The VIBRAN section supports its own I/O commands. Commands to
read, write and print coordinates (see :doc:`corman`)
and internal coordinates (see :doc:`intcor`) are
identical with those in the main program. This section can read and
write the Normal Mode data structure and write out the second
derivative matrix in several ways (for external use) and normal mode
trajectories (for the movie programs). Also, useful information about
normal modes may be printed using the PRINT NORM command.

::

   READ    { NORMal-modes  [CARD] unit-spec  [APPEnd] [mode-spec]  }

   WRITe   { NORMal-modes  [CARD] [mode-spec]                      } unit-spec
           { SECOnd-derivatives  [RAISe] [MASS] -
                   [CARD [FINIt [STEP val] [TOL val] [atom-selection] ]] }
           { TRAJectory  mode-spec magnitude-spec [PHAS real] [SHAKe] -
                   [SEQUential-files] [NCYC integer] [STEP real]
                   [SUPErpose] [RANDom] [ISEEd integer] }

   PRINt   NORMal-modes  [mode-spec] [magnitude-spec] [INTDer [FINIt]] -
                   [VECTors] [DOTProducts] [DIPOles] [STATistics] -
                                   [atom-selection]

By default, normal mode vectors are read and written in binary (FILE)
format, but ascii (CARD) can be specified. When writing, a unit must be 
specified, and a contiguous subset of modes may be specified using the 
mode-spec. When reading modes, all selected modes (from mode-spec) of the 
modes in the normal mode file will be read (assuming there is enough space). 
Note, when modes are selected on input, the selection is relative to the 
ordering of modes in the input file and does NOT correspond to the destination 
of these modes.

If the available space is exhausted, a warning is issued, and further
reading stops. Existing modes will be deleted when the READ NORM command
is executed unless the append option is used, in whichcase, the new
modes are added sequentially at the end. No modification of modes is
done upon reading (i.e. normalization, or orthogonalization).

When printing Normal Modes, a variety of options may be specified.
A contiguous subset of modes may be specified (the default is all modes), and
an appropriate magnitude may be specified (see *note modes:Normal Modes.).
For each specified mode, the frequency (cm**-1), eigenvalue (Kcal/gram/A**2),
force projection (Kcal/mole/A), % translation-rotation, and magnitude
information will be printed.  In addition, internal derivatives (or optionally
by finite differences with the FINIte keyword) of the internal coordinate
data structure (INTD keyword), displacements in coordinate space (VECTors
keyword), dipole derivatives (DIPOles keyword), dotproducts with other modes
can be printed (DOTProducts keyword), and some vector statistics (STATistics
keywords).

The Second Derivative matrix may be written out in binary
or card format for input to other programs. In particular, there is a
program that extracts normal modes from a large secular equation (DIAGIT).
There is also a restartable diagonalization program for medium sized systems.
The RAISe keyword will shift the translation-rotation modes to high
frequency (currently 5000 cm-1). This option only work with the card output
format.  The MASS keyword will calculate and write out M**(-1/2) * H * M**(-1/2)
which is the matrix used in the eigensolver with the DIAG command.
Finally, this command may be used to test the second derivative
energy surface determination by a comparison of finite derivatives and
analytic derivatives. This is done with the FINIte and CARD keywords.
Two values are also looked for, the finite difference step size (STEP,
def=0.005), and the difference tolerence for printing (TOL, def=0.0001).

Trajectory files may be written out for a set of modes with
a given magnitude factor. The modes for all specified modes may be
written out in one file, or in separate files for different modes
where sequential unit numbers are used starting with the specified
unit number. The SEQU keyword will cause sequential files to be written.
By default, one cycle of 12 frames will be written. The keywords
NCYC and PHASe can be specified to alter the number of cycles and the
phase angle between frames within a cycle. The PHASe keyword should be
an integer divisor of 360.0 ( a PHASe value of 7.2 will result in 50
frames per cycle). The total number of frames is given by NCYC * 360/PHASe.
If it is desired to specify a particular time step between frames,
the STEP keyword will specify this value in picoseconds. The total
number of frames will remain unchanged will remain unchanged.
If a movie of 10 modes and displaying one picosecond of each mode is
desired with 500 frames per picosecond (about 25 seconds of film time),
then the appropriate input could be:

::

        WRIT TRAJ MODE 7 THRU 16 STEP 0.002 PHAS 3.6 NCYC 5 TEMP 2000.0

It is also possible to write out a trajectory that consists of a
superposition of a number of normal modes. The STEP parameter needs to be
specified, together with PHAS and NCYC, to set the length of the trajectory
in ps. By default, all normal modes in the trajectory have an initial phase
of 0. Random initial phases can be specified by using the keyword RANDom
and specifying a seed for the random number generator (default ISEED =
314159). A movie of 50 ps (10 frames per ps) for modes 7 to 100 with random
initial phases could be generated by:

::

        WRIT TRAJ MODE 7 THRU 100 SUPE STEP 0.1 PHAS 3.6 NCYC 5 RAND TEMP 300.0

.. _vibran_diagonalization:

Diagonalization of the second derivative matrix
-----------------------------------------------

The DIAGonalize command will generate the second derivative
matrix, mass weight this matrix, and then diagonalize to generate
normal modes. There are some needs and restrictions for this command.
They can be summarized;

1) SHAKE may not be used
2) ST2 waters may not be used
3) Periodic boundaries and images may not be used without FINIte (see crystl.doc)
4) Group electrostatics and extended electrostatics may not be used without FINIte
5) All coordinates and energy lists (nonbond, hbond) must be set up
6) The number of atoms does not exceed the limit (300)
   (unless BOMLEV is reduced)
 
Once normal modes have been generated, they may be saved (WRITE NORM...),
or analyzed directly.

::

   DIAGonalize  [NFREquencies integer] [NADD integer] [RAISe] [FINIte [STEP real] 
                [DSCF]] [ENTRopy [TEMPerature real] [SIGMa real]]

The RAISe keyword, will cause the normal modes corresponding to
translation-rotation to have a very high frequency (currently 5000 cm-1).
This option is intended for calculations where rotational coupling terms
are to be removed.

The NADD option will cause the the specified number of lowest
modes to be skipped in the evaluation. For example, if "NADD 100"
is specified, then the first mode of the result will correspond to mode
101 in the actual matrix. This option has been added, so that modes
of moderatly sized systems (200-400 atoms) can be found in groups
when the memory requrement for a full calculation are prohibitive.

The FINIte keywords causes the Hessian to be generated from the finite
differences of the forces.  This option requires 6*N+2 energy determinations,
so it should be reserved for smaller systems.  The step length for finite
difference may be specified with the STEP keyword (default 0.005 Angstroms).
A STEP value that is too large will cause errors due to anharmonicity.
A STEP value that is too small will result in inaccurate hessian elements
(force differences are divided by the step length).

Regular vibrational analysis treats Drude particles as real atoms and
thereby shows extra peaks on the IR-spectrum due to Drude particles. This
complicates CHARMM IR-spectrum comparison with QM or experimental data.
The purpose of DSCF keyword is to allow vibrational spectrum analysis
for molecules in presence of Drude particles to be conducted in such
way that Drude particles will not be explicitly present in the IR-spectrum. 
This regimen is called SCF Drudes where position of Drude particles is
instantenously adjusted to any rearrangement in position of real atoms.
In this mode CHARMM performs calculation of second derivatives by using
finite differences applied to real atoms and followed by Drude coordinate
relaxation after every change in coordinates of real atoms. This way
Drude degrees of freedom are projected onto second derivatives of real
atoms. This algorithm works via numerical differentiation only since 
analytical solution is hardly possible. DSCF keyword can be invoked when
DIAG FINIte keywords are specified, otherwise it is ignored.


Entropy calculation
-------------------

::

   DIAGonalize  [NFREquencies integer] [NADD integer] [RAISe] [FINIte [STEP real] 
                [DSCF]] [ENTRopy [TEMPerature real] [SIGMa real]]

After second derivatives are calculated the vibrational entropy term
can be evaluated. Two other entropy terms, rotational and translational ones
are also calculated (see corman.doc for details). Entropy calculation is
implemented as an extention to VIBRAN DIAG command, since entropy calculation
is based on vibrational frequencies. 

Default value for TEMPerature is 298.15 K. Default SIGMa value is 1.0.
SIGMa is symmetry number which is 1 for non-symmetric molecule and some
low symmetry groups. For symmetric molecules one should enter a correct
value for sigma (see, for example, C.J.Cramer, "Essentials of Comp.Chem.",
2002,p.327).

EXAMPLE:

::

   VIBRAN
   DIAG ENTRopy
   END

The units for entropy are cal/(mol*K). Rotational, translational and
vibrational entropy terms along with their sum can be accessed in
CHARMM input file as ?SROT, ?STRA, ?SVIB, and ?SSUM substitution
parameters after executing the entropy command.

Alternative implementation of entropy calculation is available under 
keyword THERmo (see above). This is a kind of code duplication that should be 
resolved sometime. The THERMo also includes questionable contribution from 
trivial modes and calculates a vibrational entropy term only.


.. _vibran_quasiharmonics:

Quasiharmonic dynamics from molecular dynamics
----------------------------------------------

::

   QUASIharmonics [NUNIts integer]  FIRStunit integer [NSKIp integer]
                    atom-selection [FCUToff real]
                   [NFREquencies integer] [NADD integer] TEMP real [THERmo [RESI]]

For quasiharmonic dynamics, a dynamics trajectory file(s) is
read. During this read, atom position fluctuation tensors are generated for
the selected atoms (default: all).

The reference coordinates must be present in the main coordinate arrays.
The COOR DYNA command can preceed the QUASi command if the average
coordinates are to be used as a reference.  Another good choice is to
use an energy minimized coordinates set as the reference.  The resulting
fluctuation matrix is mass weighted and diagonalized. The resulting
modes are the quasiharmonic modes of the system, and should roughly
match those from a strightforward normal mode calculation. The significant
differences of these methods are that anharmonic terms are included here.
This method can yield misleading results when there are long lived
transitions, or when there is very slow interchange of energy between
degrees of freedom. The estimated frequency for a given mode will
strongly depend on the average energy in this mode. For the case of
slow energy transfer, this may be significantly different from kT.
If internal motion is to be studied, it is strongly recommended that
translation and rotation be removed from the dynamics trajectory.
This is done with the MERGe command using all atoms and
mass weighting. Without mass weighting, the translation-rotation motion
is not projected out.

In order for this method to estimate frequencies from fluctuations,
the average temperature (TEMP real) must be specified.

Keyword THERmo evaluates entropy [kcal/mol/K], enthalpy [kcal/mol]
and heat-capacity Cv [kcal/mol/K] at the specified temperature. RESI prints
out S, H and Cv using the fluctuations of the selected atoms in each residue.
Only frequencies > FCUToff are used in the sum (default value for FCUT=0.0001)
See also COOR COVA (corman.doc).

The dynamics trajectory is read in the usual manner.
Specifications are;

::

        NUNIts integer  - number of I/O units
        FIRStu integer  - first I/O unit
        NSKIP  integer  - integration step modulo to use

As with the DIAG command, this command will also accept the keywords;
NFREquencies integer, and NADD integer.

Reduced basis quasiharmonic calculations are also possible (RBQUas).
See the section for reduced basis normal mode analysis for information
on setting up a basis set.  The method is exactly the same as above
except the fluctuation matrix is generated in the reduced basis instead
of the mass weighted cartesian basis. 

::

   RBQUasi        IUNBas int [IUNTrans int] [NFREq integer] TEMP real
                    [NUNIts integer]  FIRStunit integer [NSKIp integer]

.. _vibran_reduce:

Reduced basis normal mode analysis
----------------------------------

With the REDUce command, it is possible to do reduced basis
normal mode analysis. This can be affected to constrain certain degrees of
freedom (leave them out of the basis), or to reduce the size of a large
problem. I can also be used to remove translation and rotation from a
calculation. The same restrictions for a full normal mode calculation
apply here as well. 

By defining the keyword CMPAct the initial second derivative matrix
in Cartesian coordinates will be stored in a compressed form. This
means that only non-zero elements are stored. This option is very
useful if you deal with very large molecules, where the initial
second derivative matrix may be too large to fit into memory. For
instance, in a system with 2000 atoms, the regular Cartesian Hessian
takes up about 144 Mbytes, whereas the compressed matrix takes up
an absolute maximum of 43 Mbytes (usually less, depending on nonbonded
cutoff distances). Currently, this option allows atoms to interact with
up to 300 other atoms, on average. This means that if large cutoff
distances (>15 A) are used, the command may exit with an error message.

This command requires an external basis set. The form of the basis
set is a CHARMM normal mode file of the correct dimension. The basis will
include all of the vectors in that file. The basis is expected to span an
orthonormal space. These basis sets may be made from the BASIS command
followed by the EDIT ORTHog command if appropriate.

The reduced basis may also be a subset of normal modes. For large
systems, this is the refinement step. By repeating the calculation in a
reduced basis, most roundoff effects from the tridiagonalization are removed.

This command removes all existing vectors and replaces them with
the result of the diagonalization in the original basis (i.e. the new vectors
are linear combinations of the basis vectors). The number of vectors to
compute may be specified with the NFREquencies keyword. The default is to
compute all frequencies (for the entire reduced basis).

::

   REDUce IUNBas int [IUNTrans int] [BIG] [NFREquencies integer] [RAISe]
          [FIX] [CMPAct]

The IUNBas keyword points to the basis vector file. If the
IUNTrans value is specified, the transformation matrix (NDIMxNFREQ)
relating the basis vectors and the final normal modes will be saved to
that unit. The number of frequencies to be computed may be specified.
The default is to compute all of them. The raise option raises the
translation/rotation degrees of freedom to a high frequency. This keyword
will have no affect if the net translation/rotation degrees of freedom
are not in the basis.

The BIG keyword signifies that the entire set of basis vectors
is not to be stored in memory. This is essential for large calculations
(i.e. when HEAP allocation errors occur without this option). In this
case, VIBRAN ... NMODE 1 ... should be specified to save space, since
no space is needed to store the basis vectors or the results.
The BIG option will cause the backtransformation step to be suppressed,
thus the IUNTrans unit should be specified (or all you will get is the
frequencies). The will be a special purpose support program to backtransform
this eigenvactors (normal modes) to the original basis.

If the BIG option is not specified, there must be enough
space allocated to store the entire basis in memory (VIBRAN NMODe int).

The keyword FIX should be used when part of the system is fixed by the
CONStraint FIX command. In the modes that are calculated, the fixed atoms
will not move, but their presence is taken into account for the motions 
of the unfixed atoms. This option does not require a basis vector file,
i.e., the IUNBas parameter should not be specified. For large systems,
the use of a compressed second derivative matrix by specifying the
keyword CMPAct is recommended. The BIG option is not allowed.


.. _vibran_dimb:

Iterative Mixed-Basis Diagonalization 
-------------------------------------

With the DIMB command, one can do an iterative diagonalization,
which will consume a lot less memory space than a regular (DIAG) 
diagonalization, but will still result in exactly the same normal modes.
The method (see L.Mouawad and D.Perahia (1993), Biopolymers, 33, 599)
does repetitive reduced-basis diagonalizations. These reduced
bases are constructed partially from the not yet converged eigenvectors
and from regular Cartesian coordinates. A unit IUNMode for writing out
the eigenvectors has to be specified. For restarts from existing
eigenvector files, the unit IUNRead should be specified. The NFREq,
NADD, and RAISe keywords have the same effect as in the regular DIAG
command. The maximum size (in atoms) for a block to be diagonalized
is defined with the PARDim (default=200) keyword, and can be tailored to
the size of the available memory. This affects the 
number of blocks the system will be divided into. Make sure that the
requested number of modes (lesser of NMODEs and NFREq) is smaller
than PARDim*3. 

If no IUNRead unit is specified, initial basis vectors are calculated
by diagonalizing main diagonal blocks of the full matrix. This initial
calculation involves the use of residual (Lanczos) vectors, and will 
result in only half the number of requested modes. The full number of 
requested modes can be obtained by using the following sequence of
commands (in this example NMODes is set to 200 and PARDim to 100, but
these values can be changed):

Create the initial basis and write out on unit 20. This results
in (NMODes+6)/2 = 103 basis vectors. The usage of BIG reduces
memory requirements by writing intermediate vectors to disk
temporarily. Specification of ITERations 0 exits the routine
after the calculation of the basis vectors:

::

   OPEN WRITe FILE UNIT 20 NAME initial.bas
   VIBRan NMODes 200
   DIMB ITERations 0 PARDim 100 IUNMode 20 BIG
   END

Do the iterative diagonalization, using the precalculated initial
basis of 103 vectors. 103 converged modes will be written to unit 10:

::

   OPEN READ  FILE UNIT 20 NAME initial.bas
   OPEN WRITe FILE UNIT 10 NAME modes.mod
   VIBRan NMODes 103
   DIMB ITERations 100 TOLErance 0.05 PARDim 100 -
        IUNMode 10 IUNRead 20  DWIN
   END

The usage of the keyword BIG is only allowed when ITERations 0 is specified.
It will use unit IUNMode to temporarily store intermediate
vectors. Make sure that there is enough available disk space. The size
of the intermediate file will be about (NMODes*NAT3*8) bytes, 
where NAT3 = 3 * number of atoms.
A rough estimate of the required memory space is:

::

   { 4*(PARDim3*PARDim3) + (2*NMODes+10)*NAT3 } * 8 bytes  (without BIG)
   { 4*(PARDim3*PARDim3) + (1*NMODes+10)*NAT3 } * 8 bytes  (with BIG)

Usually convergence will be faster with larger PARDim and with
smaller number of modes. Too small a number of modes will slow
down converge due to a limited basis set of eigenvectors.
With CUTF1 (default=50.0) a cutoff frequency (cm-1) can be defined. 
Only normal modes with frequencies below CUTF1 will be calculated. 
The maximum number of iterations can be specified with the ITERations 
(default=10) parameter, and the desired accuracy of the eigenvectors 
with TOLErance (default=0.05). TOLErance values range between 0 and 1, 
where 0 is the highest accuracy.

The standard DIMB algorithm uses a single "window" of Cartesian
coordinates that is added to the basis set. When the DWINdow keyword
is specified, two different "windows" are added. The DWINdow method
is more efficient, because it uses windows with strongly interacting
atoms more often than those with weakly interacting atoms. When the
DWINdow option is chosen, the parameter STREngth (default=0.0) defines 
which atom sets are considered to be interacting at all. If the sum of
absolute values of the second derivatives between two particular
atom sets is lower than STREngth, those second derivatives will not be 
considered.

SAVFreq (default=number of blocks) defines the frequency of saving 
the eigenvectors to disk for the DWINdow option.
When DIMB calculations are carried out on a Cray vector computer,
it is advantageous to use the general EISPACK routines do
diagonalize the submatrices. This can be defined by adding the
keyword SCILib.


.. _vibran_explore:

Exploring Energy Surfaces
-------------------------

The explore command is used to search along a particular mode and
to compute the energy at regular intervals. For now, this search is limited
to one dimension, but with the EXPL COMP (using the comparison coordinates)
option, a two dimensional search is possible. This is done by filling the
comparison coordinate set with a structure perturbed along one mode, and then
exploring along another.

::

   EXPLore         [mode-spec] [magnitude-spec]  unit-spec -
                    [GRID int] [COMP] [SHAKe] -
                     [ [ADJUst] [BOLTzmann temp] [GAUSs factor] ]

For this command, a magnitude specification and grid selection
are read. The grid selection determines the number of energy evaluations.
An odd number will include the center point, and an even number will not.
The default value is 3, one energy evaluation at each extreem (determined
by the magnitude specification), and one at the center.

When the energies are computed, they are stored in a temporary
array which may be output (for plotting) by specifying a unit number.
They may also be used to adjust the eigenvalue of the particular mode
base on a least squares quardratic fit. This is done with the ADJUst
keyword. The weighting for this fitting is currently unity for every point.
The adjust feature is intended for relatively small displacements. Large
"rigid" dislacements will necessarily come up against bad close contacts.
Quadratic fittings for this type of interactions, will be misleading.
For this reason, two weighting options are available, an energy
Boltzman weighting (BOLTzmann keyword-value), and a gaussian displacment
weighting (GAUSs keyword-value) which weights nearby points with a
larger factor. When these keyword are used together, the overall weighting
is a product of individual weightings.

To avoid problems with excessive (quartic) bond stretching which
occurs with large displacements, SHAKE may be invoked with the SHAKe
keyword. In order to use this option, the SHAKe command in CHARMM must
fist be invoked specifying which bonds (and angles) are to be maintained.
A step value may be specified. By default, all shake comparisons
are done with respect to the center coordinates. This sometimes leads
to "DEVIATION IN SHAKE TOO LARGE" errors. The step option allows shake
to be invoked in small steps. The step value gives the relative step
for each intemediate SHAKE calculation. A value of 1.0 is the default, and
a value of 0.2 will cause the extreem points to be computed in 5 steps.


.. _vibran_fluctuations:

Computing fluctuations from normal modes
----------------------------------------

The fluctuation command computes atom, internal coordinate, or
user specified fluctuations from normal modes. For atom fluctuations,
the magnitude of the components for the overall fluctuation is stored in
the comparison coordinate arrays. In addition, information giving the
contribution of each mode is printed. For the IC (internal coordinate) option,
the overall fluctuations are stored in the IC table. the "IC WRITE..."
command may be used to save this data.

::

   FLUCtuations { ATOM [atom-selection] } [mode-spec] [magnitude-spec]
                { IC                    }       [QUANtum] [VERBose]
                { USER [atom-selection] }

For the user specified option, the user may supply the subroutine
USEFLU which is call once for initialization (IMODE=0), once for each mode
(IMODE=n), and once for termination (IMODE=-1). In preparing this routine,
see the existing routine for interfacing requirements. If no USEFLU routine
is provided, the equal time cross-correlation matrix for all selected atoms
will be computed.

For this command, a selection of modes may be made (default all),
a selection of atoms may be made (default all), and a magnitude
specification may be made. For this application, a temperature factor
is usually used.

The QUANtum keyword may be used if a quantum scaling factor
is desired. For this option, higher frequency modes would have a diminished
amplitude for a given temperature.

The VERBose option causes the contribution for each selected
atom with each selected mode to be printed. This is useful when a power
spectrum is to be computed or prepared for plotting.

.. _vibran_princ_axis_flucts:

Calculation of anisotropic (and isotropic) fluctuations from Normal Modes
-------------------------------------------------------------------------

The PAFL command sums over the current set of normal modes (default)
or some subset of this (given by the [mode-spec]) to create a fluctuation
matrix. The three-by-three fluctuation matrix for each atom may be
diagonalized (default) to produce the principal axis flucutions
given in terms of three mutually perpendicular unit principle axes.  In
general, there will be a different set of principle axes for each atomic
center.  Alternatively, the COORdinate option may be used, in which case
the fluctuation matricies are not diagonalized.  Instead, the fluctuations
are given in term of the coordinate axes in which the normal modes were
calculated.  Since the size of the fluctuations are temperature-dependent,
a [magnitude-spec] in terms of temperature (such as TEMP 300.0) should
be given.

::

   PAFLuctuations { ATOM         } [atom-selection] [mode-spec] [magnitude-spec]
                  { GROUP [MASS] }   [QUANtum] [VERBose] [COORdinate] [SAVE]
                  { USER         }    [CONTinue]

The summation over modes which creates the fluctuation matrix may
be controlled in a number of ways.  The easiest of these is to
specify a TFREquency.  Normal modes of frequency
lower than TFREquency are simply ignored in the summation.  By setting
TFREquency above the rotation-translation modes, but below the lowest
vibrational mode, one achieves what is the desired result (the fluctuations
due to internal modes of vibration) for most applications.  For applications
requiring more specificty (such as calculating the fluctuations due to
modes 8,9,10,13,14, and 15 only) one may use a more detailed [mode-spec]
to select the desired modes over which the summation is to be performed).

The second way of controling the summation is necessary when
the normal modes for the desired system can not all be read into CHARMM
at one time due to space allocation limitations.  In this case one uses
the SAVE and CONTinue options.  SAVE tells the routine that you don't
want to manipulate the fluctuation matrix once it is created.  CONTinue
indicates that you don't want to reset (zero) the fluctuation matrix
before the present summation.  Therefore, if the normal modes are
stored on, let us say, four files, one reads in the first normal
mode file and gives the PAFL SAVE command. One then reads in the
second normal mode file (without the APPEnd option) and gives the
PAFL SAVE CONTinue command.  One then reads in the third normal
mode file (without APPEnd) and gives the PAFL SAVE CONTinue command.
Finally, one reads in the fourth normal mode file (without APPEnd)
and gives the PAFL CONTinue command.  Each of these four PAFL commands
should have the same options (such as GROUP, ATOM, USER; same TEMP;
same atom selection, MASS or not).  Moreover, there should be no overlap
between modes in each of the four sets.  If the four (or any n) files
have overlapping modes, then a selection must be made so that
no mode is counted twice in the summation (see the test case).

For the user specified option, the user may supply the subroutine
USPAFL which is call once for initialization (IMODE=0), once for each mode
(IMODE=n), and once for termination (IMODE=-1). In preparing this routine,
see the existing routine (and also USEFLU) for interfacing requirements.

The PAFL command is quite similar to the FLUCtuations command.
VERBose prints out lots of stuff and QUANtum uses quantum scaling in
calculating the contribution of each mode to the fluctuations.

When ATOM is specified, the fluctuations are calculated individually
for each atom currently selected.  The fluctuations are given in angstroms,
and so mass-weighting would make no sense.  When GROUp is specified, the
fluctuations for the movement of all specified atoms moving as a group are
calculated.  When MASS is also specified, the fluctuation of the
center of gravity of that group of atoms is calculated.  Otherwise, it is
the fluctuation of the center of geometry which is calculated.  At present
USPAFL does nothing, and one would need to read the code to be able
to interface this routine properly with the surrounding code.


.. _vibran_projections:

Projection of normal modes onto vectors of interest
---------------------------------------------------

The PROJect command will project the selected modes onto one
of the definable vectors.

::

   PROJect          mode-definition [mode-spec] [magnitude-spec]

The information displayed is;

::

        MODE integer    - mode number
        FREQUENCY       - frequency of this mode (cm-1)
        NORMAL DOTPR    - actual dotproduct of normalized vectors
        PROJECTION      - ratio of guess vector projection to step length
        APPROX DEL E    - estimate of energy increase along this mode
        TYPE            - type of step (FACT, RMS, KCAL, TEMP)
        STEP            - step length for this step type
        
The NONOrm keyword may be used to prevent normalization of modes before
projecting.


.. _vibran_mode_definition:

Mode Definition
---------------

Several commands use a mode_definition, which allows the
specification of a vector. There is a wide variety of possible vectors
which may be specified, and there is also a user supplied routine, USERNM,
which may be used to get any other desired motion (vector). These vectors
may be used for the analysis of normal modes, or as additions to the
basis of vectors (EDIT INCLude) for further analysis.

::

   mode-definition ::=
           { { TRAN } { X }                } [NONOrm] [NOTR]
           { { ROTA } { Y }                }
           {          { Z }                }
           {                               }
           { SPHEre   { X  } [IX int]      }
           {           { Y  }  [IY int]    }
           {           { Z  }   [IZ int]   }
           {           { R  }    [IR int]  }
           {           { TX }              }
           {           { TY }              }
           {           { TZ }              }
           {                               }
           { COMP                          }
           { DIFF                          }
           { FORC                          }
           { USER integer                  }
           { BOND atom atom                }
           { ANGL atom atom atom           }
           { DIHE atom atom atom atom      }

   atom::= residue-number atom-name


The NONOrm keyword supresses the automatic normalization of
the specified vector. This is desired for some applications, but is
not normally needed.

The NOTR keyword removes any new translation/rotation from the
specified mode. This may be needed when only internal motions are to
be analysed.

The vectors which may be defined are;

::

   TRAN  X           - Translation along the X-axis
   TRAN  Y           - Translation along the Y-axis
   TRAN  Z           - Translation along the Z-axis
   ROTA  X           - Rotation along the X-axis
   ROTA  Y           - Rotation along the Y-axis
   ROTA  Z           - Rotation along the Z-axis
   COMP              - Use the comparison coordinates (as a vector)
   DIFF              - Use the difference between the MAIN and COMP coordinates
   FORC              - Use the forces from the last energy evaluation
   USER integer      - Use a user specified vector (defined in USERNM)
   BOND atom atom    - Use vector which stretches the specified bond
   ANGL 3X(atom)     - Use vector which bends the specified angle
   DIHE 4X(atom)     - Use vector which twists the specified dihedral
   SPHEre ...        - Appropriate homogeneous motion (spherical harmonics)

For the TRANslate, ROTAte, COMP, DIFF, and FORCe options, an
atom selection may be specified. Any nonspecified atoms will have a zero
values in the vector.

For the BOND, ANGLe, and DIHEdral options, the atoms specified
don't have to be bonded together or have any special connectivity, but
the first atom specified must not close back through a loop to the last
atom specified. If this is a problem, then one of the bonds in the loop
must be deleted for this calculation to work properly.

For angle terms, the first two atoms specified (and anything
they are connected to) move as one unit, and the third atom and anything
its connected to will move as a second unit.

For dihedrals, the first 2 atoms specified define one block, and
the last 2 atoms defines the second block. The axis of rotation is about
the middle 2 atoms. For example if one specifies

::

        PROJECT DIHE 1 N 1 CA 1 CB 1 CG1
        
for an isoleucine residue, the atom 1 CG2 will also rotate with CG1 because
it is connected to the third atom (CB) which is part of the second unit.

.. _vibran_rayleigh:

Rayleigh Quotients
------------------

The RAYLeigh command will compute the second derivative matrix,
and project all selected modes onto the second derivative matrix after
mass weighting. If the modes are eigenvectors, then the resulting quotients
will correspond to the eigenvalues. The quotient values are given by;

::

                   -1/2    -1/2
        Q = < v | M    H  M    | v >


::

   RAYLeigh        [mode-spec] [SAVE]

If the SAVE keyword is given, then the quotients and corresponding
frequencies are saved and become part of the data for the selected modes.
For example, a subsequent PRINT NORM command would then use these new values
in place of the original frequencies.

.. _vibran_ped:

Potential Energy Distribution
-----------------------------

This command is designed to work with one mode at a time.
For a given magnitude specification, it computes the expectation
value for the energy contribution change for each internal coordinate
term (bond, angle, dihedral, and improper dihedral) and prints that
term if the fluctuation is greated than the tolerence (default TOL 0.0001).

::

   PED        [mode-spec] [magnitude-spec] [TOL real]

For this method, it assumes that the structure is at a stationary
point and that the potential is quadratic. In addition to printing out
the individual energy terms, it also prints out the total for each class
(bonds, angles...). These values can be used to determine the type of any
given mode.

.. _vibran_edit:

Editing the set of normal modes
-------------------------------

There are several commands that can modify the normal mode
vector space. In addition to the obvious ones such as READ NORMal-modes
and DIAGonalize,  there is also an EDIT command which can be used to
modify the normal mode data structure.

::

   EDIT    { INCL  mode-definition  [ORTHog]  [TO mode]    } [atom-selection]
           { REMOve [mode-spec] mode-definition [NONOrm]   }
           { DELEte [mode-spec]                            }
           { ORTHogonalize [PURGe] [mode-spec] [TOL real]  }
           { SHAKe mode-spec                               }
           { ZERO mode-spec                                }
           { MOVE MODE n [TO m] [SCALe real]  [NONOrm]     }
           { ADD  MODE n [TO m] [SCALe real]  [NONOrm]     }
           { MULT MODE n SCALe real                        }
           { SET  MODE n SCALe real        [NONOrm]        }

The EDIT DELETE command will delete specified modes from the
date structure.

The EDIT ORTHogonalize command is used to orthogononalize
and optionally normalize a particular subset of normal modes.

Additional modes can be added  with the EDIT INCLude command.
This command will add on the defined mode to the end of the normal mode
vector, unless the destination is specified with TO n in which case it will
be placed as mode number n.

The translation and rotation options are important when setting
up guess vectors for the iterative diagonalization program, or when
individual coreolis coupling terms are needed. Single additional modes may
be added with the EDIT INCLude DIFF or the EDIT INCLude FORCe or the
EDIT INCLude USER or the EDIT INCL COMP commands. For the DIFF option,
the difference between the comparison and main coordinate sets will
be appended, for the FORCe option, the current values in the force
arrays (from the last energy evaluation) will be appended. The COMP
option will use the comparison coordinates as the appended mode. This
is option is intended for use in the case that the comparison coordinate
set is filled with displacements. There is also the ability to specify
a user vector. This is done by the inclusion of the subroutine USERNM
in your USERSB and using the USERLINK facility. For all of the EDIT
INCLude options, the defaults are to root mass weight, normalize, and
orthogonalize to the rest of the vector space. To skip any of these steps,
the NONO, and NOOR keywords must be specified.

See (:ref:`vibran_projections`) for a complete list of
definable modes.

The following commands (which all can take an optional atom-selection)
provide a simple means of manipulating, and reshuffling modes.

::

      EDIT ZERO sets the specified modes to zero.
      EDIT MOVE MODE n [TO m] [SCALe real]       mode(m)=scale*mode(n)
      EDIT ADD  MODE n [TO m] [SCALe real]       mode(m)=mode(m)+scale*mode(n)
      EDIT SET  MODE n SCALe real                mode(n)=scale
      EDIT MULT MODE n SCALe real                mode(n)=scale*mode(n)

If the destination is not specified then the result will be appended
as a new mode after the existing ones. The default value for the SCALe factor
is 1.0.

.. _vibran_basis:

Generate entire basis sets
--------------------------

The BASIs command will append requested basis vectors to an
existing (or null) set of orthonormal vectors. In this way, basis sets
for reduced diagonalizations, constrained normal modes, or other
calculations may be generated in a simple manner.

::

   BASIS   { IC  { FIRSt BOND    }       }  [NOORthonorm]
           {     { FIRSt ANGLe   }       }
           {     { DIHEdral      }       }
           {     { SECOnd ANGLe  }       }
           {     { SECOnd BOND   }       }
           {                             }
           { TR atom-selection [BALAnce] }

The BASIS command is similar to the EDIT INCLude ORTHog command,
except that many new vectors may be added. There are two modes for this
command. The IC mode will include a set of vectors based on the IC table.
The selection of which section of the table to use is required. The choices
are;

::

        FIRSt BOND
        FIRSt ANGLE
        DIHEdral
        SECOnd ANGLe
        SECOnd BOND

All valid IC table entries (i.e. all atoms of specified set defined) will
result in the addition of one vector to the normal mode basis, provided that
this vector is not a linear combination of existing vectors. See the
description for IC mode specification (:ref:`vibran_projections`) for a
description of these vectors.

The second mode is the TR (translation/rotation) mode which will
usually add the 6 translation rotation degrees of freedom for the selected
set of atoms. If the selected set of atoms is linear, then only 5 vectors
will be added. If only one atom is selected, then just the 3 translation
vectors are added (for this case all numvers but one will contain a zero).
The BALAnce keyword is suggested, and its operation it to remove the
net (whole system) rotation/translation components from the added vectors.
If all atoms are selected, the BALAnce keyword will result in NO vectors
being added.

The purpose of this command is to facilitate setting a vector
basis for a reduced normal mode calculation. Another purpose is in setting
up translation/rotation modes as initial vectors for the iterative
diagonalization procedure (no BALAnce option).

This command does normalize the new basis vectors, but it will not
modify any existing vectors. Each added vector is then orthogonalized
from all exisiting vectors (unless the NOORthonormalize keyword is specified).
If the vector has a zero norm following orthogonalization, it is rejected.
This avoids any possibility of basis interdependancies (zero determinant).
The norm after orthogonalization is saved in the eigenvalue array.
Each added vector is then normalized to form an orthonormal basis.
NOTE: this method will not work unless all exisiting vectors form an
orthonormal subbasis (i.e. before the BASIS command is invoked, the PRINT
NORM DOTProducts should give the unit matrix). The EDIT ORTHog command may
be used to orthonormalize a set of vectors.

The user is expected to process these modes with the REDUce command
after saving them, or in an external program.


.. _vibran_fill:

The Fill command
----------------

The comparison coordinates can be modified with the FILL DIFF
command. This command will copy the main coordinate set to the comparison
set, and then step along the specified mode by the specified magnitude.
When the append (APPE keyword) is used, the main coordinates are not first
copied.

::

   FILL    {  DIFF        } [mode-spec] [magnitude-spec] [APPE]
           {  COMP }

The FILL COMP command will fill the comparison coordinate
displacement arrays with the specified vector. For this option,
the comparison coordinate are zeroed before the displacement is added.
The coordinates will then contain a displacement vector. The append option
will prevent the zeroing of these arrays before stepping along the mode.


.. _vibran_second_derivatives:

Second Derivatives
------------------

The second derivatives are computed during energy
determination and they are stored in temporary arrays that are
allocated dynamically. Once obtained, they can be written out, or
diagonalized internally.

If 'QSECD=.TRUE.' the second derivatives of the energy
are returned in the array DD1. This array contains the
full upper half of the second derivative matrix, and
contains (NAT3*(NAT3+1)/2) REAL*8 elements (NAT3 = 3 * number
of atoms). Memory storage for the DD1 array may cause 
memory overflows, especially for large systems.

If one wants to get the normal modes for large structures,
it is often advantageous to use the DIMB or REDU CMPAct option.
These options will fill a DD1CMP second derivative array with only
the non-zero elements. The size of this array depends heavily
on the non-bonded cutoff distance, since increasing that
distance increases the number of non-zero elements. At this 
moment, the DD1CMP array will allow the average number of
interactions per atom to be as high as 300. This number will
allow non-bonded cutoffs up to 13A without problems, for all-
hydrogen systems. 
A comparison between allocated memory space for DD1 versus
the maximum memory space for DD1CMP is given below:

=============== ==================  =========    ===========
Structure       #atoms(explicit H)   DD1 (MB)    DD1CMP (MB)
=============== ==================  =========    ===========
Deca-alanine             66              0.2            1.4
BPTI                    568             11.6           12.3
Lysozyme               1264             57.6           27.4
Hemoglobin             5600          1,129.2          121.2
=============== ==================  =========    ===========


.. _vibran_block_normal_mode_method:

Block Normal Mode (BNM) Analysis
--------------------------------

The Block normal mode method is based on the original work of
Tama and co-workers. BNM projects the full atomic hessian into a 
subspace spanned by the eigenvectors of blocks, which reduces the size 
of the eigenvalue problem dramatically. Each block can be defined by 
the user in very flexible manners as an amino acid or a secondary 
structural element.  Currently only the T/R vectors of the blocks are 
included as basis vectors, which therefore reduces the eigenvalue 
problem from 3Nx3N to 6n*6n; N is the total number of atoms, and n is
the total number of blocks. In the future, it might be of interest
to include other low-frequency eigenvectors of the blocks into the
definition of the subspace, which would be more robust when the blocks
are large in size. 

Compared to the original RTB work of Tama et al., the current
implementation has the following improvements:

1. The projected hessian was constructed in a direct manner; 
   i.e., the full atomic hessian was never stored, which is essential for
   large systems.
   
2. A Lanczos algorithm was adapted for diagonalizing the
   projected hessian, which can be very large but sparse for 
   super-molecular assemblies (e.g., ribosome). 

3. BNM calculations can be carried out in parallel mode.

Keywords and options:

::

   BHES {SERL} GENR [TMEM int MEMO int MEMA int]
        {PARA}      NNOD int
   POST FLAG 1

SERL/PARA is the flag for serial or parallel computations. 
GENR is the keyword for using one residue per block.
TMEM is the total memory (in MB) available to BNM calculations.
MEMO is the memory (in MB) allocated for other arrays (recommended value is 20)
MEMA is the memory (in MB) allocated to contruct super blocks. 

The most important variable is MEMA.
MEMA (in MegaB) can be estimated by the following expression 
(nres is number of residues):

::

   MEMA=6*nres*(6*nres+1)*8/2/10^6

In addition, it also depends on whether the P/ARPACK diagonalizer is 
used or not. Without P/ARPACK, the standard diagonalizer in CHARMM 
is used, which limits the size of the system that can be studied; only
MEMA is important. 

With P/ARPACK, user can set TMEM and MEMA based on the available 
resource and size of the system. It is possible to save the projected 
matrix on the hard drive (if MEMA>TMEM+MEMO) for diagonalization.
Obviously the P/ARPACK library must be available for compiling the code
if one plans to use P/ARPACK.

For parallel calculations, NNOD needs to be specified, 
which gives the number of nodes. The current version has only been 
tested with the following copilation option:

::

   ./install.com gnu SIZE M mpich

In addition, the keyword "VIBPARA" has to be specified in pref.dat.
During BNM calculations, the following files will be generated 

1. projected non-zero hessian matrix elements:hfinal-a-b.dat

   where a is the index of nodes, b is the index of files on the ath node.
   In these files, the first two columns are the Row and Column indeces of
   the projected hessian matrix elements, the third value is the derivative
   itself.

2. Normal mode eigenvectors:nmeigv-a-0.dat

   where a is the index of normal modes. 
   These files can be used to perform analysis (e.g., fluctuations) with
   other functions of VIBRAN using the POST keyword. 
   In the serial version, you can use BHES followed by "POST FLAG 1" 
   directly. In the parallel version, however, POST has to be done 
   separately.

3. The frequencies are reported in the file freq.out, which lists three 
   values for each mode: the first one is the eigenvalue itself, the
   second one is the frequency that corresponds to the eigenvalue, and the
   third one is SCALED frequency (by SCALe). As 
   discussed in the literature (see below), scaling is necessary because 
   the block approximation makes the vibrations more rigid; an appropriate 
   value is 0.5882 when one residue is treated as a block.

Finally, we note that in some systems, metal ions are present. Since a
single ion does not have rotational degrees of freedom, three redundent
vectors with zero eigenvalues are assigned. For example, for a system
with three metal ions, the BNM calculation would generate 15 
zero-frequency modes - six of which correspond to the T/R of the entire 
molecule, while nine (3*3) are redundent vectors. 

References
^^^^^^^^^^

* F. Tama et al., Proteins: Struct. Funct. Genet. 41: 1-7 (2000)
* G. Li and Q. Cui, Biophys. J. 83, 2457-2474 (2002)


.. _vibran_gaussian_network_model:

Gaussian (Anisotropic) Network Model
------------------------------------

The implemented GANM is based on papers published by the Bahar
group; the ARPACK library also works with GANM with our implementation, 
so that large systems can be studied efficiently.

The GANM calculation requires a user defined file named 
CHAINS.DAT, which contains the definition of segments of a biomolecule.

Syntax
^^^^^^

::

   GANM (atom selection) UNIT int {AISO} SERL int

==== ========================================================================
UNIT specifies the paramter file used in GANM calculation (stiffness,cutoff)
AISO if specified, the Anisotropic Network Model will be used; otherwise, the
     Gaussian Network model will be used.
SERL can be 1 or 2; 1 uses the standard diagonalization method, and 2 uses 
     ARPACK.
==== ========================================================================

atom selection: atoms included in GANM calcualtions.

Example of the parameter file (specified by UNIT):

::

   1 1                 : number of atom types, number of interaction types
   CA                  : atom type(s) included in ANM
   1 1 0.95 13.0       : index of atom type, index of atom type, force constant, 
                         distance cutoff

Example of the file CHAINS.DAT:

::

   2 148               : total number of  segments, total number of residues
   1 1   144           : the first  segment , starting and ending residues
   2 145 148           : the second segment , starting and ending residues

References:

1. GNM: Tirion, M. M. 1996, Phys. Rev. Lett. 77: 1905 
2. ANM: 

   * Atilgan, A. R., Durell, S. R., Jernigan, R. L., Demirel, M. C.
   * Keskin, O., Bahar, I., 2001 Biophys. J. 80: 505
   * Doruker, P., Atligan, A. R., Bahar, I. 2000, Proteins, 40: 512


.. _vibran_vsa:

Vibrational Subsystem Analysis (VSA) Method
-------------------------------------------

see test/c35test/VSA_Butane.inp

If this method is used please site... 

::

   Woodcock, HL; Zheng, W; Shao, Y; Kong, J; Brooks, BR. Vibrational Subsystem
   Analysis: A Method for Probing Free Energies and Correlations in the Harmonic
   Limit. To be submitted, 2008. 

The vibrational subsystem analysis (VSA) method is designed for coupling
global motion to a local subsystem. This method is a  partitioning scheme that
separates (and integrates out) the motion of the environment from the user
defined subsystem (see Method Section) while still allowing the environmental
motion to perturb the local subsystem dynamics. It was originally developed
for EN models but is now extended for all-atom representations and hybrid
quantum mechanical / molecular mechanical (QM/MM) potentials. Below is a brief
list of possible uses:

1. examination of local-global motion
2. performing accurate NMA while eliminating unwanted degrees of freedom
3. eliminating excess noise from large NMA (i.e. QM/MM)
4. performing NMA while not at a stationary point with respect to all degrees
   of freedom
5. integration of light particle during NMA (i.e. application to polarizable
   models)

