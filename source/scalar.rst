.. py:module:: scalar

======================================================
SCALar : commands to manipulate scalar atom properties
======================================================

::

   SCALar keyname  {                       }        [atom-selection]
                   {  =   keyname          } ! A = B
                   { COPY keyname          } ! A = B
                   { SUM  keyname          } ! A = A + B
                   { PROD keyname          } ! A = A * B
                   { SET <real>            } ! A = <real>
                   { ADD <real>            } ! A = <real> + A
                   { MULT <real>           } ! A = <real> * A
                   { DIVI <real>           } ! A = A / <real>
                   { SIGN                  } ! A = sign ( A )
                   { INTEger               } ! A = int ( A )
                   { RECIprocal            } ! A = 1/ A
                   { LOG                   } ! A = ln ( A )
                   { EXP                   } ! A = exp ( A )
                   { ABS                   } ! A = ABS ( A )
                   { NORM                  } ! A = A / 2-norm(A)
                   { MIN <real>            } ! A = MIN (A,<real>)
                   { MAX <real>            } ! A = MAX(A,<real>)
                   { POWEr <real>          } ! A = A ** <real>
                   { POW2r                 } ! A = A * A
                   { IPOW  <real>          } ! A = A ** int(<real>), OK for neg A
                   { SQRT                  } ! A = SQRT(A)
                   { RANDom                } ! A = random
                   { HBCOunt               } ! A = #of hbonds
                   { SHOW   [SORT]         }
                   { STATistics weight_opt }
                   { STORe  store_number   } ! S(i) = A(i)
                   { RECAll store_number   } ! A(i) = S(i)
                   { +STOre store_number   } ! S(i) = S(i) + A(i)
                   { *STOre store_number   } ! S(i) = S(i) * A(i)
                   { READ <unit>           }
                   {                       }
                   {         [ ALL       ] }
                   { AVERage [ BYSEgment ] } ! S(i) = sum(S(j))/Nj
                   {         [ BYREsidue ] } !   averaged over each selected
                   {         [ BYGRoup   ] } !   atom of each item
                         weight_opt

   store_number ::= any number between 1 and 9. store_number 
   weight_opt ::= [ WEIGht store_number ] [MASS]


   keyname ::=     { X                  }
                   { Y                  } ! main coordinates
                   { Z                  }
                   { WMAIn              } ! main coordinate weights
                   { XCOMp              }
                   { YCOMp              } ! comparison coordinates
                   { ZCOMp              }
                   { WCOMp              } ! comparison coordinate weights
                   { DX          psf_no } ! forces from last energy eval
                   { DY          psf_no } ! or force difference from last EPERT
                   { DZ          psf_no } !  eval when using PSF 0
                   { ECONt              } ! Energy partition array
                   { EPCOnt             } ! Free energy difference atom partition
                   { MASS        psf_no } ! atom masses
                   { CHARge      psf_no } ! atom charges
                   { CONStraints psf_no } ! harmonic constraint constants
                   { XREF        psf_no }
                   { YREF        psf_no } ! reference coordinates
                   { ZREF        psf_no }
                   { FBETa              } ! friction coefficients
                   { MOVE        psf_no } ! rigid constraints flag
                   { TYPE        psf_no } ! atom chemical type codes
                   { IGNOre             } ! ASP flag for ignoring atoms
                   { ASPValue           } ! ASP parameter value
                   { VDWSurface         } ! ASP van der Waals surface
                   { RSCAle      psf_no } ! Radius scale facor for vdw
                   { WCAD        psf_no } ! scale facor for WCA potential
                   { ALPHa  (**) psf_no } ! atom polarizability
                   { EFFEct (**) psf_no } ! effective number of electrons
                   { RADIus (**) psf_no } ! van der Waal radii
                   { FQPRin (**)        } ! FlucQ principal quantum number
                   { FQZeta (**)        } ! FlucQ Slater orbital exponent
                   { FQCHi  (**)        } ! FlucQ electronegativity parameter
                   { FQMAss (**)        } ! FlucQ charge mass
                   { FQJZ   (**)        } ! FlucQ self-interaction
                   { FQCForce           } ! FlucQ charge force
                   { FQOLd              } ! FlucQ charges from last timestep
                   { SCAx (x::=1,2,..,9)} ! specific scalar store array
                   { ONE    (**)        } ! vector with all 1's
                   { ZERO   (**)        } ! vector with all 0's
                   { VARC               } ! Offset for the Variable LJ cutoffs
               
   psf_no::=  {  PSF 0  } ! which PSF (use only with PERT)
              { [PSF 1] }

   (**) =  For the keynames labeled (**), the array values may not
   be modified by any scalar command, but they may be used in the
   SHOW, STORe, STATistics, or as any second keyname (e.g. COPY).
   In order to change the "ALPHa", "EFFEct", or "RADIus" value, one must
   change the atom's "TYPE" value, which in turn determines these values.

All of the SCALar commands allow an atom selection.
The dimension of all vectors is equal to the number of atoms. (see the
vibrational analysis section for 3N vector manipulations).

The READ option reads from an ascii file all selected entries,
one value per line.  This file should have a valid value within the first
40 characters of each line, and should not have a CHARMM title.

The STATistics and AVERage options accept a weighting option. If
not specified, all atoms will have an equal weighting. If specified, the
weighting of each atom will be in the specified stored vector. For the
AVERage option, the averaging may be done over ALL selected atoms, or over
groups of selected atoms defined be segment, residue, or group boundaries.

For example if one wants to print out the RMS mass weighted 
fluctuations for the sidechain of each residue, the sequence should
be something like;

::

   COOR DYNA ....          ! put the rms fluctuations in the weighting array
   SCALar MASS STORe 1     ! put the masses in the first storage vector
   SCALar WMAIn IPOW 2     ! convert rms fluctuation to mean squared fluctuations.
   SCALar WMAIn AVERage BYREsidue WEIGht 1  SELE .NOT.
      ( TYPE CA .OR. TYPE C .OR. TYPE N .OR. TYPE H .OR. TYPE O ) END
           ! does the averaging for the sidechain with mass weighting
   SCALar WMAIn SQRT      ! take root to get rms of square fluctuation data.
   SCALar WMAIn SHOW SELE TYPE CB END  ! show the results for each residue

.. note::

   to get the mass weighted average rms fluctuations, remove steps 3 and 5.

The RANDom option sets a different random number for each selected
atom based on the specifications given in the RANDom command
(see :ref:`RANDom <miscom_random>`).

Example on how to use the Variable Lennard-Jones Cutoff Method (VLJCM).
VLJCM allows the user to apply different length lennard-jones cutoff's to
different atoms in the system. This is useful when using force-fields which
are parameterized for a specific LJ cutoff distance (eg EEF1) and one wants
to introduce new types of molecules which need larger LJ cutoffs. In this
situation VLJCM allows the user to keep the original LJ cutoff for a certain
set of atoms in the system, and apply a different LJ cutoff to different
atoms in the system. Here is an example of how to set up VLJCM before
dynamics is run. Assume you have a binary mixture of atoms labeled A and B:

::

  ! The nonbonded options below are part of the model
  update ctonnb 7. ctofnb 9. cutnb @cutnb cutim @cutim group rdie

  ! scale nonbond interactions
  scalar varc set 10.0 sele resname A end
  scalar varc show

In the above example the default LJ cutoff is 9 angstroms and the switching
region starts at 7 angstroms. We set up VLJCM using the 'scalar varc' command.
In this case the new LJ cutoff will be computed based on the 10.0 angstrom
argument supplied.

If there is a B-B LJ interaction the default LJ cutoff is applied. If there is 
an A-A LJ interaction the default LJ cutoff is computed as follows:

::

    ctonnb=(varc(A)+varc(A))/2
    ctofnb=(varc(A)+varc(A))/2 + cdiff

where 'varc(A)' = 10.0 angstroms, and 'cdiff' = 2 angstroms (taken from the
default size of the switching region). Thus for A-A LJ interactions ctonnb=10.
ctofnb=12.

In the case of A-B interactions we have 'varc(A)' = 10.0 angstroms and
'varc(B)' = 0.0 angstroms, which leads to ctonnb=5. ctofnb=7. 
