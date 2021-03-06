.. py:module:: lonepair

==================
Lone Pair Facility
==================

This routine parses the lone-pair command which converts existing
atoms to lone-pairs in the PSF.

Bernard R. Brooks, NIH, October, 1997

.. _lonepair_syntax:

Syntax of the Lone-Pair Command
-------------------------------

::

   LONEpair { FIXEd   atom-spec   [ xloc-yloc-zloc ]            } [MASS]
            {                                                   }
            { CENTer  atom-spec  {  atom-selection   }          }
            {                    { repeat(atom-spec) }          }
            {                                                   }
            { COLOcate { 2x(atom-selection) }                   } 
            {          { 2x(atom-spec)      }                   } 
            {                                                   }
            { { COLInear distance-spec } { 3x(atom-selection) } }
            { { CEN2                   } { 3x(atom-spec)      } }
            {                                                   }
            { { RELAtive } { 4x(atom-selection) } position-spec }
            { { BISEctor } { 4x(atom-selection) }               }
            { { CEN3     }                                      }
            {                                                   }
            { PRINt                                             }
            { CLEAr  [MASS real]                                }

   atom-spec::= { residue-number atom-name }
                { segid  resid atom-name   }
                { BYNUm  atom-number       }

   atom-selection ::= see *note select:(chmdoc/select.doc).

   xloc-yloc-zloc::= three real numbers identifying the new position

   distance-spec::= [DISTance real] [SCAled real]
                      (def 0.0)       (def 0.0)

   postition-spec::= [DISTance real] [ANGLe real] [DIHEdral real]
                        (def 0.0)      (def 0.0)    (def 0.0)


.. _lonepair_description:

Note on the lone-pair command
-----------------------------

1. The "LONEpair FIXEd" command places atoms that are fixed in
   the unit cell fractional coordinates.  If running constant pressure,
   these atoms will move in response to the changes in the box size/shape.
   Thus, this is different than "CONStraint FIX".  The specified position
   is in cartesian space (Angstroms).

2. The "LONEpair CENTer" command places a single particle at the
   weighted center of all the subsequently specified atoms (either by
   atom-selection or atom-spec).

3. For all commands employing multiple atom selections, each atom
   selection MUST CONTAIN THE SAME NUMBER OF ATOMS.  The atoms are
   then matched off in sequential order.  The first atom selection is
   the list of lonepair atoms.  This is intended to make it easy
   to create a large number of TIP4P water molecules with one command:

   ::
   
      LONEpair BISEctor DIST 0.15 ANGLE 0.0 DIHE 0.0 -
                           SELE ATOM SOLV * OM  END - 
                           SELE ATOM SOLV * OH2 END - 
                           SELE ATOM SOLV * H1  END - 
                           SELE ATOM SOLV * H2  END

   This assumes that all of the residues in the segment SOLV are TIP4P
   types.  It places the atom OM (with zero mass) at a point 0.15 A
   from the atom OH2 in the direction of the H1-H2 bisector.

4. The MASS option is used in the CEN* commands to determine how the
   center position is computed.  The default is to use the center of
   geometry (unit weights for all atoms).

5. The CLEAr command will remove all lone-pairs from the PSF and
   resets the lone-pair facility.  WARNING:  The masses of old lone-pairs
   will still be zero.  These may be modified "SCALAR MASS SET..."
   before running further dynamics, or by using the MASS keyword value
   whereby the mass of removed lonepairs are assigned to this value.

6. Lone-pair atoms may have other lone-pair atoms as a host, provided
   that the host is not already a lonepair.  In other words, you can define
   the postition of lone-pair B based on the postition of lone-pair A
   ONLY if lonepair B is define before lonepair A.  ORDER IS IMPORTANT!!
   
7. The LONEpair command sets the MASS to zero of all selected
   lone-pairs.  This MAY change the total mass of the system.  The
   lost mass is NOT added to any other atom.

8. For the BISEctor option, the dihedral is based on: I,J,(K+L)/2,L
   where I,J,K,L are the coordinate vectors of the specified atoms.

9. For the COLInear option, the DISTance value is a signed value
   of the distance from the first host AWAY from the second host.
   The SCALed value is a relative distance from the first atoms AWAY
   from the second atom (a SCALed factor of -1.0 will put the lonepair
   at the position of the second atom).  For example, the following
   two commands do exactly the same thing:
   
   ::
   
      LONEpair COLINEAR DIST 0.0 SCALE -0.5  -
               SELE type HB  END SELE type H1  END SELE type H2  END
      LONEpair CEN2  - 
               SELE type HB  END SELE type H1  END SELE type H2  END

10. When running CHARMM in parallel, the lonepair atom should be in
    the same group as ALL of its hosts in order to ensure that these
    atoms are all in the same parallel partition.  THIS IS NOT CHECKED!

11. Lonepair data is considered to be part of the PSF.  When a PSF
    with lonepairs is read from a file, the lonepair facility is also
    read (or appended).

12. When using lonepairs with PERT, both PSF's MUST have the same
    lonepair data.  In other words, you cannot perturb a lonepair into
    a non-lonepair.

13. Lone pairs may now be set in the RTF using a syntax similar to the
    standard lonepair command (lonepair.doc).  All options for the
    generation of lonepairs (FIXed, CENTer, COLOcate etc.) as specified
    in the lonepair documentation may be used although the atom selection
    specification is simplied as shown in the following example.

    ::
    
       LONEPAIR relative LPA O  C  CL distance 0.3 angle 91.0  dihe 0.0
       LONEPAIR relative LPB O  C  N  distance 0.3 angle 91.0  dihe 0.0

14. Lone pairs with undefined coordinates can be built by COOR SHAKE.
