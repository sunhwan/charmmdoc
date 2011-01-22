.. py:module:: monitor

=================================================================
Monitor commands: Commands to monitor various dynamics properties
=================================================================

.. _monitor_syntax:

Syntax of the MONItor commands
------------------------------

::

   MONItor {DIHEdral} [SHOW] FIRSt unit-number NUNIt integer BEGIn integer -
                      STOP integer SKIP integer [SELEct atom-selection]

======= =================================================================
FIRSt   the unit number of the first file of dynamics coordinate sets
        from which the property is to be calculated.

NUNIt   the number of units of dynamics coordinate files.  Fortran unit
        numbers must be assigned to the files consecutively from FIRST.

BEGIn   the first step number for the coordinate set from which
        the property will be calculated.

STOP    the last step number for the coordinate set from which
        the property will be calculated.

SKIP    the time increment between the step numbers of the coordinates.

SELEct  selected atoms for which the property is to be monitored.  At
        this time, atoms may be selected only by the atom-selection
        keywords (e.g. RESID,TYPE,ATOM,RESN,SEGID) and NOT by
        tag-selections.  (see :doc:`select`)

DIHE    Property: monitor the dihedral transitions.

SHOW    for monitoring dihedral transitions, print out the step number,
        the cumulative number of transitions, the dihedral name, the
        current dihedral angle, and the old and new minimum well
        positions each time a transition is found.

ALL     Lots of printout.

UNIT    Unit number to write results (default: outu)
======= =================================================================


.. _monitor_properties:

Properties monitored using the MONItor commands
-----------------------------------------------

* DIHE: Dihedral transitions are monitored for any dihedral angle
  which can be made from the atoms selected.  A transition is defined as a
  change in the dihedral angle which results in going from one well of the
  torsion potential to another well, AND which involves crossing at least 30
  degrees beyond the barrier at the potential maximum.  That is, for rotation
  about a bond between tetrahedral carbons, the minima are at +60, 180 and -60,
  while the maxima are at 0, +120 and -120.  For an initial angle of +45, a
  transition is counted if the angle becomes > +150 or < -30.  The old minimum
  was +60, and the new minima would be 180 or -60, respectively.  The angle can
  change by as much as 120 degrees or as little as 60 degrees in going from one
  well to the next using this algorithm.

  For bonded atoms which both have trigonal geometry, the minima are
  +90 and -90, and a transition requires crossing 0 +- 30, or 180 +- 30 degrees.
  Only transitions for dihedrals with either 2 or 3 periodicity can be counted
  with the MONIt command.

  A word of caution:  the above algorithm for counting transitions is
  by no means fool proof, therefore one should always look at the dihedral time
  series to obtain a more precise number of transitions.  This is particularly
  true for mainchain phi and psi dihedrals which frequently have average
  positions which are not close to the minima for a tetrahedral atom. Large
  fluctuations can therefore be mistakenly (in a classical butane-type
  transition) counted as transitions.
