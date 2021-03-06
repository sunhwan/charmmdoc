.. py:module:: torque

============================
Manipulating torques: TORQUE
============================

::

    TORQue { SET   } { WATEr } [ ALL       ]  [atom-selection]
           { ADD   } { BODY  } [ BYSEgment ]
           { CLEAr }           [ BYREsidue ]
                               [ BYGRoup   ]

Meaning of individual keywords:

* Actions:

  ======== ========================================================
  SET      Sets (defines) torque centers for the select atoms.

  ADD      Adds new torque centers for the selected atoms, leaving
           existing torque centers intact.

  CLEAr    Deletes all existing torque centers.
  ======== ========================================================

* Body types:

  ======== =========================================================
  WATEr    Indicates that the selected atom(s) are waters.

  BODY     Indicates that selected atom(s) are arbitrary types and
           the rigid body may consist of any number of atoms in any
           shape.
  ======== =========================================================

* Body construction:

  ========= ========================================================
  ALL       Each atom in the selection should be given its own
            torque center (located at the atomic coordinates).

  BYSEgment Each segment is considered its own rigid body.

  BYREsidue Each residue is considered its own rigid body.

  BYGRoup   Each group is considered its own rigid body.
  ========= ========================================================


.. _torque_examples:

Examples
--------

::

    torque set water byres sele segid bwat end

Defines torque centers for all water molecules of the segment labeled
BWAT.

::

    torque add water bytes sele segid wat2 end

Defines torque centers for the waters in the wat2 segment (retaining any
previously defined torque centers).

::

    torque clear

Clears all torque centers.


.. _torque_notes:

Notes
-----

Currently, the TORQue command can only be used to define torque centers
for use with the TORQue option to MSCALe (see mscale.doc for details).
This allows the 3x3 rotation matrix of each torque center to be passed
to slaves and the 1x3 torques to be returned to the master process.

The only body type currently allowed is the WATEr type; selecting
"BODY" will produce an error. Likewise, only the BYREsidue body construction
has been tested (others may work).

Currently, coordinates are not assigned explicitly to torque centers, only
internal rotations and torques are stored (as this is all that is required
by the MSCALE implementation).
