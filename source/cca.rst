.. py:module:: cca

Common Component Architecture
-----------------------------

by Milan Hodoscek, and others... (milan@helix.nih.gov,milan@cmm.ki.si)

CCA (Common Component Architecture) specification started by
the need of interfacing a variety of computational chemistry codes, ie
GAMESS and CHARMM. For details see J. P. Kenny, et al, J. Comp. Chem.,
25, 1717-1725, 2004.


.. _cca_description:

* See J. P. Kenny, et al, J. Comp. Chem., 25, 1717-1725, 2004.


.. _cca_using:

Nothing here yet....

.. _cca_installation:

Installation
------------

Nothing yet

.. _ccs_status:

{GAMESS,GAMESS-UK,Q-CHEM}/CHARMM interface status (November 2004)

The CCA is more general effort then just interfacing some QM
program with the CHARMM. For example user of CCA would like to combine
some python code of his own and use CHARMM for setup and analysis of
the system while performing PME dynamics in parallel with
NAMD. However the good starting point can be CCA-ing one of the
current ab inito packages with the CHARMM program.

Currently GAMESS and GAMESS-UK programs communicate via
fortran COMMON blocks with CHARMM, while Q-Chem is interfaced using
external data files. Here are some details:

List of COMMON blocks that are shared in GAMESS and CHARMM:
-----------------------------------------------------------

I.  GAMESS's COMMON blocks that are currently used in CHARMM (not for
    DIESEL or GAMESS-UK):

    ::
    
       COMMON /FUNCT/  (energy and forces)
       COMMON /INFOA/  (atomic numbers, etc)
       COMMON /COORDN/ (coordinates)

II. CHARMM's COMMON blocks (defined and documented in gamess.fcm) that
    are in use in GAMESS code:

    ::
    
       COMMON /GAMESL/ (control program flags)
       COMMON /CHMGMS/ (MM atoms for GAMESS)


All the data structures for the 3 interfaces are in fcm/gamess.fcm
