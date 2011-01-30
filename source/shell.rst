.. py:module::shell 

===================
Shell Decomposition
===================

It is often desirable to decompose solvent around a solute
into shells. This module allows such a decomposition based on a distance
criterion.


.. _shell_syntax:

Syntax for the SHELL command
----------------------------

::

   SHELL  [ NSHL int ]  [ NOIMages ]  [ ATOM ] -
          [ SOLUte atom-selection ]  [ SOLVent atom-selection ] -
          [ SHTH real | SHBO real ... ]  

          [ UPDAte ] 

          [ DEFIne ] [ name ]  [ SHELL int ]
                               [ BULK ]
                               [ SOLUte ]
                               [ SOLVent ]

          [ STATistics ]

          [ OFF ]


         atom-selection::= see *note Selection:(chmdoc/select.doc)

.. _shell_overview:

General overview
----------------

The SHELL module forms the basis of any shell based analysis
of solvent properties. It partitions all given SOLVent residues /
atoms into shells depending on the distance to the closest SOLUte atom
in the solute selection. The atom numbers of all atoms in a shell are
stored in a list which other modules, such as CORRel (see :doc:`correl`)
can then use to compute properties of interest.

If images are set up SHELL will use SOLUte and SOLVent image
atoms for the distance criterion as well (although this can be switched
off explicitly by using the NOIMages keyword).

By default SHELL marks all atoms in a residue (even if
not all are included in the SOLVent selection) as being in the same
shell (if more than one atom of one residue is in the SOLVent selection
the lowest shell "wins"). This behavior can be altered by the ATOM
keyword. If present, only the atoms which are explicityl in the selection 
will be partitioned.


.. _shell_setup:

SHELL decomposition setup
-------------------------

SHELL must be called once to set up the parameters before any
analysis can take place. The following parameters can be
used:

===============  ==================================================================
NSHL <int>       number of shells: into how many shells should the solvent
                 be divided (where every solvent atom beyond the last shell
                 is considered 'bulk')
                 default: 1

SHTH <real>      shell thickness: how thick (in A) is one shell (all shells
                 have the same thickness and an atom which is exactly (i *
                 SHTH) A from a solute atom is considered to be in the i+1th
                 shell).
                 
	              default: 4.75A

SHBO <real...>   an alternative way to define the shell thickness: here
                 a real number must be given for each shell. Each number is
                 taken as the boundary of one shell to the next (so the numbers
                 should be in increasing order). This allows for shells of
                 different thickness.

ATOM             (optional) don't put whole residues (like a whole TIP3 water)
                 into the same shell, but only assign the atoms present in
                 the SOLVent selection to a shell ignoring residue
                 connectivity.
                 
	              default: OFF

NOIMage          don't use image atoms in the decomposition process even if
                 images are set up (this will usually alter only shells at
                 the rim of the coordinates but if the solute itself is
                 located at the rim (where it can drift to during a
                 simulation) this may even alter 'inner' shells.
                 
                 default: ON if images are set up and the image structure is
                 filled (e.g. by a previous call to 'ENERgy')

SOLUte           all atoms specified in the following selection are considered
                 as solute and will determine into which shell a solvent atom
                 is placed. (no default)
                 
SOLVent          these atoms will be partitioned into the specified number of 
                 shells. By default, a whole residue is in the lowest (innermost)
                 shell in which one of its atoms specified in this selection is
                 located. This can be altered by the ATOM keyword. (no default)
===============  ==================================================================

.. _shell_update:

Updating SHELL
--------------

Once the shell parameters are set up, every time new
coordinates are read into the main coordinate set, the shell lists must
be updated by calling the 'SHELL UPDAte' command. The decomposition
is performed using a cubing based distance algorithm with a cube size of
NSHL * SHTH. 

.. _shell_define:

Define
------

The individual shells as well as bulk or the solute or solvent
atoms can be converted to named selections which can then be used as
atom selections (analogue to the 'DEFIne <name> select ... end'
command). This allows output of single shell coordinates or general
access to the shells via the atom selection method.


.. _shell_statistics:

Statistics
----------

The STATistics command prints the general shell parameter
currently in use and information about the shells.

At print level 6 or higher the atom numbers for all shells will
be printed and at print level 7 also the solute, solvent and bulk atom
numbers will be shown.


.. _shell_off:

Off
---

The :chm:`SHELL OFF` command  turns off the shell
module and frees all allocated space.


.. _shell_correl:

CORREL series using SHELL
-------------------------

Currently two timeseries using shell are implemented:

- ::
 
      SDIP: syntax: ENTER <name> SDIP [SHELL int]
                                      [BULK]
                                      
  This generates a series with 4 entries: the X/Y/Z component of the
  dipole moment and the number of atoms in this shell.

- ::

     SATM: syntax: ENTER <name> SATM [SHELL int] atom-spec
                                     [BULK]
                                     
  Generates a one-dimensional timeseries containing zero or one
  indicating whether the atom is in the given shell or in the bulk at a
  given timestep.


.. _shell_caveats:

Caveats
-------

- Currently only main coordinates are considered - there is no COMP
  keyword. Since SHELL is mainly devised to be used with CORRel and
  while reading trajectories this should not pose much of a problem
  (although implementation of a 'SHELL UPDAte COMP' command should not 
  require many changes).


.. _shell_efficiency:

Efficiency
----------

The core of the SHELL UPDAte command is a cubing based algorithm 
(adapted from M. Crowley's images/nbndgcm.src)
which is aware of SOLUte and SOLVent and thus calculates only distances
between SOLVent and SOLUte atoms. The efficiency of this approach is
influenced largely by the following parameters:

- SHTH and NSHL: the size of one cube is defined by SHTH * NSHL. If this
  distance is large cubing will start to be inefficient and a brute
  force approach may be more efficient.

- The number of atoms in the SOLVent and SOLUte selections: if only one
  atom of a residue is considered 'decisive' for a whole residue the
  number of pairs to sample can be greatly reduced. E.g., if only the
  oxygen of TIP3 water is selected as SOLVent the number of SOLVent
  atoms is reduced by a factor 3 while by default all water atoms are marked
  as belonging to a shell (unless ATOM is set; then it is necessary to select
  all water atoms, see examples).


.. _shell_examples:

Examples
--------

::

   SHELL NSHL 4 SHTH 3.5 SOLU SELE RESNAME ALA END -
                         SOLV SELE RESNAME TIP3 .AND. TYPE OH2 END

Puts all TIP3 waters into 4 shells - each 3.5A thick - around all
ALAnines. Only the TIP3 oxygen (OH2) are considered for the solvent so
each water will be in the shell its oxygen is in.

::

   SHELL NSHL 4 SHTH 3.5 ATOM SOLU SELE RESNAME ALA END -
                              SOLV SELE RESNAME TIP3 END

the same as above but every atom is considered individually so at a shell
interface one or both hydrogens can be in one shell while the oxygen is
in the other.
