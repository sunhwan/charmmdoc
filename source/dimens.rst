.. py:module:: dimens

=============================================
How to set CHARMM to run with any size system
=============================================

There are two ways to change the run-time size of charmm
arrays to accomodate a system of any size: from command line, and from
the charmm input (input file or standard input). Currently only the
maximum number of atoms can be set, and other sizes like maximum
number of residues, bonds, angles, segments, and so on, are set by an
approximation based on number of atoms specified.

Command line:

::

    $ charmm -chsize nnnnnnnn

where ``nnnnnnnn`` is the maximum number of atoms in your system, for instance 200000.

Input file:

Use the charmm command dimension (only dime needs to be specified) immediately after the title:

::

    * title
    *
    dimension chsize 200000

In addition to the overall size, the DIMENSION command can also be used to set
the sizes of specific subarrays in the following list:

Data Structure Size

=================== ====================================================
chsize              This is a master size that proportions all CHARMM
                    data structures
maxai  (chsize)     This controls the maximum number of atoms
maxb   (chsize)     Maximum number of bonds
maxt   (chsize*2)   Maximum number of angles
maxp   (chsize*3)   Maximum number of proper dihedral angles
maximp (chsize/2)   Maximum number of improper dihedral angles
maxnb  (chsize/4)   Maximum number of nonbond fixes
maxpad (chsize)     Maximum number of acceptors and donors
maxres (chsize/3)   Maximum number of residues
maxseg (chsize/8)   Maximum numebr of segments
maxcrt (chsize/3)   Maximum number of CMAP dihedrals
maxshk (chsize)     Maximum number of SHAKE constraints
maxaim (chsize*2)   Maximum number of atoms including images
maxgrp (chsize*2/3) Maximum number of groups
=================== ====================================================

Example:

::

    * title
    *
    dimension maxa 10000 maxp 10 maximp 10

