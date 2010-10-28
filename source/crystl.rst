.. py:module:: crystl

####################################
Calculation on Crystals using CHARMM
####################################

The crystal section within CHARMM allows calculations on crystals to be performed. 
It is possible to build a crystal with any space group symmetry, to optimize its 
lattice parameters and molecular doordinates and to carry out a vibrational analysis 
using the options.

.. _syntax:

Syntax
------

.. py:class:: CRYStal

   .. py:method:: BUILd_crystal [CUTOff <real>] [NOPErations <int>]
   .. py:method:: DEFIne xtltype
   .. py:method:: PHONon

   :param xtltype: crystal type
   
.. _function:

Function
--------

.. _crystal_build:

Crystal Build
^^^^^^^^^^^^^

A crystal of any desired symmetry can be constructed by repeatedly applying a small number of transformations to an asymmetric collection of atoms (called here the primary atoms). The transformations include the primitive lattice translations *A*, *B* and *C* which are common to all crystals and a set of additional transformations, *T*, which determines the space group symmetry.

The Build command will generate, given *T*, a data structure of all those transformations which produce images lying within a user-specified cutoff distance of the primary atoms. The data structure can then be used by CHARMM to represent the complete crystal of the system in subsequent calculations. The symmetry operations, *T*, are read from the lines following the **CRYStal BUILd** command.

The syntax of the commmand is :

::

  CRYStal BUILd CUTOff <real> NOPErations <int>
  ! <int> lines defining the symmmetry operations.

The **CUTOff** parameter is used to determine the images which are included in the transformation list. All those images which are within the cutoff distance are included in the list.

.. Note::
   The distance test is done based on the atoms that are currently present and their symmetric representation.
   
To generate a crystal file from a box with a single atom at the center, the cutoff value will nee to be larger than the box dimensions.  If the box is filled with water and only nearest neighbor cells are desired, then the cutoff distance should be comparable to the **CUTIm** value (see :ref:`Image Updates <image_update>` ) or the CUTNB value (see :ref:`NBONDS Syntax <nbonds_syntax>`). There is no limit to the number of transformations included in the lists as they are allocated dynamically, but having too many will slow the image update step.

The crystal symmetry operations are input in standard crystallographic notation. The identity is assumed to be present so that (X, Y, Z) need not be specified (in fact, it is an error to do so). For example, a P1 crystal is defined by the identity operation and so the input would be

:: 

   CRYStal BUILd .... NOPEr 0

whilst a P21 crystal would need the following input lines :
                          
::

   CRYStal BUILd .... NOPEr 1
   (-X,Y+1/2,-Z)

A P212121 crystal is specified by NOPEr 3

::

   CRYStal BUILd .... NOPEr 3
   (X+1/2,-Y+1/2,-Z)
   (-X,Y+1/2,-Z+1/2)
   (-X+1/2,-Y,Z+1/2)

It should be noted that in those cases where the atoms in the asymmetric unit have internal symmetry or in which a molecule is sited upon a symmetry point within the unit cell not all symmetry transformations for the crystal need to be input. Some will be redundant. It is up to the user to check for these cases and modify the input accordingly.

.. _examples:

Examples
--------

.. _implementation:

Implementation
--------------
