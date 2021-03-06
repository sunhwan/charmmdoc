.. py:module:: images

######
Images
######

Original implementation by  Bernard R. Brooks, 1983

CHARMM has a general image support system that allows
the simulation of almost any crystal and also finite point groups
(such as dimers and tetramers...). There is also a facility to introduce
bond linkages (with additional energy terms including angles, dihedrals
and improper dihedrals) between the primary atoms and image atoms.
This allows infinite polymers, such as DNA to be studied.
For infinite systems, an asymmetric unit may be studied because
rotations and reflections are allowed transformations.

The IMAGE facility is invoked by reading an image transformation
file. From this point, the images of the primary atoms will be included
in any energy and force determinations for the remainder of the calculation.
A null image file with the :chm:`INIT` keyword will disable this facility.

The simple periodic boundary code is underdevelopment by
Charles L. Brooks, III at the Scripps Research Institute as of Spring
1995.

.. index:: image; read
.. _images_read:

Image Transformation File
-------------------------

The image transformation file contains all of the information needed to define
the position and orientation of all symmetric images of the primary atoms.
The file is read in free field and may be passed parameters.

::

   READ IMAGe  [UNIT integer] [INIT] [PRINT]

File Structure:

::

   (title)

   SCALE  xfactor  yfactor  zfactor
   IMAGE  image_name                       ! start of a new transformation
   DEFIne repeat( [INVErse] other_image )  ! define transformation from others
   ROTAte xdir  ydir   zdir   angle        ! specify an axis and angle
   TRANslate  xdir  ydir  zdir [distance]  ! specify a displacement
   NEGAte                                  ! invert through origin
   END                                     ! terminates transformation file

When an image file is read, any existing image transformations
are discarded, but not any information regarding image patching.
The :chm:`INIT` keyword on reading will remove ALL existing image data first.
Rereading an image transformation file without the :chm:`INIT` keyword is
useful when crystal parameters are to be modified, but the patching
data is to be retained. The :chm:`PRINT` option, prints all data as it is read.

The image file starts with a standard CHARMM title. The remaining
commands are processed sequentially.

The :chm:`SCALE` command gives three values which multiply all
subsequent transformation specifications. The default values are unity.

The :chm:`IMAGE` command initiates a new image transformation. The
transformation matrix is set to unity (no rotation, no translation).
This transformation matrix is then modified by subsequent commands until
another :chm:`IMAGE` commmand or the :chm:`END` command is found.

The :chm:`DEFIne` command multiplies the current transformation matrix
by any of the previously defined images. The :keyword:`INVErse` keyword proceeding
the other transformation name, uses the inverse of this transformation.
Any number of previous transformations may be listed with this command,
but they are processed sequentially (just in case they don't commute).

The :chm:`ROTAte` command causes the current transformation matrix
to be operated with a rotation. Four real numbers must follow which
define a rotation axis and an angle about this axis (in degrees).

The :chm:`TRANslation` command will translate the current transformation
matrix. If three values are specified, then this is used as the translation
vector. If four values are given, then the first three define a direction,
and the fourth value defines a distance. Before operating on the
transformation matrix, the elements of this vector are multiplied by
the current scale factors (from last :chm:`SCALe` command).

The :chm:`NEGAte` command projects the current transformation through
the origin. This operation changes the chirality of the system, and
is not appropriate for macromolecules. This operation is required
for simulations with glide planes.

The :chm:`END` command is used to terminate the IMAGE file. This is
required if the file is read from the input stream.

.. warning::
   One restriction on the transformations is that every
   transformation MUST have an inverse. There is a serious warning
   if this restriction is violated. This requirement is needed in defining the
   energy of the system. When computing the energy, the Hamiltonian is
   assumed to be symmetric, and only the lower half is generated. The result
   of having an image without a transformation is to remove the symmetry of
   the Hamiltonian. The considerations of program efficiency and memory
   requirements make this necessary. There may be examples where this is
   desired, such as for cases where no energy calculations are needed
   or for structural analysis. A transformation may be its own inverse
   as is the case when a transformation consists of only a 180 degree rotation.

The maximum number of allowed transformations is 100. This limit
can easily be increased.


.. index:: image; write
.. _images_write:

Image Writing and Printing
--------------------------

Several different types of image data may be written or printed.
These are used for analysis and to check the operation of the program.

::

   WRITE IMAGes  { TRANsformations   }  [UNIT integer]
                 { PSF               }
                 { FORCes            }

The :chm:`TRANsformation` option will list all image transformation
matricies as well as what the inverse transformations are.
For each transformation, there is given a 3 X 3 rotation matrix followed
by the translation vector. For the use in this program the translation
is done **AFTER** the rotation has been made.

The :chm:`PSF` option, lists information about the image atoms as
well as list all primary-image internal coordinates (bonds, angles,
dihedrals, and improper dihedrals).

The :chm:`FORCe` option, lists the total force and torque each image
transformation applies on the primary atoms. This data may be used to
estimate the pressure of a system, or to check if minimization is
complete. At the end, the total force (vector sum) and torques are
listed.

.. index:: image; update
.. _images_update:

Image Updating
--------------

.. index:: keyword; ixtfrq, keyword; imgfrq

The crystal build procedure has to be done prior to image, nonbond or
hydrogen bond updating with :chm:`IXTFrq` frequency - it can be also done manually
by issuing :ref:`crystl_build` command. Actually first time is has to be done
manually to read parameters associated with crystal build command.

The image update procedure has several functions. This updating
is done prior to any nonbond or hydrogen bond updating, because its
results may affect those updates.

::

   {                                                           } ! no change
   { IXTFrq int IMGFrq int [CUTIm real] [IMALl  ] [INVErse  ]  }
   {                                    [IMBRief] [NOINverse]  }
   { IMGFrq  0                                      } ! suppress image updating


The absence of the :chm:`IMGFrq` keyword, maintains the current status
of image updating. Specifying an :chm:`IMGFrq` value of zero, suppresses all
image update functions, but does not modify the image lists in any way.

The :chm:`IMGFrq` integer value gives the frequency of image updating
to use during dynamics or minimization. For setting up a single image
update, any positive value may be used. When negative value was used,
e.g., ``INBFRQ = -1``, all lists are updated when necessary (heuristic test).

.. index:: keyword; cutim, keyword; imbrief, keyword; imall

The :chm:`CUTIm` value gives the maximum
allowable distance of any group to be included in the image atom lists.
Normally, a group is included only if it belongs to a transformation
whose inverse transformation is of a higher index than its own.
This is because only the lower triangle of the Hamiltonian is computed
and any image interaction between primary atoms and image atoms of a
higher inverse index will already be computed. This efficiency consideration
greatly reduces the required number of image atoms and the size of the
image nonbond lists. This reduction is activated by the use of the
:chm:`IMBRief` option (default). If on the other hand, one desires these groups for
the purpose of analysis or for displays, the IMALL keyword may be used
to generate all image atoms within the :chm:`CUTIM` distance.

The sequence of events in this update are;

1. Save existing image atom lists (from the previous update).
2. Process image centering if requested to replace far off
   groups of atoms by a closer image.
3. Generate appropriate image atoms within the cutoff distance
   of the primary atoms.
4. Remap internal coordinate energy list if the new image atom
   list differs from the previous one. Also remap the IC table
   and image exclusion lists.

The INVErse and NOINverse options are internal and neither
should be specified under normal circumstances.

.. index:: image; patch
.. _images_patching:

Image Structure File Patching
-----------------------------

This command introduces bonding linkages between primary
atoms and image atoms. This allows the simulation of infinite
(or cyclic) polymers.

:: 

   IMPAtch patch_residue repeat( image_name segid resid ) [SETUp] [WARN]

The patch_residue must be present in the topology file and
the syntax of this patch residue is identical to ordinary patching
(see :ref:`struct_patch`), with the restrictions
that the :chm:`ATOM`, :chm:`DONOr`, and :chm:`ACCEptor` specifications may not be used.
Atom characteristics may not be modified with this command. The donor
and acceptor status of any image atom must match that of the
corresponding primary atom. The patching specifications that are
recognized are; :chm:`BOND`, :chm:`ANGLe`, :chm:`DIHEdral`, :chm:`IMPHi`, and :chm:`IC` (internal
coordinates)

A residue specification is required for each used in the :chm:`PRES`.
These are specified by three names, (1) the image name (for primary
atoms the name :chm:`PRIM` must be used), (2) the segid, and (3) the resid.

The :chm:`SETUp` keyword causes all :chm:`PRES IC` table entries to be added to
the current IC table.

The :chm:`WARN` makes all errors nonfatal and lists errors.

.. index:: image; center
.. _images_centering:

Image Centering
---------------

There is a set of commands that allow for the centering
of selected part of the PSF during an image update. This is primarily
designed for solvent, but may be used in many ways.

::

   IMAGE  { FIXEd        }  [ XCEN real ] [ YCEN real ] [ ZCEN real ]
          { BYSEgments   }                                 atom-selection
          { BYREsidues   }
          { BYGRoups     }
          { BYAToms      }

During dynamics, a particular water may become far from
the rest of the primary structure. The centering features allows one of its
image (the one closest to the primary space) to become the primary water.

It is also useful when setting up a crystal calculation. With a single
update, the "best" image choice of all solvent molecules may be made.
One example of this is the netropsin crystal where one of the published
sulfate groups is quite far from the primary netropsin. This command is
required for a pure solvent simulation where solvent can freely diffuse.

The execution of this command only sets up data used during the image
update. There is only one value each for :chm:`XCEN`, :chm:`YCEN`, and :chm:`ZCEN`. If these
values are not specified in any :chm:`IMAGE` command, then they are not modified
(default 0.0).

For each atom, there is a flag specifying the manner of image
centering to be used. Each invocation of the IMAGE command may modify these
flags. The default is :chm:`FIXEd` (don't center this atom). The :chm:`BYSEgment` option
will center an entire segment as a group (providing it has no FIXED atoms).
The remaining commands will allow certain other groups of atoms to be
centered as a group. It wouldn't work well if only one part of a molecule
was centered (there is no checking for this!).

The command;

::

   IMAGE FIXED SELE ALL END   - will turn off all centering
   IMAGE BYRES SELE RESNAME ST2 END - will allow centering of all ST2's
   IMAGE BYATOM SELE ALL END - will not work if there are any bonds

.. warning::
   Image centering should be set after the PSF is completed.
   Any modification to the PSF with image centering active will
   nullify the image centering function and a warning message is issued. 
   The IMAGE command will need to be reissued if image centering is desired.


.. index:: image; operation
.. _image_operation:

Image Operation
---------------

The IMAGE routines in CHARMM can be classified into five sections.

These catagories are :

* Set up images -IMREAD, REIMAG, INIMAG, IMPATC, IMATOM, IMSPEC
* Update image arrays - UPIMAG, IMCENT, MKIMAT, IMMAP, MKIMNB
* Set up energy lists - IMHBON, NEWHBL, IMHBFX, NBONDM
* Compute image energy - EIMAGE, TRANSO, TRANSI
* Print out - IMWRIT, IMPSFW

The first category involves reading the image file (IMREAD) and
setting up the data structure (REIMAG, INIMAG). In the section are also
the routines involving image patching and setting up the centering options.

The second category concerns itself with the selection of
image groups are to be kept. This selection process is repeated each
image update. Also done, is the centering, PSF remapping (if the atom list
has changed), and the generation of the image exclusion lists.

The third category in addition to finding the energy terms codes,
also generates the nonbond and hydrogen bond lists between primary
and image atoms.

The fourth category is concerned with the computation of energy
terms. For the actual computation of energy, standard routines are used
(ENBOND, EHBOND, ENST2) with a modified calling sequence. The procedure
used is:

1. Compute coordinates for all image atoms
2. Set up arrays for self energy terms (atom with its own image)
3. Compute self terms, divide energy by 2, zero out image forces
4. Compute remaining terms including forces on image atoms
5. Transform forces on image atoms back into the primary space

Using a procedure where the forces on image atoms is kept, allows for
a substantial reduction in the number of necessary image atoms. This
results in the necessity that all transformations have an inverse.
This procedure has the drawback that the self energy terms must be
treated specially and that all Hbonds between image and primary atoms
must be computed and then trimmed of any repeats.

Since there is no treatment of the second derivative of the
energy for image atoms, The procedures involving Newton-Raphson
minimizations and vibrational analysis should be avoided (see :doc:`energy`).).

.. index:: image; bound
.. _images_mipb:

Simple periodic boundaries
--------------------------

This code was developed by William A. Shirley and Charles L.
Brooks, III, Department of Molecular Biology, The Scripps Research
Institute, during the spring of 1995.

Its purpose is to speed up calculations which use the periodic
boundry conditions for A TRUNCATED OCTAHEDRAL BOX; A RHOMBIC DODECAHEDRAL
BOX; A TWO AND THREE DIMENSIONAL RHOMBOIDAL BOX.

Reference:

M. P. Allen and D. J. Tildesley, Computer Simulation of Liquids, Ch. 1.

::

   BOUNd { CUBOUNdary } { BOXL <real>                            }  CUTNB <real>
         { TOBOUNdary } { XSIZe <real> YSIZe <real> ZSIZe <real> }
         { RDBOUNdary }
         { RHBOUNdary }
 
   where
           TOBOUN =        TruncatedOctahedralBOUNd
           RDBOUN =        RhombicDodecahedralBOUNd
           RHBOUN = Two dimensional RHomboidalBOUNd
           CUBOUN =                      CUbicBOUNd

and ``BOXL <real>`` (or :chm:`XSIZe`) and ``CUTNB <real>`` are required.
:chm:``XSIZe``, :chm:``YSIZe``, and :chm:``ZSIZe`` are the edgelengths in the x, y, and z-directions.
If :chm:`YSIZe` is not specified it is assumed to equal :chm:`XSIZe` (or :chm:`BOXL`), and
if :chm:`ZSIZe` is not specified it is assumed to equal :chm:`YSIZe`.

.. warning::
   Image centering and definition of the image transformations
   through the READ IMAGE and IMAGE centering must be used in conjunction
   with this command.  Moving the commands to be available in other fast
   implementations, i.e., in vector codes, will follow shortly.

The basic code is implemented following the outline:
 
1. Parse the control variables based upon the keyword :chm:`BOUNd`, which set
   the control flags to be set.  These flags turn off the existing
   image generation code during run time.  Control is passed back to
   the CHARMM level, and the program continues.

2. Do not generate image atoms to add to the PSF.  When :chm:`UPIMAGE` 
   is called, skip :chm:`MKIMAT`.  Set the completion flag in :chm:`NBONDM`.
 
3. During the generation of the nonbonded list use the MI
   in :chm:`NBONDG` to get correct pairs onto the list.

4. During the calculation of the nonbonded energies in 
   :chm:`ENBFAST`, adjust the distances suing MI to get the correct
   VDW and Elec energies.
