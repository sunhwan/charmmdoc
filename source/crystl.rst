.. py:module:: crystal

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

::

   CRYStal  [BUILd_crystal] [CUTOff real] [NOPErations int]
   
            [DEFIne xtltype a b c alpha beta gamma]
   
            [FREE]
   
            [PHONon] [NKPOints int] 
                     [KVECtor real real real TO real real real]
   
            [VIBRation]
   
            [READ] [CARD UNIT int]
                   [PHONons UNIT int]
   
            [PRINt]
            [PRINt] [PHONons] [FACT real] [MODE int THRU int] 
                                          [KPTS int TO int]
   
            [WRITe] [CARD UNIT int]
                    [PHONons UNIT int]
                    [VIBRations] [MODE int THRU int] [UNIT int]
   
   
   xtltype ::=    { CUBIc          }
                  { TETRagonal     }
                  { ORTHorhombic   }
                  { MONOclinic     }
                  { TRIClinic      }
                  { HEXAgonal      }
                  { RHOMohedral    }
                  { OCTAhedral/trnc}
                  { RHDO           }

   a b c alpha beta gamma ::= (six real numbers)
   
The crystal module is an extension of the image facility
within the CHARMM program.  All crystal commands are invoked by the
keyword CRYStal.  The next word on the command line can be one of the
following :

* Build: builds a crystal.
* Define: defines the lattice type and constants of the crystal to be studied.
* Free: clear the crystal and image facility.
* Phonon: calculates the crystal frequencies for a single value or a range of values of the wave vector, KVEC.
* Print: prints various crystal information.
* Read: reads the crystal image file.
* Vibration: calculates the harmonic crystal frequencies when the wave vector is the zero vector.
* Write: writes out to file various crystal information.

.. _function:

Function
--------

.. _crystal_build:

Crystal Build
^^^^^^^^^^^^^

A crystal of any desired symmetry can be constructed by repeatedly applying a
small number of transformations to an asymmetric collection of atoms (called here
the primary atoms). The transformations include the primitive lattice translations 
*A*, *B* and *C* which are common to all crystals and a set of additional transformations, 
*T*, which determines the space group symmetry.

The Build command will generate, given *T*, a data structure of all those transformations 
which produce images lying within a user-specified cutoff distance of the primary atoms. 
The data structure can then be used by CHARMM to represent the complete crystal of the 
system in subsequent calculations. The symmetry operations, *T*, are read from the lines 
following the ``CRYStal BUILd`` command.

The syntax of the commmand is :

::

  CRYStal BUILd CUTOff <real> NOPErations <int>
  ! <int> lines defining the symmmetry operations.

The ``CUTOff`` parameter is used to determine the images which are included in the 
transformation list. All those images which are within the cutoff distance are included in the list.

.. Note::
   The distance test is done based on the atoms that are currently present and their symmetric representation.
   
To generate a crystal file from a box with a single atom at the center, the cutoff value 
will nee to be larger than the box dimensions.  If the box is filled with water and only 
nearest neighbor cells are desired, then the cutoff distance should be comparable to the 
``CUTIm`` value (see :ref:`Image Updates <image_update>` ) or the ``CUTNB`` value (see 
:ref:`NBONDS Syntax <nbonds_syntax>`). There is no limit to the number of transformations 
included in the lists as they are allocated dynamically, but having too many will slow the 
image update step.

The crystal symmetry operations are input in standard crystallographic notation. The 
identity is assumed to be present so that (X, Y, Z) need not be specified (in fact, 
it is an error to do so). For example, a P1 crystal is defined by the identity operation 
and so the input would be

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

It should be noted that in those cases where the atoms in the asymmetric unit have 
internal symmetry or in which a molecule is sited upon a symmetry point within the unit 
cell not all symmetry transformations for the crystal need to be input. Some will be 
redundant. It is up to the user to check for these cases and modify the input accordingly.

.. _crystl_define:

Crystal Define
^^^^^^^^^^^^^^

The define command defines the crystal-type on which calculations
are to be performed. It is usually the first crystal command that is
specified in any job using the crystal facility.  It has the format:

::

   DEFIne xtlype a  b  c  alpha beta gamma

The input lattice parameters are checked against the lattice-type to
ensure that they are compatible. Nine lattice types are permitted. They
are listed below along with any restrictions on the lattice parameters:

   CUBIc
      | a = b = c and alpha = beta = gamma = 90.0 degrees.
      | (example:  50.0 50.0 50.0 90.0 90.0 90.0 )
      | (volume = a**3)
      | (degrees of freedom = 1)
   
   TETRagonal
      | a = b and alpha = beta = gamma = 90.0 degrees.
      | (example:  50.0 50.0 40.0 90.0 90.0 90.0 )
      | (volume = c*a**2)
      | (degrees of freedom = 2)

   ORTHorhombic
      | alpha = beta = gamma = 90.0 degrees.
      | (example:  50.0 40.0 30.0 90.0 90.0 90.0 )
      | (volume = c*b*a)
      | (degrees of freedom = 3)

   MONOclinic
      | alpha = gamma = 90.0 degrees.
      | (example:  50.0 40.0 30.0 90.0 70.0 90.0 )
      | (volume = c*b*a*sin(beta) )
      | (degrees of freedom = 4)

   TRIClinic
      | no restrictions on a, b, c, alpha, beta or gamma.
      | (example:  50.0 40.0 30.0 60.0 70.0 80.0 )
      | (volume = c*b*a*sqrt(1.0 - cos(alpha)**2 - cos(beta)**2 -
      |     cos(gamma)**2 + 2.0*cos(alpha)*cos(beta)*cos(gamma)) )
      | (degrees of freedom = 6)

   HEXAgonal
      | a = b,  alpha = beta = 90.0 degrees and gamma = 120.0
      | (example:  40.0 40.0 120.0 90.0 90.0 120.0 )
      | (volume = sqrt(0.75)*c*a**2 )
      | (degrees of freedom = 2)

   RHOMbohedral
      | a = b = c ; alpha=beta=gamma<120  (trigonal)
      | (example:  40.0 40.0 40.0 67.0 67.0 67.0 )
      | (volume = a**3*(1.0-cos(alpha))*sqrt(1.0+2.0*cos(alpha)) )
      | (degrees of freedom = 2)

   OCTAhedral (a.k.a truncated octahedron)
      | a = b = c, alpha = beta = gamma = 109.4712206344907  
      | (example:  40.0 40.0 40.0 109.471220634 109.471220634 109.471220634 )
      | (volume = 4*sqrt(3))/9 * a**3 )
      | (truncated cube length = a * sqrt(4/3) )
      | (degrees of freedom = 1)

   RHDO (Rhombic Dodecahedron)
      | a = b = c, alpha = gamma = 60.0 and beta = 90.0
      | (example:  40.0 40.0 40.0 60.0 90.0 60.0 )
      | (volume = sqrt(0.5) * a**3 )
      | (truncated cube length = a * sqrt(2) )
      | (degrees of freedom = 1)

It is up to the user to ensure that the lattice parameters have the
desired values for the system at all times. The values are stored
by the program but, at present, the only way to transmit this information
between jobs is with binary coordinate, trajectory, or restart files.
For example, if the lattice parameters have been changed during a
lattice optimization then the new parameters, which are printed out at
the end of the minimization, must be input at the beginning of
the next CHARMM run, or transferred using the FILE option on coordinate
writing and reading.  Lattice parameters are stored in binary coordinate,
dynamic trajectory, and restart files only.

.. _crystl_phonon:

Crystal Phonon
^^^^^^^^^^^^^^

Phonon calculates the dispersion curves for a crystal. Any value
of the wavevector can be used (although, in practice, each component
of ``KVECector`` is normally limited to the range -0.5 to +0.5). The dynamical
matrix and normal mode eigenvectors determined in the phonon
calculation are complex although the eigenvalues remain real.

The syntax for the command is :

::

   CRYStal PHONon NKPOints <integer> KVECtor <real> <real> <real> TO <real> <real> <real>

``NKPOints`` tells the program the number of points at which the derivative
matrices must be built and diagonalised whilst the  ``KVECtor ... TO ...``
clause determines the values of KVECtor for each calculation. Thus,

::

   KVECtor 0.0 0.0 0.0   TO 0.5 0.5 0.5   NKPOints 3

would solve for the crystal frequencies at the points, KVEC=(0.0,0.0,0.0),
(0.25,0.25,0.25) and (0.5,0.5,0.5). If it is desirable, point calculations
can be carried out by omitting the  To statement and putting  Nkpoints 1.
For single calculations at KVEC=(0.0,0.0,0.0) the :ref:`crystl_vibration` command
is faster.

The eigenvalues and eigenvectors at each value of the wave vector
from the phonon calculation are saved and they can be written out to a
file using the ``Crystal Write Phonon`` command. No analysis facilities
exist within CHARMM for the phonon data structure as the eigenvectors
are complex.

It is to be noted that phonon and vibration calculations can only
be performed on crystals of P1 symmetry. No information about the
symmetry operations is used when generating the dynamical matrix.


.. _crystl_print:

Crystal Print
^^^^^^^^^^^^^

Two options exist with the ``Print`` command. If no keyword is given
then the crystal image file is printed out.

The ``Crystal Print Phonon`` command performs a similar function to the
``Print Normal_Modes`` command in the vibrational analysis facility. Selected
frequencies and eigenvectors for a range of values of the wave vector can
be printed out. The syntax is:

::

   CRYStal PRINt PHONon KPOInts <i> TO <i> MODEs <i> THRU <i> FACTor <f>

The ``Kpoints .. To ..`` clause determines the wave-vectors at which the
modes are to be printed, the ``Modes .. Thru ..`` gives the range of the
eigenvectors and the Factor command gives the scale factor to multiply
each normal mode by.


.. _crystl_read:

Crystal Read
^^^^^^^^^^^^

The :ref:`crystl_read` command reads in a crystal image file. The file
has the same output as produced by the :ref:`crystl_print` or :ref:`crystl_write`
commands.  The command is useful if a crystal image file was produced
using the :ref:`crystl_build` command and saved using the :ref:`crystal_write`
command in a previous job and it is desired to reuse the same
transformation file for analysis or comparison purposes. The command
can also be used to read in limited sets of transformations if
specific crystal interactions need to be investigated. The
transformation file is formatted so the ``Card`` keyword needs to be
specified and the unit number must be given after the ``Unit`` keyword.


.. _crystl_vibration:

Crystal Vibration
^^^^^^^^^^^^^^^^^

For a free molecule with N atoms the dynamical equations have 3N-6
non-zero eigenvalues. This is no longer so for a crystal. If a crystal
is made up of L unit cells each containing Z molecules with N atoms,
the dynamical equations would have a dimension of 3NZL. However, using
the symmetry properties of the lattice it is possible to factor the
equations into L sets each with a dimension of 3NZ and each depending
upon a vector, KVEC, which labels the irreducible representation of the
translation group to which the set belongs. The force constant matrix
is complex. Its form may be found in the references given at the end of
the documentation.

Vibration solves the dynamical equations for the case where the wave-vector
is zero, i.e. when the equations are real. The procedure is invoked by the
:ref:`crystl_vibration` command. The syntax is :

::

   Crystal Vibration


.. _crystl_write:

Crystal Write
^^^^^^^^^^^^^

There are three ``Crystal Write`` options. If no keyword is given the
crystal image file is written out, in card format, to the specified
unit. The ``CARD`` and ``UNIT`` keywords are required.

The ``Crystal Write Phonon`` command writes out the phonons from a
phonon calculation. All the eigenvalues and eigenvectors for all
values of the wavevector that are stored are written automatically.

The ``Crystal Write Vibration`` command writes out the eigenvalues and
eigenvectors from a vibration calculation. The modes to be written are
given by the ``Mode .. Thru ..`` clause. 

All ``Write`` commands require that the Fortran stream number be given
after the Unit keyword and a CHARMM title may be specified on the
following lines. 

The structure of the phonon and vibration files for a crystal may
be found by looking at the routines ``WRITDC`` and ``XFRQW2`` respectively
in the file ``[.IMAGE]XTLFRQ.SRC``. The vibration modes are written
in the same form as a for :doc:`vibran` normal mode file and may be read
in using the appropriate :doc:`vibran` commands. Unfortunately no analysis
facilities exist for complex eigenvectors within CHARMM and so users
will have to write their own if they want to perform phonon
calculations.


.. _crystl_minimization:

Crystal Minimization
^^^^^^^^^^^^^^^^^^^^

It is possible to perform a lattice minimization using the normal
CHARMM :doc:`MINImize <minimiz>` command and the :doc:`ABNR <abnr>` minimizer. Two extra keywords
have been introduced. If none of them is present then a coordinate
minimization is performed as usual. If ``LATTICE`` is specified then
the ``LATTice`` parameters and the atomic coordinates are minimized
together. If ``NOCOoordinates`` is given with the keyword ``LATTice`` then
only the lattice parameters are optimized. Specifying ``NOCOordinates``
by itself is an error.

It should be noted that when the lattice is being optimised the
crystal symmetry is maintained. A cubic crystal will remain cubic, etc.

.. _examples:

Examples
--------

Examples of input may be found in the test directory. All crystal
files are prefixed by the string :file:`xtl_{*}`. All the jobs involve
L-Alanine. Briefly the jobs are:

1. :download:`XTL_ALA1.INP </_downloads/testcases/xtlala1.inp>`

   The crystallographic fractional coordinates are
   read in and converted to real space coordinates
   using the CHARMM ``COORdinate CONVert`` command and
   the experimental values for the lattice parameters.

2. :download:`XTL_ALA2.INP </_downloads/testcases/xtlala1.inp>`

   A crystal image file is generated for the crystal
   using a value of 10.0 Angstroms for the crystal
   cutoff.

3. :download:`XTL_ALA3.INP </_downloads/testcases/xtlala1.inp>`

   A coordinate and lattice minimization are performed
   for the crystal. The crystal image file from the
   previous job is used and the optimised coordinates
   are saved. The main point to note is that before
   using the crystal package for energy calculations
   and other manipulations that involve the image
   non-bond lists an image update must be performed.
   For safety always do an update after building or
   reading in the crystal. Note too that the new,
   optimised lattice parameters are used in the all
   the subsequent input files.

4. :download:`XTL_ALA4.INP </_downloads/testcases/xtlala1.inp>`

   For subsequent calculations a coordinate file that
   contains the coordinates of all atoms (four
   molecules of L-alanine) is generated. A crystal
   image file suitable to do this is read in directly
   from the input stream. It contains 6 transformations
   (not 3 as might be expected) because the CHARMM
   image facility requires that the inverses of all
   transformations be present. The first three are the
   ones needed and the last three are their inverses.
   An update is needed after reading the file to make
   known to the program the coordinates of the atoms
   in the first transformation of all the inverse pairs
   in the image list. The ``Print Coor Image`` file will
   then print out the coordinates of the atoms in the
   original asymmetric unit and the first three of the
   images. If the coordinates of the atoms in all the
   images are required then the keyword ``NOINV`` in the
   ``UPDATE`` command must be used (check :doc:`IMAGE.DOC <image>`).

5. :download:`XTL_ALA5.INP </_downloads/testcases/xtlala2.inp>`

   The same job as the second except that the crystal
   is generated for a whole unit cell (i.e. the system
   generated in the fourth job). The same value of the
   crystal cutoff is used. An energy is calculated too.
   The energy and its RMS coordinate derivative should
   be exactly four times (apart from a small round-off
   error) the value obtained for an energy calculation
   on a single asymmetric unit with the same lattice
   parameters and crystal cutoff (see job 3).

6. :download:`XTL_ALA6.INP </_downloads/testcases/xtlala2.inp>`

   Peform a crystal vibration and phonon calculation
   for the optimised structure of the L-alanine
   crystal. The vibrational and phonon modes are
   written out to files and components of the first 24
   phonon normal modes for the three values of the
   wavevector that were calculated are printed. To
   do the same for the vibrations it would be necessary
   to use the appropriate :doc:`VIBRAN <vibran>` commands in another
   job.

.. _implementation:

Implementation
--------------

Background and Implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Crystal options and their commands were described above. The present
section discusses relevant background material and briefly reviews the
methods used in the implementation. Some technical points are also made.

The crystal option is an extension to the CHARMM program.  The source
code is in the directory ``[.IMAGE]`` whilst the crystal data structure is in
the file :file:`IMAGE.FCM`. Two additional source code files have been added -
:file:`CRYSTL.SRC` and :file:`XTLFRQ.SRC`. Small modifications have been made to the
files :file:`ENERGY.SRC` and :file:`EIMAGE.SRC`.

CHARMM Images and the Crystal Image Data Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As outlined above a crystal structure can be specified entirely
by the action of the primitive translations A, B and C, and a small set of
transformations, *T* (which themselves are functions of A, B and C), on an
asymmetric group of atoms. In CHARMM the calculation of the energy assumes
that there exists a cutoff distance beyond which all interactions between
particles are neglected so that when performing calculations on
supposedly infinite crystals only a limited portion of that crystal, i.e.
that portion containing those atoms within the cutoff distance of the
primary atoms, need be considered.

The CHARMM image option, of course, already enables the energies of
crystals to be calculated but the input required to use it to do so is
cumbersome and time consuming. It is a great simplification to include an
extra data structure that defines the crystal in terms of A, B and C and
*T*.

There are a number of advantages:

1. A crystal is regular so that its generation can be automated. All that
   needs to be done is to systematically transform the primary atoms by
   one of the set *T* and a linear combination of A, B and C.
   The result is obviously best stored in terms of A, B, and C
   rather than as absolute numerical values of the transformations.

2. It is essential to define a CHARMM crystal by A, B and C and *T* if the
   lattice parameters a, b, c, alpha, beta and gamma are to be varied
   because the coordinates of all the image atoms within the crystal will
   change during successive cycles of the optimization as a, b, c, alpha,
   beta and gamma themselves change.

3. When constructing the dynamical matrix for a non-zero wave-vector it is
   necessary to know the unit cell to which a particular atom belongs in
   order to evaluate the exponential factor in the expression.

Although the crystal data structure and the values of the lattice
parameters define the crystal the individual transformations have to be
worked out explicitly in order to determine energies, harmonic frequencies
and so on. In the present version of the program the :doc:`IMAGE <image>` facility is
used, so that a new set of :doc:`IMAGE <image>` transformations are calculated from the
crystal data structure as soon as a crystal is built or every time the
lattice parameters are changed. The use of the :doc:`IMAGE <image>` facility means that
the number of transformations that can be used is determined by the
dimension of the :doc:`IMAGE <image>` arrays (``MAXTRN`` in ``DIMENS.FCM``).


Crystal and Image Patching
^^^^^^^^^^^^^^^^^^^^^^^^^^

Crystal image patching is unavailable in the present version of the
program so that bonds between images are not permitted. Similarly
hydrogen-bond interactions described by an explicit hydrogen-bond function
are also forbidden. The only forces that can be calculated between primary
and image atoms are non-bonded ones.


The Lattice Coordinate System
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning::

   If your system is not properly rotated, there will usually be
   bad contacts.  If you have bad contacts, check the alignment.

The convention used by CHARMM for orientating the crystal in real space involves
the use of a symmetric transformation (h) matrix.  For non-orthorhombic systems,
these coordinates are different (rotated) from the aligned conventioned used by
PDB and others.  The conversion is performed by the ``COOR CONVert`` command,
see :ref:`corman_syntax <Corman Syntax>`.


The Structure of the Crystal File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The crystal file is divided into three parts.

   A standard CHARMM title.

   A symmetry operation declaration section headed by the word Symmetry
   and terminated by an End. The transformations are written in the same
   way as for the :ref:`crystl_build` command except that the identity
   transformation has to be explicitly listed.

   An image section headed by Images and terminated by an End. Here the
   images are defined in terms of the symmetry transformations and the
   lattice translations A, B and C. The comment line shows the column
   labeling.

Sometimes it is useful to write one's own crystal files without recourse
to the :ref:`crystl_build` option. In this case the symmetry and image blocks
can be put in any order (although only one of each is allowed) and there
is no restriction on the positioning of blank and comment lines.

Two examples of a crystal file are:

::

   * Crystal file for a P1bar crystal.
   *

   Symmetry
   (X,Y,Z)
   (-X,-Y,-Z)
   End

   Images
   ! Operation       a    b    c
             1       0    0   -1
             1       0    0    1
             2       0    0    0
   End


:: 

   * Crystal file for a P212121 crystal.
   *

  Symmetry
  (X,Y,Z)
  (X+1/2,-Y+1/2,-Z)
  (-X,Y+1/2,-Z+1/2)
  (-X+1/2,-Y,Z+1/2)
  End

  Images
  ! Operation       a    b    c
            2       0    0    0
            3       0    0    0
            4       0    0    0
            2      -1    0    0
            3       0   -1    0
            4       0    0   -1
  End
  

Second Derivative Calculations and the Use of Symmetry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Consider a crystal with a unit cell in which there is more than one
asymmetric unit (i.e. all space groups other than P1). The dynamical
matrix then takes a blocked form, with Z**2 blocks if Z is the number
of asymmetric units. Each block is of dimension 3N x 3N and contains
the sum over all unit cells of the second derivative interaction
elements between the *M*th and Nth asymmetric units. It is possible to
calculate only the *Z* blocks (11), (12), ..., (1M), ..., (1Z) and then
transform them to produce the full matrix. In the present program,
however, it is necessary to perform vibration calculations on entire
unit cells.

It should be emphasized that while this symmetry transformation can be
used for calculations of the normal mode eigenvectors and frequencies
for the zero wavevector it does not hold at other values for all additional
values. Therefore, simple symmetry arguments such as these do not hold
for phonon calculations.

Symmetry can also be used to block the dynamical matrix into several
smaller matrices each corresponding to a different symmetry species,
thereby greatly reducing the time needed for diagonalisation and
automatically helping to identify the normal modes. Symmetry blocking
is not coded at the moment.


References
^^^^^^^^^^

"Lattice Dynamics of Molecular Crystals", Lecture Notes in Chemistry 26,
S.Califano, V.Schettino and N.Neto (1981), Springer-Verlag, Berlin,
Heidelberg and New York. A comprehensive monograph with good sections
on the theory of lattice vibrations and normal mode symmetries.

A.Warshel and S.Lifson, J.Chem.Phys. (1970), 53, 582. The original CFF
paper on crystal calculations. It describes the theory behind crystal
optimisations and vibrational calculations.

E.Huler and A.Warshel, Acta Cryst. (1974), B30, 1822. An extension of
the work in reference 2.

"Infrared and Raman Spectra of Crystals", G.Turrell (1972), Academic
Press, London and New York. A nice clear introduction to the subject.
