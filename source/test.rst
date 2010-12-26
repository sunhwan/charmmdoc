.. py:module:: test

============================================================
Test commands: Commands to test various conditions in CHARMM
============================================================

Syntax of the TEST commands:

::

   TEST FIRSt [TOL real] [STEP real] [UNIT int] [MASS int] [atom-selection]
                 (0.005)    (0.0001)        (6)        (0)

                     [ CRYStal  [ HOMOgeneous ] ]

   TEST SECOnd  [TOL real] [STEP real] [UNIT int] [2xatom-selection]
                  (0.005)    (0.0001)  (OUTU (6)) 

   TEST COORdinates  [COMP]

   TEST CONNectvity [SUBSet atom-selection] [COMMon atom-selection] [PRINt]

   TEST PSF

   TEST PARAmeter TRIGonometry  {  DIHEdral   }
                                {  CDIHedral  }
                                {  VDIHedral  }
                                {  IMPRoper   }

   TEST HEAP

   TEST STACk

   GETHeap integer

   TEST INITialization

   RESET

   TEST NOCOmmunication  { READ   UNIT int  STEP int  MEMO int } 
                         { WRITE  UNIT int  STEP int  MEMO int }
                         { CLOSE                               }

   TEST STAMp LEVEl int

The TEST FIRSt command, tests the first derivative of the
energy by finite differences. It uses the GETE subroutine, so that before
this command is invoked, the UPDAte command must be invoked. Since
two energy evaluations are done for each degree of freedom, an atom
selection should be used for large systems. The TOL keyword may be
used to list only those terms which differ by more than a particular
tolerance factor. The STEP keyword gives the finite difference step
size in angstroms. The MASS integer value determines how mass weighting
will be used with regards to rigid or SHAKE constraints (0= no mass
weighting, 1= mass weighting, -1=ignore water hydrogen masses).
This command is analogous to the second derivative test command under
the vibrational analysis (VIBRAN:: WRITE SECOnd derivatives  CARD
FINIte).

The CRYStal keyword also causes the forces along the crystal
degrees of freedom to be checked.  The HOMOgeneous keyword causes the
causes all of the atoms to be scaled in a homogeneous fashion when
the box size/shape changes.

The TEST SECOnd command, tests the second derivative of the
energy by finite difference of the forces. It is thus assumed that the
first derivative is right. UPDAte should be invoked before use. Selection is 
allowed. The elements of the 2nd derivative matrix computed are the first
selection by the second selection in size. If only one selection is specified
then the the components of the 2nd derivative matrix corresponding to this
selection is constructed. If no selection is given, then the whole 
3Nx(3N+1)/2 second derivative matrix is created by analytic and 
finite difference methods. So don't use a too big system (a few hundred 
atoms depending on your computer). TOL and STEP have the same meaning than 
for TEST FIRSt.

The TEST COORdinates command will test to insure that all
coordinates are in range and defined. The comparison coordinates may be
tested.

The TEST CONNectivity command tests the structure for proper
connectivity between the selected set of atoms. This is to facilitate
finding loose pieces of a molecule after a DELEte ATOM command.  The
algorithm checks to make sure that some connective bonding, bonds,
angles or torsions, exist between atoms in the selection SUBSet and
atoms in the COMMon set selection.  For example, if a sequence of
amino acid residues, say 1-5, were selected as the SUBSet and only residue
5 was selected for the COMMon set the CONNectivity between the two sets
is true.  If, say, residue 10 was selected for the common set then the
test would prove false and the message "One disconnected segment found"
would be printed.

The TEST PSF command tests most of the data in the structure file.
The data is sorted and duplications are found as well.

The TEST PARAmeter command write out the tables used by energy
routines. DIHEdral, CDIHedral and IMPRoper TRIGonometry tables are the
parameters listed in the order they are stocked. Number,
Force-constant, Periodicity, Reference Angle, Cosine and sine are
listed. The VDIHedral option lists the dihedral energy term series as
they are used in the PARVect dihedral energy routines. Term-number,
Parameter-number (cf DIHEdral),  Force-constant, Periodicity,
Reference Angle, Cosine and sine are listed. TEST PARAmeter is mostly
for developing purposes. An unknown option will be ignored.

The TEST HEAP command causes the heap data to be printed
between non-miscellaneous commands. This will continue until
job termination. This makes it easy to see where space allocation
went wrong during program development.  The TEST STACk command performs
the same function for stack space usage.

The GETHeap routine allocates and then frees a block of memory
specified in words.  This command is used to prevent fragmentation of
memory, or to preallocate disk space for a long run that will need
space later.

The TEST INITialization command will initialize just about
everything. Its use is in testing the program when making major
changes or additions.

The RESEt command deletes all atoms and calls the initialization
routine.  It should be called before TEST INIT to reset everything.

The TEST NOCOmmunication performs parallel run reading data
from memory instead. TEST NOCO WRITE command has to be performed on a
single CPU and writes all the necessary data from memory to a file
specified by UNIT keyword. The STEP keyword specifies the number of
steps for which the information will be saved to a file. MEMOry
keyword reserves the memory needed for this operation. TEST NOCO READ
has the same parameters as WRITE, but the complete information is
stored to the memory from the file for the number of steps specified
with the STEP keyword. TEST NOCO CLOSe has to be specified after WRITE
command and before STOP.

The TEST STAMp command outputs a time stamp in microsecond
precision of the following events:

::

        level= 0 - OFF
              -1 - from input script
               1 - report from dynamics
               2 - report from paral1 routines
               4 - report from paral3 routines

.. note::

   It is possible to combine reports, ie level=6 is the same as
   level.eq.4 or level.eq.2

