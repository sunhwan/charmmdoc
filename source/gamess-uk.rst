.. py:module:: gamessuk

=======================================================================================
Combined Quantum Mechanical and Molecular Mechanics Method Based on GAMESS-UK in CHARMM
=======================================================================================

By Paul Sherwood (p.sherwood@dl.ac.uk)

based on the GAMESS(US) interface from Milan Hodoscek (milan@par10.mgsl.dcrt.nih.gov,milan@kihp6.ki.si)

Ab initio program GAMESS-UK (General Atomic and Molecular Electronic
Structure System, UK version) is connected to CHARMM program in a QM/MM
method.  This method is based on the interface to the GAMESS (US version),
the latter being an extension of the QUANTUM code which is
described in J. Comp. Chem., Vol. 11, No. 6, 700-733 (1990).

.. _gamessuk_description:

The GAMESS QM potential is initialized with the GAMEss command
--------------------------------------------------------------

::

   GAMEss   [REMOve] [EXGRoup] [QINPut] [BLURred RECAll [INT]] (atom selection)

   REMOve:  Classical energies within QM atoms are removed.

   EXGRoup: QM/MM Electrostatics for link host groups removed.

   QINPut:  Charges are taken from PSF for the QM atoms. Charges
            may be non integer numbers. Use this with the REMOve!

   BLURred: MM charges are scaled by a gaussian function (equivalent to ECP)
            Width of the gaussian function is specified in WMAIN array 
            (usually by SCALar command)
            The value for charge is taken from PSF. Some values of WMAIN have
            special meaning: 

            WMAIN.GT.999.0 ignore this atom from the QM/MM interaction
            WMAIN.EQ.  0.0 treat this atom as point charge in the QM/MM potential

   RECAll:  Use the RECAll array (as specified in scalar.doc) to set BLUR
            widths instead of the main WMAIN array. This is necessary when
            using Gaussian BLURred MM charges with the Replica Path or NEB
            methods as these make use of the WMAIN array. See QM-MM_DGMM.inp
            in the test directory for an example.

            The atoms in selection will be treated as QM atoms.

           Link atom may be added between an QM and MM atoms with the
   following command:


   ADDLinkatom  link-atom-name  QM-atom-spec  MM-atom-spec

         link-atom-name ::= a four character descriptor starting with QQ.

         atom-spec::= {residue-number atom-name}
                      { segid  resid atom-name }
                      { BYNUm  atom-number     }

           When using link atoms to break a bond between QM and MM
   regions bond and angle parameters have to be added to parameter file
   or better use READ PARAm APPEnd command.

           If define is used for selection of QM region put it after all
   ADDLink commands so the numbers of atoms in the selections are not
   changed. Link atoms are always selected as QM atoms.

.. _gamessuk_usage:

Usage
-----

CHARMM input scripts are the same as before except the addition of
ENVIronment commands and the GAMEss command itself. GAMESS-UK commands are in a
separate file call gamess.in, (or with an alternative name indicated by the
"gamess.in" environment variable. The GAMESS-UK input file has the same structure
as it would have for a normal GAMESS-UK run, except that the specification of
the geometry is omitted.

Names of the files for GAMESS-UK are specefied with environment
variables as follows. It is essential to provide a routing for ed3 to 
ensure it is available to hold information between GAMESS-UK calls,
other file specifications are optional.

::

     use ENVIronment command inside CHARMM
     
     envi "ed2"   "/scratch/user/test.ed2" ! quotes needed for lowercase names
     envi "ed3"   "/scratch/user/test.ed3"
     
     or use (t)csh
     
     setenv ed2     /scratch/user/test.ed2
     setenv ed3     /scratch/user/test.ed3
     
     or ksh,sh,bash
     
     export ed2=test.ed2
     export ed3=test.ed3

     or within GAMESS-UK, use the file predirective
     file ed3 /scratch/user/test.ed3
     file ed2 /scratch/user/test.ed2

You can use the "gamess.out" environment variable to 
control the routing of the GAMESS-UK output, or you can define
it as stdout as follows (csh version):

::

	setenv gamess.out stdout

in which case the GAMESS-UK output will be mixed with the charmm
output.  (Note these don't seem to work with the bash shell, as
the export command doesn't accept variable names containing a
period (.), we will have to change this part of the code.


GAMESS-UK input file directives
-------------------------------

The GAMESS-UK data is provided in a separate file, which follows GAMESS-UK
format except that there is no coordinates section.

In addition to those directives needed to switch on the reuiqred energy
expresion, you should

1) include "runtype gradient" to force the computation of both energy
   and forces (unless you are sure you are not going to invoke any calculation
   that requires a gradient from your script).

2) It is advised that the GAMESS-UK directives
   
   ::
   
	   noprint dist anal

   be included as these diagnostic calculations
   don't contribute to the charmm job but use a lot of memory
   when there are a lot of classical atoms.

3) Make sure the gamess input contains a generous time
   card, since the GAMESS calculation will be skipped if
   it thinks it has run out of time.

4) If you see that the quantum part of the energy goes to
   zero, it may reflect the timeout condition above, or some
   other non-fatal problem in GAMESS-UK. Check the GAMESS-UK log
   file.


5) The directive "chm" may be added to set CHARMM-specific options
   as follows:
   
      ===========    ===============================================================
      chm noatom     request that GAMESS-UK output items that list all
                     the atoms be suppressed This is important for macromolecular
                     systems.
                  
      chm append     request that all outputs be concatenated (the default
                     is that you will only save the last one)
                  
      chm offset     provide an energy offset to be added to all the QM
                     energies. This can help ensure the values print within
                     the fields expected in CHARMM.
                  
      chm debug      diagnostic print (not recommended unless developing)
      ===========    ===============================================================


Example:
--------

GAMESS commands have to be in a separate file. Example for the GAMESS input
follows:

::

   core 5000000
   chm append
   chm noatompr
   chm offset 100
   title
   qm region for charmm 
   charge 0
   adapt off
   nosym
   noprint distance analysis
   basis 6-31g
   scftype rhf
   runtype gradient
   vectors atoms
   enter 1

The above is for 6-31G calculation of any neutral molecule.  

.. note::
   For more examples look at test/c28test/cquantumtest/

For complete information about GAMESS input see the CFS web site
http://www.dl.ac.uk/CFS.

For further information and updates on CHARMM/GAMESS-UK interface
see http://www.cse.clrc.ac.uk/qcg/chmguk

.. _gamessuk_installtion:

Installation
------------

Installation itself cannot be fully automated yet so one has to
follow this procedure (if there are any problems ask p.sherwood@dl.ac.uk):

1. Unpack the GAMESS-UK distribution as a subdirectory of gukint:
   source/gukint/GAMESS-UK

2. install.com <machine-type> <size-keyword> U <other-options>

The build procedure works by executing a configuration script within
the GAMESS-UK source tree, (GAMESS-UK/utilities/charmm_configure).
Assuming GAMESS-UK has not already been ported to the target platform, it
is this file that will generally need modification on plaforms
for which the CHARMM/GAMESS-UK interface has not been tested.

The following is a summary of the status (c28 release)

=============  ==============  ==========  =======
Architecture   CHARMM host     Parallel    Status
               keyword         Options
=============  ==============  ==========  =======
SGI R4400      sgi             NA          OK
Pentium/Linux  gnu             NA          OK
"              gnu             mpich       OK
Compaq         alpha           NA          OK
=============  ==============  ==========  =======

Porting Notes
-------------

It is necessary to ensure that charmm_configure processes
the GAMESS-UK Makefiles with a valid set of keywords, the most
important on being the machine type. Unfortunately there isn't
a one-to-one mapping between CHARMM host types and GAMESS-UK
machine types which complicates the charmm_configure script.

Similarly, changes to charmm_configure may be needed to
request the required GAMESS-UK configuration options for
a parallel build. Probably the best bet for a simple 
parallel GAMESS-UK/CHARMM code is to select MPI with
static load balancing options, for which the "mpi"
keyword needs to be passed to configure. 

On some platforms the Global Array port of GAMESS-UK 
can be used, but this is not supported yet by the standard
distribution. See the web site  http://www.cse.clrc.ac.uk/Activity/CHMGMS
for more details of the current status.

NB The GAMESS-UK distribution can only support a single architecture
(there are no architecture dependent directories). When moving the
code from one platform to another, be sure to clear out the
object and library files

::

    % cd source/gukint/GAMESS-UK/m4
    % make clean

When building the parallel code the additional, manual steps
will be needed

For  the MPI code

- in install.com, set the environment variables 

  * MPI_LIB      - the directory holding libmpi.a (or similar)
  * MPI_INCLUDE  - the directory holding mpif.h (etc)

  you may also need to edit GAMESS-UK/m4/Makefile.in as these 
  directories are specified for most platforms.

- Some changes may be needed to build/UNX/Makefile_<arch> to support
  loading with the parallel libraries  (e.g. using MPILD)


.. _gamessuk_status:

GAMESS-UK/CHARMM interface status (July 2000)
---------------------------------------------

- Parallel version is functional, but for most platforms it 
  will require changes to install.com Makefile_$chmhost 
  and/or charmm_configure to activate

- All CHARMM testcases are still OK when CHARMM is compiled
  with GAMESS-UK inside.

- GAMESS, GAMESSUK, CADPAC and QUANTUM keywords cannot coexist
  in pref.dat

- GAMESS-UK recognizes atoms by their masses as specified in the 
  RTF file
