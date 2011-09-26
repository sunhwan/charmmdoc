.. py:module:: charmm_gen

=====================
charmm_gen.com script
=====================

The script  charmm_gen.com  was designed at NIH for easy maintenance of
multiple executables in an active research environment.  Multiple versions
versions can be derived from the same source code, incorporating different
features and maximum atom limits.  It is assumed that install.com has already
been run, and any porting or compiling issues resolved before charmm_gen.com
is used.  In fact, charmm_gen.com simply calls install.com after doing a
little creative copying and renaming.

The script is interactive; it asks a few questions, does a lot of checking,
and then proceeds to make up to nine different versions in one operation with
no further human intervention required.  A "test" or development version can
also be prepared, and is in fact the "path of least resistance", i.e. the
accepting of all the defaults to each prompt.

Since simply starting up a LARGE version of CHARMM with most of the available
feature sets can easily require 100 Mbyte of memory, we recognized the
need to have multiple executables available.  Our choice was to create 3
principal versions: "full", with most major modules included; "lite", a
version without most of the high memory usage or rarely used modules; and
"am1", which adds the QUNATUM QM/MM code and few other features to the "lite"
feature set.  Each is available in 3 sizes, small, medium, and large.

We also use a "cover" script in /usr/local/bin to run CHARMM, after parsing
feature set and size keywords, and stripping them from the command line.  An
example is included at the end of this description.

Currently, ten different sets of object libraries are maintained as well;
this does require a bit of disk space, but allows rapid re-building of all
versions when bugfixes are made.

.. _charmm_gen_configuration:

Configuration
-------------

To use charmm_gen.com, the following additional files are *required* in
build/mach, where mach = hpux in this case:

::

   Makefile_hpux.gamess
   Makefile_hpux.nogamess
   Makefile_hpux.test
   Makefile_hpux.test.list
   pref.dat.large.am1
   pref.dat.large.full
   pref.dat.large.lite
   pref.dat.medium.am1
   pref.dat.medium.full
   pref.dat.medium.lite
   pref.dat.small.am1
   pref.dat.small.full
   pref.dat.small.lite
   pref.dat.test

Makefile_hpux.gamess defines additional tools, directories, and
libraries needed to compile the GAMESS code for QM/MM calculations,
while the .nogamess version is a typical generic version of
Makefile_hpux.  Makefile_hpux.test is configured for rapid re-compiling
of CHARMM during development, while Makefile_hpux.test.list produces
cross-referenced source listings by changing the definition of the
variable FC at the top of the makefile.  The remaining files represent
the ten different possible versions; at NIH, the keywords in
pref.dat.test are usually the same as pref.dat.medium.full, but that
doesn't have to be the case.

The following listing shows the pref.dat keywords we chose for the 2
different feature sets at NIH:

* Feature set "am1"

::

   HPUX        ! machine type
   UNIX
   PARALLEL    ! multiple processors/workstations
   PARAFULL    ! req'd for parallel
   SYNCHRON    ! req'd for parallel
   SOCKET      ! req'd for parallel
   MEDIUM      ! size directive          = 25120 atom limit
   SCALAR      ! machine characteristics = default for scalar machines
   VECTOR      ! feature directive *     = Vectorized routines
   PARVECT     ! Parallel vector code (multi processor vector machines)
   CRAYVEC     ! Fast vector code (standard vector code)
   SAVEFCM     ! Include all SAVE statements
   PUTFCM
   FCMDIR=fcm
   XDISPLAY
   ASPENER     ! feature directive *     = Atomic Solvation Parameter energy term
   MOLVIB      ! feature directive       = MOLVIB vibrational analysis code
   NIH         ! feature directive *     = NIH default specs code
   OLDDYN      ! feature directive       = Old dynamics integrator
   PERT        ! feature directive *     = NIH free energy code
   QUANTUM     ! feature directiver      = include AM1 semi-empirical code
   REPLICA     ! feature directive       = Replica code 
   RISM        ! feature directive       = RISM solvation code
   RXNCOR      ! feature directive *     = RXNCOR code
   TRAVEL      ! feature directive *     = PATH and TRAVEL code
   DIMB        ! feature directive
   FMA         ! feature directive
   ZTBL        ! feature directive
   END         ! end

* Feature set "full"

::

   HPUX        ! machine type
   UNIX
   PARALLEL    ! multiple processors/workstations
   PARAFULL    ! req'd for parallel
   SYNCHRON    ! req'd for parallel
   SOCKET      ! req'd for parallel
   MEDIUM      ! size directive          = 25120 atom limit
   SCALAR      ! machine characteristics = default for scalar machines
   VECTOR      ! feature directive *     = Vectorized routines
   PARVECT     ! Parallel vector code (multi processor vector machines)
   CRAYVEC     ! Fast vector code (standard vector code)
   SAVEFCM     ! Include all SAVE statements
   PUTFCM
   FCMDIR=fcm
   XDISPLAY    ! X11 graphics display
   ASPENER     ! feature directive *     = Atomic Solvation Parameter energy term
   BLOCK       ! feature directive *     = Energy partition and free energy code
   MOLVIB      ! feature directive       = MOLVIB vibrational analysis code
   MTS         ! feature directive       = Multiple time step code
   NIH         ! feature directive *     = NIH default specs code
   OLDDYN      ! feature directive       = Old dynamics integrator
   PERT        ! feature directive *     = NIH free energy code
   GAMESS      ! GAMESS ab initio interface for QM/MM
   REPLICA     ! feature directive       = Replica code 
   RISM        ! feature directive       = RISM solvation code
   RXNCOR      ! feature directive *     = RXNCOR code
   TNPACK      ! Truncated Newton
   TRAVEL      ! feature directive *     = PATH and TRAVEL code
   TSM         ! feature directive       = TSM and ICPERT code
   DIMB        ! feature directive
   FMA         ! feature directive
   FOURD       ! feature directive
   PRIMSH      ! feature directive
   PBOUND      ! simple Periodic BOUNDary (min image)
   SHAPES      ! feature directive       = SHAPE descriptor code
   ZTBL        ! feature directive
   END         ! end

* Feature set "lite"

::

   HPUX        ! machine type
   UNIX
   PARALLEL    ! multiple processors/workstations
   PARAFULL    ! req'd for parallel
   SYNCHRON    ! req'd for parallel
   SOCKET      ! req'd for parallel
   MEDIUM      ! size directive          = 25120 atom limit
   SCALAR      ! machine characteristics = default for scalar machines
   VECTOR      ! feature directive *     = Vectorized routines
   PARVECT     ! Parallel vector code (multi processor vector machines)
   CRAYVEC     ! Fast vector code (standard vector code)
   SAVEFCM     ! Include all SAVE statements
   PUTFCM
   FCMDIR=fcm
   XDISPLAY
   ASPENER     ! feature directive *     = Atomic Solvation Parameter energy term
   MOLVIB      ! feature directive       = MOLVIB vibrational analysis code
   NIH         ! feature directive *     = NIH default specs code
   PERT        ! feature directive *     = NIH free energy code
   REPLICA     ! feature directive       = Replica code 
   RXNCOR      ! feature directive *     = RXNCOR code
   TRAVEL      ! feature directive *     = PATH and TRAVEL code
   END         ! end



.. _charmm_gen_cover:

Cover
-----

Finally, to make live easy for the end users, we use the following
script to run CHARMM on a routine basis:

::

   #! /bin/csh
   # INITIALIZE VARIABLES
   set n = $#argv
   setenv HOST `hostname | cut -d. -f1`
   set chmsiz = small
   set i = 1
   set cleanup = date
   set chmopt = lite
   # CHECK FOR OPTIONAL KEYWORDS
   while ( $i <= $n )
    switch ( $argv[$i] )
     case small:
      set chmsiz = small
      breaksw
     case large:
      set chmsiz = large
      breaksw
     case medium:
      set chmsiz = medium
      breaksw
     case test:
      set chmopt = test
      breaksw
     case lite:
      set chmopt = lite
      breaksw
     case full:
      set chmopt = full
      breaksw
    endsw
    @ i = $i + 1
   end
   # STRIP KEYWORDS FROM ARGUMENT STRING
   set t = `echo $* | sed -e 's/small//' -e 's/medium//' -e 's/test//' \
      -e 's/full//' -e 's/large//' -e 's/lite//' -e 's/am1//'`
   # CHECK FOR DESIGNATED PARALLEL HOSTS
   switch ( $HOST )
    case par0:
     set cleanup = 'qpara_clean par0 bypass'
     setenv NODE0 par0f
     setenv NODE1 par1f
     setenv NODE2 par2f
     setenv NODE3 par3f
     echo "Parallel; $NODE0 $NODE1 $NODE2 $NODE3"
     breaksw
    case par11:
     set cleanup = 'qpara_clean par11 bypass'
     setenv NODE0 par11f
     setenv NODE1 par12f
     setenv NODE2 par13f
     setenv NODE3 par14f
     echo "Parallel; $NODE0 $NODE1 $NODE2 $NODE3"
     breaksw
    default:
     echo "Single processor; $HOST"
     breaksw
   endsw
   # ECHO WORKING DIRECTORY AND CHARMM VERSION W. TIMESTAMP
   if ( $?PWD ) then
    echo $PWD
   else
    echo $cwd
    echo "Warning: env var PWD not defined; required for parallel CHARMM"
   endif
   # SET THE VERSION TO BE RUN
   if ( $chmopt == test ) then
    set exe = $chmopt
   else
    set exe = $chmsiz.$chmopt
   endif
   # VERIFY THE ACTUAL EXECUTABLE; RUN AT REDUCED PRIORITY
   ls -o ~charmm/c24n4/exec/hpux/charmm.$exe | cut -c33-
   if { /bin/nice -5 ~charmm/c24n4/exec/hpux/charmm.$exe $t } then
    echo ''
    $cleanup
   else
    echo '(charmm) ABNORMAL EXIT'
    $cleanup
    exit(1)
   endif


