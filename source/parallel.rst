.. py:module:: parallel

=================================
Parallel Implementation of CHARMM
=================================

CHARMM has been modified to allow computationally intensive simulations
to be run on multi-machines using a replicated data model.  This
version, though employing a full communication scheme, uses an efficient
divide-and-conquer algorithm for global sums and broadcasts.

Curently the following hardware platforms are supported:

  1. Cray T3D/T3E
  2. Cray C90, J90
  3. SGI Power Challenge
  4. Convex SPP-1000 Exemplar
  5. Intel iPSC/860 gamma
  6. Intel Delta machine
  7. Intel Paragon machine
 8. Thinking Machines CM-5
 9. IBM SP1/SP2 machines
10. Parallel Virtual Machine (PVM)
11. Workstation clusters (SOCKET)
12. Alpha Servers (SMP machines, PVMC)
 13. TERRA 2000
 15. Convex SPP-2000
 17. LoBoS (any Beowulf)
14. HP SMP machines
16. SGI Origin
18. IBM Power4 using GNU/Linux system


::

  * Menu:


  * Syntax::        Syntax for PARAllel command
  * Installation::  Installing CHARMM on parallel systems
  * Running::       Running CHARMM on parallel systems
  * PARAllel::      Command PARAllel controls parallel communication
  * Status::        Parallel Code Status (as of September 1998)
  * Using PVM::     Parallel Code implemented with PVM
  * Implementation:: Description of implementation of parallel code

  
  File: parallel.doc,  Node: Syntax,  Next: Installation,  Prev: Top,  Up: Top


  PARAllel command parser for controlling parallel execution
  Syntax:

  PARAllel CONCurrent <int> ...

      CONCurrent <int>   specify how many concurrent jobs
                         to run in the system

  PARAllel FIFO <int>    specify FIFO scheduler in LoBoS with
                         static priority <int>

  PARAllel BUFF <int>    specify buffer size for send/receive
                         calls. <int> is in REAL*8 units

  PARAllel INFO          Prints the hostname information for each process
                         Also fills arrays PARHOST, PARHLEN in parallel.fcm

  
  File: parallel.doc,  Node: Installation,  Next: Running,  Prev: Syntax,  Up: Top


  For support of many parallel comunication libraries the CMPI keyword
  was added. In order to get the old communication routines always
  specify CMPI otherwise MPI is the default choice (see recommended
  keyword combination for each specific platform). On some platforms
  recommended preflx directives prepare the code which does the
  communication much faster, eg on 128 nodes T3E CMPI is 4 times faster
  than MPI. For spatial decomposition method PARAFULL or PARASCAL must
  be replaced by SPACDEC pref.dat keyword

  This is a complete list of supported combinations for message passing
  libraries implemented in the parallel CHARMM

  Combinations of pref.dat keywords for MPI library (can be specified on
  any platform that support MPI):

  1. < no extra keywords > (Calls to MPI collective routines)
  2. CMPI MPI (non-blocking cube topology using send/receive from MPI)
  3. CMPI MPI GENCOMM (non-blocking ring topology, MPI send/receive)
  4. CMPI MPI SYNCHRON (blocking cube topology, MPI send/receive)
  5. CMPI MPI GENCOMM SYNCHRON (blocking ring topology, MPI send/receive)

  NOTE: using GENCOMM is slower then without it. GENCOMM is mostly used
        for QM/MM replica path method where the scaling is almost
        perfect anyway.

  Additionally there is a pref.dat keyword PARINFNTY, which simulates
  the infinitively fast network. In other words there is no communication
  involved during the dynamics after the parallel run is setup. Needles
  to say the results of such calculations are meaningless. Also in order
  to get a few 1000 of steps of dynamics one need to use very small
  timesteps, eg 0.000001. The purpose of this keyword is for testing
  CHARMM performance and also to compare a variety of parallel system
  setups. It works in combination with CMPI keyword. For example one
  should specify CMPI MPI PARAFULL PARINFNTY.



  Native library options

  6. CMPI DELTA (for Intel Paragon)
  7. CMPI IBMSP (for IBM SP2)
  8. TERRA (for TERRA 2000)
  9. CMPI CM5 (For CM5)
  10. CSPP (Convex version of MPI)

  Workstation clusters using SOCKET

  11. CMPI SOCKET SYNCRON (blocking cube topology)
  12. CMPI SOCKET SYNCRON GENCOMM (blocking ring topology)

  PVM library

  13. CMPI PVMC SYNCHRON (blocking cube, PVM send/receive)
  14. CMPI PVMC GENCOMM SYNCHRON (blocking ring, PVM send/receive)


  Combination 1., 8. and 10. are currently implemented in
  machdep/paral1.src so there is no need for paral2.src and paral3.src
  files, which will eventually become unnecessary. Efficiency of
  different topologies also varies with the number of nodes.


  Also on some platforms EXPAND keyword is recommended in the combination
  of the fastest FAST option in the CHARMM input script, eg for IBMSP:
  EXPAND (fast parvect)



  The installation script now installs default configuration for any
  parallel platform. If one of X,G,P,M,1,2,64,Q,S is specified size
  keyword must be specified too. Run install.com without parameters
  for current set of options.

  Installation command for parallel machines with relevant options:

  1. Cray T3E

  install.com t3e [size] [Q] [P] or [M]

  2. Cray T3D

  install.com t3d [size] [Q] [P] or [M]

  3. Cray C90, J90

  install.com t3d [size]

  4. SGI Power Origin

  install.com sgi64 size M [Q] [X]

  uname -a : IRIX64 icpsg1 6.2 03131016 IP25

  5. SGI Power Chellenge

  install.com sgi size P 64 [Q] [X]

  uname -a : IRIX64 icpsg1 6.2 03131016 IP25

  4a. SGI Origin

  install.com sgi64 size M 64
  /usr/include
  /usr/lib64

  uname -a : IRIX64 atlas 6.5 04131233 IP27

  6. Convex SPP-1000 or SPP-2000

  install.com cspp size P or M [Q]

  7. Intel Paragon machine

  install.com intel

  uname -a : Paragon OSF/1 timewarp 1.0.4 R1_4 paragon

  8. IBM SP1/SP2 machines

  install.com ibmsp size [Q]

  uname -a: AIX f1n3 1 4 000104697000

  8a. IBM SP3 machines

  install.com ibmsp3 size [Q]


  9. Generic Parallel Virtual Machine (PVM)

  install.com machine size P

  10. TERRA 2000

  install.com terra size

  11. Workstation clusters

  install.com machine size S [Q] [X]

  12. Alpha Servers (SMP)

  install.com alphamp size M

  13. Cluster of PCs using GNU/Linux OS - Beowulf class of machines

  0. Compile for most 64-bit machines using multilevel parallelism
     in CHARMM (Spring 2009):

     install.com gnu xxlarge x86_64 M +REPDSTR +MSCALE +ASYNC_PME +ALTIX_MPI

  A. Compile for AMD-64 machines in parallel using MPICH-1
     =====================================================

     ./install.com gnu size M MPICH AMD64

     Note, that MPICH must be compiled with -mcmodel=medium


  B. Using RedHat-6.0:
     =================
     Get and instal the official LAM MPI rpm package from
     rpm -i http://www.mpi.nd.edu/downloads/lam/lam-6.31b1-tcp.1.i386.rpm

     install.com gnu size M [Q] [X] # this asks 2 question - answers are:
     /usr/local/lam-6.3-b1/include
     /usr/local/lam-6.3-b1/lib

  C. Using Debian-potato:
     ====================
     One can use g77 with either lam or mpich (preferred)

     install.com gnu size M [Q] [X] # this asks 2 question - answers are:
     /usr/include/lam
     /usr/lib/lam/lib

     or

     install.com gnu size M mpich [Q] [X] # this asks 2 question - answers are:
     /usr/lib/mpich/build/LINUX/ch_p4/include
     /usr/lib/mpich/build/LINUX/ch_p4/lib

  This small performance table executed on a single processor Pentium
  II/450MHz machine might help you to decide which system/compiler is
  best for your needs:

  B1 = 50 steps of MbCO dynamcs + water with spherical cutoffs
  B2 = 25 steps of MbCO dynamcs + water with PM Ewald
  B3 = 10 steps of minimization of QM/MM for alanine

  All timing in seconds of elapsed time on empty machines using the
  above install procedure. (This table was made July 31, 99).

  Benchmark  |  g77/RH-6.0 | g77*/Debian| f2c/Debian | pgf77   | f77/Absoft
  =========================================================================
      B1     |    290.6 s  |   197.6 s  |  211.1 s   | 189.5 s |  196.0 s
  -------------------------------------------------------------------------
      B2     |    223.3 s  |   193.7 s  |  234.6 s   | 199.2 s |  211.3 s
  -------------------------------------------------------------------------
      B3     |     70.5 s  |    64.3 s  |   74.3 s   |  59.8 s |not working
  =========================================================================

  g77*/Debian is the newest g77-2.95 compiler from July 31, 1999. pgf77
  and f77/Absoft are also the most recent versions.

  [NOTE: pgf77 and MPI don't work out of the box. One has to recompile
  MPI library with explicit pgf77 support. Also, these are the findings
  running testcases (July 1999):


         compiler   |  f2c  | g77  | pgf77
         ---------------------------------
         NORMAL T.  |  152  | 152  |  133
         ---------------------------------
         ABNORMAL T.|   26  |  26  |   23
         ---------------------------------
         segm. fault|    4  |   4  |   23
         ---------------------------------
         total #    |  186  | 186  |  186
         ---------------------------------

  The difference between the total and the sum of other numbers is in
  the problems of CHARMM testcases suite.
  ]


  14. IBM Power4 using GNU/Linux system:
  ---------------------------------------

  Obtain the MPICH library from http://ppclinux.ncsa.uiuc.edu

  install.com imblnxmp size M

  Then depending on the system version you might get some "relocation
  truncated..." errors. If this happens, run:

  tool/ibmlnxmp_fixlibs
  cp build/UNX/Makefile_ibmlnxmp_so build/ibmlnxmp/Makefile_ibmlnxmp

  This procedure should produce an executable in exec/ibmlnxmp/charmm

  Additional note:

  Also it is needed to change INTEGER statements in mpif.h file into INTEGER*4


                       -----

  The following keywords in pref.dat are used for parallel CHARMM:

  Machine independent keywords:

  PARALLEL        Needed for parallel version
  SOCKET          If TCP/IP sockets
  PVM             If using PVM library
  PVMC            If using PVM library on some platforms (see below).
  PARAFULL        Currently the only one which works
                  (must be specified)
  PARASCAL        For force decomposition scheme
                  (not ready for general use yet.)
  SPACDEC         For spatial decomposition scheme
                  based on BYCC (BYCC must be specified in nonbond
                  options)
  SYNCHRON        Most of the machines don't do
                  receive and send at the same time
  GENCOMM         Different communication arcitecture.
                  Can run any number of nodes
  MPI             If using MPI parallel library.
                  (point-to-point routines only)
  CMPI            CHARMM implementation of the MPI library.
                  Enables all the old functionality plus some
                  combinations of libraries on the same platform.
  ASYNC_MPI       using CMPI library routines vs MPI in PME.


  Machine specific keywords:

  TERRA
  CM5
  CSPP
  DELTA
  INTEL
  PARAGON
  SHMEM
  CSPPMPI
  T3D
  T3E
  IBMSP
  ALPHAMP
  SGIMP
  ALTIX_MPI   ! also used in generic x86_64 compiles


  
  File: parallel.doc,  Node: Running,  Next: PARAllel,  Prev: Installation Top,  Up: Top

  Running CHARMM on parallel systems

  General note for MPI systems.
  Most MPI systems do not allow rewind of stdin which means charmm input files
  containing "goto" statements would not work if invoked directly
  (this example uses MPICH):
  ~charmm/exec/gnu/charmm -p4wd . -p4pg file < my.inp > my.out [charmm options]

  The workaround is simple:
  ~charmm/exec/gnu/charmm -p4wd . -p4pg file < my.stdin > my.out ZZZ=my.inp [charmm options]

  where the file my.stdin just streams to the real inputfile:
  * Stream to real file given as ZZZ=filename on commandline. Note that the filename
  * cannot consist of a mixture of upper- and lower-case letters.
  *
  stream @ZZZ
  stop

    1. Cray  T3D (Cray-PVM)

            ~charmm/exec/t3d/charmm24 -npes 256 < input_file > output_file &

       The same command may be used in a batch script but without `&'.
       Example using batch:

            #QSUB -lM 16Mw
            #QSUB -lT 600:00
            #QSUB -mb -me
            #QSUB -l mpp_p=32
            #QSUB -l mpp_t=600:00
            #QSUB -q mpp
            setenv MPP_NPES 32
            ~charmm/exec/t3d/charmm24 < Input_file > output_file

       Preflx directives required: T3D UNIX PARALLEL PARAFULL
       Additional preflx directives recommended: PVM or MPI

    2. Cray  T3E (Cray-PVM)

       CHARMM can be run on either a single processor or in parallel on the T3E.
       Single processor runs are useful for small analysis jobs and other tasks
       that are not amenable to parallel processing. The syntax for a single
       pe run is:
           charmm24 < filename.inp >& filename.out [&]
       Large CHARMM jobs should be run in parallel using the queue system.
       The syntax for a parallel run is:

          mpprun -n# charmm24 < filename.inp >& filename.out [&]
          (here # is the desired number of pe's)

       The same command may be used in a batch script but without `&'.

       Example using batch:
            #QSUB -lM 16Mw
            #QSUB -lT 600:00
            #QSUB -mb -me
            #QSUB -l mpp_p=32
            #QSUB -q mpp
            mpprun -n 32 charmm24 < Input_file > output_file

       Preflx directives required: T3E UNIX PARALLEL PARAFULL
       Additional preflx directives recommended: EXPAND(fast off)
                                                 and either PVM or MPI

       Optimization Notes:
       T3E users should use the PBOUND command for simulations of periodic
       systems.  The pbound command optimizes non-bonded list-generation and
       computations on parallel machines such as the T3E, giving significantly
       better performance for parallel applications using simple perodic
       boundary conditions. Note that the pbound command is currently
       implemented only for scalar architectures such as the T3D and T3E.


    3. Cray C90, J90 (Cray-PVM)

       No info yet

    4. SGI Power Challenge (PVM)

            pvm
            quit

            setenv NTPVM 16 (or NTPVM=16 ; export NTPVM)
            ~charmm/exe/sgi/charmm24 <input_file >output_file &

       Preflx directives required: SGI UNIX PARALLEL PARAFULL CMPI PVMC SGIMP
       Additional preflx directives recommended: EXPAND(fast off)
       Alternative, but not tested yet: SGI UNIX PARALLEL PARAFULL

       [NOTE: This is old: MPI is preffered over this. Installation
              similar to Linux, see above]

    5. Convex SPP-1000 Exemplar

       With PVM
            (see below for information setting up a PVM Hostfile)
            mpa -sc <name_of_subcomplex> /bin/csh
            setenv PVM_ROOT /usr/convex/pvm
            /usr/lib/pvm/pvm
            quit

            ~/pvm3/bin/CSPP/charmm24 -n 16  <input_file >output_file &
            ~charmm/exe/cspp/charmm24 <input_file >output_file &

       Which subcomplexes are available check with the scm utility.

       (For information on how to set up a PVM hostfile see *note 1: Using PVM.)
       Preflx directives required: CSPP UNIX PARALLEL PARAFULL PVM HPUX
       SYNCHRON (GENCOMM)

       Note: The first time that you build CHARMM with PVM specify the P option
       with install.com.  You will be asked for the location of the PVM include
       files and libraries. If these do not change and you do not reconstruct the
       Makefiles, you do not have to specify this option each time you run
       install.com.

       With MPI

            mpa -DATA -STACK -sc <name_of_subcomplex> \
            ~charmm/exe/cspp/charmm24 -np <n> <input_file >output_file &
       Where <n> is the number of processors to use.
       There are two environmanet variables that can be set:
            setenv MPI_GLOBMEMSIZE  <m>
       where <m> is the size of the shared memory region on each hypernode
       in bytes.  The default is 16MB.
       And:
            setenv MPI_TOPOLOGY <i>,<j>,<k>,<l>,...
       where <i>, <j>, <k>, <l>, ... are the number of tasks on each hypernode.
       The sum must equal the number of processors specified with -np on the
       command line.  This is optional the default behavior is generally what
       you want.  If you are using a sub-complex with more than one hypernode,
       use may want to include '-node 0' after mpa to keep the 0th process
       on the 0th hypernode of the sub-complex.

       Preflx directives required: CSPP UNIX PARALLEL PARAFULL HPUX
       MPI CSPPMPI

       The CSPPMPI directive specifies the use of extensions in the Convex
       MPI implementation. This directive is optional. Use of the MPI
       directive alone will result in a fully MPI Standard compliant program,
       albeit with a loss of performance.

       Note: The first time that you build CHARMM with MPI specify the M option
       with install.com.  You will be asked for the location of the MPI include
       files and libraries. If these do not change and you do not reconstruct the
       Makefiles, you do not have to specify this option each time you run
       install.com.

    6. Intel gamma

       Because the fortran compiler on the Intel gamma does not know how
       to rewind the redirected input file the program uses charmm.inp
       file name from current working directory. The script for running
       CHARMM should look like the following example:

            cp input_file.inp charmm.inp
            getcube -t128 > output_file
            load ~charmm/exec/intel/charmm24
            waitcube

       Preflx directives required: INTEL UNIX PARALLEL PARAFULL

    7. Intel Delta

            mexec "-t(32,16)" ~charmm/exec/intel/charmm23<input_file>output_file&

       Preflx directives required: INTEL UNIX DELTA PARALLEL PARAFULL

    8. Intel Paragon

            ~charmm/exec/intel/charmm23 -sz 64 <input_file >output_file &

       Preflx directives required: INTEL UNIX PARAGON PARALLEL PARAFULL

    9. CM-5

            ~charmm/exec/cm5/charmm23 <input_file >output_file &

       Preflx directives required:CM5 UNIX PARALLEL PARAFULL

   10. IBM SP2 or SP1

            setenv MP_RESD yes
            setenv MP_PULSE 0
            setenv MP_RMPOOL 1
            setenv MP_EUILIB us
            setenv MP_INFOLEVEL  0
            poe ~charmm/exec/ibmsp/charmm24 -hfile nodes -procs 64 <input >output

       See `man poe'  for details.

       Preflx directives required:IBMSP UNIX PARALLEL PARAFULL
       Additional preflx directives recommended: EXPAND(fast parvect)

   11. PVM

            pvm
            add host host1
            add host host2
            quit
            setenv NTPVM 3
            ~/pvm3/bin/SGI5/charmm24 <input_file >output_file&

       Preflx directives required: machine_type UNIX PARALLEL CMPI PVM
       PARAFULL SYNCHRON

   12. Linux clusters (Beowulf)

       MPICH: (MPICH doesn't need to be installed on compute nodes)

       ~charmm/exec/gnu/charmm -p4wd . -p4pg file < input > output [charmm options]

       where file is:
       host1 0
       host2 1 ~charmm/exec/gnu/charmm
       host3 1 ~charmm/exec/gnu/charmm
       etc.

              [NOTE: host1 can be the same as host2, host3, etc. for
                     SMP]

       LAM: (Every node has to have LAM installed!!)

       lamboot -v hostfile
       mpirun -O -c2c -w schema < input >output

       where schema is a file:
       ~charmm/exec/gnu/charmm n0 -- [charmm options]
       ~charmm/exec/gnu/charmm n1 -- [charmm options]
       ~charmm/exec/gnu/charmm n2 -- [charmm options]
       etc.

       and hostfile is:
       host1
       host2
       host3
       etc.

   13. PARALLEL VERSION OF CHARMM23 ON WORKSTATION CLUSTERS

       Preflx directives required: machine_type UNIX PARALLEL CMPI SOCKET
       PARAFULL SYNCHRON

       Currently the code runs on HP, DEC alpha, and IBM RS/6000
       machines. This has been tested.  The rest of UNIX world should run
       too without any changes as long as the following is true:

       Assumptions for cluster environment:

       Before you run CHARMM with SOCKET library you have to define some
       environment variables.  If you define nothing then CHARMM will
       run in a scalar mode, i.e.  default is one node run.

       PWD

       The program supports three shells: bash (Bourne Again SHell), ksh
       (Korn Shell) and tcsh, which is available from anonymous ftp. The
       only difference from csh on which CHARMM makes assumption is
       definition of variable PWD. This variable is correctly defined in
       all of the above three shells by default, while using csh it has
       to be defined by the user. Variable PWD points to the current
       working directory. If some other directory is requested the PWD
       environment variable may be changed appropriately. The program
       can figure out current working directory by itself but there are
       problems in some NFS environments, because home directory names
       can vary on different machines.( PWD is always defined correctly
       by shell which supports it ) So csh may sometimes cause
       problems. Using csh the cd command may be redefined so that it
       always defines also PWD. This is done with something like: alias
       cd 'chdir \!*; setenv PWD $cwd ' in the ~/.cshrc file.

       If you get an error which looks something like nonexistent
       directory then define PWD variable directly.

       [NIH specific (for HPUX):
       If you want to use tcsh as your login shell you may run the
       following command:
            runall chsh username /usr/local/bin/tcsh

       runall is a script which runs the command on the whole cluster of
       machines it is on /usr/local/bin at NIH.  ]

       NODEx

       In order to run CHARMM on more then one node environment variables
       NODE0, NODE1, ...,  NODEn have to be defined.

       Example for a 4 node run:

            setenv NODE0 par0
            setenv NODE1 par1
            setenv NODE2 par2
            setenv NODE3 par4

            charmm < input_file > output_file 1:parameter1 2:parameter2 ...

       "par0,par1,par2,.." are the names of the machines in the local
       network.  There is no requirement that all machines should be of
       the same type. There is nothing in the program to adjust for
       unequal load balance so all nodes will follow the slowest one. In
       near future we may implement dynamic load balance method based on
       actual time required.

       The assumption here is that the node from where CHARMM program is
       started is always NODE0!

       Setup for your login environment

       In order to run CHARMM in parallel you have to be able to rlogin to
       any of the nodes defined in NODEx environment variables. Before you
       run CHARMM check this out:

       rlogin $NODE1

       if it doesn't ask you for Password then you are OK. If it asks for
       Password then put a line like this:
            machine_name user_name

            in your ~/.rhosts file, with 600 permission.

       [NIH specific:
       How to submit job to HP.

       Currently we have assigned machines par0, par1, par2, and par4 to
       work in parallel. You may use script
       /usr/local/bin/charmm23.parallel and submit it to par0. Example:

       submit par0 charmm23.parallel <input_file >output_file ^D

       To construct your own parallel scripts look at
       /usr/local/bin/charmm23.parallel ]

       In the input scripts

       Everything should work, but avoid usage of IOLEV and PRNLEV in your
       parallel scripts.



  
  File: parallel.doc,  Node: PARAllel,  Next: Status,  Prev: Running,  Up: Top

  Syntax:

  PARAllel { FIFO       int                            }
           { BUFFer     int                            }
           { CONCurrent int  [ COUNT int  MAXI int ]   }


  Description:

  FIFO specifies priority for the Linux kernel FIFO scheduling
  scheme. Larger number means higher priority. Zero is for the default
  scheduling scheme.

  BUFFer specifies the size of the sending and receiving buffer for the
  MPI send/receive calls. It is in Real*8 units.

  CONCurrent specifies the number of independent CHARMM jobs within a
  single parallel run. If COUNt=0 it turns on the interleaving
  communication between the 2 groups, ie one group is performing
  communication while the other is doing calculation at the same
  time. Interleaving stops after MAXI steps of dynamics.

  Example:

  The following example performs interleaving between 2 jobs. The total
  number of nodes allocated has to be even. The input for job 1 has to
  be in the file with the name 1.input and for job 2 in 2.input.

  * This input script runs 2 separate jobs
  *

  paral conc 2 count 0 maxi 102 ! 1.input & 2.input are currently
                                ! hardwired into paral1.src


  
  File: parallel.doc,  Node: Status,  Prev: PARAllel,  Up: Top, Next: Using PVM

  Parallel Code Status (as of July 2003)

  NOTE: c31a1 test directory contains 276 testcases. Out of these 22
  cannot stop the execution by themself. 8 tests end with the ABNORMAL
  termination and 246 with NORMAL termination, which of course this
  doesn't guarantee that the method is working in parallel.

  The following table is the result of this testing.



  The symbol ++ indicates that parallel code development is underway.

  -----------------------------------------------------

  Fully parallel and functional features:

       Energy evaluation

       ENERgy, GETE, SKIPE, ENERgy ACE

       MINImization (CONJ,NRPH,ABNR,POWEL,TN)

       DYNAmics (leap frog integrator)

       HBOND

       BLOCK

       CRYSTAL (all)

       IMAGES

       INTEraction energy

       CONStraints (SHAKE,HARM,IC,DIHEdral,FIX,NOE,RESD,LONEPAIR)

       ANAL (energy partition)

       NBONds (generic)

       EWALD

       PME

       PERT

       GAMESS (ab initio part)

       TEST FIRST, SECOND

       REPLICA

       TREK

       EEF1

       IMCUBES (bycb)

       FSSHK   (fast non-vector shake)

       GENBORN

       GBBLOCK

       GRID

       HMCM

       BYCC

       TSM

       TMD

       GRAPE

       HQBM

       PSSP

       ADUMB

       MTS

       SSBP

       DRUDE

       VV2

       LONEPAIR

       QCHEM

       GAMESSUK

       RPATH

       QUB

       FACTS

  -----------------------------------------------------

  Functional, but nonparallel code in the parallel version (no speedup):
  ( ** indicates that these can be very computationally intensive and are
  not recommended on parallel systems)

       VIBRAN  **

       CORREL **(Except for the energy time series evaluation, which is
       parallel)

       READ, WRITE, and PRINT (I/O in general)

             NOTE:
             always protect prnlev ...
             with
             if ?mynode .eq. 0 then prnlev ...

       CORMAN commands
              COPY, ORIENT, CONVERT, SURFACE,
              CONTACT, VOLUME, LSQP, RGYR

       HBUIld **

       IC (internal coordinate commands)

       SCALar commands

       CONStraints (setup, DROPlet, SBOUnd)

       Miscellaneous commands

       GENErate, PATCh, DELEte, JOIN, RENAme, IMPAtch (all PSF
       modification commands)

       MERGE

       QUANtum **  ++

       QUICk

       REWInd (not fully supported on the Intel)

       SOLANA

       SELECT

       DEFINE

       MONITOR

       TEST

       CMDPAR and flow control

       PATH

       RXNCOR

       Commandline parameters (where supported by compiler)

       RISM

       ZMAT

       AUTOGEN

       CALC

       BOUND

       HELIX

       WHAM

       GRAPHICS

       UMBRELLA

       SBOUNDARY

       PBEQ  ++

       GSBP

  -----------------------------------------------------

  Nonfunctional code in parallel version:

       ANAL (table generation)

       DYNAmics (old integrator, NOSE integrator, 4D)

       MMFP

       TRAVEL

       VIBRAN (quasi, crystal)

       BLOCK FREE

       COOR COVARIANCE

       ST2 waters

       NMR

       DIMB

       ECONT

       PULL

       CFTI

       LUP

       GALGOR

       BYCU

       MC

       4D DYNA

       SCPISM

  -----------------------------------------------------

  Untested Features (we don't know if it works or not):
       ANALysis

       MOLVIB  (minor problems with I/O - hangs the job)

       PRESsure (the command)

       RMSD

       MBOND

       MMFF

       SHAPES

       CLUSTER


  
  File: parallel.doc, Node: Using PVM, Prev: Status, Up: Top, Next: Implementation


  Note:   Currently one should specify the absolute path to the pvm include
          files and the pvm library files.  This is done because PVM installation
          is not currently standard.  During installation, through use of
          install.com, you are asked to specify these paths.


  Convex PVM

  This version runs using PVM (Parallel Virtual Machine) versions 3.2.6 and
  higher. To run:

      1. create hostfile - as in the example below:

         #host file
         puma0 dx=/usr/lib/pvm/pvmd3 ep=/chem/sfleisch/c24a2/exec/cspp

         The first field (puma0) is the hostname of the machine.  The dx= field
         is the absolute path to the PVM daemon, pvmd3. This includes the
         filename, pvmd3.  The last field, ep= is the search path for find the
         executable when the tasks are spawned. This can be a colon (:) separated
         string for searching multiple directories. The PVM system can be
         monitored using the console program, pvm.  It has some useful commands:

             conf   list machines in the virtual machine.
             ps -a  list the tasks that are running.
             help   list the commands.
             quit   exit the console program without killing the daemon.
             halt   kill everything that is running and the daemon and exit
                    the console program.


      2. Run the PVM daemon, pvmd3:

             pvmd3 hostfile &

      3. Run the program e.g.:

         /chem/sfleisch/c24a2/exec/cspp/charmm -n <ncpu> <input_file >output_file
  &

         where -n <ncpu> indicates how many pvm controlled processes to run

      4. Halt the daemon. See above.

  The Convex Exemplar PVM implementation uses shared memory via the System V
  IPC routines, shmget and shemat.

  Generic PARALLEL PVM version for workstation clusters

  Preflx directives required: <MACHTYPE> UNIX SCALAR CMPI PVM PARALLEL
                                                         PARAFULL SYNCHRON

  Where <MACHTYPE> is the workstation you are compiling on, e.g.,
  HPUX, ALPHA, etc.

  Note:   Currently one must specify the absolute path to the pvm include
          files and the pvm library files.  This is done because PVM installation
          is not currently standard.  During installation, through use of
          install.com, you are asked to spceify these paths.

  This version runs using PVM (Parallel Virtual Machine) versions 3.2.6 and
  higher. To run:

    1. create hostfile - as in the example below:

       #host file
       boa0 dx=/usr/lib/pvm/pvmd3 ep=/cb/manet1/c24a2/exec/hpux
       boa1 dx=/usr/lib/pvm/pvmd3 ep=/cb/manet1/c24a2/exec/hpux
       boa2 dx=/usr/lib/pvm/pvmd3 ep=/cb/manet1/c24a2/exec/hpux
       boa3 dx=/usr/lib/pvm/pvmd3 ep=/cb/manet1/c24a2/exec/hpux

       The first field (boa0, etc) is the hostname of the machine. The dx= field
       is the absolute path to the PVM daemon, pvmd3. This includes the
       filename, pvmd3.  The last field, ep= is the search path for find the
       executable when the tasks are spawned. This can be a colon (:) separated
       string for searching multiple directories. The PVM system can be
       monitored using the console program, pvm.  It has some useful commands:

             conf   list machines in the virtual machine.
             ps -a  list the tasks that are running.
             help   list the commands.
             quit   exit the console program without killing the daemon.
             halt   kill everything that is running and the daemon and exit
                    the console program.


    2. Run the PVM daemon, pvmd3:

             pvmd3 hostfile &

    3. Run the program e.g.:

         /cb/manet1/c24a2/exec/hpux/charmm -n <ncpu> <input_file >output_file &

         where -n <ncpu> indicates how many pvm controlled processes to run

    4. Halt the daemon. See above.

  
  File: parallel.doc, Node: Implementation, Prev: Using PVM, Up: Top, Next: Top

  Implementation notes.
  =====================

  Currently the support for parallel machines in CHARMM is implemented
  in three levels. The topmost level is the collection of subroutines
  which are called from CHARMM itself. These subroutines are implemented
  in paral1.src. They are:

  VDGSUM  - vector distributed global sum [MPI_REDUCE_SCATTER]
  VDGBR   - vector distributed global broadcast [MPI_ALLGATHERV]
  GCOMB   - Global combine (sum) [MPI_ALLREDUCE]
  VDGBRE  - vector distributed global broadcast (one vector only) [MPI_ALLGATHERV]
  PSNDC   - Broadcast character array from node 0. [MPI_BROADCAST]
  PSND4   - Broadcast integer array from node 0. [MPI_BROADCAST]
  PSND8   - Broadcast real*8 array from node 0. [MPI_BROADCAST]
  PSYNC   - Barrier [MPI_BARRIER]
  PARFIN  - Close the parallel setup [MPI_Finalize]
  PARSTRT - Start and setup for parallel
  PARCMD  - PARAllel command parser

  The above routines then by default call the MPI equivalents as
  indicated above. Since the current status of MPI implementations is
  not efficient on most of the parallel platforms we still maintain the
  CHARMM implementation of MPI chosen by CMPI preflx keyword in pref.dat
  file. Besides the choice of standard MPI library and CMPI there are
  other choices available in paral1.src for the vendor specific
  libraries which have similar functionality as MPI library. Currently
  these are CSPP and TERRA options. So in short paral1.src is a place
  where one decides which library will be used for global parallel
  communication, such as global sum and others. It may also decide on
  machine specific libraries if they differ from MPI, but provide the
  same functionality (TERRA example).

  For the users of MPI library there are always two possibilities:

  1. Don't specify anything except PARALLEL PARAFULL in pref.dat and use
     global communication as implemented in MPI.

  2. Specify PARALLEL PARAFULL CMPI MPI and use the efficient global
     communication algorithms implemented the paral2.src and paral3.src,
     where only two primitive MPI calls are used: send and recieve. This
     choice is currently the preferred one on most of the systems
     especially for users of MPICH and its derivatives.

  Once CMPI keyword is specified the routines in paral1.src call
  another set of routines implemented in the paral2.src source file. The
  purpose of routines in this layer is to decide on which topology will
  be chosen for a given parallel system. Possible choices are:

  1. recursive halving sutable for hypercube or switched networks. This
     is the default selection.

  2. ring topology suitable for ring networks or any other where non
     power of two number of processors is selected. This is selected at
     compile time with the keyword GENCOMM in pref.dat.

  3. mesh topology for two dimensional mesh network connection, also
     sometimes works the best with FAT tree topology. Selected by
     DELTA in pref.dat.

  4. Each of the topology is by default implemented using send/receive
     routine which is capable of receiving data from the other processor
     while sending to it at the same time. If this is not supported by
     the hardware one can choose SYNCHRON keyword in pref.dat.

  All of the above topologies are then implemented in paral3.src file
  for a variety of parallel systems.

  I/O requirements for the new code
  =================================

  Each fortran WRITE statement has to be protected by PRNLEV, for
  example:

        IF(PRNLEV.GT.2) WRITE(OUTU,55) CALLNAME,N,INBLOX(NATOM)

  instead of just simply:

        WRITE(OUTU,55) CALLNAME,N,INBLOX(NATOM)


  READ statements are a little bit more complicated and they are
  controled by IOLEV. Example:

        IF(IOLEV.GT.0) THEN
           READ(UNIT)(X(I),Y(I),Z(I),I=1,NATOM)
        ENDIF
  ##IF PARALLEL
        CALL PSEND8(X,NATOM)
        CALL PSEND8(Y,NATOM)
        CALL PSEND8(Z,NATOM)
  ##ENDIF

  Any further information can be obtained from milan@cmm.ki.si.
  See also the current parallel performance table at:
  http://arg.cmm.ki.si/parallel/summary.html

