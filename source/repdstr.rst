.. py:module:: repdstr

================================
The Parallel Distributed Replica
================================

By Paul Maragakis and Milan Hodoscek, 2005

Parallel ditributed replica allows independent replicated systems
over specified number of processors. It mainly works with CMPI
pref.dat keyword (YMMV). REPDSTR is still not the default pref.dat
keyword so the recommended way to compile CHARMM is the following:

::

  install.com gnu xxlarge M mpif90 +REPDSTR +ASYNC_PME +GENCOMM [+MSCALE]

for ``install.com em64t`` add ``+CMPI`` to the above list.

Also MSCALE is not really needed for pure REPDSTR runs, but it is
needed for triple parallel CHARMM.

For one of the examples of REPDSTR usage see this reference:

Jiang, W; Hodoscek, M; Roux, B; "Computation of Absolute Hydration and
Binding Free Energy with Free Energy Perturbation Distributed
Replica-Exchange Molecular Dynamics", J. Chem. Theo. and Comp., 2009,
Vol. 5, 2583-2588.

For the reservoir replica exchange code, please cite the following
references:

* Boltzmann reservoir REX--

  Okur A., Roe D., Cui G., Hornak V., Simmerling C. J. Chem Theo. Comput. 3, 557-568 (2007).

* Non-boltzmann reservoir REX--

  Roitberg A., Okur. A., Simmerling C. J. Phys. Chem. B. 111, 2415-2418 (2007).


.. _repdstr_syntax:

Syntax
------

::

   REPDstr NREP <int>

Replicates the system <int> times. Optional NATRep limits the
number of atoms to be included in the path calculations (RPATH
commands). It also reduces the size of arrays that need to be
transfered between replicas in the RPATH calculations.


::

    REPDstr NREP <int>  [REPEat <int>] { EXCHange replica-exchange-spec }
                                       { PHREx ph-replica-exchange-spec }

    replica-exchange-spec::= FREQuency <int> temperature-spec
                             [UNIT <int>] [SUMP]
                             [NREP <int>]
                             [SGLD] sgld-replica-exchange-spec
                             [TIGEr] tiger-spec
                             [RSRV] reservoir-spec

    temperature-spec::= [ [ TEMP <real> ] TEMP <real> ... ]
                        [STEMperature <real> DTEMperature <real> MTEMperature <real>]

    tiger-spec::= [ITER <int> ] [NEQU <int>] [NMIN <int>] [TOLG <real>]

    reservoir-spec ::=  RESHigh  { BOLTzmann   }  RHTEmp <real> RLTEmp <real>
                        RESLow   { NOBOltzmann }  NINRes <int> RHUNit <int> RLUNit <int>

    ph-replica-exchange-spec ::= FREQuency <int> ph-spec [UNIT <int>]

    ph-spec ::= [ [PHVAl <real] PHVAl <real> ... ]

    sgld-replica-exchange-spec::=  [ [ SGTT <real> ] SGTT <real> ... ]
                                   [SGTE <real> DSGT <real> MSGT <real>]
                                   [SGFT <real> ] DSGF <real> ... ]


This is for the replica exchange method (see details in
c34test/rexc.inp, c36test/rexcpt.inp, c36test/rexsgld.inp) Currently
it works so that when exchange occurs all the coordinates and
velocities are exchanged, thus the lowest temperature is always on the
first replica. This also implies that the number of atoms in the
replicas have to be the same. The other method which exchange only the
temperature will be implemented later. With just temperature exchange
the replicas do not need to be the same anymore.

Current implementation of replica exchange methods has its own
temperature control independent of the CHARMM's one. So in the case of
exchanging the coordinates and velocities also the appropriate
temperature scaling is perfomed. Perhaps it is best to turn CHARMM's
own temperature controls off, but one can also combine the two. To get
both temperature control mechanisms at the same time one need to define
different temperature for each replica. This can be accomplished by
the following commands in the CHARMM input script:

::

    set st 300
    set dt 10

    repd nrep @nreps EXCHange FREQuency 50 STEMp @st DTEMp @dt sump unit 17

    mult dt by ?myrep
    incr st by @dt

    dyna cpt start nstep 1000 timestep 0.001 -
    ....
        hoover reft @st tmass 2000.0 tbath @st -
    ....

As of CHARMM version c37a2, replica exchange in pH space is also
supported. This uses the formalism described in the reference given in
:doc:`consph`. The :chm:`CONSPH` key word must be in pref.dat along with :chm:`REPDSTR`
to activate this functionality.

Details about each keyword:

============ =====================================================================
UNIT <int>   Optional keyword for exchange output file. Default for
             int is OUTU. Must be opened after repd: each replica
             writes to its own file. If no open statement all
             replicas write to the same file with the default
             fortran file name for this unit. Open before the repd
             command is not very useful. Can be also the same unit
             as on the OUTU command so exchange info is written to
             the same output files.

SUMPrint     Summary printout. On replica zero the summary from
             all other replicas is printed to UNIT, and on the
             rest of the repicas just their own data. This is
             flaged since it requires extra communication just for
             printouts.

FREQuency    when to exchange

REPEat       Number of times to repeat an exchange attempt every FREQ
             steps.

STEM         Starting-temperature.

DTEM         temperature-increase

MTEM         Top temperature. When MTEM>0, DTEM is ignored and temperatures
             of replicas are expoenentially spaced.

TEMPerature  Temperature of each replica if STEM is not set. It must be
             repeated NREP times.

PHVAl        pH value of each replica when replica exchange in pH space is
             used.

SGLD         Flag to do RXSGLD with the self-guiding
             temperature. It maybe used with the standard replica
             exchange or one can specify all the temperature the
             same, most convenient with the STEM <temp> DTEM 0.0.

SGTE    0    The self-guiding temperature for the first replica.

DSGT    0    Increment for the self-guiding temperature

MSGT    0    The top self-guiding temperature. If MSGT>0, DSGT is ignored and
             the guiding temperatures of replicas are expoenentially spaced.

SGTT         Self-guiding temperature of each replica if SGTE is not set.
             It must be repeated NREP times.

SGFT    0    The guiding factor of the first replica. When the self-guiding
             temperatures are set with SGTE..., SGFT will be adjusted
             automatically during simulation.

DSGF    0    Guiding factor increment.

TIGEr        Flag to start TIGER replica exchange.

ITER         number of iteration steps for minimization and
             equlibration procedures before the exchange.
             Default: 1

NEQU         number of steps in the equlibration process.
             Default: 1000

NMIN         number of steps in the minimization process.
             Default: 100

TOLG         gradient tolerance in the minimization step.
             Default:0.0

PHMD         Flag to allow exchange of theta variables from CPHMD
             along with spatial coordinates. Thus, replicas can be run
             at different pH (Hamiltonian replica exchange) or temperture.
============ =====================================================================

It is also possible to couple the top or bottom replicas (or
both) to a reservoir of structures. To do so, the RSVR keyword is
used. When RSVR is used, at least one of the following keywords
must also be used with the corresponding unit numbers to tell CHARMM
which replica(s) should be coupled to the reservoir.

============ =====================================================================
RESH         couple the top replica to a reservoir. RHUN must
             be specified.

RESL         couple the bottom replica to a reservoir. RLUN
             must be specified.
============ =====================================================================

RHUN and RLUN must be units that point to open files in a simplified
trajectory format. This format is a standard CHARMM binary trajectory
file, but has the header and crystal information stripped out. A utility,
simpletraj.py, is provided in the support/programs directory to convert a standard
CHARMM trajectory into a stripped down version suitable for use.

Two exchange schemes have been implemented to govern coupling of the
reservoir with its neighboring replica.

============ =====================================================================
BOLTzmann    The standard Boltzmann temperature replica exchange
             criterion is used. Use of this keyword implies that
             the reservoir is a sample from a Boltzmann distribution.
             If this option is used, RHTEmp and/or RLTEmp must
             be used to specify the temperatures of the high and
             low reservoirs, respectively.

NOBOltzmann  This option allows for a non-Boltzmann weighted
             reservoir, using a slightly different exchange
             criterion.
============ =====================================================================

In both cases, the :chm:`NINRes` parameter is required. NINR tells CHARMM how
many structures are initially in the reservoir.

::

    REPDstr RESEt [ SYNC ] [ PONE ]

Resets the run to a normal parallel run. SYNC does the global
sync before that. PONE is making for everybody NUMNOD=1. As of March
2010 RESET is still not fully supported.

::

    REPDstr IORES

Sometimes within the REPDstr run one wants to access the files
created by other replicas. After this command is executed the names in
the open command do not get _irepdstr appended!

::

    REPDstr IOSET

Sets the appending of the replica number back to original nameing
scheme in REPEDstr.

::

    REPDstr NREP <int> EXLM [EXPT NRPT <int>] FREQuency <int>

This is for Hamiltonian exchange method. Currently it works so that
when exchange occurs all the coordinates are exchanged and new nonbond list
are generated. To guarantee stable md run after exchange, velocities also
are exchanged once an exchange attempt is accepted. The present Hamiltonian-
exchange scheme works for all integrators, including VV2 integrator for
Drude oscillator model.

=============== =========================================================================
EXLM            Keyword invoking Hamiltonian exchange. Currently it can be used
                in Free Energy perturbation and umbrella sampling.

EXPT NRPT <int> Optional keyword introducing parallel tempering into
                Hamiltonian exchange. With this keyword, the replica-exchange
                consists of two alternative stages: parallel tempering and Halmiltonian
                exchange. In parallel tempering stage, the number of replicas
                participating exchange is NREP, while in Hamiltonian exchange stage,
                the number of replica is NREP/NRPT. Currently it can be used to
                accelerate the relaxation of internal degrees of freedom, such as
                sidechain dynamics and backbone dynamics
=============== =========================================================================

::

    REPD NREP <int> EXLM EX2D NRPX <int> FREQ <int>

This is a new 2 Dimensional Hamiltonian Replica exchange scheme.
Hamiltonian-Exchange is extended to PBC systems for either NVT and NPT
simulation. This new feature is especially useful to enhance samplings of
umbrella sampling that involve multiple reaction coordinates.

=============== ===========================================================================
EX2D NRPX <int> Optional keyword introducing 2D replica exchange.
                With this keyword, the number of replicas along X (one reaction coordinate)
                is NREPX, then the number of replicas along the other reaction coordinate
                is NREP/NREPX.
=============== ===========================================================================

This is for replica exchange method (see details in c34test/rexc.inp)
Currently it works so that when exchange occurs all the coordinates
and velocities are exchanged, thus the lowest temperature is always on
the first replica. The other method which exchange only the
temperature will be implemented later.

============ =====================================================================
UNIT <int>   Optional keyword for exchange output file. Default for
             int is 6. Must be opened after repd: each replica
             writes to its own file. If no open statement all
             replicas write to the same file with the default
             fortran file name for this unit. Open before the repd
             command is not very useful.

FREQuency    when to exchange

TEMPerature  the temperature of each system. It must be repeated NREP times
============ =====================================================================

Instead of repeated TEMP comands for equidistant temperature
intervals one can use:

::

    STEM <starting-temperature> DTEM <temperature-increase>


.. _repdstr_io:

I/O
---

Once REPDstr command is activated the I/O capabilities of CHARMM
are expanded. In standard parallel mode CHARMM deals with I/O only on
the first process. The rest of processes get their data through
network or memory communication. So all I/O statements that are in the
script before REPDstr command are valid only on first process. In
distributed replica mode each replica needs its own and independent
I/O which is enabled after the REPDstr keyword in the input script.
Two substitution parameters are defined after the REPD command is
specified in the input script: ?NREP (number of replicas) and ?MYREP
(current executing replica).

As of May 2009 the following is working:

1. OPEN

   The command open read|write unit 1 card name somefile will open
   somefile_0 for replica 0, somefile_1 for replica 1, etc

2. READ/WRITE

   writes to individual files one for each replica. It works for all
   I/O operations.

3. STRE stream

   This will open stream_0 for replica 0, stream_1 for replica 1, etc
   It allows CHARMM to run different input files for each replica (or
   group of processors)

4. OUTU unit

   Will stream output to individual files as specified in the open
   command for particular unit. This command should precede STRE
   command if one wants both input and output files for each group of
   processors

5. ``IF ?MYREP .EQ. n THEN ....``

   Works, too. Output only for the processor zero, unless OUTU is
   specified.

6. All the above works in parallel/parallel mode, ie each replica can
   be a parallel job in itself. The numeration of input and output
   files follows the replica numbers.

   The output is written only on a local process 0 for each replica,
   and similar is true also for stream command.  The limitation is
   that the number of replicas must divide the number of processes
   allocated for parallel. Otherwise it bombs out with the level -5.


.. _repdstr_examples:

Examples
--------

.. note::

   If you are using mpich-1.2.X then you need to use -p4wd with the
   absolute path or -p4wd `pwd`


Example 1
^^^^^^^^^

::

   read psf
   read coor
   repd nrep 4

This will replicate PSF and coordinates, so after nrep 4 there are
four independent runs with the same coordinates

Example 2
^^^^^^^^^

::

   read psf
   repd nrep 4
   read coor name system.crd

This will replicate PSF but the coordinates will be read from 4
separate files: system.crd_0, system.crd_1, etc

Example 3
^^^^^^^^^

::

   repd nrep 4
   stre inp

This will run for independent CHARMM jobs. Each inp_0, inp_1, inp_3,
and inp_4 can be different input files, with different PSFs,
parameters, etc

Example 4
^^^^^^^^^

::

   open write unit 1 card name out

   repd nrep 4
   outu 1
   stre inp

The same as example 3 but now also output files out_0, out_1, ... will
be written. Note that OUTU must precede STREam command.

Example 5: RXSGLD
^^^^^^^^^^^^^^^^^

::

    read psf
    read coor name system.crd

    !All stages have the same temperature of 300 K but have TEMPSG from 300 K to 500 K.
    repd nrep 8 EXCHange FREQuency 1000 STEMp 300 DTEMp 0  -
      SGLD  SGTE 300 MSGT 500  DSGF 0.2

    SCAL FBETA SET 1.0 SELE ALL END

    !Perform SGLD with SGFT set to 0 to allow above RXSGLD setting in control
    DYNA LANG SGLD SGFT 0


.. _repdstr_output:

Output
------

The replica exchange printout is written to the unit specified in
the command after the UNIT keyword. The output of the current results
is labeled by either REX> for temperature based replica exchange or
RXSG> for self-guding replica exchange (RXSGLD). The RXSG> line contains the
following fields:

::

    RXSG> Exchanges DynSteps StagID NeighborID ReplicaID Ep EpNeighbor
    TempScale TSGScale AcceptRatio Acceptance

A summary from all other
replicas is labeled by REXSUM>. The labels in the output are
shortened, end the meaning of some of them is as the following:

========= =====================================================================
Epot      potential energy (current)

Tscale    temperature scaling for exchange [Tscale=sqrt(Temp/NewTemp)]

Sratio    success ratio [Srate=#-of-successful-exchanges/#-of-tried-exchanges]

NewTemp   new temperature after the exchange

CurrTemp  current temperature

PROB      probability to perform exchange P=exp(-Delta(1/kT)*Delta(Epot))

Rand      random number used for exchange condition PROB>Rand => Success=T

NEIGHBOR  current neighbor with which the exchange occurs (or not)
========= =====================================================================



