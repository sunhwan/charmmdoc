.. py:module:: openmm

===========================================
OpenMM GPU acceleration interface to CHARMM
===========================================

This module describes the interface of CHARMM with the OpenMM
development platform for GPU accelerated simulations. The current
interface supports molecular dynamics on CUDA or OPENCL supported
graphical processing units (GPUs). For a full list of hardware on
which the OpenMM libraries should run, see the OpenMM website
(https://simtk.org/home/openmm). The OpenMM libraries are free and
available as pre-compiled libraries or source form. In addition, one
needs the latest NVIDIA drivers and CUDA toolkit installed on the
machine - please see OpenMM documentation for the basic procedures to
set-up and install these components.

The CHARMM/OpenMM interface is under continuing development with
new CHARMM features being added all of the time. The current
implementation supports dynamics and energy calculations for periodic
and non-periodic systems using cutoffs, nocutoffs (for finite
systems), and PME/Ewald and cutoffs for periodic systems.  Periodic
systems supported are only orthorhombic (a,b,c, alpha=beta=gamma=90.
Only Leapfrog Verlet integration and Langevin dynamics are
supported. Constant temperature molecular dynamics is also supported
through Andersen heatbath method in the OpenMM module. Additionally,
constant pressure, constant temperature dynamics are available using
MC sampled barostat. Finally, we have provided access to the variable
timestep Verlet (Leapfrog) and Langevin integrators implemented in the
OpenMM module.  SHAKE is supported as are all of the CHARMM
forcefields, e.g., CMAP.

Special Notice: The CHARMM/OpenMM interface is an evolving interface
with the OpenMM accelerated dynamics engine for GPU accelerated
molecular dynamics (See the News at www.charmm.org for a discussion of
the benchmarks and their performance). The functionality present
through the current CHARMM interface has been released prior to
"aging" in the CHARMM developmental version for a year because of the
important performance enhancements in provides through GPU
acceleration. The interface and associated modules have been well
tested, but are likely to contain yet undiscovered limitations,
compared with the full functionality in CHARMM. Additionally, we note
that this code operates in single precision, which may not be
acceptable for all applications. At present, the CHARMM/OpenMM
interface accommodates molecular dynamics with and without periodic
boundary conditions, using all of the current CHARMM force fields, and
in NVE, NVT and NPT ensembles - although not using the same methods of
achieving these as in the rest of CHARMM. Users are forewarned to
carry out some pre-testing on their system prior to initiating long
runs on GPUs. As new features and methods are added to the
CHARMM/OpenMM interface they will be described in the updated
documentation.

.. _openmm_setup:

SETTING UP AND BUILDING CHARMM/OPENMM
=====================================

To build CHARMM with the OpenMM interface to enable GPU accelerated molecular
dynamics one needs to first install the appropriate GPU drivers and software
support, e.g., NVIDIA drivers and CUDA toolkit. Additionally, the OpenMM
libraries need to be installed. Please see the OpenMM web pages and
documentation to accomplish this procedure (https://simtk.org/home/openmm).

Setup necesssary environment variables for load library and OpenMM path:
The following environment variables need to be setup OPENMM_PLUGIN_DIR
and the library path. Assuming OpenMM has been installed in it's default
directory (/usr/local/openmm), then set the following environment variable:

Mac OSX or Linux

bash shell: ```export OPENMM_PLUGIN_DIR=/usr/local/openmm/lib/plugins```
csh shell: ```setenv OPENMM_PLUGIN_DIR /usr/local/openmm/lib/plugins```

These should be added to your .bashrc or .cshrc to ensure they are always
setup.

Aditionally, one needs to tell the loader where the OpenMM libraries are
installed. This differs on Linux versus Mac OSX systems because of static
versus dynamic load library uses on these two OSs.

Linux (bash): ```export LD_LIBRARY_PATH=/usr/local/openmm/lib:$LD_LIBRARY_PATH```
Linux (csh): ```setenv LD_LIBRARY_PATH /usr/local/openmm/lib:$LD_LIBRARY_PATH```

Mac OSX (bash): ```export DYLD_LIBRARY_PATH=/usr/local/openmm/lib:$DYLD_LIBRARY_PATH```
Max OSX (csh): ```setenv DYLD_LIBRARY_PATH /usr/local/openmm/lib:$DYLD_LIBRARY_PATH```

_openmm_setup:

INSTALLING CHARMM/OpenMM
========================

Installing CHARMM with the OpenMM interface is straightforward on Linux and
Mac OSX:

Linux / Intel compilers

::

  install.com em64t <size> openmm

Linux / GCC compilers

::

  install.com gnu <size> openmm

Mac OSX / Intel compilers

::

  install.com osx <size> ifort openmm

Mac OSX / GCC compilers

::

  install.com osx <size> gfortran openmm

_openmm_usage:

USAGE and IMPLEMENTATION
========================

USAGE: add keyword omm to dynamics command *note dynamc:(chmdoc/dynamc.doc), or to energy call
(via the energy or gete commands *note energy:(chmdoc/energy.doc) . For dynamics, this
gives you the default Verlet Leapfrog integrator with timestep specified in the dynamics
command. For energy/gete calls, you get all energy terms that are active with the non-bonded and
bonded terms being computed on the GPU. One can also include the various options noted below:

SUMMARY OF OPENMM COMMANDS
==========================

::

  OMM [ openmm-control-spec ]

  openmm-control-spec
        on                 Sets omm_active to true and tells CHARMM all subsequent calls to energy,
                           dynamics or minimization will use OpenMM interface for calculation of
                           supported energies and forces. OpenMM context will be created later as needed
       off                 Sets omm_active to false but retains any OpenMM context alread created
       clear               Sets omm_ative to false and destroys the OpenMM Context

Dynamics keyword options in CHARMM/OpenMM interface

================= ==================== =================================================
keyword           default                          action
================= ==================== =================================================
omm                false               dynamics keyword to access openMM interface
langevin           false               dynamics keyword to turn on Langevin integration
andersen           false               dynamics keyword to turn on Andersen heatbath
prmc               false               dynamics keyword to turn on MC barostat
variable           false               dynamics keyword to use variable timestep md
gamma               5.0                Langevin friction coefficient in ps^-1
colfrq             1000                Andersen heatbath coupling constant ps^-1
pref                1.0                MC barostat reference pressure in atmospheres
iprsfrq             25                 MC barostat sampling frequency
vtol               1e-3                Variable timestep error tolerence
================= ==================== =================================================

.. note:: Coordinates, velocities and restart files can be written every NSAVC, NSAVV, ISVFRQ
   timesteps to files specified by IUNCRD, IUNVEL and IUNWRI. Restarts can be used specifying
   RESTSRT in the dynamics command with IUNREA also specified, like normal CHARMM runs.

.. warning:: At present the energy file is not written, since OpenMM only returns the total energies
   (TOTE, TOTKE, EPOT and TEMP) and VOLUME.

Constant Temperature Dynamics
=============================

::

  omm langevin gamma <real>      - runs Langevin dynamics with friction coefficient gamma (ps^1)
                     <5.0>         at a temperature given by finalt in dynamics command.

  omm andersen colfreq <integer> - runs constant T with Andersen collision frequency colfrq
                       <1000>      at temperature given by finalt in dynamics command

Constant Pressure/Constant Temperature Dynamics
===============================================

Using either of the integrators noted above, one can run MC barostat-ed molecular dynamics by
adding:

::

  omm langevin gamma <real> mcpr pref <real> iprsfrq <integer> - runs Langevin dynamics with
                      <5.0>            <1.0>          <25>       barostat with a reference
                                                                 pressure of pref atmospheres
                                                                 and MC volume move attempted
                                                                 every iprsfq steps.

Variable Timestep Molecular Dynamics
====================================

OpenMM has implemented a bounded error estimate driven variable timestep integration scheme
in which the size of the timestep is bounded by a specified error that would be associated with
the explicit Euler integrator. The timestep is chosen to satisfy the following relationship

::

  error = dt^2 Sum_i ( |f_i|/m_i ),

where error is the desired maximum error in the step, given the current forces. From the
user-supplied error, the timestep follows from

::

  dt = sqrt( error / Sum_i ( |f_i|/m_i ) )

Adding variable_timestep vtol <real> (default 1.0e-3) uses a variable timestep version of the
above integrators (Langevin or Leapfrog). One can run NVE dynamics with Leapfrog as well, but
this may be not useful.

One can also use the variable timestep algorithms with the barostat.

Energy Computations
===================

Energy terms supported through the CHARMM/OpenMM interface for computation
on the GPU include: BOND ANGL UREY DIHE IMPR VDW ELEC IMNB IMEL EWKS EWSE EWEX
and HARM. However, these are returned from the CHARMM/OpenMM interface as just ENER,
i.e., the sum of the components. One can evaluate the individual components through
use of the SKIPE commands.

The CHARMM/OpenMM interface supports a subset of the CONS HARM harmonic restraints.
Speciically, the default ABSOLUTE restraints with XSCALE=YSCALE=ZSCALE=1 are
supported. The COMP, WEIGHT and MASS keywords associated with this restraint are
also supported (see *note cons:(chmdoc/cons.doc) and testcase
c37test/3ala_openmm_restraints.inp)

The CHARMM/OpenMM interface can carry-out energy calculations that combine
the forces for energy terms computed on the GPU (non-bonded (VDW/ELEC) and
bonded (BOND, ANGL, DIHE, IMPHI) with those from other CHARMM functionality.
At present, aside from doing a static energy/force evaluation, one cannot use
CHARMM's minimizers or dynamics methods together with these forces (although
it is planned that we will support this functionality in the future).

.. note:: Due to single precision arithmetic on the GPU, long NVE simulations may
   have an energy drift on the order of 10^-2 * KBOLTZ * T / NDEGF per nanosecond.

As noted in the overview above, the CHARMM/OpenMM supports "no frills"
molecular dynamics for periodic and non-periodic systems. For non-periodic
systems cutoffs and no-cutoffs are supported. For cutoff based methods
a reaction field is utilized. This is also true for periodic systems that
don't employ PME/Ewald methods. The cutoff method is keyed to the value of
the energy-related cutoff CTOFNB. If CTOFNB > 99 it's assumed that no
cutoffs are to be used and OpenMM computes all interactions for non-periodic
systems. If CTOFNB < 99 then the solvent reaction field is used with a
cutoff switch such that the electrostatic energy for atom pair ij, u_ij, is
given by:

::

           q_i*q_j    /  1                    \
  u_ij = ---------- .|  ---- + k_rf*r^2 -c_rf  |
         4*pi*eps_0   \ r_ij                  /

  k_rf = (eps_solvwnt - 1)/(2*eps_solvent+1)/(r_cutoff)^3

  c_rf = (3*eps_solvent)/(2*eps_solvent+1)/(r_cutoff)

where r_cutoff is the cutoff distance (CTOFNB) and eps_solvent is the dielectric
constant of the solvent. If eps_solvent >> 1, this causes the forces to go to
zero at the cutoff.

Ewald and PME-based Ewald are both implemented. With PME-based Ewald the
OpenMM interface employs the cutoff (CUTNB), the box length, and an
estimated error desired for the long-range electrostatic forces to
determine the number of number of grid points for the PME calculations,
FFTX, FFTY and FFTZ. However, to maintain consistency with CHARMM, the
CHARMM/OpenMM interface takes CUTNB, FFTX(Y,Z) and Box_x(y,z) to determine
the error estimate and KAPPA. Thus, KAPPA as set in the CHARMM energy/nonbond
or dynamics command may be over-ridden to ensure that FFTX(Y,Z) is maintained
as requested. The error estimate is defined as delta and is related to KAPPA
via the relaitonship

::

  delta = exp[-(KAPPA*CUTNB)^2)                                   (1)

In the current implementation, this relationship in the form

::

  KAPPA = Sqrt(-ln(2*delta))/CUTNB                                (2)

is combined with

::

  FFTX(Y,Z) = 2*KAPPA*box_x(y,z)/(3*delta^(1/5))                  (3)

to eliminate KAPPA and is solved (via bisection) for a delta value that will
yield the user provided values of FFTX(Y,Z) and KAPPA is determined from
the relationship (2) above.

.. _openmm_multigpu:

Control over the number of GPU devices one uses and the platform for GPU-based
computations is available through environment variables

==================== ====================== =================================
Environment variable Setting                Effect
==================== ====================== =================================
OPENMM_DEVICE        0/1/0,1                Use device 0/1/0 and 1 (parallel)
OPENMM_PLATFORM      OpenCL/Cuda/Reference  Use OpenMM platform based on
                                            OpenCL/Cuda/Reference(CPU)
==================== ====================== =================================

.. note:: OpenMM chooses a default platform based on a guess for best performance if not is
   specified with the environment variable. Also, the platform Reference is a cpu-based
   platform for testing/validation purposes.

Example (C-shell) ```setenv OPENMM_DEVICE 0,1    # Use both GPU devices```

.. note:: Support for parallel calculations (0,1) are only supported for platform OpenCL.

.. openmm_examples:

EXAMPLES
========

Molecular dynamics using NVE with PME in a cubic system (from JACS Benchmark):

::

  set nsteps = 1000
  set cutoff = 11
  set ctofnb = 8
  set ctonnb = 7.5
  set kappa = 0.3308  ! Consistent with cutofnb and fftx,y,z
  calc cutim = @cutoff

  ! Dimension of a box
  set size 62.23
  set  theta = 90.0
  ! Dimension of a box
  Crystal define cubic @size @size @size @theta @theta @theta
  crystal build cutoff @cutim noper 0

  image byseg xcen 0.0 ycen 0.0 zcen 0.0 select segid 5dfr end
  image byres xcen 0.0 ycen 0.0 zcen 0.0 select segid wat end

  !  turn on faster options and set-up SHAKE
  faster on
  energy  eps 1.0 cutnb @cutoff cutim @cutim -
          ctofnb @ctofnb ctonnb @ctonnb vswi -
          ewald kappa @kappa pme order 4 fftx 64 ffty 64 fftz 64

  shake fast bonh tol 1.0e-8 para

  set echeck = echeck -1

  open unit 20 write form restart.res

  ! Run NVE dynamics, write restart file
  calc nwrite = int ( @nsteps / 10 )
  ! Run dynamics in periodic box
  dynamics leap start timestep 0.002 -
       nstep @nsteps nprint @nwrite iprfrq @nwrite isvfrq @nsteps iunwri 20 -
       firstt 298 finalt 298  -
       ichecw 0 ihtfrq 0 ieqfrq 0 -
       iasors 1 iasvel 1 iscvel 0  -
       ilbfrq 0 inbfrq -1 imgfrq -1 @echeck bycb -
       eps 1.0 cutnb @cutoff cutim @cutim ctofnb @ctofnb ctonnb @ctonnb vswi -
       ewald kappa @kappa pme order 4 fftx 64 ffty 64 fftz 64 ntrfq @nsteps - !PME
       omm ! Just turn on openMM, get Leapfrog Verlet, NVE

  ! Restart dynamics from current file
  ! Run dynamics in periodic box
  dynamics leap restart timestep 0.002 -
       nstep @nsteps nprint @nwrite iprfrq @nwrite isvfrq @nsteps iunwri 20 iunrea 20 -
       firstt 298 finalt 298  -
       ichecw 0 ihtfrq 0 ieqfrq 0 -
       iasors 1 iasvel 1 iscvel 0  -
       ilbfrq 0 inbfrq -1 imgfrq -1 @echeck bycb -
       eps 1.0 cutnb @cutoff cutim @cutim ctofnb @ctofnb ctonnb @ctonnb vswi -
       ewald kappa @kappa pme order 4 fftx 64 ffty 64 fftz 64 ntrfq @nsteps - !PME
       omm ! Just turn on openMM, get Leapfrog Verlet, NVE

  !!!!!!!!!!!!!!!!!!!LANGEVIN HEATBATH NVT!!!!!!!!!!!!!!!!!!!!!
  ! Run NVT dynamics with Langevin heatbath, gamma = 10 ps^-1
  ! Run dynamics in periodic box
  dynamics leap start timestep 0.002 -
       nstep @nsteps nprint @nwrite iprfrq @nwrite isvfrq @nsteps iunwri 20 -
       firstt 298 finalt 298  -
       ichecw 0 ihtfrq 0 ieqfrq 0 -
       iasors 1 iasvel 1 iscvel 0  -
       ilbfrq 0 inbfrq -1 imgfrq -1 @echeck bycb -
       eps 1.0 cutnb @cutoff cutim @cutim ctofnb @ctofnb ctonnb @ctonnb vswi -
       ewald kappa @kappa pme order 4 fftx 64 ffty 64 fftz 64 ntrfq @nsteps - !PME
       omm langevin gamma 10 ! turn on openmm, set-up Langevin

  ! Run variable timestep Langevin dynamics with error tolerance of 3e-3
  ! Run dynamics in periodic box
  dynamics leap start timestep 0.002 -
       nstep @nsteps nprint @nwrite iprfrq @nwrite isvfrq @nsteps iunwri 20 -
       firstt 298 finalt 298  -
       ichecw 0 ihtfrq 0 ieqfrq 0 -
       iasors 1 iasvel 1 iscvel 0  -
       ilbfrq 0 inbfrq -1 imgfrq -1 @echeck bycb -
       eps 1.0 cutnb @cutoff cutim @cutim ctofnb @ctofnb ctonnb @ctonnb vswi -
       ewald kappa @kappa pme order 4 fftx 64 ffty 64 fftz 64 ntrfq @nsteps - !PME
       omm langevin gamma 10 variable vtol 3e-3 ! turn on openmm, set-up variable
                                                ! timestep Langevin dynamics

  !!!!!!!!!!!!!!!!!!!LANGEVIN HEATBATH/MC BAROSTAT NPT!!!!!!!!!!!!!!!!!!!!!
  ! Run NPT dynamics with Langevin heatbath, gamma = 10 ps^-1
  ! Run dynamics in periodic box
  dynamics leap start timestep 0.002 -
       nstep @nsteps nprint @nwrite iprfrq @nwrite isvfrq @nsteps iunwri 20 -
       firstt 298 finalt 298  -
       ichecw 0 ihtfrq 0 ieqfrq 0 -
       iasors 1 iasvel 1 iscvel 0  -
       ilbfrq 0 inbfrq -1 imgfrq -1 @echeck bycb -
       eps 1.0 cutnb @cutoff cutim @cutim ctofnb @ctofnb ctonnb @ctonnb vswi -
       ewald kappa @kappa pme order 4 fftx 64 ffty 64 fftz 64 ntrfq @nsteps - !PME
       omm langevin gamma 10 - ! turn on openmm, set-up Langevin
       mcpr pref 1 iprsfrq 25  ! set-up MC barostat at 1 atm, move attempt / 25 steps

  ! Run variable timestep Langevin dynamics with error tolerance of 3e-3
  ! Run dynamics in periodic box
  dynamics leap start timestep 0.002 -
       nstep @nsteps nprint @nwrite iprfrq @nwrite isvfrq @nsteps iunwri 20 -
       firstt 298 finalt 298  -
       ichecw 0 ihtfrq 0 ieqfrq 0 -
       iasors 1 iasvel 1 iscvel 0  -
       ilbfrq 0 inbfrq -1 imgfrq -1 @echeck bycb -
       eps 1.0 cutnb @cutoff cutim @cutim ctofnb @ctofnb ctonnb @ctonnb vswi -
       ewald kappa @kappa pme order 4 fftx 64 ffty 64 fftz 64 ntrfq @nsteps - !PME
       omm langevin gamma 10 variable vtol 3e-3 - ! turn on openmm, set-up variable
       -                                          ! timestep Langevin dynamics
       mcpr pref 1 iprsfrq 25                     ! set-up MC barostat at 1 atm, move attempt / 25 steps

  !!!!!!!!!!!!!!!!!!!ANDERSEN HEATBATH NVT!!!!!!!!!!!!!!!!!!!!!
  ! Run NVT dynamics with Andersen heatbath, collision frequency = 250
  ! Run dynamics in periodic box
  dynamics leap start timestep 0.002 -
       nstep @nsteps nprint @nwrite iprfrq @nwrite isvfrq @nsteps iunwri 20 -
       firstt 298 finalt 298  -
       ichecw 0 ihtfrq 0 ieqfrq 0 -
       iasors 1 iasvel 1 iscvel 0  -
       ilbfrq 0 inbfrq -1 imgfrq -1 @echeck bycb -
       eps 1.0 cutnb @cutoff cutim @cutim ctofnb @ctofnb ctonnb @ctonnb vswi -
       ewald kappa @kappa pme order 4 fftx 64 ffty 64 fftz 64 ntrfq @nsteps - !PME
       omm andersen colfrq 250 - ! turn on openmm, set-up Andersen
       mcpr pref 1 iprsfrq 25    ! set-up MC barostat at 1 atm, move attempt / 25 steps

  ! Run variable timestep Leapfrog w/ Andersen heatbath and error tolerance of 2e-3
  ! Run dynamics in periodic box
  dynamics leap start timestep 0.002 -
       nstep @nsteps nprint @nwrite iprfrq @nwrite isvfrq @nsteps iunwri 20 -
       firstt 298 finalt 298  -
       ichecw 0 ihtfrq 0 ieqfrq 0 -
       iasors 1 iasvel 1 iscvel 0  -
       ilbfrq 0 inbfrq -1 imgfrq -1 @echeck bycb -
       eps 1.0 cutnb @cutoff cutim @cutim ctofnb @ctofnb ctonnb @ctonnb vswi -
       ewald kappa @kappa pme order 4 fftx 64 ffty 64 fftz 64 ntrfq @nsteps - !PME
       omm andersen colfrq 250 variable vtol 3e-3 ! turn on openmm, set-up variable
                                                  ! timestep Langevin dynamics

  !!!!!!!!!!!!!!!!!!!ANDERSEN HEATBATH/MC BAROSTAT NPT!!!!!!!!!!!!!!!!!!!!!
  ! Run NPT dynamics with Andersen heatbath, collision frequency = 250
  ! Run dynamics in periodic box
  dynamics leap start timestep 0.002 -
       nstep @nsteps nprint @nwrite iprfrq @nwrite isvfrq @nsteps iunwri 20 -
       firstt 298 finalt 298  -
       ichecw 0 ihtfrq 0 ieqfrq 0 -
       iasors 1 iasvel 1 iscvel 0  -
       ilbfrq 0 inbfrq -1 imgfrq -1 @echeck bycb -
       eps 1.0 cutnb @cutoff cutim @cutim ctofnb @ctofnb ctonnb @ctonnb vswi -
       ewald kappa @kappa pme order 4 fftx 64 ffty 64 fftz 64 ntrfq @nsteps - !PME
       omm andersen colfrq 250 - ! turn on openmm, set-up Andersen
       mcpr pref 1 iprsfrq 25    ! set-up MC barostat at 1 atm, move attempt / 25 steps

  ! Run variable timestep Leapfrog w/ Andersen heatbath and error tolerance of 2e-3
  ! Run dynamics in periodic box
  dynamics leap start timestep 0.002 -
       nstep @nsteps nprint @nwrite iprfrq @nwrite isvfrq @nsteps iunwri 20 -
       firstt 298 finalt 298  -
       ichecw 0 ihtfrq 0 ieqfrq 0 -
       iasors 1 iasvel 1 iscvel 0  -
       ilbfrq 0 inbfrq -1 imgfrq -1 @echeck bycb -
       eps 1.0 cutnb @cutoff cutim @cutim ctofnb @ctofnb ctonnb @ctonnb vswi -
       ewald kappa @kappa pme order 4 fftx 64 ffty 64 fftz 64 ntrfq @nsteps - !PME
       omm andersen colfrq 250 variable vtol 3e-3 - ! turn on openmm, set-up variable
       -                                            ! timestep Andersen dynamics
       mcpr pref 1 iprsfrq 25                       ! set-up MC barostat at 1 atm, move attempt / 25 steps

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!EXAMPLE ENERGY CALCULATIONS!!!!!!!!!!!!!!!!!!!!!!!!
  ! Use omm on/off/clear to set-up and carry-out energy calculations using CPU and/or GPU

  ! Energy calculation for periodic system use PME on CPU
  energy eps 1.0 cutnb @cutoff cutim @cutim ctofnb @ctofnb ctonnb @ctonnb vswi -
         ewald kappa @kappa pme order 4 fftx 64 ffty 64 fftz 64

  ! Same calculation using GPU throuhg CHARMM/OpenMM interface
  energy eps 1.0 cutnb @cutoff cutim @cutim ctofnb @ctofnb ctonnb @ctonnb vswi -
         ewald kappa @kappa pme order 4 fftx 64 ffty 64 fftz 64 -
         omm

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!EXAMPLE II!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! Energy calculation for periodic system use PME on CPU
  energy eps 1.0 cutnb @cutoff cutim @cutim ctofnb @ctofnb ctonnb @ctonnb vswi -
         ewald kappa @kappa pme order 4 fftx 64 ffty 64 fftz 64

  omm on  ! subsequent invocations of energy will use CHARMM/OpenMM interface
  ! Same calculation using GPU through CHARMM/OpenMM interface
  energy eps 1.0 cutnb @cutoff cutim @cutim ctofnb @ctofnb ctonnb @ctonnb vswi -
         ewald kappa @kappa pme order 4 fftx 64 ffty 64 fftz 64

  omm off  ! turn off use of GPU calculation but leave OpenMM "Context" intact
  ! Energy calculation for periodic system use PME on CPU
  energy eps 1.0 cutnb @cutoff cutim @cutim ctofnb @ctofnb ctonnb @ctonnb vswi -
         ewald kappa @kappa pme order 4 fftx 64 ffty 64 fftz 64

  omm clear   ! Deactivate (until next omm on) calculations using GPU

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TEST CASES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  The relevant test cases for the CHARMM/OpenMM functionality are:

      Test case               Purpose
  omm_nonperiodic.inp  Test and compare CHARMM and CHARMM/OpenMM energy and forces for vacuum system
  omm_periodic.inp     Test and compare CHARMM and CHARMM/OpenMM energy and forces for solvated system
  omm_dynamics.inp     Test CHARMM/OpenMM dynamics with various integrators
  omm_restraints.inp   Test restraint methods between CHARMM and CHARMM/OpenMM for vacuum system
  omm_modpsf.inp       Test whether CHARMM/OpenMM senses psf changes and rebuilds OpenMM context
  omm_nbexcl.inp       Test whether CHARMM/OpenMM handles nb exclusions correctly
  omm_nbfix.inp        Test whether CHARMM/OpenMM handles nbfixes correctly


