.. py:module: eds

=======================================
Enveloping Distribution Sampling Method
=======================================

| Method development by Christ and van Gunsteren
| Implemented in CHARMM by Tim Miller, Gerhard Koenig, and Bernard R. Brooks

If you use this code, please cite the following works:

for EDS--

* Christ CD, van Gunsteren WF. J. Chem. Phys. 126: 184110 (2007).

* Christ CD, van Gunsteren WF. J. Comp. Chem. 30: 1664-1679 (2009).

* Riniker S, Christ CD, Hansen N, Mark AE, Nair PC, van Gunsteren WF.
  J. Chem. Phys. 135: 24105 (2011).

for MSCALE--

* Woodcock HL, Miller BT, Hodoscek M, Okur A, Larkin JD, Ponder JW,
  Brooks BR. J. Chem. Theo. and Comp., 2011, Vol 7, 1208-1219.

.. _eds_syntax:

Syntax
------

EDS allows the calculation of free energy differences between multiple end
states from a single molecular dynamics simulation of a reference state
R. The reference state is designed to contain the important parts of phase
space for all end states, thus improving the convergence of one step free
energy calculations.

EDS is designed to work with the MSCALE command in CHARMM (see mscale.doc for
further information). The TERM key words in the EDS command must refer to
previously defined MSCALE subsystems. The various MSCALE subsystems represent
end points for the EDS free energy calculation. The free energy differences
must be calculated from the individual energies of the subsystems and the
total EDS energy. These can be stored conveniently in an ASCII formatted file
if the UPEN argument to MSCALE is used.

To run EDS, the MSCALE subsystems should be set up, followed by the EDS and
DYNAmics commands. An example is given below.

::

  EDS TEMP <real> NEDS <int> eds-term-spec

  eds-term-spec ::= TERM <string> <float> ... TERM <string> <float>

  The meaning of the key words is as follows:


====================== =======================================================
TEMP <real>            The temperature at which the EDS simulation is
                       being run. Currently, constant temperature dynamics
                       should be used with EDS. Future enhancements will
                       allow for NVE ensembles to be used.
NEDS <int>             The number of EDS end points that will be used. This
                       should correspond to the number of MSCALE subsystems.

TERM <string> <float>  One "TERM" keyword must be given for each end point.
                       The <string> is the name of the MSCALE subsystem
                       corresponding to that end point, and the <float> is
                       the EDS free energy offset to be used for that
                       end point (refer to the EDS equation in the reference
                       for more details).
====================== =======================================================


.. _eds_examples:

Examples
--------

The following shows how to set up an EDS calculation with two endpoints:

::

  ! first set up both subsystems using MSCALE, the subsystem energies and
  ! the EDS energy will be stored in the file named energies.dat.

  open unit 50 write form name energies.dat

  mscale nsubs 2 upen 50

  subs one coef 1.0 prog "./charmm" -
       outp "subsys/s1.out" inpu "subsys/s1.inp" -
       nproc 1 sele all end

  subs two coef 1.0 prog "./charmm" -
       outp "subsys/s2.out" inpu "subsys/s2.inp" -
       nproc 1 sele all end

  end

  ! now use the EDS command to turn on EDS. In this case, the
  ! first end point will have an energy offset of 1, and the
  ! second has an energy offset of 2 (adjust as needed for
  ! your system).
  eds temp 300. neds 2 term one 1.0 term two 2.0

  ! now run dynamics for stastics collection and sampling.
  dyna leap langevin strt -
    timestep 0.001 nstep 10000 nprint 100 -
    iseed 314159 314159 314159 314159 -
    iasors 1 iasvel 1 iscvel 0 -
    firstt 300. finalt 300. tbath 300. tstruct 300. -
    echeck 500.

.. _eds_notes:

Notes
-----

If the energy differences between the end states are large, the end
states with low energies will be sampled preferentially. In extreme
cases, only one end state might be sampled during an EDS simulation. To
avoid that problem, adequate energy offset parameters have to be chosen
in the TERM command. In an ideal setting, the energy offsets should be
good estimates of free energy difference between the corresponding end
state and the EDS reference state. For optimal performance, the
parameters should be updated periodically. A parameter update scheme for
EDS can be found in J. Chem. Phys. 135, 024105 (2011).

In order to avoid numerical instabilities of EDS, the fluctuations of
the potential energies in each subsystem should be minimized. This can
be achieved by calculating the energy contributions of the common
environment in the main script, while calculating the contributions of
the mutation sites in the subsystems.
