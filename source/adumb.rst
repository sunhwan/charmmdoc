.. py:module:: adumb

=================================
Adaptive Umbrella Sampling Module 
=================================

Setting up of adaptive umbrella potentials. Currently supported types
of umbrella potentials are functions of dihedral angles and functions of the 
potential energy of the system (energy sampling).

.. warning::
   The module is still being developed and some details are likely
   to change in future versions.
   Please report problems to Christian Bartels at cb@brel.u-strasbg.fr

REFERENCES:
* C. Bartels & M. Karplus, J. Comp. Chem. 18 (1997) 1450-
* C. Bartels & M. Karplus, J. Phys. Chem. 102 (1998) 865-
* M. Schaefer, C. Bartels, & M. Karplus, J. Mol. Biol. (1998)


.. _adumb_syntax:

Syntax
------

::

   ADUMb CORR    DIST  UNIT int  SELE...END SELE...END (atom selection x 2)

         CORR    RMSD  COR1 
                       COR2

         CORR    RMSD  SETUp NATOms int NSTRuctures int 

         CORR    RMSD  UNT1 1 int  UNT2 int (atom selection x 3) -
   	       ORIEnt SYMMetry 4X(atom-spec) FOLD int
	      
   ADUMb DIHE    NRES int  TRIG int  POLY int  4X(atom-spec)

   ADUMb ENER    NRES int   TRIG int   POLY int
                 MAXE real  MINE real  [MAXT real] [MINT real]

   ADUMb INIT    NSIM int  [UPDA int] [EQUI int]  [TEMP real] 
                 [AGIN real] [NEXT int] [THRE real]
                 [UCUN int]  [WUNI int] [RUNI int] [FREQ int] 
                 [WCUN int]  [RCUN int] [CPFR int]

   ADUMb PROB    UCUN int [TEMP real] [PUNI int] [TUNI int]

   ADUMb STON

   ADUMb STOFf

   where:      atom-spec ::= { segid resid iupac }
                             { resnumber iupac   }


.. _adumb_function:

Introduction
------------

The module provides commands to define degrees of freedom along which 
adaptive umbrella potentials are applied in molecular dynamics 
simulations. Statistics on the sampling of the degrees of freedom are recorded
during the md simulations and periodically used to update the umbrella
potential such that uniform sampling of the degrees of freedom can be
expected. Currently, dihedral angles and the potential energy are supported
as degrees of freedom.

If several degrees of freedom are defined, multidimensional adaptive
umbrella sampling is performed.

Two sorts of input/output files are used by the module. The "umbrella" 
files contain the umbrella potentials that were used in the simulations 
together with the statistics of the sampling of the bins during the 
simulations. Based on this information the potential of mean force can be 
calculated and the umbrella potential expected to lead to uniform sampling 
can be determined. The second sort of files contains the values of the 
umbrella coordinates (=degree of freedom for adaptive umbrella sampling) 
for each time step in which coordinates were saved to the trajectory files. 
The umbrella coordinates are normalized to the range 0 to 1, independent of 
the degrees of freedom used. From the umbrella coordinates saved, weighting 
factors can be calculated which are needed to calculate average properties of 
the unbiased system.

ADUMb DIHE 
^^^^^^^^^^

Define a dihedral angle as degree of freedom for adaptive umbrella
sampling. To record the statistics the degree of freedom is partitioned
into NRES bins. The umbrella potentials are represented as a linear combination
of two times TRIG trigonometric functions and polynomial functions of degree 
0 to POLY - 1. Repeating the command results in a multidimensional adaptive
umbrella potential.

The coordinates written to the umbrella coordinates file are normalized 
to the range 0 to 1 with 0 corresponding to -180 degrees and 1 corresponding
to +180 degrees.

ADUMb ENER
^^^^^^^^^^

Define the potential energy as degree of freedom for adaptive umbrella
sampling. NRES, TRIG and POLY have the same meaning as in ADUMb DIHE.
MINE and MAXE specify the potential energy range: Statistics on the sampling 
are recorded in the range MINE-0.5*(MAXE-MINE) to MAXE+0.5*(MAXE-MINE). In 
the range outside of MINE to MAXE the umbrella potential is kept constant 
to prevent the system from leaving the range in which statistics are recorded. 
MINT and MAXT (default values: 273 K and 1000 K, respectively) are minimal 
and maximal temperatures to restrict sampling in the relevant temperature 
range. To set up a system, get a rough estimate of the potential energy of the 
system at the desired TMIN and TMAX (from short unbiased simulations at 
TMIN and TMAX).  Set EMIN and EMAX to the values determined minus/plus a 
small tolerance, respectively.

The coordinates written to the umbrella coordinates file are normalized 
to the range 0 to 1 with 0 corresponding to MINE-0.5*(MAXE-MINE) and 1
corresponding to MAXE+0.5*(MAXE-MINE).


ADUMb INIT
^^^^^^^^^^

Defines or redefines the parameters for adaptive umbrella sampling and
initializes the umbrella potential. The umbrella potential is updated
every UPDAte steps. After each update, no statistics are recorded for
EQUI steps. For the remaining UPDA - EQUI steps, statistics on the sampling
of the umbrella coordinates are recorded and stored separately from previous
statistics and together with the umbrella potential active when recording the 
statistics. NSIM separate statistics can be kept in memory. If the number of 
updates performed in a run exceeds NSIM, the oldest statistics are discarded 
to make space for the most recent statistics.

After each update the umbrella potential and the statistics are written
to standard output (the log file). The written table contains, from left
to right, the number of the bin, the number of integration time steps in
which the system was in the bin since the last update, the potential of mean
force calculated with the WHAM equations, the negative of the updated
umbrella potential (potential of mean force modified to restrict sampling if 
necessary and fitted to the set of trigonometric and polynomial functions),
the total number of times the bin was visited in the entire simulation, and
the umbrella coordinates of the center of the bin.

The temperature TEMP should be set to the temperature used in 
the simulations. It is used to calculate the umbrella potentials from
the sampling statistics and to restrict sampling if potential energy 
sampling is performed. 

Umbrella coordinates are written to unit UCUN. At each update,
the statistics are written to unit WUNI together with the umbrella potential 
active when recording the statistics. Statistics from previous runs can
be read from unit RUNI. The statistics read must be from adaptive umbrella
sampling simulations with the same parameters as the present one, in 
particular, the same degrees of freedom have to be used as umbrella 
coordinates. If adaptive umbrella sampling of the potential energy is 
used, umbrella potentials from runs at different temperatures can
be read by repeating the ADUMb INIT command with RUNI set to the unit
containing the statistics of each of the runs and TEMP set to the temperature
of the run.

To define the umbrella potential of bins for which no statistics 
have been acquired so far, the umbrella potential has to be extrapolated.
In the current implementation (might change in future implementations),
the umbrella potential of the bins that were not sampled is set to 
the same value (ext-cons). To determine ext-cons, the potential of the bins
that were sampled is linearly extrapolated for NEXT bins, and the maximal
value (max-extrapolated) of the linearly extrapolated potentials is 
determined. Then, the minimal value (min-sampled) of the potentials of the
bins that were sampled is determined and ext-cons is set to min-sampled
or max-extrapolated whatever value is smaller.

A few statistics that differ significantly from the rest of the 
statistics can be due to problems with the convergence caused by the 
extrapolation or due to the occurrence of rare events. In the former case, 
outliers should occur only in the first few simulations and it is advantageous 
to eliminate them. By default, the module eliminates statistics that 
differ from the averaged statistics by THRE times the average deviation. If
one wants to prevent statistics from being eliminated THRE has to be set to 
a value larger than NSIM. At each update, the deviations of the statistics
from the averaged statistics is printed to standard output (log file), e.g.,

::

 0 Deviation of simulation     1 :    0.955    
 0 Deviation of simulation     2 :    0.513E-01
 0 Deviation of simulation     3 :    0.787E-01
 0 Deviation of simulation     4 :    0.292    
 0 Deviation of simulation     5 :    0.170    
 0 Deviation of simulation     6 :    0.201    
 0 Deviation of simulation     7 :    0.933    
 0 Deviation of simulation     8 :    0.208    
 0 Deviation of simulation     9 :    0.270    
 0 Deviation of simulation    10 :    0.131    
 0 Deviation of simulation    11 :    0.394    
 0 Deviation of simulation    12 :     1.52    
 0 Deviation of simulation    13 :    0.969    
 0 Deviation of simulation    14 :    0.502    
 0 Deviation of simulation    15 :     1.47    
 0 Deviation of simulation    16 :     2.97    
-1 Deviation of simulation    17 :     210.    
 0 Deviation of simulation    18 :    0.695E-01
 0 Deviation of simulation    19 :    0.160    
 0 Deviation of simulation    20 :    0.450    

The 0 or -1 on each line indicates whether the statistics of a particular 
simulation are used (0) or were discarded (-1) based on the THRE criterion.

For complex systems, there might exist no umbrella potential that
enables the system to diffuse rapidly along the umbrella coordinate. In 
such cases it has been found to be advantageous to give a higher weight 
to the most recent statistics. This is implemented using the AGINg factor. 
For an umbrella potential calculated from n statistics, the i'th statistics
(i=1,2,..,n) are weighted by AGINg**(n-i).

The FREQ keyword specifies the frequency with which the dynamics 
trajectories are sampled for compilation of the umbrella potential 
statistics.  FREQ 2 for example means that every other point along the trajectory
is sampled.

WCUNit and RCUNit specify the units to which the accumulators for the
correlated structural variables are to be written and read, for the purposes
of restarting trajectories.  The accumulators, along with the updated
average results, will be written every CPFRequency updates of the umbrella
potential (see also ADUMb CORR).  If WCUNit and RCUNit are omitted, no 
writing of the accumulator statistics will be done.

ADUMb PROB
^^^^^^^^^^

Average properties of the unbiased system can be obtained by weighting
the conformations of an adaptive umbrella sampling run by appropriate
factors. The ADUMb PROB command calculates these weighting factors from
the umbrella coordinates read from unit UCUN and writes them to unit PUNI. 
For the command to work the umbrella potentials and statistics from the
run must have been read with the ADUMb INIT command. If the potential
energy was used as umbrella coordinate, the TEMP specifies the temperature
at which properties of the unbiased system should be calculated.



ADUMb [ STON | STOFF ]
^^^^^^^^^^^^^^^^^^^^^^

By default statistics on the sampling of the umbrella coordinates are
recorded in each call to the energy routines. The ADUMb STOFf command
prevents that statistics are recorded. This might be useful when doing
a minimization or running a md simulation with an umbrella potential
that should not change during the simulation.


ADUMb CORRelations
^^^^^^^^^^^^^^^^^^

The CORRelations keyword allows for the running calculation of the
average values of specified structural variables over the course of the 
trajectores as a function of the reaction coordinates.  It is intended 
as a tool for examining correlations between the reaction coordinates
and various other structural variables in the system.  It is currently 
implemented for interatomic distances and substructure rmsd's. The average
values for the specified variables (distances or rmsd's) are written to
a file (or to standard output) every CPFR times the umbrella potential is
updated, where CPFR is a keyword specified in the UMBR INIT command.

Correlated Distances
....................

The UMBR CORR DIST command sets up the calculation of an average inter-
atomic distance, between atoms specified with a double atom selection. 
Only one atom may be specified for each atom selection. The UMBR CORR DIST
command must be given once for each interatomic distance to be calcu-
lated.  The UNIT keyword is followed by the unit number to which the
results are to be written.  If no unit number is specified,  the results
for the correlated distances will be written to standard output. If a
unit number is specified for any distance, they must be specified for
all distances.  The average distance results will be written every 
CPFRequency updates of the umbrella potential (see ADUMb INIT).
Up to 100 distances can be specified.

EXAMPLE:

::

   umbrella corr dist unit 17 sele atom1 end -
     sele  atom2 end

This will result in the calculation of the running average of the
distance between atom1 and atom2.  (The selection of less than or
greater than exactly 2 atoms will result in an error.) 

The output is formatted as follows:

::

   Average vals of distance fr     17 to      6 at step            500
          2        1    -1.00000000     7.25430918
          2        2    -1.00000000     6.89725628
          2        3    -1.00000000     6.69046274
          2        4     6.38194491     6.41586493
          2        5     5.92699253     5.84622204

The first line describes the variable
The first column gives the assigned number of the distance variable.
The second column gives the position of the reaction coordinate
(same as in free energy output).  The third column gives the average
value of the distance over the last trajectory.  The fourth column gives
the cumulative average over all trajectories.  A "-1" value indicates
that the reaction coordinate position has not been visited.

Correlated Substructure RMSD's
..............................

The UMBR CORR RMSD commands allow for the calculation of the running 
average of the rmsd's, as a function of the reaction coordinates,
for specified parts of the system relative to 2 reference structures.

::

   UMBR CORR RMSD COR1 !saves the current coords as reference structure #1.

   UMBR CORR RMSD COR2 !saves the current coords as reference structure #2.

   UMBR CORR RMSD SETUp NATOms int NSTRuctures int WCUNit int RCUNit int 

This command gives the memory specifications, where NATOms is the
total number of atoms that will be selected for all UMBR CORR RMSD
calculations, and NSTRuctures is the number of sets of substructures 
for which RMSD calculations are to be carried out.  WCUNit and RCUNit
specify the units to which the accumulators are to be written/read for
the purposes of restarting trajectories.  The accumulator values will
be written every CPFRequency updates of the umbrella potential (see also
ADUMb INIT).

The above three commands must be invoked prior to the last set of commands
(UMBR CORR RMSD SUBStructure), which specifies the atoms involved in
the rmsd calculations:

::

   UMBR CORR RMSD SUBStructure UNT1 1 int  UNT2 int (atom selection x 3) -
   	       ORIEnt SYMMetry 4X(atom-spec) FOLD int

UNT1 and UNT2 are the unit numbers for the output (average rmsd's 
relative to reference structures 1 and 2, respectively).  If no unit 
numbers are specified, the results are written to standard output.  Unit
numbers must be specified for either all UMBR CORR RMSD SUBS commands or 
none of them (i.e. either all results are written to files or all are 
written to standard output).  

The three atom selections specify the following (in order): 

1) the atoms whose rmsd is to be calculated 
2) the atoms relative to which a reorientation of the
   system is to take place prior to calculation of the rmsd
3) the atoms involved in a symmetry operation that will be 
   done prior to the calculation of the rmsd.

The ORIEnt keyword invokes a reorientation of the system.
The SYMMetry keyword invokes the symmetry operation, which is a dihedral
angle rotation specified by 4 atoms.  The FOLD keyword specifies the 
multiplicity of the symmetry. The final rmsd will be the lowest one
calculated for any of the symmetric positions. (Only 1 symmetry oper-
ation is allowed per rmsd calculation, currently).  Since the positions 
of atoms in this (3rd) selection will be initialized and rebuilt according
to the internal coordinates of the initialized fragment and the cartesian
coordinates of the rest of the structure, care must be taken in the
selection so as to ensure the initialized fragment is not too large.
If the ORIEnt keyword is specified and only one atom selection is given,
the reorientation (as well as the rmsd calculation) will be done relative
to this selection.  If only one or two atom selections are given, no
symmetry operation will occur (irrespective of the presence or absence of
reorientation).  

The UMBR CORR RMSD SUBStructure command must be invoked once for each set of
rmsd substructure calculations to be done during the dynamics.

Example

::

   UMBRELLA CORR RMSD SUBS UNT1 27 UNT2 28 SELE (phe residue) END -
      SELE (phe backbone) END -
      SELE (phe sidechain) END -
      ORIE SYMM 2 CA 2 CB 2 CG 2 CD1 FOLD 2

This command specifies that the rmsd will be calculated relative to 
the "phe residue" atoms.  Reorientation will be done relative to the
"phe backbone" atoms prior to the rmsd calculation.  A 2-fold symmetry 
operation will be carried out involving the "phe sidechain" atoms and
a rotation about the dihedral defined  by 2 CA 2 CB 2 CG 2 CD1. The
rmsd relative to reference structure 1 will be written to unit 27 and 
that relative to reference structure 2 will be written to unit 28.

Example

::

   UMBRELLA CORR RMSD SUBS UNT1 27 UNT2 28 SELE (phe residue) END -
    SELE (phe backbone) END - 
    ORIE 

This will result in the same calculation as above, absent the symmetry
operation.

The output is formatted as follows:

::

   Average RMSDs from Ref #1 for set      35 at step         500000
         35        1    -1.00000000     1.11056063
         35        2    -1.00000000     1.09866706
         35        3    -1.00000000     1.05065449
         35        4    -1.00000000     1.09327534
         35        5    -1.00000000     1.07153876
         35        6    -1.00000000    -1.00000000
         35        7    -1.00000000    -1.00000000

The first column gives the number of the substructure (numbered
serially from 1 with each UMBR CORR RMSD command). The second column
gives the reaction coordinate gridpoint.   The third column gives 
the average rmsd over the last trajectory.  The fourth column gives
the average rmsd over all trajectories.

NOTE that specification of any correlated variables must be followed
by an UMBRella INIT command, prior to the start of dynamics.
In addition, the specification of any correlated variables will reset
the umbrella potential, causing the previously accumulated statistics
to be discarded.  This is to ensure exact correspondence between the
statistical ensembles that are sampled for the free energy surface
and the structural variables.

In adaptive umbrella sampling without structural correlations,
trajectories (sampling runs between updates) that deviate more than a
specified tolerance from the average trajectory are removed from the
statistics.  This filtering feature is disabled when structural
correlations are invoked, due to the large memory requirements.

The "aging" option, whereby older trajectories may be weighted by
the user less heavily than more recent trajectories, is preserved
for the free energy surfaces when structural correlations are invoked,
but the feature is not implemented for the structural correlations,
themselves, again because of large memory requirements. Hence aging
the trajectories may result in free energy surfaces and structural
correlations that are derived from different statistical distributions.
 

.. _adumb_examples:

Examples
--------

This examples are meant to be a partial guide in setting up
an input file for ADUMB. There are three test files, adumb-phichi.inp,
adumb-enum.inp and ace2.inp.

Example (1) 
^^^^^^^^^^^

Set up and run an adaptive umbrella sampling simulation using two dihedral
angles as umbrella coordinates.

::

   ! define the phi and chi1 dihedral angle as the two umbrella coordinates
   umbrella dihe nresol 36 trig  6 poly 1 pept 1 N  pept 1 CA pept 1 CB pept 1 OG1
   umbrella dihe nresol 36 trig  6 poly 1 pept 1 CY  pept 1 N pept 1 CA pept 1 C

   umbrella init nsim 100 update 10000 equi 1000 thresh 10 temp 300 -
                 ucun 10 wuni 11

   ! perform adaptive umbrella sampling md simulation
   dynamics nose tref 300 qref 20 start -
                nstep 20000 timestep 0.001 -
                ihbfrq 0 inbfrq 10  ilbfrq 5 -
                iseed 12 -
                nprint 1000  iprfreq 1000 -
                isvfrq 1000  iunwrite -1 iunread -1 -
                wmin 1.2

Example(2)
^^^^^^^^^^

Set up and run an adaptive umbrella sampling simulation using the potential
energy as umbrella coordinate (=energy sampling, multicanonical simulation,
entropic sampling).

::

   ! set up umbrella; the range of relevant potential energies is assumed to
   ! extend form -50 kcal/mol to 100 kcal/mol. 
   umbrella ener nresol 200 trig 20 poly 5 mine -50 maxe 100.0 mint 280 maxt 2000

   open write formatted   unit 9 name  @9enum.umb
   open write formatted   unit 10 name @9enum.uco
   open write unformatted unit 11 name  @9enum.cor 

   umbrella init nsim 100 update 10000 equi 1000 temp 1000 thres 100 -
                 wuni 9 ucun 10

   ! energy sampling simulation 
   dynamics langevin start -
                nstep 50000 timestep 0.001 -
                inbfrq 10  ilbfrq 10 rbuffer 0.0 tbath 1000 -
                iseed 12 -
                nprint 1000  iprfreq 1000 -
                isvfrq 1000  iunwrite -1 iunread -1 -
                nsavc  100 iuncrd 11  -
                wmin 1.2

Example(3)
^^^^^^^^^^

Determine the weighting factors to calculate properties of the
unbiased system.

::

   ! define the umbrella coordinates
   umbrella ener nresol 200 trig 20 poly 5 mine -50 maxe 100.0 mint 280 maxt 2000
         
   open read formatted    unit 10 name ../scr/@n.umb
   umbrella init nsim 100 update 10000 equi 1000 runi 10 temp 1000 thres 200

   ! translate umbrella coordinates into probability factors at 300K
   open read formatted    unit 11 name ../scr/@n.uco
   open write formatted   unit 12 name ../scr/@nT300K.pfa

   umbrella prob ucun 11 puni 12 temp 300

   ! translate umbrella coordinates into probability factors at 1000K
   open read formatted    unit 11 name ../scr/@n.uco
   open write formatted   unit 12 name ../scr/@nT1000K.pfa

   umbrella prob ucun 11 puni 12 temp 1000

