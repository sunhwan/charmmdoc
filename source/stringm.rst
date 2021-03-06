.. py:module:: stringm

=================
0-K String Method
=================

V. Ovchinnikov (ovchinnv@MIT.edu),
B.L. Trout &
M. Karplus

.. _stringm_description:

Introduction
------------

The string method is an algorithm for finding paths between different
configurations of a molecular system.  It is a `chain-of-states' method,
similar in principle to the Replica Path and Nudged Elastic Band (NEB)
methods, in which a continuous transition path is `discretized'
into a finite collection of system replicas under the constraint
of approximately equal spacing (in an appropriate RMSD metric) between
adjacent replicas.  The `string' is essentially a discretization strategy
that, in principle, can be applied to finding paths on any energy landscape,
both at zero [1] and at a finite temperature [2].

The current implementation in CHARMM is limited to the string method at
zero temperature (0-K String), although further development is on-going.

Algorithm
^^^^^^^^^

The 0-K string is implemented within the framework of the ensemble module
(:doc:`ensemble`) according to reference [1].

For each replica, the 0-K string algorithm performs steepest descent (SD)
evolution on the potential energy landscape (i.e. minimization) followed
by a collective `reparameterization' to maintain equal spacing between
adjacent replicas.  The reparameterization step guarantees that the
replicas do not slide down into the nearest energy minimum.  The procedure
is iterated until the string no longer moves. At this point it has
converged to the minimum energy path (MEP) to within the discretization
approximation [1].


Remarks on the implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Because the string algorithm requires inter-replica communication only during
the reparameterization step, it is naturally compatible with all of the
features of CHARMM available in the ensemble module, and in particular
with all solvent models (although currently only EEF1 and ACE have been tested)

The number of replicas in the discretized string, however, is equal to the number
of processors in the ensemble (?nensem).  Thus, this implementation
only allows one processor per replica.  The string functionality is accessed
through the ensemble module, as described in the next section.

To build the string method use the E and STRINGM keywords, e.g.

::

  $> ./install.com gnu medium g77 E MPICH STRINGM

.. _stringm_syntax:

Command syntax
--------------

::

  ENSEmble STRIng { STATistics }     [ ENER {<energy_terms>} [ENAM <chracter*>] -
                                     [ RMSD [RNAM <character*>] [RAPP] ] -
                                     [ DELS [DNAM <character*>] [DAPP] ] -
                                     [ ARCL [ANAM <character*>] [AAPP] ]
  		{ REPArameterize } [ ZERO ] [ITER <int>] [DEFI <real>]-
  		                   [ LINEar|CSPLines|BSPLines] [MASS] -
  			           [ ORIE] (atom-selection) [MASS]
                  { MINImize }       {REPF <int>} mini-spec
                                     (*note minmiz:(chmdoc/minmiz.doc))

.. _stringm_function:

Description of commands and options
-----------------------------------

::

  ENSEmble STRIng { STATistics }     [ COUNt  <int> ] -
                                     [ ENER {<energy_terms>} [ENAM <chracter*>] ] END -
                                     [ RMSD [RNAM <character*>] [RAPP] ] -
                                     [ DELS [DNAM <character*>] [DAPP] ] -
                                     [ ARCL [ANAM <character*>] [AAPP] ]

This command sets up options for the output of statistics.  When called with no
arguments after STAT, it causes a instance of statistics output to be written
out.

``[COUNt <int>]`` specifies iteration counted for statistics output.  The default
value is 0.  This number is incremented after every call to statistics, and
is printed in the output files corresponding to RMSD, DELS, and ARCL. It is
also appended to the base file name specified by the ``[ENAM <character*>]``

``[ ENER {<energy_terms>} [ENAM <character*>] ]`` sets up output of energy values
along the path.  At each statistics call, the energy values corresponding
to the terms in {<energy_terms>} will be listed in one file per iteration,
with one line of output per replica. The list {<energy_terms>} can contain
one or more energy substitution keywords described in energy.doc
(*note Energy:(chmdoc/energy.doc)).  The current iteration will be appended
to the base energy file name specified with [ENAM <character*>], followed by
extension '.DAT'.  If ENAM is omitted, energy output is directed to the
output stream.  This subcommand must be terminated by END.

``[ RMSD [RNAM <character*>] [RAPP] ]`` sets up output of RMSD values between the
replica coordinates at the current iteration and the coordinates present in
the comparison set (usually the initial string).  At each call to string
statistics, a line will be added to the file with the name specified in
[RNAM <character*>], with the columns in the file corresponding to different
replicas.  If RAPP is specified, output will be appended to the file.
If RNAM is omitted, RMSD output is directed to the output stream.
This subcommand is useful for gauging convergence of the string.

``[ DELS [DNAM <character*>] [DAPP] ]`` sets up output of RMSD values between the
replica coordinates at the current iteration and those at the previous iteration.
At each call to string statistics, a line will be added to the file with the
name specified with [DNAM <character*>].  If DAPP is specified, output will be
appended to the file. If DNAM is omitted, RMSD output is directed to the output
stream.  This subcommand is useful for gauging the convergence of the string.

``[ ARCL [ANAM <character*>] [AAPP] ]`` sets up output of the distance between
the adjacent replicas (such that their sum yields the string length)
at the current iteration.  At each call to string statistics, a line
will be added to the file with the name specified with [ANAM <character*>].
If AAPP is specified, output will be appended to the file. If ANAM is omitted,
RMSD output is directed to the output stream.  This subcommand is useful for
gauging the convergence of the string.

::

  ENSEmble STRIng { REPArameterize } [ ZERO ] [ITER <int>] [DEFI <real>]-
  		                   [ LINEar|CSPLines|BSPLines] [MASS] -
  			           [ ORIE] (atom-selection) [MASS]

This command sets up options for string reparameterization, which ensures
that the string replicas remain equidistant (in the sense of equal RMSD between
adjacent replicas).  When called with no arguments after REPA, it causes the
string to be reparameterized using the options specified in the most recent
REPA call.

``[ZERO]`` indicates that the the reparameterization will use all Cartesian
coordinates (in correspondence with the zero-temperature method).  This
option is assumed since the current implementation of the string method only
contains the 0-K string method.

``[ITER <int>]`` specifies the maximum number of iterations in the
reparameterization call. The reparameterization algorithm is iterative,
and after each iteration the value d=max RMSD[i,i+1]/RMSD[i,i-1] moves
closer to unity.  The default value is 10.

``[DEFI <real>]`` specifies the maximum allowed RMSD error in the distance between
adjacent replicas. The default value is 1.1.  Thus, by default,
reparameterization iterations will continue until
(1) d=max RMSD[i,i+1]/RMSD[i,i-1] >1.1, or
(2) the maximum number of iterations specified by [iter <int>] is exceeded.

``[ LINEar|CSPLines|BSPLines] [MASS]`` specifies the reparameterization algorithm.
Linear, cubic splines and B-splines are supported.  The default is CSPLines.
Specifying BSPLines will usually produce smoother paths. The preferred
method is LINEar. [MASS] causes the atom coordinates to be mass-weighted
in when the string length is computed.

``[ ORIE] (atom-selection) [MASS]`` specifies that the adjacent string replicas
are to be RMSD-aligned based on the atom selection.  This option should
almost always be present.  [MASS] specifies that the orientation is to
use mass-weighting.  Prior to string reparameterization, replica i
will be rotated/translated such that the RMSD between the orientation sets
of atoms of replicas i and i-1 is minimized.

::

  ENSEmble STRIng { MINImize } {REPF <int>} -
                    mini-spec (*note minmiz:(chmdoc/minmiz.doc))

This command calls string steepest descent dynamics (using the SD minimizer).
{REPF <int>} specifies the reparameterization frequency, i.e. the string is to
be reparameterized after every <int> minimization iterations.  The REPA
command must be called prior to string minimization to set up reparameterization
options.  Statistics output will be called before each reparameterization.
The STAT command must be called prior to string minimization to set up
statistics output.

Any options to the minimizer are passed via mini-spec.  Note that only the SD
minimizer is supported.

.. _stringm_examples:

Usage example (assumes that each replica has a defined set of coordinates)
--------------------------------------------------------------------------

::

  ! 1) set up statistics
  ensemble string ener ener bonds angles dihe impr vdw elec enam @outdir/string end -
                  rmsd rname @outdir/rmsd.dat dels dname @outdir/dsdt.dat arcl aname @outdir/arc.dat

  ! 2) setup reparameterization
  ensemble string repa zero iter 10 defi 1.021 linear orie select all end

  ! 3) minimize string -- 2019 iterations with reparameterization after every 20
  ! iterations
  ensemble string mini repf 20 nstep 2019
  ! decrease the minimization step and minimize again
  ensemble string mini repf 20 nstep 2019 step 0.001

Testcase: c35test/zts_ens.inp

.. _stringm_references:

References
----------

[1] E, W., Ren, W. & Vanden-Eijnden, E. 2007.
    Simplified and improved string method for computing the minimum energy paths
    in barrier-crossing events. J. Chem. Phys. 126, 164103-164103-8

[2] Maragliano, L., Fischer, A., Vanden-Eijnden & Ciccotti, G. 2006
    String method in collective variables: Minimum free energy paths and
    isocommittor surfaces. J. Chem. Phys. 125, 024106
