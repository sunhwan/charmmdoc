.. py:module:: olap

==================================
Dynamic Importance Sampling (DIMS)
==================================

Dynamic Importance Sampling (DIMS) is a method that generates transitions
between a given initial and final state. Typically, those states are
experimental structures in two different functional states. What sets DIMS
apart from other methods is that no reaction coordinate needs to be defined in
advance and that the quality of a transition can be assessed with a score
during the simulation. The theory of the method is described in the articles
by Woolf et al (:ref:`dims_references`).


.. _dims_syntax:

Syntax of DIMS Commands
-----------------------

The main command is DIMS but a few other commands also have DIMS-related
options. Those options are documented here extensively and their main
documentation links here.


Main DIMS command:

::

  DIMS { DBNM  [HARD]         }
               [DCAR]         }
               [TMD]          }  dims-specs bhes-specs atom-selection
               [HALT]         }     [bnm-atom-selection]
  DIMS { HARD                 }
       { DCARtesian           }  dims-common-specs     atom-selection


  DBNM|HARD|DCAR biasing mode ('flavor')


  bhes-specs             Block Normal Mode parameters
                         (see *note BNM:(chmdoc/vibran.doc)Block normal mode
                          method.)

  dims-specs         ::= dims-common-specs [dims-nm-specs]

  atom-selection         atoms to which the DIMS bias is applied;
                         (see *note atom selections:(chmdoc/select.doc).)

  bnm-atom-selection ::= atom-selection
                         (atoms for which the normal modes are computed. This
                         MUST contain all atoms of the protein)


  dims-common-specs  ::= COFF real SCAV int DSUNit int

                        (see *note Common DIMS parameters::)


  dims-nm-specs     ::= DOFF real  DSCAle real
                        NBIAs int   SKIP int    BSKIp int
                        NBESt int  COMBinations int  NWINdow int
                        ORIEnt int
                        NMPRint int UNIT unit
                        TMD HARD DCAR
                        DFIX
                        MTRA int   NMUNit unit

                        (see *note DIMS-NM parameters::)


Extensions to other commands:

::

  ------------------------------------------------------------------------
  READ NM UNIT unit            See also *note READ:(chmdoc/io.doc)Other files.

  ------------------------------------------------------------------------
  COORdinates COPY DIMS        See *note simple coordinate
                               commands:(chmdoc/corman.doc)Simple. for a
                               description of the new DIMS set.

  ------------------------------------------------------------------------
  DYNAMICS ... OMSC            Compute the *note Onsager-Machlup score::
                               during dynamics.

  ?OMSCORE                     After the DIMS run, this energy variable
                               holds the trajectory's *note
                               Onsager-Machlup score::

  ?DSCORE                      DIMS score *note DIMS score::

  ?LDSCORE                     logarithmic DIMS score

  ------------------------------------------------------------------------
  ENTEr <name> OMSC ETA <friction>
                               *note Onsager-Machlup score:: score as a CORREL
                               timeseries. See also
                               *note CORREL OMSC:(chmdoc/correl.doc)Enter.

.. _dims_description:

Description of DIMS Commands
----------------------------

DIMS has a rather large number of options implemented and it can access
several other parameters from other functions that it uses. DIMS can use
different 'flavors' to bias the transition. (We use 'flavors' in favor of
'modes' in order to avoid confusion with normal modes.):

* DBNM

  DIMS-NM or DIMS-Block Normal Modes goes from the origin toward the target
  by displacing atoms on the conformational space, using collective motion
  information---the normal modes---as bias
  (:ref:`Perilla 2007 <dims_references>`).
  This algorithm produces the best transitions but close to the target the
  bias may not be strong enough to reach the target. At this point one can
  have DIMS use a different algorithm to reach the target.

  The algorithm employs the :ref:`vibran_block_normal_mode_method` by :ref:`Li and
  Cui <dims_references>`.

  DIMS defaults to DBNM.

* HARD

  Atoms are pulled from the origin towards the target based on the distance
  and the remaining time steps (:ref:`Woolf 1998 <dims_references>`). It is
  guaranteed to reach the target state in a given number of steps but the
  transitions can necessarily become forced with rather low quality scores.

* DCARtesian

  DIMS-Cartesian accepts moves that go toward the target
  (see :ref:`Zuckerman 2002 <dims_references>`) and uses bias
  moves towards the target with an acceptance function of the form:

  ::

                                            2
    P  (DeltaPsi) = exp[ |DeltaPsi/DeltaPsi|  ]  if DeltaPsi < 0,
     acc                                         otherwise P  (DeltaPsi)=1
                                                            acc

  where P_acc(DeltaPsi) is the selected order parameter (:ref:`dims_rmsd_score`
  or, for instance, the :ref:`dims_interatomic_distance` score.)

  This gently moves the system toward the target without a restriction on
  the total time.  If the barrier height is not high enough under some
  conditions this algorithm will not converge. When the barriers to
  conformation change is small this approach will converge with a better
  DIMS or OM score.

.. _dims_common_dims_parameters:

Common DIMS parameters
^^^^^^^^^^^^^^^^^^^^^^

The common parameters are:

* ``SCAV int``

  Number of dims scores to include in the computation of the scaling
  factor, default: 5. The scaling factor is computed as the average
  of the first SCAV scores, afterwards each new score is multiplied
  by this scaling factor.

* ``DSUNit int``

  Unit number to store DIMS score for the current run. default: -1.
  NOTE: The unit must be open before calling DIMS or Charmm will crash.

* ``atom-selection``

  The bias is applied to the :ref:`selection of atoms <select_syntax>`.
  A sensible choice is all heavy atoms or
  the back bone. Note that an atom selection must be provided or Charmm
  bombs.


.. _dims_dims_nm_parameters:

DIMS-NM parameters
^^^^^^^^^^^^^^^^^^

DIMS-NM is signified by the DBNM keyword. In this mode, the transition is
biased by using a combination of normal modes. The normal modes are computed
using the :ref:`vibran_block_normal_mode_method`.

Near the target configuration the energy landscape becomes rather soft and
normal modes are often not sufficient to drive the transition to the exact
target configuration. For this case DIMS-NM includes the 'Last-Mile-Hard'
option which allows it to exactly reach the target by using the DIMS-hard
mode. The switch to DIMS-hard occurs once the progress threshold COFF has
been reached.

The default score to measure the transition progress is RMS distance to the
target (:ref:`dims_rmsd_score`), although in principle it can use any :ref:`dims_progress_score`.

In order to increase the variety in an ensemble of transition trajectories
DIMS-NM can use the :ref:`dims_mode_self_avoidance` algorithm.

The syntax of the DIMS-NM command is

::

   DIMS DBNM { [HARD] }
             { [DCAR] } dims-specs bhes-specs atom-selection
             { [TMD]  }            [bnm-atom-selection]
             { [HALT] }

The :ref:`vibran_block_normal_mode_method` defines the set of atoms for which the block normal modes are
computed. It defaults to the first selection. Its main use is when
simulations with explicit solvent are performed. In this case the
normal modes should only be computed for the protein (although the
Hessian is built from all interactions, including the solvent). If the
DIMS bias should only be applied to, say, the backbone, then the
bnm-atom-selection MUST contain the whole protein, including all
hydrogens as otherwise the normal modes would be calculated wrongly.


DIMS-NM supports the following options:

* COFF real

  This option tells DIMS-NM when to stop biasing based on the proximity to
  the target (measured by the progress score, which is by default the RMS
  distance from the target in Angstrom. Default value: 1.0

  Depending on the options given, DIMS uses different approaches for the
  remaining steps after the COFF threshold has been reached:

 (no keyword, the default)
    After COFF has been reached, the remaining NSTEP steps will be run
    with unbiased MD. The trajectory is not guaranteed to exactly reach the
    target.

* HALT

  After COFF has been reached the run stops.

* DCARtesian

  DCARtesian accepts moves that go toward the target and uses bias
  moves towards the target with an acceptance function.  Biasing
  on the cartesian coordinates is being done using 'soft-ratcheting;
  it is not guaranteed to reach the target.

* HARD

  The "Last-Mile-Hard" version is used (which is equivalent to running
  'DIMS HARD') and the target will be reached by forcing the atoms to go
  towards the target during the remaining steps. This can be important for
  trajectory annealing schemes.

* TMD

  If targeted molecular dynamics (TMD) is enabled in Charmm and
  the TMD flag is set then a last-mile TMD approach will be
  run. This is equivalent to stopping your simulation at a given
  cutoff and then running regular TMD from that final state
  towards the target. TMD has to be configured via the regular TMD
  commands (see :doc:`TMD`) prior to the DIMS
  call as DIMS does not handle any of the TMD parameters in any
  way except for the target array orientation (if enabled).

* DOFF

  After the cutoff has been reached sometimes the structure tends to go back
  thus increasing the order parameter. 'DOFF ("DynCutOff") prevents this by
  re-computing the collective motions and forcing the structure to stay
  within certain distance to the target. This option must be used with
  caution as it might lead to undesired impulses in the dynamics.

* DSCAle real

  This is the NM-vector scaling. The force of the bias highly depends on
  this parameter. The bias is applied for NBIAs steps. It is gradually
  switched on with a sigmoidal function (over 1/3 NBIAS), set to a constant
  DSCAle for 1/3 NBIAS, and switched off gradually over the remaining 1/3
  NBIAs steps. Reasonable values range from 2.5*10^-2 to 2.5*10^-3

* NBIAs int

  The bias is applied for NBIAs steps.

* SKIP int

  Recompute the normal modes every SKIP steps. This is computationally
  expensive so it is prudent to use a large SKIP value and a small BSKIP
  value (see BSKIp). In this case, SKIP should be a multiple of BSKIp.

* BSKIp int

  The bias is applied every BSKIp steps for the next NBIAs steps. The
  default value of BSKIP is the value of SKIP, which means that by default
  the normal modes are recomputed every BSKIP steps. However, it is more
  efficient (and seems to lead to more natural transitions) to only
  recompute the normal modes every few thousand steps and reuse the same set
  of normal modes for many cycles of biasing and relaxation. For example, if
  SKIP 5000, BSKIP 40, and NBIAS 21, then every 5000 steps the Hessian is
  diagonalized and the normal modes are recomputed. Every 40 steps, the bias
  is applied for 21 steps, then for 19 steps the system evolves without
  bias.

* NBESt int

  make a list of the NBESt "best" normal modes, where "best" means that
  moving the system along this mode improves the progress score (by default
  the RMSD) in the direction of the target structure.

* COMB int

  From the NBEST modes build combinations of up to COMB modes and evaluate
  those combinations. E.g. if COMB 3 then singlets, doublets and triplets of
  modes will be evaluated and ranked.

* ORIEnt int

  Re-orient the target every ORIEnt steps. If set to -1 then the structure
  is not reoriented. Reorientation of the target does not need to be done
  very frequently unless large changes happen quickly. A value of the order
  of 1000...10,000 is probably appropriate.

* NMPR int

  Write the selected normal modes to the unit defined by UNIT every
  NMPR steps.

* UNIT unit

  UNIT number to write the normal modes.

* NWIND int

  In order to generate variety in the transition, avoid the same combination
  of modes as a bias within a window of +/-NWIND steps around the current
  time step. The sequence of modes used must have been saved within a
  previous run using the NMUN keyword and then read with READ NM UNIT unit.

* MTRA int

  MTRA is the number of NM bias-sequences stored in the file read with READ
  NM UNIT unit.

* NMUN unit

  Unit to write the sequence of normal modes used as bias. This is used in
  subsequent runs to avoid re-using the same modes
  ("self-avoidance"). Setting NMUN -1 disables writing of normal mode
  combinations.

* DFIX

  This option enables DynFix which automatically sets to zero the
  contribution to the motion from regular-MD for the steps in which the bias
  from the collective motions is included. The system evolves exclusively
  along the normal modes chosen as bias. NOT RECOMMENDED FOR STANDARD USE.
  Default setting: OFF.

The command

::

   READ NM UNIT unit

reads the sequence of normal mode combinations that were used in previous
DIMS-NM runs. It is used in conjunction with the MTRA, NWIND, and NMUN
keywords to compute an ensemble of trajectories with
:ref:`dims_mode_self_avoidance`.

Also note that DIMS makes use of the Block Normal Mode subroutine implemented
by Dr. Guohui Li. Convergence also depends on those parameters; for further
information please refer to :ref:`BNM <vibran_block_normal_mode_method>`
method and his paper on BNM (:ref:`Li and Cui 2002 <dims_references>`)

The main features of DIMS-NM are described in more detail in their
own entries:


.. _dims_mode_self_avoidance:

Mode self-avoidance in DIMS
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To estimate transition rates a diverse ensemble of trajectories is
required. In order to increase diversity, one can calculate trajectories
sequentially and use information from the previous runs to avoid recreating
very similar trajectories.

We employ an approach inspired by self-avoiding random walks: DIMS-NM can
ignore modes or mode combination that were used in previous run at a given
time step (or window around a time step). Here the assumption is that modes
with the same mode number are the same mode and hence ignoring a given mode
forces the system to evolve in a different direction. Of course, this
assumption is not strictly true. In practice we found that this approach does
lead to an increased spread in trajectory space.

The 'mode self-avoidance' algorithm requires the use of a new files. This file
is used to store the modes that were used during a run. On subsequent runs if
is read with READ NM, and new modes are appended to it. The modes must be read
using the READ NM command,

::

   open read unit 1 card name nmavoid.dat
   read nm unit 1
   close unit 1
   open append card unit 10 name nmavoid.dat
   dims ... nmun 10 ... mtraj 2

The write unit must be passed to DIMS using the NMUN keyword. If a unit other
than -1 is specified then the self-avoidance feature is active.

The MTRA parameter specifies how many trajectories are included in the
file.

DIMS can also avoid normal modes previously used within an specified window
with the NWIND keyword: If the mode (or mode combination) already occurred in
a previous trajectory within +/-NWIND steps of the current step then those
modes are ignored.



.. _dims_mode_combinatorics:

Mode combinatorics in DIMS
^^^^^^^^^^^^^^^^^^^^^^^^^^

By default DIMS-NM will only use one normal mode to bias the
transition. From all NMOD modes it uses the one which results in the
largest change towards the target (measured by the progress score).

However, it can be beneficial to combine normal modes for the biasing
step, say a combination of three modes. This is implemented as
'combinatorial normal mode DIMS' and signified with the parameter COMB
having a value larger than 1. COMB gives the maximum number of modes
to be combined. For instance, if COMB = 3, then at each biasing step
DIMS looks for the best singlet, doublet, or triplet of modes to use as
a bias. To speed up the combinatorical search it is prudent to
restrict the initial mode space from NMOD to the NBESt singlet modes.

Our tests show that the combinatorial version (COMB>1)
gives better transitions than the singlet version (COMB=1).

.. _dims_mode_self_avoidance_file:

The mode self-avoidance file
............................

When using self avoidance DIMS in combination with the combinatorial
version, DIMS will avoid a number COMB of combinations of normal
modes. One can not directly mix normal modes obtained for a simulation
with different COMB values. If this is desired then one will have to
pre-process the input NM file to match the file format expected for
the new COMB value.

For example, suppose a simulation was run with COMB=1 and self
avoidance so the normal modes at each NMPR steps were written to
NMUNit. Since COMB=1 is being used, DIMS will save just one mode for
each step. For the second simulation we want to increase the
combinatorics to three, i.e. COMB=3, but still avoid the modes
previously used.  This is not a straightforward procedure as DIMS will
be expecting to see a file with triplets of nodes instead of singlets
from the first simulation, thus the NM file must be pre-processed
externally. Two skip modes must be added for each step in the
modes file. A skip mode is symbolized by -1. Examples should
make this clear:

Example file for COMB 1 (the ####### symbolizes a new trajectory):

::

   ** TITLE
   ** My normal modes singlets
   *
   2
   33
   21
   #######

Example file for COMB 3 but based on a previous COMB 1 run:

::

   ** TITLE
   ** My modified normal modes -> triplets
   *
   8
   -1
   -1
   33
   -1
   -1
   21
   ########
   -1
   -1

.. _dims_dims_nm_algorithm:

Outline of the DIMS-NM algorithm with mode combinations and mode self-avoidance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DIMS-NM uses normal modes to bias the transition. Two additional
features increase the quality and diversity of trajectories: linear
combinations of normal modes and mode self-avoidance. The two main
ideas are:

(1) To not only use the best normal mode but combinations of modes.
(2) To increase diversity in an ensemble of trajectories by remembering
    which modes were used at or near the same time step in previous
    trajectories and then choosing different modes to bias the dynamics.

The basic algorithm for DIMS-NM that incorporates (1) and (2) reads thus:

1. Every SKIP steps, diagonalize the Hessian and compute the first
   NMOD normal modes.

2. Rank those modes individually by how much they move the current
   structure towards the target; the "best" mode is the one that
   moves the structure closest towards the target as measured by
   the progress variable and is ranked first.

3. Choose the top NBES best modes.

4. From those best modes, compute the possible change in progress
   variable for all combinations of 2, 3, ... up to COMB modes.
   Select that combination of modes that reduces the progress
   variable most.

5. If self-avoidance is selected check if this combination has already been
   used at this step in a previous run (or within a window of +/-NWIN steps
   around the current step). If this is the case forget this combination
   and try the next one (step 4).

6. Apply the bias (scaled by the factor DSCAle) along the mode(s)
   for the next  NBIAS steps during the dynamics.

7. Run unbiased dynamics for the remaining SKIP - NBIAS steps.

8. Check if the RMSD to the target has reached the cut-off distance
   COFF (in Angstrom).

   * If this is the case, switch to the hard DIMS version or TMD-DIMS
     to move directly to the target.


.. _dims_progress_score:

Measures of progress of the transition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To evaluate the conformational transition a progress variable must
be specified. For our current purposes we select an order parameter
based on RMS. Any other variable that attributes the transition
progress can be used.  An example is a set of pairwise distances or
a vector RMS between two structures.

Other measures can be implemented in Charmm. As an example,
Interatomic distance was also implemented.


.. _dims_rmsd_score:

Root mean square distance score
...............................

Root mean square deviation in configuration space (RMSD) is the default
score. It is measured in Angstrom. Small values mean that a given configuration
is close to the target.

For the calculation of the RMSD, the target configuration must be repeatedly
reoriented with respect to the evolving configuration of the molecule in the
trajectory. This is accomplished by setting the ORIEnt parameter to the DIMS
command.

.. _dims_interatomic_distance:

Interatomic distance score
..........................

NOTE: currently not enabled by default, requires ANNLIB.
see :ref:`dims_developer`.

The interatomic distance score is based on the distances between a
predetermined set of atoms during the simulation. Two set of atoms are
required for this approach: The first set, called the data set is
composed of the atoms/points whose distance to a second set, the query
set, is computed. These pair of sets can be equal or have
common elements, there are not restrictions as on what to put on
them. Note that both set of points dynamically evolve during the
simulation, the only things that are kept fixed are the atoms/points
indexes. Let C0 be the set of the k-nearest-neighbors in the data set
for every point in the query set for the target structure. Let I0 be
the atom indexes for the atoms in C0. The atomic distances between the
elements of C0 and the ones from the query set are then stored in D0
for the target structure. After a simulation step new coordinates are
obtained and both the data set and the query set are updated. Now the
set C is the set of atoms that belong to the data set that are indexed
by I0. The distances between the atoms in C and the ones in the query
set are then stored in D. The progress score is then defined as
ABS(D-D0). Therefore the score is positive definite and a lower score
means a closer match and 0 for identical. This approach is called
using the keyword IATD and an integer (k) and the two selection of
atoms:

::

  dims ... IATD 3 sele type N end sele type CA end ! Compute the score from
                      ! the distance between the three nearest N and each CA

The second approach also uses the data and query sets concept, however it
doesn't use the k-NNs and is slightly different from the previous one.
Two molecular descriptors are needed to compute the score. Each descriptor
contains the first three moment of the distance distribution for each query
point. The score is computed as the Manhattan distance between two molecular
descriptors. To use this metric we use the keyword adist plust two selection
of atoms:

::

 dims ... ADIST sele all end sele type CA end ! all-Atom distribution about CAs


.. _dims_trajectory_scores:

Trajectory scores
^^^^^^^^^^^^^^^^^

Trajectory scores are used to rank trajectories in an
ensemble. Higher-ranked ('better') trajectories are the ones which are
more likely to occur without bias.


.. _dims_onsager_machlup_score:

Onsager-Machlup Score
.....................

The Onsager-Machlup score (OM score) is an action, computed along the
whole trajectory. The smaller this number is the more the (biased)
transition resembles a transition that could have naturally occurred.

The step score s(t) is the Onsager-Machlup action for the given time
step,

::

          N_atom / x_i(t) - x_i(t-dt)     F_i  \ 2
  s(t) :=  Sum   |------------------- - ------- |
           i=1   \        dt            m_i*eta/

The cumulative OM score S(t) is

::

           t
  S(t) := Sum   dt * s(t')
          t'=0

The OM score S_OM of a trajectory of length t_traj is the cumulative
OM score of the last frame,

::

  S_OM = S(t_traj)

After the DIMS run, the energy variable ?OMSCORE holds the
trajectory's Onsager-Machlup score. If the trajectory is continued
with another DIMS run, one can simply add the two scores for the score
of the combined trajectory.

With the OMSC keyword to DYNAMics (:ref:`Langevin dynamics <dynamc_syntax>`)
the OM score is computed during the
simulation and printed out. It is the sum of all the step scores over
the whole trajectory so far.

The normalized-cumulative-score

::

  S*(t) := S(t)/s(0)

is computed with the :ref:`OMSC time series <correl_enter>`.

.. _dims_dims_score:

DIMS score
..........

The simplified DIMS score is described in
:ref:`Jang 2006 <dims_references>`. It approximates a rigorous score
based on transition probabilities as a ratio of Boltzmann probabilities of the
unbiased and the biased movement of the system.

The single step score R_i^s at step i is:

::

    s   exp(-DeltaE_Q/kT)
   R  = -----------------
    i   exp(-DeltaE_D/kT)

where DeltaE_Q is the change in total energy if the system evolves unbiased
(according to the distribution Q) and DeltaE_D is the change in energy under
the bias D.

The DIMS score is the product along the trajectory

::

         N      S
   S = Product R
        i=1     i

A trajectory close to a naturally occurring one should have single step scores
close to 1 and a DIMS score close to 1. In practice, the DIMS score almost
always becomes very small and eventually underflows.

The logarithmic DIMS score ln S is slightly more robust in this respect.

In order to calculate reaction rates only relative scores between trajectories
are important so it is only necessary to record a non-vanishing score. Thus,
it is not a problem per se if the score is much smaller than 1.

The final DIMS score is available in the energy variable ?DSCORE and is also
written to the file designated by DSUNIT. The logarithmic DIMS score is
provided in ?LDSCORE. During a run, DIMS writes those score to the standard
output as well (see :ref:`dims_output`).


.. _dims_restriction:

Restrictions when using DIMS
----------------------------

The DIMS code is currently (December 2007) in alpha release. Feedback is very
welcome. See :ref:`dims_references` for contact details.

DIMS

* Mostly tested with DYNAMICS Langevin.
* No good default values - most parameters are probably system dependent.
* DIHEdral bias currently not implemented.
  (see :ref:`Jang 2006 <dims_references>`)
* PSF for initial and final state must be the same.
* DCAR runs in parallel as an MPI version but DBNM or HARD cannot be run
  in parallel.
* The interatomic distances progress score is not enabled by
  default because it requires an additional library.

CORREL OMSC

* OMSC is incompatible with RMS because they both use the same reference
  array for different things (RMS stores the comparison structure, OMSC
  the previous frame to compute velocities X(t) - X(t-1).) To be safe,
  only use a single ENTER name OMSC .. time series per CORREL command.

.. _dims_examples:

DIMS examples
-------------

It is recommended to read through the examples in sequence as important setup
steps are only shown in the first example. Further options are then added in
the other examples.

(See also the scripts in the test directory.)

The basic work flow with DIMS is simple:

1. Load structures; the target is stored in the DIMS coordinates set using
   COOR COPY DIMS

2. Setup DIMS using the DIMS command.

3. Run dynamics.

If the 'mode self-avoidance' is used, also save the modes (see the second
example below) and repeat 3 to generate an ensemble or
trajectories. Alternatively one can also use an ensemble of initial and final
structures (e.g. from short MD) and/or different random seeds for initial
velocity assignment and Langevin random forces.

.. dims_nm_dims_example:

NM DIMS example
^^^^^^^^^^^^^^^

* Example 1a: Normal Mode-biased DIMS

  DIMS with normal mode bias (mode keyword DBNM) tends to produce better
  transitions than the cartesian or dihedral-based schemes. In this example we
  combine up to 3 normal modes (COMB 3) out of the 15 best normal modes (NBES
  15).

  The following Charmm script fragment is not complete but shows the most
  important steps.

  ::

    ! generate psf (same for start and target structure)
    OPEN READ CARD UNIT 1 NAME start.crd
    READ SEQUENCE COOR UNIT 1
    GENERATE prot SETUP FIRST NTER LAST CTER      ! PSF for initial and final
                                                  ! states must be the same

    ! target configuration
    OPEN READ UNIT 1 CARD NAME target.crd
    READ COOR CARD UNIT 1
    CLOSE UNIT 1

    COOR COPY DIMS         ! target must be copied to DIMS set

    ! starting configuration
    OPEN READ UNIT 1 CARD NAME start.crd
    READ COOR CARD UNIT 1
    CLOSE UNIT 1


    DEFINE mydims SELECT all END      ! atoms to apply bias to

    DIMS DBNM  -                      ! set up NM-DIMS
         DSCALE 8e-3 SKIP 500 BSKIP 50 NBIAS 27 SCAV 10 - ! dims options
         SERL GENR SCAL 0.5882 TMEM 420 MEMO 20 MEMA 350 NMOD 50  - ! BNM options
         COFF 0.8 HARD ORIENT 400                         - ! Lastmile / Orient
         COMB 3 NBES 15                                   - ! Combinatorics
         MTRA 0  NMUN -1 NWIND 0 DSUNit 11                - ! Self Avoidance
         SELE mydims END                                    ! DIMS atom selection

    SCALAR fbeta SET 50.0             ! friction coefficient in 1/ps

    ! run Langevin dynamics (in implicit solvent)
    DYNA LEAP LANGEVIN START ... -
         OMSC ... -                   ! compute *note Onsager-Machlup score::


.. _dims_dcar_dims_example:

DCAR DIMS example
^^^^^^^^^^^^^^^^^

* Example 1b: DIMS-CARTESIAN-biased DIMS

  ::

    open read....                 ! Read input similar to that of NM DIMS example

    DEFINE mydims SELECT all END  ! atoms to apply bias to

    DIMS DCAR  1e-7  -            ! set up DCAR-DIMS
         orient 10                ! re-orient at every 10 nstep
         SELE mydims END          ! DIMS atom selection
         COFF 0.6 halt            ! stop biasing when target is 0.6 A from target

    SCALAR fbeta SET 25.0         ! friction coefficient in 1/ps

    ! run Langevin dynamics (in implicit solvent)
    DYNA LEAP LANGEVIN START ... -
         OMSC ... -               ! compute *note Onsager-Machlup score::


.. _dims_nm_dims_self_avoidance_example:

Normal Mode-biased DIMS with self avoidance of modes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DIMS is capable of avoiding previously used normal modes. On the first run,
save the normal modes:

::

  open write card unit 10 name nmavoid.dat
  dims ... nmun 10 ... mtraj 2

For subseqent runs, load the modes and append new ones by
using the READ NM command (:doc:`io`):

::

   open read unit 1 card name nmavoid.dat
   read nm unit 1
   close unit 1
   open append card unit 10 name nmavoid.dat
   dims ... nmun 10 ... mtraj 2

The write unit must be passed to DIMS using the NMUN keyword; using NMUN=-1
does not save the modes. One can also specify how many trajectories are
included in the file with the MTRA keyword. DIMS can also avoid normal modes
previously used within an specified window with the NWIND keyword.

::

  DEFINE mydims SELECT all END    ! atoms to apply bias to

  DIMS DBNM -
    DSCALe 0.1 SKIP 500 BSKIP 50 NBIAs 27                 - ! dims options
    SERL GENR SCAL 0.5882 TMEM 420 MEMO 20 MEMA 400 NMOD 30 - ! BNM options
    COFF 2.0 HARD                                         - ! NM Hard Cutoff
    ORIEnt 20                                             - ! Orient every 20th
    SELE mydims END                                       - ! selection
    COMB 3 NBES 15     -  !  Store 15 best modes and then group them
                       -  !     by 1, 2 and 3 for each try.
    NWINDow 12         -  !  Avoid normal modes previously used within
                       -  !     a window of 12 modes
    MTRA @I NMUNit 10  -  !  @i trajectories per file, write output to unit 10
    DSUNit 11          -  !  Use unit 11 to store dims score


.. _dims_explicit_solvent_example:

Explicit solvent example
^^^^^^^^^^^^^^^^^^^^^^^^

* Example 3: Normal Mode-biased DIMS with explicit solvent

In this example, we will use DIMS on a peptide in water. The system is
simulated using the Spherical Solvent Boundary Potential :doc:`SSBP <mmfp>`.

::

  ! read psf (peptide + water)
  ...

  DEFINE solute SELECT SEGID pept END ! only the peptide
  DEFINE mydims SELECT solute END     ! DIMS on the whole peptide but not water

  ! read target coordinates
  ...
  COOR COPY DIMS  SELECT mydims END ! setup DIMS target structure (no solvent!)

  ! read starting structure
  ...


  ! non bonded interactions
  NBONDS  EXTEND  GRAD   QUAD   GROUP   SWITCH  CDIE      EPS 1.0 -
          VDW          VSWITCH      -
          CUTNB  12.0  CTOFNB 12.0  CTONNB 12.0  WMIN 1.2     WRNMXD 1.2

  !------------------------------------------------------------
  ! SSBP & restraints
  !------------------------------------------------------------
  COOR STAT SELECT solute END
  SET xcen = ?XAVE
  SET ycen = ?YAVE
  SET zcen = ?ZAVE

  MMFP
    ! Use Stochastic Boundary Potential to constrain water
    ! (flexible boundary which adjusts shape). Leave out KIRKWOOD for
    ! faster, less accurate simulations
    SSBP  KIRKWOOD ANGU HSR EMPI CAVITY   select type OH2 end

    ! keep peptide at centre; otherwise it may diffuse to the
    ! SSBP interface
    GEO RCM SPHERE XREF @xcen YREF @ycen ZREF @zcen -
        HARMONIC FORCE 5.0 -
        SELECT solute END
  END

  !------------------------------------------------------------
  ! brief minimization
  !------------------------------------------------------------
  ! should be on water only
  CONS HARM FORCE 50.0 SELECT solute END
  MINI SD NSTEP 100
  CONS HARM FORCE 0 SELECT all END       ! free solute


  !------------------------------------------------------------
  ! DIMS
  !------------------------------------------------------------
  DIMS DBNM DSCALE 8e-3 SKIP 1000 BSKIP 40 NBIAS 21  -  ! dims options
     SERL GENR SCAL 0.5882 TMEM 420 MEMO 20 MEMA 350 NMOD 50  -  ! NM options
     COFF 0.5 HARD ORIENT 100 COMB 3 NBES 15  SCAV 10         -
     MTRA 0 NMUN -1 NWIND 0  DSUN -1                          -
     SELECT mydims END                            -  ! DIMS selection
     SELECT mydims END                               ! BNM selection (optional)

  !------------------------------------------------------------
  ! dynamics
  !------------------------------------------------------------
  SCALAR fbeta SET 5.0 SELECT .NOT. TYPE H* END
  SHAKE BONH PARA

  OPEN WRITE FILE UNIT 52 NAME dims.dcd     ! trajectory

  set TEMPERATURE = 300
  set Nstep = 10000

  ! Run with Langevin dynamics
  ! (frequent output for testing)
  PRNLEV 4    ! verbose DIMS output
  DYNAMICS  start           nstep   @Nstep  timestp   0.002  iprfrq    100  -
            nprint   100    echeck  10000.0  -
            iasvel       1  -
            firstt @TEMPERATURE  finalt  @TEMPERATURE  tstruc  @TEMPERATURE  -
            langevin        tbath  @TEMPERATURE -
            inbfrq      10  imgfrq      -1  ihbfrq        0  ilbfrq       0  -
            nsavcrd     50  isvfrq       0  -
            iunread     -1  -
            iunwrite    -1  iuncrd      52 -
            omsc                                ! DIMS Onsager-Machlup score


.. _dims_tmd_last_mile_example:

NM-DIMS, using Targeted MD for the 'last mile'
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Close to the target the normal modes may not be specific enough to drive the
transition completely to the target. If this is desired, other methods can be
used to close the transition. The simplest approach is to use targeted MD once
the RMSD becomes smaller than COFF = 0.9 Angstrom.

::

   COOR COPY dims      ! for DIMS
   COOR COPY targ      ! import for TMD calculations

   TMDInitialize sele all end sele all end inrt 10 dincr 0.00025 ! Initialize TMD

   DIMS DBNM DSCALE 0.08 SKIP 500 BSKIP 50 NBIAS 27  - ! dims options
      SERL GENR SCAL 0.5882 TMEM 900 MEMO 30 MEMA 700 NMOD 50 - ! BNM options
      COFF 0.9 TMD                                   - ! NM TMD Cutoff
      ORIENT 400                                     - ! Orient every 400th
      COMB 3 NBES 15  MTRA 0 NMUN -1 NWIND 0 DSUN 11 - ! NM comb options
      SELE mydims END                                  ! selection

TMD pulls the structure linearly. It is also possible to employ the HARD
flavor of DIMS for a more natural transition. Note that CHARMM must be
compiled with TMD support (not enabled by default).


.. _dims_choosing_parameters:

Choosing parameters for DIMS-NM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Initial testing shows that a suitable choice of parameters can be rather
system dependent. The following is meant as a rough guideline to find
appropriate parameters.

* Use DBNM (the Block Normal Modes version) in conjunction with a
  last-mile option.  Typical DSCALE values are on the order of 1.0*10^-2
  to 1.0*10^-5.

* Optimize BSKIp, SKIP, and NBIAs so that a transition close to the target
  is achieved, say within COFF 0.5 Å (just using DBNM).

* Start with a reasonable number of dynamics steps, e.g. NSTEp 100,000 and
  a time step of 1-2 fs (SHAKE can be used).

* Look at the Onsager-Machlup score (should rise as slowly as possible)
  and the per-step DIMS score (should be as close to 1 as possible).

  * The total DIMS score will eventually decay to 0; that's a current
    limitation. A logarithmic DIMS score is displayed as well.
  * To get a feeling for the OM score, compute it for your system when
    running free Langevin MD (use the new OMSC time series). Note that
    the OM score depends on the step size (ie your trajectory SKIP
    value).
  * The aim is to produce a variety of trajectories with low OM
    scores. This may require running the dynamics for longer (larger
    NSTEp) to bias the system more gradually.

* Analyze the transition:

  * Plot RMSD (or more generally, the progress score) over steps.
  * Plot total potential energy over steps.
  * Plot the OM-score over time. (This is a cumulative measure and the
    last frame's score is the score for the trajectory. It is also
    available in the energy variable ?OMSCORE.)


.. _dims_dims_output:

Explanation of the DIMS output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

(1) Regular DIMS output

    ::

      DIMS> DIMS Score:   0. -28164.3312 -28159.1454  2.34808744  0.00559546636 3 2
                          |  |            |           |           |             | |
        DIMS_averaged-----+  |            |           |           |             | |
        Energy_unbiased------+            |           |           |             | |
        Energy_biased---------------------+           |           |             | |
        Normalization_constant------------------------+           |             | |
        current_move_score----------------------------------------+             | |
        dims_counter------------------------------------------------------------+ |
        number of steps over which normalization constant was computed------------+

      DIMS> DIMS LogScore:  0.28724E+04
                            |
        Log(DIMS)-----------+

    If the store goes to zero then the  best way to monitor the dims score is to
    save the scores to a file (using DSUN) and post-process them.

(2) At prnlev 3

    When running NM-DIMS the Progress Score is printed ; by default the progress
    score is the RMSD to the target (lower is better).

    ::

      DIMS> Progress Score0:   1.58859248
       DIMS> Move accepted:
         1.57390
              8.01332  36
              3.88053  30
             11.78028  46

      Progress Score0
          RMSD from target before DIMS: previous DIMS and MD

      Move accepted
          DIMS found a move along normal modes. After the move the new RMSD is given
          (here: 1.57390). The modes are listed below in the format

             frequency   mode_number

          Note that this triple (because COMB 3) of modes is the "best" combination
          of modes out of all NBES 15 modes in the list. (Note: all combinations of
          modes are checked and the best one is used, which may be a single mode or
          only a combination of two, even if COMB 3.)

(3) At prnlev 4

    * score of every combination

(4) At prnlev 6

    * lots of output - only useful for debugging


.. _dims_references:

References for DIMS
-------------------

1.  T. B. Woolf, Path corrected functionals of stochastic
    trajectories: Towards relative free energy and reaction
    coordinate calculations, Chemical Physics Letters 289(5-6)
    (1998) 433-441.

2.  D.M. Zuckerman and T. B. Woolf, Dynamic reaction paths and
    rates through importance-sampled stochastic dynamics, J Chem
    Phys 111 (1999) 9475-9484.

3.  D. M. Zuckerman, T.B. Woolf, Rapid Determination of Multiple Reaction
    Pathways in Molecular Systems: The Soft-Ratcheting Algorithm.
    arxiv:physics/0209098 (2002)

4.  H. Jang and T. B. Woolf, Multiple pathways in conformational
    transitions of the alanine dipeptide: An application of dynamic
    importance sampling, Journal of Computational Chemistry 27(11)
    (2006) 1136-1141.

5.  J. R. Perilla, A. Nagarajan, E. J. Denning, J. M. Johnston, O. Beckstein,
    T.B. Woolf, Sampling macromolecular transitions with dynamic
    importance sampling. (in preparation)

NM-DIMS uses the Block Normal Mode  routines in Charmm so
please also cite

6. Li G and Cui Q. A coarse-grained normal mode approach for
   macromolecules: an efficient implementation and application to
   Ca(2+)-ATPase. Biophys J 2002 Nov; 83(5) 2457-74.

Ref 6 describes in more detail how to choose blocks for the BNM approach:

7. Tama F, Gadea FX, Marques O, and Sanejouand YH. Building-block approach
   for determining low-frequency normal modes of macromolecules. Proteins
   2000 Oct 1; 41(1) 1-7

Targeted Molecular Dynamics (TMD option):

1. J. Schlitter, M. Engels, P. Krueger, E. Jacoby and A. Wollmer,
   Targeted Molecular Dynamics Simulation of Conformational
   Change-Application to the T --> R Transition in Insulin
   Mol. Sim. 10 (1993) 291-308.

Authors and Contact details:

Please direct feedback (questions, bug reports, suggestions) to

|   Tom Woolf <twoolf@jhmi.edu>
|   Juan Roberto Perilla <jrperillaj@jhu.edu>
|   Oliver Beckstein <orbeckst@jhmi.edu>


.. _dims_developer:

Remarks (for developers)
------------------------

Necessary flags for compiling with DIMS:

::

  VIBBLOCK
  DIMS

Optional flags:

::

  TMD        # targeted MD for the 'last mile'
  ARPK       # see ARPACK section below
  ANNLIB     # see ANNLIB section below. Required to activate IATD and ADIST.

External libraries:

::

  ANNLIB  -       A Library for Approximate Nearest Neighbor Searching.
                  http://www.cs.umd.edu/~mount/ANN/
  ARPK    -       ARPACK is designed to solve large scale eigenvalue problems.
                  http://www.caam.rice.edu/software/ARPACK/

Testcase files for DIMS can be found in the c35test directory:

::

  dims.inp
  dims_nm.inp

(Email Thomas B. Woolf <twoolf@jhmi.edu> or Juan Roberto Perilla
<jrperillaj@jhu.edu> with questions or feature requests.)



