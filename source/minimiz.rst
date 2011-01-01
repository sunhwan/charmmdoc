.. py:module:: minimiz

===============================================
Energy Manipulations: Minimization and Dynamics
===============================================

One can minimize the energy by adjusting the coordinates
of all the atoms in order to reduce its value. Several minimization
algorithms are provided. They include:

*        Steepest Descent (SD)
*        Conjugate Gradient (CONJ)
*        Adopted Basis Newton-Raphson (ABNR)
*        Newton-Raphson (NRAP)
*        Powell (POWE)
*        Truncated Newton Method (TNPACK)

.. _minimiz_syntax:

Syntax for Energy Manipulation Commands
---------------------------------------

::

   MINI     { SD     steepd-spec  } [ nonbond-spec ] [ hbond-spec ] -
            { CONJ   conj-spec    } [   INBFrq 0   ] [  IHBFrq 0  ] [NOUPdate]
            { ABNR   abnr-spec    }
            { NRAP   nrap-spec    }
            { POWEll powell_spec  }
            { TN     tnpack-spec  }

                   [STEP real] [GRADient] [NUMErical]
                      [ frequency-spec ] [ tolerance-spec ] [ io-spec ] }

               [ CHEQ [CGMD int] [CGIN] [CGFC] [PBEQ] [QPOL [IPOL int] ] ]

   hbond-spec::=     *Note Hbonds:(chmdoc/hbonds.doc).
   nonbond-spec::=   *Note Nbonds:(chmdoc/nbonds.doc).

   frequency-spec::= [NSTEP int] [IHBFrq int] [INBFrq int] [NPRInt int]

   tolerance-spec::= [TOLENR real] [TOLGRD real] [TOLITR int] [TOLSTP real]

   io-spec::= [DEBUG] [IUNCrd int [NSAVC int] ] [IUNXyz [NSAVX int] [MXYZ int] ]

   conj-spec::= [NCGCyc int] [PCUT real] [PRTMin int]
                     [LATTice] [NOCOordinates]

   powell-spec::=    [LATTice] [NOCOordinates]

   steepd-spce::=  [NOENergy ] [LATTice] [NOCOordinates]

   abnr-spec::= [EIGRng real] [MINDim int] [STPLim real] -
                   [STRIct real] [ MASS ] [PSTRct real]
                     [LATTice] [NOCOordinates] [FMEM real]

   nrap-spec::=  [TFREq real]

   tnpack-spec::= [NCGCyc int] 
                   [PREC or NOPR] [USER or OURH]  [REST or QUAT]  [NOSC or SCHE]
                    [DEFS or SEAR] [LOWP or HIGP]  [IORD or NOOR] [PERM or NOPM]  


.. _minimiz_description:

Options common to minimization and dynamics
-------------------------------------------

The following table describes the keywords which apply to all
minimization methods.

+-------------+-------+----------------------------------------------------------+
|Keyword      |Default|Purpose                                                   |
+-------------+-------+----------------------------------------------------------+
|NSTEP        |100    |The number of steps to be taken. This                     |
|             |       |is the number of cycles of minimization, not the number   |
|             |       |of energy evaluations.                                    |
+-------------+-------+----------------------------------------------------------+
|INBFRQ       |50     |The frequency of regenerating the non-bonded list.        |
|             |       |The list is regenerated if the current step number        |
|             |       |modulo INBFRQ is zero and if INBFRQ is non-zero.          |
|             |       |Specifying zero prevents the non-bonded list from being   |
|             |       |regenerated at all.                                       |
+-------------+-------+----------------------------------------------------------+
|IHBFRQ       |50     |The frequency of regenerating the hydrogen bond list.     |
|             |       |Analogous to INBFRQ                                       |
+-------------+-------+----------------------------------------------------------+
|non-bond-spec|       |The specifications for generating the non-bonded list.    |
|             |       |See doc:`nbonds`.                                         |
+-------------+-------+----------------------------------------------------------+
|hbond-spec   |       |The specifications for generating the hydrogen bond list. |
|             |       |See doc:`hbonds`                                          |
+-------------+-------+----------------------------------------------------------+
|NPRINT       | 1     |The frequency with which energies are printed during      |
|             |       |the course of dynamics or minimization.                   |
+-------------+-------+----------------------------------------------------------+
|GRADient     |       |Minimize the magnitude of the gradient of the energy      |
|             |       |instead of the energy.                                    |
+-------------+-------+----------------------------------------------------------+
|NUMErical    |       |Forces will be determined by finite differences           |
+-------------+-------+----------------------------------------------------------+
|IUNCrd       |-1     |Unit to write out a "trajectory" file for the minimization|
+-------------+-------+----------------------------------------------------------+
|NSAVC        | 1     |Frequency for writing out frames (only with IUNCrd)       |
+-------------+-------+----------------------------------------------------------+
|DEBUg        |       |Extra print for debug purposes                            |
+-------------+-------+----------------------------------------------------------+

In the table which follows, keywords enclosed in square brackets
means that one can choose one in the set. Such enclosed keywords do not
expect a value after them. All other keywords are used for specifying
values, see :ref:`minimiz_syntax`. The method column shows which method the
keyword affects.

+--------+--------+----------+------------------------------------------------------+
|Keyword |Default |Method    |Purpose                                               |
+--------+--------+----------+------------------------------------------------------+
|[ CONJ ]| CONJ   |          | Do conjugate gradient minimization.                  |
+--------+        |          +------------------------------------------------------+
|[ SD   ]|        |          | Do steepest descent minimization.                    |
+--------+        |          +------------------------------------------------------+
|[ NRAP ]|        |          | Do Newton-Raphson minimization.                      |
+--------+        |          +------------------------------------------------------+
|[ ABNR ]|        |          | Do Adopted Basis Newton-Raphson minimization,        |
+--------+        |          +------------------------------------------------------+
|[ TN   ]|        |          | Do Truncated-Newton minimization.                    |
+--------+--------+----------+------------------------------------------------------+
|[MASS]  |        |          | with mass weighted forces if specified.              |
+--------+--------+----------+------------------------------------------------------+
|STEP    | .02    | ALL      | Initial step size for the minimization algorithms.   |
|        |        | except TN| Reasonable values for the various methods are best   |
|        |        |          | determined by trial and error.                       |
+--------+--------+----------+------------------------------------------------------+
|LATTice |        |   ABNR   | With the CRYSTAL facility, also optimize the unit    |
|        |        |          | cell box size and/or shape.                          |
+--------+--------+----------+------------------------------------------------------+
|NOCOords|        |   ABNR   | With the CRYSTAL facility, only optimize the unit    |
|        |        |          | cell. This leaves coordinates unchanged, but         |
|        |        |          | allows the box size and/or shape to change.          |
+--------+--------+----------+------------------------------------------------------+
|PRTMIN  |   1    |   CONJ   | A flag indicating how much to print during           |
|        |        |          | minimization.                                        |
|        |        |          | If less than 2, the energy is printed only once      |
|        |        |          | each cycle. A setting of 2 shows the energy for      |
|        |        |          | each evaluation plus variables used in the method.   |
+--------+--------+----------+------------------------------------------------------+
|NCGCYC  |  100   |    CONJ  | The number of conjugate gradient cycles executed     |
|        |        |          | before the algorithm restarts.                       |
+--------+--------+----------+------------------------------------------------------+
|PCUT    | .9999  |    CONJ  | If the cosine of the angle between the old and new   |
|        |        |          | P vector is greater than PCUT, the algorithm will be |
|        |        |          | restarted. This prevents the algorithm from plodding |
|        |        |          | down the same path repeatedly. If PRTMIN is less     |
|        |        |          | than 2, one effect of the restart is that the step   |
|        |        |          | size will go to its initial value. If this happens   |
|        |        |          | many times, you've converged.                        |
+--------+--------+----------+------------------------------------------------------+
|EIGRNG  | .0005  |   ABNR   | The smallest eigenvalue (relative to the largest)    |
|        |        |          | that will be considered nonsingular.                 |
+--------+--------+----------+------------------------------------------------------+
|MINDIM  |   5    |   ABNR   | The dimension of the basis set stored.               |
+--------+--------+----------+------------------------------------------------------+
|STPLIM  |   1.0  |   ABNR   | The maximum Newton Raphson step that will            |
|        |        |          | be allowed.                                          |
+--------+--------+----------+------------------------------------------------------+
|STRICT  |   0.1  |   ABNR   | The strictness of descent.  The energy of a new step |
|        |        |          | must be less than the previous best energy + STRICT  |
|        |        |          | for the new step to be accepted.                     |
+--------+--------+----------+------------------------------------------------------+
|MASS    | false  |   ABNR   | Use unweighted forces by default or mass weighted    |
|        |        |          | if specified.  Mass weights converge more slowly but |
|        |        |          | allow association with normal mode frequencies.      |
+--------+--------+----------+------------------------------------------------------+
|TFREQ   |   1.0  |   NRAP   | The smallest eigenvalue that is considered to be     |
|        |        |          | non-negative (i.e. do cubic fitting on all           |
|        |        |          | eigenvalues smaller than this).                      |
+--------+--------+----------+------------------------------------------------------+
|TOLENR  |   0.0  |   ABNR   | A tolerance applied to the change in total energy    |
|        |        |          | change during a cycle of minimization (NCYCLE steps).|
|        |        |          | If the energy change is less than or equal to        |
|        |        |          | TOLENR, the minimization routine will exit.          |
+--------+--------+----------+------------------------------------------------------+
|TOLGRD  |   0.0  |   ABNR   | A tolerance applied to the average gradient during   |
|        |        |          | a cycle of minimization.  If the average gradient    |
|        |        |          | is less than or equal to TOLGRD, the routine         |
|        |        |          | will exit.                                           |
|        +--------+----------+------------------------------------------------------+
|        |   1.0  |    TN    | A parameter which determines the desired accuracy    |
|        |        |          | of the computed solution. The following four         |
|        |        |          | convergence tests are checked:                       |
|        |        |          | T1) f(x_{k-1})-f(x_k) < tolgrd (1+|f(x_k)|)          |
|        |        |          | T2) ||x_{k-1} - x_k|| < sqrt(tolgrd) (1+||x_k||)     |
|        |        |          | T3) ||g(x_k)|| < tolgrd^(1/3) (1+ ||f(x_k)||)/100    |
|        |        |          | T4) ||g(x_k)|| < eg (1+ ||f(x_k)||)                  |
|        |        |          |                                                      |
|        |        |          | If TOLGRD is equal to 0. in the input file, TOLGRD   |
|        |        |          | set to 10^(-8) in the calculation. If it is equal    |
|        |        |          | to 1., it is set to 10^(-12).                        |
|        |        |          | eg is the square root of machine precision.          |
|        |        |          |                                                      |
|        |        |          | The routine will exit when either (T1,T2, and T3)    |
|        |        |          | are satisfied or (T4). (T4) is a useful test at      |
|        |        |          | the first Newton iteration or for comparison with    |
|        |        |          | other methods (see TNPACK paper).                    |
+--------+--------+----------+------------------------------------------------------+
|TOLITR  |  100   |   ABNR   | The maximum number of energy evaluations allowed     |
|        |        |   CONJ   | for a single step of minimization.                   |
+--------+--------+----------+------------------------------------------------------+
|TOLSTP  |  0.0   |   ABNR   | A tolerance applied to the average step size during  |
|        |        |          | a cycle of minimization.  If the average step size   |
|        |        |          | is less than or equal to TOLSTP, the routine         |
|        |        |          | will exit.                                           |
+--------+--------+----------+------------------------------------------------------+
|FMEM    |  0.0   |    ABNR  | Memory factor. It is used to compute average         |
|        |        |          | gradient and step size according to the formula :    |
|        |        |          |                                                      |
|        |        |          | AVERAGE_VALUE = FMEM * AVERAGE_VALUE                 |
|        |        |          | + (1-FMEM) * CURRENT_VALUE.                          |
|        |        |          |                                                      |
|        |        |          | FMEM=0 means no memory (i.e current value is used)   |
|        |        |          | and FMEM=1 means infinitely long memory (i.e.        |
|        |        |          | initial value will be used all the time).            |
+--------+--------+----------+------------------------------------------------------+
|NOUP    | false  |   ALL    | Turns off the list updates.                          |
+--------+--------+----------+------------------------------------------------------+
|PREC or | NOPR   |   TN     | selects preconditioning (PREC) or no preconditioning |
|NOPR    |        |          | (NOPR).                                              |
+--------+--------+----------+------------------------------------------------------+
|ANAL or | ANAL   |   TN     | selects option for Hd multiplication:                |
|FDIF    |        |          | ANAL for analytic version,                           |
|        |        |          | FDIF for the finite-difference formula.              |
+--------+--------+----------+------------------------------------------------------+
|REST or | REST   |  TN      | specifies choice of PCG truncation test:             |
|QUAT    |        |          | residual (REST) or quadratic (QUAT).                 |
|        |        |          |                                                      |
+--------+--------+----------+------------------------------------------------------+
|SCOF or | SCOF   |   TN     | specifies whether the scheduling subroutine is used  |
|SCON    |        |          | (SCON for on, SCOF for off). The subroutine turns    |
|        |        |          | on preconditioning (if chosen) when the gradient is  |
|        |        |          | smaller than some tolerance, and uses steepest       |
|        |        |          | descent steps beforehand.                            |
+--------+--------+----------+------------------------------------------------------+
|SRON or | SROF   |   TN     | specifies whether the optimal search-vector subrou-  |
|SROF    |        |          | tine is turned on (SRON) or off (SROF). This subrou- |
|        |        |          | tine considers more than one possible descent        |
|        |        |          | directions at a Newton iteration and chooses the     |
|        |        |          | one that leads to greatest energy reduction.         |
|        |        |          | Additional energy + gradient evaluations are         |
|        |        |          | required.                                            |
+--------+--------+----------+------------------------------------------------------+
|IORD or | NOOR   |   TN     | specifies whether a reordering of M will be performed|
|NOOR    |        |          | to minimize fill-in (IORD) or not (NOOR). This might |
|        |        |          | be useful if M is very large and sparse. The         |
|        |        |          | reordering is done only once, but the savings are    |
|        |        |          | reflected in each-inner loop iteration where a       |
|        |        |          | linear system involving M is solved.                 |
+--------+--------+----------+------------------------------------------------------+
|PERM or | NOPM   |   TN     | determines if the permutation array for reordering   |
|NOPM    |        |          | M is known when the current TNPACK call is made      |
|        |        |          | (PERM - known, NOPM - unknown).                      |
|        |        |          |                                                      |
+--------+--------+----------+------------------------------------------------------+
|NOEN    | FALSE  |   SD     | only use the information of force to minimize        |
|        |        |          | a system. implemented for the case of minimizing     |
|        |        |          | a reaction path using the eudged elastic band method.|
+--------+--------+----------+------------------------------------------------------+
|NSADD   | 0      |   NRAP   | sets the order of saddle point you want to find.     |
|        |        |          | NSADD=1 will search in the opposite direction of     |
|        |        |          | the most negative eigenvector (i.e. uphill) until    |
|        |        |          | a stationary point is located (i.e. transition state |
|        |        |          | at NSADD=1).                                         |
+--------+--------+----------+------------------------------------------------------+

Note that the following commands are equivalent:

* ANAL = USER
* FDIF = OURH
* SCOF = NOSC
* SCON = SCHE
* SRON = DEFS
* SROF = SEAR


.. _minimize_discussion:

Discussion of the various minimization methods
----------------------------------------------

There are several different algorithms for minimizing the energy
of the system. They all involve calculating the derivative of the
potential energy, and possibly the second derivative, and using that
information to adjust the coordinates in order to find a lower energy.
In the descriptions of the algorithms below, a capitalized keyword at
the beginning of each paragraph is the key word used to invoke that
method. After the descriptions is a listing of all keywords associated
with minimization.

The simplest minimization algorithm is steepest descent (SD).
In each step of this iterative procedure, we adjust the coordinates in
the negative direction of the gradient. It has one adjustable parameter,
the step size, which determines how far to shift the coordinates at each
step. The step size is adjusted depending on whether a step results in a
lower energy. I.e., if the energy drops, we increase the step size by
20% to accelerate the convergence. If the energy rises, we overshot a
minimum, so the step size is halved. Steepest descent does not converge
in general, but it will rapidly improve a very poor conformation.

A second method is the conjugate gradient technique (CONJ) which has
better convergence characteristics (R. Fletcher & C.M. Reeves, The Computer
Journal 7:149 (1964)). The method is iterative and makes use of the previous
history of minimization steps as well as the current gradient to determine the
next step. It can be shown that the method converges to the minimum energy in
N steps for a quadratic energy surface where N is the number of degrees of
freedom in the energy. Since several terms in the potential are quadratic, it
requires less evaluations of the energy and gradient to achieve the same
reduction in energy in comparison to steepest descent. Its main drawback is
that with very poor conformations, it is more likely to generate numerical
overflows than steepest descent. The algorithm used in CHARMM has a slightly
better interpolation scheme and automatic step size selection.

A third method is the conjugate gradient powell method minimizer
(POWE) (A. Brunger). Its efficiency is much improved over the Fletcher
and Reeves method. The use of the POWELL minimizer is encouraged whenever
ABNR is not possible. The POWELL minimizer has no INBFRQ or IHBFRQ
feature. The use of CHARMM loops to mimic this important feature is
encouraged. The CHARMM loop facilities allow harmonic constrained
minimizations with periodic updates. In case of bad contacts or
unlikely conformations the SHAKE method might give an error when the
displacements become to large. Using a harmonic constraint setup
with periodic updates solves this problem.

A fourth method involves solving the Newton-Raphson minimization
equations iteratively (NRAP). This procedure requires the computation of
the derivative of the gradient which is a matrix of size O(n**2). The
procedure here is to find a point where the gradient will be zero
(hopefully a minimum in energy) assuming that the potential is
quadratic. The Newton-Raphson equations can be solved by a number of
means, but the method adopted for CHARMM involves diagonalizing the
second derivative matrix and then finding the optimum step size along
each eigenvector. When there are one or more negative eigenvalues, a
blind application of the equations will find a saddle point in the
potential. To overcome this problem, a single additional energy and
gradient determination is performed along the eigenvector displacement
for each small or negative eigenvalue. From this additional data, the
energy function is approximated by a cubic potential and the step size
that minimizes this function is adopted. The advantages of this method
are that the geometry cannot remain at a saddle point, as sometimes
occurs with the previous procedures, and that the procedure converges
rapidly when the potential is nearly quadratic (or cubic). The major
disadvantage is that this procedure can require excessive storage
requirements, O(n**2), and computation time, O(n**3), for large
molecules. Thus we are currently restricted to systems with about 200
atoms or less for this method. IMAGES and SHAKE are currently unavailable
with this method.

The fifth method available is an adopted basis Newton-Raphson
method (ABNR) (D. J. States). This routine performs energy minimization
using a Newton-Raphson algorithm applied to a subspace of the coordinate
vector spanned by the displacement coordinates of the last (mindim)
positions. The second derivative matrix is constructed numerically from
the change in the gradient vectors, and is inverted by an eigenvector
analysis allowing the routine to recognize and avoid saddle points in
the energy surface. At each step the residual gradient vector is
calculated and used to add a steepest descent step onto the
Newton-Raphson step, incorporating new direction into the basis set.
This method is the best for most circumstances.
SHAKE is currently unavailable with this method.
 	
The sixth method is the truncated-Newton (TN) minimization 
package TNPACK, developed by T. Schlick and A. Fogelson.  TNPACK is
based on the preconditioned linear conjugate-gradient technique for
solving the Newton equations.  The structure of the problem ---
sparsity of the Hessian --- is exploited for preconditioning.
Thorough experience with the new version of TNPACK in CHARMM has been
described in the paper: Journal of Computational Chemistry: 15,
532--552, 1994.  Through comparisons among the minimization algorithms
available in CHARMM, we find that TNPACK compares favorably to ABNR in
terms of CPU time when curvature information is calculated by a
finite-difference of gradients (the "numeric" option of TNPACK).
With the analytic option, TNPACK can converge more rapidly than ABNR
for small and medium systems (up to 400 atoms) as well as large
molecules that have reasonably good starting conformations; for large
systems that are poorly relaxed (i.e., the initial Brookhaven Protein
Data Bank structures are poor approximations to the minimum), TNPACK
performs similarly to ABNR.  SHAKE is currently unavailable with this
method.


Barriers and Minima
^^^^^^^^^^^^^^^^^^^

The GRADient option causes the minimizers to find a zero of the
target function (grad(V))^2.  The square of the gradient replaces the
energy in the minimizers.  Depending on the initial condition (initial
set of coordinates), the search can either be terminated in a minimum
or in a saddle point of the potential energy function (a barrier). If
the second derivative of the initial condition is negative BARI will
look for a saddle point; if it is positive it will stop at a minimum.
The second derivative matrix is employed to calculate first derivatives
of the target function. As a result it is much slower compared to
ABNR and NRAP in reaching a minimum.  For minimum energy calculations:
DO NOT USE THE GRADient OPTION.

The NSADD keyword turns on special code that follows positive
eigenvectors thus searching for a saddle point. Care must be taken
when choosing the starting structure for this code (i.e. you should
not start the search from a true minima as the code can get confused
about which eigenvector to follow). The best suggestion is to
slightly perturb your structure in the direction you believe that
transition state (or higher order saddle point) of interest to be. 
