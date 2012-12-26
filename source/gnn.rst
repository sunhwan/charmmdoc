.. py:module:: olap

======================
Genetic Neural Network
======================

A genetic neural network (GNN) method is provided for efficient determination
of quantitative structure property relationships. See the references given
below for a description of the GNN and its application. Some details specific
to the CHARMM implementation follow.

The GNN keyword must be included in pref.dat for the code to be compiled.

The input and output vectors of the data set are internally scaled to take
values between 0.1 and 0.9. The format of the data file is described in the
examples section.

Steepest descent back-propagation neural network is used to evaluate model
predictive quality. Jackknife cross-validation residual rms errors are
reported if no test data are specified. Only one hidden layer is employed.

Exhaustive enumeration and two genetic algorithm variants, genetic function
approximation (GFA) and evolutionary programming (EP), are available for
selecting models (sets of descriptors). The stochastic reminder method and
elitism are included for GFA reproduction.


.. _gnn_syntax:

Syntax required to invoke GNN
-----------------------------

::

    GNN [ data-spec ] [ nn-spec ] [ ga-spec ]

    data-spec ::= [ NDATa 1 ] [ NPROd 0 ] [ NPARa 1 ] [ UNIT -1 ] [ SEED 123 ]

    nn-spec ::= [ NDES 1 ] [ NHIDden 2 ] [ NTARg 1 ] [ NSWEep 100 ] [ MU 0.5 ] [ ETA 0.5 ]

    ga-spec ::= [ EXHAust ] [ GFA ] [ EP ] [ NPOPu 500 ] [ NGEN 200 ] [ FITNess 5.0 ]

.. _gnn_description:

Description of GNN specific keywords
------------------------------------

======== ========================================================================
NDATa    Number of data points in the training set.

NPROd    Number of data points in the test set.

NPARa    Number of candidate descriptors.

UNIT     Unit number from which data are imported.

SEED     Seed for random number generator.

NDES     Number of descriptors for the neural network.

NHIDden  Number of nodes in the hidden layer.

NTARg    Number of target parameters to predict.

NSWEep   Number of sweeps training.

MU       Momentum constant.

ETA      Learning rate.

EXHAust  Exhuastive enumeration.

GFA      Genetic function approximation.

EP       Evolutionary programming.

NPOPu    Number of individual models in reproduction pool.

NGEN     Number of generations to reproduce.

FITNess  Average fitness of models in the reproduction pool before terminating
         genetic algorithms. Fitness is defined as the reciprocal of the
         residual rms error.
======== ========================================================================

.. _gnn_examples:

Examples
--------

Data File (represented here symbolically)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    P1(1) P2(1) P3(1) P4(1) P5(1)
    P1(2) P2(2) P3(2) P4(2) P5(2)
    P1(3) P2(3) P3(3) P4(3) P5(3)
    P1(4) P2(4) P3(4) P4(4) P5(4)
    P1(5) P2(5) P3(5) P4(5) P5(5)

Note: Descriptors P1, P2 and P3, target parameters P4 and P5, training set
with data points (1), (2), and (3), and test set with data points (4) and (5).

Input
^^^^^

For example, file gnn.dat has 53 lines and 6 columns.

open read card unit 18 name gnn.dat

1. Exhaustive enumeration + Cross-validation + 1-Descriptor network

   ::

     gnn ndata 53 nprod 0 npara 5 unit 18 seed 123 -
         ndes 1 nhidden 2 ntarg 1 nsweep 100 mu 0.5 eta 0.5 -
         exhaust

2. GFA + Test set residual rms error evaluated + 2-Descriptor network

   ::

     gnn ndata 30 nprod 23 npara 4 unit 18 seed 123 -
         ndes 2 nhidden 3 ntarg 2 nsweep 100 mu 0.5 eta 0.5 -
         gfa npopu 2 ngen 10 fitness 5.0

   Note: ndata + nprod = 53 (number of lines),
   npara + ntarg = 6 (number of columns).

.. _gnn_references:

References

The GNN method was originally introduced in:

* Sung-Sau So and Martin Karplus, Evolutionary optimization in quantitative
  structure-activity relationship: An application of genetic neural networks,
  J. Med. Chem., 39:1521-1530 (1996).

* Sung-Sau So and Martin Karplus, Genetic Neural Networks for Quantitative
  Structure-Activity Relationships:  Improvements and application of
  benzodiazepine affinity for benzodiazepine/GABA_A receptors, J. Med. Chem.,
  39:5246-5256 (1996).

* Jie Hu and Aaron Dinner implemented the version in CHARMM.  It differs from
  the HIPPO program of So and Karplus primarily in that steepest descents rather
  than scaled conjugate gradients optimization is used to train the neural
  networks.  Performance is comparable for the data in:

  Jie Hu, Ao Ma, and Aaron R. Dinner, A two-step nucleotide flipping mechanism
  enable kinetic discrimination of DNA lesions by AGT, Proc. Natl. Acad. Sci.
  USA, in press (2008).

In addition to the above studies, papers using the CHARMM GNN method should
cite the introduction of the use of the GNN (and, more generally, informatic
methods) for determination of reaction coordinates:

* Ao Ma and Aaron R. Dinner, Automatic method for identifying reaction
  coordinates in complex systems, J. Phys. Chem. B, 109:6769-6779 (2005).

For a review of the GNN in other contexts, see:

* Aaron R. Dinner, Sung-Sau So, and Martin Karplus, Statistical analysis of
  protein folding kinetics, Adv. Chem. Phys., 120:1-34 (2002).

For a more detailed discussion of the neural network component used in the
CHARMM implementation, see the book:

* Jure Zupan and Johann Gasteiger, Neural Networks for Chemists:  An
  Introduction, VCH, New York (1993).

