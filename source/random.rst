.. py:module:: random

============================================
Random Number Generator Controlling Commands
============================================

The commands described in this section are for a control of
the random number generators in CHARMM.

.. _random_syntax:

Syntax of RANDom commands
-------------------------

Current Syntax:

::

    RANDom specifications:


    RANDom specifications:

    RANDom  { [CLCG] { [TIME]          } } [UNIForm]        [ASIN] [SCALe real] [OFFSet real] [TEST]
            {        { ISEEd 4X(int    } } [GAUSsian real ] [ACOS]
            {                            }
            { SYSTem { [TIME]          } }
            {        { ISEEd seed-specs} }
            {                            }
            { USER   { [TIME]          } }
            {        { ISEEd seed-specs} }
            {                            }
            { OLDRandom [ISEEd int]      }

    seed-specs::= repeat integer ?NSEED number of times. For the default
                  CLCG RNG 4 integer numbers need to be specified. For the
                  standard fortran randum_number() routine one need to
                  specify ?NSEED integer numbers. This number is compiler
                  dependent!

Integer random number generator IRANdom specifications:

::

    IRANdom                 [SERIes int]  [SETUp]  [BEGInt int]  [ENDInt int]
                            [SEED int]


.. _random_function:

1) RANDom command.

   The expression ?RAND will have a random number substituted for it
   during command line evaluation.  The default is to provide a number
   from a uniform distribution, between 0.0 and 1.0; the RANDom command
   allows modification of the distribution type and specification of
   other factors.  The only required keyword is the distribution type,
   which must be second; for a GAUSsian distribution, a value for sigma
   is required; the default mean is 0.0.


   There is a variety of random number generators (RNG) available in
   CHARMM. They are:

   ========= ===========================================================
   OLDRandom legacy RNG. It is not appropriate to use for production
             simulations anymore, but it can be used for simple testing
             and for comparison with older results and test cases. "Use
             The Source Luke" to find out the CHARMM command line
             parameters to switch back to either old random or old
             CLCG. It is very convenient to put these flags into
             test.com script and compare the test results with older
             CHARMM versions. The alternative would be to put random
             command in test/datadir.def file and in some 40 or so test
             cases input scripts which don't stream datadir.def.

   CLCG      new RNG that supports 100 series and uses 4 seeds. It is a
             modern RNG and is the default choice.

   SYSTem    whatever is provided by random_number() routine in
             fortran. It supports only one series but uses different
             numbers of seeds, depending on compilers and integer(8)
             vs integer(4) compilation. Use ?NRAND value to query

   USER      a user_random() function is provided in
             charmm/usersb.src. This routine can be replaced by users
             to test their own RNG.
   ========= ===========================================================

   Additional keywords:

   ================= ===========================================================
   UNIForm           uniform distribution - default

   GAUSSian sigma    Gaussian distribution. Value of sigma must be specified

   SCALe  scale      multiply the number by scale

   OFFSet offset     add offset to the number

   ACOS              treat the number as a cosine and return the angle (deg)

   ASIN              treat the number as a sine and return the angle (deg)

   ISEEd iseed       specify a new random seed(s) (integer(s)). Use
                     ?NSEED parameter to query how many seeds are
                     needed for the current random number generator.
                     NOTE: No big numbers for CLCG (iseed < 2 giga).
                     No limit for SYSTem random generator.

   TIME              seeds are assigned from the current system time.
                     The default.

   TEST              this command will test the random number generator for
                     its poriodicity.

   PARAllel          stores the seeds from all the processors in
                     the restart file, so parallel run can be
                     restored if needed. Without the parallel
                     keyword every processors has its own random
                     number series initialized from system time.
   ================= ===========================================================

   Note that OLDRandom sub-command sets OLDRNG, which runs "old" random
   number generator instead of "new" CLCG method. CLCG unsets OLDRNG, and
   runs the CLCG random number generator.

   About the seeds:

   The use of fixed seeds is discouraged for production runs so
   by default the system clock provides an initial seeds for the
   RNGs. One can specify seeds on the RANDom command or in some other
   commands, eg DYNA. The seed number stored in the trajectory file is
   for legacy RNG (OLDRandom) and its use is deprecated. However the full
   functionality of seed storage is implemented for the restart file. It
   is implemented so that the restart files from older versions of CHARMM
   can still be used, only the seed is ignored when running from the old
   restart file. The new restart file saves all the necessary seeds.

   Examples:

   ::

      RANDOM GAUSS 0.2 SCALE 10.0   !     gaussian  mean of 0.0 with a sigma of 2.
      RANDOM UNIFORM SCALE 360.     !     uniform   0. to 360
      RANDOM UNIFORM ACOS SCALE .5  !     uniform   angles with cosines from 0. to .5
      RAND GAUS 5. OFFS 60.         !     gaussian  mean of 60. with a sigma of 5.
      RAND UNIF ISEED 7734          !     uniform   new random seed

   Subsequent use of ?RAND will substitute a number from the appropriate
   distribution.


2) IRANdom command.

   This command is designed to generate series of random integers taken
   from uniform distributions between user-specified limits.  Each series
   or distribution must first be set up with the ``IRANdom SETUp`` command,
   in which the lower and upper limits of the distribution, the series
   number, and an integer seed are specified. E.g.

   ::

       IRAND SERIES 1 SETUp BEGI 1 ENDI 18  SEED 2346
       IRAND SERIES 2 SETUp BEGI 1 ENDI 402 SEED 4028987

   The random integers for each series are then generated with the commands

   ::

       IRAND SERIES 1
       IRAND SERIES 2
       etc.

   The ?iran expression accesses the last random integer generated.

   The purpose of the multiple series feature is at least two-fold. First, it
   allows users to generate random numbers easily from many different
   distributions during the same charmm run (e.g. for use in different parts
   of the same calculation). Second, it may help the user avoid correlations
   between random numbers generated for different parts of a calculation.
   An internal counter, corresponding initially to the seed, is incremented
   by several units with each instance of the IRANdom command; by separating
   the seeds of the various distributions sufficiently, the user can thus avoid
   cross-series correlations.  The use of multiple seeds for a given series
   should be unnecessary and is discouraged.  The IRANdom function has an
   overall period of no less than 10^12 for distribution widths of 10^10 or less.
   IRANdom can also be used to effectively generate random real numbers, through
   a division of the generated integers by a constant, with the use of the
   CALC command.
