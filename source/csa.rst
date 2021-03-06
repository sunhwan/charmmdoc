.. py:module:: csa

========================================================
Distributed CSA (Conformational Space Annealing) Command
========================================================

1. The distributed CSA commands will be based on the recent
   MSCALE commands in CHARMM.  The MSCALE commands allow diverse components
   of a single hamiltonian to be calculated on additional processors.

2. The CSA module will distribute a repetitive workload of many processors.
   It will consist of a series of commands, some on the the master and
   some on the slaves.  Slaves can be setup as CHARMM scripts, or run
   as other separate utility programs.

.. _csa_syntax:

Suggested Syntax
----------------

ON THE MASTER:

::

    MASTer -                          ! Use this processor as a master
           [ NSUBsystems integer ] -  ! How many slaves to generate (typ. 50)
           [ PROGram filename ]   -   ! What program to run the slaves with (typ. CHARMM)
           [ NPROC integer ] -        ! How many processors each slave will use (def. 1)
           atom-selection -           ! ?? Which atoms will participate (def. all)
           [ INPUt filename ] -       ! Input script for each slave
           [ OUTPut filename ]        ! Output file from each slave

    CSA  -          ! The Conformational Space Annealing command
           [  ] -   ! Details to be worked out and implemented
           [  ] -
           [  ]

    ETRAjectory   ! Compute the frame by frame energy of a trajectory file.

    SYSDisplay   !  Display the info about the whole setup

ON THE SLAVES

::

    RECEive -                 ! Receive data and coordinates from the Master
             [ IFDOne word ] - ! Where to go if exit condition is received
             [ repeat(flags) ]

    TRANsmit -                ! Send back energies to the Master
             [ COORdinates ] - ! Send back coordinates too
             [ FORCes ] -      ! Send back forces too
             [ DATA ]          ! ??  Send back other data


.. _csa_example:

Example
-------

Sample CHARMM Slave Script

::

    ...
    (setup PSF and stuff)
    ...
    LABEL LOOP
    RECEIVE  IFDONE NEXT  flag1 flag2
    IF @flag1 .eq. 1 THEN MINI SD 50  ...
    MINI ABNR 200 ...
    if @flag2 .eq. 1 THEN print energy unit 22
    TRANSMIT COORDINATES
    GOTO LOOP

    LABEL NEXT
    ...
    ...


.. _csa_notes:

Additional Notes
----------------

Slaves should send a semaphore when ready to receive.
Communication should be non-blocking in general.  The master should
not assume that all slaves will take the same amount of time.

Details for the high level FORTRAN CALLS for master-slave
communication need to be worked out (a job for Milan).

Suggestion for development:
Develop the ETRAJ command first.  This will have all of the basic
communication calls, and can be used as a base by the CSA developers.
The ETRAJ command may proves useful as a stand alone tool...

Check syntax above for possible conflicts with other developments

