.. py:module:: gpu

====================================
Installing and running CHARMM on GPU
====================================

Using GPUs for general computing tasks is becoming
increasingly popular. CHARMM is capable to use the GPU for
calculations of non-bonded part of the energy and forces. The
code in CHARMM is based on the GRAPE/CHARMM interface mostly
implemented in the nbonds/grape.src file.

.. _gpu_syntax:

::

  [SYNTAX NBONDs][SYNTAX ENERgy][SYNTAX DYNAmics]

    [NBONds]   GRAPe <int> restricted-nonbond-spec
    [ENERgy]
    [DYNAmics]

  restricted-nonbond-spec::= ... shift vshift ...

Any nonbond keyword and value may be specified, except that only shift
cuttof method is supported on the GPU.

As of November 2011 only shift cutoff method is currently
supported. But support for others is coming soon. The rest of
nonbonded specifications are as in standard CHARMM.

When grape <int> is specified all the non-bond calculations are performed
on the GPU so there is no need for the non-bond list updates on the host
CPU. This means that the inbfrq and imgfrq can be set to either 0 or a
large number. For an example see: test/cbenchtest/mbco/gputest.inp file.

Useful grape <int> flag values represent the following:

========= =================================================================
grape 11  the fastes GPU calculations. It is limited only for single
          CPU/GPU runs. CHARMM may be compiled with parallel, so it
          is still usable in parallel/parallel setups.

grape 13  calculates pressure, and it also supports parallel
          calculations: both multiple boards in one box as multiple
          boxes with arbitrary number of boards.
          There is an extra performance penalty for this
          calculation, ie it is slower than grape 11. In this
          case the image nonbond list must be updated like in the
          host calculations, which slows the calculation even further.
          Currently works with BYCB only!
========= =================================================================

The rest are more for debugging purposes

========= =================================================================
grape 0   separate calculation of VDW and electrostatic forces, so it is
          slower than 1. PME is calculated on the host.

grape 1   PME is calculated on the host (usefull for some Macs, which
          lack the hardware to perform PME calculations on the GPU).

grwpe 2   everything is calculated on the host, but with the mr3
          library, so it is extremly slow. Currently not working properly.

grape 3   pressure (general image) calculation with PME calculated on the host.

grape 10   separate calculation of VDW and electrostatic forces, so it
           is slower than 11

grwpe 12   Real space is calculated on the host, but with the mr3
           library, so it is extremly slow. Currently not working properly.
========= =================================================================

Also note that the printout of individual energy terms differ when
using different grape flags. Notably when using methods other than 3
or 13 there is no separation for image and main atom energies.

The gpu code is using cell index method for cutoff calculations and
one can specify the number of cells in the input scripts by the
following environment commands. By default these are now calculated for
both grape 11 or grape 13. But if specified in the input script the
values from the script are taken for the setup of the cell structure.

::

  envi VG_LDIMX 7
  envi VG_LDIMY 7
  envi VG_LDIMZ 7

Higher number gives better performance, but there is a limit
which can be estimated from a trial calculations. Once the number is
too high the results are either wrong or in some cases the execution
of the program may stop.

For "grape 11" the code supports only orthorombic symmetry. In the
case the system is non-cubic one must specify the cell dimensions in
the input script by the following commands:

::

  envi vg_cellsizex 90.0
  envi vg_cellsizey 90.0
  envi vg_cellsizez 84.0

For "grape 13" there is no limit in the crystal symmetry. Whatever
CHARMM can build gpu will calculate.

There is also a memory limitation on PME calculations with cell index
method. The code automatically switches to no cell index when the
FFTX*FFTY*FFTZ > 128*128*128. But this is also the system size
dependent so for larger systems at 128*128*128 cell index might not
work. In this case please use

::

  envi VG_NOCELLINDEX_WAVE 1

in your input script.

There is also a limit for FFTX*FFTY*FFTZ > 384*384*384 even when using
no cell index method.


.. _gpu_installation:

Installation
============

Installing GPU version of CHARMM is performed by the following command:

::

  install.com gpu {charmm-size} [ifort] [M] [mpif90]

Before this command is issued several libraries have to be installed
in the system. Currently it is tested for 64-bit linux systems
only. The following is the list of the files and libraries that are
needed to compile and run GPU CHARMM.

1. Install CUDA library on the system:

   - from NVIDIA web site this is:

     * CUDA driver  (normal X11 graphics driver, needed to run CHARMM)
     * CUDA Toolkit (library needed to compile and run CHARMM)

     As of August 2011 the following packages are available from
     NVIDIA website:

     * NVIDIA-Linux-x86_64-260.19.36
     * cudatoolkit_3.2.16_linux_64

     or for CUDA version 4

     * NVIDIA-Linux-x86_64-275.21.run
     * cudatoolkit_4.0.17_linux_64

     At this point it is maybe worth trying some of the demo programs
     from NVIDA SDK package. At least matrixMul may tell you if your
     hardware is working properly.

2. set additional environment variables to compile CUDA code with CHARMM:

   For example if the CUDA libraries from NVIDIA were installed in
   /opt/cuda directory then you need to do the following:

   ::

      export PATH=/opt/cuda/bin:$PATH
      export LD_LIBRARY_PATH=/opt/cuda/lib64:$LD_LIBRARY_PATH

   The following is guessed from the path above but in the case of any
   problem you can set it, too:

   ::

      export CUDATK=/opt/cuda

   or for csh like shells:

   ::

      set path = ( /opt/cuda $path )
      setenv LD_LIBRARY_PATH /opt/cuda/lib64:$LD_LIBRARY_PATH
      setenv CUDATK /opt/cuda

3. Then run the above install.com command. If you need to recompile
   everyhting then also the tool/gpu/rs and tool/gpu/pme need to be
   cleaned up: make clean in both drectories will do it.

4. For parallel compile add M mpif90 flags on install.com command as
   described above. Also GRAPE 11 in the dyna command needs to be
   changed to GRAPE 13 IMALl. To use different GPU card for each
   process in the same box these commands can be specified in the
   input script:

   ::

       envi vg_deviceid 0
       if ?mynode .eq. 1 then envi vg_deviceid 1

   The following is now just for the debugging purposes. The table files
   r1.g80emu and rsqrt.g80emu are now produced during the initialization
   process at the start of the run. The tables are not generated by
   default anymore.

5. When running charmm, don't forget to link to the binary tables
   tool/gpu/r1.g80emu, and tool/gpu/rsqrt.g80emu to a current
   directory! One can also use VG_EMUDIR environment variable to
   specify the location of the *.g80emu file. This can be done also
   inside the CHARMM input script with envi vg_emudir "directory-name"
   These tables may be generated by running ./emutestvg w in the
   tool/gpu/rs directory.

.. _gpu_parallel:

Parallel
========

When running GPU compiled CHARMM in parallel it works the same
as a standard parallel CHARMM. The calculation takes the same number
of GPUs as is the number of CPUs. However in some setups it is
desirable to split the calculations in a variety of parts. To setup
this one need to specify some parameters with the PARAllel command in
the CHARMM input script. The syntax is:

::

  PARAllel GPUSplit  [ NGR <int> ] [ NGK <int> ]
                       [ NCH <int> ] [ NCE <int> ] ..

work in progress....

::

    NGR   Number of Gpu for Real space
             how many GPUs for the RS calculaiton
    NGK   Number of Gpus for K-space calculation (PME)
             how many GPUs for the KSPACE calculaiton.
             Obviously only 1 or 0 make sense here
             In case of 0 host will do the KSPACE
    NCH   Number of Cpus for Host calculations
             How many CPUs are used for the non GPU calculations
             Maybe this one is not needed, since normally it
             would be 0. If NHG=0 then this is used, but
             basically we don't need the check below, since
             we can use the formula: NCH=numnodg-NGK-NGR
    NCE   Number of Cpus for Exclusion calculaitons
             We don't know yet about this... But possible extension


For example if we have 8 multicore boxes, with one GPU each, one can
use the following command to run CHARMM:

::

  mpirun -n 40 -hostfile $PBS_NODEFILE charmm -i input-script > output

But in the input script one then needs to specify the following:

::

  PARAllel GPUSplit ngr 8 nch 32

This will setup the parallel in the following order:

1. the first 8 CPUs will be controling GPUs
2. The 32 additional cores will do the rest of CHARMM calculations
   including PME

By adding one more box into the system one can perform PME on the
extra GPU:

::

  mpirun -n 41 -hostfile $PBS_NODEFILE charmm -i input-script > output

In the input script one then needs to specify the following:

::

  PARAllel GPUSplit ngr 8 ngk 1 nch 32

NOTE: one needs to be sure that pbs submit command reserves exactly 8
or 9 boxes in the above examples. To check this one can use the:

::

  PARAllel INFO

command in the begining of CHARMM input script. It tells which box is
executing which part of the calculation.




