.. py:module:: grid

Grid: A general facility to implment grid-based potentials for docking
======================================================================

# Charles L. Brooks III, TSRI. December 2000.

This document node describes the implementation, commands and syntax
associated with an implementation of grid-based potentials to be used
in ligand-docking studies, or when an additional set of potentials are
to be added to augment. It can be used with dynamics as well as the
GA/MC module.

.. _grid_implementation:

Grid-based potentials:  Description and Discussion
--------------------------------------------------


CHARMM modules involved: misc/grid.src, fcm/grid.fcm, fcm/energy.fcm
energy/energy.src, energy/eutil.src, energy/intere.src,
energy/printe.src, charmm/iniall.src, charmm/charmm_main.src

This module provides code to 1) generate a set of van der Waals and
electrostatic grid-based potentials and to 2) use these potentials in
dynamics, minimization and GA/MC-based searching algorithms.

Generation of the grid-based van der Waals potentials is accomplished
by establishing a series of vdW radius based potential surfaces over a
limited spatial extent specified by the user. This set of potentials
is built for radii of a series of test particles of unit epsilon
parameter. The general idea is to use radii that span the range of
radii used in the force field of interest, either on a discrete grid
or at particular values. In utilizing these grids for energy and force
calculations, the vdW radius of the atoms in the target molecule are
mapped to the nearest probe radius, with a warning being given if the
radii differ by more than 0.1 A, and the overall energy is scaled by
the square-root of the specific atoms vdW epsilon. This simplification
provides a means to minimize the number of 3D grids that must be
generated to represent the potential for a complex system. However, if
memory is not an issue, in principle grid-based potentials may be
generated for all vdW-based atom radii.

The electrostatic-based potential is the electrostatic potential
associated with a test charge throughout the user specified space.

The potential energy is computed as a 3D, 8-point linear interpolation
with forces computed from the analytic gradient of this interpolation
formula. The potential energy and forces beyond the grid edges is
constructed as a quadric potential away from the grid edge.


.. _grid_syntax:

Syntax for the Grid-based potentials
------------------------------------

Generation:

::

    grid generate select <atom selection> end   -
      [xcen <real>] [ycen <real>] [zcen <real>] -
      [xmax <real>] [ymax <real>] [zmax <real>] -
      [dgrid <real>] [Force <real>]             -
      [OutUnit <integer>] [Formatted] [Print]
     
Initiailization:

::

    grid read select <atom selection> end -
      Unit <integer> [Formatted] [Print]
 
    grid on select <atom selection> end 
 
    grid off
 
    grid clear


.. _grid_description:

There are two basic parts to utilizing grid-based potentials in CHARMM:

Description of the basic key words for grid-based potentials:
 
The following is the description of the setup commands for setting up the
system

===============   =======   ===============================================
Keyword/Syntax    Default   Purpose
===============   =======   ===============================================
GENErate                    Setting up the data structure and calculating
                            the potential grids.

READ                        Keyword to read a set of grid-based potentials
                            and set-up grid-based energy calculations.

ON                          Keyword to activate grid-based potential 
                            calculations.

OFF                         Keyword to de-activate grid-based potential 
                            calculations.

CLEAr                       Keyword to clear all grid-based potential
                            data structures from heap and stack.

XCEN              0.0 (A)   X-position for center of grid-based potentials.

YCEN              0.0 (A)   Y-position for center of grid-based potentials.

ZCEN              0.0 (A)   Z-position for center of grid-based potentials.

XMAX              0.0 (A)   X-direction extent of the potential grid.

YMAX              0.0 (A)   Y-direction extent of the potential grid.

ZMAX              0.0 (A)   Z-direction extent of the potential grid.

DGRId             0.5 (A)   Spacing between consequetive points in 
                            potential grids.

FORCe             300.0     Force constant for quadratic extention of 
                            potential beyond grid edges, in kcal/mol/A^2.

OUTUnit           Std Out   Unit to write grid-based potential file.

UNIT              Std Out   Unit from which to read grid-based potential file.

FORMatted         .false.   Logical to set reading/writing of grid-based
                            potentials in ascii format.

PRINt             .false.   Logical to set whether grid-based potentials
                            will be printed to standard out.
===============   =======   ===============================================

.. _grid_restrictions:

Restrictions
------------

This module is in alpha release and subject to change. All aspects
should work but this energy term has not been implemented in all other
CHARMM modules, e.g., it cannot be used with the free energy modules
pert, tsm or block, or with the MC module of Arron Dinner.
   

.. _grid_examples:

Supplementary examples
----------------------

::

   Generate and test grid for simple example of test atom.
   probes.RTF:
   * ...
   *
      22    0
   MASS   301 P1    1.00    P1 ! 
   MASS   302 P2    1.00    P2 ! 
   MASS   303 P3    1.00    P3 ! 
   MASS   304 P4    1.00    P4 ! 
   MASS   305 P5    1.00    P5 ! 
   MASS   306 P6    1.00    P6 ! 
   MASS   307 P7    1.00    P7 ! 
   MASS   308 P8    1.00    P8 ! 
   MASS   309 P9    1.00    P9 ! 
   MASS   310 P10   1.00    P10 ! 
   MASS   311 P11   1.00    P11 ! 
   MASS   312 P12   1.00    P12 ! 
   MASS   313 P13   1.00    P13 ! 
   MASS   314 P14   1.00    P14 ! 
   MASS   315 P15   1.00    P15 ! 
   MASS   316 P16   1.00    P16 ! 
   MASS   317 P17   1.00    P17 ! 
   MASS   318 P18   1.00    P18 ! 
   MASS   319 P19   1.00    P19 ! 
   MASS   320 P20   1.00    P20 ! 

   RESI PROB 20.000
   ATOM P1 P1 1.0 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20
   ATOM P2 P2 1.0 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20
   ATOM P3 P3 1.0 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20
   ATOM P4 P4 1.0 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20
   ATOM P5 P5 1.0 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20
   ATOM P6 P6 1.0 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20
   ATOM P7 P7 1.0 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20
   ATOM P8 P8 1.0 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20
   ATOM P9 P9 1.0 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20
   ATOM P10 P10 1.0 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20
   ATOM P11 P11 1.0 p12 p13 p14 p15 p16 p17 p18 p19 p20
   ATOM P12 P12 1.0 p13 p14 p15 p16 p17 p18 p19 p20
   ATOM P13 P13 1.0 p14 p15 p16 p17 p18 p19 p20
   ATOM P14 P14 1.0 p15 p16 p17 p18 p19 p20
   ATOM P15 P15 1.0 p16 p17 p18 p19 p20
   ATOM P16 P16 1.0 p17 p18 p19 p20
   ATOM P17 P17 1.0 p18 p19 p20
   ATOM P18 P18 1.0 p19 p20
   ATOM P19 P19 1.0 p20
   ATOM P20 P20 1.0
   END

   probes.prm:
   * Test probes for grid potential set-up
   *

   NBONDED  NBXMOD 5  ATOM RDIEL SWITCH VATOM VDISTANCE VSWITCH -
   CUTNB 999 CTOFNB 999 CTONNB 999 EPS 3 E14FAC 0.5 WMIN 1.5
   !
   !                 EMIN         Rmin     These columns used for
   !              (kcal/mol)      (A)      1-4 interactions
   !
   P1       0.00     -1.0000       0.65
   P2       0.00     -1.0000       0.75
   P3       0.00     -1.0000       0.85
   P4       0.00     -1.0000       0.95
   P5       0.00     -1.0000       1.05
   P6       0.00     -1.0000       1.15
   P7       0.00     -1.0000       1.25
   P8       0.00     -1.0000       1.35
   P9       0.00     -1.0000       1.45
   P10      0.00     -1.0000       1.55
   P11      0.00     -1.0000       1.75
   P12      0.00     -1.0000       1.85
   P13      0.00     -1.0000       1.95
   P14      0.00     -1.0000       2.05
   P15      0.00     -1.0000       2.15
   P16      0.00     -1.0000       2.25
   P17      0.00     -1.0000       2.35
   P18      0.00     -1.0000       2.45
   P19      0.00     -1.0000       2.55
   P20      0.00     -1.0000       2.65

   END



   * GRIDTEST.INP
   * This test-case demonstrates features of the grid-based potentials.
   * It utilizes the MSI CHARMm (Momany & Rone) force field and the
   * trypsin/benzamidine receptor/ligand pair.
   * Required files: MASSES.RTF, probes.RTF, AMINO.RTF, PARM.PRM, probes.prm
   *  3ptb_complex.psf, 3ptb_complex.pdb
   *

   open unit 1 read card name "MASSES.RTF"
   read rtf card unit 1

   open unit 1 read card name  "probes.RTF"
   read rtf card unit 1 append

   open unit 1 read card name "AMINO.RTF"
   read rtf card unit 1 append

   open unit 3 read card name "PARM.PRM"
   read param card unit 3


   open unit 1 read card name "probes.prm"
   read param  card unit 1 append

   open unit 1 read form name "3ptb_complex.psf"
   read psf card unit 1

   open unit 1 read form name "3ptb_complex.pdb"
   read coor pdb unit 1

   ! Find the center of the binding site
   coor stat select resname ptb end
   set xcen = ?xave
   set ycen = ?yave
   set zcen = ?zave

   ! Remove "real" ligand
   delete atom select resname ptb end

   ! Generate test probe atoms
   read sequ card
   * title
   *
   1
   prob
   generate  prob  setup

   ! Delete all atoms but single representative for first grid test
   delete atom select .not. ( type p15 .or. segid seg1 ) end

   ! Set-up position of test atom
   scalar x set @xcen select segid prob end
   scalar y set @ycen select segid prob end
   scalar z set @zcen select segid prob end

   ! Fix receptor atoms
   cons fix select segid seg1 end
   energy

   open unit 3 write form name grid.ascii
   title
   * Test grid for system
   *

   grid generate xmax 1 ymax 1 zmax 1 xcen @xcen ycen @ycen zcen @zcen -
        force 300 dgrid 0.5 select segid prob end outu 3 formatted print

   grid clear

   open unit 3 write unform name grid.bin
   title
   * Test grid for system
   *

   grid generate xmax 1 ymax 1 zmax 1 xcen @xcen ycen @ycen zcen @zcen -
        force 300 dgrid 0.5 select segid prob end outu 3 print

   grid clear

   open unit 3 read form name grid.ascii
   grid read unit 3 formatted select type p15 end print
   close unit 3

   grid clear

   open unit 3 read unform name grid.bin
   grid read unit 3 select type p15 end print
   close unit 3

   ! Generate positions on grid, vdW and elec should match grid terms
   energy inbfrq 0
   Calc Xmax = @Xcen + .5
   Calc Ymax = @ycen + .5
   Calc zmax = @zcen + .5
   Calc Xmin = @Xcen - .5
   Calc Ymin = @ycen - .5
   Calc zmin = @zcen - .5

   set x = @xmax
   label ix
     set y = @ymax
     label iy
        set z = @zmax
        label iz

          scalar x set @x select type p15 end
          scalar y set @y select type p15 end
          scalar z set @z select type p15 end
          energy
          Calc dvdW = ( ?vdW - ?Grvd ) / ?vdw
          Calc delec = ( ?elec - ?Grel ) / ?elec
   write title unit 12
   * ?Grvd ?Grel ?vdW ?elec @dvdw @delec
   *
          Calc z = @z - 0.5
        if z ge @zmin goto iz
        Calc y = @y - 0.5
     if y ge @ymin goto iy
     Calc x = @x - 0.5
   if x ge @xmin goto ix

   ! Test on/off components of grid energy terms
   grid off
   energy

   grid on select type p15 end
   energy

   skipe all excl grvd grel
   energy

   ! Generate energy curve along diagonal of cube to demonstrate interpolation
   ! and extrapolation.

   label dodiagonal
   Calc xlow = @Xmin - 0.5
   Calc x = @xmax+0.5
   Calc y = @ymax+0.5
   Calc z = @zmax+0.5
   set cnt = 0
   skipe all excl elec vdw grel grvd
   label diagonal
     scalar x set @x select type p15 end
     scalar y set @y select type p15 end
     scalar z set @z select type p15 end
     energy
     incr cnt by 1

   write title unit 13
   * @cnt ?Grvd ?vdW ?Grel ?elec
   *
     Calc z = @z - 0.1
     Calc y = @y - 0.1
     Calc x = @x - 0.1

   if x ge @xlow goto diagonal
   
   grid clear

   stop	

Example 2
^^^^^^^^^

An exploration of grid-based potential versus full molecular
potential for benzamidine-trypsin pair.

::

   * GRID_2.INP
   * This test-case demonstrates features of the grid-based potentials.
   * It utilizes the MSI CHARMm (Momany & Rone) force field and the
   * trypsin/benzamidine receptor/ligand pair.
   * Required files: MASSES.RTF, probes.RTF, AMINO.RTF, PARM.PRM, probes.prm
   *  3ptb_complex.psf, 3ptb_complex.pdb
   *

   open unit 1 read card name "MASSES.RTF"
   read rtf card unit 1

   open unit 1 read card name  "probes.RTF"
   read rtf card unit 1 append

   open unit 1 read card name "AMINO.RTF"
   read rtf card unit 1 append

   open unit 3 read card name "PARM.PRM"
   read param card unit 3

   open unit 1 read card name "probes.prm"
   read param  card unit 1 append

   open unit 1 read form name "3ptb_complex.psf"
   read psf card unit 1

   open unit 1 read form name "3ptb_complex.pdb"
   read coor pdb unit 1

   ! Define dimensions of volume for docking
   coor stat select resname ptb end
   set xcen = ?xave
   set ycen = ?yave
   set zcen = ?zave

   ! Set dimensions of grid as maximum extent of ligand + 4 A
   Calc Xmax = ?xmax - ?xmin + 4
   Calc Ymax = ?ymax - ?ymin + 4
   Calc Zmax = ?zmax - ?zmin + 4
   Let Xmax = Max @Xmax @Ymax 
   Let Xmax = Max @Xmax @Zmax

   ! If we have already generated the grid potentials go to final part.
   ! Uncomment after grid generation and run again.
   !goto alreadygener

   ! Remove ligand and generate probe atoms.
   delete atom select resname ptb end

   read sequ card
   * title
   *
   1
   prob
   generate  prob  setup

   ! Set positions for all probe atoms
   scalar x set @xcen select segid prob end
   scalar y set @ycen select segid prob end
   scalar z set @zcen select segid prob end

   ! Fix position of receptor.
   cons fix select segid seg1 end
   skipe all excl vdw elec

   energy 

   open unit 3 write unform name grid_3ptb.bin
   title
   * Test grid for system
   *

   ! Generate grid-based potentials for 20 probe atoms + electrostatic
   ! using default grid spacing of 0.5 A and default harmonic potential
   ! beyond grid edges (300 kcal/mol/A^2).
   grid generate xmax @xmax ymax @xmax zmax @xmax -
        xcen @xcen ycen @ycen zcen @zcen -
        select segid prob end outu 3

   grid clear

   stop

   ! Begin here after grid potentials have been generated
   label alreadygener

   ! Fiex receptor atoms for "rigid"-receptor docking
   cons fix select segid seg1 end

   ! Read grid and set-up for ligand (seg2)
   open unit 3 read unform name grid_3ptb.bin
   grid read unit 3 select segid seg2 end 
   close unit 3

   ! Randomly rotate ligand about its center and minimize
   Calc phi = ?rand * 30

   coor rota xdir @xcen ydir @ycen zdir @zcen phi @phi select segid seg2 end

   ! Turn off grid potential and minimize using "true" receptor.
   grid off
   energy inbfrq 1
   coor copy compare
   mini sd nstep 200 inbfrq 0
   coor rms select segid seg2 end

   ! Turn on grid potential, restore coordinates of ligand and remove receptor
   ! then minimize using grid-based potential only.
   grid on select segid seg2 end
   coor swap
   coor translate xdir 10000 select segid seg1 end
   energy inbfrq 1
   mini sd nstep 200 inbfrq 0

   ! Check rmsd between ligand minimized in actual receptor and in grid-based
   ! receptor.
   coor rms select segid seg2 end

   stop

