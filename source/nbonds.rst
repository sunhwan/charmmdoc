.. py:module:: nbonds

#####################################
Generation of Non-bonded Interactions
#####################################

Nonbonded interactions (frequently abreviated "nbond") refer to
van der Waals terms and the electrostatic terms between all atom pairs
that are not specifically excluded from nonbond calculations
as for example are directly bonded atoms :ref:`nbx <struct_nbx>`.
These terms are defined on atom pairs and to a first aproximation would
require the number of atoms squared calculations. To avoid this burden
various truncation and approximation schemes can be employed in the
program, breaking the nonbonded calculation into two parts,
initialization and actual energy calculation.

The method of approximation, cutoffs, and other relevant
parameters can be entered any time the nbond specification parser is
invoked. See the syntax section for a list of all commands that invoke this
parser.

Simple Ewald function is modified so it works with any shape of the 
simulation box. Comments/suggestions at hkamberaj@asu.edu

.. _nbonds_syntax:

Syntax
------

::

   { NBONds       }   { [INBFrq integer] nonbond-spec  }
   { UPDAte ...   }   {                                }
   { ENERgy ...   }   {                                }
   { MINImize ... }   {                                }
   { DYNAmics ... }   {                                }
   { VIBRAN ...   }   {                                }
   { CORRel ...   }   {                                }

.. note::
   The INBFrq value is remembered. If its value is zero,
   no interpretation of [nonbond-spec] will be made, as well as no
   modifications of the nonbond lists. It's default value is -1 .

In all cases as many keywords and values as desired may be specified.
The keywords are:

:: 

   nonbond-spec::= [method-spec] [distances-spec] [misc-specs] [INIT] [RESET]

   method-spec::= [ ELEC electrostatics-spec ] [ VDW vdw-spec  ] [ BYCUbes ]
                  [ NOELectrostatics         ] [ NOVDwaals     ] [ BYGRoup ]
                  [ GRAPe grape-spec         ] [ LIST          ] [ BYCBim  ]
                  [ NOGRape                  ] [ NOLIst        ] [ BYCC    ]
                  [ LRC                      ]
                  [ IPS  ips-spec            ]
           
   electrostatics-spec::= [ EWALD ewald-spec ] [ FMA fma-spec ]  elec-opt-spec 
                          [ NOEWald          ] [ NOFMA        ] [ ACE ace-spec ]
                          [ EIPS          ] 

   elec-opt-spec::= 
      [ ATOM                                     ] [ CDIElec ] [ SHIFted  ]
      [ GROUp [ EXTEnded [GRADients] [QUADrip] ] ] [ RDIElec ] [ SWITched ]
      [       [          [NOGRad   ] [NOQUads] ] ]             [ FSWItch  ]
      [       [ NOEXtended                     ] ]             [ FSHIft   ]
                                                               [ MSHIft   ] ! MMFF
                                                               [ TRUNcate ] ! MMFF

   vdw-spec::=  [ VGROUP [ VSWITched ]           ]
                [ VATOM  [ VSHIfted           ]  ]
                [        [ VSWItched          ]  ]
                [        [ VFSWitch           ]  ]
                [        [ VTRUnc [ CTVT real ]  ]  ! MMFF only
                [        [ VIPS               ]  ]  ! NBIPS

   ips-spec::=  [RAIPS real] [RIPS real] [NIPSFRQ int] [DVBIPS real] -
                [ [MIPSX int] [MIPSY int] [MIPSZ int] [MIPSO int]] -
                    [ [PXYZ]                    ]
                    [ [PXY|PYZ|PZX|PYX|PZY|PXZ] ] 
	                  [ [PX|PY|PZ]                ]    
	         

   distances-spec::= [general-dist] [warning-dist]

   general-dist::=
           [ CUTNB  real ] [ CTONNB real ] [ CTOFNB real ] [ CTEXNB real ]
           SOFT [EMAX  real ] [MINE real] [MAXE real] VDWE ELEE  ! SOFTVDW 

   warning-dist::=
           [ WMIN   real ] [ WRNMXD real ]

   misc-specs::= [ EPS real ] [ E14Factor real ] [ NBXM integer] -
                    [ NBSCale real] [IMSCale real] [EXOForce]  -

                        [ NORXN                ]
                        [ RXNFLD  rxnfld-spec  ]
                        [ RXNNB   rxnfld-spec  ]

   rxnfld-spec::= [ EPSEXT real ] [ ORDER integer ] [ SHELL real ]

   fma-spec::= [LEVEL int] [TERMS int]

   ewald-spec::=  (see ewald::(chmdoc/ewald.doc)syntax )
   ace-spec::=    (see ace::(chmdoc/ace.doc)syntax )


.. _nbonds_defaults:

Defaults
--------

The defaults for the nonbond specification reside with the
parameter file. The defaults are specified at the begining of the van der
Waal section. These defaults are the recommended options.

The following command contains all defaults for one of the
older protein parameter files, and is equvalent to the command :chm:`NBONDS INIT`
in it usage when this parameter file is present.

::

   NBONDS ELEC ATOM NOEX NOGR NOQU SWIT RDIE  VATOM VDW VSWI -
          CUTNB 8.0  CTEXNB 999.0 CTOFNB 7.5 CTONNB 6.5 WMIN 1.5 WRNMXD 0.5 -
          EPS 1.0 NORXN EPSEXT 80.0 ORDER 10 SHELL 2.0  CTVTRN 10.0 -
          E14FAC 1.0 NBSCAL 1.0 IMSCAL 1.0 NBXMOD 5 NOFMA -
          NOEWALD NOPME KAPPA 1.0 KMAX 5 KMAXSQ 27 ERFMOD -1

MMFF specific defaults: :chm:`VTRUNC MSHFT E14FAC 0.75 CTVTRN 8.0`

SOFTVDW specific defaults:

::

   SOFT EMAX 30.0/EPS MINE -300.0/EPS MAXE -2.*MINE ! CDIE
   SOFT EMAX 15.0/EPS MINE -120.0/EPS MAXE -2.*MINE ! RDIE

Values do not change unless explicitly specified, except for
the ON/OFF values which cascade when the cutoff values are changed as;

::

   CTOFNB=CUTNB-0.5
   CTONNB=CTOFNB-1.0

.. warning::

   These old defaults have been shown to be detrimental to protein
   behavior.  It is generally better to use the defaults in the parameter sets.
   
The presence of soft core nonbonded terms is recommended for
calculations in a dense system (docking, loop refinement, NMR refinement).
The initial value of RMIN (switching distance for the soft core potential)
is recommended to be 0.885  
 

RECOMMENDED:

Presented here is a suggested list of options. Where specifications are
missing, substitute the defaults:

Use Isotropic Periodic Sum (IPS) calculation for either
finite (in vacuum) or periodic systems:

NBONDS  IPS  CUTNB 12.0  CTOFNB 10.0   EPS 1.0 CDIE

or use IPS with fully homogenouse assumption for either
finite (in vacuume) or periodic systems::

NBONDS  IPS PXYZ CUTNB 12.0  CTOFNB 10.0   EPS 1.0

or use IPS with 2D homogenouse assumption for interfacrial membrane systems:

NBONDS  IPS PXY CUTNB 12.0  CTOFNB 10.0   EPS 1.0

or use VIPS for vdw and Ewald for charge interaction for periodic systems:

NBONDS  VIPS   -
        ATOM  EWALD PMEWALD KAPPA 0.32  FFTX 32 FFTY 32 FFTZ 32 ORDER 6 -
        CUTNB 12.0  CTOFNB 10.0   EPS 1.0 CDIE  



For no-cutoff periodic systems:
                                           **      **      **
NBONDS ATOM EWALD PMEWALD KAPPA 0.32  FFTX 32 FFTY 32 FFTZ 32 ORDER 6 -
        CUTNB 12.0  CTOFNB 10.0  VDW  VSHIFT  

(** system size dependent - use about 0.8-1.0 grids per angstrom)

For atom based cutoffs:

NBONDS  ATOM  FSHIFT CDIE  VDW VSHIFT  -
        CUTNB 13.0  CTOFNB 12.0 CTONNB 8.0  WMIN 1.5  EPS 1.0

or (perhaps better for longer cutoff distances, but more expensive)

NBONDS  ATOM  FSWITCH CDIE  VDW VSHIFT  -
        CUTNB 13.0  CTOFNB 12.0 CTONNB 8.0  WMIN 1.5  EPS 1.0

For group based cutoffs (doesn't vectorize well):

NBONDS  GROUP  FSWITCH CDIE  VDW VSWITCH  -
        CUTNB 13.0  CTOFNB 12.0 CTONNB 8.0  WMIN 1.5  EPS 1.0

For extended electrostatics :

NBONDS  GROUP  SWITCH CDIE  VDW VSWI  EXTEND GRAD QUAD -
        CUTNB 13.0  CTOFNB 12.0 CTONNB 8.0  WMIN 1.5  EPS 1.0

For a better description of these methods and how they perform, see:
 P.J. Steinbach, B.R. Brooks: "New Spherical-Cutoff Methods for Long-Range
 Forces in Macromolecular Simulation,"  J. Comp. Chem. 15, 667-683 (1994).

OPTIONS THAT ARE NOT RECOMMENDED (OR REALLY BAD):

[ ATOM  ] [ CDIElec ] [ SHIFted  ]   no (obsolete, but used in the past)
[ ATOM  ] [ CDIElec ] [ SWITched ]   NO!  Very bad - do not use!
[ GROUp ] [ CDIElec ] [ SHIFted  ]   no (obsolete)
[ GROUp ] [ CDIElec ] [ SWITched ]   NO!  Very bad with non-neutral groups!
[ ATOM  ] [ RDIElec ] [ SHIFted  ]   yes, but do you really want RDIE??
[ ATOM  ] [ RDIElec ] [ SWITched ]   no. switch is bad here.


File: Nbonds, Node: Function, Up: Top, Previous: Defaults, Next: Tables

NBSCale & IMSCale
=================

   The first time that the primary or image non-bond list is generated, an
estimate is made, based on the number of atoms, of how much memory will be
needed to store the pair list.  If too large an estimate is made, memory will
be wasted.  If too small an estimate is made, a second (and larger) estimate
will be made and the memory allocated on the first attempt is wasted.  NBSCale
is a correction factor to the intial estimate allowing better control of 
memory allocation.  For example NBSCale 1.5 allocated 50% more memory than
CHARMM would usually guess and NBSCale 0.8 allocated 20% less.  IMSCale does
the same thing when the image pair list is generated.  The values of NBSCale
and IMSCale must be determined empirically, but they can generate huge
memory savings on large systems.

These keywords are valid wherever nonbond options may appear, e.g. ENERgy,
DYNAmics, MINImiz, and UPDAte.  Note that NBSCale must be used in the first 
statement which generates a nonbond list; an UPDAte without NBSCale followed 
by DYNAmics with NBSCale is ineffective. 

For a system of about 17,000 atoms, a value of NBSCALE 1.5 was effective in
providing about 25 MB of memory reduction compared to the the default NBSCale
value (1.0); while for a system of about 10,000 atoms, the optimum of 1.4 gave
a reduction of about 12 MB.  For example:

mini abnr nstep 500 nprint 10 tolenr 1.e-6 cutnb 14. ctofnb 12. ctonnb 10. -
 fshift cdie eps 1.0 vshift inbfrq 20 imgfrq 20 cutim 14. nbscale 1.5

Likewise for DNYAmics, ENERgy, or the UPDAte commands; the latter is useful in
doing some trial and error probes to determine the optimum NBSCale value,
with fixed sizes for CUTNB and CUTIM:

update cutnb 14. ctofnb 12. ctonnb 10. fshift cdie eps 1.0 vshift -
  inbfrq 20 imgfrq 20 cutim 14. nbscale @1 -

where the csh or tcsh command line might be something like:

% charmm medium 1:1.5 < tst_nbscal.inp >& nbscal_1.5 &

The above is based on single processor calculations; the same general idea
applies to parallel calculations, but the optimum value for NBSCale will
be less than 1.0, perhaps 0.7 to 0.8 for systems in the 10K to 17K atom
range.  Memory usage can be further reduced using the IMSCale keyword; 
some experimentation will be required depending on the number of atoms
and the cutoffs being used.

INBFrq n
========

   Update frequency for the non-bonded list. Used in the subroutine ENERGY()
to decide whether to update the non-bond list. When set to :

     0 --> no updates of the list will be done.

    +n --> an update is done every time  MOD(ECALLS,n).EQ.0  . This is the old
           frequency-scheme, where an update is done every n steps of dynamics
           or minimization.

    -1 --> heuristic testing is performed every time ENERGY() is called and
           a list update is done if necessary. This is the default, because
           it is both safer and more economical than frequency-updating.

Description of the heuristic testing algorythm :
-----------------------------------------------
   Every time the energy is called, the distance is computed each atom moved
since the last list-update.
   If any atom moved by more than (CUTNB - CTOFNB)/2 since the last list-update
was done, then it is possible that some atom pairs in which the two atoms are
now separated by less than CTOFNB are not in the pairs-list. So a list update
is done.
   If all atoms moved by less than (CUTNB - CTOFNB)/2 , then all atom
pairs within the CTOFNB distance are already accounted for in the non-bond list
and no update is necessary.

Description of the code for the heuristic testing :
--------------------------------------------------
   This section describes how programmers can control the list-updating
behavior when their routines call the ENERGY() subroutine.
   All list-updating decisions, whether they are frequency based or heuristic
based, are made in the subroutine  UPDECI(ECALLS) , which is called from only
one place : at the very beginning of ENERGY().
   UPDECI(ECALLS) can be controled through INBFRQ (via the CONTRL.FCM common
block) and ECALLS (via the ENERGY.FCM common block) as follows :

     If INBFRQ = +n --> non-bond list is performed when MOD(ECALLS,n).EQ.0 .
                        Image and H-bond lists are updated according to IMGFRQ
                        and IHBFRQ.

     If INBFRQ =  0 --> non-bond list update is not performed. Image and H-bond
                        lists are updated according to IMGFRQ and IHBFRQ.

     If INBFRQ = -1 --> all lists are updated when necessary (heuristic test).

     (note that ECALLS is incremented by ICALL every time ENERGY(,,,ICAL) is
      called. In most cases, ICALL=1 )

   The current implementation of UPDECI() will work (without modifications) to
decide whether the image/crystal non-bond lists need updating, provided the
periodicity parameters don't change (i.e. constant Volume).
    UPDECI() is easily adapted to variable Volume dynamics/minimizations. This
is described in comments of the routine itself.

Further computational economy in update-testing :
-------------------------------------------------
   A programmer can sometimes skip the heuristic test itself, making the
decision whether to do list-updating even more economical.
   This option is only available if the size of the step taken since the last
call to ENERGY() is known. For an example of usage, see the subroutine ENERG()
in TRAVEL.


NON-BOND ENERGY TERMS.
======================
        The electrostatic options are separate from the van der Waal
options, though some keywords are shared between them. The following is
a description of all options and keywords.

1) Electrostatics

        The ELEC keyword (default) invokes electrostatics. The NOELec
keyword supresses all electrostatic energy terms and options.
There are two basic methods for electrostatics, GROUp and ATOM.  A 
model based on the GROUp method is the extended electrostatics model 
which approximates the full electrostatic interaction and eliminates
the need for a cutoff function. This model is based on the partitioning
of the electrostatic term into two contributions.  One comes from the
interaction between particles which are spatially close and is treated
by conventional pairwise summation. The second contribution
comes from interactions between particles which are spatially distant from
one  another and is treated by a multipole moment approximation.
[The original model was described in B. R. Brooks, R. E. Bruccoleri,
B. D.  Olafson, D. J. States, S. Swaminathan, M. Karplus.  J. Comp. Chem.,
4, 187, (1983) and more recently in R.H. Stote, D.J. States and M. Karplus,
J. Chimie Physique (1991)]


        A) Functional Forms

Atom electrostatics indicates that interactions are computed on an
atom-atom pair basis. There are two options that specify the radial
energy functional form. The keywords CDIE and RDIE select the basic
functional form. The SWIT, SHIF, FSWI, and FSHI keywords determine the
long-range truncation option.

CDIE - Constant dielectric. Energy is proportional to 1/R.
RDIE - Distance dielectric. Energy is proportional to 1/(R-squared)

SWIT - Switching function used from CTONNB to CTOFNB values.
SHIF - Shifted potential acting to CTOFNB and zero beyond.
FSWI - Switching function acting on force only.  Energy is integral of force.
FSHI - Classical force shift method for CDIE (force has a constant offset).

EIPS - Isotropic periodic sum method for CDIE or RDIE electrostatic energy

        B) Atom electrostatics

Atom electrostatics indicates that interactions are computed on an
atom-atom pair basis. There are two options that specify the radial
energy functional form. The keywords CDIE and RDIE select the basic
functional form. The SWIT and SHIF keywords determine the long-range
truncation option.

[ ATOM ] [ CDIElec ] [ SHIFted  ]
         [ RDIElec ] [ SWITched ]
                     [ FSWItch  ]
                     [ FSHIft   ]
                     [ EIPS     ]

CDIE - Constant dielectric. Energy is proportional to 1/R.
RDIE - Distance dielectric. Energy is proportional to 1/(R-squared)

SWIT - Switching function used from CTONNB to CTOFNB values.
SHIF - Shifted potential acting to CTOFNB and zero beyond.
FSWI - Switching function acting on force only.  Energy is integral of force.
FSHI - Classical force shift method for CDIE (force has a constant offset)

EIPS - Isotropic Periodic Sum for charge interaction.

        C) Group Electrostatics
electrostatics-spec::=
                                                              
[ GROUp [ EXTEnded [ GRADients ] [ QUADrip ] ] ] [ CDIElec ] [ SWITched ]
[       [          [ NOGRad    ] [ NOQUads ] ] ] [ RDIElec ] [ FSWItch  ]
[       [ NOEXtended                         ] ]                         
[       [ EIPS                               ] ]                         

SWIT - Switching function used from CTONNB to CTOFNB values.
FSWI - Switching function, but QiQj/Rcut is added before switching.
         (FSWI has no effect on neutral groups).
EIPS - Isotropic Periodic Sum using CTOFNB as the local region radius.


        D) Electrostatic Distances
electrostatic-dist::=
        [ CUTNB  real ] [ CTEXNB real ]        [ CTONNB real ] [ CTOFNB real ]
        [ EMAX   real ] [ MINE real ] [ MAXE real ] ! SOFTVDW only

EMAX - Twice the VDW energy from which the soft core potential becomes active
       (if SOFT key word is employed).
       It has a linear form for r<rc (E(rc)=EMAX/2) :
       Esoftl=EMAX/2+alfa*(r-rc)
       unless the VDWE key word turns on the exponential form :
       EsoftE=EMAX+alfa*r**beta  
       For the exponential form E(0)=EMAX, so EMAX is also the VDW energy
       at r=0 for the exponential form  

MINE - Twice the energy from which the electrostatic attractive soft  potential
       begins. ELEE turns on the exponential form.

MAXE - Twice the energy from which the electrostatic repulsive soft potential
       begins.  ELEE turns on the exponential form.  

       The soft core potential is currently implemented
       only in fast energy routine, so fast option has to be used to
       allow it. The form of the soft potential for the electrostatics is
       the same as for the VDW, the defaults are different for VDW,
       electrostatic repulsion and attraction.  The defaults need to be
       modified for some cases (i.e. for the spc water model) to prevent
       shifting of energy minima.
       The SOFT key word turns on both VDW and electrostatic soft core. 
       To turn off the soft core set EMAX = 0 with the SOFT keyword (in energy,
       minimization or nbonds call).  

CTEXNB - defines the cutoff distance beyond which interaction pairs are
         excluded from the Extended Electrostatics calculation.


        E) Extended (group) Electrostatics
electrostatics-spec::=
[ ATOM                                         ] [ CDIElec ] [ SHIFted  ]
[ GROUp [ EXTEnded [ GRADients ] [ QUADrip ] ] ] [ RDIElec ] [ SWITched ]
[       [          [ NOGRad    ] [ NOQUads ] ] ]
[       [ NOEXtended                         ] ]

EXTE - invokes the extended electrostatics command for calculating long
       range electrostatic interactions.
NOEX - suppress the extended calculation.
GRAD - keyword flags the inclusion of the field of the extended gradient in
       calculating the force on each atom,i.e. include first and second
       derivatives.
QUAD - flags the inclusion of the quadrupole in the multipole expansion.


        F) Reaction Fields
misc-specs::= [ EPS real ] [ E14Factor real ] [ NORXN                ]
                                              [ RXNFLD  rxnfld-spec  ]
                                              [ RXNNB   rxnfld-spec  ]

rxnfld-spec::= [ EPSEXT real ] [ ORDER integer ] [ SHELL real ]

        G) Isotropic Periodic Sum
electrostatics-spec::={EIPS] [ ATOM   ] [ CDIElec ]
                             [ GROUp  ] [ RDIElec ] 


2) Van Der Waal Interactions

The VDW keyword (default) invokes the van der waal energy term.
To supress this term, the NOVDw keyword may be used.

        A) Distance specified van der Waal Function
             [ VGROUP [ VSWITched ]        ]
             [        [ VIPS      ]        ]
vdw-spec::=  [ VATOM  [ VSHIfted  ]        ]
             [        [ VSWItched ]        ]
             [        [ VFSWitch  ]        ]
             [        [ VIPS      ]        ]

VIPS - Isotropic Periodic Sum for VDW interaction.

3) Miscellaneous options and keywords

    A) Dielectric specification
misc-specs::= [ EPS real ] [ E14Factor real ]


    B) Warning Distance Specifications
warning-dist::=
        [ WMIN   real ] [ WRNMXD real ]

 WRNMXD - keyword defines a warning cutoff for maximum atom displacement from
          the last close contact list update (used in EXTEnded)


    C) Initialization

In all cases as many keywords and values as desired may be specified.
The key words, their functions, and defaults are:

        1) The method to be used

        2) Distance cutoff in generating the list of pairs

                CUTNB value (default 8.0)

        3) Distance cut at which the switching function eliminates
                all contributions from a pair in calculating energies.
                Once specified, This value is not reset unless respecified.

                CTOFNB value (default CUTNB-0.5)

        4) Distance cut at which the smoothing function begins to reduce
                a pair's contribution. This value is not used with SHFT.
                Once specified, This value is not reset unless respecified.

                CTONNB value (default CTOFNB-1.0)

        6) Dielectric constant for the extened electrostatics routines
                (RDIE option sets the dielectric equal to r
                times the EPS value)

                EPS value (default 1.0 for r dielectric)
                EPS 0.0  or  NOELec  (zero elecrostatic energy)

        7) Warning cutoff for minimum atom to atom distance. Pairs are
                checked during close contact list compilation.

                WMIN value (default 1.5)

        8) Warning cutoff for maximum atom displacement from the last
                close contact list update (used only in EXTEnded)

                WRNMXD value (default 0.5)
        9) The presence of soft core nonbonded interactions (turned on
           only if SOFT key word is present)


    D) Exclusion Lists

By default, vdw and electrostatic interactions between two bonded
1-2 interactions) and two atoms bonded to a common atom (1-3 interactions)
atoms are excluded from the calculation of energy and forces.  Also,
special vdw parameters and an electrostatic scale factor (E14FAC) can be
applied to atom pairs separated by 3 bonds (1-4 interactions).  The
control of the exclusion list is by the integer variable, NBXMod.

   NBXMod =     0        Preserve the existing exclusion lists
   NBXMod = +/- 1        Add nothing to the exclusion list
   NBXMod = +/- 2        Add only 1-2 (bond) interactions
   NBXMod = +/- 3        Add 1-2 and 1-3 (bond and angle)
   NBXMod = +/- 4        Add 1-2 1-3 and 1-4'S
   NBXMod = +/- 5        Add 1-2 1-3 and special 1-4 interactions

A positive NBXMod value causes the explicit exclusions in the PSF (inb array)
to be added to the exclusion list.  A negative value causes the use of only
the bond connectivity data to construct the exclusion list (thus, ignoring
the PSF data).

   E) LRC 
          Long range correction to cutoff van der Waal's energy and
          its contribution to the pressure.

   F)  Isotropic Periodic Sum (IPS)
       This is a newly developed method to calculate electrostatic and/or
VDW interaction accurately and efficiently. Using EIPS, VIPS to setup IPS 
calculation for electrostatic and vdw interactions, respectively. Both atom
and group nonbonded list are supported. Alternative, one can set IPS for both
electrostatic and vdw interactions. Also, one can use VIPS for vdw and use
Ewald for charge interaction. The IPS method can be applied to both periodic
and finit(in vacuume) systems. The IPS method calculate long-range
interactions in two steps. The first step calculates the interaction with
the local region defined by CTOFNB or RIPS. The second step calculates
the difference of an anisotropic system (defined by radius RAIPS, which is
set by default the diagonal distance of the PBC box. The first step is done
the same way as the cutoff methods by summing over local atom pairs.
The second step is done through the convolution thereom which can be
efficiently calculated using FFT technique. This method can be used for both
homogenouse and hetergenous systems. By setting PXYZ, the second step will
be turned off by assuming the system is fully homogenouse. Setting
PXY|PYZ|PXZ or PX|PY|PZ will use 1+2D IPS for the second step calculation.
Other than default values determined by system sizes, grid numbers for FFT
can be set through MIPSX, MIPSY,MIPSZ, and Bspline order by MIPSO. If use
VIPS with EWALD (PME), grid numbers and bspline order will be defined by PME
input. For constant pressure simulation, the IPS energies at grid points are
updated according to the updating frequency, NIPSFRQ (default 1), and the
volume change ratio, DVBIPS(default 10e-9). Increase NIPSFRQ or DVBIPS can
slightly speed up simulation.
     For finite systems( such as in vacuume), the calculation is done by
assuming a PBC box that is large enough (twice the size of the system) so
that molecule will not see any images within the interaction range.
     The original description of the IPS method can be found at:
"Xiongwu Wu and Bernard R. Brooks, Isotropic Periodic Sum: A method for the 
calculation of long-range interactions. J. Chem. Phys., Vol.122, No.4, 
article 044107 (2005)" (http://link.aip.org/link/?JCP/122/044107/1)

Here are some examples to using IPS.  Using IPS for any simulation system:

DYNA LEAP CPT STRT  NSTE 100 TIME 0.001 -
   EIPS VIPS  -
   CUTNB 12 CTOFNB 10 imgfrq 10 inbfrq 10

or use IPS for vdw and Ewald for charge interaction:

DYNA LEAP CPT STRT  NSTE 100 TIME 0.001 -
   VIPS  -
   ATOM  EWALD PMEWALD KAPPA 0.32  FFTX 32 FFTY 32 FFTZ 32 ORDER 6 -
   CUTNB 12 CTOFNB 10 IMGFRQ 10 INBFRQ 10 NTRFRQ 100

or for fully homogenous systems:

DYNA LEAP CPT STRT  NSTE 100 TIME 0.001 -
   EIPS VIPS PXYZ -
   CUTNB 12 CTOFNB 10 IMGFRQ 10 INBFRQ 10

or for interfacial systems using 1+2D IPS:

DYNA LEAP CPT STRT  NSTE 100 TIME 0.001 -
   EIPS VIPS PXY -
   CUTNB 12 CTOFNB 10 IMGFRQ 10 INBFRQ 10

___________________________________________________________________


ALGORITHMS

        There are four algorithms used in calulating the nonbonded
energies, each making different approximations in an attempt to speed
the calulation. Electrostatic interactions are the most difficult to
deal with for two reasons. They do not fall off quickly with distance
(so it is inappropriate to simply ignore all interactions beyond some
cutoff), and they depend on odd powers of r necessitating expensive
square root caluculations for each pair evaluated. The approximations
used to make the electrostatics calculation more tractable are setting
the dielectric constant equal to r or using a constant dielectric but
only calculating distant interactions periodically (and storing the
value in between).

        Setting the dielectric constant equal to the atom atom distance
times a constant factor ( determined by the EPS keyword value )
makes the computation easier by eliminating the need to calculate square
roots and by making the calculated contribution fall off more quickly.
It also introduces problems. The force calculated using an r dependent
dielectric will be larger than the force from a constant dielectric at
short distances (5.0 angstroms or less by comparsion to a constant
dielectric of 2.5). In addition, the electrostatic contribution still
falls off relatively slowly and large distance cutoffs are needed. As
the number of atom pairs included will be proportional to the cutoff
cubed, this is a significant disadvantage.

        The SHIFt option is similar to SWITch except, the potential:

        E = (QI*QJ/(EPS*R)) * ( 1 - (R/CTOFNB)**2 )**2

is used when ( R < CTOFNB ) and zero otherwise. This potential
and it first derivative approach zero as R becomes CTOFNB,
without the messy computation of switching functions and steep
forces at large R. 

CDIE uses a constant dielectric everywhere.  This requires a square root
to be calculated in the inner loop of ENBOND, slowing things down a bit,
but it is physically more reasonable and widely employed by other groups
doing empirical energy modelling.  The short range forces are
identical to those calculated with the other options, reflecting
the decrease in dielectric shielding at short ranges.

        The constant dielectric routines compile the close contact list
using the same two stage minimum rectangle box search that is described
above. In this way the efficiency of a residue by residue search is
exploited while being certain that all necessary pairs are included. For
close residue pairs an atom by atom search is then performed. Atom pairs
are either included in the list of close contacts or their electrostatic
interactions are calculated and stored.

Description of the Extended Electrostatics method
-------------------------------------------------
For the long range forces there is effectively no cutoff in the
electrostatic energy when using the Extended Electrostatics model.
The Extended Electrostatics model approximates the full electrostatic
interaction by partitioning the electric potential and the resulting forces
at any point ri into a near and extended contribution. The near contribution
arises from the charged particles which are spatially close to ri while the
extended contribution arises from the particles which are spatially distant
from ri. The total electrostatic potential can be written as a sum of the
two.  The near region is defined in terms of a radial distance, CUTNB, for
each atom.  Interactions between atoms separated by a distance greater than
CUTNB are calculated using a time saving multipole approximation when the 
nbond list is updated.  These interactions are stored together with their 
first (NOGRad) or first and second (GRADients) derivatives.  Interactions
between particles within CUTNB are calculated by the conventional pairwise 
additive scheme. (For a more complete development of the model, see
R.H. Stote, D.J. States and M. Karplus, J. Chimie Physique Vol. 11/12, 1991).
The energy is calculated by explicitly evaluating pairs in the list and using 
the stored potentials, fields, and gradients to approximate the distant
pairs. In essence the routines assume that for distant pairs the atom
movements will be small enough that the changes in their electrostatic
interactions can be accurately calculated using local expansions.
In using this model the GROUp method for constructing the nonbond list must
be used.  The interactions between particles within CUTNB are truncated 
rather than having a SHIFt or SWITch function applied.  Additionally, as
one is calculating all electrostatic interactions in the system, the 
dielectric constant should be set to 1.0.

Not Available at this time:
        An option is offered to increase the accuracy of residue residue
interactions by using a multipole expansion of one residue evaluated for
each atom of the other. This cutoff for this treatment is CUTMP. For
residue pairs outside of CUTMP only a single multipole evaluation is
made and second order polynomial expansion is used to extrapolate to
each atom.  Ordinarily this is sufficient and CUTMP is set to 0.0.

IMPLEMENTATION AND DATA STRUCTURES

        The initialization and list compilation is performed by the
subroutine NBONDS which in turn will call a lower level routine that
will do all of the work.  It functions by guessing how much space will be
needed to store the close contact list, allocating that space (and space
for electrostatic potentials and gradients if necessary) on the heap,
and calling the appropriate subroutine to actually compile the nonbonded
list (NBONDG,...). If sufficient space was not available 1.5 times as
much is allocated and another attempt is made.

        ENBOND evaluates the nonbonded energy, calling EEXEL to evaluate
the stored electric potentials and fields. Double precision is used for
all arithmatic.

        All of the nonbonded cutoffs and lists are stored on the heap.
BNBND is the descriptor array passed through most of the program (in
some of the analysis routines an additional array BNBNDC is used for the
comparision data structure). BNBND holds heap adresses and LNBND holds
the lengths of the elements in the data structure. To actually access
the data it is necessary to include INBND.FCM (an index common
block) and specify HEAP(BNBND(xxx)) where xxx is the desired element
name in INBND.FCM. This is arrangement has the advantage of allowing
dynamic storage allocation and easy modification of the types of
information passed from routine to routine.

The nonbonded data structure is described in: source/fcm/inbnd.fcm


FAST VECTOR/SHARED-MEMORY PARALLEL ROUTINES

      There are 5 sets of standard routines used to compute nonbond energies
and forces. The can be summarized:

      FAST CRAYVEC - optimized vector code for a Cray (array processor)
      FAST PARVECT - optimized vector code for an SMP machine
      FAST VECTOR  - general optimized vector code
      FAST SCALAR  - general optimized scalar code
      FAST OFF     - the generic - support everything routine

These routines are processed in the order listed.  The highest gaining priority
based on what options have been compiled and what the user requests.
If the user does not specify a fast option and all code is compiled,
then the CRAYVEC code will probably be used, unless this routine does
not support the requested options (in which case, the next routine is tried).
For example, the PARVECT code supports PME Ewald, but CRAYVEC does not, so
a calculation with PME will run with PARVECT (unless otherwise specified).
To determine which routine is actually doing your calculation, use "PRNLEV 6"
to list energy routines as they are processed.


OTHER

The option EXOForce, forces update of exclusion lists.
The option SOFT will turn on the soft core nonbonded potential.

The option GRAPe will perform all the nonbond interactions in the
specialized molecular dynamics hardware, called GRavity PipE.
There is an environment (M2_ON) variable to specify which board to
use:
Example for MDGRAPE2:
envi m2_on 0 (use 1st board)
envi m2_on 1 (use 2nd board)
envi m2_on "0,1" (use both: 1st and 2nd)
If you don't specify this environment variable then it will use all
available boards.
This variable is changed for MDGRAPE3 to MR3_BOARD.
grap-spec is an integer variable (default -1 = don't use GRAPE)
If zero then normal usage.
NOTE for GRAPE: Because of the energy calculation is done only whe it
is printed in the output, default ECHECK is too low. Please increase!

The option NOLIst will perform all the no cuttof nonbond interactions
without the use of nonbond list, the same as on GRAPE machine.

  

File: Nbonds, Node: Tables, Up: Top, Previous: Function, Next: Cubes

There are two independent implementations of table lookups. LOOKup uses fast
lookups into internally generated tables (using standard CHARMM FAST energy
routines), and ETABLE reads in tables from external files.

LOOKup
Fast non-bonded force and energy calculation for standard MD simulations,
in particular simulations in explicit solvent (water). Not working with
(or not tested with) BLOCK, TSM, PERT, the various implicit solvent models.
Speedup depends on the size of the system and the number of water molecules
but typically a 2-fold speedup may be obtained compared to the standard fast
expanded routines. 

Reference: Nilsson, L. "Fast lookup tables for pairwise interactions" (2007),
in preparation.

Non-bonded interactions within a user specifed selection of solvent molecules 
(any three site water model with the O first, and followed by two identical
hydrogens in the PSF should work) are removed from the regular non-bonded
lists and are instead handled by a dedicated routine. Speedup is achieved
through the use of a table lookup of energies and forces. It is easy to extend
to other solvent models, although a table lookup may be inefficient if there
are more than a few atomtypes in the model. 

Similar tables are used for interactions between the selected solvent molecules
and the rest of the system ("solvent-solute") and for the solute-solute
interactions, execpt 1-4 interactions which are sent on to the standard
routine.

The code works with PBC (images/crystal) and in parallel, with all cutoff
methods implemented in the ENBAEXP routine (eg, FSHIFT,FSWITCH,VSHIFT and the
real space part of PME).

Works with nonbond-list methods that generate an atom based non-bond list
(BYGRoup, BYCBim).

No second derivatives.

Use:

LOOKup { RESEt                                                           }
       { atom-selection [[NO]INTerpolate] [TABIncr <int>] [[NO]ENERgy] -
         [NOVU] [NOUU]                                                   }  

        
RESEt         The regular routines will be used

INTerpolate   Linear interpolation will be used in the table lookup
NOINterpolate

TABIncrement  Determines resolution and size (=TABINCREMENT*CTOFNB**2) of
              lookup tables, which are indexed using Rij**2. For instance
              the energy of two water oxygens at a distance of R is found
              in EOO(R*R*TABINCREMENT).

NOENergy      Energies will only be evaluated when non-bond list is updated
              This may result in apparent ECHEck violations; if so you can
              increase ECHEck (or turn off the check: ECHEck -1.0)
ENERgy        Energies will always be evaluated

NOVU          Do not use lookup fopr the solVent-solUte interactions 
NOUU          Do not use lookup fopr the solUte-solUte interactions.
              NB! 1-4 interactions are always handled by standard routine 

atom-selection  Select the waters to be included. This is mandatory.

DEFAULTS:   
  NOENergy; INTErpolation; TABIncrement 20; no selection; use solvent-solute 
  and solute-solute lookup routines.

Increasing TABI, using the ENERGY and INTERPOLATE options will increase
the size of the lookup tables, which may reduce the speed of the calculation
if it leads to cache misses - experimentation is adviced.
With INTERpolate you can use a smaller TABI and get good accuracy.

Everything is reset on each invocation.

Usage example (see also test/c34test/lookup.inp):
! When the system is completely defined (PSF,cooordinateds, IMAGES/CRYSTAL,,,)
! first a call to energy to fill interaction coefficient arrays
ENERGY FSHIFT VSHIFT CDIE CUTNB 14.0 CTOFNB 12.0 
LOOKUP  SELE SEGI WAT END

ENERGY / MINIMIZE / DYNAMICS   - but do NOT change the cutoff distances 
                                 or options!
 
Implementation:
Code is protected by ##LOOKUP pref.dat keyword
LOOKUP     Identifies the water molecules to get special treatment 
           (or turns the whole thing off). 
           Fills the lookup tables  using calls to ENBAEXP
           - so it is important that an energy
           evaluation has already been performed in order to have
           all coefficient arrays properly filled.
           
           When the non-bond list has been generated in NBONDS, all atom pairs
           belonging to this set are removed from the regular and image lists,
           and placed in their own lists (one for water-water oxygen-oxygen
           pairs, and one each for solvent-solute and solute-solute (non 1-4)
           atom pairs).

           This scheme should work transparently for all methods as
           long as the atom based lists are generated in the first place. There
           is a slight penalty for this second pass through the list, but the
           implementation is very non-intrusive.
           
           Routine ELOOKUP is called from ENERGY (and EIMAGE) to compute the
           water-water and water-solute interactions as specified. 
           Since there is only one table for solvent-solvent energies all the
           solvent-solvent nonbond energy (vdW+Coulomb) computed in EWWNB is
           reported as electrostatic, with zero returned for EvdW; 
           the solvent-solute part separates Coulomb and vdW.
 
           Specifying NOENenergy skips the energy lookups, except at nonbond 
           list updates. This can give  a modest (5-10%) speedup, at the risk
           of getting caught by ECHECK; ECHECK can be turned off by
           ECHECK -1.0.
 
           The tables are stored in single precision and some intermediate
           variables are also in single precision. For typical systems the
           relative error in total energy  and GRMS is of the order of 10^-4.
           Using the INTErpolation option and/or increasing TABIncrement can
           reduce this somewhat, but the error seems to be dominated by the
           single precision noise rather quickly.

 -----------------------------------------------------
        The nonbond energy terms may also be specified with a user supplied
binary lookup table. The command

                READ TABLE UNIT int

will invoke this feature and disable all other energy term options.
The nonbond list specifiers will still be used (cutoff distances...).

This is completely separate from the NBSOlv method outlined above.

        This feature is not designed for casual users, and is not supported
with test cases. Also, in version17, there is an uncorrected bug in
the second derivative determination.

        To use this feature, first read the common file ETABLE.FCM
for a description of variables, and then create a file the the routine
REATBL (consult the source) can read. The sources for this option
are contained in the file ETABLE.FLX.


File: Nbonds, Node: Cubes, Up: Top, Next: Top, Previous: Tables

SUMMARY

     The purpose of using finite cutoffs in energy calculations is to 
reduce, from O(N*N) to O(N), the number of nonbonded interaction 
terms that actually need to be computed.  CHARMM has three ways of 
building the nonbonded interaction list:  BYGRoups, BYCUbes, and 
'By-Clusters-in-Cubes,' or BYCC.  The BYCC method is a combination  
of the earlier BYGR and BYCU methods.  For a given set of atoms
and a given cut-off distance, all three algorithms should generate
the same non-bonded list.

     BYCBim extends the BYCUbes method to systems with images or periodic
boundaries, but is more restricted in other options that are compatible
with it. It will not work with replica, extended electrostatics, or constant
pressure dynamics since no group list is generated. However, BYCBim does work
in parallel mode; BYCUbes does not.

     The differences between BYCUbes, BYGRoups, and BYCC are in speed and
memory requirements. This section gives a description of the algorithms, 
followed by a description of some of the unique aspects of BYCC.  Finally, 
a few general guidelines for choosing between the earlier BYCUbes and
BYGRoups methods are presented.

ALGORITHMS

   Atom-based calculations:

	The basis of the efficiency of the BYGRoups algorithm over a 
brute-force comparison of each atomic position with all others in the 
system is that it clusters atoms into chemical groups, initially ignoring 
the individual atoms, themselves.  This significantly reduces the 
number of pairs of particles that need to be examined.  Effectively, BYGR
speeds up the calculation by reducing the particle density, which it
does by simply redefining the particles. Once a list of group-group
pairs satisfying the initial distance criterion is made, only atoms 
from this relatively short list are then considered for further atom-
atom distance testing.

	In contrast, the efficiency of the BYCUbes algorithm is 
based on the partitioning of the system into small cubical regions.
A synopsis of how it works:

(1) Find a rectangular parallelepiped that bounds the system (with
    margins), is aligned with the Cartesian axes, and has sides that
    are integral multiples of the cutoff distance.  Divide this
    parallelepiped, or box, into cubes whose sides are equal in
    length to the cutoff distance.

(2) Identify the cube that contains each atom.

(3) For each cube C, loop over each atom A contained in C.  For each
    such A, loop over each atom A' contained in the 27-cube surrounding 
    region, which is the (3 cube x 3 cube x 3 cube) region that contains 
    the central cube. If the A--A' pair falls within the cutoff distance,  
    check various other inclusion criteria; e.g. that the pair is not a
    1-2 or 1-3 excluded pair.

Hence, BYCU speeds up the generation of the non-bonded list by distance-
testing only pairs of atoms that are in nearby regions.

	The 'By-Clusters-in-Cubes', or BYCC, algorithm 
incorporates both the partitioning technique of BYCUbes and also
the atomic grouping technique of BYGRoups.  BYCC divides the system
into a cubical grid and compares particles only in adjacent cubes,
as BYCU does.  However, the particles it compares are clusters of
atoms, in the spirit of BYGR, as opposed to individual atoms.
The use of both techniques in BYCC allows for greater efficiency
than is possible using either technique alone, and for this reason
BYCC is generally faster than the other algorithms.  In addition, 
because it does the final atom-atom distance calculations, 
the exclusions, and the formation of the non-bonded list all in one 
final loop, the routine needs to store only a cluster-cluster pair-
list as a work array.  Thus, the memory requirements are reduced 
relative to BYCU, which stores a much longer atom-atom pairlist
(essentially a second non-bonded list) internally.

	Like BYCU, the computational time for BYCC increases 
linearly with the number of atoms and sigmoidally as a function of
cut-off distance.  When (cut-off distance) << (radius of system) the 
time dependence on cut-off distance is essentially third order, but
depends on the system, the machine type, and the actual cutoff value.
The computational time for BYGR increases quadratically with the size 
of the system. NB: The energy calculations are always O(N) if a non-
bonded list is being used, regardless of which method is used to 
generate the list.  However, generation of the list itself can be
either O(N*N) (e.g. BYGR) or O(N) (e.g. BYCU or BYCC).

   Group-based calculations: 

	The BYCC option, like the others, supports calculations 
based on chemical groups instead of atoms, and an additional option 
for extended electrostatics.  However, clusters are not used in 
group-based calculations, since their role is in effect subsumed by 
the groups.  The relative speed advantage of BYCC over BYCU is 
therefore diminished in group-based computations.  
	For the BYCC or BYCU options, when extended 
electrostatics are requested, CPU time for nonbonded list 
generation and calculation of extended electrostatics will depend 
on the extended electrostatics cut-off distance, CTEXNB. 
Hence, smaller CTEXNB values will significantly speed up the
calculations. For BYGR, CPU time is independent of CTEXNB.
	Note that for group-based calculations with either BYCC, 
BYCU, or BYGR, the groups are treated as point-like.  This means 
that, unlike the case for cluster size in atom-based calculations, 
there is no dependence of CPU time on group size in the generation
of the non-bonded list.
	
BYCC:  CLUSTERS

	With the BYCC option, a requirement for the atom-based 
calculations is that clusters of atoms need to be created.  A cluster is 
a set of atoms that have mutually close connectivity relationships. 
Generally, small, dense, uniformly-sized clusters yield the most 
efficient non-bonded list generation. The use of a separate clustering 
scheme instead of the usual chemical grouping arrangement used 
elsewhere in CHARMM (e.g. in the BYGRoups algorithm) provides a
separate spatial organization of the system that can be manipulated 
and optimized largely independently of the other CHARMM functions.
Because the criteria for the grouping of atoms are chemical, 
whereas those for the clustering of atoms are spatial, the optimal
arrangement of atoms for CHARMM groups is not necessarily (and, 
indeed, not usually) the optimal arrangement for clusters.
       Clusters are generated by default just before the non-bonded list
is updated for the first time with BYCC invoked.  Thereafter, they are
regenerated when a change in the connectivity or the topology of the system
is detected.

	Clusters can also be (re)made with the command:

	MKCLuster  [CADL] [CMARgin] [EXCL] 

        This command is no longer required for the use of BYCC
(since its equivalent is called by default), but can be used in the
middle of a CHARMM run (for example to change the CADL or CMAR
parameters) or to ensure reformation of atom clusters.
 1) Keyword CADL 
    Consecutive atom distance limit -- distance between consecutive
     atoms (as defined by order in RTF/PSF) in a cluster beyond which
     the cluster will be split.  Helps prevent the generation of overly large
     clusters. Default : 4.0 Angstroms.
 2) Keyword CMAR
    Limiting cluster margin width.  The cluster margin is the calculated
     distance that is added to the non-bonded cutoff distance,
     when partitioning the clusters in the system, to ensure that all      
     atom-atom pairs are in fact included in the distance testing.  This 
     margin depends on the largest cluster size and will affect CPU time  
     (larger margin means slower non-bonded list generation). CMAR sets 
     the upper limit for this added margin. Default: 5.0 Angstroms.
     Note that while the calculated cluster margin is the width that 
     GUARANTEES all-atom testing, in practice CMAR can often be 
     set below the calculated value (speeding up the calculations) without   
     changing the results. This is particularly true when there is a large 
     spread in the cluster sizes or when there are statistical outliers.
     If CMAR < calculated cluster margin, a warning is issued.
 3) Keyword EXCL 
    The MKCL command will result in a re-generation of the exclusion
     table (from the current exclusion list) if the "EXCLusions" 
     keyword is specified.  This may be necessary
     if the connectivity of the system is altered during a CHARMM run
     or if the exclusion list otherwise changes.

BYCC: ACTIVE ATOMS

An additional feature of the BYCC algorithm is its ability to handle 
"active" atom selections, which are specified with the command 
	
	NBACtive [atom selection]

	The purpose of this feature is to allow the user to focus 
calculations on regions of interest without having to alter the psf or 
coordinate files. It allows the non-bonded list generation routine 
(BYCC) to ignore completely atoms that are not defined as active, so 
that only active atoms appear in the non-bonded list. This speeds up 
non-bonded updates and energy calculations, and, in addition, it 
allows for more selectivity in energy calculations-- in dynamics 
and minimization, particularly.
	An example of the ideal use of this command is the study 
of a single subunit in a protein containing several.  Another is the 
study of a single side chain, or a group of side chains, on a fixed 
protein backbone.
        If the energy calculations on the active portion of the 
system are to be "consistent" with the inactive portion of the 
system, a buffer region (ideally of width CUTNB) is required 
surrounding the region of interest. This buffer region should 
be active and fixed.
	Active clusters and groups are defined automatically on 
the basis of the active atom selection.  If no active atom selection
is given, it is assumed that all atoms are active.
	In its current implementation, the active atom selection does 
not affect the bonded energies (note to developers:  this should be 
rectified in the future).  This has at least two implications:  1) It
is currently recommended that in conjunction with this command, 
'inactive' atoms be fixed with the CONS FIX command, since otherwise
they will contribute to the bonded terms but not to the non-bonded
terms.  2) The connectivity between inactive and active regions, if
it exists in the original system, is not broken by the NBACtive 
command.  This means that regions defined as active must remain
'tethered' to the inactive regions.


BYCC: EXAMPLE OF ADDITIONAL COMMANDS

The additional commands in a typical CHARMM script using the 
BYCC option for non-bonded list generation would be as follows:
The BYCC keyword in the NBOND command is necessary for using
the BYCC option.

The other commands are optional.

! make clusters
 MKCL           !no longer necessary 

! define active atom set and buffer subset (optional)
DEFINE ACTIVE SELE (active atom selection) END  
DEFINE BUFFER SELE (buffer atom selection) END  !subset of ACTIVE

! "activate" the selected atoms (optional)
 NBACtive SELE ACTIVE END   !default is all atoms

! fix inactive atoms (optional)
 CONS FIX SELE ((.NOT.ACTIVE) .OR. BUFFER) END

! call non-bonded list generation routine
 NBONd -
  BYCC           !necessary 


PAYOFF THRESHOLD:  Comparison of BYCUbes and BYGroups
        
	Given any pair of functions F1(N) and F2(N), which are O(N)
and O(N*N) respectively, there will always be some constant N0 such
that F2(N) > F1(N) for all N > N0.  This constant N0 may be referred
to as the "payoff threshold," since it is the system size above 
which the BYCUbes algorithm will be faster than BYGRoups for a 
given cut-off distance, particle density, machine type, and set of
options.

    The following are some properties of the payoff threshold that 
    can be used as guidelines for choosing between BYGRoups and 
    BYCUbes:

(1) The payoff threshold is smaller on a vector machine than on a
    scalar machine.  That is, BYCUbes is more vectorizable than
    BYGRoups.

(2) The shape of the system does not have a big impact on the payoff
    threshold.  BYCUbes operates by drawing a rectangular box around the
    system and dividing it into cubes, but on a serial machine the 
    empty cubes take little time.

(3) The payoff threshold usually grows steeply with cutoff distance.
    This is because the speed advantage of the BYCUbes algorithm 
    is based upon the compartmentalization of the system.  The smaller
    the compartments (cubes), the more advantage.  Since the compartment
    size is based on the cutoff distance, as the cutoff distance 
    increases, there is less advantage for BYCUbes over BYGRoups.

(4) BYCUbes trades memory for time.  Its memory requirements are 
    significantly higher than those of BYGR and they may be 
    prohibitive for large systems.

