.. py:module::valbond

=============
VALBOND-TRANS
=============

by  Ivan Tubert-Brohman, Maurus Schmid
and Markus Meuwly (m.meuwly@unibas.ch)

(Based on work by Landis and co-workers; see references.)

This is an implementation of a modified VALBOND force field in CHARMM.
The force field supports hypervalent molecules (such as SF6) and
transition metal complexes. For more information about VALBOND, see
the references at the end.


.. _valbond_syntax:

Syntax for the VALBOND-TRANS commands
-------------------------------------

VALBOND COMMANDS

::

    VALB DONE
    VALB E <atom-name> <e>
    VALB LP <atom-name> <lp-count>
    VALB HYBR <atom-name-1> <atom-name-2> <p> <d> <ishyp>
    VALB HYBA <atom-name> <p> <d> <ishyp>
    VALB PARA <param-name> <z1> <z2> <value>
    VALB SKIP <atom-selection>
    VALB INCL <atom-selection>
    VALB PRINT

DONE

::

    VALB DONE

Calculates all the hybridizations and other terms needed for computing
the VALBOND energy. It must be called AFTER the molecule is generated
but before energies (optimizations, dynamics, etc.) are required.

----

::

    VALB E <atom-name> <e>

Sets the formal number of metal d-electrons. This is the group
number minus the oxidation number. For example, for Ir3+, E = 6.

This will be used for assigning the hybridization of the metal
automatically.  However, it can be overridden manually using HYBR or
HYBA (see below).

----

::

    VALB LP <atom-name> <n>

Set the number of lone pairs for an atom manually. For example,

----

::

    VALB LP OH2 2

will assign two lone pairs to atoms named OH2. This is the "iupac"
atom name; note that this will affect atoms with that name in all
residues.  LP must be called before DONE.

----

::

    VALB HYBR <atom-name-1> <atom-name-2> <p> <d> <ishyp>

Set the hybridization for an orbital manually.

This is similar to the LP command, but requires two atoms, to set
the hybridization of the bond between atom1 and atom2. Note that
this is not commutative; the hybridization between atom2 and atom1
can be different.  p and d are real numbers; for example it is
possible to have sp3.591.

The <ishyp> term should be non-zero if atom1 is hypervalent.

----

::

    VALB HYBA <atom-name> <p> <d> <ishyp>

Similar to HYBR, but sets all the hybridizations for the bonds
centered on a given atom.

----

::

    VALB PARA <param-name> <z1> <z2> <value>

Override the default value of a VALBOND parameter. Available names are:

One-element parameters:

::

    KHV - hypervalent force constant
    LP  - lone pair weight of element
    U2R - UFF nonbond radius
    U2S - UFF nonbond scaling
    U2E - UFF nonbond energy
    VAL - valence
    EN  - electronegativity
    BLI - bond lengthening intensity
    BLS - bond lenghtening sensitivity

Two-element parameters:

::

    K - force constant for element I bonded to J
    WT- weight of element I bonded to J
    TR- trans energy offset for bond I trans to bond to J (symmetric)

For parameters depending on only one element I, the value of z2 is
ignored. Note that z1 and z2 are atomic numbers (that is H = 1, C =
6).

----

::

    SKIP <atom-selection>

Skips the atoms in <atom-selection> for the calculation in valbond.
For example:

::

    VALB SKIP SELE RESN TIP3 END

----

::

    INCL <atom-selection>

Includes the atoms in <atom-selection> for the calculation in valbond.
Forces override, use with caution to not calculate energies
in valbond and classical CHARMM.
For example:

::

    VALB INCL SELE RESN LIG END


----

::

    VALB PRINT

Produce a verbose output the next time the energy is evaluated. By
default, this verbose output is not printed out. Using this
command, one can force it to be printed.


.. _valbond_rtf:

Notes about the RTF
-------------------

The RTF defining the residues that will be treated with valbond should not
include AUTO ANGLES. That is,

::

    AUTO DIHE

is OK, but

::

    AUTO ANGLES DIHE

is not. This will ensure that the bending energies are not computed
using the standard CHARMM harmonic bending terms.

It is possible to have a residue that mixes VALBOND angles with CHARMM
angles.  This is achieved by adding the CHARMM angles manually to the
RTF. VALBOND will skip all angles centered on an atom that is in the
middle of any CHARMM angle.  For example, consider the following RTF
definition for methanol:

::

    ATOM C1
    ATOM C2
    ATOM O3
    ATOM H4
    ATOM H5
    ATOM H6
    ATOM H7
    BOND C1 C2
    BOND C2 O3
    BOND C1 H4
    BOND C1 H5
    BOND C1 H6
    BOND O3 H7
    ANGLE H4 C1 H5
    ANGLE H4 C1 H6
    ANGLE H6 C1 H5
    ANGLE H4 C1 C2
    ANGLE H5 C1 C2
    ANGLE H6 C1 C2

Here all the angles on the carbon are defined explicitly, but no angle
on the oxygen is defined. Therefore, the C-O-H angle will be treated by
VALBOND. Using the SKIP command it is possible to skip atoms manually.
where some angles use VALBOND and some angles use CHARMM. For example,
on a tetrahedral atom, which has six angles, all six angles use CHARMM
or all six angles use VALBOND.

.. _valbond_hybridization:

Hybridization for transition metals
-----------------------------------

According to the VALBOND model, transition metal complexes have an sd^n
hybridization, with no p orbital participation.

* For non-hypervalent compounds, n = N - 1
* For hypervalent compounds, n = N - 1 - H

where N is the number of ligands, and H is the number of 3c4e bonds
H = (E - 12)/2, where E is the electron count including the ligands.

Therefore, for a hypervalent complex n = ((12 - e)/2) - 1, where e is
the formal d-electron count for the metal itself. For example, Ir3+ with
six 2-e ligands, e = 6, N = 6, E = 18, n = 2, H = 3.


.. _valbond_caveats:

Caveats
-------

* Second derivatives are NOT implemented. Therefore some optimization
  methods and frequency calculations are not likely to work.
* The method has been tested with geometry optimizations and simple
  molecular dynamics. Compatibility with FEP and other modules should
  be carefully checked.
* Valbond is compatible with parallel execution, but will use only one CPU.


.. _valbond_references:

References for VALBOND-TRANS
----------------------------

0. Tubert-Brohman, I.; Schmid, M.; Meuwly, M. A molecular mechanics force field
   for octahedral organometallic compounds with inclusion of the trans influence.
   J. Chem. Theory Comput. 2009, 5, 530-539.

1. Root, D. M.; Landis, C. R.; Cleveland, T. Valence Bond Concepts Applied to
   the Molecular Mechanics Description of Molecular Shapes. 1. Application to
   Nonhypervalent Molecules of the P-Block. J. Am. Chem. Soc. 1993, 115,
   4201-4209.

2. Cleveland, T.; Landis, C. R. Valence Bond Concepts Applied to the
   Molecular Mechanics Description of Molecular Shapes. 2. Application to
   Hypervalent Molecules of the P-Block. J. Am. Chem. Soc. 1996, 118, 6020-6030.
   doi:10.1021/ja9506521

3. Landis, C. R.; Cleveland, T.; Firman; T. K. Valence Bond Concepts Applied
   to the Molecular Mechanics Description of Molecular Shapes. 3. Application to
   Transition Metal Alkyls and Hydrides. J. Am. Chem. Soc. 1998, 120, 2641-2649.
   doi:10.1021/ja9734859

