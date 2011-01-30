.. py:module:: sbound

=======================================================
Method and implementation of deformable boundary forces
=======================================================

The use of deformable boundary forces is in studying small
localized regions of solvent, say around an active site.  The boundary
forces are applied to the atoms in the solvent and serve to contain
the reaction zone.

Generally the boundary forces are computed from the deformable
boundary method of C. L. Brooks III and M. Karplus, J. Chem. Phys., 79,
6312(1983).  Following generation one must;

i) Generate the corresponding boundary potential

ii) Read the tabulated boundary potential into CHARMM

iii) Set-up the mapping CHARMM uses to connect table entries with
     boundary constrained atoms

Steps ii) and iii) ** MUST ** be done everytime the boundary forces are to be
used.  For example, during the initial stages of a dynamics simulation and
at all subsequent restarts.

The syntax for generating the potential and reading and setting up
the table structure is given in the following mode.

.. _sbound_syntax:

Syntax
------

::

   SBOUnd POTEntial INPUt <integer> OUTPut <integer>

Integrates forces to get potential and generates cubic spline
approximation of the potential.

::

   SBOUnd SET XREF <real> YREF <real> ZREF <real>
          ASSIgn <table number> <selection-syntax>

Solvent boundary routine to set boundary geometry
and specify the atoms referring to the tables.  Note the
assign must be repeated for each table entry.

::

   SBOUnd READ UNIT <integer>

This routine reads in the boundary potential splines
for the boundary constraints.
