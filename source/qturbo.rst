.. py:module:: qturbo

==========================================================================================
A Quantum Mechanical / Molecular Mechanical (QM/MM) Interface Between TURMOBOLE and CHARMM
==========================================================================================

Christopher N. Rowley (cnrowley@mun.ca)
based on the Q-Chem interface from H. Lee Woodcock (hlw@mail.usf.edu)
which was
based on the GAMESS(US) interface from Milan Hodoscek
(milan@par10.mgsl.dcrt.nih.gov,milan@kihp6.ki.si)
and
the GAMESS(UK) interface from Paul Sherwood
(p.sherwood@dl.ac.uk)

A QM/MM interface between CHARMM and the ab initio/DFT code
TURBOLE is implemented through this module.


.. _qturbo_description:

Description
-----------

The TURBOMOLE QM potential is initialized with the QTURBO command.

::

  qturbo    [REMOve]

  REMOve:  The force field terms involving the atoms designated as QM are removed.
           This is needed for any QM/(MM) calculation as it it how atoms are
           designated as QM atoms.

           QTURBO REMOve SELEct RESId 1 END
           ENERgy

.. _qturbo_usage:

Usage
-----

This interface follows the same conventions as the CHARMM/GAMESS and
CHARMM/Q-Chem interfaces, where a standard CHARMM MM topology and coordinates
configuration is prepared and then some atoms are designated to be modeled using
QM. The location of the TURBOMOLE input files (qturboinpath), the directory where
CHARMM and TURBOMOLE should write files for communication between them
(qturbooutpath), and the Python wrapper script  that performs the execution of
TURBOMOLE (qturboexe) are specified through ENVI commands. Alternatively, these
could be set in the shell environment that CHARMM is executed in

::

   envi qturboinpath "data/qturbo/"
   envi qturbooutpath "cquantumtest/turbodir/"
   envi qturboexe "../support/programs/turbo.py"

For every QM/MM energy/gradient evaulation, CHARMM writes a coord file
containing the coordinates of the QM atoms. The TURBOMOLE execution script
is then run by CHARMM, which copies this file to a temporary directory
along with the neccessary TURBOMOLE input files from qturboinpath (control,
mos, mos.bin, basis, auxbasis...). If there are MM atoms present, the location
and magnitude of these charges will be written to the file pc_charges.
Depending on whether the jobtype is a DFT calculation or an ab initio method,
ridft/rdgrad or dscf/ricc2 are then executed by Python wrapper script,
qturboexe. TURBOMOLE will generate a file containing the forces on the QM
atoms (gradient) and the gradient on the point charges due to the QM region
(pc_gradient). These files are copied back to qturbooutpath by the wrapper
script. The wrapper script then exits and CHARMM will read in the gradient
files and continue the calculation.

The files in qturboinpath must be configured by the user using the TURBOMOLE
define program before the execution of CHARMM. See the TURBOMOLE manual
for instructions on how to do this.

.. _qturbo_installation:

To use the CHARMM TURBOMOLE interface, CHARMM must be compiled with
the QTURBO keyword, which is invoked by the command:

::

    install.com <machine-type> <CHARMM size> QT <other CHARMM options>
