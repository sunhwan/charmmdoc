.. py:module:: analysis

=================
Analysis Commands
=================

The ANALysis command is an energy and structure analysis
facility that has been developed to examine both static and dynamic
properties.  The current code allows energy partition analysis and
energy contribution analysis from free energy simulations.  It also
can produce a detailed printout of structural and energy term
contributions for selected atoms


.. _analysis_description:

Description of the ANALysis Command
-----------------------------------

::

   Syntax:

   ANALys { ON                                                         }
          { TERM  { [ALL] }  { NONBond     } [UNIT int] atom-selection }
          {       {  ANY  }  { [NONOnbond] }                           }
          { OFF                                                        }


   ON    Enable energy partition analysis and disable FAST routines.

   OFF   Disable analysis and restore FAST option defaults.

   TERM  Setup energy term print data and disable FAST routines.

   ALL (default)       Print energy terms involving only selected atoms
   ANY                 Print energy terms when any of the atoms is selected

   NONBond             In addition to internal terms, also print nonbond terms
   NONOnbond (default) Do not print electrostatic and vdw energy data

   UNIT integer        Write the energy term printout data to a formatted file
                       Otherwise, write data to the output file.


.. _analysis_energy:

Energy term option of the ANALysis Command
------------------------------------------

The ANALysis ON command enables energy partition analysis and disables the
FAST routines.  This will slow the calculation (especially on vector machines),
but allows a detailed, atom by atom, energy analysis.  Everytime the energy
routine is invoked, the energy for each atom is stored in the ECONT array.
During PERT dynamics, the EPCONT is filled with the time average energy
difference on a atom by atom basis including every step of dynamics.  This
allows the free energy differences to be analyzed based on atom contributions.

The energy partition array can be accessed with the SCALar ECONt
commands.  :ref:`Econt <scalar>`.  The sum of all of the elements
of the ECONT array is usually the total energy, but some energy terms, such
as extended electrostatics, will not be included.  The command:

::

        SCALar ECONT STATistics
        
can be used to check the total energy and the command

::

        SCALar EPCONT ....
        
can be used to examine atom contributions to energy differences for PERT.

The ANALysis TERM command will cause all selected energy terms to be
printed to the specified output unit (default: standard CHARMM output).
The ALL keyword (default) will list elements where all atoms are
selected. The ANY keyword will cause terms including any selected atom
to be listed.  The NONOnbond keyword (default) suppresses the listing
of vdw and electrostatic energy terms.  The NONBonded keywords also
allows the analysis of vdw and electrostatic interactions.

The ANALys OFF command enables the FAST routines and disables the resetting
of the ECONT array (i.e. the ECONT array will not change, but may still be
accessed.  This command also disables the energy term analysis.
