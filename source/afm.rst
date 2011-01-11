.. py:module:: afm
.. include:: <isonum.txt>

========================
The AFM Module of CHARMM
========================

By Emanuele Paci, 1997/2000


AFM is an external perturbation designed to pull macromolecules
mimicking single molecule experiments (AFM or LOT).
There are three possible way of simulating the pulling; all consist in
applying a suitable force to two atoms.  The force is identical in
magnitude for the two atoms, parallel to the two atoms and directed in
the direction of increasing distance.  The difference in the forces
applied concerns their dependence on the time. The three methods
currently implemented are:

1) constant force.
2) steered molecular dynamics
3) biased molecular dynamics (see HQBM.DOC)

References:

(1) Paci & Karplus, PNAS, 97, 6521-6526, (2000),
(2) Paci et al., JMB, 314, 589-605, (2001).
 

.. _afm_syntax:

Syntax
------

::

   [INPUT AFM command]

     AFM  METHOD ALPHA real [BETA real] two-atom-selection -
                [IUNJ integer] [XIMAX real]

     AFM RESEt


.. _afm_function:
====== =====================================================================
METHOD one of CF, BMD, SMD
ALPHA  force constant (in pN (CF) or pN/Å (SMD/BMD)) 
BETA   cantilever speed (in |Aring|/ps only for SMD)
IUNJ   write the output 

       ::
       
         istep rc[Å] rc_ref[Å] force[pN] energy[kcal/mol]
      
        on unit IUNJ
XIMAX  stops the run when the distance between the two selected atoms
       exceeds XIMAX
RESEt  resets the logicals and turns off the afm module.
====== =====================================================================