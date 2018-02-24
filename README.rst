================
LAMMPS_TO_CUBE
================

Description
-------------------

Converts a LAMMPS 3d histogram to a cube file for visualization in VMD. The current way I compute histograms in LAMMPS is with the following
commands in the input file:

:: 
 
    compute         cc1 TIP4P_Water chunk/atom bin/3d x 0.0 0.01 y 0.0 0.01 z 0.0 0.01 units reduced
    fix             7 TIP4P_Water ave/chunk ${nevery} ${nrepeat} ${nfreq} cc1 density/number norm all file TIP4P.profile ave running overwrite

To understand these commands, please consult the LAMMPS website. Anything defined with the dollar sign followed by curly braces is a variable that can be 
defined earlier in the lammps input script.

**note: the above command can change between different LAMMPS versions. This was working with LAMMPS-17Nov16.**


:Authors: Peter Boyd
:Date: 24/2/2018
