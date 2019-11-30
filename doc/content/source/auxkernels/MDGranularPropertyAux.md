# MDGranularPropertyAux

!syntax description /AuxKernels/MDGranularPropertyAux

Injects properties collected for MD particles from MDRunBase object user_object into auxiliary variable.
The properties are obtained from an `MDRunBase` object. The list of available properties is:
velocity in $x, y, z$ `vel_x`, `vel_y`, `vel_z`; force in $x, y, z$ `force_x`, `force_y`, `force_z`;
charge, and radius.

!syntax parameters /AuxKernels/MDGranularPropertyAux

!syntax inputs /AuxKernels/MDGranularPropertyAux

!syntax children /AuxKernels/MDGranularPropertyAux

!bibtex bibliography
