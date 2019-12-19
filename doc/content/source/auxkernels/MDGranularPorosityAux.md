# MDGranularPorosityAux

!syntax description /AuxKernels/MDGranularPorosityAux

Computes porosity $\phi$ or packing fraction $1 - \phi$ and injects it into an aux variable.
The packing fraction is obtained from an `MDRunBase` object. The `MDRunBase` object needs
to operate in `granular` mode which is triggered by providing the `radius` property
in the `md_particle_properties` parameter.

!syntax parameters /AuxKernels/MDGranularPorosityAux

!syntax inputs /AuxKernels/MDGranularPorosityAux

!syntax children /AuxKernels/MDGranularPorosityAux

!bibtex bibliography
