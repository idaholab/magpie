# RadialGreensSource

!syntax description /Kernels/RadialGreensSource

Source term that taked the the difference between the original variable field
and the convolution result and applies it as a source term. The $\gamma$
(`gamma`) rate parameter together with the current timestep size $dt$ is
multiplied as a prefactor onto the data from the user object (`convolution`).

!syntax parameters /Kernels/RadialGreensSource

!syntax inputs /Kernels/RadialGreensSource

!syntax children /Kernels/RadialGreensSource

!bibtex bibliography
