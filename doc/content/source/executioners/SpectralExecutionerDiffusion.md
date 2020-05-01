# SpectralExecutionerDiffusion

!syntax description /Executioner/SpectralExecutionerDiffusion

## Diffusion spectral solve

`SpectralExecutionerDiffusion` inherits from `SpectralExecutionerBase` and implements problem-specific methods. For example, it computes the diffusion equation's Green's function for a time step `_dt` through the method `getGreensFunction`. This function is used in the `execute` method to explicitly obtain the diffusion solution step by step. The Green's function is transformed to the reciprocal space where it is point-wise multiplied by the diffused variable in the Fourier. The transformation of its result into the real space carries the proper grid-related scale factor through the `FFTBufferBase` method `backward`.

As shown below, this diffusion execution does not need to iterate due to its linear nature. Spectral executioners are, however, not restricted to linear solves.

!listing src/executioners/SpectralExecutionerDiffusion.C
         re=void\sSpectralExecutionerDiffusion::execute.*}

!syntax parameters /Executioner/SpectralExecutionerDiffusion

!syntax inputs /Executioner/SpectralExecutionerDiffusion

!syntax children /Executioner/SpectralExecutionerDiffusion

!bibtex bibliography
