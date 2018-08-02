# FourierTransform

!syntax description /UserObjects/FourierTransform

This user object computes a Fourier transform of a selected variable field using
the FFT library FFTW3 (needs to be installed prior to building Magpie).
Internally a _Real to Real_ transform is computed. The fransform is computed as
a discrete Fourier transform (DFT)_ on a regular orthogonal grid. This grid is
either specified using the `grid` parameter or can be auto-detected when a
[MyTRIMMesh](/MyTRIMMesh.md) is used in the simulation.

The fourier transform data is used by the following objects:

- [FourierLengthScale](/FourierLengthScale.md) Postprocessor - Extracts a scalar
  average lenghth scale from the transform

- [FourierPowerSpectrum](/FourierPowerSpectrum.md) VectorPostprocessor -
  Extracts a radially (in frequency space) averaged power spectrum from the transform


!syntax parameters /UserObjects/FourierTransform

## Installing FFTW3

### Ubuntu

```
sudo apt install libfftw3-dev pkg-config
```

### macOS

Using homebrew do

```
brew install fftw
```

!syntax inputs /UserObjects/FourierTransform

!syntax children /UserObjects/FourierTransform

!bibtex bibliography
