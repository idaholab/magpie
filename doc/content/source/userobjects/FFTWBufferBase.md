# FFTWBufferBase

## Buffer for FFTW spectral executioners

`FFTWBufferBase` is a templated buffer that contains specific interleaved data buffer. Its implementation is linked to the external library FFTW. On construction, `FFTWBufferBase` creates `fftw_plan`s for forward and inverse discrete Fourier transform (DFTs). The methods `forward` and `backward` are equipped with profiling information and call FFTW to perform DFT. `FFTWBufferBase` takes care of deleting FFTW pointers by calling `fftw_destroy_plan` upon destruction. 

Use of `FFTWBufferBase`s can be observed in `SpectralExecutioners`.

!syntax parameters /UserObjects/FFTWBufferBase

!syntax inputs /UserObjects/FFTWBufferBase

!syntax children /UserObjects/FFTWBufferBase

!bibtex bibliography
