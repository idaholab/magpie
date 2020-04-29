# FFTBufferBase


## Buffer for spectral executioners

`FFTBufferBase` is an `ElementUserObject` that holds information on domain variables for FFT solves. 
Buffers receive their initial value from Moose variables. Domain size, corner coordinates, dimensions, 
and  k vector components can be accessed through methods such as `getBoxSize`, `getMinCorner`, 
`dim` and `kTable`, respectively. `kTable` provides a useful factor $\frac{2 \pi i}{L}k$ to 
perform computations in the reciprocal space. This is used, for example, in `SpectralExecutionerBase`
 to obtain spatial derivatives. `FFTBufferBase` holds pointers to the start of the buffer data, 
which has real and complex counterparts. Real space and reciprocal space data are member variables 
in the form of templated `FFTData`.  

Real space data can be accessed by location, i.e. `Point`, through parenthesis operator methods. 
`forwardRaw` and `backwardRaw` are pure virtual methods which must be implemented in derived classes.
 `FFTWBufferBase` is an example of such classes, where direct and inverse discrete Fourier 
transforms are done via the external library FFTW (Fastest Fourier Transform in the West).

!bibtex bibliography
