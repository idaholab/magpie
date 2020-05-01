# SpectralExecutionerBase

## Spectral executioners

In general, executioners determine the type of solve for finite element problems. In spectral solvers, the executioner is also responsible for **actually** solving the system. `SpectralExecutionerBase` is an base class which holds a pointer to an `FFTProblem` and initializes an `FEProblem`.  

The executioner has helper methods such as `kVectorMultiply` which may be called, e.g., to obtain spatial derivatives in the reciprocal space by child classes.

The spectral solve application developer will need to inherit from this class and provide appropriate explicit or implicit algorithm to obtain each step's solution.

## ::execute for computing spatial derivatives

The `execute` method is pure virtual and handles the numerical time stepping, the forward and inverse discrete Fourier transforms through `FFTWBufferBase`.   


!syntax parameters /Executioner/SpectralExecutionerBase

!syntax children /Executioner/SpectralExecutionerBase

!bibtex bibliography
