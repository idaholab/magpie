# FFTProblem

`FFTProblem` is an augmented `FEProblem` that supports FFT buffers as variables. It skips the nonlinear system check of base classes as the system is solved in the `SpectralExecutioner`. 


!syntax parameters /Problem/FFTProblem

!syntax inputs /Problem/FFTProblem

!syntax children /Problem/FFTProblem

!bibtex bibliography
