# FourierLengthScale

!syntax description /Postprocessors/FourierLengthScale

This Postprocessor acceses a [FourierTransform](/FourierTransform.md) user
object and computes an average frequency weighted by the square of the spectral
intensity. It reports back the reciprocal of this frequence, which corresponds
to an average length scale in length units.

This Postprocessor can for example be used to monitor coarsening in an evolving
micristructure, such as a spinodal decomposition.

!media media/fourierlengthscale.jpg
       caption=Fourier lengthscale plot of an evolving microstructure under spinodal decomposition and coarsening. Insets show the actual microstructure.
       style=width:90%;padding:20px;

!syntax parameters /Postprocessors/FourierLengthScale

!syntax inputs /Postprocessors/FourierLengthScale

!syntax children /Postprocessors/FourierLengthScale

!bibtex bibliography
