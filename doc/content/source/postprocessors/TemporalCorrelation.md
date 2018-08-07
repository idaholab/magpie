# TemporalCorrelation

!syntax description /Postprocessors/TemporalCorrelation

\begin{equation}
  \Delta U = \sqrt{\frac1\Omega\int_\Omega |\dot u|^2 dr},
\end{equation}

where $\Omega$ is the simulation cell volume, and $\Delta U$ can be taken as a
measure of the correlation between timesteps. Small values indicate slow,
gradual change, large values indicate strong fluctuations.

!syntax parameters /Postprocessors/TemporalCorrelation

!syntax inputs /Postprocessors/TemporalCorrelation

!syntax children /Postprocessors/TemporalCorrelation

!bibtex bibliography
