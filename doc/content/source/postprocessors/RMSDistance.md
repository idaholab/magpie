# RMSDistance

!syntax description /Postprocessors/RMSDistance

The RMSDistance postprocessor computes the RMS distance of the variable parameter $v$ from a fixed point $\vec{r}_0$:
\begin{equation}
  \text{RMS} = \sqrt{\frac{\int\limits_{V} \| \vec{r} - \vec{r}_0\|^2 v(\vec{r})dV}{\int\limits_{V} v(\vec{r})dV}}.
\end{equation}

!syntax parameters /Postprocessors/RMSDistance

!syntax inputs /Postprocessors/RMSDistance

!syntax children /Postprocessors/RMSDistance

!bibtex bibliography
