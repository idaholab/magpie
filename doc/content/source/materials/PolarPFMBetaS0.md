# PolarPFMPhi

This material implements the gradient energy coefficient $\beta^{S0}(\vartheta)$
from [!cite](Momeni2014).

This object can be used with the [MatDiffusion](/MatDiffusion.md) and
[PolarPFMDerivative](/PolarPFMDerivative.md) kernels to construct the highlighted
terms in eqs. (14) and (15)

\begin{equation}
\gray{\frac1{L_\Upsilon}\frac{\partial\Upsilon}{\partial t} =
-\frac{\partial\psi^l}{\partial\Upsilon}
-\frac{\beta^{21}}2\frac{\partial\phi(\Upsilon,a_\phi,a_0)}{\partial\Upsilon}|\nabla\vartheta|^2}
+\nabla\cdot\left[\red{\beta^{S0}(\vartheta)}\nabla\Upsilon\right]
\end{equation}

\begin{equation}
\gray{\frac1{L_\vartheta}\frac{\partial\vartheta}{\partial t} =
-\frac{\partial\psi^l}{\partial\vartheta}}
-\frac12\frac{\partial\red{\beta^{S0}(\vartheta)}}{\partial\vartheta}|\nabla\Upsilon|^2
\gray{+\nabla\cdot\left[\beta^{21}\phi(\Upsilon,a_\phi,a_0)\nabla\vartheta\right]}
\end{equation}

!syntax parameters /Materials/PolarPFMPhi

!syntax inputs /Materials/PolarPFMPhi

!syntax children /Materials/PolarPFMPhi

!bibtex bibliography
