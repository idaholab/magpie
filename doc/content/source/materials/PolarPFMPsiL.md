# PolarPFMPsiL

The Bulk Helmholtz energy $\psi^l$ of the [polar phase field model](/PolarPhaseField/index.md)
from [!cite](Momeni2014). $\psi^l$ is a function
of $\Upsilon$ and $\vartheta$.

\begin{equation}
\gray{\frac1{L_\Upsilon}\frac{\partial\Upsilon}{\partial t} =}
-\frac{\partial\red{\psi^l}}{\partial\Upsilon}
\gray{-\frac{\beta^{21}}2\frac{\partial\phi(\Upsilon,a_\phi,a_0)}{\partial\Upsilon}|\nabla\vartheta|^2
+\nabla\cdot\left[\beta^{S0}(\vartheta)\nabla\Upsilon\right]}
\end{equation}

\begin{equation}
\gray{\frac1{L_\vartheta}\frac{\partial\vartheta}{\partial t} =}
-\frac{\partial\red{\psi^l}}{\partial\vartheta}
\gray{-\frac12\frac{\partial\beta^{S0}(\vartheta)}{\partial\vartheta}|\nabla\Upsilon|^2
+\nabla\cdot\left[\beta^{21}\phi(\Upsilon,a_\phi,a_0)\nabla\vartheta\right]}
\end{equation}

!syntax parameters /Materials/PolarPFMPsiL

!syntax inputs /Materials/PolarPFMPsiL

!syntax children /Materials/PolarPFMPsiL

!bibtex bibliography
