# Polar Phase Field model

With the help of the [PolarPhaseFieldAction](/PolarPhaseFieldAction.md) users can
set up the three phase polar phase field model from [!cite](Momeni2014).

## Ginzburg-Landau equations

The model consists of the following main equations (eqs. (14) and (15) in the 2014 paper)

\begin{equation}
\underbrace{\frac1{L_\Upsilon}\frac{\partial\Upsilon}{\partial t} }_{\text{SusceptibilityTimeDerivative}} =
\underbrace{-\frac{\partial\psi^l}{\partial\Upsilon}}_{K_1}\,
\underbrace{-\frac{\beta^{21}}2\frac{\partial\phi(\Upsilon,a_\phi,a_0)}{\partial\Upsilon}|\nabla\vartheta|^2}_{K_2}
+ \underbrace{\nabla\cdot\left[\beta^{SO}(\vartheta)\nabla\Upsilon\right]}_{\text{MatDiffusion}}
\end{equation}

\begin{equation}
\underbrace{\frac1{L_\vartheta}\frac{\partial\vartheta}{\partial t} }_{\text{SusceptibilityTimeDerivative}} =
\underbrace{-\frac{\partial\psi^l}{\partial\vartheta}}_{K_1}\,
\underbrace{-\frac12\frac{\partial\beta^{SO}(\vartheta)}{\partial\vartheta}|\nabla\Upsilon|^2}_{K_2}
+ \underbrace{\nabla\cdot\left[\beta^{21}\phi(\Upsilon,a_\phi,a_0)\nabla\vartheta\right]}_{\text{MatDiffusion}}
\end{equation}

This model requires the implementation of two new kernels $K_1$ and $K_2$ for
the terms indicated. The remaining terms can be modeled using the existing
[SusceptibilityTimeDerivative](/SusceptibilityTimeDerivative.md) and
[MatDiffusion](/MatDiffusion.md) kernels.

!bibtex bibliography
