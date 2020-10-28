# Polar Phase Field model

With the help of the [PolarPhaseFieldAction](/PolarPhaseFieldAction.md) users can
set up the three phase polar phase field model from [!cite](Momeni2014).

## Ginzburg-Landau equations

The model consists of the following main equations (eqs. (14) and (15) in the 2014 paper)

\begin{equation}
\underbrace{\frac1{L_\Upsilon}\frac{\partial\Upsilon}{\partial t} }_{\text{SusceptibilityTimeDerivative}} =
\underbrace{-\frac{\partial\psi^l}{\partial\Upsilon}}_{\text{PolarPFMDerivative}}
\underbrace{-\frac{\beta^{21}}2\frac{\partial\phi(\Upsilon,a_\phi,a_0)}{\partial\Upsilon}|\nabla\vartheta|^2}_{\text{PolarPFMGradient}}
+ \underbrace{\nabla\cdot\left[\beta^{S0}(\vartheta)\nabla\Upsilon\right]}_{\text{MatDiffusion}}
\end{equation}

\begin{equation}
\underbrace{\frac1{L_\vartheta}\frac{\partial\vartheta}{\partial t} }_{\text{SusceptibilityTimeDerivative}} =
\underbrace{-\frac{\partial\psi^l}{\partial\vartheta}}_{\text{PolarPFMDerivative}}
\underbrace{-\frac12\frac{\partial\beta^{S0}(\vartheta)}{\partial\vartheta}|\nabla\Upsilon|^2}_{\text{PolarPFMGradient}}
+ \underbrace{\nabla\cdot\left[\beta^{21}\phi(\Upsilon,a_\phi,a_0)\nabla\vartheta\right]}_{\text{MatDiffusion}}
\end{equation}

This model requires the implementation of two new kernels
[PolarPFMDerivative](/PolarPFMDerivative.md) and [PolarPFMGradient](/PolarPFMGradient.md)
for the terms indicated. The remaining terms can be modeled using the existing
[SusceptibilityTimeDerivative](/SusceptibilityTimeDerivative.md) and
[MatDiffusion](/MatDiffusion.md) kernels.

Sub-terms appearing the equations above are implemented as materials

|Material  |  Term |
| - | - |
| [PolarPFMBetaS0](/PolarPFMBetaS0.md) | $\beta^{S0}(\vartheta)$ |
| [PolarPFMPhi](/PolarPFMPhi.md) | $\beta^{21}\phi(\Upsilon,a_\phi,a_0)$ |
| [PolarPFMPsiL](/PolarPFMPsiL.md) | $\psi^l$ |

## Physical meaning

The model implements a two sold phase system with a third interfacial melt
phase. The $\vartheta$ order parameter switched between solid1 ($\vartheta=0$)
and solid2 ($\vartheta=0$). The $\Upsilon$ order parameter switches between
solid ($\vartheta=0$ - governed by $\vartheta$) and interfacial melt
($\vartheta=1$). The indices $s$ and $1,2$ refer to an unspecified and specified
solid phase respectively, while the index $0$ refers to the interfacial melt
phase.

Both order parameters $\vartheta$ and $\Upsilon$ are non-conserved
\emph{Allen-Cahn like} phase order parameters.

The parameters in the test and example files are using the parameterization for
HMX from the original publication.

!bibtex bibliography
