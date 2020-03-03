# PolarPFMGradient

!syntax description /Kernels/PolarPFMGradient

This kernel implements the term

\begin{equation}
-\frac12 \frac{\partial F}{\partial u} |\nabla v|^2
\end{equation}

in the polar phase field model from [!cite](Momeni2014). This term is used in the
evolution equations of both order parameters.

| $u$ (`variable`)| $v$ (`v`)| $F$ (`F`) |
| - | - | - |
| $\Upsilon$ | $\vartheta$ | [PolarPFMPhi](/PolarPFMPhi.md) |
| $\vartheta$ | $\Upsilon$ | [PolarPFMBetaS0](/PolarPFMBetaS0.md) |

!syntax parameters /Kernels/PolarPFMGradient

!syntax inputs /Kernels/PolarPFMGradient

!syntax children /Kernels/PolarPFMGradient

!bibtex bibliography
