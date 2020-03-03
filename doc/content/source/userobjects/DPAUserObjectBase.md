# DPAUserObjectBase

## Description

A base class for computing the displacement per atom (DPA) caused by neutron irradiation. The DPA calculation requires uses the following
information:

- the scalar flux $\phi(E)$, where $E$ denotes neutron energy.
- the number density of isotope $i$ denoted by $N_i$. The index $i$ runs over all isotopes $i=1,..,I$.
- the microscopic cross section of isotope $i$ for reaction type $x$ (currently $x$ can stand for elastic and inelastic scattering).
- the atomic number $Z_i$ and mass number $A_i$ of isotope $i$.
- the number of displacements of type $j$ that a primary knock-on atom (PKA) of type $i$ and energy $T$ causes is denoted by $\nu_{i,j}(T)$.

The increase in dose $\Delta D$ measured in dpa after a time interval $\Delta t$ is given by:
\begin{equation}
   \Delta D = \Delta t \sum\limits_x \sum\limits_{i=1}^I
   \sum\limits_{j=1}^J N_i
  \int_{0}^{\infty} \sigma_{i,x}(E)  \phi(E)
  \left[ \int_{T_{\text{min}, i, x}}^{T_{\text{max}, i, x}}  \nu_{i,j}(T) H_{i,x}(E,T) dT \right] dE,
\end{equation}
where $H_x(E,T)$ is the probability density function that a reaction of type $x$ causes a PKA of type $i$ with energies in the interval $dT$ about $T$; $T_{\text{min}, i, x}$ and $T_{\text{max}, i, x}$ are the maximum and minimum recoil energy that a neutron in reaction channel $x$ can cause in a collision with isotope $i$.

It is convenient to define the parameter $\gamma_i$:
\begin{equation}
  \gamma_i = \frac{4 A_i}{(A_i + 1)^2}
\end{equation}

For inelastic reactions, the reduced energy is defined:
\begin{equation}
   \delta_i = \frac{Q_i}{E} \frac{A_i + 1}{A_i}
\end{equation}

For elastic and inelastic scattering:

- Elastic scattering: we assume isotopic scattering in the center of mass frame; then we have $H_{i,x}(E,T) = \frac{1}{\gamma_i E}$, $T_{\text{min}, i, x} = 0$, $T_{\text{max}, i, x} = \gamma_i E$.
- Inelastic scattering: $H_{i,x}(E,T) = \frac{1}{\gamma_i E\sqrt(1 - \delta_i)}$,

\begin{equation}
\begin{aligned}
  T_{\text{min}, i, x} &=\frac{\gamma_i E}{2}\left(1 - \sqrt{1 - \delta_i}\right)-\frac{1}{2}\delta_i \\
  T_{\text{max}, i, x} &=\frac{\gamma_i E}{2}\left(1 + \sqrt{1 - \delta_i}\right)-\frac{1}{2}\delta_i.  
\end{aligned}
\end{equation}

We observe that in both cases $H_{i,x}(E,T) = H_{i,x}(E)$.

From deterministic transport theory, usually the group flux $\phi_g$ is provided. It is defined by
\begin{equation}
  \phi_g = \int_{E_{g+1}}^{E_{g}} \phi(E) dE,
\end{equation}
where $E_g$ and $E_{g+1}$ are energy group boundaries ordered from high to low, i.e. $E_{g+1} < E_g$. We make the approximation:
\begin{equation}
  \phi(E) \approx \frac{\phi_g}{E_{g} - E_{g+1}} \text{ for } E_{g+1} \le E < E_{g}.
\end{equation}
The group-averaged cross section is defined by:
\begin{equation}
   \sigma_{g,i,x} = \frac{1}{\phi_g} \int_{E_{g+1}}^{E_{g}} \sigma_{i,x}(E) \phi(E) dE.
\end{equation}

Then the dpa increment $\Delta D$ can be approximated by:
\begin{equation}
  D = \delta t \Delta t \sum_{g} \sum\limits_x \sum\limits_{i=1}^I
  \sum\limits_{j=1}^J \frac{N_i \sigma_{g,i,x} \phi_g}{E_{g} - E_{g+1}} \Gamma_{i,j,g,x},
\end{equation}
where
\begin{equation}
  \Gamma_{i,j,g,x} = \int_{E_{g+1}}^{E_{g}} H_{i,x}(E) \left[\int_{T_{\text{min}, i, x}}^{T_{\text{max}, i, x}}  \nu_{i,j}(T)  dT  \right]  dE.
\end{equation}

The approximation of using constant $\phi(E)$ within each energy group is expected to lead to large errors if the fast energy range is represented by too few groups. An improvement is needed that accounts for the variation of $\phi(E)$ within each energy group.
