# PKAFissionFragmentNeutronics

`PKAFissionFragmentNeutronics` generates primary knock-on atoms (PKA) originating
from fission reactions. The difference to PKAFissionFragmentEmpirical is that
it uses isotopic fission rates computed by neutronics calculations and ENDF data
for sampling fission product species.

The partial fission rate for isotope $i$ is denoted by $F_i(\vec{r}, \vec{\rho})$.
It is computed by:

\begin{equation}
F_i(\vec{r}, \vec{\rho}) = \sum\limits_{g=1}^G N_i(\vec{r},\vec{\rho})\sigma_{f,g,i} \phi_g(\vec{r}),
\end{equation}

where $\vec{r}$ is the location in the macroscopic, neutronics domain and
$\vec{\rho}$ is the location in the microscopic domain. These two locations are
separated because their scale is significantly separated by 3-4 orders of
magnitude. Changes with $\vec{r}$ are understood to be smooth changes of the
average composition of the microscopic domain (orders of centimeters), while
changes with $\vec{\rho}$ captures compositional changes between grains on the
microscopic domain, i.e. on the orders of micro-meters. The notation of
separating the spatial dependence into two scales follows homogenization theory.
$N_i$ is the number density of isotope $i$, $\sigma_{f,g,i}$ is the microscopic
fission cross section of isotope $i$ in group $g$, and $\phi_g$ is the scalar
flux in group $g$.

Denoting the microscopic domain as $\Theta_m$ located at $\vec{r}_m$, we can define a slowly varying average of the number densities:

\begin{equation}
  N_i(\vec{r}) = \frac{1}{\Theta_m} \int_{\Theta_m} N_i(\vec{r},\vec{\rho}) d\vec{\rho}.
\end{equation}

$N_i(\vec{r})$ is the number density provided to neutronics calculations. The
nuclide fission rates computed by the neutronics calculation is consequently the
slowly varying average given by:

\begin{equation}
F_i (\vec{r}) = \sum\limits_{g=1}^G N_i(\vec{r})\sigma_{f,g,i}.
\end{equation}

The `PKAFissionFragmentNeutronics` accepts the values of $F_i (\vec{r}_m)$ as
the `partial_reaction_rates` parameter. The fission rate density in the
microscopic domain is computed by:

\begin{equation}
  F_i (\vec{r}_m, \vec{\rho}) = \frac{N_i(\vec{r}_m, \vec{\rho})}{N_i(\vec{r}_m)} F_i(\vec{r}_m),
\end{equation}  

where $N_i(\vec{r}_m, \vec{\rho}) = \gamma_i(\vec{\rho}) / \Omega(\vec{\rho})$,
$\Omega(\vec{\rho})$ is the site volume that is provided to the
`MyTRIMRasterizer` and $\gamma_i$ is the $i$-th variable provided to the
`MyTRIMRasterizer's` `var` parameter. Hence, we have:

\begin{equation}
  F_i (\vec{r}_m, \vec{\rho}) = \frac{\gamma_i(\vec{\rho})}{\Omega(\vec{r}) N_i(\vec{r}_m)} F_i(\vec{r}_m),
\end{equation}  

`PKAFissionFragmentNeutronics` accepts the values of $N_i(\vec{r}_m)$ in the
`averaged_number_densities` parameter.

The expected number of fissions in a mesh element with index $j$ and volume
$V_j$ is given by:

\begin{equation}
C_i = \Delta t \int_{V_j} F_i (\vec{r}_m, \vec{\rho})  d\vec{\rho}.
\end{equation}

Non-integer results are rounded up with a probability of the $C_i -
\text{int}(C_i)$; otherwise rounded down. Two PKAs are created from each fission
event. The algorithm for sampling their type, energy, and direction of motion is
described in [cite:SchunertTREATHeatSource].

!syntax parameters /UserObjects/PKAFissionFragmentNeutronics

!syntax inputs /UserObjects/PKAFissionFragmentNeutronics

!syntax children /UserObjects/PKAFissionFragmentNeutronics

!bibtex bibliography
