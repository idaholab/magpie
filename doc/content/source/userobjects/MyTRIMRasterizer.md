# MyTRIMRasterizer

## Providing isotopic species and concentrations

`MyTRIMRasterizer` uses information of the isotopic composition in the domain
provided by MooseVariables specified in parameter `var` to prepare a rasterized
domain on which binary collision Monte Carlo is performed. Rasterization is the
process of averaging smoothly varying compositional data over each FEM mesh element.

The entries in of the `var`, `M`, and `Z` parameters are paired; we index them
as $\gamma_j(\vec{r})$, $M_j$, and $Z_j$, respectively. Loosely speaking
$\gamma_j(\vec{r})$ represents the number density of isotope $j$ at point
$\vec{r}$, $M_j$ is the atomic mass, and $Z_j$ is the charge in units of $e$.

The physical meaning of the `var` argument depends on what the user provides for
the `site_volume` parameter. `site_volume` is a material property and can vary
with space; it is denoted by $\Omega(\vec{r})$. MyTRIM computes the mass density
$\rho$ of the material by:

\begin{equation}
  \rho = \sum\limits_{j=1}^J \frac{\gamma_j(\vec{r}) M_j}{N_A \Omega(\vec{r})},
\end{equation}

where $N_A$ is Avogadro's number. We identify $\gamma_j / \Omega$ as the number
density $N_j$ defined as the number of atoms of species $j$ per unit volume. The
following choices for $\Omega$ may be usage of the `site_volume` parameter:

$\Omega = 1$: $\gamma_j = N_j$ so that the user should store number densities in
the variables specified in `var`.

$\Omega \equiv \text{compound volume, e.g. volume of UO}_2$: $\gamma_j$ should
be the stoichiometric content of species $j$. Using the $\text{UO}_2$ example:
$\gamma_U = 1$, $\gamma_O=2$.

$\Omega \equiv \text{average volume per atoms, i.e. }
\frac{\text{volume}}{\text{number of atoms}}$: $\gamma_j$ should be the number
fraction. In this case $\sum\limits_{j=1}^J \gamma_j=1$.

## Computing the PKA list

Before performing damage cascades, the rasterizer computes the number of primary
knock-on atoms (PKA) $n_{\text{PKA,l, M, Z}$ for each mesh element $l$, mass $M$
and atomic number $Z$ by:

\begin{equation}
 n_{\text{PKA},l,M,Z} = \Delta t \sum\limits_{k=1}^K  p_k(A,Z) \int\limits_{V_l} dV R_k(\vec{r}),
\end{equation}

where $\Delta t$ is the current time step, $k$ indexes the `PKAGenerator`
objects provided in the `pka_generator` parameter, $p_k(M,Z)$ is the probability
that `PKAGenerator` $k$ produces a PKA with mass $M$ and atomic number $Z$, and
$R_k(\vec{r})$ is the recoil rate of `PKAGenerator` $k$ that may depend on
space. $n_{\text{PKA},l,M,Z}$ is an integer number so it is rounded to the next
higher value with a probability proportional to its truncated integer complement
($1.4$ would be $2$ with a probability of $0.4$, and $1$ with a probability of
$0.6$). The total number of PKAs is given by $n_{\text{PKA}} =
\sum\limits_{l,M,Z} n_{\text{PKA},l,M,Z}$.

## Scaling PKA rates

The rasterizer provides two parameters for modifying the number of simulated
PKAs without biasing the results. The parameter `max_pka_count` limits the
number of PKAs roughly to the provided integer number denoted by $n$. If the
number of computed PKAs is smaller than or equal to `max_pka_count`, the
parameter is ignored. If the number of PKAs is larger than `max_pka_count`, PKAs
are accepted with a probability of $n/n_{\text{PKA}}$. The actual number of PKAs
is usually not identical to $n$ and is denoted by $n'$. Any tallied result is,
e.g. vacancy or interstitial densities or energy deposition, are multiplied by
$n_{\text{PKA}/n'$.

The parameter `recoil_rate_scaling` is a real number that is multiplied to the
recoil rate $R_k(\vec{r})$. All tallied results are divided by
`recoil_rate_scaling`.

## Recombination

A simple model for recombination *within* damage cascades is currently
implemented. A vacancy-interstitial pair is annihilated if their distance is
less than `r_rec`. Note that `r_rec` must be given in Angstrom. Recombination
within damage cascades typically occurs within thermal spike areas where
sufficient energy has been deposited to displace virtually all contained atoms.
In thermal spike areas the atoms are effectively in a liquid state.
Recrystallization takes place when energy has been dissipated from the thermal
spike area and most Frenkel pairs annihilate. Surviving Frenkel pairs have
usually been transported out of the melted zone through collision sequences
along closely packed directions. The current recombination model does not
accurately describe this physics.

## Periodic variables

If all variables in the `var` parameters are auxiliary variables, the rasterizer
has to be made aware of periodic boundaries. The user needs to provide a
nonlinear variable that the desired periodic variable apply to in the
`periodic_var` parameter.

!syntax parameters /UserObjects/MyTRIMRasterizer

!syntax inputs /UserObjects/MyTRIMRasterizer

!syntax children /UserObjects/MyTRIMRasterizer

!bibtex bibliography
