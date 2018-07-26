# RadialGreensConvolution

!syntax description /UserObjects/RadialGreensConvolution

Given a Green's function $G(r)$ and a variable field $c(\vec r)$ this userobject
calculates a convolution $c'(\vec r)$ such that

\begin{equation}
c'(\vec r) = \frac1{4\pi}\int_{\Omega(\vec{r})} \frac{G(|\vec r - \vec{r'}|)}{|\vec r - \vec{r'}|^2} c(\vec{r'})\,d\vec{r'},
\end{equation}

where the integration domain $\Omega$ is defined as the ball with radius $r_{cut}$ (`r_cut`)

\begin{equation}
\Omega(\vec r) = \left\{\vec{r'} \in \mathbb{R}^3: |\vec r - \vec{r'}| < r_{cut}\right\}.
\end{equation}

!alert note
Note that the geometrical attenuation $\frac1{4\pi |\vec r - \vec{r'}|^2}$ is not part of $G$.

!alert note
The convolution will always be computed in 3D, no matter what the mesh dimension is.

## Applications

$G(r)$ could be a probability distribution for ballistically knocking an atom in
the sample for a distance $r$. Such a function could be obtained using a binary
collision Monte Carlo code.

Applying this convolution using a [RadialGreensSource](/RadialGreensSource.md)
kernel would then correspond to applying ballistic mixing to the sample. The
probability of hitting an atom at all would be cntrolled using the `gamma`
parameter in the kernel.

## Design

The RadialGreensFunction user object is derived from `ElementUserObject` and
works in two stages.

1. In the element loop (in the `execute()` method) a list of all quadrature
   points, their locations, indices, and selected variable value is compiled.

2. In the `finalize()` method

    1. the list is communicated in parallel

    2. a KD-tree is filled witha ll quadrature point entries (utilizing the
        [nanoflann](https://github.com/jlblancoc/nanoflann) library bundled with
        libMesh)

    3. a loop over all QPs is performed and at each QP a $r_{cut}$ (`r_cut`)
        radius search in the KD-tree is performed

    4. the results from the search are used to spatially integrate the selected
        variable field multiplied with the Green's Function $G$ and the geometric
        attenuation $\left(4\pi |\vec r - \vec{r'}|^2\right)^{-1}$.

    5. The integral is normalized to ensure mass conservation

!syntax parameters /UserObjects/RadialGreensConvolution

!syntax inputs /UserObjects/RadialGreensConvolution

!syntax children /UserObjects/RadialGreensConvolution

!bibtex bibliography
