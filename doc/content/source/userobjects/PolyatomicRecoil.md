# PolyatomicRecoil

!syntax description /UserObjects/PolyatomicRecoil

# Description

PolyatomicRecoil allows computation of total and net displacement functions described in [!cite](PC1981) and damage energy functions defined in [!cite](PC1980). PolyatomicRecoil will
print output in csv format. In the future the preparation of displacement cross sections, described in
[!cite](HG1993), will be implemented. PolyatomicRecoil also implements the computation of first derivatives of the net displacement function with respect to number densities.

PolyatomicRecoil allows the displacement and capture energies
of different species at a particular lattice to be different, e.g. the capture energy of species i at an
i-site may be different than the capture energy of species j at an i-site. Integration of the
integro-differential equations uses GSL's ODE and numerical integration capabilities; an in-depth discussion
of the utilized algorithms is contained in "doc/parkin_coulter/parkin_coulter.pdf".

Energy spacing is logarithmic (constant ratio) except for energies $$$E < E_t$$$, where $$$E_t$$$ is set via the uniform_energy_spacing_threshold parameter.

Stopping powers are computed using mytrim routines; Thomas-Fermi potential and Lindhard's expression for the scattering cross section are used.

!syntax parameters /UserObjects/PolyatomicRecoil

!bibtex bibliography
