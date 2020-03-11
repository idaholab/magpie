# ParkinCoulterDPAUserObject

This userobject computes the displacement per atoms
caused by neutron irradiation for displacement functions $\nu_{i,j}(T)$
computed by a Parkin-Coulter displacement function. In these MOOSE Background information on DPA calculation can be found [here](/DPAUserObjectBase.md).

## Example Input File Syntax

!listing tests/userobjects/dpa/huang_ghoniem_92.i 
block=UserObjects

!syntax description /UserObjects/ParkinCoulterDPAUserObject

!syntax parameters /UserObjects/ParkinCoulterDPAUserObject

!syntax inputs /UserObjects/ParkinCoulterDPAUserObject

!syntax children /UserObjects/ParkinCoulterDPAUserObject

!bibtex bibliography
