# FunctionDPAUserObject

This userobject computes the displacement per atoms
caused by neutron irradiation for displacement functions $\nu_{i,j}(T)$
provided by functions. In these MOOSE functions the time argument
is reinterpreted as the energy ($T$) slot. Background information can be found [here](/DPAUserObjectBase.md).

## Example Input File Syntax

!listing tests/userobjects/dpa/verification_1.i
block=UserObjects

!syntax description /UserObjects/FunctionDPAUserObject

!syntax parameters /UserObjects/FunctionDPAUserObject

!syntax inputs /UserObjects/FunctionDPAUserObject

!syntax children /UserObjects/FunctionDPAUserObject

!bibtex bibliography
