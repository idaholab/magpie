# MyTRIMMesh

!syntax description /Mesh/MyTRIMMesh

Restricted regular orthogonal [GeneratedMesh](/GeneratedMesh.md) with equal size
level 0 elements (no bias) and EDGE2, QUAD4, or HEX8 elements only. These
restrictions allow the use of a custom super fast point locator to find elements
based on spatial locations (required for fast material property lookup in the
Binary Collision Monte Carlo stage). This mesh is also recommended for use with
the FourierTransform user object.

The use of adaptivity is fully supported.

!syntax parameters /Mesh/MyTRIMMesh

!syntax inputs /Mesh/MyTRIMMesh

!syntax children /Mesh/MyTRIMMesh

!bibtex bibliography
