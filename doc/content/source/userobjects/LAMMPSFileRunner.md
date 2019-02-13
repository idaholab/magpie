# LAMMPSFileRunner

`LAMMPSFileRunner` reads in LAMMPS dumps file and maps them on an FEM mesh.
The `time_sequence` parameter determines if a sequence of LAMMPS files is read
or a single file is loaded and used throughout. If a sequence of files is used,
the `lammps_file` parameter should contain the file base, e.g. `path/to/file/filebase`
if the files are called `path/to/file/filebase.<time>.xyz`.

Time conversion from FEM time to LAMMPS time is facilitated by the `time_conversion`
parameter that accepts a MOOSE function object. Denoting FEM time by `t` and LAMMPS time
by `T`, `time_conversion` is a function `T = F(t)`.

!syntax description /UserObjects/LAMMPSFileRunner

!syntax parameters /UserObjects/LAMMPSFileRunner

!syntax inputs /UserObjects/LAMMPSFileRunner

!syntax children /UserObjects/LAMMPSFileRunner

!bibtex bibliography
