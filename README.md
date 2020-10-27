# Magpie

[![Build status](https://civet.inl.gov/idaholab/magpie/master/branch_status.svg)](https://civet.inl.gov/repo/idaholab/magpie/)


![Magpie Logo](http://idaholab.github.io/img/magpie/magpie_logo_small.png)

**Magpie** (_Mesoscale Atomistic Glue Program for Integrated Execution_) is a glue application to link various atomistic codes to applications based on the [MOOSE Finite Element Framework](http://mooseframework.org).

Magpie provides coupling modules for

* [SPPARKS](http://spparks.sandia.gov/) - the kinetic Monte Carlo code by Sandia National Laboratory
* [MyTRIM](http://github.com/idaholab/mytrim/) - a binary collision Monte Carlo code for ion transport in materials

Support is planned for
* [LAMMPS](http://lammps.sandia.gov/) - the molecular dynamics code by Sandia National Laboratory

and possibly others.

### Other Software

Idaho National Laboratory is a cutting edge research facility which is a constantly producing high quality research and software. Feel free to take a look at our other software and scientific offerings at:

[Primary Technology Offerings Page](https://www.inl.gov/inl-initiatives/technology-deployment)

[Supported Open Source Software](https://github.com/idaholab)

[Raw Experiment Open Source Software](https://github.com/IdahoLabResearch)

[Unsupported Open Source Software](https://github.com/IdahoLabCuttingBoard)

### Contrib

Part of this code uses the exact calculation of the [overlap](contrib/overlap/README.md) volume of spheres
and mesh elements routines copyright (C) 2015-2018 by Severin Strobl
<git@severin-strobl.de> (http://dx.doi.org/10.1016/j.jcp.2016.02.003).

### License

Copyright 2017 Battelle Energy Alliance, LLC

This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
