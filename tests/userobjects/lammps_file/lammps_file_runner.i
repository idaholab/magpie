[Mesh]
  type = GeneratedMesh
  dim = 3
  xmin = -1
  xmax = 2
  nx = 3
  ymin = -2
  ymax = 1
  ny = 3
  zmin = 0
  zmax = 2
  nz = 2
[]

[Problem]
  kernel_coverage_check = false
[]

[Variables]
  [./d]
  [../]
[]

[AuxVariables]
  [./nparticles]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[AuxKernels]
  [./np]
    type = MDNParticleAux
    variable = nparticles
    user_object = lammps_runner
  [../]
[]

[UserObjects]
  [./lammps_runner]
    type = LAMMPSFileRunner
    lammps_file = 'simple.0.xyz'
    xyz_columns = '0 1 2'
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
[]

[Outputs]
  exodus = true
[]
