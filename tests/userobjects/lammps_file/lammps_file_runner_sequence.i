[Mesh]
  type = GeneratedMesh
  dim = 3
  xmin = 0
  xmax = 2
  nx = 2
  ymin = 0
  ymax = 2
  ny = 2
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
    lammps_file = 'sequence/simple'
    time_sequence = true
    xyz_columns = '0 1 2'
    time_conversion = tc
  [../]
[]

[Functions]
  [./tc]
    type = ParsedFunction
    value = '2 * t + 5'
  [../]
[]

[Executioner]
  type = Transient
  end_time = 3
  dt = 0.1
[]

[Outputs]
  exodus = true
[]
