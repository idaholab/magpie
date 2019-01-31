[Mesh]
  type = GeneratedMesh
  dim = 3
  xmin = 0
  xmax = 1
  nx = 1
  ymin = 0
  ymax = 1
  ny = 1
  zmin = 0
  zmax = 1
  nz = 1
  elem_type = TET4
[]

[Problem]
  kernel_coverage_check = false
[]

[Variables]
  [./d]
  [../]
[]

[AuxVariables]
  [./packing]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[AuxKernels]
  [./packing_fraction]
    type = MDGranularPorosityAux
    variable = packing
    compute_packing_fraction = true
    user_object = lammps_runner
  [../]
[]

[UserObjects]
  [./lammps_runner]
    type = LAMMPSFileRunner
    lammps_file = 'single.0.xyz'
    xyz_columns = '0 1 2'
    md_particle_properties = 'radius'
    property_columns = '3'
  [../]
[]

[Postprocessors]
  [./sphere_volume]
    type = ElementIntegralVariablePostprocessor
    variable = packing
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
  dt = 1
[]

[Outputs]
  exodus = true
[]
