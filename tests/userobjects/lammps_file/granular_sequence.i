[Mesh]
  type = GeneratedMesh
  dim = 3
  xmin = 0
  xmax = 1
  nx = 2
  ymin = 0
  ymax = 1
  ny = 2
  zmin = 0
  zmax = 1
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
  [./packing]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[AuxKernels]
  [./packing_fraction]
    type = MDGranularPropertyAux
    variable = packing
    average_type = granular_density
    md_particle_property = charge
    user_object = lammps_runner
  [../]
[]

[UserObjects]
  [./lammps_runner]
    type = LAMMPSFileRunner
    lammps_file = 'granular'
    time_sequence = true
    xyz_columns = '0 1 2'
    md_particle_properties = 'radius charge'
    property_columns = '3 4'
  [../]
[]

[Postprocessors]
  [./mid1]
    type = PointValue
    variable = packing
    point = '0.25 0.25 0.25'
  [../]
  [./mid2]
    type = PointValue
    variable = packing
    point = '0.75 0.25 0.25'
  [../]
  [./mid3]
    type = PointValue
    variable = packing
    point = '0.75 0.75 0.25'
  [../]
  [./mid4]
    type = PointValue
    variable = packing
    point = '0.25 0.75 0.25'
  [../]
  [./mid5]
    type = PointValue
    variable = packing
    point = '0.25 0.25 0.75'
  [../]
  [./mid6]
    type = PointValue
    variable = packing
    point = '0.75 0.25 0.75'
  [../]
  [./mid7]
    type = PointValue
    variable = packing
    point = '0.75 0.75 0.75'
  [../]
  [./mid8]
    type = PointValue
    variable = packing
    point = '0.25 0.75 0.75'
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 2
  dt = 0.5
[]

[Outputs]
  exodus = true
[]
