[Mesh]
  type = GeneratedMesh
  dim = 2
[]

[Variables]
  [./c]
    initial_condition = 1.0
  [../]
[]

[Problem]
  kernel_coverage_check = false
[]

[UserObjects]
  [./thermal_fission]
    type = PKAFissionFragmentEmpirical
    relative_density = 1
    fission_rate = 1e5
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = c
    M = 40
    Z = 20
    site_volume = 0.0404 # nm^3 per UO2 unit
    pka_generator = thermal_fission
  [../]
[]

[VectorPostprocessors]
  [./mass]
    type = MyTRIMPKAStatistics
    rasterizer = rasterizer
    execute_on = timestep_end
    value_type = MASS
  [../]
  [./zaid]
    type = MyTRIMPKAStatistics
    rasterizer = rasterizer
    execute_on = timestep_end
    value_type = ZAID
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
[]

[Outputs]
  csv = true
  execute_on = 'TIMESTEP_END'
[]
