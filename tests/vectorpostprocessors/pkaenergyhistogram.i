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
  [./histo]
    type = MyTRIMPKAEnergyHistogram
    rasterizer = rasterizer
    execute_on = timestep_end
    channel_number = 80
    channel_width = 2.5e+06
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
