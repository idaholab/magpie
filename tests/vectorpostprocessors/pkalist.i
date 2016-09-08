[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 60
  ny = 60
  xmin = -100
  xmax = 100
  ymin = 0
  ymax = 200
[]

[AuxVariables]
  [./c]
    initial_condition = 1.0
  [../]
[]

[Variables]
  [./u]
  [../]
[]

[Problem]
  kernel_coverage_check = false
[]

[UserObjects]
  [./thermal_fission]
    type = PKAConstant
    E = 200000
    Z = 40
    m = 80
    pka_rate = 5e-4
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
  [./pkas]
    type = PKAList
    rasterizer = rasterizer
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
[]

[Outputs]
  execute_on = 'TIMESTEP_END'
  csv = true
[]
