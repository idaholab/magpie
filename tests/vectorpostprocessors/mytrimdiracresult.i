[Mesh]
  type = MyTRIMMesh
  dim = 2
  nx = 60
  ny = 60
  xmin = -1000
  xmax = 1000
  ymin = 0
  ymax = 2000
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
    E = 200
    Z = 40
    m = 80
    pka_rate = 2e-3
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = c
    M = 40
    Z = 20
    site_volume = 0.0404 # nm^3 per UO2 unit
    pka_generator = thermal_fission
  [../]
  [./runner]
    type = MyTRIMDiracRun
    rasterizer = rasterizer
  [../]
[]

[VectorPostprocessors]
  [./int]
    type = MyTRIMDiracResult
    runner = runner
    defect = INT
    ivar = 0
  [../]
  [./vac]
    type = MyTRIMDiracResult
    runner = runner
    defect = VAC
    ivar = 0
  [../]
  [./pka]
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
