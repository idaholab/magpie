[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50
  ny = 50
  xmin = -100
  xmax = 100
  ymin = 0
  ymax = 200
[]

[Variables]
  [./T]
  [../]
[]

[Kernels]
  [./source]
    type = MyTRIMElementHeatSource
    variable = T
    runner = runner
  [../]
  [./conduction]
    type = Diffusion
    variable = T
  [../]
  [./time]
    type = TimeDerivative
    variable = T
  [../]
[]

[AuxVariables]
  [./edep]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./c]
    initial_condition = 1.0
  [../]
[]

[AuxKernels]
  [./edep]
    variable = edep
    type = MyTRIMElementEnergyAux
    runner = runner
  [../]
[]

[UserObjects]
  [./thermal_fission]
    type = PKAFissionFragmentEmpirical
    relative_density = 1
    fission_rate = 0.01
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = c
    M = 40
    Z = 20
    site_volume = 0.0404 # nm^3 per UO2 unit
    pka_generator = thermal_fission
    trim_module = ENERGY_DEPOSITION
  [../]
  [./runner]
    type = MyTRIMElementRun
    rasterizer = rasterizer
  [../]
[]

[Postprocessors]
  [./E]
    type = ElementIntegralVariablePostprocessor
    variable = edep
    execute_on = timestep_end
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 3
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv = true
  hide = c
[]
