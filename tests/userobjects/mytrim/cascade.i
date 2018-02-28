[Mesh]
  type = MyTRIMMesh
  dim = 2
  nx = 20
  ny = 20
  xmin = -100
  xmax = 100
  ymin = 0
  ymax = 200
[]

[Variables]
  [./c]
    initial_condition = 1.0
  [../]
[]

[AuxVariables]
  [./int]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vac]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Problem]
  kernel_coverage_check = false
[]

[AuxKernels]
  [./int]
    variable = int
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 0
    defect = INT
  [../]
  [./vac]
    variable = vac
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 0
    defect = VAC
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
  [../]
  [./runner]
    type = MyTRIMElementRun
    rasterizer = rasterizer
  [../]
[]

[Postprocessors]
  [./int]
    type = ElementIntegralVariablePostprocessor
    variable = int
    execute_on = timestep_end
  [../]
  [./vac]
    type = ElementIntegralVariablePostprocessor
    variable = vac
    execute_on = timestep_end
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv = true
  hide = c
[]
