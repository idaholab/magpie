[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 20
  xmin = -100
  xmax = 100
  ymin = -100
  ymax = 100
[]

[Variables]
  [./c_U]
    initial_condition = 1.0
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      variable = c_U
      auto_direction = 'x y'
    [../]
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
  [./c_O]
    initial_condition = 2.0
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
  [./func_pka]
    type = PKAFunction
    m = 238
    Z = 54
    E = '1000*10^t'
    point = '0 0 0'
    pka_rate = (0.6-t)*0.1
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = 'c_U  c_O'
    M   = '235  16 '
    Z   = '92   8  '
    site_volume = 0.0404 # nm^3 per UO2 unit
    pka_generator = func_pka
    r_rec = 0
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
  [./npka]
    type = MyTRIMPKAInfo
    rasterizer = rasterizer
    value_type = TOTAL_NUMBER
  [../]
  [./E]
    type = FunctionValuePostprocessor
    function = '1000*10^t'
  [../]
[]

[Executioner]
  type = Transient
  dt = 0.1
  num_steps = 5
  nl_abs_tol = 1e-10
[]

[Outputs]
  csv = true
[]
