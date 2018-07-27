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
  [./c]
    initial_condition = 1.0
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      variable = c
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
    type = PKAFixedPointGenerator
    m = 131
    Z = 54
    E = 100000 # 100keV
    point = '0 0 0'
    num_pkas = 10
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = c
    M = 40
    Z = 20
    site_volume = 0.0404 # nm^3 per UO2 unit
    pka_generator = thermal_fission
    r_rec = 5
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
