[Mesh]
  type = MyTRIMMesh
  dim = 3
  nx = 2
  ny = 2
  nz = 2
  xmin = 0
  xmax = 100
  ymin = 0
  ymax = 100
  zmin = 0
  zmax = 100
  elem_type = HEX8
[]

[Problem]
  kernel_coverage_check = false
[]

[Variables]
  [./dummy]
  [../]
[]

[AuxVariables]
  [./cC]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1
  [../]

  [./vac]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./vac]
    variable = vac
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 0
    defect = VAC
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'x y z'
    [../]
  [../]
[]

[UserObjects]
  [./constant]
    type = PKAConstant
    E = 1000
    Z = 6
    m = 12
    pka_rate = 0.005
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = 'cC'
    M = '12'
    Z = '6'
    Edisp = '16.3'
    site_volume = 0.0404
    pka_generator = constant
    periodic_var = dummy
  [../]
  [./runner]
    type = MyTRIMElementRun
    rasterizer = rasterizer
  [../]
[]

[Postprocessors]
  [./total_vacancies]
    type = ElementIntegralVariablePostprocessor
    variable = vac
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
[]
