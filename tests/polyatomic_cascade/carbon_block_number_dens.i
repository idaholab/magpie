[Mesh]
  type = MyTRIMMesh
  dim = 3
  nx = 5
  ny = 5
  nz = 5
  xmin = 0
  xmax = 1e5
  ymin = 0
  ymax = 1e5
  zmin = 0
  zmax = 1e5
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

[UserObjects]
  [./constant]
    type = PKAFixedPointGenerator
    E = 1000
    Z = 6
    m = 12
    point = '5e4 5e4 5e4'
    num_pkas = 1000
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = 'cC'
    M = '12'
    Z = '6'
    Edisp = '28'
    pka_generator = constant
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
