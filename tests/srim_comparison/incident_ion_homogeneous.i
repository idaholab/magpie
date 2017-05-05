#
# ions incident on a homogeneous slab
#
[Mesh]
  type = MyTRIMMesh
  dim = 3
  nx = 11
  ny = 11
  nz = 11
  xmin = 0
  xmax = 100000
  ymin = -50000
  ymax =  50000
  zmin = -50000
  zmax =  50000
[]

[Problem]
  kernel_coverage_check = false
[]

[Variables]
  [./dummys]
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

[UserObjects]
  [./constant]
    type = PKAGun
    E = 1e3
    Z = 14
    m = 28
    num_pkas = 100
    point = '0 0 0'
    direction = '1 0 0'
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = 'cC'
    M = '12'
    Z = '6'
    Edisp = '16.3'
    site_volume = 0.0404
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
