[Mesh]
  type = MyTRIMMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  uniform_refine = 2
[]

[Variables]
  [./c]
    [./InitialCondition]
      type = RandomIC
      min = 0.1
      max = 1
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
    Z = 60
    m = 120
    pka_rate = 10
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = c
    M = 40
    Z = 20
    site_volume = 0.0404 # nm^3 per UO2 unit
    pka_generator = constant
  [../]
  [./runner]
    type = MyTRIMElementRun
    rasterizer = rasterizer
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
  nl_abs_tol = 1e-10
[]

[Problem]
  kernel_coverage_check = false
[]

[Outputs]
  exodus = true
[]
