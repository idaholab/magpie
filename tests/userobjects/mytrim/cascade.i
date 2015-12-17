[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 100
  xmin = -15
  xmax = 15
  ymin = -5
  ymax = 25
[]

[Variables]
  [./c]
    initial_condition = 1.0
  [../]
[]

[Kernels]
  [./dt]
    type = TimeDerivative
    variable = c
  [../]
  [./diff]
    type = Diffusion
    variable = c
  [../]
[]

[AuxVariables]
  [./int]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./int]
    variable = int
    type = MyTRIMAux
    runner = runner
    ivar = 0
    defect = INT
  [../]
[]

[UserObjects]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = c
    M   = 40
    Z   = 20
  [../]
  [./runner]
    type = MyTRIMRun
    rasterizer = rasterizer
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
