[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = -10
  ymin = -10
  xmax = 10
  ymax = 10
  nx = 60
  ny = 60
[]

[Variables]
  [./c]
    [./InitialCondition]
      type = FunctionIC
      function = 'if(x^2+y^2<16,1,0)'
    [../]
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

[UserObjects]
  [./green1]
    type = RadialGreensConvolution
    execute_on = TIMESTEP_BEGIN
    v = c
    r_cut = 4
    function = 'if(x<=1e-9,1.0,1e-5*exp(-x))'
    normalize = true
  [../]
  [./green2]
    type = RadialGreensConvolution
    execute_on = TIMESTEP_BEGIN
    v = c
    r_cut = 6.28
    function = 'if(x<=1e-9,1.0,1e-5*sin(x))'
    normalize = true
  [../]
[]

[AuxVariables]
  [./cc1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./cc2]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./cc1]
    type = RadialGreensAux
    variable = cc1
    convolution = green1

  [../]
  [./cc2]
    type = RadialGreensAux
    variable = cc2
    convolution = green2
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
  dt = 1e-6
[]

[Outputs]
  exodus = true
[]
