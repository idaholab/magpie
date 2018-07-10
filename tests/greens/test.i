[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = -10
  ymin = -10
  xmax = 10
  ymax = 10
  nx = 20
  ny = 20
[]

[Variables]
  [./c]
    [./InitialCondition]
      type = FunctionIC
      function = 'if((x+6)^2+(y-5)^2<16,1,0)'
    [../]
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'x y'
    [../]
  [../]
[]

[Kernels]
  [./dt]
    type = TimeDerivative
    variable = c
  [../]
  [./cc1]
    type = RadialGreensSource
    variable = c
    gamma = 1
    convolution = green1
  [../]
  [./cc2]
    type = RadialGreensSource
    variable = c
    gamma = 1e-1
    convolution = green2
  [../]
[]

[UserObjects]
  [./green1]
    type = RadialGreensConvolution
    execute_on = TIMESTEP_BEGIN
    v = c
    r_cut = 8
    function = 'exp(-x/2)'
    normalize = true
  [../]
  [./green2]
    type = RadialGreensConvolution
    execute_on = TIMESTEP_BEGIN
    v = c
    r_cut = 6.28
    function = 'sin(x)'
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

[Postprocessors]
  [./C]
    type = ElementIntegralVariablePostprocessor
    execute_on = 'INITIAL TIMESTEP_END'
    variable = c
  [../]
  [./CC1]
    type = ElementIntegralVariablePostprocessor
    execute_on = 'INITIAL TIMESTEP_END'
    variable = cc1
  [../]
  [./CC2]
    type = ElementIntegralVariablePostprocessor
    execute_on = 'INITIAL TIMESTEP_END'
    variable = cc2
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
[]

[Outputs]
  exodus = true
  execute_on = FINAL
[]
