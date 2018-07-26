[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = 0
  ymin = 0
  xmax = 100
  ymax = 100
  nx = 100
  ny = 100
  parallel_type = DISTRIBUTED
[]

[AuxVariables]
  [./rank]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[Variables]
  [./c]
    initial_condition = 1
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
    convolution = green
  [../]
[]

[AuxKernels]
  [./rank]
    type = ProcessorIDAux
    variable = rank
  [../]
[]

[UserObjects]
  [./green]
    type = RadialGreensConvolution
    execute_on = TIMESTEP_BEGIN
    v = c
    r_cut = 5.0001
    function = '5.0001-x'
    normalize = true
  [../]
[]


[Postprocessors]
  [./c_min]
    type = ElementExtremeValue
    value_type = min
    execute_on = 'INITIAL TIMESTEP_END'
    variable = c
  [../]
  [./c_max]
    type = ElementExtremeValue
    value_type = max
    execute_on = 'INITIAL TIMESTEP_END'
    variable = c
  [../]
[]

[Executioner]
  type = Transient
  nl_abs_tol = 1e-10
  num_steps = 1
[]

[Outputs]
  csv = true
  execute_on = FINAL
  [./perf]
    type = PerfGraphOutput
    level = 3
  [../]
[]
