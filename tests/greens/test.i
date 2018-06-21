[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 2
  ny = 2
[]

[AuxVariables]
  [./c]
    [./InitialCondition]
      type = FunctionIC
      function = 'x + 10.0 * y'
    [../]
  [../]
[]

[Variables]
  [./dummy]
  [../]
[]

[UserObjects]
  [./green1]
    type = RadialGreensGather
    v = c
    r_cut = 1
    function = exp(-10*x)
  [../]
[]

[Problem]
  kernel_coverage_check = false
  solve = false
[]

[Executioner]
  type = Transient
  num_steps = 1
[]
