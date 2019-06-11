[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 20
[]

[AuxVariables]
  [./T]
    [./InitialCondition]
      type = FunctionIC
      function = y/2+0.1
    [../]
  [../]
  [./c]
    [./InitialCondition]
      type = FunctionIC
      function = x
    [../]
  [../]
[]

[Materials]
  [./Fn]
    type = NeuralNetFreeEnergy
    file_name = weights_biases.txt
    file_format = MAGPIE
    inputs     = 'T c'
    prop_names = 'F'
    outputs = exodus
  [../]
[]

[Problem]
  kernel_coverage_check = false
  solve = false
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
[]
