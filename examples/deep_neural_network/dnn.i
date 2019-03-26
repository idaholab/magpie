[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
[]

[Variables]
  [./T]
    [./InitialCondition]
      type = FunctionIC
      function = x
    [../]
  [../]
  [./c]
    [./InitialCondition]
      type = FunctionIC
      function = y
    [../]
  [../]
[]

[Problem]
  kernel_coverage_check = false
  solve = false
[]

[Materials]
  [./F]
    type = DeepNeuralNetFreeEnergy
    filename = weights_biases.txt
    inputs     = 'T c'
    prop_names = 'F'
    outputs = exodus
  [../]
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
[]

# ./maagpie-opt -i dnn.i | grep '^DNN' | cut -c4- > out
# gnuplot:
# spl 'out' w d, y*(1-y)+0.001*(x*300+400)*(y*log(y)+(1-y)*log(1-y))
