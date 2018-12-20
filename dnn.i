[Mesh]
  type = GeneratedMesh
  dim = 1
[]

[Variables]
  [./T]
  [../]
  [./c]
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
    inputs = 'T c'
  [../]
[]

[Executioner]
  type = Steady
[]


# ./maagpie-opt -i dnn.i | grep '^DNN' | cut -c4- > out
# gnuplot:
# spl 'out' w d, y*(1-y)+0.001*(x*300+400)*(y*log(y)+(1-y)*log(1-y))
