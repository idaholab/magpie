[Mesh]
  type = MyTRIMMesh
  dim = 2
  xmax = 100
  ymax = 100
  nx = 100
  ny = 100
[]

[Problem]
  type = FFTProblem
[]

[AuxVariables]
  [./c_aux]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = FunctionIC
      function = 'cos(x/100*2*pi*4)*cos(y/100*2*pi*3)'
    [../]
  [../]
[]

[Materials]
  [./test]
    type = ParsedMaterial
    args = c
    function = c*c
  [../]
[]

[UserObjects]
  # Buffers
  [./c]
    type = RealFFTWBuffer
    moose_variable = c_aux
  [../]
  [./R]
    type = RankTwoTensorFFTWBuffer
  [../]

  # Solver
  # ...
[]

[AuxKernels]
  [./c_aux]
    type = FFTBufferAux
    variable = c_aux
    fft_buffer = c
    execute_on = FINAL
  [../]
[]

[Executioner]
  type = SpectralExecutionerBase
[]

[Outputs]
  exodus = true
  execute_on = 'INITIAL FINAL'
  perf_graph = true
[]
