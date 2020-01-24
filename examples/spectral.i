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

  [./R0_aux]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = FunctionIC
      function = 'cos(x/100*2*pi*2)*cos(y/100*2*pi*4)'
    [../]
  [../]
  [./R1_aux]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = FunctionIC
      function = 'cos(x/100*2*pi*3)*cos(y/100*2*pi*3)'
    [../]
  [../]
  [./R2_aux]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = FunctionIC
      function = 'cos(x/100*2*pi*4)*cos(y/100*2*pi*2)'
    [../]
  [../]


  [./u_aux]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = SmoothCircleIC
      x1 = 50
      y1 = 50
      radius = 30
      int_width = 20
      invalue = 1
      outvalue = 0
    [../]
  [../]
  [./grad_u0_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./grad_u1_aux]
    order = CONSTANT
    family = MONOMIAL
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
    type = RealVectorValueFFTWBuffer
    moose_variable = 'R0_aux R1_aux R2_aux'
  [../]

  [./u]
    type = RealFFTWBuffer
    moose_variable = u_aux
  [../]
  [./grad_u]
    type = RealVectorValueFFTWBuffer
  [../]
[]

[AuxKernels]
  [./c_aux]
    type = FFTBufferAux
    variable = c_aux
    fft_buffer = c
    execute_on = FINAL
  [../]

  [./R0_aux]
    type = FFTBufferAux
    variable = R0_aux
    fft_buffer = R
    component = 0
    execute_on = FINAL
  [../]
  [./R1_aux]
    type = FFTBufferAux
    variable = R1_aux
    fft_buffer = R
    component = 1
    execute_on = FINAL
  [../]
  [./R2_aux]
    type = FFTBufferAux
    variable = R2_aux
    fft_buffer = R
    component = 2
    execute_on = FINAL
  [../]

  [./u_aux]
    type = FFTBufferAux
    variable = u_aux
    fft_buffer = u
    execute_on = FINAL
  [../]

  [./grad_u0_aux]
    type = FFTBufferAux
    variable = grad_u0_aux
    fft_buffer = grad_u
    component = 0
    execute_on = FINAL
  [../]
  [./grad_u1_aux]
    type = FFTBufferAux
    variable = grad_u1_aux
    fft_buffer = grad_u
    component = 1
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
