[Mesh]
  type = MyTRIMMesh
  dim = 2
  xmax = 100
  ymin = -25
  ymax = 25
  nx = 100
  ny = 50
[]

[Problem]
  type = FFTProblem
[]

[AuxVariables]
  [./c_fft]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = FunctionIC
      function = 'cos(x/100*2*pi*4)*cos(y/100*2*pi*3)'
    [../]
  [../]

  [./R0_fft]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = FunctionIC
      function = 'cos(x/100*2*pi*2)*cos(y/100*2*pi*4)'
    [../]
  [../]
  [./R1_fft]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = FunctionIC
      function = 'cos(x/100*2*pi*3)*cos(y/100*2*pi*3)'
    [../]
  [../]
  [./R2_fft]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = FunctionIC
      function = 'cos(x/100*2*pi*4)*cos(y/100*2*pi*2)'
    [../]
  [../]


  [./u_fft]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = SmoothSuperellipsoidIC
      x1 = 50
      y1 = 0
      a = 30
      b = 15
      int_width = 15
      invalue = 1
      outvalue = 0
      n = 2
    [../]
  [../]
  [./v_fft]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = FunctionIC
      function = if(y<0,y,-y)
    [../]
  [../]
  [./grad_u0_fft]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./grad_u1_fft]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./grad_v1_fft]
    order = CONSTANT
    family = MONOMIAL
  [../]
  
  [./u_fem]
    [./InitialCondition]
      type = SmoothSuperellipsoidIC
      x1 = 50
      y1 = 0
      a = 30
      b = 15
      int_width = 15
      invalue = 1
      outvalue = 0
      n = 2
    [../]
  [../]
  [./v_fem]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = FunctionIC
      function = if(y<0,y,-y)
    [../]
  [../]
  [./grad_u0_fem]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./grad_u1_fem]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./grad_v1_fem]
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
    moose_variable = c_fft
  [../]
  [./R]
    type = RealVectorValueFFTWBuffer
    moose_variable = 'R0_fft R1_fft R2_fft'
  [../]

  [./u]
    type = RealFFTWBuffer
    moose_variable = u_fft
  [../]
  [./grad_u]
    type = RealVectorValueFFTWBuffer
  [../]
  [./v]
    type = RealFFTWBuffer
    moose_variable = v_fft
  [../]
  [./grad_v]
    type = RealVectorValueFFTWBuffer
  [../]
[]

[AuxKernels]
  [./c_fft]
    type = FFTBufferAux
    variable = c_fft
    fft_buffer = c
    execute_on = FINAL
  [../]

  [./R0_fft]
    type = FFTBufferAux
    variable = R0_fft
    fft_buffer = R
    component = 0
    execute_on = FINAL
  [../]
  [./R1_fft]
    type = FFTBufferAux
    variable = R1_fft
    fft_buffer = R
    component = 1
    execute_on = FINAL
  [../]
  [./R2_fft]
    type = FFTBufferAux
    variable = R2_fft
    fft_buffer = R
    component = 2
    execute_on = FINAL
  [../]

  [./u_fft]
    type = FFTBufferAux
    variable = u_fft
    fft_buffer = u
    execute_on = FINAL
  [../]
  [./v_fft]
    type = FFTBufferAux
    variable = v_fft
    fft_buffer = v
    execute_on = FINAL
  [../]

  [./grad_u0_fft]
    type = FFTBufferAux
    variable = grad_u0_fft
    fft_buffer = grad_u
    component = 0
    execute_on = FINAL
  [../]
  [./grad_u1_fft]
    type = FFTBufferAux
    variable = grad_u1_fft
    fft_buffer = grad_u
    component = 1
    execute_on = FINAL
  [../]
  [./grad_v1_fft]
    type = FFTBufferAux
    variable = grad_v1_fft
    fft_buffer = grad_v
    component = 1
    execute_on = FINAL
  [../]

  [./grad_u0_fem]
    type = GradientComponentAux
    variable = grad_u0_fem
    v = u_fem
    component = 0
    execute_on = FINAL
  [../]
  [./grad_u1_fem]
    type = GradientComponentAux
    variable = grad_u1_fem
    v = u_fem
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
