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
  [./u_aux]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = SmoothCircleIC
      x1 = 50
      y1 = 50
      radius = 10
      int_width = 1
      invalue = 1
      outvalue = 0.0001
    [../]
  [../]
[]

[UserObjects]
  # Buffers
  [./u]
    type = RealFFTWBuffer
    moose_variable = u_aux
  [../]
  [./grad_u]
    type = RealVectorValueFFTWBuffer
  [../]
[]

[AuxKernels]
  [./u_aux]
    type = FFTBufferAux
    variable = u_aux
    fft_buffer = u
    execute_on = final
  [../]
[]

[Executioner]
  type = SpectralExecutionerDiffusion

  # Parameters for diffusion executioner
  diffusion_coefficient = 400.0 # 1000.0
  time_step = 0.0005
  number_steps = 200
[]

[VectorPostprocessors]
  [./linevalue]
    type = LineValueSampler
    variable = 'u_aux'
    start_point = '0 0 0'
    end_point = '99.999999999 0 0'
    num_points = 101
    sort_by = x
    execute_on = final
  [../]
[]

[Outputs]
  exodus = true
  execute_on = 'INITIAL FINAL'
  perf_graph = true
  [./comp]
    type = CSV
  [../]
[]
