[Mesh]
  type = MyTRIMMesh
  dim = 2
  nx = 30
  ny = 20
[]

[Variables]
  [./c]
    [./InitialCondition]
      type = FunctionIC
      function = 'sin(x*3*(2*pi))*sin(y*2*(2*pi))'
    [../]
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 0
[]

[Problem]
  kernel_coverage_check = false
[]

[UserObjects]
  [./fft]
    type = FourierTransform
    variable = c
    execute_on = 'INITIAL'
  [../]
[]

[VectorPostprocessors]
  [./power]
    type = FourierPowerSpectrum
    fourier_transform = fft
    execute_on = 'INITIAL'
  [../]
[]

[Outputs]
  csv = true
  exodus = true
  execute_on = INITIAL
[]
