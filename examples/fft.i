[Mesh]
  type = MyTRIMMesh
  dim = 2
  xmax = 100
  ymax = 100
  nx = 100
  ny = 100
[]

[Variables]
  [./c]
    [./InitialCondition]
      type = RandomIC
      min = 0.49
      max = 0.51
    [../]
  [../]
  [./w]
  [../]
[]

[Kernels]
  [./cres]
    type = SplitCHParsed
    variable = c
    f_name = fbulk
    kappa_name = kappa_c
    w = w
  [../]
  [./wres]
    type = SplitCHWRes
    variable = w
    mob_name = M
  [../]
  [./time]
    type = CoupledTimeDerivative
    variable = w
    v = c
  [../]
[]

[Materials]
  [./const]
    type = GenericConstantMaterial
    prop_names =  'kappa_c M'
    prop_values = '8       1'
  [../]
  [./F]
    type = DerivativeParsedMaterial
    args = c
    f_name = fbulk
    constant_names       = 'A B'
    constant_expressions = '1 0.5'
    function = '-A * (2*c-1)^2 + B * (2*c-1)^4'
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'x y'
    [../]
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 100
  l_max_its = 100
  nl_max_its = 15

  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 7
    iteration_window = 2
    growth_factor = 1.5
    cutback_factor = 0.8
    dt = 0.1
  [../]
[]

[UserObjects]
  [./fft]
    type = FourierTransform
    variable = c
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
[]

[VectorPostprocessors]
  [./power]
    type = FourierPowerSpectrum
    fourier_transform = fft
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
[]

[Postprocessors]
  [./length]
    type = FourierLengthScale
    fourier_transform = fft
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
  [./perf]
    type = PerfGraphOutput
    level = 2
  [../]
[]
