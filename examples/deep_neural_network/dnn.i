[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 100
  xmax = 70
  ymax = 70
[]

[AuxVariables]
  [./T]
    initial_condition = 0.2
  [../]
  [./diff]
  [../]
[]

[AuxKernels]
  [./diff]
    type = ParsedAux
    variable = diff
    function = 'cn - cp'
    args = 'cn cp'
  [../]
[]

[Variables]
  # active = 'cp wp'
  [./cn]
    [./InitialCondition]
      type = RandomIC
      min = 0.45
      max = 0.55
      seed = 1234
    [../]
  [../]
  [./wn]
  [../]

  [./cp]
    [./InitialCondition]
      type = RandomIC
      min = 0.45
      max = 0.55
      seed = 1234
    [../]
  [../]
  [./wp]
  [../]
[]

[Kernels]
  # active = 'cp_res wp_res timep'
  [./cn_res]
    type = ADSplitCHParsed
    variable = cn
    f_name = Fn
    kappa_name = kappa_c
    w = wn
  [../]
  [./wn_res]
    type = ADSplitCHWRes
    variable = wn
    mob_name = M
  [../]
  [./timen]
    type = ADCoupledTimeDerivative
    variable = wn
    v = cn
  [../]

  [./cp_res]
    type = SplitCHParsed
    variable = cp
    f_name = Fp
    kappa_name = kappa_c
    w = wp
  [../]
  [./wp_res]
    type = SplitCHWRes
    variable = wp
    mob_name = M
  [../]
  [./timep]
    type = CoupledTimeDerivative
    variable = wp
    v = cp
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'x y'
    [../]
  [../]
[]

[Materials]
  # active = 'Fp pfmobility'
  [./Fn]
    type = DeepNeuralNetFreeEnergy
    filename = weights_biases.txt
    inputs     = 'T cn'
    prop_names = 'Fn'
  [../]
  [./Fp]
    type = DerivativeParsedMaterial
    f_name = 'Fp'
    function = 'cp*(1-cp)+0.001*(T*300+400)*(cp*log(cp)+(1-cp)*log(1-cp))'
    args = 'cp T'
  [../]
  [./pfmobility]
    type = GenericConstantMaterial
    prop_names  = 'M kappa_c'
    prop_values = '1 0.25'
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Problem]
  kernel_coverage_check = false
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'

  nl_abs_tol = 1e-12
  num_steps = 150
  dt = 5
[]

[Outputs]
  exodus = true
  print_linear_residuals = false
[]
