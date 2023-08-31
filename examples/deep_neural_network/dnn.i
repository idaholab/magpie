[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 100
  xmax = 70
  ymax = 70
[]

[AuxVariables]
  [T]
    initial_condition = 0.2
  []
  [diff]
  []
[]

[AuxKernels]
  [diff]
    type = ParsedAux
    variable = diff
    function = 'cn - cp'
    coupled_variables = 'cn cp'
  []
[]

[Variables]
  # active = 'cp wp'
  [cn]
    [InitialCondition]
      type = RandomIC
      min = 0.45
      max = 0.55
      seed = 1234
      #type = FunctionIC
      #function = cos(x/70*2*pi)/2+0.5
    []
  []
  [wn]
  []

  [cp]
    [InitialCondition]
      type = RandomIC
      min = 0.45
      max = 0.55
      seed = 1234
      #type = FunctionIC
      #function = cos(x/70*2*pi)/2+0.5
    []
  []
  [wp]
  []
[]

[Kernels]
  # active = 'cp_res wp_res timep'
  [cn_res]
    type = ADSplitCHParsed
    variable = cn
    property_name = Fn
    kappa_name = ad_kappa_c
    w = wn
  []
  [wn_res]
    type = ADSplitCHWRes
    variable = wn
    mob_name = ad_M
  []
  [timen]
    type = ADCoupledTimeDerivative
    variable = wn
    v = cn
  []

  [cp_res]
    type = SplitCHParsed
    variable = cp
    property_name = Fp
    kappa_name = kappa_c
    w = wp
  []
  [wp_res]
    type = SplitCHWRes
    variable = wp
    mob_name = M
  []
  [timep]
    type = CoupledTimeDerivative
    variable = wp
    v = cp
  []
[]

[BCs]
  [Periodic]
    [all]
      auto_direction = 'x y'
    []
  []
[]

[Materials]
  # active = 'Fp pfmobility'
  [Fn]
    type = NeuralNetFreeEnergy
    file_name = weights_biases.txt
    file_format = MAGPIE
    activation_function = SOFTSIGN
    inputs = 'T cn'
    prop_names = 'Fn'
    outputs = exodus
    debug = true
  []
  [Fp]
    type = DerivativeParsedMaterial
    property_name = 'Fp'
    function = 'cp*(1-cp)+0.001*(T*300+400)*(cp*log(cp)+(1-cp)*log(1-cp))'
    coupled_variables = 'cp T'
    derivative_order = 2
    outputs = exodus
  []
  [pfmobility]
    type = GenericConstantMaterial
    prop_names = 'M kappa_c'
    prop_values = '1 0.25'
  []
  [adpfmobility]
    type = ADGenericConstantMaterial
    prop_names = 'ad_M ad_kappa_c'
    prop_values = '1 0.25'
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Problem]
  kernel_coverage_check = false
[]

[Postprocessors]
  [diff]
    type = ElementL2Difference
    variable = cn
    other_variable = cp
  []
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
  gnuplot = true
  print_linear_residuals = false
[]
