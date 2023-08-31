#
# Demonstarte how to use adaptive mesh coarsening to generate an output showing
# defect tracks at a desired detail level withthe least amount of mesh elements
#
[Mesh]
  type = MyTRIMMesh
  dim = 2
  nx = 8
  ny = 8
  xmin = -2e5
  xmax = 2e5
  ymin = -2e5
  ymax = 2e5
  zmin = -2e5
  zmax = 2e5
  uniform_refine = 4
[]

[Variables]
  [v]
    initial_condition = 0.0
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxVariables]
  [int]
    order = CONSTANT
    family = MONOMIAL
  []
  [vac]
    order = CONSTANT
    family = MONOMIAL
  []
  [c_U]
    order = CONSTANT
    family = MONOMIAL
    [InitialCondition]
      type = FunctionIC
      function = func_c_U
    []
  []
  [c_O]
    order = CONSTANT
    family = MONOMIAL
    [InitialCondition]
      type = FunctionIC
      function = func_c_O
    []
  []
  [c_C]
    order = CONSTANT
    family = MONOMIAL
    [InitialCondition]
      type = FunctionIC
      function = func_c_C
    []
  []
  [c]
    initial_condition = 1.0
  []
[]

[Kernels]
  [dt]
    type = TimeDerivative
    variable = v
  []
  [diff]
    type = Diffusion
    variable = v
  []
  [src]
    type = MyTRIMElementSource
    variable = v
    runner = runner
    ivar = 2
    defect = VAC
  []
[]

[AuxKernels]
  [int]
    variable = int
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 2
    defect = INT
  []
  [vac]
    variable = vac
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 2
    defect = VAC
  []
[]

[UserObjects]
  [fision]
    type = PKAFissionFragmentEmpirical
    relative_density = 1
    fission_rate = 4.0e-11
  []
  [rasterizer]
    type = MyTRIMRasterizer
    var = 'c_U  c_O  c_C'
    M = '235  16    12'
    Z = '92   8     6'
    site_volume = 0.0404 # nm^3 per UO2 unit
    pka_generator = fision
    print_pka_statistics = true
    interval = 2
  []
  [runner]
    type = MyTRIMElementRun
    rasterizer = rasterizer
  []
[]

[Adaptivity]
  [Markers]
    [marker]
      type = ValueThresholdMarker
      variable = v
      coarsen = 1e-9
      refine = 1e9
    []
  []

  marker = marker
  max_h_level = 4
  start_time = 0.1
  cycles_per_step = 4
[]

[Functions]
  [func_c_U]
    type = ParsedFunction
    expression = 'r := sqrt(x*x+y*y); rho := if(r/rf<1.0e-10,1.0e-10,r/rf); th := tanh(a*(rho-1/rho)); 0.5 * (1+th) * ev + 0.5 * (1-th) * iv'
    symbol_names = 'rf    a ev iv'
    symbol_values = '2.0e4 2  0 0.3333333'
  []
  [func_c_O]
    type = ParsedFunction
    expression = 'r := sqrt(x*x+y*y); rho := if(r/rf<1.0e-10,1.0e-10,r/rf); th := tanh(a*(rho-1/rho)); 0.5 * (1+th) * ev + 0.5 * (1-th) * iv'
    symbol_names = 'rf    a ev iv'
    symbol_values = '2.0e4 2  0 0.66666666'
  []
  [func_c_C]
    # ev computed as follows:
    # rho_G = 1.72 g/cc (see Kun MO's paper)
    # rho_UO2 = 10.963 g/cc
    # MUO2 = 235 + 2*16 = 267
    # MG = 12
    # xG = rho_G / rho_UO2 * MUO2 / MG = 3.4908
    type = ParsedFunction
    expression = 'r := sqrt(x*x+y*y); rho := if(r/rf<1.0e-10,1.0e-10,r/rf); th := tanh(a*(rho-1/rho)); 0.5 * (1+th) * ev + 0.5 * (1-th) * iv'
    symbol_names = 'rf    a ev     iv'
    symbol_values = '2.0e4 2 3.4908 0'
  []
[]

[Executioner]
  type = Transient
  num_steps = 2
  nl_abs_tol = 1e-10

  [TimeStepper]
    type = TimeSequenceStepper
    time_sequence = '1 1.0001'
  []
[]

[Outputs]
  exodus = true
  csv = true
[]
