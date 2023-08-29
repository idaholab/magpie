# dimensions are in Angstrom (A = 1e-10 m)
# typical treat fuel particle is rf = 20 mu-m = 2e5 A
# typical distance between fuel particles is D = 0.893 * rf / cube_root(1/2571)
# D = 244 mu-m = 2.44e6 A.
# the typical fuel particle distance is much larger than the mfp so we cut the domain
# at 4 * rf = 80 mu-m = 8e5 A
#
# _NOTE_ this calculation is actually for 3D but this is a 2D setup to save
# runtime
[Mesh]
  type = MyTRIMMesh
  dim = 2
  nx = 200
  ny = 200
  xmin = -8.0e5
  xmax = 8.0e5
  ymin = -8.0e5
  ymax = 8.0e5
  elem_type = QUAD4
  uniform_refine = 0
[]

[Problem]
  kernel_coverage_check = false
[]

[Variables]
  [c]
    [InitialCondition]
      type = ConstantIC
      value = 1
    []
  []
[]

[AuxVariables]
  [int_U]
    order = CONSTANT
    family = MONOMIAL
  []
  [vac_U]
    order = CONSTANT
    family = MONOMIAL
  []
  [int_O]
    order = CONSTANT
    family = MONOMIAL
  []
  [vac_O]
    order = CONSTANT
    family = MONOMIAL
  []
  [int_C]
    order = CONSTANT
    family = MONOMIAL
  []
  [vac_C]
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
  [pka_count]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [int_U_aux]
    variable = int_U
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 0
    defect = INT
  []
  [vac_U_aux]
    variable = vac_U
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 0
    defect = VAC
  []
  [int_O_aux]
    variable = int_O
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 1
    defect = INT
  []
  [vac_O_aux]
    variable = vac_O
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 1
    defect = VAC
  []
  [int_C_aux]
    variable = int_C
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 2
    defect = INT
  []
  [vac_C_aux]
    variable = vac_C
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 2
    defect = VAC
  []
[]

[UserObjects]
  [neutronics_fission_generator]
    type = PKAFissionFragmentNeutronics
    partial_reaction_rates = '2.0e-11 0 0'
  []
  [rasterizer]
    type = MyTRIMRasterizer
    var = 'c_U  c_O  c_C'
    M = '235  16    12'
    Z = '92   8     6'
    site_volume = 0.0404 # nm^3 per UO2 unit
    pka_generator = neutronics_fission_generator
  []
  [runner]
    type = MyTRIMElementRun
    rasterizer = rasterizer
  []
[]

[Postprocessors]
  #[./pka_info]
  #  type = MyTRIMPKAInfo
  #  rasterizer = rasterizer
  #  value_type = TOTAL_ENERGY
  #[../]
  [int_U_pp]
    type = ElementIntegralVariablePostprocessor
    variable = int_U
    execute_on = timestep_end
  []
  [vac_U_pp]
    type = ElementIntegralVariablePostprocessor
    variable = vac_U
    execute_on = timestep_end
  []
  [int_O_pp]
    type = ElementIntegralVariablePostprocessor
    variable = int_O
    execute_on = timestep_end
  []
  [vac_O_pp]
    type = ElementIntegralVariablePostprocessor
    variable = vac_O
    execute_on = timestep_end
  []
  [int_C_pp]
    type = ElementIntegralVariablePostprocessor
    variable = int_C
    execute_on = timestep_end
  []
  [vac_C_pp]
    type = ElementIntegralVariablePostprocessor
    variable = vac_C
    execute_on = timestep_end
  []
  [n_pkas]
    type = ElementIntegralVariablePostprocessor
    variable = pka_count
    execute_on = timestep_end
  []
[]

[Functions]
  [func_c_U]
    type = ParsedFunction
    expression = 'r := sqrt(x*x+y*y); rho := if(r/rf<1.0e-10,1.0e-10,r/rf); th := tanh(a*(rho-1/rho)); 0.5 * (1+th) * ev + 0.5 * (1-th) * iv'
    symbol_names = 'rf    a ev iv'
    symbol_values = '2.0e5 8  0 0.3333333'
  []
  [func_c_O]
    type = ParsedFunction
    expression = 'r := sqrt(x*x+y*y); rho := if(r/rf<1.0e-10,1.0e-10,r/rf); th := tanh(a*(rho-1/rho)); 0.5 * (1+th) * ev + 0.5 * (1-th) * iv'
    symbol_names = 'rf    a ev iv'
    symbol_values = '2.0e5 8  0 0.66666666'
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
    symbol_values = '2.0e5 8 3.4908 0'
  []
[]

[Executioner]
  type = Transient
  num_steps = 1
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv = true
[]
