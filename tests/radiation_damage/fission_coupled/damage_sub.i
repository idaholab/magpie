# dimensions are in nm
# typical treat fuel particle is rf = 20 mu-m = 2e4 nm
# typical distance between fuel particles is D = 0.893 * rf / cube_root(1/2571)
# D = 244 mu-m
# the typical fuel particle distance is equal to the domain size
#
# _NOTE_ this calculation is actually for 3D but this is a 2D setup to save
# runtime
[Mesh]
  type = MyTRIMMesh
  dim = 2
  nx = 100
  ny = 100
  xmin = -1.22e5
  xmax = 1.22e5
  ymin = -1.22e5
  ymax = 1.22e5
  elem_type = QUAD4
  uniform_refine = 0
[]

[Variables]
  [./c]
    initial_condition = 1.0
  [../]
[]

[AuxVariables]
  [./int]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vac]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./c_U]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = FunctionIC
      function = func_c_U
    [../]
  [../]
  [./c_O]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = FunctionIC
      function = func_c_O
    [../]
  [../]
  [./c_C]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = FunctionIC
      function = func_c_C
    [../]
  [../]
  [./pka_count]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./dt]
    type = TimeDerivative
    variable = c
  [../]
  [./diff]
    type = Diffusion
    variable = c
  [../]
[]

[AuxKernels]
  [./int]
    variable = int
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 0
    defect = INT
  [../]
  [./vac]
    variable = vac
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 0
    defect = VAC
  [../]
  [./count_pkas]
    type = PKAAux
    variable = pka_count
    rasterizer = rasterizer
  [../]
[]

[UserObjects]
  [./neutronics_fission_generator]
    type = PKAFissionFragmentNeutronics
    relative_density = 1
    fission_rate = 3.2498e-6
    #fission_rate = 3.25e-6
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = 'c_U  c_O  c_C'
    M   = '235  16    12'
    Z   = '92   8     6'
    site_volume = 0.0404 # nm^3 per UO2 unit
    pka_generator = neutronics_fission_generator
  [../]
  [./runner]
    type = MyTRIMElementRun
    rasterizer = rasterizer
  [../]
[]

[Postprocessors]
  [./int]
    type = ElementIntegralVariablePostprocessor
    variable = int
    execute_on = timestep_end
  [../]
  [./vac]
    type = ElementIntegralVariablePostprocessor
    variable = vac
    execute_on = timestep_end
  [../]
  [./n_pkas]
    type = ElementIntegralVariablePostprocessor
    variable = pka_count
    execute_on = timestep_end
  [../]
[]

[Functions]
  [./func_c_U]
    type = ParsedFunction
    value = 'r := sqrt(x*x+y*y); rho := if(r/rf<1.0e-10,1.0e-10,r/rf); th := tanh(a*(rho-1/rho)); 0.5 * (1+th) * ev + 0.5 * (1-th) * iv'
    vars = 'rf    a ev iv'
    vals = '2.0e4 2  0 0.3333333'
  [../]
  [./func_c_O]
    type = ParsedFunction
    value = 'r := sqrt(x*x+y*y); rho := if(r/rf<1.0e-10,1.0e-10,r/rf); th := tanh(a*(rho-1/rho)); 0.5 * (1+th) * ev + 0.5 * (1-th) * iv'
    vars = 'rf    a ev iv'
    vals = '2.0e4 2  0 0.66666666'
  [../]
  [./func_c_C]
    # ev computed as follows:
    # rho_G = 1.72 g/cc (see Kun MO's paper)
    # rho_UO2 = 10.963 g/cc
    # MUO2 = 235 + 2*16 = 267
    # MG = 12
    # xG = rho_G / rho_UO2 * MUO2 / MG = 3.4908
    type = ParsedFunction
    value = 'r := sqrt(x*x+y*y); rho := if(r/rf<1.0e-10,1.0e-10,r/rf); th := tanh(a*(rho-1/rho)); 0.5 * (1+th) * ev + 0.5 * (1-th) * iv'
    vars = 'rf    a ev     iv'
    vals = '2.0e4 2 3.4908 0'
  [../]
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
