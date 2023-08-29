#
# Unit is Angstrom, N_Ta and N_O are number densities # / volume and the volume
# is in A^3. Total number density N = 0.048
#

[Mesh]
  type = MyTRIMMesh
  dim = 3
  nx = 5
  ny = 5
  nz = 5
  xmin = 0
  xmax = 1e5
  ymin = 0
  ymax = 1e5
  zmin = 0
  zmax = 1e5
[]

[Problem]
  kernel_coverage_check = false
[]

[Variables]
  [dummy]
  []
[]

[AuxVariables]
  [N_Ta]
    order = CONSTANT
    family = MONOMIAL
    [InitialCondition]
      type = FunctionIC
      function = numd_Ta
    []
  []

  [N_O]
    order = CONSTANT
    family = MONOMIAL
    [InitialCondition]
      type = FunctionIC
      function = numd_O
    []
  []

  [vac_Ta]
    order = CONSTANT
    family = MONOMIAL
  []

  [vac_O]
    order = CONSTANT
    family = MONOMIAL
  []

  [inter_Ta]
    order = CONSTANT
    family = MONOMIAL
  []

  [inter_O]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Functions]
  [numd_Ta]
    type = ParsedFunction
    expression = 'N := 0.048; L := 1e5; N * x / L'
  []

  [numd_O]
    type = ParsedFunction
    expression = 'N := 0.048; L := 1e5; N * (1 - x / L)'
  []
[]

[AuxKernels]
  [vac_Ta]
    variable = vac_Ta
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 0
    defect = VAC
  []

  [vac_O]
    variable = vac_O
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 1
    defect = VAC
  []

  [inter_Ta]
    variable = inter_Ta
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 0
    defect = INT
  []

  [inter_O]
    variable = inter_O
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 1
    defect = INT
  []
[]

[UserObjects]
  [constant]
    type = PKAGun
    Z = 8
    m = 16
    num_pkas = 200
    E = 5e5
    point = '0.001 5e4 5e4'
    direction = '1 0 0'
  []

  [rasterizer]
    type = MyTRIMRasterizer
    var = 'N_Ta N_O'
    M = '181 16'
    Z = '73  8'
    Ebind = '3  5'
    Edisp = '60 40'
    site_volume = 1
    pka_generator = constant

    # control NRT
    max_nrt_difference = 0.4
    analytical_energy_cutoff = 500
  []

  [runner]
    type = MyTRIMElementRun
    rasterizer = rasterizer
  []
[]

[Postprocessors]
  [integral_vac_Ta]
    type = ElementIntegralVariablePostprocessor
    variable = vac_Ta
  []

  [integral_vac_O]
    type = ElementIntegralVariablePostprocessor
    variable = vac_O
  []

  [integral_inter_Ta]
    type = ElementIntegralVariablePostprocessor
    variable = inter_Ta
  []

  [integral_inter_O]
    type = ElementIntegralVariablePostprocessor
    variable = inter_O
  []
[]

[Executioner]
  type = Transient
  num_steps = 1
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
[]
