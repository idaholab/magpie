# dimensions are in Angstrom (A = 1e-10 m)
# dimensions here are 1 x 1 x 1 micro-meter
#
# _NOTE_ this calculation is actually for 3D but this is a 2D setup to save
# runtime

# This input is a good example on how to set periodicity when all concentration
# variables are AuxVariables

[Mesh]
  type = MyTRIMMesh
  dim = 3
  nx = 20
  ny = 20
  nz = 20
  xmin = 0
  xmax = 1.0e4
  ymin = 0
  ymax = 1.0e4
  zmin = 0
  zmax = 1.0e4
[]

[Problem]
  kernel_coverage_check = false
[]

[Variables]
  [./c]
    [./InitialCondition]
      type = ConstantIC
      value = 1
    [../]
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'x y z'
    [../]
  [../]
[]

[AuxVariables]
  [./c_U238]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0.6
  [../]
  [./c_O]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0.2
  [../]
  [./c_Th222]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0.02
  [../]
  [./c_He4]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0.0
  [../]
  [./int_Th222]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vac_Th222]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./int_alpha]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vac_alpha]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./int_U_aux]
    variable = int_Th222
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 1
    defect = INT
  [../]

  [./vac_U_aux]
    variable = vac_Th222
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 1
    defect = VAC
  [../]

  [./int_alpha]
    variable = int_alpha
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 3
    defect = INT
  [../]

  [./vac_alpha]
    variable = vac_alpha
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 3
    defect = INT
  [../]
[]

[UserObjects]
  [./pka_alpha_generator]
    type = PKAGeneratorAlphaDecay
    ZAID = '922380 902220 80160 20040'
    time_unit = microsecond
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = 'c_U238 c_Th222 c_O c_He4'
    M   = '238    222     16  4'
    Z   = '92     90      8   2'
    site_volume = 0.0404 # nm^3 per UO2 unit
    periodic_var = c
    pka_generator = pka_alpha_generator
  [../]
  [./runner]
    type = MyTRIMElementRun
    rasterizer = rasterizer
  [../]
[]

[Postprocessors]
  [./vac_Th222]
    type = ElementIntegralVariablePostprocessor
    variable = vac_Th222
    execute_on = timestep_end
  [../]
  [./int_Th222]
    type = ElementIntegralVariablePostprocessor
    variable = int_Th222
    execute_on = timestep_end
  [../]
  [./int_alpha]
    type = ElementIntegralVariablePostprocessor
    variable = int_alpha
    execute_on = timestep_end
  [../]
  [./vac_alpha]
    type = ElementIntegralVariablePostprocessor
    variable = vac_alpha
    execute_on = timestep_end
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
  dt = 1.0e-5
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
[]
