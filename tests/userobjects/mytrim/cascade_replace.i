[Mesh]
  type = MyTRIMMesh
  dim = 2
  nx = 20
  ny = 20
  xmin = -100
  xmax = 100
  ymin = 0
  ymax = 200
[]

[Variables]
  [./c0]
    initial_condition = 0.5
  [../]
  [./c1]
    initial_condition = 0.5
  [../]
[]

[AuxVariables]
  [./rep0_in]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rep0_out]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rep1_in]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rep1_out]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Problem]
  kernel_coverage_check = false
[]

[AuxKernels]
  [./rep0_in]
    variable = rep0_in
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 0
    defect = REPLACEMENT_IN
  [../]
  [./rep0_out]
    variable = rep0_out
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 0
    defect = REPLACEMENT_OUT
  [../]
  [./rep1_in]
    variable = rep1_in
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 1
    defect = REPLACEMENT_IN
  [../]
  [./rep1_out]
    variable = rep1_out
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 1
    defect = REPLACEMENT_OUT
  [../]
[]

[UserObjects]
  [./thermal_fission]
    type = PKAFissionFragmentEmpirical
    relative_density = 1
    fission_rate = 0.01
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = 'c0 c1'
    M =   '40 40'
    Z =   '20 20'
    site_volume = 0.0404 # nm^3 per UO2 unit
    pka_generator = thermal_fission
  [../]
  [./runner]
    type = MyTRIMElementRun
    rasterizer = rasterizer
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  hide = 'c0 c1'
[]
