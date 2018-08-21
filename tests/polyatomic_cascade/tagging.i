[Mesh]
  type = MyTRIMMesh
  dim = 3
  nx = 2
  ny = 2
  nz = 2
  xmin = 0
  xmax = 100
  ymin = 0
  ymax = 100
  zmin = 0
  zmax = 100
[]

[Problem]
  kernel_coverage_check = false
[]

[Variables]
  [./dummy]
  [../]
[]

[AuxVariables]
  [./cC]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0.5
  [../]
  [./cSi]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0.5
  [../]
  [./cXe]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0.0
  [../]

  [./int]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./int]
    variable = int
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 2
    defect = INT
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'x y z'
    [../]
  [../]
[]

[UserObjects]
  [./ffe]
    type = PKAFissionFragmentEmpirical
    fission_rate = 0.00001
    relative_density = 1
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = 'cC cSi cXe'
    M =   '12 28  136'
    Z =   '6  14  54'
    # MTol = '0.5 0.5 4'
    site_volume = 0.0404
    pka_generator = ffe
    periodic_var = dummy
    analytical_energy_cutoff = 1e5
  [../]
  [./runner]
    type = MyTRIMElementRun
    rasterizer = rasterizer
  [../]
[]

[Postprocessors]
  [./total_vacancies]
    type = ElementIntegralVariablePostprocessor
    variable = int
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
[]
