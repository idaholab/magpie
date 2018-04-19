# Unit of length = Angstrom
[Mesh]
  type = MyTRIMMesh
  dim = 3
  nx = 8
  ny = 8
  nz = 8
  xmin = 0
  xmax = 1e3
  ymin = 0
  ymax = 1e3
  zmin = 0
  zmax = 1e3
[]

[Problem]
  kernel_coverage_check = false
[]

[Variables]
  [./dummy]
  [../]
[]

[AuxVariables]
  [./Cu]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'z'
    [../]
  [../]
[]

[UserObjects]
  [./constant]
    type = PKAGun
    E = 1.0e4
    Z = 29
    m = 63.546
    point = '0 500 500'
    direction = '0.6 0 0.8'
    num_pkas = 500
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = 'Cu'
    M = '63.546'
    Z = '29'
    Edisp = '33'
    site_volume = 0.0118
    pka_generator = constant
    periodic_var = dummy
  [../]
  [./runner]
    type = MyTRIMDiracRun
    rasterizer = rasterizer
  [../]
[]

[VectorPostprocessors]
  [./vacancy_sites]
    type = MyTRIMDiracResult
    runner = runner
    ivar = 0
    defect = VAC
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
  nl_abs_tol = 1e-10
[]

[Outputs]
  csv = true
[]
