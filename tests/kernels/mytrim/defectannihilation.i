[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  xmin = 0
  xmax = 100
  ymin = 0
  ymax = 100
[]

[Variables]
  [./c_vac]
    [./InitialCondition]
      type = FunctionIC
      function = sin(x/100*2*pi)/2+0.5
    [../]
  [../]
  [./c_int]
    [./InitialCondition]
      type = FunctionIC
      function = sin(y/100*2*pi)/2+0.5
    [../]
  [../]
[]

[Kernels]
  [./vac_dt]
    type = TimeDerivative
    variable = c_vac
  [../]
  [./vac_react]
    type = DefectAnnihilation
    variable = c_vac
    v = c_int
  [../]
  [./vac_diff]
    type = Diffusion
    variable = c_vac
  [../]

  [./int_dt]
    type = TimeDerivative
    variable = c_int
  [../]
  [./int_react]
    type = DefectAnnihilation
    variable = c_int
    v = c_vac
  [../]
  [./int_diff]
    type = Diffusion
    variable = c_int
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      variable = 'c_int c_vac'
      auto_direction = 'x y'
    [../]
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    coupled_groups = 'c_vac,c_int'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  num_steps = 1
  nl_abs_tol = 1e-10
  dt = 0.5
[]

[Outputs]
  exodus = true
  print_linear_residuals = false
[]
