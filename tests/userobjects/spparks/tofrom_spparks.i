[GlobalParams]
  xmin = 0
  ymin = 0
  zmin = 0
  xmax = 8
  ymax = 8
  zmax = 4
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 8
  ny = 8
  nz = 4
[]

[Functions]
  [conc]
    type = ParsedFunction
    expression = 0.5
  []
  [linear]
    type = ParsedFunction
    expression = t+7
  []
[]

[Variables]
  [c]
    order = FIRST
    family = LAGRANGE
  []
[]

[ICs]
  [c]
    type = FunctionIC
    function = conc
    variable = c
  []
[]

[Kernels]
  [diff]
    type = Diffusion
    variable = c
  []
  [time]
    type = TimeDerivative
    variable = c
  []
[]

[BCs]
  [Periodic]
    [all]
      auto_direction = 'x y'
    []
  []
[]

[Executioner]
  type = Transient

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm      101                 preonly       lu           1'

  l_max_its = 30
  l_tol = 1.0e-3

  nl_max_its = 50
  nl_rel_tol = 1e-8 #1.0e-10
  nl_abs_tol = 1e-10

  num_steps = 1
[]

[Outputs]
  exodus = true
  interval = 1
  execute_on = 'TIMESTEP_END'
  [console]
    type = Console
    output_linear = true
  []
[]

[UserObjects]
  [spparks]
    type = SPPARKSUserObject
    file = in.rpv.fs
    one_time_run = true
    from_ivar = '1 2'
    from_dvar = 5
    sol_vars = c
    init_spparks = false
    execute_on = 'timestep_end'
    time_spparks_time_ratio = 0.000001
  []
[]
