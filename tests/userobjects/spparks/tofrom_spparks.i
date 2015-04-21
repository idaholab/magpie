[GlobalParams]
  xmin = 0
  ymin = 0
  zmin = 0 
  xmax = 5
  ymax = 5
  zmax = 5 
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 5
  ny = 5
  nz = 5
[]

[Functions]
  [./conc]
    type = ParsedFunction
    value = 0.5
  [../]
  [./linear]
    type = ParsedFunction
    value = t+7
  [../]
[]

[AuxVariables]
  [./q_alpha]
    initial_condition = 0
  [../]
  [./spin]
    initial_condition = 50
  [../]
[]

[AuxKernels]
#  [./get_spin]
#    type = SPPARKSAux
#    variable = spin
#    ivar = 1
#    user_object = spparks
#    execute_on = timestep_begin
#  [../]
#  [./get_phase]
#    type = SPPARKSAux
#    variable = q_alpha
#    ivar = 2
#    user_object = spparks
#    execute_on = timestep_begin
#  [../]
[]

[Variables]
  [./c]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[ICs]
  [./c]
    type = FunctionIC
    function = conc
    variable = c
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = c
  [../]
  [./time]
    type = TimeDerivative
    variable = c
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'x y'
    [../]
  [../]
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

  dt = 1 #200
  start_time = 0
  end_time = 1 #80000 #1600
[]

[Outputs]
  output_initial = true
  exodus = true
  interval = 1
  [./console]
    type = Console
    output_linear = true
  [../]
[]

[UserObjects]
  [./spparks]
    type = SPPARKSUserObject
    file = in.rpv.fs
    #spparks_only = true
    one_time_run = true
    #to_dvar = 1
    #double_vars = c
    from_ivar = '1 2'
    from_dvar = 5
    sol_vars = c
    init_spparks = false
    execute_on = timestep_end
    time_spparks_time_ratio = 10
  [../]
[]
