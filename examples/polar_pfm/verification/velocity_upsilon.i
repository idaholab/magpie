#
# HMX example from page 8 of Physical Review B 89, 184102 (2014)
# delta-phase is solid 1
# beta-phase is solid 2
#

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 500
  xmax = 100
  xmin = -100
[]

[Variables]
  [./Upsilon]
  [../]
  [./theta]
    initial_condition = 1
  [../]
[]

# i) disappearing melt
# KE = 1.933
KE = 2.2
KD = 1   # energy ratios

# ii) barrierless formation of melt
#kE = 3.39

# simulation temperature θ
Temp_ = 432 # equilibrium temperature
P4 = 2.415  # scale factor for interface width
ru0 = ${units 1848.78 kg/m^3 -> g/nm^3} # density conversion
z_mb = ${units 505.98 J/(kg*K) -> eV/(g*K)}
z_md = ${units 429.358 J/(kg*K) -> eV/(g*K)}
z_ss = ${fparse z_mb-z_md}
E_ss = ${units 1 J/m^2 -> eV/nm^2} # length conversion
E_md = ${fparse E_ss/KE}
E_mb = ${fparse E_ss/KE}
Delta_ss = ${units 3e-9 m -> nm} # length conversion
Delta_md = ${fparse Delta_ss/KD}
Delta_mb = ${fparse Delta_ss/KD}
Theta_e_ss = ${units 432 K}
Theta_e_md = ${units 550 K}
Theta_e_mb = ${units 532.14 K}
Theta_c_ss = ${fparse Theta_e_ss-(E_ss*P4)/(z_ss*Delta_ss*ru0)}
Theta_c_md = ${fparse Theta_e_md-(E_md*P4)/(z_md*Delta_md*ru0)}
Theta_c_mb = ${fparse Theta_e_mb-(E_mb*P4)/(z_mb*Delta_mb*ru0)}
A1_t = 0
A2_t = 0
Ab_t = 0
A0_md = ${fparse 3*z_md*ru0}
A1    = ${fparse A1_t+A0_md*(Temp_-Theta_c_md)}
A0_mb = ${fparse 3*z_mb*ru0}
A2    = ${fparse A2_t+A0_mb*(Temp_-Theta_c_mb)}
A0_ss = ${fparse 3*z_ss*ru0}
Ab = ${fparse Ab_t+A0_ss*(Temp_-Theta_c_ss)}
G1 = ${fparse ${units 429.39 J/(kg*K) -> eV/(g*K)}*ru0*(Temp_-Theta_e_md)}
G2 = ${fparse ${units 505.97 J/(kg*K) -> eV/(g*K)}*ru0*(Temp_-Theta_e_mb)}
# lambda_s = ${fparse ${units 2596.5 m^2/(N*s)}/2}
# lambda_m = ${units 2596.5 m^2/(N*s)}
beta1   = ${fparse 6*E_md*Delta_md/P4}
beta2   = ${fparse 6*E_mb*Delta_mb/P4}
beta_ss = ${fparse 6*E_ss*Delta_ss/P4}
a0 = 1e-4 #0.01
# x0 = 0
L_theta = ${units 1298.3 m*s/kg -> nm*s/g}
L_upsilon = ${units 2596.5 m*s/kg -> nm*s/g}

[ICs]
  [./upsilon]
    type = FunctionIC
    variable = Upsilon
    function = '1/(1 + exp(-${P4}*2*(x - 0)/${Delta_mb}))'
    #function = 'x/40'
  [../]
[]

[Functions]
  [vanalytic]
    type = ParsedFunction
    # theta = 1 -> probing L20
    value = 'dG20:=${G2}; 6*${L_upsilon}*${Delta_mb}*dG20/${P4}'
  []
[]

[Magpie]
  [./PolarPhaseField]
    [./all]
      a0       = ${a0}
      a_A      = 3
      a_theta  = 3
      a_phi    = 3
      a_beta   = 3 # "It is also assumed that all a = 3."
      beta10 = ${beta1} # m-delta  g(kE,kdelta) (nJ/m)
      beta20 = ${beta2} # m-beta (nJ/m)
      beta21 = ${beta_ss} # delta-beta (nJ/m)
      G0       = 0
      DeltaG10 = ${G1}
      DeltaG20 = ${G2}
      A10      = ${A1}
      A20      = ${A2}
      A21      = ${Ab}
      L_theta = ${L_theta}
      L_upsilon = ${L_upsilon}
    [../]
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  nl_abs_tol = 1e-12

  nl_max_its = 30
  l_max_its = 100

  # [./TimeStepper]
  #   type = IterationAdaptiveDT
  #   optimal_iterations = 6
  #   iteration_window = 2
    dt = 1e-11
  # [../]

  [./Predictor]
    type = AdamsPredictor
    scale = 0.5
  [../]

  num_steps = 200

  steady_state_detection = true
  steady_state_tolerance = 1e-08
  automatic_scaling = true
[]

[VectorPostprocessors]
  [op]
    type = LineValueSampler
    variable = 'theta Upsilon'
    start_point = '-100 0 0'
    end_point = '100 0 0'
    sort_by = x
    num_points = 501
  []
[]

[Postprocessors]
  [./pos]
    type = FindValueOnLine
    v = Upsilon
    start_point = ' -100 0 0'
    end_point   = '100 0 0'
    target = 0.5
    tol = 1e-6
    error_if_not_found = true
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./dpos]
    type = ChangeOverTimestepPostprocessor
    postprocessor = pos
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [vanalytic]
    type = FunctionValuePostprocessor
    function = vanalytic
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[UserObjects]
  [./terminate]
    type = Terminator
    expression = 'pos < 40 | pos > 160'
  [../]
[]

[Outputs]
  print_linear_residuals = false
  file_base = 'out/velocity_upsilon_${KE}_${KD}_${Temp_}'
  csv = true
  execute_on = 'TIMESTEP_END'
[]
