#
# HMX example from page 8 of Physical Review B 89, 184102 (2014)
# delta-phase is solid 1
# beta-phase is solid 2
#

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 200
  xmax = 40
[]

[Variables]
  [./Upsilon]
    initial_condition = 1
  [../]
  [./theta]
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
a0 = 0.01
x0 = 20

[Functions]
  [theta]
    type = ParsedFunction
    value = '1/(1 + exp(-${P4}*(x - ${x0})/${Delta_ss}))'
  []
[]

[ICs]
  [./theta]
    type = FunctionIC
    variable = theta
    #function = '1/(1 + exp(-${P4}*2*(x - ${x0})/${Delta_ss}))'
    function = 'x/40'
  [../]
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
  nl_abs_tol = 1e-12

  nl_max_its = 30
  l_max_its = 100

  [./TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 8
    iteration_window = 2
    dt = 1e-7
  [../]

  steady_state_detection = true
  steady_state_tolerance = 1e-06
  automatic_scaling = true
[]

[VectorPostprocessors]
  [./sol]
    type = LineValueSampler
    variable = 'Upsilon theta'
    start_point = ' 0 0 0'
    end_point   = '40 0 0'
    sort_by = x
    num_points = 801
  [../]
  [./ana]
    type = LineFunctionSampler
    functions = theta
    start_point = ' 0 0 0'
    end_point   = '40 0 0'
    sort_by = x
    num_points = 801
  [../]
[]

[Outputs]
  print_linear_residuals = false
  file_base = 'out/steady_theta_${KE}_${KD}_${Temp_}'
  csv = true
[]
