#
# HMX example from page 8 of Physical Review B 89, 184102 (2014)
# delta-phase is solid 1
# beta-phase is solid 2
#

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 200
  xmax = 200
[]

[Variables]
  [Upsilon]
  []
  [theta]
  []
[]

[Kernels]
  [dtu]
    # general form would be SusceptibilityTimeDerivative with 1/L^Upsilon
    type = TimeDerivative
    variable = Upsilon
  []
  [dLu]
    type = PolarPFMDerivative
    variable = Upsilon
    F = psiL
  []
  [gradu]
    type = PolarPFMGradient
    variable = Upsilon
    v = theta
    F = beta21phi
  []
  [diffu] # 0 if Upsilon is constant
    type = MatDiffusion
    variable = Upsilon
    diffusivity = betaS0
    args = 'theta'
  []

  [dtt]
    # general form would be SusceptibilityTimeDerivative with 1/L^theta
    type = TimeDerivative
    variable = theta
  []
  [dLt]
    type = PolarPFMDerivative
    variable = theta
    F = psiL
  []
  [gradt] # 0 if Upsilon is constant
    type = PolarPFMGradient
    variable = theta
    v = Upsilon
    F = betaS0
  []
  [difft]
    type = MatDiffusion
    variable = theta
    diffusivity = beta21phi
    args = 'Upsilon'
  []
[]

# i) disappearing melt
# KE = 1.933
KE = 4
KD = 1

# ii) barrierless formation of melt
#kE = 3.39

# simulation temperature θ
Temp_ = 535
P4 = 2.415 # scale factor for interface width
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
A1 = ${fparse A1_t+A0_md*(Temp_-Theta_c_md)}
A0_mb = ${fparse 3*z_mb*ru0}
A2 = ${fparse A2_t+A0_mb*(Temp_-Theta_c_mb)}
A0_ss = ${fparse 3*z_ss*ru0}
Ab = ${fparse Ab_t+A0_ss*(Temp_-Theta_c_ss)}
G1 = ${fparse ${units 429.39 J/(kg*K) -> eV/(g*K)}*ru0*(Temp_-Theta_e_md)}
G2 = ${fparse ${units 505.97 J/(kg*K) -> eV/(g*K)}*ru0*(Temp_-Theta_e_mb)}
beta1 = ${fparse 6*E_md*Delta_md/P4}
beta2 = ${fparse 6*E_mb*Delta_mb/P4}
beta_ss = ${fparse 6*E_ss*Delta_ss/P4}
a0 = 0.01

[ICs]
  [theta]
    type = FunctionIC
    variable = theta
    function = 'x/200'
  []

  [Upsilon_i]
    type = FunctionIC
    variable = Upsilon
    function = 'abs(x-100)/100'
  []
[]

[Materials]
  [bs0]
    type = PolarPFMBetaS0
    property_name = betaS0
    theta = theta
    a_beta = 3 # "It is also assumed that all a = 3."
    beta10 = ${beta1} # m-delta  g(kE,kdelta) (nJ/m)
    beta20 = ${beta2} # m-beta (nJ/m)
    derivative_order = 2
  []
  [beta21phi]
    type = PolarPFMPhi
    property_name = beta21phi
    a0 = ${a0}
    a_phi = 3
    beta21 = ${beta_ss} # delta-beta (nJ/m)
    upsilon = Upsilon
    derivative_order = 2
  []
  [psiL]
    type = PolarPFMPsiL
    property_name = psiL
    DeltaG10 = ${G1}
    DeltaG20 = ${G2}
    G0 = 0
    a_A = 3
    a_theta = 3

    A10 = ${A1}
    A20 = ${A2}
    A21 = ${Ab}
    theta = theta
    upsilon = Upsilon
    derivative_order = 2
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  nl_abs_tol = 1e-12

  nl_max_its = 30
  l_max_its = 100
  dt = 1e-1
  num_steps = 5
  automatic_scaling = true
[]

[Outputs]
  print_linear_residuals = false
  exodus = true
[]
