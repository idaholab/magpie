#
# HMX example from page 8 of Physical Review B 89, 184102 (2014)
# delta-phase is solid 1
# beta-phase is solid 2
#

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 40
  ny = 40
[]

[Variables]
  [./Upsilon]
  [../]
  [./theta]
  [../]
[]

[Kernels]
  [./dtu]
    # general form would be SusceptibilityTimeDerivative with 1/L^Upsilon
    type = TimeDerivative
    variable = Upsilon
  [../]
  [./dLu]
    type = PolarPFMDerivative
    variable = Upsilon
    F = psiL
  [../]
  [./gradu]
    type = PolarPFMGradient
    variable = Upsilon
    v = theta
    F = beta21phi
  [../]
  [./diffu]
    type = MatDiffusion
    variable = Upsilon
    F = betaS0
    args = 'theta'
  [../]

  [./dtt]
    # general form would be SusceptibilityTimeDerivative with 1/L^theta
    type = TimeDerivative
    variable = theta
  [../]
  [./dLt]
    type = PolarPFMDerivative
    variable = theta
    F = psiL
  [../]
  [./gradt]
    type = PolarPFMGradient
    variable = theta
    v = Upsilon
    F = betaS0
  [../]
  [./difft]
    type = MatDiffusion
    variable = theta
    F = beta21phi
    args = 'Upsilon'
  [../]
[]

# i) disappearing melt
kE = 1.933

# ii) barrierless formation of melt
#kE = 3.39

# simulation temperature θ
T = 532 # (K)

kdelta = 1.0
beta21 = 2.4845

T10e = 550    # (K)
T20e = 532.14 # (K)
T21e = 432    # (K)

DeltaS10 = -793.79 # (kJ/m^3 K)
DeltaS20 = -935.45 # (kJ/m^3 K)
DeltaS21 = ${fparser DeltaS20 - DeltaS10} # (kJ/m^3 K)

T21c = -16616 # (K)
T10c = ${fparser kdelta/kE * DeltaS21/DeltaS10 * (T21c - T21e) + T10e}
T20c = ${fparser kdelta/kE * DeltaS21/DeltaS20 * (T21c - T21e) + T20e}

A21tilde = 0
A10c = ${fparse -3*DeltaS10}
A20c = ${fparse -3*DeltaS20}
A21c = ${fparse -3*DeltaS21}

[ICs]
  active = 'theta upsilon_i'

  [./theta]
    type = FunctionIC
    variable = theta
  [../]

  [./Upsilon_i]
    type = PolarPFMInterfaceIC
    variable = Upsilon
    a_beta = 3 # "It is also assumed that all a = 3."
    beta10 = ${fparser beta21/(kE*kdelta)} # m-delta  g(kE,kdelta) (nJ/m)
    beta20 = ${fparser beta21/(kE*kdelta)} # m-beta (nJ/m)
  [../]

  [./Upsilon_ii]
    type = ConstantIC
    variable = Upsilon
    value = 0.99
  [../]
[]

[Materials]
  [./bs0]
    type = PolarPFMBetaS0
    f_name = betaS0
    theta = theta
    a_beta = 3 # "It is also assumed that all a = 3."
    beta10 = ${fparser beta21/(kE*kdelta)} # m-delta  g(kE,kdelta) (nJ/m)
    beta20 = ${fparser beta21/(kE*kdelta)} # m-beta (nJ/m)
  [../]
  [./beta21phi]
    type = PolarPFMPhi
    f_name = beta21phi
    a0     = 3
    beta21 = ${beta21} # delta-beta (nJ/m)
    upsilon = Upsilon
  [../]
  [./psiL]
    type = PolarPFMPsiL
    f_name = psiL
    DeltaG10 = ${fparser -DeltaS10 * (T - T10e)}
    DeltaG21 = ${fparser -DeltaS21 * (T - T21e)}
    G0       = ?
    a_A      = 3
    a_theta  = 3

    A10      = ${fparser A10c*(T-T10c)}
    A20      = ${fparser A20c*(T-T20c)}
    A21      = ${fparser A21tilde + A21c*(T-T21c)} # A21tilde_c?
    theta = theta
    upsilon = Upsilon
  [../]
[]

[Executioner]
  type = Transient
[]

[Outputs]
  exodus = true
[]
