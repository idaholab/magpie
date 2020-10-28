[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 200
  xmax = 200
[]

[Variables]
  [theta]
  []
  [Upsilon]
  []
[]

# i) disappearing melt
# KE = 1.933
KE = 4
KD = 1

# simulation temperature θ
P4 = 2.415 # scale factor for interface width
E_ss = ${units 1 J/m^2 -> eV/nm^2} # length conversion
E_md = ${fparse E_ss/KE}
E_mb = ${fparse E_ss/KE}
Delta_ss = ${units 3e-9 m -> nm} # length conversion
Delta_md = ${fparse Delta_ss/KD}
Delta_mb = ${fparse Delta_ss/KD}
beta1 = ${fparse 6*E_md*Delta_md/P4}
beta2 = ${fparse 6*E_mb*Delta_mb/P4}
x0 = 100

[ICs]
  [theta]
    type = FunctionIC
    variable = theta
    function = '1/(1 + exp(-${P4}*(x - ${x0})/${Delta_ss}))'
  []
  [Upsilon]
    type = PolarPFMInterfaceIC
    variable = Upsilon
    theta = theta
    a_beta = 3 # "It is also assumed that all a = 3."
    beta10 = ${beta1} # m-delta  g(kE,kdelta) (eV/nm)
    beta20 = ${beta2} # m-beta (eV/nm)
    p = ${P4}
  []
[]

[Problem]
  solve = false
[]

[Executioner]
  type = Steady
[]

[Outputs]
  print_linear_residuals = false
  exodus = true
[]
