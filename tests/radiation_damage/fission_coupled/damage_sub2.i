[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50
  ny = 50
  xmin = -100
  xmax = 100
  ymin = -100
  ymax = 100
  elem_type = QUAD4
  uniform_refine = 1
[]

[Problem]
 type = FEProblem
 kernel_coverage_check = false
[]

[Variables]
 [./u]
 [../]
[]

#[BCs]
#  [./Periodic]
#    [./c_U]
#      variable = c_U
#      auto_direction = 'x y'
#    [../]
#  [../]
#[]

[AuxVariables]
  [./U_int]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./U_vac]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./O_int]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./O_vac]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./C_int]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./C_vac]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./c_U]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = FunctionIC
      function = fuel_conc
    [../]
  [../]
  [./c_O2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./c_C]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rho]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

#[Kernels]
#  [./dt]
#    type = TimeDerivative
#    variable = c_U
#  [../]
#  [./diff]
#    type = Diffusion
#    variable = c_U
#  [../]
#[]

[AuxKernels]
  [./U_int]
    type = MyTRIMElementResultAux
    variable = U_int
    runner = runner
    ivar = 0
    defect = INT
  [../]
  [./U_vac]
    type = MyTRIMElementResultAux
    variable = U_vac
    runner = runner
    ivar = 0
    defect = VAC
  [../]

  [./O_int]
    type = MyTRIMElementResultAux
    variable = O_int
    runner = runner
    ivar = 1
    defect = INT
  [../]
  [./O_vac]
    type = MyTRIMElementResultAux
    variable = O_vac
    runner = runner
    ivar = 1
    defect = VAC
  [../]

  [./C_int]
    type = MyTRIMElementResultAux
    variable = C_int
    runner = runner
    ivar = 2
    defect = INT
  [../]
  [./C_vac]
    type = MyTRIMElementResultAux
    variable = C_vac
    runner = runner
    ivar = 2
    defect = VAC
  [../]


  [./c_O2]
    type = ParsedAux
    variable = c_O2
    args = c_U
    function = 0#2*c_U
    execute_on = 'initial LINEAR'
  [../]
  [./c_C]
    type = ParsedAux
    variable = c_C
    args = c_U
    function = 1-c_U#1-3*c_U
    execute_on = 'initial LINEAR'
  [../]

  [./rho]
    type = MyTRIMDensityAux
    variable = rho
    rasterizer = rasterizer
  [../]
[]


[UserObjects]
  [./neutronics_fission_generator]
    type = PKAFissionFragmentNeutronics
    relative_density = c_U
    fission_rate = 1e-4
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = 'c_U  c_O2  c_C'
    M   = '235  16    12'
    Z   = '92   8     6'
    site_volume = 0.0435 # nm^3 per UO2 unit
    pka_generator = neutronics_fission_generator
  [../]
  [./runner]
    type = MyTRIMElementRun
    rasterizer = rasterizer
  [../]
[]

[Functions]
  [./fuel_conc]
    type = ParsedFunction
    value = 'r:=sqrt(x*x+y*y); Rlow:=low*R; Rhigh:=high*R; Z:= (r-high*R)/(r-low*R); if(r<Rlow,1,
                                                                                   if(r>Rhigh,0,
                                                                                     -(exp(2*Z)-1)/(exp(2*Z)+1)))'
                                                                                    #-(1/(high-low))*(r/R-1)+0.5))'
    #value = 'if((abs(x+y)+abs(x-y))<=sqrt(200),1,0)'
    vars = 'R low high'
    vals = '20 0.9 1.1'
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 2
  nl_abs_tol = 1e-10
  dt = 1
[]

[Outputs]
  exodus = true
  csv = true
[]
