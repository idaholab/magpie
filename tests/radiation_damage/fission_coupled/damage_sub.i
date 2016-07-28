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
  uniform_refine = 0
[]

[Variables]
  [./c]
    initial_condition = 1.0
  [../]
[]

[AuxVariables]
  [./int]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vac]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fuel_conc]
    order = CONSTANT
    family = MONOMIAL
    [./InitialCondition]
      type = FunctionIC
      function = U02_conc
    [../]
  [../]
[]

[Kernels]
  [./dt]
    type = TimeDerivative
    variable = c
  [../]
  [./diff]
    type = Diffusion
    variable = c
  [../]
[]

[AuxKernels]
  [./int]
    variable = int
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 0
    defect = INT
  [../]
  [./vac]
    variable = vac
    type = MyTRIMElementResultAux
    runner = runner
    ivar = 0
    defect = VAC
  [../]
[]

[UserObjects]
  [./neutronics_fission_generator]
    type = PKAFissionFragmentNeutronics
    relative_density = 1
    fission_rate = 0.001
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = fuel_conc
    M = 40
    Z = 20
    site_volume = 0.0404 # nm^3 per UO2 unit
    pka_generator = neutronics_fission_generator
  [../]
  [./runner]
    type = MyTRIMElementRun
    rasterizer = rasterizer
  [../]
[]

[Postprocessors]
  [./int]
    type = ElementIntegralVariablePostprocessor
    variable = int
    execute_on = timestep_end
  [../]
  [./vac]
    type = ElementIntegralVariablePostprocessor
    variable = vac
    execute_on = timestep_end
  [../]
[]

[Functions]
  [./U02_conc]
    type = ParsedFunction
    value = 'r:=sqrt(x*x+y*y); Rlow:=low*R; Rhigh:=high*R; Z:= (r-high*R)/(r-low*R); if(r<Rlow,1,
                                                                                     if(r>Rhigh,0,
                                                                                     -(exp(2*Z)-1)/(exp(2*Z)+1)))'
                                                                                     #-(1/(high-low))*(r/R-1)+0.5))'
    vars = 'R low high'
    vals = '20 0.9 1.1'
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv = true
[]
