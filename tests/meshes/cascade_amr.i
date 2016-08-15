[Mesh]
  type = MyTRIMMesh
  dim = 2
  nx = 10
  ny = 10
  xmin = -100
  xmax = 100
  ymin = -100
  ymax = 100
  elem_type = QUAD4
[]

[Variables]
  [./c]
    [./InitialCondition]
      type = FunctionIC
      function = 'r:=sqrt(x^2+y^2)/40-1;if(r<0,1,if(r>1,0.5,cos(r*pi)/4+3/4))'
    [../]
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

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'x y'
    [../]
  [../]
[]

[UserObjects]
  [./constant]
    type = PKAConstant
    E = 1000
    Z = 60
    m = 120
    pka_rate = 0.01
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = c
    M = 40
    Z = 20
    site_volume = 0.0404 # nm^3 per UO2 unit
    pka_generator = constant
  [../]
  [./runner]
    type = MyTRIMElementRun
    rasterizer = rasterizer
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
  nl_abs_tol = 1e-10

  [./Adaptivity]
    initial_adaptivity = 2
    refine_fraction = 0.4
    coarsen_fraction = 0.8
    max_h_level = 2
  [../]
[]

[Outputs]
  exodus = true
[]
