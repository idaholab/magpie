[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 60
  ny = 60
  xmin = -100
  xmax = 100
  ymin = 0
  ymax = 200
  elem_type = TRI3
[]

[Variables]
  [./c]
    initial_condition = 1.0
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

[AuxKernels]
  [./int]
    variable = int
    type = MyTRIMAux
    runner = runner
    ivar = 0
    defect = INT
  [../]
  [./vac]
    variable = vac
    type = MyTRIMAux
    runner = runner
    ivar = 0
    defect = VAC
  [../]
[]

[UserObjects]
  [./thermal_fission]
    type = PKAFissionFragmentEmpirical
    relative_density = 1
    fission_rate = 0.001
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = c
    M   = 40
    Z   = 20
    site_volume = 0.0404 # nm^3 per UO2 unit
    pka_generator = thermal_fission
  [../]
  [./runner]
    type = MyTRIMRun
    rasterizer = rasterizer
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
[]
