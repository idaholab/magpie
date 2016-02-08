[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
[]

[Variables]
  [./dummy]
  [../]
[]

[Kernels]
  [./dummy]
    type = Diffusion
    variable = dummy
  [../]
[]

[AuxVariables]
  [./c_U]
    [./InitialCondition]
      type = FunctionIC
      function = x
    [../]
  [../]
  [./c_O2]
    [./InitialCondition]
      type = FunctionIC
      function = 2*y
    [../]
  [../]
  [./rho]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./rho]
    type = MyTRIMDensityAux
    variable = rho
    rasterizer = rasterizer
  [../]
[]

[UserObjects]
  [./dummy]
    type = PKAConstant
    Z = 0
    m = 0
    E = 0
    pka_rate = 0
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = 'c_U   c_O2'
    M   = '253 16'
    Z   = '92  8'
    site_volume = 0.0404 # nm^3 per UO2 unit
    pka_generator = dummy
  [../]
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
[]
