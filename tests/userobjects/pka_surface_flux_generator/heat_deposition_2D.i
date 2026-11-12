#
# ions depositing energy and heat conduction in 2D
#

[Mesh]
  type = MyTRIMMesh
  dim = 3
  nx = 20
  ny = 20
  xmin = -100
  xmax = 100
  ymin = 0
  ymax = 200
[]

[Variables]
  [T]
  []
[]

[Kernels]
  [source]
    type = MyTRIMElementHeatSource
    variable = T
    runner = runner
  []
  [conduction]
    type = ADHeatConduction
    variable = T
    thermal_conductivity = 1 # dummy value
  []
  [time]
    type = TimeDerivative
    variable = T
  []
[]

[AuxVariables]
  [energy_dep]
    order = CONSTANT
    family = MONOMIAL
  []

  [./cC]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1
  [../]
[]

[AuxKernels]
  [energy_dep]
    variable = energy_dep
    type = MyTRIMElementEnergyAux
    runner = runner
  []
[]

[UserObjects]
  [./constant]
    type = PKASurfaceFluxGenerator
    E = 1e3
    Z = 14
    m = 28
    direction = '1 0 0'
    boundary = 'bottom'
    flux = 1e3
    dt = 1
    boundary_surface_area = 1
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = 'cC'
    M = '12'
    Z = '6'
    Edisp = '16.3'
    site_volume = 0.0404
    pka_generator = constant
    trim_module = ENERGY_DEPOSITION
  [../]
  [./runner]
    type = MyTRIMElementRun
    rasterizer = rasterizer
  [../]
[]

[Postprocessors]
  [E]
    type = ElementIntegralVariablePostprocessor
    variable = energy_dep
    execute_on = timestep_end
  []
  [pka_total_E]
    type = MyTRIMPKAInfo
    rasterizer = rasterizer
    value_type = TOTAL_ENERGY
  []
  [pka_total_Z]
    type = MyTRIMPKAInfo
    rasterizer = rasterizer
    value_type = TOTAL_CHARGE
  []
  [pka_total_m]
    type = MyTRIMPKAInfo
    rasterizer = rasterizer
    value_type = TOTAL_MASS
  []
  [pka_total_num]
    type = MyTRIMPKAInfo
    rasterizer = rasterizer
    value_type = TOTAL_NUMBER
  []
[]

[Executioner]
  type = Transient
  num_steps = 3
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
[]
