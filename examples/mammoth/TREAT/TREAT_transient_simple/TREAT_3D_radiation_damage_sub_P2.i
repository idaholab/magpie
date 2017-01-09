# dimensions are in micometer
# typical treat fuel particle is rf = 20 mu-m
# domain size is 400^3 mu-m^3

#
# units:
# length: mu-m
# temperature: K
# mass: g
# energy: mu-J
# time: s
#

[Mesh]
  type = MyTRIMMesh
  dim = 3
  nx = 50
  ny = 50
  nz = 50
  xmin = -200.0
  xmax =  200.0
  ymin = -200.0
  ymax =  200.0
  zmin = -200.0
  zmax =  200.0
  uniform_refine = 0
[]

[GlobalParams]
  order = CONSTANT
  family = MONOMIAL
[]

[Variables]
  [./temperature]
    order = FIRST
    family = LAGRANGE
    initial_condition = 300
  [../]
[]

[Kernels]
  [./HeatConduction]
    type = HeatConduction
    variable = temperature
    diffusion_coefficient_dT_name = dthermal_conductivity/dtemperature
  [../]
  [./HeatStorage]
    type = HeatCapacityConductionTimeDerivative
    heat_capacity = cp_rho
    variable = temperature
  [../]
  [./HeatSource]
    type = CoupledForce
    v = heat_source
    variable = temperature
  [../]
[]

[AuxVariables]
  [./vac_U]
  [../]
  [./vac_O]
  [../]
  [./vac_C]
  [../]
  [./vacancies]
  [../]
  [./scaled_vacancies]
  [../]
  [./arho]
  [../]
  [./dpa]
  [../]

  ## energy deposition: Edep is correctly scaled!
  [./Edep_old]
  [../]
  [./Edep_current]
  [../]
  [./Edep]
  [../]

  # the heat source is in units of W / A^3
  [./heat_source]
  [../]
  # necessary to save the local power density pp as variable
  [./qflat]
  [../]

  ## variables separating in and out
  [./dpa_in]
  [../]
  [./dpa_out]
  [../]
  [./Edep_in]
  [../]
  [./Edep_out]
  [../]
  [./temperature_in]
  [../]
  [./temperature_out]
  [../]

  ## this fuel particle is an ellispoid with the same volume as the spherical fuel particle
  ## with r = 20 mu-m
  ## the semiaxes are a = b = kappa * c
  ## c = r / (kappa)**(2/3)
  ## a = b = kappa**(1/3) * r
  ## kappa = 10
  [./c_U]
    [./InitialCondition]
      type = SmoothSuperellipsoidIC
      x1 = 0
      y1 = 0
      z1 = 0
      a = 20
      b = 20
      c = 20
      n = 2
      invalue = 0.3333333
      outvalue = 0
      int_width = 0
    [../]
  [../]
  [./c_O]
    [./InitialCondition]
      type = SmoothSuperellipsoidIC
      x1 = 0
      y1 = 0
      z1 = 0
      a = 20
      b = 20
      c = 20
      n = 2
      invalue = 0.66666666
      outvalue = 0
      int_width = 0
    [../]
  [../]
  [./c_C]
    [./InitialCondition]
      type = SmoothSuperellipsoidIC
      x1 = 0
      y1 = 0
      z1 = 0
      a = 20
      b = 20
      c = 20
      n = 2
      invalue = 0
      outvalue = 3.4908
      int_width = 0
    [../]
  [../]

  [./rho_UO2]
    [./InitialCondition]
      type = SmoothSuperellipsoidIC
      x1 = 0
      y1 = 0
      z1 = 0
      a = 20
      b = 20
      c = 20
      n = 2
      invalue = 1
      outvalue = 0
      int_width = 0
    [../]
  [../]
[]

[AuxKernels]
  [./vacancies_U]
    type = MyTRIMElementResultAux
    variable = vac_U
    runner = runner
    ivar = 0
    defect = VAC
    execute_on = timestep_end
  [../]
  [./vacancies_O]
    type = MyTRIMElementResultAux
    variable = vac_O
    runner = runner
    ivar = 1
    defect = VAC
    execute_on = timestep_end
  [../]
  [./vacancies_C]
    type = MyTRIMElementResultAux
    variable = vac_C
    runner = runner
    ivar = 2
    defect = VAC
    execute_on = timestep_end
  [../]

  [./vacancies]
    type = AccumulateAux
    variable = vacancies
    accumulate_from_variable = 'vac_U vac_O vac_C'
    execute_on = timestep_end
  [../]

  [./scale_vacancies]
    type = ParsedAux
    variable = scaled_vacancies
    # r is the scaling factor accouting for the 300 times reduction in
    # fission rate
    constant_names =        'r'
    constant_expressions =  '10'
    function = 'r * vacancies'
    args = 'vacancies'
  [../]

  [./arho]
    type = AtomicDensityAux
    #
    # Check the units of this!!! Is it in Ang^-3 ????!!!!
    #
    variable = arho
    rasterizer = rasterizer
    execute_on = timestep_end
  [../]

  [./dpa]
    type = ParsedAux
    variable = dpa
    function = 'scaled_vacancies / arho'
    args = 'arho scaled_vacancies'
    execute_on = timestep_end
  [../]

  ## energy deposition aux kernels
  [./save_old]
    type = ParsedAux
    variable = Edep_old
    args = Edep
    function = Edep
    execute_on = timestep_begin
  [../]
  [./current_edep]
    type = MyTRIMElementEnergyAux
    variable = Edep_current
    runner = runner
    execute_on = timestep_end
  [../]
  [./accumulated_edep]
    type = ParsedAux
    variable = Edep
    # r is the scaling factor accouting for the 300 times reduction in
    # fission rate
    constant_names =        'r'
    constant_expressions =  '10'
    args = 'Edep_current Edep_old'
    function = 'r * Edep_current + Edep_old'
    execute_on = timestep_end
  [../]

  ## unfortunately need to set qflat through a function
  [./qflat]
    type = FunctionAux
    variable = qflat
    function = local_power_density
  [../]

  ## heat source aux:
  ## the heat source is in mu-W / mu-m^3
  ## r = the scaling factor accouting for the 300 times reduction in
  ## fission rate
  ## dt = timestep
  ## tau = fraction of energy released as non-fission products = 0.2
  ## ev_to_muJ: conversion from EV to micro Joule
  ## vol_conv: conversion from W/cm3 to mu-W/mu-m^3
  ##
  [./heat_source]
    type = ParsedAux
    variable = heat_source
    constant_names =        'r   dt   tau   ev_to_muJ   vol_conv'
    constant_expressions =  '10  0.01 0.132 1.602e-13   1e-6'
    args = 'Edep_current qflat'
    function = 'r * Edep_current / dt * ev_to_muJ + tau * qflat * vol_conv'
    execute_on = timestep_end
  [../]

  ## auxkernels for variables that are defined in/out
  [./dpa_in]
    type = ParsedAux
    variable = dpa_in
    function = 'dpa * rho_UO2'
    args = 'dpa rho_UO2'
    execute_on = timestep_end
  [../]

  [./dpa_out]
    type = ParsedAux
    variable = dpa_out
    function = 'dpa * (1 - rho_UO2)'
    args = 'dpa rho_UO2'
    execute_on = timestep_end
  [../]

  [./Edep_in]
    type = ParsedAux
    variable = Edep_in
    function = 'Edep * rho_UO2'
    args = 'Edep rho_UO2'
    execute_on = timestep_end
  [../]

  [./Edep_out]
    type = ParsedAux
    variable = Edep_out
    function = 'Edep * (1 - rho_UO2)'
    args = 'Edep rho_UO2'
    execute_on = timestep_end
  [../]

  [./T_in]
    type = ParsedAux
    variable = temperature_in
    # note division by volume of fuel to ensure we get an average
    function = 'temperature * rho_UO2 / 3.350400e+04'
    args = 'temperature rho_UO2'
    execute_on = timestep_end
  [../]

  [./T_out]
    type = ParsedAux
    variable = temperature_out
    # note division by volume of graphite to ensure we get an average
    function = 'temperature * (1 - rho_UO2) / 63966496.0'
    args = 'temperature rho_UO2'
    execute_on = timestep_end
  [../]
[]

[Materials]
  #
  # Thermal conducitivity and heat capacity from "Heat transfer simulation of the UO2 particle-graphite system in TREAT fuel"
  # by K. Mo et al., Nucl. Eng. Design, 293 (2015)
  #

  # thermal conductivity in units of: [k] = mu-W / (mu-m * K) == W / (m * K)
  [./thermalconductivity_mat]
    type = DerivativeParsedMaterial
    constant_names =       'r'
    constant_expressions = '1000.0'
    function = 'k_fuel := 100 / (6.548 + 23.533 * temperature / r) + 6400 / pow(temperature / r, 2.5) * exp(-16.35 * r / temperature);
                k_graphite := 1.046;
                rho_UO2 * k_fuel + (1 - rho_UO2) * k_graphite'
    f_name = thermal_conductivity
    args = 'temperature rho_UO2'
    derivative_order = 1
  [../]

  # heat capacity * rho in units of: [cp * rho] = mu-J / (mu-m^3 * K)
  #
  # conversion factors:
  # fuel: density 10.963 g / cm3, molar mass 267 g / mol, cp in J / (mol K) ==> 1e-6
  # graphite: density 1.72 g / cm3, cp is in J / kg / K; g / cm3 = 1e-12 * g / mu-m3, J / kg / K = 1e6 mu-J / (1000 g K) ==> 1e-9
  [./cp_rho_mat]
    type = ParsedMaterial
    constant_names =       'rho_fuel M_fuel rho_graphite'
    constant_expressions = '10.963   267    1.72'
    function = 'rT := temperature / 1000;
                cp_f := 52.1743 + 87.951 * rT - 82.2411 * pow(rT, 2) + 31.542 * pow(rT, 3) - 2.6334 * pow(rT, 4) - 0.7139 * pow(rT, 4);
                cp_rho_fuel := 1.0e-6 * rho_fuel / M_fuel * cp_f;
                cp_rho_graphite := 1e-9 * rho_graphite / (11.07 * pow(temperature, -1.644) + 0.0003688 * pow(temperature, 0.02191));
                rho_UO2 * cp_rho_fuel + (1 - rho_UO2) * cp_rho_graphite'
    f_name = cp_rho
    args = 'temperature rho_UO2'
  [../]

  # it is more convenient to merge the rho * cp statement to cancel out some
  # scaling issues
  [./density_mat]
    type = GenericConstantMaterial
    prop_names  = 'density'
    prop_values = '1'
  [../]
[]

[Functions]
  [./local_power_density]
    type = ParsedFunction
    vars = 'local_power_density'
    vals = 'local_power_density'
    value = 'local_power_density'
  [../]
[]

[UserObjects]
  [./neutronics_fission_generator]
    type = PKAFissionFragmentNeutronics
    relative_density = 1
    fission_rate = fission_rate
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = 'c_U  c_O  c_C'
    M   = '235  16    12'
    Z   = '92   8     6'
    site_volume = 0.0404 # nm^3 per UO2 unit
    pka_generator = neutronics_fission_generator
    interval = 1
    length_unit = MICROMETER
    trim_module = ENERGY_DEPOSITION
    analytical_energy_cutoff = 1.0e4
  [../]
  [./runner]
    type = MyTRIMElementRun
    rasterizer = rasterizer
  [../]
[]

[Postprocessors]
  # unit is W/cm3
  [./local_power_density]
    type = Receiver
    default = 1.274404e-01
  [../]
  [./fission_rate]
    type = ScalePostprocessor
    value = local_power_density
    # scaling factor: f = p * V / (Q * Vf)
    #
    # p: local power density in W / cm^3
    # V / Vf: ratio of total and fuel volume = 2571
    # Q = 192.9 MeV
    # conversion from MeV to J: 1.60218e-13
    # Conversion to mu-m multiply by 1e-12
    # r = 10
    # The 3 is because C_U is only 1/3
    # f = p [J/s/cm^3] * 2571 / (192.9 * 1.60218e-13 [J]) * 1e-12 [cm^3 / mu-m^3] / r * 3
    scaling_factor = 24.956277010364971
    execute_on = 'initial timestep_begin'
  [../]

  [./total_heat_source]
    type = ElementIntegralVariablePostprocessor
    variable = heat_source
  [../]

  [./grain_volume]
    type = ElementIntegralVariablePostprocessor
    variable = rho_UO2
  [../]

  ## pps for variables that are defined in/out
  [./total_dpa_in]
    type = ElementIntegralVariablePostprocessor
    variable = dpa_in
    execute_on = timestep_end
  [../]
  [./total_dpa_out]
    type = ElementIntegralVariablePostprocessor
    variable = dpa_out
    execute_on = timestep_end
  [../]
  [./total_Edep_in]
    type = ElementIntegralVariablePostprocessor
    variable = Edep_in
    execute_on = timestep_end
  [../]
  [./total_Edep_out]
    type = ElementIntegralVariablePostprocessor
    variable = Edep_out
    execute_on = timestep_end
  [../]

  [./average_fuel_temperature]
    type = ElementIntegralVariablePostprocessor
    variable = temperature_in
    execute_on = timestep_end
  [../]
  [./average_graphite_temperature]
    type = ElementIntegralVariablePostprocessor
    variable = temperature_out
    execute_on = timestep_end
  [../]
[]

[VectorPostprocessors]
  active = ''
  [./profile]
    type = SphericalAverage
    variable = 'dpa'
    radius = 40
    bin_number = 8
  [../]
[]

[Adaptivity]
  [./Markers]
    [./fuel]
      type = ValueThresholdMarker
      variable = rho_UO2
      coarsen = 0
      refine = 0.5
    [../]
  [../]
  initial_marker = fuel
  initial_steps = 3
[]

[Executioner]
  type = Transient
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-7
  l_tol = 1.0e-3
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_max_iter -pc_hypre_boomeramg_tol'
  petsc_options_value = 'hypre boomeramg 100 20 1.0e-6'
  [./TimeStepper]
    type = ConstantDT
    dt = 1.0e-2
    growth_factor = 1.5
  [../]
[]

[Outputs]
  exodus = true
  csv = true
[]
