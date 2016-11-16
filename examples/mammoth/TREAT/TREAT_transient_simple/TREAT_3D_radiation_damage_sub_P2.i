# dimensions are in Angstrom (A = 1e-10 m)
# typical treat fuel particle is rf = 20 mu-m = 2e5 A
# typical distance between fuel particles is D = 0.893 * rf / cube_root(1/2571)
# D = 244 mu-m = 2.44e6 A.
# the typical fuel particle distance is much larger than the mfp so we cut the domain
# at 4 * rf = 80 mu-m = 8e5 A
[Mesh]
  type = MyTRIMMesh
  dim = 3
  nx = 20
  ny = 20
  nz = 20
  xmin = -4.0e5
  xmax = 4.0e5
  ymin = -4.0e5
  ymax = 4.0e5
  zmin = -4.0e5
  zmax = 4.0e5
  uniform_refine = 0
[]

[GlobalParams]
  order = CONSTANT
  family = MONOMIAL
[]

[Problem]
  kernel_coverage_check = false
[]

[Variables]
  [./c]
    initial_condition = 1.0
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

  [./c_U]
    [./InitialCondition]
      type = SmoothCircleIC
      x1 = 0
      y1 = 0
      z1 = 0
      radius = 2.0e5
      invalue = 0.3333333
      outvalue = 0
      int_width = 5e4
    [../]
  [../]
  [./c_O]
    [./InitialCondition]
      type = SmoothCircleIC
      x1 = 0
      y1 = 0
      z1 = 0
      radius = 2.0e5
      invalue = 0.66666666
      outvalue = 0
      int_width = 5e4
    [../]
  [../]
  [./c_C]
    [./InitialCondition]
      type = SmoothCircleIC
      x1 = 0
      y1 = 0
      z1 = 0
      radius = 2.0e5
      invalue = 0
      outvalue = 3.4908
      int_width = 5e4
    [../]
  [../]

  [./rho_UO2]
    [./InitialCondition]
      type = SmoothCircleIC
      x1 = 0
      y1 = 0
      z1 = 0
      radius = 2.0e5
      invalue = 1
      outvalue = 0
      int_width = 5e4
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
    function = '300.0 * vacancies'
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
[]

[UserObjects]
  [./neutronics_fission_generator]
    type = PKAFissionFragmentNeutronics
    relative_density = 1
    #fission_rate = 3.2498e-6
    # 1.0e-12 yields 2e4 PKAs
    fission_rate = fission_rate
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = 'c_U  c_O  c_C'
    M   = '235  16    12'
    Z   = '92   8     6'
    site_volume = 0.0404 # nm^3 per UO2 unit
    pka_generator = neutronics_fission_generator
    print_pka_statistics = true
  [../]
  [./runner]
    type = MyTRIMElementRun
    rasterizer = rasterizer
  [../]
[]

[Postprocessors]
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
    # Q = 200 MeV = 3.204e-11 J
    # Conversion to A multiply by 1e-24
    # f = p [J/s/cm^3] * 2571 / 3.204e-11 [J] * 1e-24 [cm^3 / A]
    #
    # To save execution time we reduce the number of PKAs by factor of 30
    # scaling_factor = 8.0243e-11 / 300
    scaling_factor = 2.6748e-13
    execute_on = 'initial timestep_begin'
  [../]
[]

[VectorPostprocessors]
  [./profile]
    type = SphericalAverage
    variable = 'dpa'
    radius = 4e5
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
[]

[Outputs]
  exodus = true
  csv = true
[]
