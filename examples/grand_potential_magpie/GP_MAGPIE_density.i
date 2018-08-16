# This is based on Larry's grandpotential.i file in marmot-problems/Magpie/examples
# The domain is 2D 400 nm x 400 nm with a xenon bubble in the middle of UO2 matrix, radius = 44 nm
# The variables are chemical potentials of U vacancies and interstitials, Xe atoms
# This is solving the rate theory equation for wi, wv, wg
# It is using a source term from MAGPIE based on BCMC

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 30
  ny = 30
  xmin = 0
  xmax = 400
  ymin = 0
  ymax = 400
  # uniform_refine = 2
[]

[GlobalParams]
  op_num = 2
  var_name_base = etam
[]

[Variables]
  # chemical potentials for vacancies, interstitials and gas atoms
  [./wv]
  [../]
  [./wi]
  [../]
  [./wg]
  [../]

  # order parameters: etab0 is the bubble, etam0 is the grain, etam1 is currently not used
  [./etab0]
  [../]
  [./etam0]
  [../]
  [./etam1]
  [../]
[]

[AuxVariables]
  # grain boundaries
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]

  # concentration of Xe
  [./c_Xe]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  [../]


  # number density variables created for the rasterizer in Magpie
  [./rho_U235]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rho_U238]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rho_O]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rhog_var]
    order = CONSTANT
    family = MONOMIAL
  [../]


  # concentrations of (other) defects
  # vacancies
  [./cv]
    order = CONSTANT
    family = MONOMIAL
  [../]

  # interstitials
  [./ci]
    order = CONSTANT
    family = MONOMIAL
  [../]

  # number density of interstitials for the rasterizer
  [./rhoi_var]
    order = CONSTANT
    family = MONOMIAL
  [../]

  # number of defects generated each time step
  [./interstitial_rate_Xe]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vacancy_rate_U]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./interstitial_rate_U]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./interstitial_rate_i]
    order = CONSTANT
    family = MONOMIAL
  [../]

  # burnup in GWd/tHM
  [./burnup_var]
    order = CONSTANT
    family = MONOMIAL
  [../]

  # temperature
  [./T]
    initial_condition = 1200
  [../]

  # rate of absorbed defects due to dislocations
  [./absorbed_int]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./absorbed_vac]
    order = CONSTANT
    family = MONOMIAL
  [../]

  # Replacement atoms
  [./rep_in_235]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./rep_out_235]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./rep_in_238]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./rep_out_238]
    order = CONSTANT
    family = MONOMIAL
  [../]

  # Total number of uranium atoms lost to dislocation sinks
  [./total_U_lost]
    order = CONSTANT
    family = MONOMIAL
  [../]

  # UO2 matrix density
  [./UO2_density]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./relative_density]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1.0
  [../]
[]

[ICs]
  [./IC_etab0]
    type = FunctionIC
    variable = etab0
    function = ic_func_etab0
  [../]
  [./IC_etam0]
    type = FunctionIC
    variable = etam0
    function = ic_func_etam0
  [../]
  [./IC_etam1]
    type = FunctionIC
    variable = etam1
    function = ic_func_etam1
  [../]
  [./IC_wv]
    type = ConstantIC
    value = 0.0
    variable = wv
  [../]
  [./IC_wi]
    type = ConstantIC
    value = 0.0
    variable = wi
  [../]
  [./IC_wg]
    type = ConstantIC
    value = 0.0
    variable = wg
  [../]
[]

[Functions]
  # bubble IC
  [./ic_func_etab0]
    type = ParsedFunction
    vars = 'kappa   mu'
    vals = '0.5273  0.004688'
    value = 'r:=sqrt((x-200)^2+(y-200)^2);0.5*(1.0-tanh((r-44)*sqrt(mu/2.0/kappa)))'
  [../]

  # grain IC
  [./ic_func_etam0]
    type = ParsedFunction
    vars = 'kappa   mu'
    vals = '0.5273  0.004688'
    value = 'r:=sqrt((x-200)^2+(y-200)^2);0.5*(1.0+tanh((r-44)*sqrt(mu/2.0/kappa)))'
  [../]
  [./ic_func_etam1]
    type = ParsedFunction
    vars = 'kappa   mu'
    vals = '0.5273  0.004688'
    value = '0'
  [../]

  # burnup function
  [./fburnup] # in GWD/tHM
    type = ParsedFunction
    value = 'q*sigmaf_u235*flux*t*950'
    vars = 'q sigmaf_u235 flux'
    vals = '0.045 5.906e-22 1.5e13'
  [../]

  # function to define the total number of lost U atoms to dislocations with a postprocessor, in the auxvariable total_U_lost, to use it in the material property
  [./lost_U_func]
    type = ParsedFunction
    vals = total_U_lost_pp
    vars = total_U_lost_pp
    value = 'total_U_lost_pp'
  [../]
[]



[BCs]
  [./Periodic]
    [./All]
      auto_direction = 'x y' # 2D
    [../]
  [../]
[]

[Kernels]
# Order parameter eta_b0 for bubble phase
  [./ACb0_bulk]
    type = ACGrGrMulti
    variable = etab0
    v =           'etam0 etam1'
    gamma_names = 'gmm   gmm'
  [../]
  [./ACb0_sw]
    type = ACSwitching
    variable = etab0
    Fj_names  = 'omegab omegam'
    hj_names  = 'hb     hm'
    args = 'etam0 etam1 wv wi wg'
  [../]
  [./ACb0_int]
    type = ACInterface
    variable = etab0
    kappa_name = kappa
  [../]
  [./eb0_dot]
    type = TimeDerivative
    variable = etab0
  [../]
# Order parameter eta_m0 for matrix grain 1
  [./ACm0_bulk]
    type = ACGrGrMulti
    variable = etam0
    v =           'etab0 etam1'
    gamma_names = 'gmb   gmm'
  [../]
  [./ACm0_sw]
    type = ACSwitching
    variable = etam0
    Fj_names  = 'omegab omegam'
    hj_names  = 'hb     hm'
    args = 'etab0 etam1 wv wi wg'
  [../]
  [./ACm0_int]
    type = ACInterface
    variable = etam0
    kappa_name = kappa
  [../]
  [./em0_dot]
    type = TimeDerivative
    variable = etam0
  [../]
# Order parameter eta_m1 for matrix grain 2
  [./ACm1_bulk]
    type = ACGrGrMulti
    variable = etam1
    v =           'etab0 etam0'
    gamma_names = 'gmb   gmm'
  [../]
  [./ACm1_sw]
    type = ACSwitching
    variable = etam1
    Fj_names  = 'omegab omegam'
    hj_names  = 'hb     hm'
    args = 'etab0 etam0 wv wi wg'
  [../]
  [./ACm1_int]
    type = ACInterface
    variable = etam1
    kappa_name = kappa
  [../]
  [./em1_dot]
    type = TimeDerivative
    variable = etam1
  [../]
#Chemical potential for vacancies
  [./wv_dot]
    type = SusceptibilityTimeDerivative
    variable = wv
    f_name = chiv
    args = '' # in this case chi (the susceptibility) is simply a constant
  [../]
  [./Diffusion_v]
    type = MatDiffusion
    variable = wv
    D_name = Dchiv
    args = ''
  [../]
  [./coupled_v_etab0dot]
    type = CoupledSwitchingTimeDerivative
    variable = wv
    v = etab0
    Fj_names = 'rhovbub rhovmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1'
  [../]
  [./coupled_v_etam0dot]
    type = CoupledSwitchingTimeDerivative
    variable = wv
    v = etam0
    Fj_names = 'rhovbub rhovmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1'
  [../]
  [./coupled_v_etam1dot]
    type = CoupledSwitchingTimeDerivative
    variable = wv
    v = etam1
    Fj_names = 'rhovbub rhovmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1'
  [../]
  [./vac_source_U238] # Magpie source term for U238
    type = MyTRIMElementSource
    variable = wv
    runner = runner
    ivar = 1
    defect = VAC
    prefactor = 0.04092
  [../]
  [./vac_source_U235] # Magpie source term for U235
    type = MyTRIMElementSource
    variable = wv
    runner = runner
    ivar = 0
    defect = VAC
    prefactor = 0.04092
  [../]
  [./recombination_v] # Recombination term = Recombination rate * rhoi * rhov
    type = GrandPotentialRecombination
    variable = wv
    rho = rhov
    rho_r = rhoi
    value = 1
    D = Di
    omega = 0.04092
    a0 = 0.25
    Z = 50
    hm = hm
    args = 'wi etab0 etam0'
  [../]
  [./dislocation_sink_v] # Dislocation sink = sink strength * Diffusion coefficient * rho
     type = GrandPotentialSink
     variable = wv
     rho = rhov
     D = D
     mask = hm
     value = 1
     sink_strength = dislocation_density
     args = 'etab0 etam0 etam1'
 [../]

# Chemical potential for interstitials (same as vacancy)
  [./wi_dot]
    type = SusceptibilityTimeDerivative
    variable = wi
    f_name = chii
    args = ''
  [../]
  [./Diffusion_i]
    type = MatDiffusion
    variable = wi
    D_name = Dchii
    args = ''
  [../]
  [./coupled_i_etab0dot]
    type = CoupledSwitchingTimeDerivative
    variable = wi
    v = etab0
    Fj_names = 'rhoibub rhoimatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1'
  [../]
  [./coupled_i_etam0dot]
    type = CoupledSwitchingTimeDerivative
    variable = wi
    v = etam0
    Fj_names = 'rhoibub rhoimatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1'
  [../]
  [./coupled_i_etam1dot]
    type = CoupledSwitchingTimeDerivative
    variable = wi
    v = etam1
    Fj_names = 'rhoibub rhoimatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1'
  [../]
  [./source_i_U238]
    type = MyTRIMElementSource
    variable = wi
    runner = runner
    ivar = 1
    defect = INT
    prefactor = 0.04092
  [../]
  [./source_i_U235]
    type = MyTRIMElementSource
    variable = wi
    runner = runner
    ivar = 0
    defect = INT
    prefactor = 0.04092
  [../]
  [./source_i_i]
    type = MyTRIMElementSource
    variable = wi
    runner = runner
    ivar = 4
    defect = INT
    prefactor = 0.04092
  [../]
  [./source_i_v]
    type = MyTRIMElementSource
    variable = wi
    runner = runner
    ivar = 4
    defect = VAC
    prefactor = -0.04092
  [../]
  [./recombination_i]
    type = GrandPotentialRecombination
    variable = wi
    rho = rhoi
    rho_r = rhov
    value = 1
    D = Di
    omega = 0.04092
    a0 = 0.25
    Z = 50
    hm = hm
    args = 'wv etab0 etam0'
  [../]
  [./dislocation_sink_i]
     type = GrandPotentialSink
     variable = wi
     rho = rhoi
     mask = hm
     D = Di
     value = 1.02
     sink_strength = dislocation_density
     args = 'etab0 etam0 etam1'
   [../]

#Chemical potential for gas atoms
  [./wg_dot]
    type = SusceptibilityTimeDerivative
    variable = wg
    f_name = chig
    args = '' # in this case chi (the susceptibility) is simply a constant
  [../]
  [./Diffusion_g]
    type = MatDiffusion
    variable = wg
    D_name = Dchig
    args = ''
  [../]
  [./coupled_g_etab0dot]
    type = CoupledSwitchingTimeDerivative
    variable = wg
    v = etab0
    Fj_names = 'rhogbub rhogmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1'
  [../]
  [./coupled_g_etam0dot]
    type = CoupledSwitchingTimeDerivative
    variable = wg
    v = etam0
    Fj_names = 'rhogbub rhogmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1'
  [../]
  [./coupled_g_etam1dot]
    type = CoupledSwitchingTimeDerivative
    variable = wg
    v = etam1
    Fj_names = 'rhogbub rhogmatrix'
    hj_names = 'hb      hm'
    args = 'etab0 etam0 etam1'
  [../]
  [./source_g]
    type = MyTRIMElementSource
    variable = wg
    runner = runner
    ivar = 3
    defect = INT
    prefactor = 0.04092
  [../]
[]

[AuxKernels]
  # Numbers of defects created each time step
  [./interstitial_rate_Xe]
    type = MyTRIMElementResultAux
    runner = runner
    defect = INT
    variable = interstitial_rate_Xe
    ivar = 3
    execute_on = timestep_end
  [../]
  [./vacancy_rate_U]
    type = MyTRIMElementResultAux
    runner = runner
    defect = VAC
    variable = vacancy_rate_U
    ivar = 1
    execute_on = timestep_end
  [../]
  [./interstitial_rate_U]
    type = MyTRIMElementResultAux
    runner = runner
    defect = INT
    variable = interstitial_rate_U
    ivar = 1
    execute_on = timestep_end
  [../]
  [./interstitial_rate_i]
    type = MyTRIMElementResultAux
    runner = runner
    defect = INT
    variable = interstitial_rate_i
    ivar = 4
    execute_on = timestep_end
  [../]

  [./rep_in_rate_235]
    type = MyTRIMElementResultAux
    runner = runner
    defect = REPLACEMENT_IN
    variable = rep_in_235
    ivar = 0
    execute_on = timestep_end
  [../]
  [./rep_out_rate_235]
    type = MyTRIMElementResultAux
    runner = runner
    defect = REPLACEMENT_OUT
    variable = rep_out_235
    ivar = 0
    execute_on = timestep_end
  [../]
  [./rep_in_rate_238]
    type = MyTRIMElementResultAux
    runner = runner
    defect = REPLACEMENT_IN
    variable = rep_in_238
    ivar = 1
    execute_on = timestep_end
  [../]
  [./rep_out_rate_238]
    type = MyTRIMElementResultAux
    runner = runner
    defect = REPLACEMENT_OUT
    variable = rep_out_238
    ivar = 1
    execute_on = timestep_end
  [../]

  # Grain boundaries
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]

  # Number densities for different concentrations (U,O,Xe,Uint)
  [./rho_U235_aux] # U235
    type = MaterialRealAux
    property = rho_U235_matl
    variable = rho_U235
    execute_on = 'initial timestep_begin'
  [../]
  [./rho_U238_aux] # U238
    type = MaterialRealAux
    property = rho_U238_matl
    variable = rho_U238
    execute_on = 'initial timestep_begin'
  [../]
  [./rho_O] # Oxygen
    type = MaterialRealAux
    property = rho_O_matl
    variable = rho_O
    execute_on = 'initial timestep_begin'
  [../]
  [./rhog_aux] # Xenon
    type = MaterialRealAux
    property = rhog
    variable = rhog_var
    execute_on = 'initial timestep_begin'
  [../]
  [./rhoi_aux] # U interstitial
    type = MaterialRealAux
    property = rhoi
    variable = rhoi_var
    execute_on = 'initial timestep_begin'
  [../]

  # Number fractions of point defects
  [./cv_aux]
    type = MaterialRealAux
    property = cv_mat
    variable = cv
    execute_on = 'initial timestep_begin'
  [../]
  [./ci_aux]
    type = MaterialRealAux
    property = ci_mat
    variable = ci
    execute_on = 'initial timestep_begin'
  [../]
  [./c_Xe]
    type = MaterialRealAux
    property = c_Xe_matl
    variable = c_Xe
    execute_on = 'initial timestep_begin'
  [../]

  # U lost to dislocations
  [./lost_U_aux]
    type = FunctionAux
    variable = total_U_lost
    function = lost_U_func
    execute_on = timestep_end
  [../]

  # Auxkernels for sink rates
  # The value of absorbed_int/vac is integrated over the domain to get the total number of atoms absorbed (see postprocessor)
  # It is also multiplied by the time step because the absorption rate is per second
  [./disloc_absorption_int]
    type = SinkAbsorptionRate
    variable = absorbed_int
    rho = rhoi
    D = Di
    mask = hm
    sink_strength = disloc_density_dt
  [../]
  [./disloc_absorption_vac]
    type = SinkAbsorptionRate
    variable = absorbed_vac
    rho = rhov
    D = D
    mask = hm
    sink_strength = disloc_density_dt
  [../]

  # Burnup auxkernel defined by the function fburnup
  [./burnup_aux]
    type = FunctionAux
    function = fburnup
    variable = burnup_var
  [../]

  # UO2 density aukernel that defines the density from the material property
  # Density is not required as an auxvariable, but it could make coupling easier
  [./UO2_density_aux]
    type = MaterialRealAux
    variable = UO2_density
    property = UO2_density_matl
  [../]
  [./relative_density_aux]
    type = MaterialRealAux
    variable = relative_density
    property = relative_density_matl
  [../]
[]

[Materials]

  # Switching functions for bubble and matrix
  [./hb]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = hb
    all_etas = 'etab0 etam0 etam1'
    phase_etas = 'etab0'
    # outputs = exodus
  [../]
  [./hm]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = hm
    all_etas = 'etab0 etam0 etam1'
    phase_etas = 'etam0 etam1'
    outputs = exodus
  [../]

  # Grand potential densities for bubble and matrix
  [./omegab]
    type = DerivativeParsedMaterial
    args = 'wv wi wg'
    f_name = omegab
    material_property_names = 'Va kvbub cvbubeq kibub cibubeq kgbub cgbubeq f0'
    function = '-0.5*wv^2/Va^2/kvbub-wv/Va*cvbubeq-0.5*wi^2/Va^2/kibub-wi/Va*cibubeq-0.5*wg^2/Va^2/kgbub-wg/Va*cgbubeq+f0'
    derivative_order = 2
    #outputs = exodus
  [../]
  [./omegam]
    type = DerivativeParsedMaterial
    args = 'wv wi wg'
    f_name = omegam
    material_property_names = 'Va kvmatrix cvmatrixeq kimatrix cimatrixeq kgmatrix cgmatrixeq'
    function = '-0.5*wv^2/Va^2/kvmatrix-wv/Va*cvmatrixeq-0.5*wi^2/Va^2/kimatrix-wi/Va*cimatrixeq-0.5*wg^2/Va^2/kgmatrix-wg/Va*cgmatrixeq'
    derivative_order = 2
    #outputs = exodus
  [../]

  # Number densities defind in each phase for each defect (vacancy, interstitial, gas atoms)
  [./rhovbub]
    type = DerivativeParsedMaterial
    args = 'wv'
    f_name = rhovbub
    material_property_names = 'Va kvbub cvbubeq'
    function = 'wv/Va^2/kvbub + cvbubeq/Va'
    derivative_order = 2
    # outputs = exodus
  [../]
  [./rhovmatrix]
    type = DerivativeParsedMaterial
    args = 'wv'
    f_name = rhovmatrix
    material_property_names = 'Va kvmatrix cvmatrixeq'
    function = 'wv/Va^2/kvmatrix + cvmatrixeq/Va'
    derivative_order = 2
    # outputs = exodus
  [../]
  [./rhoibub]
    type = DerivativeParsedMaterial
    args = 'wi'
    f_name = rhoibub
    material_property_names = 'Va kibub cibubeq'
    function = 'wi/Va^2/kibub + cibubeq/Va'
    derivative_order = 2
    #outputs = exodus
  [../]
  [./rhoimatrix]
    type = DerivativeParsedMaterial
    args = 'wi'
    f_name = rhoimatrix
    material_property_names = 'Va kimatrix cimatrixeq'
    function = 'wi/Va^2/kimatrix + cimatrixeq/Va'
    derivative_order = 2
    #outputs = exodus
  [../]
  [./rhogbub]
    type = DerivativeParsedMaterial
    args = 'wg'
    f_name = rhogbub
    material_property_names = 'Va kgbub cgbubeq'
    function = 'wg/Va^2/kgbub + cgbubeq/Va'
    derivative_order = 2
    #outputs = exodus
  [../]
  [./rhogmatrix]
    type = DerivativeParsedMaterial
    args = 'wg'
    f_name = rhogmatrix
    material_property_names = 'Va kgmatrix cgmatrixeq'
    function = 'wg/Va^2/kgmatrix + cgmatrixeq/Va'
    derivative_order = 2
    #outputs = exodus
  [../]

  # Total number densities in the whole domain for each defect
  [./rhov]
    type = DerivativeParsedMaterial
    args = 'wv'
    f_name = rhov
    material_property_names = 'hm hb rhovmatrix(wv) rhovbub(wv)'
    function = 'hm * rhovmatrix + hb * rhovbub'
    derivative_order = 2
    # outputs = exodus
  [../]
  [./rhoi]
    type = DerivativeParsedMaterial
    args = 'wi'
    f_name = rhoi
    material_property_names = 'hm hb rhoimatrix(wi) rhoibub(wi)'
    function = 'hm * rhoimatrix + hb * rhoibub'
    derivative_order = 2
    # outputs = exodus
  [../]
  [./rhog]
    type = DerivativeParsedMaterial
    args = 'wg'
    f_name = rhog
    material_property_names = 'hm hb rhogmatrix(wg) rhogbub(wg)'
    function = 'hm * rhogmatrix + hb * rhogbub'
    derivative_order = 2
    # outputs = exodus
  [../]

  # U and O number densities
  [./rho_U235_matl]
    type = ParsedMaterial
    material_property_names = c_U235_matl
    f_name = rho_U235_matl
    constant_names = Va
    constant_expressions = 0.04092
    function = 'c_U235_matl/Va'
    # outputs = exodus
  [../]
  [./rho_U238_matl]
    type = ParsedMaterial
    material_property_names = c_U238_matl
    f_name = rho_U238_matl
    constant_names = Va
    constant_expressions = 0.04092
    function = 'c_U238_matl/Va'
    # outputs = exodus
  [../]
  [./rho_O_matl]
    type = ParsedMaterial
    material_property_names = 'c_O_matl'
    f_name = rho_O_matl
    constant_names = Va
    constant_expressions = 0.04092
    function = 'c_O_matl/Va/2'
    # outputs = exodus
  [../]

  # Constants
  [./const]
    type = GenericConstantMaterial
    prop_names =  'kappa   mu       L   D    Va      cvbubeq cgbubeq cibubeq  kgbub  kvbub kibub gmb     gmm T    Efvbar    Efgbar    kTbar     f0     tgrad_corr_mult  kappa_c kappa_op gamma_asymm Di'
    prop_values = '0.5273  0.004688 0.1 0.01 0.04092 0.5459  0.4541  0.0      1.41   1.41  1.41  0.9218 1.5 1200 7.505e-3  7.505e-3  2.588e-4  0.0    0.0              1.0     0.5273   1.5         1 '
  [../]

  # Equilibrium concentrations of defects
  [./cvmatrixeq]    #For values, see Li et al., Nuc. Inst. Methods in Phys. Res. B, 303, 62-27 (2013).
    type = ParsedMaterial
    f_name = cvmatrixeq
    material_property_names = 'T'
    constant_names        = 'kB           Efv'
    constant_expressions  = '8.6173324e-5 3.0'
    function = 'exp(-Efv/(kB*T))'
  [../]
  [./cimatrixeq]
    type = ParsedMaterial
    f_name = cimatrixeq
    material_property_names = 'T'
    constant_names        = 'kB           Efi'
    constant_expressions  = '8.6173324e-5 3.0'
    function = 'exp(-Efi/(kB*T))'
  [../]
  [./cgmatrixeq]
    type = ParsedMaterial
    f_name = cgmatrixeq
    material_property_names = 'T'
    constant_names        = 'kB           Efg'
    constant_expressions  = '8.6173324e-5 3.0'
    function = 'exp(-Efg/(kB*T))'
  [../]

  # See free energy expression
  [./kvmatrix_parabola]
    type = ParsedMaterial
    f_name = kvmatrix
    material_property_names = 'T  cvmatrixeq'
    constant_names        = 'c0v  c0g  a1                                               a2'
    constant_expressions  = '0.01 0.01 0.178605-0.0030782*log(1-c0v)+0.0030782*log(c0v) 0.178605-0.00923461*log(1-c0v)+0.00923461*log(c0v)'
    function = '((-a2+3*a1)/(4*(c0v-cvmatrixeq))+(a2-a1)/(2400*(c0v-cvmatrixeq))*T)'
    # outputs = exodus
  [../]
  [./kimatrix_parabola]
    type = ParsedMaterial
    f_name = kimatrix
    material_property_names = 'kvmatrix'
    function = 'kvmatrix'
  [../]
  [./kgmatrix_parabola]
    type = ParsedMaterial
    f_name = kgmatrix
    material_property_names = 'kvmatrix'
    function = 'kvmatrix'
  [../]

  # Material properties to define the number fractions of U, O, Xe, U vacancy and interstitial
  [./c_U235_matl]
    type = ParsedMaterial
    f_name = c_U235_matl
    material_property_names = 'hm'
    function = 'hm*0.01485'
    outputs = exodus
  [../]
  [./c_U238_matl]
    type = ParsedMaterial
    f_name = c_U238_matl
    material_property_names = 'hm'
    function = 'hm*0.31515'
    outputs = exodus
  [../]
  [./c_O_matl]
    type = ParsedMaterial
    f_name = c_O_matl
    material_property_names = 'hm'
    function = '0.66*hm'
    outputs = exodus
  [../]
  [./c_Xe_matl]
    type = ParsedMaterial
    f_name = c_Xe_matl
    material_property_names = 'hm hb Va rhogmatrix rhogbub'
    function = 'Va * (hm * rhogmatrix + hb * rhogbub)'
  [../]
  [./cv_mat]
    type = ParsedMaterial
    f_name = cv_mat
    material_property_names = 'hm hb Va rhovmatrix rhovbub'
    function = 'Va * (hm * rhovmatrix + hb * rhovbub)'
    outputs = exodus
  [../]
  [./ci_mat]
    type = ParsedMaterial
    f_name = ci_mat
    material_property_names = 'hm hb Va rhoimatrix rhoibub'
    function = 'Va * (hm * rhoimatrix + hb * rhoibub)'
    outputs = exodus
  [../]

  # Mobilities for Grand Potential Model
  [./Mobility_v]
    type = DerivativeParsedMaterial
    f_name = Dchiv
    material_property_names = 'D chiv'
    function = 'D*chiv'
    derivative_order = 2
    #outputs = exodus
  [../]
  [./Mobility_i]
    type = DerivativeParsedMaterial
    f_name = Dchii
    material_property_names = 'Di chii'
    function = 'Di*chii'
    derivative_order = 2
    #outputs = exodus
  [../]
  [./Mobility_g]
    type = DerivativeParsedMaterial
    f_name = Dchig
    material_property_names = 'Dtot chig'
    function = 'Dtot*chig'
    derivative_order = 2
    #outputs = exodus
  [../]

  # Susceptibilities
  [./chiv]
    type = DerivativeParsedMaterial
    f_name = chiv
    material_property_names = 'Va hb kvbub hm kvmatrix '
    function = '(hm/kvmatrix + hb/kvbub) / Va^2'
    derivative_order = 2
    # outputs = exodus
  [../]
  [./chii]
    type = DerivativeParsedMaterial
    f_name = chii
    material_property_names = 'Va hb kibub hm kimatrix '
    function = '(hm/kimatrix + hb/kibub) / Va^2'
    derivative_order = 2
    # outputs = exodus
  [../]
  [./chig]
    type = DerivativeParsedMaterial
    f_name = chig
    material_property_names = 'Va hb kgbub hm kgmatrix '
    function = '(hm/kgmatrix + hb/kgbub) / Va^2'
    derivative_order = 2
    # outputs = exodus
  [../]

  # Dislocation density per nm2
  [./dislocation_density]
    type = ParsedMaterial
    f_name = dislocation_density
    args = 'burnup_var'
    function = '(10^(2.2e-2*burnup_var+13.8))*1e-18' #This is an empirical expression  https://www.sciencedirect.com/science/article/pii/S002231151200637X
    # outputs = exodus
  [../]
  # Dislocation density per nm2 * time step in sec
  [./disloc_dens_dt]
    type = ParsedMaterial
    f_name = disloc_density_dt
    material_property_names = 'dislocation_density dt'
    function = 'dislocation_density * dt'
    # outputs = exodus
  [../]

  # Computes the difference between absorbed interstitials and vacancies only in the matrix (not in the interface nor the bubble phase)
  [./U_lost_mat]
    type = ParsedMaterial
    args = 'absorbed_int absorbed_vac'
    f_name = U_lost_mat
    function = 'if(absorbed_int > absorbed_vac, absorbed_int - absorbed_vac, 0)'
  [../]

  # Computes the change in the UO2 density in the matrix
  [./UO2_density_matl]
    type = ParsedMaterial
    args = 'total_U_lost'
    f_name = UO2_density_matl
    material_property_names = 'hm'
    constant_names = 'theo_dens initial_vol Va' # Va is the atomic volume
    constant_expressions = '10.97  1.537735e+05 0.04092'
    function = 'initial_vol * theo_dens/(initial_vol + Va*(total_U_lost))'
    # outputs = exodus
  [../]

  [./relative_density_matl]
    type = ParsedMaterial
    f_name = relative_density_matl
    material_property_names = 'UO2_density_matl hm'
    constant_names = 'theo_dens' # Va is the atomic volume
    constant_expressions = '10.97'
    function = 'if(hm > 0.95, UO2_density_matl/theo_dens, 0.0)'
  [../]

  # Time step material
  [./dt]
    type = TimeStepMaterial
  [../]

[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  # total number of defects in each time step
  [./total_Xe_interstitial_production]
    type = ElementIntegralVariablePostprocessor
    variable = interstitial_rate_Xe
  [../]

  [./total_U_vacancy_production]
    type = ElementIntegralVariablePostprocessor
    variable = vacancy_rate_U
  [../]

  [./total_U_interstitial_production]
    type = ElementIntegralVariablePostprocessor
    variable = interstitial_rate_U
  [../]

  [./total_i_interstitial_production]
    type = ElementIntegralVariablePostprocessor
    variable = interstitial_rate_i
  [../]
  #
  # [./total_rep_in_235_production]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = rep_in_235
  # [../]
  #
  # [./total_rep_out_235_production]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = rep_out_235
  # [../]
  #
  # [./total_rep_in_238_production]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = rep_in_238
  # [../]
  #
  # [./total_rep_out_238_production]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = rep_out_238
  # [../]

 # average number fraction of defect
  [./cv_average]
    type = ElementAverageValue
    variable = cv
  [../]
  [./ci_average]
    type = ElementAverageValue
    variable = ci
  [../]
  [./cg_average]
    type = ElementAverageValue
    variable = c_Xe
  [../]

 # initial matrix volume
  [./UO2_volume]
    type = ElementIntegralMaterialProperty
    mat_prop = hm
    execute_on = initial
  [../]

 # time step
  [./dt]
    type = TimestepSize
  [../]

 # Total number of absorbed U atoms in each time step
  [./lost_U_pp]
    type = ElementIntegralMaterialProperty
    mat_prop = U_lost_mat
  [../]

 # Cumulative number of absorbed u atoms at dislocations
  [./total_U_lost_pp]
    type = CumulativeValuePostprocessor
    postprocessor = lost_U_pp
  [../]

 # Average UO2 density
  [./my_density]
    type = ElementAverageValue
    variable = UO2_density
  [../]

[]

[Executioner]
  type = Transient
  nl_max_its = 15
  scheme = bdf2
  solve_type = NEWTON
  # solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type'
  petsc_options_value = 'asm       1               ilu'
  l_max_its = 15
  l_tol = 1.0e-4
  nl_rel_tol = 1.0e-8
  start_time = 0.0
  num_steps = 1500
  dtmax = 1e6
  nl_abs_tol = 1e-10
  [./TimeStepper]
    type = IterationAdaptiveDT
        dt = 1000
        optimal_iterations = 8
        iteration_window = 2
  [../]
[]

[UserObjects]
  active = 'rasterizer runner neutronics_fission_generator'
  [./neutronics_fission_generator]
    type = PKAFissionFragmentEmpirical
    relative_density = 'relative_density'
    fission_rate = 1.5e-08
  [../]
  [./xenon]
    type = PKAConstant
    pka_rate = 3e-8
    m = 131
    Z = 54
    E = 70e6
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = 'rho_U235  rho_U238   rho_O  rhog_var rhoi_var'
    M   = '235       238        16      135     238'
    Z   = '92        92         8       54      92'
    site_volume = 1 # nm^3 per UO2 unit
    periodic_var = wv
    pka_generator = neutronics_fission_generator
    length_unit = NANOMETER
    max_pka_count = 1000
    recoil_rate_scaling = 1e2
    r_rec = 5.45
  [../]
  [./runner]
    type = MyTRIMElementRun
    rasterizer = rasterizer
  [../]
[]

[Adaptivity]
 marker = errorfrac
 max_h_level = 2
 [./Indicators]
   [./error]
     type = GradientJumpIndicator
     variable = bnds
   [../]
 [../]
 [./Markers]
   [./errorfrac]
     type = ErrorFractionMarker
     coarsen = 0.1
     indicator = error
     refine = 0.7
   [../]
 [../]
[]

[Outputs]
  exodus = true
  csv = true
  file_base = density
[]
