box_half_length = 50

[Mesh]
  type = MyTRIMMesh
  dim = 3
  nx = 10
  ny = 10
  nz = 10
  xmin = -${box_half_length}
  xmax =  ${box_half_length}
  ymin = -${box_half_length}
  ymax =  ${box_half_length}
  zmin = -${box_half_length}
  zmax =  ${box_half_length}
[]

[Problem]
  kernel_coverage_check = false
[]

[Variables]
  [./T]
  [../]
[]

[AuxVariables]
  # unit is atoms per nm^3
  [./c]
    initial_condition = 84.91232
  [../]
[]

[UserObjects]
  [./pka_gun]
    type = PKAGun
    m = 63.55
    Z = 29
    E = 50e3
    direction = '1 0 0'
    point = '0 0 0'
    num_pkas = 1
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = c
    M = '63.55'
    Z = '29'
    pka_generator = pka_gun
    length_unit = NANOMETER
    var_physical_meaning = NUMBER_DENSITY
    trim_module = ENERGY_DEPOSITION
  [../]
  [./runner]
    type = MyTRIMDiracRun
    rasterizer = rasterizer
  [../]
[]

[VectorPostprocessors]
  [./dirac_energy]
    type = MyTRIMDiracEnergyResult
    runner = runner
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
  hide = c
[]
