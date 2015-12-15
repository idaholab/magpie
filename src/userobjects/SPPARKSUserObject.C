#include "SPPARKSUserObject.h"

#include "MooseRandom.h"
#include "libmesh/mesh_tools.h"
#include <numeric>

template<>
InputParameters validParams<SPPARKSUserObject>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addParam<std::string>("file", "", "SPPARKS input file");
  params.addParam<bool>("spparks_only", false, "Whether to run SPPARKS independently of MOOSE");

  params.addParam<std::vector<unsigned> >("from_ivar", std::vector<unsigned>(), "Index into SPPARKS iarray.  This data will be extracted from SPPARKS.");
  params.addParam<std::vector<unsigned> >("from_dvar", std::vector<unsigned>(), "Index into SPPARKS darray.  This data will be extracted from SPPARKS.");
  params.addParam<std::vector<unsigned> >("to_ivar", std::vector<unsigned>(), "Index into SPPARKS iarray.  This data will be copied to SPPARKS.");
  params.addParam<std::vector<unsigned> >("to_dvar", std::vector<unsigned>(), "Index into SPPARKS darray.  This data will be copied to SPPARKS.");

  params.addParam<std::vector<std::string> >("int_vars", "Integer Vars to send to SPPARKS.");
  params.addParam<std::vector<std::string> >("double_vars", "Double Vars to send to SPPARKS.");
  params.addParam<std::vector<std::string> >("sol_vars", "Double solving-for Vars obtained from SPPARKS.");

  params.addParam<Real>("xmin", 0.0, "Lower X Coordinate of the generated mesh");
  params.addParam<Real>("ymin", 0.0, "Lower Y Coordinate of the generated mesh");
  params.addParam<Real>("zmin", 0.0, "Lower Z Coordinate of the generated mesh");
  params.addParam<Real>("xmax", 1.0, "Upper X Coordinate of the generated mesh");
  params.addParam<Real>("ymax", 1.0, "Upper Y Coordinate of the generated mesh");
  params.addParam<Real>("zmax", 1.0, "Upper Z Coordinate of the generated mesh");

  params.addParam<bool>("init_spparks", false, "Set values of SPPARKS variables at time zero.");
  params.addParam<bool>("one_time_run", false, "Run SPPARKS just once at time zero."); // added by YF
  params.addParam<Real>("time_spparks_time_ratio", 1, "Ratio of time to SPPARKS time.");

  // Hide from input file dump
  // params.addPrivateParam<std::string>("built_by_action", "");
  return params;
}

SPPARKSUserObject::SPPARKSUserObject(const InputParameters & params) :
    GeneralUserObject(params),
    _spparks(NULL),
    _file(getParam<std::string>("file")),
    _spparks_only(getParam<bool>("spparks_only")),
    _from_ivar(getParam<std::vector<unsigned> >("from_ivar")),
    _from_dvar(getParam<std::vector<unsigned> >("from_dvar")),
    _to_ivar(getParam<std::vector<unsigned> >("to_ivar")),
    _to_dvar(getParam<std::vector<unsigned> >("to_dvar")),
    _xmin(getParam<Real>("xmin")),
    _ymin(getParam<Real>("ymin")),
    _zmin(getParam<Real>("zmin")),
    _xmax(getParam<Real>("xmax")),
    _ymax(getParam<Real>("ymax")),
    _zmax(getParam<Real>("zmax")),
    _init_spparks(getParam<bool>("init_spparks")),
    _time_ratio(getParam<Real>("time_spparks_time_ratio")),
    _initialized(false),
    _one_time_run(getParam<bool>("one_time_run")),
    _n_spparks_run(0),
    _int_vars(),
    _double_vars(),
    _sol_vars(),
    _last_time(std::numeric_limits<Real>::min())
{
  if (isParamValid("int_vars"))
  {
    const std::vector<std::string> & names = getParam<std::vector<std::string> >("int_vars");
    for (unsigned i = 0; i < names.size(); ++i)
      _int_vars.push_back(&_fe_problem.getVariable(0, names[i]));
  }
  if (_int_vars.size() != _to_ivar.size())
    mooseError("Mismatch with int_vars and to_ivar");

  if (isParamValid("double_vars"))
  {
    const std::vector<std::string> & names = getParam<std::vector<std::string> >("double_vars");
    for (unsigned i = 0; i < names.size(); ++i)
      _double_vars.push_back(&_fe_problem.getVariable(0, names[i]));
  }
  if (_double_vars.size() != _to_dvar.size())
    mooseError("Mismatch with double_vars and to_dvar");

  if (isParamValid("sol_vars"))
  {
    const std::vector<std::string> & names = getParam<std::vector<std::string> >("sol_vars");
    for (unsigned i = 0; i < names.size(); ++i)
      _sol_vars.push_back(&_fe_problem.getVariable(0, names[i]));
  }

  _console << "\n>>>> STARTING SPPARKS <<<<\n";
  spparks_open(0, NULL, _communicator.get(), &_spparks);
  if (!_spparks)
    mooseError("Error initializing SPPARKS");

  _console << "\n>>>> RUNNING SPPARKS FILE " << _file << " <<<<\n";
  spparks_file(_spparks, _file.c_str());

  int * iptr;
  double * dptr;
  double ** ddptr;

  // Extract and print information about the SPPARKS internals
  getSPPARKSPointer(iptr, "dimension");
  _dim = *iptr;
  _console << "\n>>>> SPPARKS DIMENSION: " << _dim << " <<<<\n";

  getSPPARKSPointer(dptr, "boxxlo");
  double xlo = *dptr;
  _console << "\n>>>> SPPARKS BOXXLO: " << *dptr << " <<<<\n";

  getSPPARKSPointer(iptr, "nlocal");
  int nlcl = *iptr;
  _console << "\n>>>> SPPARKS NLOCAL: " << *iptr << " <<<<\n";

  getSPPARKSPointer(iptr, "id");
  _console << "n>>>> SPPARKS ID array pointer: " << iptr << " <<<<\n";

  getSPPARKSPointer(ddptr, "xyz");
  _console << "\n>>>> SPPARKS XYZ array pointer: " << ddptr << " <<<<\n";

  // for (unsigned i = 0; i < nlcl; ++i)
  // {
  //   std::cout << id_array[i] << "\t"
  //             << xyz_array[i][0] << " "
  //             << xyz_array[i][1] << " "
  //             << xyz_array[i][2] << std::endl;
  // }
}

SPPARKSUserObject::~SPPARKSUserObject()
{
  spparks_close(_spparks);
}

void
SPPARKSUserObject::initialize()
{
  if (_spparks_only) return;

  // getSPPARKSData();
  // setSPPARKSData();
}

int
SPPARKSUserObject::getIntValue(unsigned fem_node_id, unsigned index) const
{
  return _spparks_only ? 0 : getValue(_int_data_for_fem, fem_node_id, index);
}

Real
SPPARKSUserObject::getDoubleValue(unsigned int fem_node_id, unsigned index) const
{
  return _spparks_only ? 0.0 : getValue(_double_data_for_fem, fem_node_id, index);
}

char *
SPPARKSUserObject::runSPPARKSCommand(const std::string & cmd)
{
  return spparks_command(_spparks, cmd.c_str());
}

void
SPPARKSUserObject::getSPPARKSData()
{
  // _console << "\ngetSPPARKSData\n";
  // Update the integer data
  for (unsigned i = 0; i < _from_ivar.size(); ++i)
    getSPPARKSData(_int_data_for_fem[_from_ivar[i]], "iarray", _from_ivar[i]);

  // Update the double data
  for (unsigned i = 0; i < _from_dvar.size(); ++i)
    getSPPARKSData(_double_data_for_fem[_from_dvar[i]], "darray", _from_dvar[i]);
}

void
SPPARKSUserObject::setSPPARKSData()
{
  // _console << "\nsetSPPARKSData\n";
  // Update the integer data
  int * ip = NULL;
  for (unsigned i = 0; i < _to_ivar.size(); ++i)
    setSPPARKSData(ip, "iarray", _to_ivar[i], *_int_vars[i]);

  // Update the double data
  double * dp = NULL;
  for (unsigned i = 0; i < _to_dvar.size(); ++i)
    setSPPARKSData(dp, "darray", _to_dvar[i], *_double_vars[i]);
}

void
SPPARKSUserObject::setFEMData()
{
  // _console << "\nsetFEMData\n";
  // Update the double solving variables using the data from SPPARKS, added by YF
  for (unsigned int i = 0; i < _sol_vars.size(); ++i)
    setFEMData<double>("darray", _from_dvar[i], *_sol_vars[i]);
}

void
SPPARKSUserObject::execute()
{
  if (_one_time_run && _n_spparks_run > 0) return;
  if (_spparks_only) return;

  if (_init_spparks)
  {
    initSPPARKS();
    _init_spparks = false;
  }

  if (_t != _last_time)
  {
    _last_time = _t;

    // set variables in SPPARKS defined by to_ivar and to_dvar
    setSPPARKSData();

    // Run SPPARKS over a certain time defined by sp_time
    const Real sp_time = getSPPARKSTime(_dt);
    std::stringstream cmd;
    cmd << "run ";
    cmd << sp_time;
    cmd << " pre no" << std::endl;
    runSPPARKSCommand(cmd.str());

    // times that SPPARKS has been called
    _n_spparks_run ++;

    // obtain data from SPPARKS
    // getSPPARKSData update the auxvariables defined by from_ivar and from_dvar
    // setFEMData update the MOOSEvariables defined by sol__vars
    getSPPARKSData();
    setFEMData();
  }
}

void
SPPARKSUserObject::initSPPARKS()
{
  // Set SPPARKS values based on 3.3.2 in Veena's paper
  // This is needed only to couple POTTS model with MARMOT

  // Loop over local FEM nodes
  // For each node, pick a random spin
  // Half alpha phase, half beta phase
  // Composition for alpha phase is 0.25; for beta phase, 0.75.

  // _console << "\ninitSPPARKS\n";
  if (_from_ivar.size() != 2)
    mooseError("Must have two integer variables from SPPARKS");

  if (_double_vars.size() != 1)
    mooseError("Must have one double variable to send to SPPARKS");

  SystemBase & sys = _double_vars[0]->sys();
  NumericVector<Number> & solution = sys.solution();

  std::map<libMesh::dof_id_type, int> ints[2];
  std::map<libMesh::dof_id_type, Real> comp;
  ConstNodeRange & node_range = *_fe_problem.mesh().getLocalNodeRange();
  for (ConstNodeRange::const_iterator i = node_range.begin(); i < node_range.end(); ++i)
  {
    int value = int(100*MooseRandom::rand()) + 1;
    ints[0][(*i)->id()] = value;
    ints[1][(*i)->id()] = value < 51;

    Real comp = 0.75;
    if (value < 51)
      comp = 0.25;

    // Set data
    solution.set((*i)->dof_number(sys.number(), _double_vars[0]->number(), 0), comp);
  }
  solution.close();

  int * pint = NULL;
  char iarray[] = "iarray";
  for (unsigned i(0); i < 2; ++i)
  {
    getSPPARKSDataPointer(pint, iarray, _from_ivar[i]);

    // Index into data is SPPARKS node id.
    for (std::multimap<FEMID, SPPARKSID>::const_iterator it = _fem_to_spparks.begin(); it != _fem_to_spparks.end(); ++it)
      pint[it->second.id] = ints[i][it->first.id];

    // Copy data across processors
    if (n_processors() > 1)
      sendRecvFEMData(ints[i], pint);
  }
}

void
SPPARKSUserObject::initialSetup()
{
  if (_spparks_only || _initialized)
    return;

  _initialized = true;

  // Initialize communication maps
  // _console << "\ninitialSetup: begin\n";

  // 1. Get on-processor map from SPPARKS ID to FEM ID
  int * iptr;
  getSPPARKSPointer(iptr, "nlocal");
  int nlcl = *iptr;

  double ** xyz_array;
  getSPPARKSPointer(xyz_array, "xyz");

  std::set<SPPARKSID> spparks_id; // Local SPPARKS nodes
  for (int i = 0; i < nlcl; ++i)
  {
    Point p(xyz_array[i][0],
            _dim > 1 ? xyz_array[i][1] : 0,
            _dim > 2 ? xyz_array[i][2] : 0);
    spparks_id.insert(SPPARKSID(i, p));
  }
  _num_local_spparks_nodes = spparks_id.size();

  // MOOSE FEM nodes are set up with periodic bcs (nodes at each end of periodicity),
  // but SPPARKS nodes are not.  We may have two or more MOOSE FEM nodes at SPPARKS
  // node locations.
  std::multiset<FEMID> fem_id; // Local MOOSE FEM nodes
  ConstNodeRange & node_range = *_fe_problem.mesh().getLocalNodeRange();
  for (ConstNodeRange::const_iterator i = node_range.begin(); i < node_range.end(); ++i)
  {
    Point coor(**i);

    if (coor(0) == _xmax)
      coor(0) = _xmin;
    if (coor(1) == _ymax)
      coor(1) = _ymin;
    if (coor(2) == _zmax)
      coor(2) = _zmin;

    fem_id.insert(FEMID((*i)->id(), coor));
  }
  _num_local_fem_nodes = fem_id.size();

  // SPPARKS nodes not found on this processor
  std::set<SPPARKSID> unmatched_spparks;

  for (std::set<SPPARKSID>::iterator i = spparks_id.begin(); i != spparks_id.end(); ++i)
  {
    std::multiset<FEMID>::iterator fem_iter = fem_id.find(*i);
    if (fem_iter != fem_id.end())
      _spparks_to_fem.insert(std::pair<SPPARKSID, FEMID>(*i, *fem_iter)); // SPPARKSID to FEMID
    else
      unmatched_spparks.insert(*i);
  }

  for (std::multiset<FEMID>::iterator i = fem_id.begin(); i != fem_id.end(); ++i)
  {
    std::set<SPPARKSID>::iterator spparks_iter = spparks_id.find(*i);
    if (spparks_iter != spparks_id.end())
      _fem_to_spparks.insert(std::pair<FEMID, SPPARKSID>(*i, *spparks_iter)); // FEMID to SPPARKSID
    // else
    //   _spparks_to_proc.insert(std::pair<SPPARKSID, unsigned>(*i, -1));
  }

  const unsigned num_procs = n_processors();

  if (num_procs == 1)
  {
    if (spparks_id.size() != _spparks_to_fem.size())
    {
      // Error skipped to allow unmatched meshes
      // mooseError("Did not find MOOSE FEM node for each SPPARKS node, " << spparks_id.size()
      //            << ", " << fem_id.size() << ", " << _spparks_to_fem.size());

      _console << "\n  spparks size  "<< spparks_id.size() << "not equal to fem size  " << fem_id.size()  << '\n';
    }

    return;
  }

  // _console << "\ninitialSetup: 1\n";

  // 2. Get send map (spparks id -> proc id)

  // A. Get local SPPARKS bounding box
  // TODO: Expand bounding boxes by half cell width?
  //       May not be necessary.

  unsigned proc_id = processor_id();

  std::vector<Real> spparks_bounds(num_procs*6, 0.0);
  unsigned offset = proc_id*6;
  Real * s_bounds = &spparks_bounds[0] + offset;
  s_bounds[0] = s_bounds[1] = s_bounds[2] = std::numeric_limits<Real>::max();
  s_bounds[3] = s_bounds[4] = s_bounds[5] = std::numeric_limits<Real>::min();
  for (std::set<SPPARKSID>::const_iterator i = unmatched_spparks.begin(); i != unmatched_spparks.end(); ++i)
  {
    s_bounds[0] = std::min(s_bounds[0], i->coor(0));
    s_bounds[1] = std::min(s_bounds[1], i->coor(1));
    s_bounds[2] = std::min(s_bounds[2], i->coor(2));
    s_bounds[3] = std::max(s_bounds[3], i->coor(0));
    s_bounds[4] = std::max(s_bounds[4], i->coor(1));
    s_bounds[5] = std::max(s_bounds[5], i->coor(2));
  }
  MeshTools::BoundingBox spparks_bb(Point(s_bounds[0], s_bounds[1], s_bounds[2]),
                                     Point(s_bounds[3], s_bounds[4], s_bounds[5]));
  _communicator.sum(spparks_bounds);

  // _console << "\ninitialSetup: 2A\n";

  //
  // B: Get MOOSE bounding boxes
  //

  //
  // TODO: These bboxes could/should be based on non-matched fem nodes.
  //
  std::vector<Real> fem_bounds(num_procs*6, 0.0);
  Real * e_bounds = &fem_bounds[0] + offset;
  e_bounds[0] = e_bounds[1] = e_bounds[2] = std::numeric_limits<Real>::max();
  e_bounds[3] = e_bounds[4] = e_bounds[5] = std::numeric_limits<Real>::min();
  for (std::multiset<FEMID>::const_iterator i = fem_id.begin(); i != fem_id.end(); ++i)
  {
    e_bounds[0] = std::min(e_bounds[0], i->coor(0));
    e_bounds[1] = std::min(e_bounds[1], i->coor(1));
    e_bounds[2] = std::min(e_bounds[2], i->coor(2));
    e_bounds[3] = std::max(e_bounds[3], i->coor(0));
    e_bounds[4] = std::max(e_bounds[4], i->coor(1));
    e_bounds[5] = std::max(e_bounds[5], i->coor(2));
  }
  MeshTools::BoundingBox fem_bb(Point(e_bounds[0], e_bounds[1], e_bounds[2]),
                                 Point(e_bounds[3], e_bounds[4], e_bounds[5]));
  _communicator.sum(fem_bounds);

  // _console << "\ninitialSetup: 2B\n";

  //
  // C: Get number of processors that overlap my SPPARKS and MOOSE domains
  //

  std::vector<unsigned> procs_overlapping_spparks_domain;
  std::vector<unsigned> procs_overlapping_fem_domain;
  for (unsigned i = 0; i < num_procs; ++i)
  {
    if (i == proc_id)
      continue;

    offset = i * 6;
    e_bounds = &fem_bounds[0] + offset;
    MeshTools::BoundingBox e_box(Point(e_bounds[0], e_bounds[1], e_bounds[2]),
                                  Point(e_bounds[3], e_bounds[4], e_bounds[5]));
    if (spparks_bb.intersect(e_box))
      procs_overlapping_spparks_domain.push_back(i);

    s_bounds = &spparks_bounds[0] + offset;
    MeshTools::BoundingBox s_box(Point(s_bounds[0], s_bounds[1], s_bounds[2]),
                                  Point(s_bounds[3], s_bounds[4], s_bounds[5]));
    if (fem_bb.intersect(s_box))
      procs_overlapping_fem_domain.push_back(i);
  }

  // _console << "\ninitialSetup: 2C\n";
  // _console << "\noverlapping spparks domain:  "<< procs_overlapping_spparks_domain.size() << '\n';
  // _console << "\noverlapping fem domain:  "<< procs_overlapping_fem_domain.size() << '\n';

  // D: Communicate number of MOOSE FEM nodes, number of SPPARKS nodes

  std::vector<unsigned> num_fem_nodes(procs_overlapping_spparks_domain.size(), 0); // Number remote fem nodes for each proc overlapping my spparks
  std::vector<unsigned> num_spparks_nodes(procs_overlapping_fem_domain.size(), 0);

  std::vector<MPI_Request> recv_request1(std::max(procs_overlapping_spparks_domain.size(), procs_overlapping_fem_domain.size()));
  std::vector<MPI_Request> recv_request2(std::max(procs_overlapping_spparks_domain.size(), procs_overlapping_fem_domain.size()));
  int comm_tag = 100;

  // _console << "\ninitialSetup: 2D1\n";

  for (unsigned i = 0; i < procs_overlapping_spparks_domain.size(); ++i)
  //  if (num_fem_nodes.size() && procs_overlapping_spparks_domain.size())
       MPI_Irecv(&num_fem_nodes[i], 1, MPI_UNSIGNED, procs_overlapping_spparks_domain[i], comm_tag, _communicator.get(), &recv_request1[i]);

  // _console << "\ninitialSetup: 2D2\n";

  for (unsigned i = 0; i < procs_overlapping_fem_domain.size(); ++i)
  //  if (num_spparks_nodes.size() && procs_overlapping_fem_domain.size())
       MPI_Irecv(&num_spparks_nodes[i], 1, MPI_UNSIGNED, procs_overlapping_fem_domain[i], comm_tag+11, _communicator.get(), &recv_request2[i]);

  // _console << "\ninitialSetup: 2D3\n";

  for (unsigned i = 0; i < procs_overlapping_fem_domain.size(); ++i)
  //  if (procs_overlapping_fem_domain.size())
      MPI_Send(&_num_local_fem_nodes, 1, MPI_UNSIGNED, procs_overlapping_fem_domain[i], comm_tag, _communicator.get());

  // _console << "\ninitialSetup: 2D4\n";

  for (unsigned i = 0; i < procs_overlapping_spparks_domain.size(); ++i)
  //  if (procs_overlapping_spparks_domain.size())
      MPI_Send(&_num_local_spparks_nodes, 1, MPI_UNSIGNED, procs_overlapping_spparks_domain[i], comm_tag+11, _communicator.get());

  // _console << "\ninitialSetup: 2D5\n";

  std::vector<MPI_Status> recv_status1(std::max(procs_overlapping_spparks_domain.size(), procs_overlapping_fem_domain.size()));
  std::vector<MPI_Status> recv_status2(std::max(procs_overlapping_spparks_domain.size(), procs_overlapping_fem_domain.size()));
  MPI_Waitall(procs_overlapping_spparks_domain.size(), &recv_request1[0], &recv_status1[0]);
  MPI_Waitall(procs_overlapping_fem_domain.size(), &recv_request2[0], &recv_status2[0]);

  // _console << "\ninitialSetup: 2D\n";

  //
  // E: Communicate MOOSE FEM nodes, SPPARKS nodes
  //

  comm_tag = 200;
  int comm_tag_double = comm_tag + 1;
  const unsigned num_fem_nodes_total = std::accumulate(&num_fem_nodes[0], &num_fem_nodes[0]+num_fem_nodes.size(), 0);
  const unsigned num_spparks_nodes_total = std::accumulate(&num_spparks_nodes[0], &num_spparks_nodes[0]+num_spparks_nodes.size(), 0);
  std::vector<libMesh::dof_id_type> remote_fem_nodes(num_fem_nodes_total);
  std::vector<Real> remote_fem_coords(num_fem_nodes_total*3);
  std::vector<unsigned> remote_spparks_nodes(num_spparks_nodes_total);
  std::vector<Real> remote_spparks_coords(num_spparks_nodes_total*3);
  offset = 0;

  // Receive MOOSE FEM nodes
  std::vector<MPI_Request> recv_request_coor1(procs_overlapping_spparks_domain.size());
  for (unsigned i = 0; i < procs_overlapping_spparks_domain.size(); ++i)
  {
    MPI_Irecv(&remote_fem_nodes[offset], num_fem_nodes[i], MPI_UNSIGNED, procs_overlapping_spparks_domain[i], comm_tag, _communicator.get(), &recv_request1[i]);
    MPI_Irecv(&remote_fem_coords[offset*3], num_fem_nodes[i]*3, MPI_DOUBLE, procs_overlapping_spparks_domain[i], comm_tag_double, _communicator.get(), &recv_request_coor1[i]);
    offset += num_fem_nodes[i];
  }
  // Receive SPPARKS nodes
  std::vector<MPI_Request> recv_request_coor2(procs_overlapping_fem_domain.size());
  offset = 0;
  for (unsigned i = 0; i < procs_overlapping_fem_domain.size(); ++i)
  {
    MPI_Irecv(&remote_spparks_nodes[offset], num_spparks_nodes[i], MPI_UNSIGNED, procs_overlapping_fem_domain[i], comm_tag+11, _communicator.get(), &recv_request2[i]);
    MPI_Irecv(&remote_spparks_coords[offset*3], num_spparks_nodes[i]*3, MPI_DOUBLE, procs_overlapping_fem_domain[i], comm_tag_double+11, _communicator.get(), &recv_request_coor2[i]);
    offset += num_spparks_nodes[i];
  }


  // Prepare vectors of MOOSE FEM ids, coordinates
  std::vector<libMesh::dof_id_type> fem_ids(_num_local_fem_nodes);
  std::vector<Real> fem_coords(3*_num_local_fem_nodes);
  offset = 0;
  for (std::multiset<FEMID>::const_iterator i = fem_id.begin(); i != fem_id.end(); ++i)
  {
    fem_ids[offset] = i->id;
    fem_coords[offset*3+0] = i->coor(0);
    fem_coords[offset*3+1] = i->coor(1);
    fem_coords[offset*3+2] = i->coor(2);
    ++offset;
  }

  // Send MOOSE FEM ids, coordinates
  offset = 0;
  for (unsigned i = 0; i < procs_overlapping_fem_domain.size(); ++i) {
    MPI_Send(&fem_ids[0], _num_local_fem_nodes, MPI_UNSIGNED, procs_overlapping_fem_domain[i], comm_tag, _communicator.get());
    MPI_Send(&fem_coords[0], _num_local_fem_nodes*3, MPI_DOUBLE, procs_overlapping_fem_domain[i], comm_tag_double, _communicator.get());
    ++offset;
  }

  // Prepare SPPARKS ids, coordinates
  std::vector<unsigned> spparks_ids(_num_local_spparks_nodes);
  std::vector<Real> spparks_coords(3*_num_local_spparks_nodes);
  offset = 0;
  for (std::set<SPPARKSID>::const_iterator i = spparks_id.begin(); i != spparks_id.end(); ++i)
  {
    spparks_ids[offset] = i->id;
    spparks_coords[offset*3+0] = i->coor(0);
    spparks_coords[offset*3+1] = i->coor(1);
    spparks_coords[offset*3+2] = i->coor(2);
    ++offset;
  }

  // Send SPPARKS ids, coordinates
  offset = 0;
  for (unsigned i = 0; i < procs_overlapping_spparks_domain.size(); ++i) {
    MPI_Send(&spparks_ids[0], _num_local_spparks_nodes, MPI_UNSIGNED, procs_overlapping_spparks_domain[i], comm_tag+11, _communicator.get());
    MPI_Send(&spparks_coords[0], _num_local_spparks_nodes*3, MPI_DOUBLE, procs_overlapping_spparks_domain[i], comm_tag_double+11, _communicator.get());
    ++offset;
  }

  MPI_Waitall(procs_overlapping_spparks_domain.size(), &recv_request1[0], &recv_status1[0]);
  MPI_Waitall(procs_overlapping_spparks_domain.size(), &recv_request_coor1[0], &recv_status1[0]);
  MPI_Waitall(procs_overlapping_fem_domain.size(), &recv_request2[0], &recv_status2[0]);
  MPI_Waitall(procs_overlapping_fem_domain.size(), &recv_request_coor2[0], &recv_status2[0]);

  // _console << "\ninitialSetup: 2E \n";

  //
  // F: Count matching nodes for each proc that sent MOOSE FEM nodes, SPPARKS nodes
  //

  comm_tag = 300;
  // Receive number of MOOSE FEM nodes on this processor that match SPPARKS nodes on another
  std::vector<unsigned> num_remote_fem_matches(procs_overlapping_fem_domain.size());
  for (unsigned i = 0; i < procs_overlapping_fem_domain.size(); ++i)
  {
    MPI_Irecv(&num_remote_fem_matches[i], 1, MPI_UNSIGNED, procs_overlapping_fem_domain[i], comm_tag, _communicator.get(), &recv_request1[i]);
  }
  // Receive number of SPPARKS nodes on this processor that match MOOSE FEM nodes on another
  std::vector<unsigned> num_remote_spparks_matches(procs_overlapping_spparks_domain.size());
  for (unsigned i = 0; i < procs_overlapping_spparks_domain.size(); ++i)
  {
    MPI_Irecv(&num_remote_spparks_matches[i], 1, MPI_UNSIGNED, procs_overlapping_spparks_domain[i], comm_tag+11, _communicator.get(), &recv_request2[i]);
  }

  // Count number of remote MOOSE FEM nodes that match the processor's SPPARKS nodes
  // Store the matches
  std::vector<std::vector<unsigned> > spparks_matches(procs_overlapping_spparks_domain.size());
  offset = 0;
  for (unsigned i = 0; i < procs_overlapping_spparks_domain.size(); ++i)
  {
    for (unsigned j = 0; j < num_fem_nodes[i]; ++j)
    {
      SPPARKSID tmp(remote_fem_nodes[offset],
                     Point(remote_fem_coords[offset*3+0],
                           remote_fem_coords[offset*3+1],
                           remote_fem_coords[offset*3+2]));
      std::set<SPPARKSID>::iterator iter = spparks_id.find(tmp);
      if (iter != spparks_id.end())
      {
        spparks_matches[i].push_back(tmp.id);

        _spparks_to_proc[*iter].push_back(procs_overlapping_spparks_domain[i]);
      }
      ++offset;
    }
  }
  // Count number of remote SPPARKS nodes that match this processor's MOOSE FEM nodes
  // Store the matches
  std::vector<std::vector<unsigned> > fem_matches(procs_overlapping_fem_domain.size());
  offset = 0;
  for (unsigned i = 0; i < procs_overlapping_fem_domain.size(); ++i)
  {
    for (unsigned j = 0; j < num_spparks_nodes[i]; ++j)
    {
      FEMID tmp(remote_spparks_nodes[offset],
                 Point(remote_spparks_coords[offset*3+0],
                       remote_spparks_coords[offset*3+1],
                       remote_spparks_coords[offset*3+2]));
      std::multiset<FEMID>::iterator iter = fem_id.find(tmp);
      if (iter != fem_id.end())
      {
        fem_matches[i].push_back(tmp.id);

        _fem_to_proc[*iter].push_back(procs_overlapping_fem_domain[i]);
      }
      ++offset;
    }
  }

  // Send number of MOOSE FEM nodes that match this processor's SPPARKS nodes
  std::vector<unsigned> spparks_sizes(procs_overlapping_spparks_domain.size());
  for (unsigned i = 0; i < procs_overlapping_spparks_domain.size(); ++i)
  {
    spparks_sizes[i] = spparks_matches[i].size();
    MPI_Send(&spparks_sizes[i], 1, MPI_UNSIGNED, procs_overlapping_spparks_domain[i], comm_tag, _communicator.get());
  }
  // Send number of SPPARKS nodes that match this processor's MOOSE FEM nodes
  std::vector<unsigned> fem_sizes(procs_overlapping_fem_domain.size());
  for (unsigned i = 0; i < procs_overlapping_fem_domain.size(); ++i)
  {
    fem_sizes[i] = fem_matches[i].size();
    MPI_Send(&fem_sizes[i], 1, MPI_UNSIGNED, procs_overlapping_fem_domain[i], comm_tag+11, _communicator.get());
  }

  MPI_Waitall(procs_overlapping_fem_domain.size(), &recv_request1[0], &recv_status1[0]);
  MPI_Waitall(procs_overlapping_spparks_domain.size(), &recv_request2[0], &recv_status2[0]);

  // _console << "\ninitialSetup: 2F \n";

  //
  // G: Communicate matching nodes
  //

  comm_tag = 400;
  // Receive MOOSE DEM ids that match SPPARKS nodes on another processor
  std::vector<std::vector<unsigned> > matched_fem_ids(procs_overlapping_fem_domain.size());
  for (unsigned i = 0; i < procs_overlapping_fem_domain.size(); ++i)
  {
    matched_fem_ids[i].resize(num_remote_fem_matches[i]);
    if (matched_fem_ids[i].size()) MPI_Irecv(&matched_fem_ids[i][0], num_remote_fem_matches[i], MPI_UNSIGNED, procs_overlapping_fem_domain[i], comm_tag, _communicator.get(), &recv_request1[i]);
  }
  // Receive SPPARKS ids that match MOOSE FEM nodes on another processor
  std::vector<std::vector<unsigned> > matched_spparks_ids(procs_overlapping_spparks_domain.size());
  for (unsigned i = 0; i < procs_overlapping_spparks_domain.size(); ++i)
  {
    matched_spparks_ids[i].resize(num_remote_spparks_matches[i]);
    if (matched_spparks_ids[i].size()) MPI_Irecv(&matched_spparks_ids[i][0], num_remote_spparks_matches[i], MPI_UNSIGNED, procs_overlapping_spparks_domain[i], comm_tag+11, _communicator.get(), &recv_request2[i]);
  }

  // Send remote MOOSE FEM ids that match SPPARKS nodes on this processor
  for (unsigned i = 0; i < procs_overlapping_spparks_domain.size(); ++i)
  {
    spparks_matches[i].resize(spparks_sizes[i]);
    if (spparks_matches[i].size()) MPI_Send(&spparks_matches[i][0], spparks_sizes[i], MPI_UNSIGNED, procs_overlapping_spparks_domain[i], comm_tag, _communicator.get());
  }
  // Send remote SPPARKS ids that match MOOSE FEM nodes on this processor
  for (unsigned i = 0; i < procs_overlapping_fem_domain.size(); ++i)
  {
    fem_matches[i].resize(fem_sizes[i]);
    if (fem_matches[i].size()) MPI_Send(&fem_matches[i][0], fem_sizes[i], MPI_UNSIGNED, procs_overlapping_fem_domain[i], comm_tag+11, _communicator.get());
  }

  MPI_Waitall(procs_overlapping_fem_domain.size(), &recv_request1[0], &recv_status1[0]);
  MPI_Waitall(procs_overlapping_spparks_domain.size(), &recv_request2[0], &recv_status2[0]);

  // _console << "\ninitialSetup: 2G \n";

  //
  // H: Generate final recv communication maps
  //

  for (unsigned i = 0; i < procs_overlapping_fem_domain.size(); ++i)
  {
    for (unsigned j = 0; j < num_remote_fem_matches[i]; ++j)
    {
      _sending_proc_to_fem_id[procs_overlapping_fem_domain[i]].push_back(matched_fem_ids[i][j]);
    }
  }
  for (unsigned i = 0; i < procs_overlapping_spparks_domain.size(); ++i)
  {
    for (unsigned j = 0; j < num_remote_spparks_matches[i]; ++j)
    {
      _sending_proc_to_spparks_id[procs_overlapping_spparks_domain[i]].push_back(matched_spparks_ids[i][j]);
    }
  }

  //  _console << "\ninitialSetup: 2H\n";

  /*
  if (_init_spparks)
  {
    initSPPARKS();
    _init_spparks = false;

    setSPPARKSData();
  } */
}
