/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "RadialGreensConvolution.h"
#include "libmesh/nanoflann.hpp"

#include <list>

registerMooseObject("MagpieApp", RadialGreensConvolution);

// serialization helper for parallel communication
template <>
void
dataStore(std::ostream & stream, RadialGreensConvolution::QPData & qpd, void * context)
{
  dataStore(stream, qpd._q_point, context);
  dataStore(stream, qpd._elem_id, context);
  dataStore(stream, qpd._qp, context);
  dataStore(stream, qpd._volume, context);
  dataStore(stream, qpd._value, context);
}

// unserialization helper for parallel communication
template <>
void
dataLoad(std::istream & stream, RadialGreensConvolution::QPData & qpd, void * context)
{
  dataLoad(stream, qpd._q_point, context);
  dataLoad(stream, qpd._elem_id, context);
  dataLoad(stream, qpd._qp, context);
  dataLoad(stream, qpd._volume, context);
  dataLoad(stream, qpd._value, context);
}

// specialization for PointListAdaptor<RadialGreensConvolution::QPData>
template <>
inline const Point &
PointListAdaptor<RadialGreensConvolution::QPData>::getPoint(const size_t idx) const
{
  return _pts[idx]._q_point;
}

template <>
InputParameters
validParams<RadialGreensConvolution>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addClassDescription("Perform a radial Green's function convolution");
  params.addCoupledVar("v", "Variable to gather");
  params.addRequiredParam<FunctionName>("function",
                                        "Green's function (distance is substituted for x) without geometrical attenuation");
  params.addRequiredParam<Real>("r_cut", "Cut-off radius for the Green's function");
  params.addParam<bool>("normalize", false, "Normalize the Green's function integral to one");

  // we run this object once at the beginning of the timestep by default
  params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;

  return params;
}

RadialGreensConvolution::RadialGreensConvolution(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _v(coupledValue("v")),
    _v_var(coupled("v")),
    _function(getFunction("function")),
    _r_cut(getParam<Real>("r_cut")),
    _normalize(getParam<bool>("normalize")),
    _dim(_mesh.dimension()),
    _correction_integral(100)
{
  // collect mesh periodicity data
  for (unsigned int i = 0; i < _dim; ++i)
  {
    _periodic[i] = _mesh.isRegularOrthogonal() && _mesh.isTranslatedPeriodic(_v_var, i);

    _periodic_min[i] = _mesh.getMinInDimension(i);
    _periodic_max[i] = _mesh.getMaxInDimension(i);
    _periodic_vector[i](i) = _mesh.dimensionWidth(i);

    // we could allow this, but then we'd have to search over more than just the nearest periodic
    // neighbors
    if (_periodic[i] && 2.0 * _r_cut > _periodic_vector[i](i))
      paramError(
          "r_cut",
          "The cut-off radius cannot be larger than half the periodic size of the simulation cell");
  }
}

void
RadialGreensConvolution::initialize()
{
  _qp_data.clear();

  /* compute the dimensional correction
   *
   *  |---__
   *  |     / .
   *  |  R/  | .
   *  | /   h|  \ r_cut
   *  |______|__|____
   *         r
   *
   * Out of plane (2D) and out of line (1D) 3D contributions to the convolution
   * are not naturally captured by the 1D or 2D numerical integration of QP neighborhood.
   * At every Qp an integration in the h direction needs to be performed.
   * We can pre-tabulate this integral.
   */
  if (_dim < 3)
  {
    const Real dr = _r_cut / _correction_integral.size();
    for (unsigned int i = 0; i < _correction_integral.size(); ++i)
    {
      const Real r = i * dr;
      const Real h = std::sqrt(_r_cut * _r_cut - r * r);
      const unsigned int hsteps = std::ceil(h/dr);
      const Real dh = h / hsteps;

      _correction_integral[i] = 0.0;

      // for r = 0 we skip the first bin of the integral!
      if (i == 0)
        _zero_dh = dh;

      for (unsigned int j = i > 0 ? 0 : 1; j < hsteps; ++j)
      {
        const Real h1 = j * dh;
        const Real h2 = h1 + dh;

        // we evaluate the Green's function at the geometric middle of the height interval
        // and assume it to be constant.
        const Real R2 = r * r + h1 * h2;
        const Real G = _function.value(_t, Point(std::sqrt(R2), 0.0, 0.0));

        // We integrate the geometric attenuation analytically over the interval [h1, h2).
        _correction_integral[i] += G * attenuationIntegral(h1, h2, r, _dim);
      }
    }
  }
}

Real
RadialGreensConvolution::attenuationIntegral(Real h1, Real h2, Real r, unsigned int dim) const
{
  if (dim == 2)
  {
    // in 2D we need to return the above and below plane contribution (hence 2.0 * ...)
    if (r > 0.0)
      // integral 1/(h^2+r^2) dh = 1/r*atan(h/r)|
      return 2.0 * 1.0/(4.0 * libMesh::pi * r) * (std::atan(h2/r) - std::atan(h1/r));
    else
      // integral 1/h^2 dh = -1/h|
      return 2.0 * 1.0/(4.0 * libMesh::pi) * (-1.0/h2 + 1.0/h1);
  }
  else
    // dim = 1, we multiply the attenuation by 2pi*h (1/2 * 2pi/4pi = 0.25)
    // integral h/(h^2+r^2) dh = 1/2 * ln(h^2+r^2)|
    return 0.25 * (std::log(h2*h2 + r*r) - std::log(h1*h1 + r*r));
}

void
RadialGreensConvolution::execute()
{
  auto id = _current_elem->id();

  // collect all QP data
  for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
    _qp_data.emplace_back(_q_point[qp], id, qp, _JxW[qp] * _coord[qp], _v[qp]);

  // make sure the result map entry for the current element is sized correctly
  auto i = _convolution.find(id);
  if (i == _convolution.end())
    _convolution.insert(std::make_pair(id, std::vector<Real>(_qrule->n_points())));
  else
    i->second.resize(_qrule->n_points());
}

void
RadialGreensConvolution::finalize()
{
  _console <<  "RadialGreensConvolution::finalize()" << std::flush;

  // the first chunk of data is always the local data - remember its size
  unsigned int local_size = _qp_data.size();

  // communicate the qp data list if n_proc > 1
  if (_app.n_processors() > 1)
  {
    // pack the local qp data into a string buffer
    std::string send_buffer;
    std::ostringstream oss;
    dataStore(oss, _qp_data, this);
    send_buffer.assign(oss.str());

    // create byte buffers for the streams received from all processors
    std::vector<std::string> recv_buffers;

    // broadcast serialized data to and receive from all processors
    _communicator.allgather(send_buffer, recv_buffers);
    mooseAssert(recv_buffers.size() == _app.n_processors(),
                "Unexpected size of recv_buffers: " << recv_buffers.size());

    // Loop over all data structures for all processors to perform the gather operation
    std::istringstream iss;
    for (unsigned int rank = 0; rank < recv_buffers.size(); ++rank)
    {
      // skip the current processor (its data is already in the structures)
      if (rank == processor_id())
        continue;

      // populate the stream with a new buffer and reset stream state
      iss.str(recv_buffers[rank]);
      iss.clear();

      // Load the communicated data into temporary structures
      std::vector<QPData> other_data;
      dataLoad(iss, other_data, this);

      // merging the qp data lists (by appending to the end to keep local data at the beginning of
      // the vector)
      _qp_data.insert(_qp_data.end(), other_data.begin(), other_data.end());
    }
  }

  _console <<  "  communiction done." << std::flush;

  // if normalization is requested we need to integrate the current variable field
  Real _source_integral = 0.0;
  if (_normalize)
    for (const auto & qpd : _qp_data)
      _source_integral += qpd._volume * qpd._value;

  // build KD-Tree using data we just allgathered
  const unsigned int max_leaf_size = 20; // slightly affects runtime
  auto point_list = PointListAdaptor<QPData>(_qp_data);
  auto kd_tree = libmesh_make_unique<KDTreeType>(
      LIBMESH_DIM, point_list, nanoflann::KDTreeSingleIndexAdaptorParams(max_leaf_size));

  mooseAssert(kd_tree != nullptr, "KDTree was not properly initialized.");
  kd_tree->buildIndex();

  _console <<  "  KD-tree built." << std::flush;

  // result map entry
  const auto end_it = _convolution.end();
  auto it = end_it;

  // integral of the convolution (used if normalization is requested)
  Real _convolution_integral = 0.0;

  // radial bin size
  const Real dr = _r_cut / _correction_integral.size();

  // iterate over the local portion of the gathered QP data (TODO: make this threaded)
  std::vector<std::pair<std::size_t, Real>> ret_matches;
  nanoflann::SearchParams search_params;

  _console <<  "  iterating over " << local_size << " QPs..." << std::flush;
  const int interval = local_size / 100;

  for (unsigned int i = 0; i < local_size; ++i)
  {
    if (_dim == 3 && i % interval == 0)
      _console <<  "  " << i/interval << "%..." << std::flush;

    const auto & local_qp = _qp_data[i];

    // Look up result map iterator only if we enter a new element. this saves a bunch
    // of map lookups because same element entries are consecutive in the _qp_data vector.
    if (it == end_it || it->first != local_qp._elem_id)
      it = _convolution.find(local_qp._elem_id);

    // initialize result entry
    mooseAssert(it != end_it, "Current element id not found in result set.");
    auto & sum = it->second[local_qp._qp];
    sum = 0.0;

    // if the variable is periodic we need to perform extra searches translated onto
    // the periodic neighbors
    std::list<Point> cell_vector = {Point()};
    for (unsigned int j = 0; j < _dim; ++j)
      if (_periodic[j])
      {
        std::list<Point> new_cell_vector;

        for (const auto & cell : cell_vector)
        {
          if (local_qp._q_point(j) + _periodic_vector[j](j) - _r_cut < _periodic_max[j])
            new_cell_vector.push_back(cell + _periodic_vector[j]);

          if (local_qp._q_point(j) - _periodic_vector[j](j) + _r_cut > _periodic_min[j])
            new_cell_vector.push_back(cell - _periodic_vector[j]);
        }

        cell_vector.insert(cell_vector.end(), new_cell_vector.begin(), new_cell_vector.end());
      }

    // perform radius search and aggregate data considering potential periodicity
    Point center;
    for (const auto & cell : cell_vector)
    {
      ret_matches.clear();
      center = local_qp._q_point + cell;
      std::size_t n_result =
          kd_tree->radiusSearch(&(center(0)), _r_cut * _r_cut, ret_matches, search_params);
      for (std::size_t j = 0; j < n_result; ++j)
      {
        const auto & other_qp = _qp_data[ret_matches[j].first];
        const Real r = std::sqrt(ret_matches[j].second);

        // R is the equivalent sphere radius for the quadrature point. The spherical integral
        // integral_0^R 1/(4pi*r^2) *4pi*r^2 = R
        // So R is the integral over the geometric attenuation at the center quadrature point
        switch (_dim)
        {
          case 1:
          {
            // correction integral is the integral over the geometric attenuation
            // times the Green's function over a 2D disc inscribed into the r_cut
            // sphere and perpendicular to the mesh.
            Real add = _correction_integral[std::floor(r/dr)] * other_qp._volume;

            if (r == 0)
            {
              const Real R = 0.5 * other_qp._volume;
              add += _function.value(_t, Point(0.0, 0.0, 0.0)) * (
                // add the center sphere attenuation integral
                R +
                // add the section missing or overlapping between _correction_integral and center sphere
                attenuationIntegral(R, _zero_dh, 0.0, 1) * other_qp._volume);
            }
            sum += add * other_qp._value;
            break;
          }

          case 2:
          {
            // correction integral is the integral over the geometric attenuation
            // times the Green's function over a 1D line segment inscribed into the r_cut
            // sphere and perpendicular to the mesh.
            Real add = _correction_integral[std::floor(r/dr)] * other_qp._volume;
            if (r == 0)
            {
              const Real R = std::sqrt(other_qp._volume / (2.0 * libMesh::pi));
              add += _function.value(_t, Point(0.0, 0.0, 0.0)) * (
                // add the center sphere attenuation integral
                R +
                // add the section missing or overlapping between _correction_integral and center sphere
                attenuationIntegral(R, _zero_dh, 0.0, 2) * other_qp._volume);
            }
            sum += add * other_qp._volume * other_qp._value;
            break;
          }

          case 3:
          {
            const Real G = _function.value(_t, Point(r, 0.0, 0.0));
            if (r == 0)
            {
              // R is the integral over the geometric attenuation in a sphere around the origin
              const Real R = std::cbrt(3.0/4.0 * other_qp._volume / libMesh::pi);
              sum += G * R * other_qp._value;
            }
            else
              sum += G * 0.25 / (libMesh::pi * r * r) * other_qp._volume * other_qp._value;
          }
        }
      }
    }

    // integrate the convolution result
    _convolution_integral += sum;
  }

  // if normalization is requested we need to communicate the convolution result
  // and normalize the result entries
  if (_normalize)
  {
    gatherSum(_convolution_integral);

    // we may not need to or may not be able to normalize
    if (_convolution_integral == 0.0)
    {
      if (_source_integral == 0.0)
        return;
      mooseError("Unable to normalize Green's function. Is it all zero?");
    }

    // normalize result entries
    for (auto & re : _convolution)
      for (auto & ri : re.second)
        ri *= _source_integral / _convolution_integral;
  }

  // make it a differential result
  it = end_it;
  for (unsigned int i = 0; i < local_size; ++i)
  {
    const auto & local_qp = _qp_data[i];

    // Look up result map iterator only if we enter a new element. this saves a bunch
    // of map lookups because same element entries are consecutive in the _qp_data vector.
    if (it == end_it || it->first != local_qp._elem_id)
      it = _convolution.find(local_qp._elem_id);

    it->second[local_qp._qp] -= local_qp._volume * local_qp._value;
  }
}

void
RadialGreensConvolution::threadJoin(const UserObject & y)
{
  const RadialGreensConvolution & uo = static_cast<const RadialGreensConvolution &>(y);
  _qp_data.insert(_qp_data.begin(), uo._qp_data.begin(), uo._qp_data.end());
  _convolution.insert(uo._convolution.begin(), uo._convolution.end());
}
