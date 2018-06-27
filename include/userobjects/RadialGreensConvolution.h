/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#ifndef RADIALGREENSCONVOLUTION_H
#define RADIALGREENSCONVOLUTION_H

#include "ElementUserObject.h"
#include "DataIO.h"

class RadialGreensConvolution;

template <>
InputParameters validParams<RadialGreensConvolution>();

/**
 * Gather and communicate a full list of all quadrature points and the values of
 * a selected variable at each point. Use a KD-Tree to integrate the weighted
 * neighborhood of each QP to obtain the convolution.
 */
class RadialGreensConvolution : public ElementUserObject
{
public:
  RadialGreensConvolution(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;
  virtual void threadJoin(const UserObject & y) override;

  /// quaddrature point data
  struct QPData
  {
    /// physical coordinates of the quadrature point
    Point _q_point;
    /// element id
    dof_id_type _elem_id;
    /// index of the quadrature point
    short _qp;
    /// current value * _JxW
    Real _integral;

    QPData() : _q_point(), _elem_id(libMesh::invalid_uint), _qp(0), _integral(0.0) {}
    QPData(const Point & q_point, dof_id_type elem_id, short qp, Real integral)
      : _q_point(q_point), _elem_id(elem_id), _qp(qp), _integral(integral)
    {
    }
  };

  using Result = std::map<dof_id_type, std::vector<Real>>;
  const Result & getConvolution() const { return _convolution; }

protected:
  /// variable field to be gathered
  const VariableValue & _v;

  /// index of field variable
  unsigned int _v_var;

  /// Green's function
  Function & _function;

  /// Green's function cut-off radius
  const Real _r_cut;

  /// Normalize the Green's function to one to make the integral of teh convolution
  /// the same as the integral of the original data.
  const bool _normalize;

  /// gathered data
  std::vector<QPData> _qp_data;

  /// convolution result
  Result _convolution;

  /// is the mesh translated periodic in a given cardinal direction
  std::array<bool, LIBMESH_DIM> _periodic;

  ///@{ periodic size per component
  std::array<Real, LIBMESH_DIM> _periodic_min;
  std::array<Real, LIBMESH_DIM> _periodic_max;
  std::array<Point, LIBMESH_DIM> _periodic_vector;
  ///@}

  using KDTreeType = nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<Real, PointListAdaptor<QPData>>,
      PointListAdaptor<QPData>,
      LIBMESH_DIM>;
};

template <>
void dataStore(std::ostream &, RadialGreensConvolution::QPData &, void *);

template <>
void dataLoad(std::istream &, RadialGreensConvolution::QPData &, void *);

#endif // RADIALGREENSCONVOLUTION_H
