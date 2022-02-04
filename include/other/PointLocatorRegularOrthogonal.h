/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#pragma once

#include "Moose.h"
#include "libmesh/point.h"
#include "libmesh/point_locator_base.h"

namespace libMesh
{
class MeshBase;
}
class PointLocatorRegularOrthogonalData;

class PointLocatorRegularOrthogonal : public PointLocatorBase
{
public:
  /**
   * Constructor.  Needs the \p mesh in which the points
   * should be located.  Optionally takes a master
   * interpolator.  This master helps in saving memory
   * by reducing the number of root element tables in use.  Only the
   * master locator holds a root element table, the others simply
   * use the master's table.
   */
  PointLocatorRegularOrthogonal(const MeshBase & mesh, const PointLocatorBase * master = nullptr);

  virtual ~PointLocatorRegularOrthogonal();

  /**
   * Clears the locator.  This function frees dynamic memory with "delete" for the master.
   */
  virtual void clear() override;

  /// Do not use this init method, use the one below.
  virtual void init() override;

  /**
   * Initializes the locator, so that the \p operator() methods can
   * be used. This method allocates dynamic memory with "new".
   */
  void init(const std::vector<unsigned int> & cell_count,
            const Point & min_corner,
            const Point & max_corner);

  /**
   * Locates the element in which the point with global coordinates
   * \p p is located, optionally restricted to a set of allowed subdomains.
   */
  virtual const Elem *
  operator()(const Point & p,
             const std::set<subdomain_id_type> * allowed_subdomains = nullptr) const override;

  /**
   * Locates a set of elements in proximity to the point with global coordinates
   * \p p  Pure virtual. Optionally allows the user to restrict the subdomains searched.
   */
  virtual void
  operator()(const Point & p,
             std::set<const Elem *> & candidate_elements,
             const std::set<subdomain_id_type> * allowed_subdomains = nullptr) const override;

  /**
   * Enables out-of-mesh mode.  In this mode, if asked to find a point
   * that is contained in no mesh at all, the point locator will
   * return a NULL pointer instead of crashing.  Per default, this
   * mode is off.
   */
  virtual void enable_out_of_mesh_mode() override { _out_of_mesh_mode = true; }

  /**
   * Disables out-of-mesh mode (default).  If asked to find a point
   * that is contained in no mesh at all, the point locator will now
   * crash.
   */
  virtual void disable_out_of_mesh_mode() override { _out_of_mesh_mode = false; }

protected:
  /// true if out-of-mesh mode is enabled
  bool _out_of_mesh_mode;

  /// internal data object, shared between master and servants
  PointLocatorRegularOrthogonalData * _data;
};
