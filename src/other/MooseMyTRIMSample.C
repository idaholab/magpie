/**********************************************************************/
/*                     DO NOT MODIFY THIS HEADER                      */
/* MAGPIE - Mesoscale Atomistic Glue Program for Integrated Execution */
/*                                                                    */
/*            Copyright 2017 Battelle Energy Alliance, LLC            */
/*                        ALL RIGHTS RESERVED                         */
/**********************************************************************/

#include "MooseMyTRIMSample.h"
#include "MooseMesh.h"

MooseMyTRIMSample::MooseMyTRIMSample(const MyTRIMRasterizer & rasterizer,
                                     const MooseMesh & mesh,
                                     MyTRIM_NS::SimconfType * simconf)
  : MyTRIM_NS::SampleBase(),
    _rasterizer(rasterizer),
    _trim_parameters(_rasterizer.getTrimParameters()),
    _nvars(_trim_parameters.nVars()),
    _mesh(mesh),
    _dim(_mesh.dimension()),
    _pl(_mesh.getPointLocator()),
    _simconf(simconf)
{
  if (_dim < 2 || _dim > 3)
    mooseError("TRIM simulation works in 2D or 3D only.");

  // permit querying points that are potentially outside the mesh
  _pl->enable_out_of_mesh_mode();
}

void
MooseMyTRIMSample::averages(const MyTRIM_NS::IonBase * pka)
{
  _current_ion = pka;

  // averages do NOT depend on ion energy! Cache averaged materials for all possible PKAs
  // this means adding comparison operators to PKAs so that we can use them as map indices
  // have a master material cache and per-PKA-type cache. Do not average stuff in the master
  // cache at all use it only as a copy construction source for the per-pka-cache
  // std::map<MyTRIM_NS::IonBase, std::map<Elem *, MooseMyTRIMMaterial> _per_pka_cache
}

MyTRIM_NS::MaterialBase *
MooseMyTRIMSample::lookupMaterial(Point & pos)
{
  // point to sample the material at
  Point p = _rasterizer.periodicPoint(pos);

  // get element containing the point
  mooseAssert(_pl != nullptr, "initialize() must be called on the MooseMyTRIMSample object.");
  const Elem * elem = (*_pl)(p);

  // no element found means we have left the mesh
  if (elem == nullptr)
    return nullptr;

  // obtain the per pka cache entry (insert if not found)
  auto & cache = _per_pka_materials_cache[*_current_ion];

  // look up current element
  auto i = cache.find(elem);

  // not averaged yet
  if (i == cache.end())
  {
    // look up if the element is prepared in the master cache
    auto j = _materials_master_cache.find(elem);

    // not prepared yet
    if (j == _materials_master_cache.end())
    {
      // prepare the material using data from the rasterizer
      const std::vector<Real> & material_data = _rasterizer.material(elem);
      MooseMyTRIMMaterial material(_simconf);

      // set elements
      MyTRIM_NS::Element element;
      for (unsigned int i = 0; i < _nvars; ++i)
      {
        element = _trim_parameters.element_prototypes[i];
        element._t = material_data[i];
        material._element.push_back(element);
      }

      // calculate the density (must happen first!)
      material.calculateDensity(_rasterizer.siteVolume(elem));

      // prepare material
      material.prepare();
      j = _materials_master_cache.insert(_materials_master_cache.begin(),
                                         std::make_pair(elem, material));
    }

    // create a copy from the master cache entry, average it and file it
    MooseMyTRIMMaterial material(j->second);
    material.average(_current_ion);

    i = cache.insert(cache.begin(), std::make_pair(elem, material));
  }

  return &i->second;
}
