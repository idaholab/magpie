#include "MooseMyTRIMSample.h"
#include "MooseMesh.h"

#include "mytrim/element.h"

MooseMyTRIMSample::MooseMyTRIMSample(const MyTRIMRasterizer & rasterizer, const MooseMesh & mesh) :
    _rasterizer(rasterizer),
    _nvars(_rasterizer.nVars()),
    _trim_mass(_rasterizer.mass()),
    _trim_charge(_rasterizer.charge()),
    _mesh(mesh)
{
}

void
MooseMyTRIMSample::initialize()
{
  _pl = _mesh.getMesh().sub_point_locator();
  _materials_cache.clear();
}

MyTRIM_NS::materialBase *
MooseMyTRIMSample::lookupMaterial(double * pos)
{
  // point to sample the material at
  Point p(pos[0], pos[1], pos[2]);

  // get element containing the point
  const Elem * elem = (*_pl)(p);

  // try to find the element in the cache
  MaterialsCache::iterator i = _materials_cache.find(elem);
  if (i != _materials_cache.end())
    return &i->second;

  // otherwise prepare the material using data from the rasterizer
  const std::vector<Real> & material_data = _rasterizer.material(elem);
  MyTRIM_NS::materialBase material(1.0);

  // set elements
  for (unsigned int i = 0; i < _nvars; ++i)
  {
    MyTRIM_NS::elementBase * element = new MyTRIM_NS::elementBase;
    element->z = _trim_mass[i];
    element->m = _trim_charge[i];
    element->t = material_data[i];
    material.element.push_back(element);
  }

  material.prepare();

  // add the material to the cache and return pointer
  MaterialsCache::iterator mi = _materials_cache.insert(_materials_cache.begin(), std::make_pair(elem, material));
  return &mi->second;
}
